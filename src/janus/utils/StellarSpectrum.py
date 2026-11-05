import numpy as np
import re
import shutil , os
import subprocess
import tempfile
import threading
from collections import deque
from scipy.interpolate import PchipInterpolator

import logging 
log = logging.getLogger("fwl."+__name__)


from .. import set_socrates_env 

PREP_SPEC_TIMEOUT = 1800.0      # [s] limit of one prep_spec call; 4096 bands take about 50 s
MAX_OUTPUT_FOLDER_BYTES = 190   # With the trailing separator
STAR_DATA_MARKERS = ("*BEGIN_DATA", "*BEGIN_BINNED_DATA")
BLOCK_HEADER = re.compile(r"\*BLOCK: TYPE =\s*(\d+)")

def PrepareStellarSpectrum(wl, fl, star_file, nbins_max=95000):
    """Write a stellar spectrum.

    This function supplements InsertStellarSpectrum by writing a stellar 
    spectrum to a file in the format that SOCRATES expects. The flux needs to 
    be provided at 1 AU. Bins the spectrum down to at most nbins_max bins.

    Parameters
    ----------
        wl : list
            Wavelength [nm]
        fl : list
            Flux at 1 AU [erg s-1 cm-2 nm-1]
        star_file : str
            Path to output file, which will contain the stellar spectrum.
        nbins_max : int
            Number of spectral bins to down-sample to (-1 for no rebinning)
            
    """

    socrates_nbins_max = int(1e5 - 1)

    # Validate
    if (len(wl) != len(fl)):
        raise Exception("Stellar wavelength and flux arrays have different lengths")
    if (len(wl) < 500):
        log.warning("Loaded stellar spectrum is very short!")

    nbins_max = min(len(wl), nbins_max)
    
    if (nbins_max > socrates_nbins_max):
        raise Exception("Too many bins requested for stellar spectrum (maximum is %d)" % socrates_nbins_max)

    # Down-sample spectrum when necessary or requested
    if len(wl) > socrates_nbins_max:
        log.info("Rebinning stellar spectrum")

        # Store old wl,fl arrays
        wl_orig = wl
        fl_orig = fl

        if (nbins_max < 500):
            log.warning("Requested number of bins is small (%d bins)" % nbins_max)

        nbins_max = min( int(socrates_nbins_max), nbins_max) # Must be fewer than 100k

        # create interpolator on log-wavelength grid
        itp = PchipInterpolator(np.log10(wl_orig), fl_orig)
        wl = np.linspace(np.log10(wl_orig[0]), np.log10(wl_orig[-1]), nbins_max)

        # interpolate new fluxes
        fl = itp(wl)
        wl = 10**wl
            
    # Convert units
    wl = np.array(wl) * 1.0e-9  # [nm] -> [m]
    fl = np.array(fl) * 1.0e6   # [erg s-1 cm-2 nm-1] -> [W m-3]

    # Store header
    content = ""
    content += "Star spectrum at 1 AU. Created using PrepareStellarSpectrum() \n"
    content += "      WAVELENGTH        IRRADIANCE\n"
    content += "          (m)               (W/m3)\n"
    content += "*BEGIN_DATA\n"

    # Store body of data
    for i in range(len(wl)):
        content += str("      %1.7e      %1.7e\n" % (wl[i],fl[i]))

    # Store footer
    content += "*END\n"
    content += " "

    # Write content to file
    with open(star_file,'w') as handle:
        handle.write(content)


def InsertStellarSpectrum(orig_file:str, star_file:str, output_folder:str):
    """Insert a stellar spectrum.

    It's nice to be able to switch out the stellar spectrum for a different one. 
    This function takes in an original spectral file, with opacity data, and 
    inserts a stellar spectrum into a copy of it.

    Parameters
    ----------
        orig_file : str
            Path to original spectral file WITHOUT stellar spectrum.
        star_file : str
            Path to file containing stellar spectrum in the SOCRATES format.
        output_folder : str
            Path to output folder, at most MAX_OUTPUT_FOLDER_BYTES bytes with its
            trailing separator.

    Raises
    ------
        ValueError
            If the output folder path is too long, or the stellar spectrum has no data.
        RuntimeError
            If prep_spec is missing, fails, passes PREP_SPEC_TIMEOUT or writes no stellar
            spectrum. No star.sf remains in the output folder after a failure in prep_spec.
    """

    # Cl_run_cdf later holds the folder and "mcica_data" in CHARACTER (LEN=200) MCICA_DATA.
    folder_bytes = len(os.fsencode(os.path.join(output_folder, "")))
    if folder_bytes > MAX_OUTPUT_FOLDER_BYTES:
        raise ValueError(f"Output folder path has {folder_bytes} bytes, above the SOCRATES "
                         f"limit of {MAX_OUTPUT_FOLDER_BYTES}: {output_folder}")

    # prep_spec reads on without end if the star file has no data marker
    has_data = False
    if os.path.isfile(star_file):
        with open(star_file, errors="replace") as handle:
            has_data = any(line.startswith(STAR_DATA_MARKERS) for line in handle)
    if not has_data:
        raise ValueError(f"Stellar spectrum has no {STAR_DATA_MARKERS[0]} line or is missing: "
                         f"{star_file}")

    prep_spec = shutil.which("prep_spec")
    if prep_spec is None:
        raise RuntimeError("Cannot find prep_spec on PATH: check the SOCRATES build in RAD_DIR")

    outp_file  = os.path.join(output_folder, "star.sf")
    outp_filek = outp_file+"_k"

    # Delete "new" files if they already exist
    if os.path.exists(outp_file):
        os.remove(outp_file)
    if os.path.exists(outp_filek):
        os.remove(outp_filek)

    # Relative names avoid the 132-character limit of prep_spec on any file name; the lock
    # file stops its interactive mode, which loops on a bad answer.
    with tempfile.TemporaryDirectory(prefix="prep_spec_", dir=output_folder) as work:
        shutil.copyfile(orig_file,      os.path.join(work, "star.sf"))
        shutil.copyfile(orig_file+"_k", os.path.join(work, "star.sf_k"))
        shutil.copyfile(star_file,      os.path.join(work, "star.txt"))
        open(os.path.join(work, "radiation_code.lock"), "w").close()

        inputs = ["star.sf",'a',                     # append existing file
                  '6','n','T','100 4000','250',      # tabulate thermal source function
                  '2','n',"star.txt",'y',            # insert stellar spectrum
                  '-1','EOF'                         # save and exit
                  ]
        answers = os.path.join(work, "answers.txt")
        with open(answers, "w") as handle:
            handle.write('\n'.join(inputs))

        # stderr holds the error line of a star that prep_spec rejects; stdout holds the menus
        # and some other error text, and is dropped
        with open(answers) as stdin, subprocess.Popen(
                [os.path.abspath(prep_spec)], cwd=work, stdin=stdin, stdout=subprocess.DEVNULL,
                stderr=subprocess.PIPE, encoding="ascii", errors="replace") as proc:
            timer = threading.Timer(PREP_SPEC_TIMEOUT, proc.kill)
            timer.daemon = True
            timer.start()
            tail = "".join(deque(proc.stderr, maxlen=20)).strip()
            expired = not timer.is_alive()
            timer.cancel()
        if expired:
            raise RuntimeError(f"prep_spec passed its time limit, {PREP_SPEC_TIMEOUT:g} s: "
                               f"{tail}")
        if proc.returncode != 0:
            raise RuntimeError(f"prep_spec returned with code {proc.returncode}: {tail}")

        # prep_spec exits with 0 when it leaves out a block after an error
        with open(os.path.join(work, "star.sf"), errors="replace") as handle:
            blocks = {int(m[1]) for m in map(BLOCK_HEADER.match, handle) if m}
        if not {2, 6} <= blocks:
            raise RuntimeError("prep_spec wrote no stellar spectrum (block 2) or thermal "
                               f"source function (block 6): {tail}")

        os.replace(os.path.join(work, "star.sf"),   outp_file)
        os.replace(os.path.join(work, "star.sf_k"), outp_filek)

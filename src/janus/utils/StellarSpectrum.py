import numpy as np
import shutil , os
import subprocess
from scipy.interpolate import PchipInterpolator

from janus.utils.logs import GetLogger
log = GetLogger()

from .. import set_socrates_env 

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
            Path to output folder
            
    """

    # k files
    orig_filek = orig_file+"_k"

    outp_file  = os.path.join(output_folder, "star.sf")
    outp_filek = outp_file+"_k"

    # Delete "new" files if they already exist
    if os.path.exists(outp_file):
        os.remove(outp_file)
    if os.path.exists(outp_filek):
        os.remove(outp_filek)

    # Copy original files to new location (retain old files)
    shutil.copyfile(orig_file,  outp_file)
    shutil.copyfile(orig_filek, outp_filek)

    # Run prep_spec
    inputs = [outp_file,'a',                    # append existing file 
              '6','n','T','100 4000','250',     # tabulate thermal source function
              '2','n',star_file,'y',            # insert stellar spectrum
              '-1','EOF'                        # save and exit
              ]
    p = subprocess.run(['prep_spec'], stdout=subprocess.PIPE, input='\n'.join(inputs), encoding='ascii')
    if (p.returncode != 0):
        log.warning("prep_spec returned with code %d" % p.returncode)
    


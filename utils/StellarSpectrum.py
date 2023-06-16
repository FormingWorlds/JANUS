#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 15 16:07:20 2023

@author: nichollsh
"""

import numpy as np
import shutil , os
import subprocess

def PrepareStellarSpectrum(wl, fl, star_file:str):
    """Write a stellar spectrum.

    This function supplements InsertStellarSpectrum by writing a stellar 
    spectrum to a file in the format that SOCRATES expects. The flux needs to 
    be provided at 1 AU.

    Parameters
    ----------
        wl : list
            Wavelength [nm]
        fl : list
            Flux at 1 AU [erg s-1 cm-2 nm-1]
        star_file : str
            Path to output file, which will contain the stellar spectrum.
            
    """

    # Convert units
    wl = np.array(wl) * 1.0e-9  # [nm] -> [m]
    fl = np.array(fl) * 1.0e6   # [erg s-1 cm-2 nm-1] -> [W m-3]

    # Validate
    if (len(wl) != len(fl)):
        raise Exception("Wavelength and Flux arrays have different lengths")

    if (len(wl) < 2000):
        print("WARNING: Stellar spectrum array is very short!")

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


def InsertStellarSpectrum(orig_file:str, star_file:str, outp_file:str):
    """Insert a stellar spectrum.

    It's nice to be able to switch out the stellar spectrum for a different one. 
    This function takes in an original spectra file, with opacity data, and 
    inserts a stellar spectrum into a copy of it.

    Parameters
    ----------
        orig_file : str
            Path to original spectral file WITHOUT stellar spectrum.
        star_file : str
            Path to file containing stellar spectrum in the SOCRATES format.
        outp_file : str
            Path to output file, containing both opacity and stellar data.
            
    """

    # k files
    orig_filek = orig_file+"_k"
    outp_filek = outp_file+"_k"

    # Delete "new" files if they already exist
    if os.path.exists(outp_file):
        os.remove(outp_file)
    if os.path.exists(outp_filek):
        os.remove(outp_filek)

    # Copy original files to new location (retain old files)
    shutil.copyfile(orig_file,  outp_file)
    shutil.copyfile(orig_filek, outp_filek)

    # Run prep_spec from SOCRATES
    inputs = [outp_file,'a','6','n','T','100 4000','100','2','n',star_file,'y','-1','EOF']
    p = subprocess.run(['prep_spec'], stdout=subprocess.PIPE, input='\n'.join(inputs), encoding='ascii')
    if (p.returncode != 0):
        print("WARNING: prep_spec returned with code %d" % p.returncode)
    


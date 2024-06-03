import os
import numpy as np

# Read band information from spectral file
def ReadBandEdges(spfile:str):

    if not os.path.exists(spfile):
        raise Exception("Canot find spectral file '%s'"%spfile)

    with open(spfile,'r') as hdl:
        lines = hdl.readlines()
    band_edgesm = []
    block_idx:int = -999999999
    for l in lines:
        # We want block 1 data
        if "Band        Lower limit         Upper limit" in l:
            block_idx = 0 
        if (block_idx > 0) and ("*END" in l):
            break 
        
        # Read bands
        if (block_idx > 0):
            s = [float(ss.strip()) for ss in l.split()[1:]]
            # First edge
            if block_idx == 1:
                band_edgesm.append(s[0])  # Note that these are in units of [metres]
            # Upper edges
            band_edgesm.append(s[1])  # [metres]

        # Block index 
        block_idx += 1
        
    return [v*1e9 for v in band_edgesm]  # convert to nm and return


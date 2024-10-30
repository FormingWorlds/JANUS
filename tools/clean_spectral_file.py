#!/usr/bin/env python3

# Script to remove NaN values from spectral file k-tables. These are sometimes
# produced in data for continuum absorption, and are replaced with zeros to
# avoid issues at runtime.

import os, shutil

kfile = 'Reach_k'

# Open original file
print("Cleaning k-table file")
old = open(kfile,'r')

# Work out the pitch of data in the file
lines = old.readlines()
sample = str(lines[10]).split() # use line 0 as a sample of the data
pitch = int(len(sample[0]) + 1)
if pitch < 8:
	raise Exception("Could not parse file '%s'" % kfile)

# Remove temp file if it already exists
temp_file = "temp_k_clean"
if os.path.isfile(temp_file):
    os.remove(temp_file)

# Replace NaN values with zero
new = open(temp_file,'w')
for l in lines:
	f = " "*(pitch-3) + "NaN" 				# search for this
	r = " 0." + "0"*(pitch-7) + "E+00" 		# replace with this
	n = str(l).replace(f,r)					# do replacement
	new.write(n)							# write to temp file
	
# Done writing
old.close()
new.close()

# Copy 'temp' file into location of 'old' file
os.remove(kfile)
shutil.copyfile(temp_file, kfile)
os.remove(temp_file)




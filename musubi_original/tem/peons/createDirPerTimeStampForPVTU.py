#!/usr/bin/env python
# Copyright (c) 2019 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this
# list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
# this list of conditions and the following disclaimer in the documentation
# and/or other materials provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
################################################################################
# AUTHOR: KANNAN MASILAMANI
################################################################################
# This python script creates pvd file from list of multiple pvtu files in 
# given directory
################################################################################
import sys
import os
import glob
import shutil
from tools.functions import *
################################################################################

in_files = None

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--files', help='Input *.pvtu files')
parser.add_argument('--out', help='Output folder')
args = parser.parse_args()

# First look for mandatory files as input.
if not args.files:
  print('''No files to work on provided.

           Please provide the files to process as arguments like:
           treelm/peons/harvest_series.py fileA fileB pre*.lua

           or in the configuration file as a list of comma separated files:
           files: fileA,fileB,pre*.lua

           Note that you can make use of globbing expansions here.

           As there are no files provided, there is nothing to do...
        ''')
  sys.exit(1)

print('Input files: {0}'.format(args.files))

# Now go on with the optional arguments, either provided in the config file,
# or the command line.
# The command line arguments take precedence over settings in the config file.
if args.out:
  out_folder = expanded_path(args.out)

print('Output folder: {0}'.format(out_folder))

# Check if the folder exists. If not -> create
if not os.path.exists(out_folder):
  if args.verbose:
    print('Output folder does not exist. Creating {0}'.format(out_folder))
  os.makedirs(out_folder)

if isinstance(args.files, str):
  infileList = getFileNames(args.files.split(','))
else:
  infileList = getFileNames(args.files)

# Sort files by timestamps, required for proper visit movie files.
# Only add each timestamp once.
filewithtime = []
seen = set()
for iFile in infileList:
  print('Found file {0}'.format(iFile))
  timestamp = float(iFile.split('_t')[1].strip('.pvtu'))
  dirName = 't'+str(timestamp)

  if timestamp not in seen:
    seen.add(timestamp)
    filewithtime.append( (timestamp, dirName, iFile) )

infileList = sorted(filewithtime)

funit = None

iFile_count = 0
# Process each file
for filetime in infileList:
  time  = filetime[0]
  dName = filetime[1]
  iFile = filetime[2]
  print('Processing file {0}'.format(iFile))

  pvtufile = os.path.basename(iFile.rstrip().split(' ')[-1])
  basepath = os.path.abspath(iFile.rstrip().split('/')[0])
  dPath = os.path.join(basepath, dName)
  if os.path.exists(dPath):
    try:
      os.rmdir(dPath)
    except OSError as error:
      print("Error: {0} : {1}".format(dPath, error.strerror))

  try:
    os.mkdir(dPath)
  except OSError as error:
    print("Error: {0} : {1}".format(dPath, error.strerror))

  timestamp = iFile.split('_t')[1].strip('.pvtu')
  # Move vtu and pvtu files matching current time stamp to dPath
  pattern = basepath + "/*_t" + timestamp + '.*'
  for tFile in glob.glob(pattern, recursive=True):
    # Extract file name form file path
    file_name = os.path.basename(tFile)
    shutil.move(tFile, dPath + '/' + file_name)

  if iFile_count == 0:
    pvdname = pvtufile.split('_t')
    pvdbase = ''.join(pvdname[0:-1]) + '.pvd'
    pvdfile = os.path.join(out_folder, pvdbase)
    if os.path.exists(pvdfile):
      try:
        os.remove(pvdfile)
      except OSError as e:
        print("Error: %s : %s" % (pvdfile, e.strerror))
    funit = open(pvdfile, 'w')
    writePVDHeader(funit)

  timestr = '{0:e}'.format(time)
  writePVDData(funit, dName+'/'+pvtufile, timestr)

  iFile_count = iFile_count + 1


# close pvd file
try:
  writePVDClose(funit)
except:
  print('Could not close pvd unit. Most likely something went wrong. Or no files were provided')


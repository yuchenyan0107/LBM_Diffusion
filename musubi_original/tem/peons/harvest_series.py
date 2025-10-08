#!/usr/bin/env python
# Copyright (c) 2016-2018 Harald Klimach <harald.klimach@uni-siegen.de>
# Copyright (c) 2016 Verena Krupp <verena.krupp@uni-siegen.de>
# Copyright (c) 2016-2017 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
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
# This python script runs harvester for time series and generates a *.pvd xml
# file to load multiple *.(p)vtu files in paraview.
# The configuration file (series.config) can be used to specify
# 1. time series of filenames as list or with linux wild cards
# 2. path to the harvesting executable
# 3. a template, in while filename and output folder is modified by this
#    script to generate the input harvester.lua for the harvesting tool.
#    These variables are named 'file' and 'out' in the template and are
#    replaced by the arguments 'files' and 'out' respectively.
#
# All options might be specified in the config file or as arguments on the
# command line.
# Except for the in_files, all settings have some sane default.
# The default configuration file is searched in series.config.
# You can provide an alternative configuration file name with the -c option.
# If no configuration file is found, only the command line arguments are used.
# Command line arguments take precedence over settings in the configuration
# file and settings in the configuration take precedence over the defaults.
################################################################################
import sys
################################################################################
# This script requires at least Python 2.7 to run.
assert sys.version_info >= (2,7), "Need at least Python 2.7 to run harvest_series!"
##----------  IMPORT LIBRARIES  ----------------------------------------------##
import re
import logging
import subprocess
import os
import shutil
import string
import fileinput

from tools.functions import *
from tools import configuration

##---- DONE Importing Libraries ----##

paramFile = None
in_files = None

args = configuration.get_args()

form = 'pvd'
stlBaseName = ''


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

# Executable to call for harvesting
print('Harvester executable: {0}'.format(expanded_path(args.harvester)))

# Command to run the executable (e.g. mpirun)
if args.run != '':
  print('Run command: {0}'.format(args.run))

# Visualize mesh only
if args.mesh: # template
  visMesh = True
else:
  visMesh = False

if visMesh:
  print('Visualizing MESHES only!')

# Template to use for harvester input
print('Template: {0}'.format(args.template))

# If template does not exist, create one
if not os.path.isfile(args.template) and not args.harvester == 'stl2vtp':
  print('Template file does not exist!')
  print('Creating a simple one in {0}'.format(args.template))
  tmp = open(args.template,'w')
  if visMesh: 
    print('Creating MESH visualization template {0}'.format(args.template))
    tmp.write("mesh = '$!file!$'\n")
  else:
    print('Creating data visualization template {0}'.format(args.template))
    tmp.write("restart = { read = '$!file!$' }\n")

  if int(args.samplelevel) > 0:
    tmp.write("ply_sampling = { nlevels = $!samplelevel!$ }\n")
  tmp.write("output_folder = '$!out!$' \n")
  tmp.write("output = { format = 'vtk', write_pvd=false }\n")
  tmp.close()
  created_template = True

else:
  created_template = False

# Get the format for the movie file
print('Movie format: {0}'.format(args.form))

### Finished option section ###


hvsfile = 'harvest_series.lua'
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
  if args.harvester == 'stl2vtp':
    timestamp = float(iFile.split('_t')[1].strip('.stl'))
  elif not visMesh:
    timestamp = float(getTimestamp(iFile,args.lua).strip())
  else:
    timestamp = 0.0

  if timestamp not in seen:
    seen.add(timestamp)
    filewithtime.append( (timestamp, iFile) )

infileList = sorted(filewithtime)

finished_pvdHeader = False
finished_trackHeader = False

# Create a sequence for the command to run in Popen
if args.run != '':
  runcmd = args.run.split()
else:
  runcmd = [] 
runcmd += [expanded_path(args.harvester), hvsfile]

regexp_track = re.compile('Writing TRACK to disk')
regexp_pvtu = re.compile('Opening P?VTU file')
funit = None
tl_unit = None

try:
  begintime = float(args.begin)
except:
  print('Your begin setting {0} needs to be a number!'.format(args.begin))
  raise

try:
  endtime = float(args.end)
except:
  print('Your end setting {0} needs to be a number!'.format(args.end))
  raise

try:
  delta_t = float(args.interval)
except:
  print('Your interval setting {0} needs to be a number!'.format(args.interval))
  raise

iFile_count = 0
nextinter = begintime
# Postprocess each file
for filetime in infileList:
  iFile = filetime[1]
  time  = filetime[0]

  if time >= begintime and time <= endtime:

    if time >= nextinter:
      print('Processing file {0}'.format(iFile))
      iFile_count = iFile_count + 1

      if not args.harvester == 'stl2vtp':
        # generate harvester.lua
        success = genscript(hvsfile, iFile, vars(args))
        if not success:
          print(' Error in generating intermediate script')
          sys.exit(1)

        # run harvester
        hvspipe = subprocess.Popen(runcmd, stdout=PIPE, stderr=PIPE)
        (hvslog, hvserr) = hvspipe.communicate()
        if sys.version_info[0] > 2:
          hvslog = hvslog.decode('ascii')
          hvserr = hvserr.decode('ascii')
        success = hvspipe.returncode
        if not args.quiet:
          print(hvslog)
        if success == 0:
          pvtufile = None
          trackfile = None
          for line in hvslog.split('\n'):
            if regexp_track.search(line):
              trackfile_withprefix = line.rstrip().split(' ')[-1]
              trackfile = os.path.basename(trackfile_withprefix)
            if regexp_pvtu.search(line):
              vtufile_withprefix = line.rstrip().split(' ')[-1]
              pvtufile = os.path.basename(vtufile_withprefix)
      else:
        pvtufile = None
        trackfile = None
        pvtufile = stl2vtp(iFile, out_folder, stlBaseName)
        if not pvtufile == None:
          success = 0
        else:
          success = 1
      
      if success == 0:
        if pvtufile:
          if form == 'visit':
            if os.path.extension(vtufile_withprefix) == '.pvtu':
              vtufiles = all_vtus(vtufile_withprefix)
            else:
              vtufiles = [vtufile_withprefix]
          if not finished_pvdHeader:
            pvdname = pvtufile.split('_t')
            pvdbase = ''.join(pvdname[0:-1]) + '.' + form
            pvdfile = os.path.join(out_folder, pvdbase)
            funit = open(pvdfile, 'w')
            if form == 'pvd':
              writePVDHeader(funit)
            elif form == 'visit':
              writeVISITHeader(funit, len(vtufiles))
            finished_pvdHeader = True

          if form == 'pvd':
            timestr = '{0:e}'.format(time)
            writePVDData(funit, pvtufile, timestr)
          elif form == 'visit':
            writeVISITData(funit, vtufiles, out_folder, iFile_count)

        if trackfile:
          if not finished_trackHeader:
            tl_file = os.path.join(out_folder, 'timeline.dat')
            tl_unit = open(tl_file, 'w')
            with open(trackfile_withprefix, 'r') as f:
              headline = f.readline()
            newhead = '# time' + headline[1:]
            tl_unit.write(newhead)
            finished_trackHeader = True

          # Get all the data from the tracking file an put it into the timeline
          with open(trackfile_withprefix, 'r') as f:
            next(f) # skip first line
            for line in f:
              tl_unit.write("{0} {1}".format(time, line))
      else:
        print('Harvester failed to run for restart file ' + iFile)
        if not args.quiet:
          print(hvserr)

      if delta_t > 0:
        nextinter = begintime + ((time - begintime)//delta_t + 1) * delta_t

    else:
      print('Skipping file {0}'.format(iFile))

  else:
    print('Skipping file {0}'.format(iFile))

# close pvd file
if finished_pvdHeader:
  if form == 'pvd':
    try:
      writePVDClose(funit)
    except:
      print('Could not close pvd unit. Most likely something went wrong. Or no files were provided')

if created_template:
  os.remove(args.template)

if finished_trackHeader:
  tl_unit.close()

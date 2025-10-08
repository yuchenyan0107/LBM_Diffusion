#!/usr/bin/env python
# Copyright (c) 2016 Harald Klimach <harald.klimach@uni-siegen.de>
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
## This file contains all functions used by harvest_series.py
from subprocess import Popen, PIPE
import os
import sys
import glob
import re

# Generate the lua script from the .template
# We use the pyratemplate to do this process
def genscript( hvsfile, infile, st_dict):
  from tools.pyratemp import pyratemp

  # Define Lua params as dict and use them in the dict
  #st_dict = {'out':out_folder+os.sep, 'samplelevel':samplelevel, 'file':infile}
  if st_dict['out'].endswith(os.sep):
    # If the out setting ends in the path separator, ensure the directory
    # exists.
    if not os.path.isdir(st_dict['out']):
      # The provided out path is not a directory, create it if it does not
      # already exist.
      if os.path.lexists(st_dict['out']):
        # The path already exists, but as it is not a directory we have a
        # conflict, raise this to the user...
        print('ERROR: Provided out prefix {0} exists already'.format(st_dict['out']))
        print('       and is not a path!')
        return False
      # Create not yet existing out directory.
      os.makedirs(st_dict['out'])

  luafile = open(hvsfile, "w")
  mytemplate = pyratemp.Template(filename=st_dict['template'], data=st_dict)
  luafile.write(mytemplate(file=infile))
  luafile.close()
  return True

def expanded_path(pathname):
  ''' Return the path with expanded environment variables and expanded user
      path. This is just a short cut for os.path.expandvars and
      os.path.expanduser applied both to the provided pathname.
  '''
  return(os.path.expandvars(os.path.expanduser(pathname)))

# find all the files match the string pattern given in infiles
def getFileNames(infiles):
  infileList = []
  for iFile in range(len(infiles)):
    infile_tmp = glob.glob(expanded_path(infiles[iFile]))
    for iFile2 in range(len(infile_tmp)):
       infileList.append(infile_tmp[iFile2])
  return infileList     

# get timestamp from restart file
def getTimestamp(infile, lua_exec):
  getTime = "dofile '"+infile+"'; print(string.format('%.7E', time_point.sim))"
  time = Popen([lua_exec,"-e", getTime], stdout=PIPE)
  stamp = time.communicate()[0]
  if sys.version_info[0] > 2:
    stamp = stamp.decode('ascii')
  return stamp

def writePVDHeader(funit):
  funit.write("<?xml version=\"1.0\"?>\n")
  funit.write("<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n")
  funit.write("<Collection>\n")

def writePVDData(funit, fname, time):
  funit.write(" <DataSet timestep=\"" + time + "\" file=\""+fname+"\"/>\n")

def writePVDClose(funit):
  funit.write("</Collection>\n")
  funit.write("</VTKFile>")

def writeVISITHeader(funit, nfiles):
  funit.write("!NBLOCKS {0}\n".format(nfiles))

def writeVISITData(funit, vtufiles, out_folder, iFile_count):
  for vtu in vtufiles:
    # old filename
    oldfile = os.path.join(out_folder, vtu)
    # new vtu filename with cycle count to generate movie in visit
    newvtu = vtu.strip('.vtu')+'_c'+str(iFile_count-1)+'.vtu'
    newfile = os.path.join(out_folder, newvtu)
    # rename oldfile with new filename
    os.rename(oldfile, newfile)
    funit.write(newvtu+'\n')

def stl2vtp(stlfile, out_folder, stlBaseName):
  import vtk
  # new vtk filename
  tmpvtpbase = stlfile.split('/')
  newvtpbase = tmpvtpbase[len(tmpvtpbase)-1]
  newvtp = newvtpbase.split('.stl')[0]+'.vtp'
  vtpfile = os.path.join(out_folder, newvtp)
  # read the stl file
  reader = vtk.vtkSTLReader()
  reader.SetFileName(stlfile)
  reader.Update()
  # catch empty stl basename
  if not stlBaseName == '':
    reader2 = vtk.vtkSTLReader()
    reader2.SetFileName(stlBaseName)
    reader2.Update()
    data1 = vtk.vtkPolyData()
    data2 = vtk.vtkPolyData()
    data1.ShallowCopy(reader.GetOutput())
    data2.ShallowCopy(reader2.GetOutput())
    # build polydata out of the two stl
    addPolyData = vtk.vtkAppendPolyData()
    addPolyData.AddInputConnection( data1.GetProducerPort() )
    addPolyData.AddInputConnection( data2.GetProducerPort() )
    addPolyData.Update()
  else:
    data1 = vtk.vtkPolyData()
    data1.ShallowCopy(reader.GetOutput())
    addPolyData = vtk.vtkAppendPolyData()
    addPolyData.AddInputConnection( data1.GetProducerPort() )
    addPolyData.Update()

  # write the vtk file
  w = vtk.vtkXMLPolyDataWriter()
  w.SetInputConnection(addPolyData.GetOutputPort())
  w.SetFileName(vtpfile)
  w.Write()
  return newvtp

def all_vtus(filename):
  piece = re.compile("<Piece Source=\"([^\"]*)\"/>")

  vtus = []

  with (open(filename)) as pf:
    for line in pf:
      pmatch = piece.search(line)
      if pmatch:
        vtus.append(pmatch.group(1))

  return vtus

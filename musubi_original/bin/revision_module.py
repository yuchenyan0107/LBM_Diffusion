def options(opt):
  """ Add an option to explicitly set the revision string for a project. """
  opt.add_option('--revision_string', action='store',
                 help='Set the revision string instead of trying to obtain it automatically.')


def configure(conf):
  """ The revision string may be set during configuration. """
  if conf.options.revision_string:
    conf.env.revision_string = conf.options.revision_string


def get_revision_of(projectdir):
  ''' Try to find a revision string of the code in projectdir.

      This routine will try to determine the git version
      in the provided project directory.
  '''
  import sys
  solver_rev = None
  try:
     import subprocess # If subprocess is available and provides check_output
     # This is working with Python > 2.7.
     # The check on check_output is necessary, because 2.6 provides subprocess
     # put without check_output.
     if getattr(subprocess, 'check_output'):
       use_subproc = True
     else:
       use_subproc = False
  except:
     import commands
     use_subproc = False

  if use_subproc:
    try:
      git_stat = 0
      git_out = subprocess.check_output(['git', 'describe', '--abbrev=12', '--always', '--dirty=+'],
                                        cwd=projectdir)
    except subprocess.CalledProcessError as e:
      git_stat = e.returncode
      git_out = e.output
    except OSError:
      git_stat = 1
  else:
    (git_stat, git_out) = commands.getstatusoutput('git describe --abbrev=12 --always --dirty=+',
                                                   cwd=projectdir)

  if git_stat == 0:
    if sys.version_info[0] > 2:
      solver_rev = git_out.split()[-1].decode('ascii')
    else:
      solver_rev = git_out.split()[-1]

  return solver_rev

def get_tag_of(projectdir):
  ''' Try to find a version string of the code in projectdir.

      This routine will try to determine the git tag
      in the provided project directory.
  '''
  import sys
  solver_tag = None
  try:
     import subprocess # If subprocess is available and provides check_output
     # This is working with Python > 2.7.
     # The check on check_output is necessary, because 2.6 provides subprocess
     # put without check_output.
     if getattr(subprocess, 'check_output'):
       use_subproc = True
     else:
       use_subproc = False
  except:
     import commands
     use_subproc = False

  if use_subproc:
    try:
      git_stat = 0
      git_out = subprocess.check_output(['git', 'describe', '--tags'], cwd=projectdir)
    except subprocess.CalledProcessError as e:
      git_stat = e.returncode
      git_out = e.output
    except OSError:
      git_stat = 1
  else:
    (git_stat, git_out) = commands.getstatusoutput('git describe --tags', cwd=projectdir)

  if git_stat == 0:
    if sys.version_info[0] > 2:
      solver_tag = git_out.split()[-1].decode('ascii')
    else:
      solver_tag = git_out.split()[-1]

  if not solver_tag:
    solver_tag = "0.0.0"

  return solver_tag

def fill_revision_string(bld, subdir='.'):
  """ Update the bld.env.revision_string information.

      If the bld.env.revision_string is not yet known, try to
      construct it from mercurial information in the provided
      subdirectory.
  """
  import os
  if not bld.env.revision_string:
    construct_revision = True
  elif bld.env.revision_string == 'UNKNOWN':
    construct_revision = True
  else:
    construct_revision = False

  bld.env.version_tag = get_tag_of(bld.path.abspath())
  if construct_revision:
    projectdir = os.path.join(bld.path.abspath(), subdir)
    bld.env.revision_string = get_revision_of(projectdir)


def revision_module_text(env):
  """ Define a revision module describing the project revision and the
      compilation information.

      The revision information is taken from env.revision_string,
      which might be set on the command line by --revision_string.
      To have the revision string properly defined before creating this
      module file, fill_revision_string should be called.
  """

  import datetime
  from waflib import Logs

  fc_name_str = "".join(env['FC_NAME'])
  fc_version_str = ".".join(env['FC_VERSION'])
  fc_flags_str = " ".join(env['FCFLAGS'])
  link_flags_str = " ".join(env['LINKFLAGS'])
  builddate = datetime.datetime.now().strftime("%Y-%m-%d")

  Logs.info('Compiler name   : {0}'.format(fc_name_str))
  Logs.info('Compiler version: {0}'.format(fc_version_str))
  Logs.info('Compiler options: {0}'.format(fc_flags_str))
  Logs.info('Project version : {0}'.format(env.version_tag))
  Logs.info('Project revision: {0}'.format(env.revision_string))


  flaglen = len(fc_flags_str)
  nFlagLines = max((flaglen+71) // 72, 1)

  modtext = """!> SOIL module for holding the revision and compilation information
!! of the executable.
!!
!! This information will be written by soi_world_init.
!! This source file is generated during compilation!
! *************************************************************************** !
! WARNING: Do NOT change this file, as it will be overwritten during
!          compilation.
!          See bin/revision_module.py for the generating script.
!          (in the apes build infrastructure project)
! *************************************************************************** !

module soi_revision_module
  implicit none
  !> The git revision of the application used for this executable.
  character(len=13), parameter :: soi_solver_revision &
    &                            = '{0}'

  !> The version of this application
  character(len=16), parameter :: soi_solver_version &
    &                            = '{5}'

  !> Name of the compiler.
  character(len=32), parameter :: soi_FC_name &
    &                            = '{1}'

  !> The compilation command that was used to build this executable.
  character(len=32), parameter :: soi_FC_command &
    &                            = '{2}'

  !> The version of the Fortran compiler used in the compilation of this
  !! executable.
  character(len=32), parameter :: soi_FC_version &
    &                            = '{3}'

  !> Number of lines needed to represent the compiler flags
  integer, parameter :: soi_FC_nFlagLines = {4}

  !> The Fortran compiler flags used to compile this executable.
  character(len=72), parameter :: soi_FC_flags(soi_FC_nFlagLines) &
""".format(env.revision_string[:13],
           fc_name_str[:32], " ".join(env.FC)[:32],
           fc_version_str[:32], nFlagLines, env.version_tag[:16])
  tempstr = fc_flags_str[0:72]
  tempstr = tempstr + ' '*(72-len(tempstr))
  modtext = modtext + """    & = [ '%s'""" % (tempstr[0:72])
  for line in range(2,nFlagLines+1):
     tempstr = fc_flags_str[(line-1)*72:line*72]
     tempstr = tempstr + ' '*(72-len(tempstr))
     modtext = modtext + """, &
    &     '{0}'""".format(tempstr)
  modtext = modtext + """ ]

  !> The date when this executable was built.
  character(len=10), parameter :: soi_build_date &
    &                            = '{0}'
end module soi_revision_module
""".format(builddate)

  return(modtext)


def revision_module_file(task):
  return(task.outputs[0].write(revision_module_text(task.env)))


from waflib import TaskGen

@TaskGen.feature('revmod')
@TaskGen.before('process_source')
def create_revmod(self):
  """ Define a revision module describing the project revision and the
      compilation information.

      The revision information is taken from self.env.revision_string,
      which might be set on the command line by --revision_string.
      To have the revision string properly defined before creating this
      module file, fill_revision_string should be called.
  """
  node = self.path.get_bld()
  node = node.make_node('soi_revision_module.f90')
  node.parent.mkdir()
  node.write(revision_module_text(self.env))
  if "fc" in self.features:
    self.source.append(node)

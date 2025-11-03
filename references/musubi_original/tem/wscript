#! /usr/bin/env python
# encoding: utf-8
# Harald Klimach 2011

APPNAME = 'tem-source'
VERSION = '1'

import sys
import os

top = '.'
out = 'build'

def options(opt):
    '''Building options provided by the treelm library.
       Remember, all options can be displayed with waf --help.'''

    opt.load('coco')

    # General options
    opt.add_option('--nowarn', action='store_true', default=False,
                   help='Do not show compiler warnings.')
    opt.add_option('--openmp', action='store_true', default=False,
                   help='Build with OpenMP support')
    opt.add_option('--mic', action='store_true', default=False,
                   help='Build for Intels MIC architecture')
    opt.add_option('--vlen', action='store',
                   help='Vector length to use')
    opt.add_option('--revision_string', action='store',
                   help='Set the revision string instead of trying to obtain it autmatically.')
    opt.add_option('--print-commands', action='store_true', default=False,
                   help='Pretty print the commands executed during build.', dest='print_cmds')
    opt.add_option('--nag_path', action='store', help='Path to NAG Fortran library')
    opt.add_option('--with-distcrc', action='store_true', default=False,
                   help='Enable usage of distributed checksums. Useful for debugging.',
                   dest='with_distcrc')
    opt.add_option('--precice_fbindpath', action='store',
                   help='Path to Fortran 2003 binding sources of PreCICE. '
                        'Can also be set via the PRECICE_FBIND environment variable. '
                        'REQUIRED if you want to enable PreCICE support. '
                        'If this option is not given, waf will not attempt '
                        'to find the PreCICE library.')
    opt.add_option('--precice_libpath', action='store', help='Path to PreCICE library')

    #Option to use Xevolver
    opt.add_option('--xevolver', action='store_true', default=False,
                   help='Preprocessing using Xevolver')

    # Get paths to libraries required by Atles
    # (if not located in default directories)
    atl_opt = opt.add_option_group('Ateles specific options')
    atl_opt.add_option('--fftw_path', action='store', help='Path to the FFTW')
    atl_opt.add_option('--fftw_libpath', action='store', help='Path to the FFTW library')
    atl_opt.add_option('--fftw_incpath', action='store', help='Path to the FFTW includes')
    atl_opt.add_option('--nfft_path', action='store', help='Path to NFFT library')

    # Get paths to libraries required by Seeder
    # (if not located in default directories)
    sdr_opt = opt.add_option_group('Seeder specific options')
    sdr_opt.add_option('--libNBC_path', action='store', help='Path to NBC library')

    # De-activate vtk output.
    opt.add_option('--no_vtk', action = 'store_true',
                   default = False,
                   help = 'Do not include vtk output files for compilation',
                   dest = 'no_vtk')


def configure(conf):
    '''Common configuration tasks for all projects using treelm.'''
    import os
    from waflib import Logs

    conf.load('revision_module')
    conf.load('coco')

    conf.setenv('')
    conf.env.ford_mainpage = 'tem_mainpage.md'
    conf.env.fordurl_tem = 'https://apes-suite.github.io/treelm/'
    conf.env.distcrc = ''

    conf.setenv('cenv')
    if conf.options.with_distcrc:
      # Checking availability of ZLIB only if MPI enabled C-Compiler available.
      conf.check_cc(header_name='mpi.h', mandatory=False)
      if conf.is_defined('HAVE_MPI_H'):
        conf.check(lib='z', uselib_store='ZLIB', mandatory=False)
        if conf.env.LIB_ZLIB:
          conf.all_envs[''].distcrc = 'distcrc.o'

    conf.setenv('')

    ### MPI ###
    conf.check_fc(fragment = '''
program check_mpi
  use mpi
end program check_mpi''',
                  msg = 'Checking if MPI provides a module',
                  errmsg = 'NO! Attempting to create a custom mpi module.',
                  mandatory=False, define_name='has_mpi')
    # If the simple mpi program failed, let's see if a manual inclusion
    # of mpif.h does the trick. If so, then we can proceed and later on
    # include the user defined mpi module, which simply contains include mpif.h
    # If not, this is fatal and we have to abort compilation
    if ( not conf.is_defined('has_mpi')):
      conf.check_fc(fragment = '''
program check_mpi
  implicit none
  include "mpif.h"
  write(*,*) MPI_COMM_WORLD
end program check_mpi''',
               msg = 'Checking for mpif.h header',
               okmsg = 'MPI environment successful',
               errmsg = '''ERROR!

Check if your FC variable points to the MPI wrapper.
(e.g.: export FC=mpif90)
               ''',
               mandatory=True)

    if conf.check_fc( features  = 'fc fcprogram',
                      msg       = 'Checking for non-blocking collectives',
                      errmsg    = ''' Non-blocking collectives NOT supported.

    The sparse alltoall implementation will not be available!
                               ''',
                      mandatory = False,
                      fragment  = '''
program check_nbc
  use mpi

  implicit none

  integer  :: req,iError
  call mpi_ibarrier(MPI_COMM_WORLD, req, iError)
end program check_nbc
'''):
      conf.env.sparse_a2a_sources = ['nbc_mpi/tem_sparse_comm_module.f90']
    else:
      conf.env.sparse_a2a_sources = ['nbc_mpi/tem_sparse_comm_dummy.f90']

    conf.env.mpi_mem_c_sources = []
    if conf.check( features    = 'fc',
                   fragment    = conf.path.find_node('mem_for_mpi/mem_for_mpi_f08_module.f90').read(),
                   msg         = 'Checking for mpi_f08 MPI_Alloc_mem',
                   mandatory   = False ):
        conf.env.mpi_mem_f_sources = ['mem_for_mpi/mem_for_mpi_f08_module.f90']

    else:
      if conf.check( features    = 'fc',
                     fragment    = conf.path.find_node('mem_for_mpi/mem_for_mpi_f_module.f90').read(),
                     msg         = '> Checking for (modern) Fortran MPI_mem Wrapper',
                     mandatory   = False ):
        conf.env.mpi_mem_f_sources = ['mem_for_mpi/mem_for_mpi_f_module.f90']

      else:
        conf.env.mpi_mem_f_sources = ['mem_for_mpi/mem_for_mpi_c_module.f90']
        if conf.check( features    = 'c',
                       fragment    = conf.path.find_node('mem_for_mpi/mem_for_mpi.c').read(),
                       msg         = '>> Checking for C MPI_mem Wrapper',
                       errmsg      = 'C compiler does not support MPI, using dummy implementation.',
                       mandatory   = False ):
          conf.env.mpi_mem_c_sources = ['mem_for_mpi/mem_for_mpi.c']
        else:
          conf.env.mpi_mem_c_sources = ['mem_for_mpi/mem_for_dummy.c']
    ### End MPI section ###

    tmpDEF = conf.env.DEFINES

    ### Check for NAG library ###
    if conf.options.nag_path:
       Logs.info('Checking path to NAG library: ' + conf.options.nag_path + '/lib')
       conf.check_fc(lib='nag_nag', libpath=(conf.options.nag_path+'/lib'),
                  uselib_store='NAG', mandatory=False)
       if conf.env.LIB_NAG:
          conf.env.INCLUDES_NAG = [conf.options.nag_path+'/nag_interface_blocks']
    if conf.env.LIB_NAG:
       try:
          conf.check_fc(fragment = '''
program check_NAG_header
  use nag_library
  implicit none
  real(kind=nag_wp) :: a_real
  real(kind=nag_wp) :: res
  read(*,*) a_real
end program check_NAG_header''',
                      msg = 'Checking for NAG library header',
                      errmsg = 'NAG headers not found, disabling the NAG library',
                      use = 'NAG')
       except Errors.ConfigurationError:
          conf.env.LIB_NAG=None
          conf.env.INCLUDES_NAG=None

    ### END NAG library ###


    ### PreCICE ###
    conf.setenv('cenv')
    conf.env.PRECICE_FBINDPATH = os.environ.get('PRECICE_FBINDPATH')
    if conf.options.precice_fbindpath:
      conf.env.PRECICE_FBINDPATH = conf.options.precice_fbindpath


    if conf.env.PRECICE_FBINDPATH:
      # Activate PreCICE only, if the fbindpath option is given.
      Logs.info('Adding PreCICE Fortran interface from: '
                + conf.env.PRECICE_FBINDPATH)
      if os.path.exists(conf.env.PRECICE_FBINDPATH):
        if conf.options.precice_libpath:
           Logs.info('Checking path to PreCICE library: ' + conf.options.precice_libpath)
           conf.check(lib='precice', libpath=conf.options.precice_libpath,
                      uselib_store='PRECICE', mandatory=True)
        else:
           conf.check(lib='precice', uselib_store='PRECICE', mandatory=True)
      else:
        # The fbindpath is given but not found
        Logs.warn('Given directory for PreCICE Fortran bindings does not exist: ',
                  conf.env.PRECICE_FBINDPATH)
        Logs.warn('!! PreCICE support DEACTIVATED !!')

    conf.all_envs[''].LIB_PRECICE = conf.env.LIB_PRECICE
    conf.all_envs[''].LIBPATH_PRECICE = conf.env.LIBPATH_PRECICE
    conf.all_envs[''].PRECICE_FBINDPATH = conf.env.PRECICE_FBINDPATH
    conf.setenv('')
    if conf.env.LIB_PRECICE:
      # The directory for the fortran bindings exists, and libprecice was found.
      # Should be fine to proceed with PreCICE support.
      conf.check_fc(lib='stdc++', uselib_store='STDCXX', mandatory=True)

      # mpi_cxx needed when using openmpi:
      conf.check_fc(lib='mpi_cxx', uselib_store='MPICXX', mandatory=False)
      conf.check_fc(lib='rt', uselib_store='RT', mandatory=False)

      # libprecice might make use of libpython, check for 2.7 or, if that
      # fails 2.6:
      conf.check_fc(lib='python2.7', uselib_store='PYLIB', mandatory=False)
      if not conf.env.LIB_PYLIB:
        conf.check_fc(lib='python2.6', uselib_store='PYLIB', mandatory=False)
      if not conf.env.LIB_PYLIB:
        Logs.warn('No Python library found to support libprecice!')

      prcf_parent = conf.path.abspath()+'/external/precice'
      precicef = prcf_parent + '/conf_fbind'
      if os.path.exists(precicef):
        #print (os.path.islink(precicef))
        if os.path.islink(precicef):
          os.unlink(precicef)
        else:
          if os.path.isdir(precicef):
            os.removedir(precicef)
          else:
            os.remove(precicef)
      else:
        if os.path.islink(precicef):
          os.unlink(precicef)

      path_to_prec = os.path.relpath(os.path.abspath(conf.env.PRECICE_FBINDPATH),
                                     prcf_parent)
      os.symlink(path_to_prec, precicef)
      conf.all_envs[''].precicef = os.path.relpath(precicef, conf.path.abspath())
    ### End PreCICE section ###

    conf.env.DEFINES = tmpDEF

    if conf.options.vlen:
      conf.env.vlen = conf.options.vlen
    else:
      conf.env.vlen = 256

    if conf.options.revision_string:
      conf.env.revision_string = conf.options.revision_string

    ### Check for Fortran features ###
    check_fortran_features(conf)


    conf.env['no_vtk'] = conf.options.no_vtk

    # Deactivate VTK output, if no sizeof available
    if not 'sizeof_source' in conf.env:
       Logs.warn('No sizeof function found, deactivating VTK output!')
       conf.env['no_vtk'] = True

    set_variant_flags(conf)


def check_fortran_features(conf):

    import fortran_language

    # Looking for is_NaN support:
    fortran_language.supports_ieee_is_NaN(conf = conf, mandatory = False)
    if conf.env['fortsupp_ieee_is_NaN']:
      conf.env['isNaN_source'] = 'source/Fortran_Features/isNaN/tem_isNaN_IEEE.f90'
    else:
      fortran_language.supports_vendor_is_NaN(conf = conf, mandatory = False)
      if conf.env['fortsupp_vendor_is_NaN'].lower() == 'isnan':
        conf.env['isNaN_source'] = 'source/Fortran_Features/isNaN/tem_isNaN_vendor.f90'
      else:
        conf.env['isNaN_source'] = 'source/Fortran_Features/isNaN/tem_isNaN_dummy.f90'

    # Do not change the DEFINES themselves, they are not needed for Fortran.
    # Copying them here to restore them later on.
    tmpDEF = conf.env.DEFINES

    ################# F2008 #######################

    # Check for Gamma function:
    fortran_language.supports_f2008_gamma(conf = conf, mandatory = False)
    if conf.env['fortsupp_f2008_gamma']:
        conf.env['gamma_source'] = 'f2008'
    else:
        # No compiler support, look for NAG gamma function #
        conf.check_fc(fragment = '''
program check_gamma
  use nag_library, only: s14aaf, nag_wp
  implicit none
  real(kind=nag_wp) :: a_real
  real(kind=nag_wp) :: res
  integer :: iFail
  read(*,*) a_real
  res = s14aaf(a_real,iFail)
end program check_gamma''',
                      msg = 'Checking for NAG specific Gamma function',
                      mandatory=False, use='NAG', define_name='has_NAG_gamma')
        if (conf.is_defined('has_NAG_gamma')):
            conf.env['gamma_source'] = 'nag'
        else:
            # Neither compiler support nor NAG implementation, use slatec #
            conf.env['gamma_source'] = 'slatec'

    # Check for Bessel function:
    fortran_language.supports_f2008_bessel(conf = conf, mandatory = False)


    # Check for sizeof function
    if conf.env.FC_NAME == 'NEC':
       conf.env['sizeof_source'] = ['wrap_sizeof/hvs_nec_sx_sizeof.f90']
    else:
      fortran_language.supports_c_sizeof(conf = conf, mandatory = False)
      if conf.env['fortsupp_c_sizeof']:
        conf.env['sizeof_source'] = ['wrap_sizeof/hvs_iso_c_sizeof.f90']
      else:
        fortran_language.supports_vendor_sizeof(conf = conf, mandatory = False)
        if conf.env['fortsupp_vendor_sizeof']:
          conf.env['sizeof_source'] = ['wrap_sizeof/hvs_vendor_sizeof.f90']


    # Restore the DEFINES again
    conf.env.DEFINES = tmpDEF


def set_variant_flags(conf):
    '''Setting the compilation flags for the different variants.
       This is separated from the configuration step above, to
       allow for more configuration settings in the individual
       project wscripts, using the treelm library, that are to
       affect all build variants.'''

    from waflib import Logs

    from fortran_compiler import set_fc_flags, fcopts
    # includes options for:
    # * 'warn': activate compile time warnings
    # * 'w2e': turn warnings into errors
    # * 'standard': check for standard compliance
    # * 'debug': activate debugging facilities
    # * 'optimize': turn optimization on
    # * 'profile': activate profiling facilities
    # * 'double': promote default reals to double precision
    # * 'openmp': enable OpenMP

    conf.env.fordurl_tem = 'https://apes-suite.github.io/treelm/'

    omp_flags = fcopts[conf.env.FC_NAME, 'noomp']
    conf.env.flag_fixform = fcopts[conf.env.FC_NAME, 'fixform']
    if conf.options.openmp:
       omp_flags = fcopts[conf.env.FC_NAME, 'openmp']

    mic_flags = []
    if conf.options.mic:
      mic_flags = ['-mmic']
      conf.env.append_value('CFLAGS', '-mmic')
      if conf.options.fc_delflags:
        conf.options.fc_delflags = delflags + ['-xHOST']
      else:
        conf.options.fc_delflags = ['-xHOST']

    mismatch_flags = []
    if conf.env.FC_NAME == 'NAG':
       mismatch_flags = ['-wmismatch=mpi_bcast,mpi_send,mpi_isend,mpi_irecv,mpi_allreduce,mpi_reduce']

    osfcflags = conf.env.FCFLAGS + omp_flags + mic_flags + mismatch_flags


    # Default Fortran compiler
    FC_def = conf.env.FC[0]

    # Default C compiler
    CC_def = conf.all_envs['cenv'].CC[0]

    if conf.options.nowarn:
      warn = []
    else:
      warn = ['warn', 'w2e']

    conf.env.fcwarn = []
    for opt in warn:
      conf.env.fcwarn += fcopts[conf.env.FC_NAME, opt]

    # Flags for the default (production) variant
    set_fc_flags(conf, ['optimize'] + warn, osfcflags)
    Logs.warn('Default flags: ' + ' '.join(conf.env['FCFLAGS']))

    # Set flags for the debug variant
    conf.setenv('debug', conf.env)
    set_fc_flags(conf, ['debug'] + warn, osfcflags)
    Logs.warn('Debug flags: ' + ' '.join(conf.env['FCFLAGS']))
    conf.env.append_value('COCOFLAGS', ['-DDEBUG'])

    conf.setenv('')
    conf.setenv('profile', conf.env)
    set_fc_flags(conf, ['optimize', 'profile'] + warn,
                 osfcflags)

    #todo#  # Set flags for VAMPIR instrumented compilation
    #todo#  conf.setenv('vampir', conf.env)
    #todo#  myflags = ['-vt:f90',FC_def,'-g'] + fcopts[conf.env.FC_NAME, 'optimize'] + addflags
    #todo#  conf.env['FCFLAGS'] = [flag for flag in myflags if not flag in delflags]
    #todo#  conf.env['CCFLAGS'] = ['-vt:cc',CC_def,'-g']
    #todo#  conf.env['FC'] = 'vtf90'
    #todo#  conf.env['CC'] = 'vtcc'
    #todo#  conf.env['LINKFLAGS_fcprogram'] = conf.env['FCFLAGS']

    # Set flags for SCALASCA instrumented compilation
    conf.setenv('')
    conf.setenv('scalasca', conf.env)
    set_fc_flags(conf, ['optimize'], [FC_def, '-g'] + osfcflags)
    conf.env['FC'] = 'skin'

    # Set flags for SCOREP instrumented compilation
    conf.setenv('')
    conf.setenv('scorep', conf.env)
    if conf.options.openmp:
      set_fc_flags(conf, ['optimize'], ['--openmp', '--user', FC_def]
                                       + osfcflags)
    else:
      set_fc_flags(conf, ['optimize'], ['--user', FC_def] + osfcflags)

    conf.env['FC'] = 'scorep'
    conf.env.append_value('COCOFLAGS', ['-DSCOREP_USER_ENABLE'])
    conf.setenv('cenv')
    conf.setenv('cenv_scorep',conf.env)

    #todo#  # Set flags for TAU instrumented compilation
    #todo#  conf.setenv('tau', conf.env)
    #todo#  myflags = ['-g'] + fcopts[conf.env.FC_NAME, 'optimize'] + addflags
    #todo#  conf.env['FCFLAGS'] = [flag for flag in myflags if not flag in delflags]
    #todo#  conf.env['FC'] = 'tau_f90.sh'
    #todo#  conf.env['CC'] = 'tau_cc.sh'
    #todo#  conf.env['LINKFLAGS_fcprogram'] = conf.env['FCFLAGS']


def build(bld):
    import revision_module
    from waflib.extras.utest_results import utests
    from waflib import Logs
    import sys

    bld.load('coco')

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
    import datetime
    from waflib import Utils
    if bld.options.print_cmds:
        import waflib.extras.print_commands

    ### Setup the soi_revision_module ###
    bld( features = "revmod fc",
         target = "rev_module" )
    ### End setup soi_revision_module ###

    ccmod = bld.path.find_or_declare('source/tem_compileconf_module.f90')
    modtext = """!> TreElm module for holding compile time configurations.
!*****************************************************************************!
! WARNING: Do NOT change this file, as it will be overwritten during
!          compilation.
!*****************************************************************************!

module tem_compileconf_module
  implicit none

  integer, parameter :: vlen = %s
end module tem_compileconf_module
""" % (bld.env.vlen)

    ccmod.write(modtext)
    ccmod.sig = Utils.h_file(ccmod.abspath())

    tem_ppsources = bld.path.ant_glob('source/*.fpp')
    tem_ppsources += bld.path.ant_glob('source/variables/*.fpp')
    tem_ppsources += bld.path.ant_glob('source/shapes/*.fpp')

    lapack_sources = bld.path.ant_glob('external/lapack/*.f')
    blas_sources = bld.path.ant_glob('external/blas/*.f')
    fnlib_sources = bld.path.ant_glob('external/fnlib/*.f')

    tem_sources = bld.path.ant_glob('source/*.f90')
    tem_sources += bld.path.ant_glob('source/control/*.f90')
    tem_sources += bld.path.ant_glob('source/faces/*.f90')
    tem_sources += bld.path.ant_glob('source/variables/*.f90')
    tem_sources += bld.path.ant_glob('source/shapes/*.f90')
    tem_sources += bld.path.ant_glob('source/libharvesting/*.f90',
                                     excl=['source/libharvesting/hvs_vtk_dummy.f90',
                                     'source/libharvesting/hvs_vtk_module.f90'])
    tem_sources += bld.env.sizeof_source
    tem_sources.append('source/space_time_functions/tem_polygon_material_module.f90')
    tem_sources += bld.path.ant_glob('external/stla_io.f90')

    # Bessel dependend features:
    if bld.env['fortsupp_f2008_bessel']:
      tem_sources.append('source/space_time_functions/tem_miescatter_module.f90')
      tem_sources.append('source/space_time_functions/tem_cylindricalWave_module.f90')
    else:
      tem_sources.append('source/space_time_functions/dummy/tem_miescatter_module.f90')
      tem_sources.append('source/space_time_functions/dummy/tem_cylindricalWave_module.f90')

    tem_sources.append('source/space_time_functions/tem_acoustic_pulse_module.f90')

    tem_sources += bld.env.mpi_mem_f_sources
    tem_sources += bld.env.sparse_a2a_sources

    tem_sources.append(ccmod)
    tem_sources.append(bld.env['isNaN_source'])
    if bld.env['gamma_source'] == 'f2008':
      tem_sources.append('source/Fortran_Features/gamma/tem_gamma_F2008.f90')
    elif bld.env['gamma_source'] == 'nag':
      tem_sources.append('source/Fortran_Features/gamma/tem_gamma_NAG.f90')
    else:
      tem_sources.append('source/Fortran_Features/gamma/tem_gamma_SLATEC.f90')
      tem_sources += bld.path.ant_glob('external/fnlib/*.f')

    if ( not bld.is_defined('has_mpi')):
        tem_sources.append('source/mpi_header/mpi_module.f90')

    if bld.env.LIB_PRECICE:
      tem_sources += bld.path.ant_glob(bld.env.precicef + '/*.f03')
      tem_sources += bld.path.ant_glob(bld.env.precicef + '/*.f90')
      tem_sources += bld.path.ant_glob('external/precice/*.f90')
    else:
      # No precice library available, use the dummy implementation instead.
      tem_sources += ['external/dummy/tem_precice_module.f90']
      tem_sources += ['external/dummy/tem_precice_interpolation_module.f90']

    if bld.all_envs['cenv'].LIB_ZLIB:
      tem_sources += ['external/dist_check/distcheck_module.f90']

    # If not vtk, include dummy implementation
    if bld.options.no_vtk or bld.env.no_vtk:
      print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
      print('        NO VTK SUPPORT')
      print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
      tem_sources += bld.path.ant_glob('source/libharvesting/hvs_vtk_dummy.f90')
    else:
      tem_sources += bld.path.ant_glob('source/libharvesting/hvs_vtk_module.f90')
      tem_ppsources.append(bld.path.find_node('source/libharvesting/hvs_base64_module.fpp'))

    tem_sources += tem_ppsources


    if bld.cmd != 'docu':
      bld.env.mpi_mem_c_obj = 'mem_for_mpi.o'
      if bld.env.mpi_mem_c_sources != '':
        bld.env.mpi_mem_c_obj = ''
        bld( features = 'c',
             source = bld.env.mpi_mem_c_sources,
             target = bld.env.mpi_mem_c_obj )

      if bld.all_envs['cenv'].LIB_ZLIB:
        bld( features = 'c',
             source = 'external/dist_check/distcrc.c',
             target = bld.env.distcrc )

      if not bld.options.no_vtk:
        bld(
            features = 'c',
            source   = 'external/base64/Base64EncodeDecode.c',
            target   = 'base64')

      bld(
        features = 'fc',
        fcflags  = bld.env.flag_fixform,
        delflags = bld.env.fcwarn,
        source   = blas_sources,
        target   = 'blas_objs')

      bld(
        features = 'fc',
        fcflags  = bld.env.flag_fixform,
        delflags = bld.env.fcwarn,
        source   = lapack_sources,
        target   = 'lapack_objs')

      bld(
        features = 'coco fc',
        source   = tem_sources,
        use      = ['blas_objs', 'rev_module', 'lapack_objs', 'aotus', 'NAG', 'base64'],
        target   = 'tem_objs')

      utest_sources = bld.path.ant_glob('utests/*_module.f90')
      lib_deps = ['lapack_objs', 'rev_module', 'aotus', 'tem_objs',
                  bld.env.mpi_mem_c_obj,
                  bld.env.distcrc,
                  'PRECICE', 'MPICXX', 'PYLIB', 'STDCXX', 'ZLIB']
      bld(
        features = 'fc',
        source   = utest_sources,
        target   = 'tem_utest_objs')
      test_deps = lib_deps + ['tem_utest_objs']

      bld(
        features = 'fc fcprogram',
        source   = 'peons/treelmesh_connectivity.f90',
        use      = lib_deps,
        target   = 'treelmesh_connectivity')

      bld(
        features = 'fc fcprogram',
        source   = 'peons/iar_bench.f90',
        use      = lib_deps,
        target   = 'iar_bench')

      utests(bld, test_deps, preprocessor='coco')

    else:
      import os
      from waflib.extras.make_fordoc import gendoc

      tpp = bld(
        features = 'includes coco',
        source   = tem_ppsources)

      tpp.post()
      tem_preprocessed = []
      for ppt in tpp.tasks:
        for f in ppt.outputs:
          tem_preprocessed.append(f)

      if not bld.env.fordonline:
        tem_preprocessed.append(bld.env.fordext_aotus)

      tgt = bld.path.get_bld().make_node('docu/modules.json')
      bld.env.fordext_tem = tgt

      bld( rule = gendoc,
           source = tem_preprocessed,
           src_paths = [bld.path.find_node('source').abspath()],
           target = tgt,
           extern = ['aoturl = {0}'.format(bld.env.fordext_aotus)],
           extern_urls = ['aoturl = {0}'.format(bld.env.fordurl_aotus)],
           mainpage = os.path.join(bld.top_dir, 'tem', 'tem_mainpage.md')
      )


from waflib.Build import BuildContext
class debug(BuildContext):
    "Build a debug executable"
    cmd = 'debug'
    variant = 'debug'

class test(BuildContext):
    "Unit Tests"
    cmd = 'test'
    variant = 'debug'

class profile(BuildContext):
    "Build a profile executable"
    cmd = 'profile'
    variant = 'profile'

class scalasca(BuildContext):
    "Build a scalasca instrumented executable"
    cmd = 'scalasca'
    variant = 'scalasca'

class scorep(BuildContext):
    "Build a scorep instrumented executable"
    cmd = 'scorep'
    variant = 'scorep'
#todo#
#todo# class vampir(BuildContext):
#todo#     "Build a vampir instrumented executable"
#todo#     cmd = 'vampir'
#todo#     variant = 'vampir'
#todo#
#todo# class tau(BuildContext):
#todo#     "Build a tau instrumented executable"
#todo#     cmd = 'tau'
#todo#     variant = 'tau'

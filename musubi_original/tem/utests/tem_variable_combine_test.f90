! Copyright (c) 2015-2016 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2016, 2019 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2016-2017 Peter Vitt <peter.vitt2@uni-siegen.de>
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are met:
!
! 1. Redistributions of source code must retain the above copyright notice, this
! list of conditions and the following disclaimer.
!
! 2. Redistributions in binary form must reproduce the above copyright notice,
! this list of conditions and the following disclaimer in the documentation
! and/or other materials provided with the distribution.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
! DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
! FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
! DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
! SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
! OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!> This unit test tests the different evaluation types for user defined
!! variables.
program tem_variable_combine_Test

  use, intrinsic :: iso_c_binding, only: C_NEW_LINE, c_ptr, c_loc, c_f_pointer

  use env_module,               only: eps, rk, solSpecLen, fin_env

  use aotus_module,             only: open_config_chunk, &
    &                                 close_config,      &
    &                                 flu_state
  use aot_table_module,         only: aot_table_open,   &
    &                                 aot_table_close,  &
    &                                 aot_table_length, &
    &                                 aot_get_val

  use tem_logging_module,       only: tem_logging_init
  use treelmesh_module,         only: treelmesh_type
  use tem_bc_prop_module,       only: tem_bc_prop_type
  use tem_general_module,       only: tem_general_type
  use tem_variable_module,      only: tem_variable_type, &
    &                                 tem_variable_load
  use tem_spacetime_fun_module, only: tem_st_fun_linkedList_type
  use tem_derived_module,       only: tem_varSys_append_luaVar
  use tem_spacetime_fun_module, only: tem_create_subTree_of_st_funList
  use tem_varSys_module,        only: tem_varSys_init,            &
    &                                 tem_varSys_type,            &
    &                                 tem_varSys_op_type,         &
    &                                 tem_varSys_append_stateVar, &
    &                                 tem_varSys_append_derVar,   &
    &                                 tem_varSys_proc_point,      &
    &                                 tem_varSys_proc_element
  use tem_dyn_array_module,     only: PositionOfVal

  use tem_utestEnv_module,      only: load_env

  !mpi!nprocs = 1

  implicit none

  character, parameter :: nl = C_NEW_LINE

  character(len=solSpecLen), parameter :: sysConf =              &
    &    'variable = {' // nl                                    &
    & // '  {' // nl                                             &
    & // '    name = "velocity",' // nl                          &
    & // '    ncomponents = 3,' // nl                            &
    & // '    vartype = "st_fun",' // nl                         &
    & // '    st_fun = {' // nl                                  &
    & // '      const = {1,2,3}' // nl                           &
    & // '    }' // nl                                           &
    & // '  },' // nl                                            &
    & // '  {' // nl                                             &
    & // '    name = "density",' // nl                           &
    & // '    ncomponents = 1,' // nl                            &
    & // '    vartype = "st_fun",' // nl                         &
    & // '    st_fun = {' // nl                                  &
    & // '      const = 4' // nl                                 &
    & // '    }' // nl                                           &
    & // '  },' // nl                                            &
    & // '  {' // nl                                             &
    & // '    name = "all",' // nl                               &
    & // '    ncomponents = 4,' // nl                            &
    & // '    vartype = "operation",' // nl                      &
    & // '    operation = {' // nl                               &
    & // '      kind = "combine",' // nl                         &
    & // '      input_varname = { "density", "velocity" }' // nl &
    & // '    }' // nl                                           &
    & // '  },' // nl                                            &
    & // '}' // nl

  type solver_type
    integer :: nDofs
    type(treelmesh_type) :: tree
    type(tem_bc_prop_type) :: boundary
    type(tem_general_type) :: general
  end type solver_type

  integer :: nComponents, nDofs, nElemsToTrack, pos
  type(solver_type), target :: solver
  type(tem_varSys_type) :: varSys
  type(tem_st_fun_linkedList_type) :: st_funList
  type(tem_variable_type), allocatable :: newVar(:)
  real(kind=rk), allocatable :: res(:), error(:)
  integer, allocatable :: vError(:)

  write(*,*) 'Hello from tem_variable_combine_Test'
  ! load utest mesh
  call load_env( tree     = solver%tree,     &
    &            boundary = solver%boundary, &
    &            general  = solver%general   )

  call tem_logging_init( level = 10, rank = 0 )
  write(*,*) 'nElems ', solver%tree%nElems

  allocate(solver%general%solver%conf(1))
  call open_config_chunk(L=solver%general%solver%conf(1), chunk=trim(sysConf))

  call tem_variable_load( me     = newVar,                        &
    &                     conf   = solver%general%solver%conf(1), &
    &                     key    = 'variable',                    &
    &                     vError = vError                         )
  call close_config(L=solver%general%solver%conf(1))

  nElemsToTrack = 1
  nComponents = 1
  nDofs = 1

  ! initialize variable system
  write(*,*) 'calling varsys init'
  call tem_varSys_init(me = varSys, systemName = 'utest')

  call tem_varSys_append_luaVar( luaVar     = newVar,    &
    &                            varSys     = varSys,    &
    &                            st_funList = st_funList )

  call tem_create_subTree_of_st_funList( &
    & me      = st_funList,              &
    & tree    = solver%tree,             &
    & bc_prop = solver%boundary          )

  pos = positionOfVal( varSys%varname, 'all' )
  if (pos==0) then
    write(*,*) 'FAILED: variable "all" not found'
    stop
  end if

  nComponents = varSys%method%val(pos)%nComponents
  allocate(res(nElemsToTrack*nComponents*nDofs))
  call varSys%method%val(pos)%get_element(     &
    & varSys  = varSys,                        &
    & elemPos = (/ 1 /),                       &
    & time    = solver%general%simControl%now, &
    & tree    = solver%tree,                   &
    & nElems  = nElemsToTrack,                 &
    & nDofs   = nDofs,                         &
    & res     = res                            )
  write(*,*) res
  allocate(error(nElemsToTrack*nComponents*nDofs))
  error = res - (/ 4, 1, 2, 3 /)
  if( any(error > eps) ) then
    write(*,*) 'FAILED: Variable "all" returned wrong value.'
    stop
  endif

  ! for the sake of the address sanitizer
  deallocate(error, res, vError, newVar)

  write(*,*) 'PASSED'

  call fin_env()

end program tem_variable_combine_Test

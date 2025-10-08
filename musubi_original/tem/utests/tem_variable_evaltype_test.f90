! Copyright (c) 2015-2018 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2016, 2018-2019 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2016 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2018 Robin Weihe <robin.weihe@student.uni-siegen.de>
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
program tem_variable_evaltype_test

  use, intrinsic :: iso_c_binding, only: C_NEW_LINE, c_ptr, c_loc, c_f_pointer

  use env_module,               only: rk, solSpecLen, fin_env

  use aotus_module,             only: open_config_chunk, &
    &                                 close_config,      &
    &                                 flu_state
  use aot_table_module,         only: aot_table_open,   &
    &                                 aot_table_close,  &
    &                                 aot_table_length, &
    &                                 aot_get_val

  use tem_logging_module,       only: tem_logging_init
  use tem_float_module,         only: operator(.fne.)
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

  character(len=solSpecLen), parameter :: sysConf =      &
    &    'function varfun(x, y, z, t)' // nl             &
    & // '  return 5' // nl                              &
    & // 'end' // nl                                     &
    & // 'firstelem = {' // nl                           &
    & // '  kind = "canoND",' // nl                      &
    & // '  object = {' // nl                            &
    & // '    origin = {0.25, 0.25, 0.25}' // nl         &
    & // '  }' // nl                                     &
    & // '}' // nl                                       &
    & // 'variable = {' // nl                            &
    & // '  {' // nl                                     &
    & // '    name = "3",' // nl                         &
    & // '    ncomponents = 1,'// nl                     &
    & // '    vartype = "st_fun",'// nl                  &
    & // '    st_fun = {' // nl                          &
    & // '      { const = {1} },' // nl                  &
    & // '      { const = {2} }' // nl                   &
    & // '    }' // nl                                   &
    & // '  },' // nl                                    &
    & // '  {' // nl                                     &
    & // '    name = "explicitly_3",' // nl              &
    & // '    ncomponents = 1,' // nl                    &
    & // '    vartype = "st_fun",' // nl                 &
    & // '    evaltype = "add",' // nl                   &
    & // '    st_fun = {' // nl                          &
    & // '      { const = {1} },' // nl                  &
    & // '      { const = {2} }' // nl                   &
    & // '    }' // nl                                   &
    & // '  },' // nl                                    &
    & // '  {' // nl                                     &
    & // '    name = "just_2",' // nl                    &
    & // '    ncomponents = 1,' // nl                    &
    & // '    vartype = "st_fun",' // nl                 &
    & // '    evaltype = "first",' // nl                 &
    & // '    st_fun = {' // nl                          &
    & // '      { const = {1} },' // nl                  &
    & // '      { const = {2} }' // nl                   &
    & // '    }' // nl                                   &
    & // '  },' // nl                                    &
    & // '  {' // nl                                     &
    & // '    name = "shaped_3",' // nl                  &
    & // '    ncomponents = 1,' // nl                    &
    & // '    vartype = "st_fun",' // nl                 &
    & // '    evaltype = "add",' // nl                   &
    & // '    st_fun = {' // nl                          &
    & // '      {' // nl                                 &
    & // '        const = {1},' // nl                    &
    & // '        shape = firstelem' // nl               &
    & // '      },' // nl                                &
    & // '      { const = {2} }' // nl                   &
    & // '    }' // nl                                   &
    & // '  },' // nl                                    &
    & // '  {' // nl                                     &
    & // '    name = "shaped_reverse_3",' // nl          &
    & // '    ncomponents = 1,' // nl                    &
    & // '    vartype = "st_fun",' // nl                 &
    & // '    evaltype = "add",' // nl                   &
    & // '    st_fun = {' // nl                          &
    & // '      { const = {2} },' // nl                  &
    & // '      {' // nl                                 &
    & // '        const = {1},' // nl                    &
    & // '        shape = firstelem' // nl               &
    & // '      }' // nl                                 &
    & // '    }' // nl                                   &
    & // '  },' // nl                                    &
    & // '  {' // nl                                     &
    & // '    name = "var_5",' // nl                     &
    & // '    ncomponents = 1,' // nl                    &
    & // '    vartype = "st_fun",' // nl                 &
    & // '    evaltype = "add",' // nl                   &
    & // '    st_fun = {' // nl                          &
    & // '      {' // nl                                 &
    & // '        fun = varfun,' // nl                   &
    & // '        shape = firstelem' // nl               &
    & // '      },' // nl                                &
    & // '      { const = 2 }' // nl                     &
    & // '    }' // nl                                   &
    & // '  },' // nl                                    &
    & // '  {' // nl                                     &
    & // '    name = "var_5_reverse",' // nl             &
    & // '    ncomponents = 1,' // nl                    &
    & // '    vartype = "st_fun",' // nl                 &
    & // '    evaltype = "add",' // nl                   &
    & // '    st_fun = {' // nl                          &
    & // '      { const = 2 },' // nl                    &
    & // '      {' // nl                                 &
    & // '        fun = varfun,' // nl                   &
    & // '        shape = firstelem' // nl               &
    & // '      }' // nl                                 &
    & // '    }' // nl                                   &
    & // '  },' // nl                                    &
    & // '  {' // nl                                     &
    & // '    name = "first_5",' // nl                   &
    & // '    ncomponents = 1,' // nl                    &
    & // '    vartype = "st_fun",' // nl                 &
    & // '    evaltype = "first",' // nl                 &
    & // '    st_fun = {' // nl                          &
    & // '      {' // nl                                 &
    & // '        fun = varfun,' // nl                   &
    & // '        shape = firstelem' // nl               &
    & // '      },' // nl                                &
    & // '      { const = 2 }' // nl                     &
    & // '    }' // nl                                   &
    & // '  }' // nl                                     &
    & // '}' // nl

  type solver_type
    integer :: nDofs
    real(kind=rk), allocatable :: state(:)
    type(treelmesh_type) :: tree
    type(tem_bc_prop_type) :: boundary
    type(tem_general_type) :: general
  end type solver_type

  integer :: nComponents, nDofs, nElemsToTrack, pos
  type(solver_type), target :: solver
  type(tem_varSys_type) :: varSys
  type(tem_st_fun_linkedList_type) :: st_funList
  type(tem_variable_type), allocatable :: newVar(:)
  real(kind=rk), allocatable :: res(:)
  integer, allocatable :: vError(:)

  write(*,*) 'Hello from tem_variable_evaltype_test'
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
  ! for the sake of the address sanitizer
  deallocate(vError)

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

  allocate(res(nElemsToTrack*nComponents*nDofs))
  pos = positionOfVal( varSys%varname, "3" )
  call varSys%method%val(pos)%get_element(     &
    & varSys  = varSys,                        &
    & elemPos = (/ 1 /),                       &
    & time    = solver%general%simControl%now, &
    & tree    = solver%tree,                   &
    & nElems  = nElemsToTrack,                 &
    & nDofs   = nDofs,                         &
    & res     = res                            )
  write(*,*) res
  if( res(1) .fne. 3._rk) then
    write(*,*) 'FAILED: Variable 3 returned wrong value.'
    stop
  endif

  pos = positionOfVal( varSys%varname, "explicitly_3" )
  call varSys%method%val(pos)%get_element(     &
    & varSys  = varSys,                        &
    & elemPos = (/ 1 /),                       &
    & time    = solver%general%simControl%now, &
    & tree    = solver%tree,                   &
    & nElems  = nElemsToTrack,                 &
    & nDofs   = nDofs,                         &
    & res     = res                            )
  write(*,*) res
  if( res(1) .fne. 3._rk ) then
    write(*,*) 'FAILED: Variable explicitly_3 returned wrong value.'
    stop
  endif

  pos = positionOfVal( varSys%varname, "just_2" )
  call varSys%method%val(pos)%get_element(     &
    & varSys  = varSys,                        &
    & elemPos = (/ 1 /),                       &
    & time    = solver%general%simControl%now, &
    & tree    = solver%tree,                   &
    & nElems  = nElemsToTrack,                 &
    & nDofs   = nDofs,                         &
    & res     = res                            )
  write(*,*) res
  if( res(1) .fne. 1._rk ) then
    write(*,*) 'FAILED: Variable just_2 returned wrong value.'
    stop
  endif

  ! Check a function where only the first element
  ! should have the value 3, all others the value
  ! 2.
  pos = positionOfVal( varSys%varname, "shaped_3" )
  call varSys%method%val(pos)%get_element(     &
    & varSys  = varSys,                        &
    & elemPos = (/ 1 /),                       &
    & time    = solver%general%simControl%now, &
    & tree    = solver%tree,                   &
    & nElems  = nElemsToTrack,                 &
    & nDofs   = nDofs,                         &
    & res     = res                            )
  write(*,*) res
  if ( res(1) .fne. 3._rk ) then
    write(*,*) 'FAILED: Variable shaped_3 returned wrong value' &
      &        //' in first element.'
    stop
  end if
  call varSys%method%val(pos)%get_element(     &
    & varSys  = varSys,                        &
    & elemPos = (/ 2 /),                       &
    & time    = solver%general%simControl%now, &
    & tree    = solver%tree,                   &
    & nElems  = nElemsToTrack,                 &
    & nDofs   = nDofs,                         &
    & res     = res                            )
  write(*,*) res
  if ( res(1) .fne. 2._rk ) then
    write(*,*) 'FAILED: Variable shaped_3 returned wrong value' &
      &        //' in second element.'
    stop
  end if

  ! Check a function where only the first element
  ! should have the value 3, all others the value
  ! 2, this time reversed ordering of ST fun definitions.
  pos = positionOfVal( varSys%varname, "shaped_reverse_3" )
  call varSys%method%val(pos)%get_element(     &
    & varSys  = varSys,                        &
    & elemPos = (/ 1 /),                       &
    & time    = solver%general%simControl%now, &
    & tree    = solver%tree,                   &
    & nElems  = nElemsToTrack,                 &
    & nDofs   = nDofs,                         &
    & res     = res                            )
  write(*,*) res
  if ( res(1) .fne. 3._rk ) then
    write(*,*) 'FAILED: Variable shaped_reverse_3 returned wrong value' &
      &        //' in first element.'
    stop
  end if
  call varSys%method%val(pos)%get_element(     &
    & varSys  = varSys,                        &
    & elemPos = (/ 2 /),                       &
    & time    = solver%general%simControl%now, &
    & tree    = solver%tree,                   &
    & nElems  = nElemsToTrack,                 &
    & nDofs   = nDofs,                         &
    & res     = res                            )
  write(*,*) res
  if ( res(1) .fne. 2._rk ) then
    write(*,*) 'FAILED: Variable shaped_reverse_3 returned wrong value' &
      &        //' in second element.'
    stop
  end if

  ! Check a function where only the first element should have the
  ! pseudo-variable value 7 (5+2), all others the value 2.
  pos = positionOfVal( varSys%varname, "var_5" )
  call varSys%method%val(pos)%get_element(     &
    & varSys  = varSys,                        &
    & elemPos = (/ 1 /),                       &
    & time    = solver%general%simControl%now, &
    & tree    = solver%tree,                   &
    & nElems  = nElemsToTrack,                 &
    & nDofs   = nDofs,                         &
    & res     = res                            )
  write(*,*) res
  if( res(1) .fne. 7._rk ) then
    write(*,*) 'FAILED: Variable var_5 returned wrong value in first element.'
    stop
  endif
  call varSys%method%val(pos)%get_element(     &
    & varSys  = varSys,                        &
    & elemPos = (/ 2 /),                       &
    & time    = solver%general%simControl%now, &
    & tree    = solver%tree,                   &
    & nElems  = nElemsToTrack,                 &
    & nDofs   = nDofs,                         &
    & res     = res                            )
  write(*,*) res
  if( res(1) .fne. 2._rk ) then
    write(*,*) 'FAILED: Variable var_5 returned wrong value in second element.'
    stop
  endif

  ! Check a function where only the first element should have the
  ! pseudo-variable value 7 (5+2), all others the value 2. This time with
  ! reversed stfuns
  pos = positionOfVal( varSys%varname, "var_5_reverse" )
  call varSys%method%val(pos)%get_element(     &
    & varSys  = varSys,                        &
    & elemPos = (/ 1 /),                       &
    & time    = solver%general%simControl%now, &
    & tree    = solver%tree,                   &
    & nElems  = nElemsToTrack,                 &
    & nDofs   = nDofs,                         &
    & res     = res                            )
  write(*,*) res
  if( res(1) .fne. 7._rk ) then
    write(*,*) 'FAILED: Variable var_5_reverse returned wrong value in first element.'
    stop
  endif
  call varSys%method%val(pos)%get_element(     &
    & varSys  = varSys,                        &
    & elemPos = (/ 2 /),                       &
    & time    = solver%general%simControl%now, &
    & tree    = solver%tree,                   &
    & nElems  = nElemsToTrack,                 &
    & nDofs   = nDofs,                         &
    & res     = res                            )
  write(*,*) res
  if( res(1) .fne. 2._rk ) then
    write(*,*) 'FAILED: Variable var_5_reverse returned wrong value in second element.'
    stop
  endif

  ! Check a function where only the first element should have the
  ! pseudo-variable value 5, all others the value 2.
  pos = positionOfVal( varSys%varname, "first_5" )
  call varSys%method%val(pos)%get_element(     &
    & varSys  = varSys,                        &
    & elemPos = (/ 1 /),                       &
    & time    = solver%general%simControl%now, &
    & tree    = solver%tree,                   &
    & nElems  = nElemsToTrack,                 &
    & nDofs   = nDofs,                         &
    & res     = res                            )
  write(*,*) res
  if( res(1) .fne. 5._rk ) then
    write(*,*) 'FAILED: Variable first_5 returned wrong value in first element.'
    stop
  endif
  call varSys%method%val(pos)%get_element(     &
    & varSys  = varSys,                        &
    & elemPos = (/ 2 /),                       &
    & time    = solver%general%simControl%now, &
    & tree    = solver%tree,                   &
    & nElems  = nElemsToTrack,                 &
    & nDofs   = nDofs,                         &
    & res     = res                            )
  write(*,*) res
  if( res(1) .fne. 2._rk ) then
    write(*,*) 'FAILED: Variable first_5 returned wrong value in second element.'
    stop
  endif

  ! for the sake of the address sanitizer
  deallocate(res)
  call close_config(L=solver%general%solver%conf(1))

  write(*,*) 'PASSED'
  call fin_env()

end program tem_variable_evaltype_test

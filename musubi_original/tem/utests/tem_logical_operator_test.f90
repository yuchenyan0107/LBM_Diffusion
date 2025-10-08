! Copyright (c) 2016 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2019 Harald Klimach <harald.klimach@uni-siegen.de>
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
program tem_logical_opertor_test

  use, intrinsic :: iso_c_binding,      only: C_NEW_LINE
  use env_module,                       only: rk, solSpecLen, fin_env

  use aotus_module,             only: open_config_chunk, &
    &                                 flu_state

  use tem_float_module,                 only: operator(.feq.), &
    &                                         operator(.fgt.), &
    &                                         operator(.flt.)
  use tem_logical_operation_var_module, only: logicalToReal,      &
    &                                         logicalToRealArray, &
    &                                         realToLogical,      &
    &                                         realToLogicalArray, &
    &                                         numTrue,            &
    &                                         numFalse
  use tem_utestEnv_module,      only: load_env
  use tem_bc_prop_module,       only: tem_bc_prop_type
  use tem_general_module,       only: tem_general_type
  use tem_varSys_module,        only: tem_varSys_init, &
    &                                 tem_varSys_type
  use treelmesh_module,         only: treelmesh_type
  use tem_variable_module,      only: tem_variable_type, &
    &                                 tem_variable_load
  use tem_dyn_array_module,     only: PositionOfVal
  use tem_derived_module,       only: tem_varSys_append_luaVar
  use tem_spacetime_fun_module, only: tem_st_fun_linkedList_type,      &
    &                                 tem_create_subTree_of_st_funList
  use tem_geometry_module,      only: tem_BaryOfId

  !mpi!nprocs = 1

  implicit none

  character, parameter :: nl = C_NEW_LINE
  character(len=solSpecLen), parameter :: sysConf =                 &
    &    'variable = {'                                       // nl &
    & // '  {'                                                // nl &
    & // '    name = "true_var",'                             // nl &
    & // '    ncomponents = 1,'                               // nl &
    & // '    vartype = "st_fun",'                            // nl &
    & // '    st_fun = 1.0'                                   // nl &
    & // '  },'                                               // nl &
    & // '  {'                                                // nl &
    & // '    name = "false_var",'                            // nl &
    & // '    ncomponents = 1,'                               // nl &
    & // '    vartype = "st_fun",'                            // nl &
    & // '    st_fun = 0.0'                                   // nl &
    & // '  },'                                               // nl &
    & // '  {'                                                // nl &
    & // '    name = "true_and_true",'                        // nl &
    & // '    ncomponents = 1,'                               // nl &
    & // '    vartype = "operation",'                         // nl &
    & // '    operation = {'                                  // nl &
    & // '      kind = "and",'                                // nl &
    & // '      input_varname = { "true_var", "true_var" }'   // nl &
    & // '    }'                                              // nl &
    & // '  },'                                               // nl &
    & // '  {'                                                // nl &
    & // '    name = "true_and_false",'                       // nl &
    & // '    ncomponents = 1,'                               // nl &
    & // '    vartype = "operation",'                         // nl &
    & // '    operation = {'                                  // nl &
    & // '      kind = "and",'                                // nl &
    & // '      input_varname = { "true_var", "false_var" }'  // nl &
    & // '    }'                                              // nl &
    & // '  },'                                               // nl &
    & // '  {'                                                // nl &
    & // '    name = "true_or_true",'                         // nl &
    & // '    ncomponents = 1,'                               // nl &
    & // '    vartype = "operation",'                         // nl &
    & // '    operation = {'                                  // nl &
    & // '      kind = "or",'                                 // nl &
    & // '      input_varname = { "true_var", "true_var" }'   // nl &
    & // '    }'                                              // nl &
    & // '  },'                                               // nl &
    & // '  {'                                                // nl &
    & // '    name = "true_or_false",'                        // nl &
    & // '    ncomponents = 1,'                               // nl &
    & // '    vartype = "operation",'                         // nl &
    & // '    operation = {'                                  // nl &
    & // '      kind = "or",'                                 // nl &
    & // '      input_varname = { "true_var", "false_var" }'  // nl &
    & // '    }'                                              // nl &
    & // '  },'                                               // nl &
    & // '  {'                                                // nl &
    & // '    name = "false_or_false",'                       // nl &
    & // '    ncomponents = 1,'                               // nl &
    & // '    vartype = "operation",'                         // nl &
    & // '    operation = {'                                  // nl &
    & // '      kind = "or",'                                 // nl &
    & // '      input_varname = { "false_var", "false_var" }' // nl &
    & // '    }'                                              // nl &
    & // '  },'                                               // nl &
    & // '}'

  logical :: res = .true.

  call check_logicalToReal(res)
  call check_realToLogical(res)
  call check_arrayRoutines(res)
  call check_variableOperations(res)

  if (res) write(*,*) 'PASSED'
  call fin_env()


contains


  subroutine checkResult( res, msg )
    logical, intent(in) :: res
    character(len=*), intent(in) :: msg

    if ( .not. res ) then
      write(*,*) msg
      stop
    end if
  end subroutine checkResult

  subroutine check_logicalToReal(res)
    logical, intent(inout) :: res

    res = res .and. (logicalToReal(.true.) .feq. numTrue)
    call checkResult( res, 'Test logicalToReal(.true.) .feq. numTrue failed' )
    res = res .and..not. (logicalToReal(.true.) .feq. numFalse)
    call checkResult( res, 'Test logicalToReal(.true.) .feq. numFalse failed' )
    res = res .and. (logicalToReal(.false.) .feq. numFalse)
    call checkResult( res, 'Test logicalToReal(.false.) .feq. numFalse failed' )
    res = res .and..not. (logicalToReal(.false.) .feq. numTrue)
    call checkResult( res, 'Test logicalToReal(.false.) .feq. numTrue failed' )

  end subroutine check_logicalToReal

  subroutine check_realToLogical(res)
    logical, intent(inout) :: res

    res = res .and..not. realToLogical(numFalse)
    call checkResult( res, 'Test realToLogical(numFalse) failed' )
    res = res .and. realToLogical(numTrue)
    call checkResult( res, 'Test realToLogical(numTrue) failed' )
    res = res .and. realToLogical(0._rk+spacing(0._rk))
    call checkResult( res, 'Test realToLogical(0._rk+spacing(0._rk)) failed' )
    res = res .and. realToLogical(huge(0._rk))
    call checkResult( res, 'Test realToLogical(huge(0._rk)) failed' )
    res = res .and. realToLogical(tiny(0._rk))
    call checkResult( res, 'Test realToLogical(tiny(0._rk)) failed' )

  end subroutine check_realToLogical

  subroutine check_arrayRoutines(res)
    logical, intent(inout) :: res
    real(kind=rk) :: logicalAsReals(5)

    logicalAsReals = logicalToRealArray(                      &
      & value = (/ .true., .true., .true., .true., .true. /), &
      & n     = 5                                             )
    res = res .and. all(                                                  &
      & pack( array = realToLogicalArray( value = logicalAsReals, n = 5), &
      &       mask  = (/ .true., .true., .true., .true., .true. /)      ) )
    call checkResult( res, 'Test array with all true failed' )


    logicalAsReals = logicalToRealArray(                      &
      & value = (/ .false., .true., .false., .true., .false. /), &
      & n     = 5                                             )
    res = res .and. all(                                                  &
      & pack( array = realToLogicalArray( value = logicalAsReals, n = 5), &
      &       mask  = (/ .false., .true., .false., .true., .false. /)   ) )
    call checkResult( res, 'Test mixed array for true failed' )
    res = res .and. .not. any(                                            &
      & pack( array = realToLogicalArray( value = logicalAsReals, n = 5), &
      &       mask  = (/ .true., .false., .true., .false., .true. /)    ) )
    call checkResult( res, 'Test mixed array for false failed' )

    logicalAsReals = logicalToRealArray(                           &
      & value = (/ .false., .false., .false., .false., .false. /), &
      & n     = 5                                                  )
    res = res .and. .not. any(                                            &
      & pack( array = realToLogicalArray( value = logicalAsReals, n = 5), &
      &       mask  = (/ .true., .true., .true., .true., .true. /)      ) )
    call checkResult( res, 'Test array with all false failed' )

  end subroutine check_arrayRoutines

  subroutine check_variableOperations(res)
    logical, intent(inout) :: res

    type(treelmesh_type) :: tree
    type(tem_bc_prop_type) :: boundary
    type(tem_general_type) :: general
    type(tem_varSys_type) :: varSys
    type(tem_st_fun_linkedList_type) :: st_funList
    type(flu_state) :: conf
    type(tem_variable_type), allocatable :: newVar(:)
    real(kind=rk) :: point(1,3)
    integer, allocatable :: indices(:,:)
    integer, allocatable :: vError(:)
    integer :: iVar

    call load_env( tree     = tree,     &
      &            boundary = boundary, &
      &            general  = general   )

    call open_config_chunk(L=conf, chunk=trim(sysConf))

    call tem_variable_load( me   = newVar,     &
      &                     conf = conf,       &
      &                     key  = 'variable', &
      &                     vError  = vError   )

    call tem_varSys_init(me = varSys, systemName = 'utest')

    call tem_varSys_append_luaVar( luaVar     = newVar,    &
      &                            varSys     = varSys,    &
      &                            st_funList = st_funList )

    call tem_create_subTree_of_st_funList( &
      & me      = st_funList,              &
      & tree    = tree,                    &
      & bc_prop = boundary                 )

    point(1,:) = tem_BaryOfId( tree   = tree,          &
      &                        treeID = tree%treeID(1) )

    allocate(indices(varSys%method%nVals,1))
    do iVar = 1, varSys%method%nVals
      call varSys%method%val(iVar)%setup_indices( &
        & varSys = varSys,                        &
        & point  = point,                         &
        & iLevel = 1,                             &
        & tree   = tree,                          &
        & nPnts  = 1,                             &
        & idx    = indices(iVar,:)                  )
    end do

    call check_variableOperations_byPoint( res, general, tree, varSys, point )
    call check_variableOperations_byElement( res, general, tree, varSys )
    call check_variableOperations_byIndex( res, general, varSys, indices )

  end subroutine check_variableOperations

  subroutine check_variableOperations_byPoint(res, general, tree, varSys, point)
    logical, intent(inout) :: res
    type(treelmesh_type), intent(in) :: tree
    type(tem_general_type), intent(in) :: general
    type(tem_varSys_type), intent(in) :: varSys
    real(kind=rk), intent(in) :: point(1,3)

    real(kind=rk) :: resVar(1)
    integer :: varPos

    varPos = positionOfVal(varSys%varname, 'true_and_true')
    call varSys%method%val(varPos)%get_point( &
      & varSys = varSys,                      &
      & point  = point,                       &
      & time   = general%simControl%now,      &
      & tree   = tree,                        &
      & nPnts  = 1,                           &
      & res    = resVar                       )
    write(*,*) resVar
    res = res .and. realToLogical(resVar(1))
    call checkResult( res, 'Logical operator: true_and_true' )

    varPos = positionOfVal(varSys%varname, 'true_and_false')
    call varSys%method%val(varPos)%get_point( &
      & varSys = varSys,                      &
      & point  = point,                       &
      & time   = general%simControl%now,      &
      & tree   = tree,                        &
      & nPnts  = 1,                           &
      & res    = resVar                       )
    write(*,*) resVar
    res = res .and..not. realToLogical(resVar(1))
    call checkResult( res, 'Logical operator: true_and_false' )

    varPos = positionOfVal(varSys%varname, 'true_or_true')
    call varSys%method%val(varPos)%get_point( &
      & varSys = varSys,                      &
      & point  = point,                       &
      & time   = general%simControl%now,      &
      & tree   = tree,                        &
      & nPnts  = 1,                           &
      & res    = resVar                       )
    write(*,*) resVar
    res = res .and. realToLogical(resVar(1))
    call checkResult( res, 'Logical operator: true_or_true' )

    varPos = positionOfVal(varSys%varname, 'true_or_false')
    call varSys%method%val(varPos)%get_point( &
      & varSys = varSys,                      &
      & point  = point,                       &
      & time   = general%simControl%now,      &
      & tree   = tree,                        &
      & nPnts  = 1,                           &
      & res    = resVar                       )
    write(*,*) resVar
    res = res .and. realToLogical(resVar(1))
    call checkResult( res, 'Logical operator: true_or_false' )

    varPos = positionOfVal(varSys%varname, 'false_or_false')
    call varSys%method%val(varPos)%get_point( &
      & varSys = varSys,                      &
      & point  = point,                       &
      & time   = general%simControl%now,      &
      & tree   = tree,                        &
      & nPnts  = 1,                           &
      & res    = resVar                       )
    write(*,*) resVar
    res = res .and..not. realToLogical(resVar(1))
    call checkResult( res, 'Logical operator: false_or_false' )

  end subroutine check_variableOperations_byPoint

  subroutine check_variableOperations_byElement(res, general, tree, varSys)
    logical, intent(inout) :: res
    type(treelmesh_type), intent(in) :: tree
    type(tem_general_type), intent(in) :: general
    type(tem_varSys_type), intent(in) :: varSys

    real(kind=rk) :: resVar(1)
    integer :: varPos

    varPos = positionOfVal(varSys%varname, 'true_and_true')
    call varSys%method%val(varPos)%get_element( &
      & varSys  = varSys,                       &
      & elemPos = (/ 1 /),                      &
      & time    = general%simControl%now,       &
      & tree    = tree,                         &
      & nElems  = 1,                            &
      & nDofs   = 1,                            &
      & res     = resVar                        )
    write(*,*) resVar
    res = res .and. realToLogical(resVar(1))
    call checkResult( res, 'Logical operator: true_and_true' )

    varPos = positionOfVal(varSys%varname, 'true_and_false')
    call varSys%method%val(varPos)%get_element( &
      & varSys  = varSys,                       &
      & elemPos = (/ 1 /),                      &
      & time    = general%simControl%now,       &
      & tree    = tree,                         &
      & nElems  = 1,                            &
      & nDofs   = 1,                            &
      & res     = resVar                        )
    write(*,*) resVar
    res = res .and..not. realToLogical(resVar(1))
    call checkResult( res, 'Logical operator: true_and_false' )

    varPos = positionOfVal(varSys%varname, 'true_or_true')
    call varSys%method%val(varPos)%get_element( &
      & varSys  = varSys,                       &
      & elemPos = (/ 1 /),                      &
      & time    = general%simControl%now,       &
      & tree    = tree,                         &
      & nElems  = 1,                            &
      & nDofs   = 1,                            &
      & res     = resVar                        )
    write(*,*) resVar
    res = res .and. realToLogical(resVar(1))
    call checkResult( res, 'Logical operator: true_or_true' )

    varPos = positionOfVal(varSys%varname, 'true_or_false')
    call varSys%method%val(varPos)%get_element( &
      & varSys  = varSys,                       &
      & elemPos = (/ 1 /),                      &
      & time    = general%simControl%now,       &
      & tree    = tree,                         &
      & nElems  = 1,                            &
      & nDofs   = 1,                            &
      & res     = resVar                        )
    write(*,*) resVar
    res = res .and. realToLogical(resVar(1))
    call checkResult( res, 'Logical operator: true_or_false' )

    varPos = positionOfVal(varSys%varname, 'false_or_false')
    call varSys%method%val(varPos)%get_element( &
      & varSys  = varSys,                       &
      & elemPos = (/ 1 /),                      &
      & time    = general%simControl%now,       &
      & tree    = tree,                         &
      & nElems  = 1,                            &
      & nDofs   = 1,                            &
      & res     = resVar                        )
    write(*,*) resVar
    res = res .and..not. realToLogical(resVar(1))
    call checkResult( res, 'Logical operator: false_or_false' )

  end subroutine check_variableOperations_byElement

  subroutine check_variableOperations_byIndex(res, general, varSys, indices)
    logical, intent(inout) :: res

    type(tem_general_type), intent(in) :: general
    type(tem_varSys_type), intent(in) :: varSys
    integer, intent(in) :: indices(:,:)

    real(kind=rk) :: resVar(1)
    integer :: varPos

    varPos = positionOfVal(varSys%varname, 'true_and_true')
    call varSys%method%val(varPos)%get_valOfIndex( &
      & varSys = varSys,                           &
      & time   = general%simControl%now,           &
      & iLevel = 1,                                &
      & idx    = indices(varPos,:),                &
      & nVals  = 1,                                &
      & res    = resVar                            )
    write(*,*) resVar
    res = res .and. realToLogical(resVar(1))
    call checkResult( res, 'Logical operator: true_and_true' )

    varPos = positionOfVal(varSys%varname, 'true_and_false')
    call varSys%method%val(varPos)%get_valOfIndex( &
      & varSys = varSys,                           &
      & time   = general%simControl%now,           &
      & iLevel = 1,                                &
      & idx    = indices(varPos,:),                &
      & nVals  = 1,                                &
      & res    = resVar                            )
    write(*,*) resVar
    res = res .and..not. realToLogical(resVar(1))
    call checkResult( res, 'Logical operator: true_and_false' )

    varPos = positionOfVal(varSys%varname, 'true_or_true')
    call varSys%method%val(varPos)%get_valOfIndex( &
      & varSys = varSys,                           &
      & time   = general%simControl%now,           &
      & iLevel = 1,                                &
      & idx    = indices(varPos,:),                &
      & nVals  = 1,                                &
      & res    = resVar                            )
    write(*,*) resVar
    res = res .and. realToLogical(resVar(1))
    call checkResult( res, 'Logical operator: true_or_true' )

    varPos = positionOfVal(varSys%varname, 'true_or_false')
    call varSys%method%val(varPos)%get_valOfIndex( &
      & varSys = varSys,                           &
      & time   = general%simControl%now,           &
      & iLevel = 1,                                &
      & idx    = indices(varPos,:),                &
      & nVals  = 1,                                &
      & res    = resVar                            )
    write(*,*) resVar
    res = res .and. realToLogical(resVar(1))
    call checkResult( res, 'Logical operator: true_or_false' )

    varPos = positionOfVal(varSys%varname, 'false_or_false')
    call varSys%method%val(varPos)%get_valOfIndex( &
      & varSys = varSys,                           &
      & time   = general%simControl%now,           &
      & iLevel = 1,                                &
      & idx    = indices(varPos,:),                &
      & nVals  = 1,                                &
      & res    = resVar                            )
    write(*,*) resVar
    res = res .and..not. realToLogical(resVar(1))
    call checkResult( res, 'Logical operator: false_or_false' )

  end subroutine check_variableOperations_byIndex


end program tem_logical_opertor_test

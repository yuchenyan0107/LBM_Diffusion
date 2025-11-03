! Copyright (c) 2014, 2016 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2014, 2016, 2019 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2014-2015, 2017 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2014-2015 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2016-2017 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
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
!> This UTEST is to test variable system.
!! load variables and access state and derive quantities from
!! c pointer methoddata
program tem_varSys_test
  use, intrinsic :: iso_c_binding, only: C_NEW_LINE, c_ptr, c_loc, c_f_pointer

  use env_module,               only: rk, solSpecLen, labelLen, eps, fin_env
  use tem_bc_prop_module,       only: tem_bc_prop_type
  use tem_general_module,       only: tem_general_type
  use treelmesh_module,         only: treelmesh_type
  use tem_time_module,          only: tem_time_type
  use tem_utestEnv_module,      only: load_env
  use tem_tools_module,         only: upper_to_lower, tem_PositionInSorted
  use tem_varSys_module,        only: tem_varSys_init, tem_varSys_type,        &
    &                                 tem_varSys_op_type,                      &
    &                                 tem_varSys_append_stateVar,              &
    &                                 tem_varSys_append_derVar,                &
    &                                 tem_varSys_proc_point,                   &
    &                                 tem_varSys_proc_element,                 &
    &                                 tem_varSys_proc_setParams,               &
    &                                 tem_varSys_proc_getParams,               &
    &                                 tem_varSys_proc_setupIndices,            &
    &                                 tem_varSys_proc_getValOfIndex
  use tem_variable_module,      only: tem_variable_type, tem_variable_load
  use tem_grow_array_module,    only: grw_intArray_type, init, append, destroy
  use tem_dyn_array_module,     only: PositionOfVal
  use tem_spacetime_fun_module, only: tem_spacetime_fun_type,                  &
    &                                 tem_load_spacetime,                      &
    &                                 tem_st_fun_listElem_type,                &
    &                                 tem_st_fun_linkedList_type, append,      &
    &                                 tem_spacetime_for
  use tem_spacetime_var_module, only:                                      &
    &                           evaluate_add_spacetime_scalarByTreeID,     &
    &                           evaluate_add_spacetime_vectorByTreeID
  use tem_subTree_module,       only: tem_create_subTree_of
  use tem_stencil_module,       only: tem_stencilHeader_type

  use aotus_module,          only: open_config_chunk, close_config, flu_state
  use aot_table_module,      only: aot_table_open, aot_table_close,            &
    &                              aot_table_length, aot_get_val

  !mpi!nprocs = 1

  implicit none

  character, parameter :: nl = C_NEW_LINE

  character(len=solSpecLen), parameter :: sysConf = &
    &    'function vel_fun(x,y,z,t)' // nl          &
    & // '  return {x,y,z}' // nl                   &
    & // 'end' // nl                                &
    & // 'variables = {' // nl                      &
    & // '  {' // nl                                &
    & // '    name = "velocity",' // nl             &
    & // '    ncomponents = 3, '// nl               &
    & // '    vartype = "st_fun",'// nl             &
    & // '    st_fun = vel_fun' // nl               &
    & // '  }' // nl                                &
    & // '}' // nl                                  &
    & // 'track_variable = {' // nl                 &
    & // '  "state",'// nl                          &
    & // '  "density",'// nl                        &
    & // '  "velocity"' // nl                       &
    & // '}' // nl

  ! write output to screen
  logical, parameter :: dumpRes = .true.
  ! number of elements to track
  integer, parameter :: nElems_track = 5
  ! position in global treeID list to track
  integer, dimension(nElems_track), parameter :: elemPos = (/ 1, 3, 5, 7, 8 /)

  ! reference values
  real(kind=rk), dimension(nElems_track,4), parameter ::     &
    & state = reshape((/ 1.0_rk,  9.0_rk, 17.0_rk, 25.0_rk, 29.0_rk, &
    &                    2.0_rk, 10.0_rk, 18.0_rk, 26.0_rk, 30.0_rk, &
    &                    3.0_rk, 11.0_rk, 19.0_rk, 27.0_rk, 31.0_rk, &
    &                    4.0_rk, 12.0_rk, 20.0_rk, 28.0_rk, 32.0_rk /), &
    &         (/nElems_track,4/) )
  real(kind=rk), dimension(nElems_track), parameter :: &
    & dens = (/ 10.0_rk,  42.0_rk, 74.0_rk, 106.0_rk, 122.0_rk /)
  real(kind=rk), dimension(nElems_track,3), parameter ::                &
    & vel = reshape((/ 0.25_rk,  0.25_rk, 0.25_rk, 0.25_rk, 0.75_rk,    &
    &                  0.25_rk,  0.75_rk, 0.25_rk, 0.75_rk, 0.75_rk,    &
    &                  0.25_rk,  0.25_rk, 0.75_rk, 0.75_rk, 0.75_rk /), &
    &       (/nElems_track,3/) )

  type solver_type
    integer :: nDofs
    real(kind=rk), allocatable :: state(:)
    real(kind=rk), allocatable :: dens(:)
    type(treelmesh_type) :: tree
    type(tem_bc_prop_type) :: boundary
    type(tem_general_type) :: general
  end type solver_type

  type tracking_type
    integer :: nRequestedVars
    type(grw_intArray_type) :: varPos
    character(len=labelLen), allocatable :: variable(:)
  end type tracking_type

  type(tem_varSys_type) :: varSys
  procedure(tem_varSys_proc_point), pointer :: get_point => NULL()
  procedure(tem_varSys_proc_element), pointer :: get_element => NULL()
  procedure(tem_varSys_proc_setParams), pointer :: set_params => null()
  procedure(tem_varSys_proc_getParams), pointer :: get_params => null()
  procedure(tem_varSys_proc_setupIndices), pointer :: &
    &                                      setup_indices => null()
  procedure(tem_varSys_proc_getValOfIndex), pointer :: &
    &                                       get_valOfIndex => null()
  type(tem_st_fun_linkedList_type) :: st_funList
  type(tem_st_fun_listElem_type), pointer :: st_fun
  type(solver_type), target :: solver
  type(tracking_type) :: tracking
  type(tem_variable_type), allocatable :: newVar(:)
  integer :: addedPos
  integer :: iElem, iComp, iDof, elem_offset, comp_offset
  integer :: nComp
  integer :: iVar, iList, iSt
  logical :: wasAdded
  real(kind=rk), allocatable :: res(:)
  integer, allocatable :: error(:)
  character(len=labelLen), allocatable :: input_varname(:)
  !> dummy variable required by creat_subTree
  type( tem_stencilHeader_type ) :: stencil
  type(tem_st_fun_listElem_type), pointer :: newElem

  write(*,*) 'Hello from tem_varSys_test'
  ! load utest mesh
  call load_env( tree     = solver%tree,     &
    &            boundary = solver%boundary, &
    &            general  = solver%general   )

  write(*,*) 'nElems ', solver%tree%nElems

  allocate(solver%general%solver%conf(1))
  call load_config( conf             = solver%general%solver%conf(1), &
    &               chunk            = trim(sysConf),                 &
    &               stfun_linkedList = st_funList,                    &
    &               tracking         = tracking,                      &
    &               newVar           = newVar                         )

  ! initialize variable system
  write(*,*) 'calling varsys init'
  call tem_varSys_init(me = varSys, systemName = 'utest')

  ! add state variable
  write(*,*) 'add state variable '
  get_element => access_state
  call tem_varSys_append_stateVar( me             = varSys,         &
    &                              varName        = 'state',        &
    &                              nComponents    = 4,              &
    &                              method_data    = c_loc(solver),  &
    &                              get_point      = get_point,      &
    &                              get_element    = get_element,    &
    &                              set_params     = set_params,     &
    &                              get_params     = get_params,     &
    &                              setup_indices  = setup_indices,  &
    &                              get_valOfIndex = get_valOfIndex, &
    &                              pos            = addedPos,       &
    &                              wasAdded       = wasAdded        )

  ! assign values to state
  solver%nDofs = 1
  allocate(solver%state(solver%tree%nElems*varSys%nScalars*solver%nDofs))
  do iVar = 1, varSys%nStateVars
    write(*,*) 'initializing state for ', trim(varSys%varname%val(iVar))
    nComp = varSys%method%val(iVar)%nComponents
    do iElem = 1, solver%tree%nElems
      write(*,*) 'iElem ', iElem
      elem_offset =  (iElem-1)*varSys%nScalars*solver%nDofs
      do iDof = 1, solver%nDofs
        comp_offset = elem_offset + (iDof-1)*varSys%nScalars
        do iComp = 1, varSys%method%val(iVar)%nComponents
          solver%state( comp_offset + varSys%method%val(iVar)%state_varPos(iComp) )  &
            &   = real( comp_offset + varSys%method%val(iVar)%state_varPos(iComp),   &
            &           kind=rk)
        end do
        write(*,*) 'iDof ', iDof , ' state ', &
          & solver%state(comp_offset + varSys%method%val(iVar)%state_varPos(1) : &
          &              comp_offset + varSys%method%val(iVar)%state_varPos(nComp))
      end do
    end do
  end do


  write(*,*)
  write(*,*) 'add density variable'
  ! add density variable
  get_element => derive_density
  allocate(input_varname(1))
  input_varname(1) = 'state'
  call tem_varSys_append_derVar( me             = varSys,         &
    &                            varName        = 'density',      &
    &                            operType       = 'state',        &
    &                            nComponents    = 1,              &
    &                            input_varname  = input_varname,  &
    &                            method_data    = c_loc(solver),  &
    &                            get_point      = get_point,      &
    &                            get_element    = get_element,    &
    &                            set_params     = set_params,     &
    &                            get_params     = get_params,     &
    &                            setup_indices  = setup_indices,  &
    &                            get_valOfIndex = get_valOfIndex, &
    &                            pos            = addedPos,       &
    &                            wasAdded       = wasAdded        )


  write(*,*)
  write(*,*) 'add new variables. nVars:', size(newVar)
  call init(tracking%varPos)
  do iVar = 1, size(newVar)
    write(*,*) 'appending variable, ',iVar,':'//trim(newVar(iVar)%label)
    select case(trim(newVar(iVar)%varType))
    case('st_fun')
      call append( st_funList, newVar(iVar)%st_fun, newElem )

      if ( newVar(iVar)%nComponents == 1 ) then
        get_element => evaluate_add_spacetime_scalarByTreeID
      else
        get_element => evaluate_add_spacetime_vectorByTreeID
      end if

      !do ifun = 1, size(newVar(iVar)%st_fun)
      call tem_varSys_append_derVar(me             = varSys,                   &
        &                           varName        = newVar(iVar)%label,       &
        &                           operType       = 'st_fun',                 &
        &                           nComponents    = newVar(iVar)%nComponents, &
        &                           method_data    = c_loc(newElem),&
        &                           get_point      = get_point,                &
        &                           get_element    = get_element,              &
        &                           set_params     = set_params,               &
        &                           get_params     = get_params,               &
        &                           setup_indices  = setup_indices,            &
        &                           get_valOfIndex = get_valOfIndex,           &
        &                           pos            = addedPos,                 &
        &                           wasAdded       = wasAdded                  )
    !case('operation')
    !  !\todo create subroutines for different operations
    case default
      write(*,*) 'varType not supported. Variable is not appended.'
      continue
    end select
  end do

  write(*,*) 'create subtree for all space time functions stored '
  write(*,*) 'in linked list of spacetime function'

  st_fun => st_funList%head
  iList = 0
  do
    if (.not. associated(st_fun)) EXIT
    iList = iList + 1
    !write(*,*) 'iList ', iList, ' st_fun nVals ', st_fun%nVals
    do iSt = 1, st_fun%nVals
      !write(*,*) 'iSt ', iSt, ' stfun kind ', trim(st_fun%val(iSt)%fun_kind)
      call tem_create_subTree_of( inTree = solver%tree, &
        &                         bc_prop = solver%boundary, &
        &                         stencil = stencil, &
        &                         subTree = st_fun%val(iSt)%subTree, &
        &                         inShape = st_fun%val(iSt)%geom )
    end do
    st_fun => st_fun%next
  end do
  write(*,*) 'Done creating subtree for all spacetime functions'

  write(*,*)

  ! create tracking variable poisition in the global varSys
  call init(tracking%varPos)
  do iVar = 1, tracking%nRequestedVars
    addedPos = PositionOfVal( me  = varSys%varname, &
      &                       val = trim(tracking%variable(iVar)) )
    if (addedPos > 0) then
      call append( me = tracking%varPos, val = addedPos )
    else
      write(*,*) 'Variable ', trim(tracking%variable(iVar)), ' not found in varSys'
    end if
  end do
  write(*,*)

  ! track variables
  allocate(error(tracking%varPos%nVals))
  write(*,*) 'Checking variable extraction '
  do iVar=1,tracking%varPos%nVals
    if (allocated(res))  deallocate(res)
    addedPos = tracking%varPos%val(iVar)
    write(*,*) 'track var ', trim(varSys%varname%val(addedPos))
    allocate(res(nElems_track*varSys%method%val(addedPos)%nComponents &
      &          *solver%nDofs))
    ! access state array for given elemPos
    call varSys%method%val(addedPos)%get_element(                       &
      &                                varSys  = varSys,                &
      &                                elemPos = elemPos,               &
      &                                time    = solver%general         &
      &                                                %simControl%now, &
      &                                tree    = solver%tree,           &
      &                                nElems  = nElems_track,          &
      &                                nDofs   = solver%nDofs,          &
      &                                res     = res                    )

    call check_res( varname = varSys%varname%val(addedPos) , &
      &             res     = res,                           &
      &             error   = error(iVar)                    )

    if (error(iVar) == -1) exit

!    do iElem = 1, nElems_track
!      write(*,*) 'iElem ', iElem, 'elemPos ', elemPos(iElem)
!      write(*,*) 'res ', &
!        & res( (iElem-1)*varSys%method%val(addedPos)%nComponents + 1 : &
!        &      iElem*varSys%method%val(addedPos)%nComponents )
!    end do
  end do

  if (any(error/=0)) then
    write(*,*) 'FAILED'
  else
    write(*,*) 'PASSED'
  end if  
  call close_config(L=solver%general%solver%conf(1))
  call fin_env()


contains


  ! ****************************************************************************!
  !> load variables defined in sysConf
  subroutine load_config( conf, chunk, stfun_linkedList, tracking, newVar )
    ! --------------------------------------------------------------------------!
    type(flu_state) :: conf
    character(len=*), intent(in) :: chunk
    type(tem_st_fun_linkedList_type), intent(out) :: stfun_linkedList
    type(tracking_type), intent(out) :: tracking
    type(tem_variable_type), allocatable, intent(out) :: newVar(:)
    ! --------------------------------------------------------------------------!
    integer :: nVars, iVar, thandle, iError
    integer, allocatable :: vError(:)
    ! --------------------------------------------------------------------------!
    call open_config_chunk(L=conf, chunk=trim(chunk))

    call tem_variable_load( me   = newVar,      &
      &                     conf = conf,        &
      &                     key  = 'variables', &
      &                     vError  = vError    )

    call aot_table_open( L       = conf,            &
      &                  thandle = thandle,         &
      &                  key     = 'track_variable' )
    nVars = aot_table_length( L = conf, thandle = thandle )
    tracking%nRequestedVars = nVars
    allocate(tracking%variable(nVars))
    write(*,*) 'track_variables '
    do iVar = 1, nVars
      ! Get the names of the variable
      call aot_get_val( L       = conf,                    &
        &               thandle = thandle,                 &
        &               val     = tracking%variable(iVar), &
        &               ErrCode = iError,                  &
        &               pos     = iVar                     )
    write(*,*) trim(tracking%variable(iVar))
    end do
    call aot_table_close( L = conf, thandle = thandle)

  end subroutine load_config
  ! ****************************************************************************!


  ! ****************************************************************************!
  !> access state variables
  subroutine access_state(fun, varsys, elempos, time, tree, n, nDofs, res)
    ! --------------------------------------------------------------------------!
    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> TreeID of the element to get the variable for.
    integer, intent(in) :: elempos(:)

    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in)  :: time

    !> global treelm mesh info
    type(treelmesh_type), intent(in) :: tree

    !> Number of values to obtain for this variable (vectorized access).
    integer, intent(in) :: n

    !> Number of degrees of freedom within an element.
    integer, intent(in) :: nDofs

    !> Resulting values for the requested variable.
    !!
    !! The first dimension are the n requested entries and the second
    !! dimension are the components of this variable.
    !! Third dimension are the degrees of freedom.
    real(kind=rk), intent(out) :: res(:)
    ! --------------------------------------------------------------------------!
    integer :: iElem, iComp, iDof, pos
    type(solver_type), pointer :: fPtr
    ! --------------------------------------------------------------------------!
    write(*,*) 'access variable label: ', trim(varSys%varname%val(fun%myPos))
    write(*,*) 'at time: ', time%sim
    write(*,*) 'tree%nElems: ', tree%nElems
    call C_F_POINTER( fun%method_Data, fPtr )

    res = 0.0_rk
    do iElem = 1, n
      ! if state array is defined level wise then use levelPointer(pos)
      ! to access state array
      pos = elemPos(iElem)
      do iDof = 1, nDofs
        do iComp = 1, fun%nComponents
          res( (iElem-1)*fun%nComponents*nDofs           &
            &  + (iDof-1)*fun%nComponents                &
            &  + iComp ) =                               &
            &  fPtr%state( (pos-1)*varSys%nScalars*nDofs &
            &  + (iDof-1)*varSys%nScalars                &
            ! position of this state variable in the state array
            &  + fun%state_varPos(iComp) )

        end do
      end do
    end do

  end subroutine access_state
  ! ****************************************************************************!


  ! ****************************************************************************!
  !> derive density variables
  subroutine derive_density(fun, varsys, elempos, time, tree, n, nDofs, res)
    ! --------------------------------------------------------------------------!
    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> TreeID of the element to get the variable for.
    integer, intent(in) :: elempos(:)

    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in)  :: time

    !> global treelm mesh info
    type(treelmesh_type), intent(in) :: tree

    !> Number of values to obtain for this variable (vectorized access).
    integer, intent(in) :: n

    !> Number of degrees of freedom within an element.
    integer, intent(in) :: nDofs

    !> Resulting values for the requested variable.
    !!
    !! The first dimension are the n requested entries and the second
    !! dimension are the components of this variable.
    !! Third dimension are the degrees of freedom.
    real(kind=rk), intent(out) :: res(:)
    ! --------------------------------------------------------------------------!
    integer :: iElem, iComp, iDof, pos, input_varPos, inVar_nComp
    real(kind=rk) :: dens
    type(solver_type), pointer :: fPtr
    ! --------------------------------------------------------------------------!
    write(*,*) 'myPos ', fun%myPos
    write(*,*) 'derive variable label: ', trim(varSys%varname%val(fun%myPos))
    write(*,*) 'at time: ', time%sim
    write(*,*) 'tree%nElems: ', tree%nElems
    call C_F_POINTER( fun%method_Data, fPtr )

    input_varPos = fun%input_varPos(1)
    inVar_nComp = varSys%method%val(input_varPos)%nComponents
    write(*,*) 'input_varPos ', input_varPos
    res = 0.0_rk

    do iElem = 1, n
      ! if state array is defined level wise then use levelPointer(pos)
      ! to access state array
      pos = elemPos(iElem)
      ! use state_varPos(iComp) if state has more than one variable
      do iDof = 1, nDofs
        dens = 0.0_rk
        do iComp = 1, inVar_nComp
          dens = dens + fPtr%state( (pos-1)*varSys%nScalars*nDofs            &
            &             + (iDof-1)*varSys%nScalars                         &
            &             + varSys%method%val(input_varPos)%state_varPos(iComp) )
        end do
        res( (iElem-1)*nDofs + iDof )  = dens
      end do
    end do

  end subroutine derive_density
  ! ****************************************************************************!

  ! ***************************************************************************!
  !> Check tracking results
  subroutine check_res( varname, res, error )
    ! ---------------------------------------------------------------------- !
    character(len=labelLen), intent(in) :: varname
    real(kind=rk), intent(in) :: res(:)
    integer, intent(out) :: error
    ! ---------------------------------------------------------------------- !
    integer :: iElem  
    real(kind=rk) :: diff_dens(nElems_track)
    real(kind=rk) :: diff_state(nElems_track,4)
    real(kind=rk) :: diff_vel(nElems_track,3)
    ! ---------------------------------------------------------------------- !
    error = 0 
    select case (trim(varname))
    case ('state')
      do iElem = 1, nElems_track
        diff_state(iElem,:) = res((iElem-1)*4 + 1 : (iElem-1)*4 + 4 ) &
          &                 - state(iElem,:)
      end do    
      if ( any(diff_state > eps) ) then
        error = -1
        write(*,*) 'Refernece value for state does not match'
      end if  
      if (dumpRes) then
        do iElem = 1, nElems_track
          write(*,*) 'reference: ', state(iElem, :)
          write(*,*) '   output: ', res( (iElem-1)*4 + 1 : (iElem-1)*4 + 4 )
        end do
      end if  
    case ('density')
      diff_dens = res(:) - dens
      if ( any(diff_dens > eps) ) then
        error = -1
        write(*,*) 'Refernece value for density does not match'
      end if  
      if (dumpRes) then
        write(*,*) 'reference: ', dens
        write(*,*) '   output: ', res
      end if  
    case ('velocity')
      do iElem = 1, nElems_track
        diff_vel(iElem,:) = res((iElem-1)*3 + 1 : (iElem-1)*3 + 3 ) &
          &               - vel(iElem,:)
      end do
      if ( any(diff_vel > eps) ) then
        error = -1
        write(*,*) 'Refernece value for velocity does not match'
      end if  
      if (dumpRes) then
        do iElem = 1, nElems_track
          write(*,*) 'reference: ', vel(iElem, :)
          write(*,*) '   output: ', res( (iElem-1)*3 + 1 : (iElem-1)*3 + 3 )
        end do  
      end if  
    case default
      write(*,*) 'Unknown variable name', trim(varname)
      error = -1
    end select
      
   end subroutine check_res
  ! ***************************************************************************!


end program tem_varSys_test


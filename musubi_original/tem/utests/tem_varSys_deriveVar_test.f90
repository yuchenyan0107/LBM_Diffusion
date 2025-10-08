! Copyright (c) 2014, 2016 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2014 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2016-2017 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
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
!> This UTEST is to test variable system.
!! load variables and access state and derive quantities from 
!! c pointer methoddata
program tem_varSys_deriveVar_test
  use, intrinsic :: iso_c_binding, only: C_NEW_LINE, c_ptr, c_loc, c_f_pointer

  use env_module,               only: rk, labelLen, fin_env
  use tem_bc_prop_module,       only: tem_bc_prop_type
  use tem_general_module,       only: tem_general_type
  use treelmesh_module,         only: treelmesh_type
  use tem_time_module,          only: tem_time_type
  use tem_utestEnv_module,      only: load_env
  use tem_tools_module,         only: upper_to_lower
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
  use tem_grow_array_module,    only: grw_intArray_type, init, append, destroy
  use tem_dyn_array_module,     only: PositionOfVal
  use tem_spacetime_fun_module, only: tem_spacetime_fun_type,                  &
    &                                 tem_load_spacetime,                      &
    &                                 tem_st_fun_listElem_type,                &
    &                                 tem_st_fun_linkedList_type, append,      &
    &                                 tem_spacetime_for

  use aotus_module,          only: open_config_chunk, close_config, flu_state
  use aot_table_module,      only: aot_table_open, aot_table_close,            &
    &                              aot_table_length, aot_get_val

  !mpi!nprocs = 1

  implicit none

  type solver_type
    integer :: nComps
    integer :: nDofs
    real(kind=rk), allocatable :: state(:)
    type(treelmesh_type) :: tree
    type(tem_bc_prop_type) :: boundary
    type(tem_general_type) :: general
  end type solver_type

  type(tem_varSys_type) :: varSys
  procedure(tem_varSys_proc_point), pointer :: get_point => NULL()
  procedure(tem_varSys_proc_element), pointer :: get_element => NULL()
  procedure(tem_varSys_proc_setParams), pointer :: set_params => null()
  procedure(tem_varSys_proc_getParams), pointer :: get_params => null()
  procedure(tem_varSys_proc_setupIndices), pointer :: &
    &                                      setup_indices => null()
  procedure(tem_varSys_proc_getValOfIndex), pointer :: &
    &                                       get_valOfIndex => null()
  type(solver_type), target :: solver
  integer :: addedPos, statePos, densityPos
  integer :: iElem, iComp, iDof, elem_offset, comp_offset, nElems_track
  ! position in global treeID list
  integer, allocatable :: elemPos(:)
  logical :: wasAdded
  real(kind=rk), allocatable :: res(:)
  character(len=labelLen), allocatable :: input_varname(:)

  write(*,*) 'Hello from tem_varSys_test'
  ! load utest mesh
  call load_env( tree     = solver%tree,     &
    &            boundary = solver%boundary, &
    &            general  = solver%general   ) 

  write(*,*) 'nElems ', solver%tree%nElems

  solver%nDofs = 1
  solver%nComps = 4
  allocate(solver%state(solver%tree%nElems*solver%nComps*solver%nDofs))
  ! assign values to state
  do iElem = 1, solver%tree%nElems
    write(*,*) 'iElem ', iElem
    elem_offset =  (iElem-1)*solver%nComps 
    do iComp = 1, solver%nComps
      comp_offset = elem_offset + (iComp-1)*solver%nDofs
      do iDof = 1, solver%nDofs
        solver%state( comp_offset + iDof ) = real(comp_offset+iDof, kind=rk)
      end do
      write(*,*) 'iComp ', iComp , ' state ', &
        &        solver%state(comp_offset + 1 : comp_offset + solver%nDofs)
    end do
  end do   


  ! initialize variable system
  write(*,*) 'calling varsys init'
  call tem_varSys_init(me = varSys, systemName = 'utest')

  ! add state variable
  write(*,*) 'add state variable '
  get_element => access_state
  call tem_varSys_append_stateVar( me             = varSys,         &
    &                              varName        = 'state',        &
    &                              nComponents    = solver%nComps,  &
    &                              method_data    = c_loc(solver),  &
    &                              get_point      = get_point,      &
    &                              get_element    = get_element,    &
    &                              set_params     = set_params,     &
    &                              get_params     = get_params,     &
    &                              setup_indices  = setup_indices,  &
    &                              get_valOfIndex = get_valOfIndex, &
    &                              pos            = addedPos,       &
    &                              wasAdded       = wasAdded        )

!KM!  write(*,*) 'addedPos ', addedPos
!KM!  write(*,*) 'added variable ', trim(varSys%varname%val(addedPos))
!KM!  write(*,*) 'state_varPos ', varSys%method%val(addedPos)%state_varPos

  statePos = PositionOfVal( varSys%varname, 'state' )
!KM!  write(*,*) 'state_pos ', statePos

  nElems_track = 5
  allocate(elemPos(nElems_track))
  elemPos = (/ 1, 3, 5, 7, 8 /)

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

  write(*,*) 'addedPos ', addedPos
  write(*,*) 'added variable ', trim(varSys%varname%val(addedPos))

  densityPos = PositionOfVal( varSys%varname, 'density' )
  write(*,*) 'density_pos ', densityPos, &
    & 'input_varPos ', varSys%method%val(densityPos)%input_varPos

  allocate(res(nElems_track*varSys%method%val(densityPos)%nComponents*solver%nDofs))
  ! access state array for given elemPos
  call varSys%method%val(densityPos)%get_element(                             &
    &                                varSys  = varSys,                        &
    &                                elemPos = elemPos,                       &
    &                                time    = solver%general%simControl%now, &
    &                                tree    = solver%tree,                   &
    &                                nElems  = nElems_track,                  &
    &                                nDofs   = solver%nDofs,                  &
    &                                res     = res                            )

  do iElem = 1, nElems_track
    write(*,*) 'iElem ', iElem, 'elemPos ', elemPos(iElem)
    ! nComp = 1
    write(*,*) 'density ', &
      & res((iElem-1)*solver%nDofs+1: iElem*solver%nDofs)
  end do  

  write(*,*) 'PASSED'
  call fin_env()


contains


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
    integer :: iElem, iComp, iDof, pos, nScalars
    type(solver_type), pointer :: fPtr
    ! --------------------------------------------------------------------------!
    write(*,*) 'access variable label: ', trim(varSys%varname%val(fun%myPos))
    write(*,*) 'at time: ', time%sim
    write(*,*) 'tree%nElems: ', tree%nElems
    call C_F_POINTER( fun%method_Data, fPtr )

    ! number of scalars in state array
    nScalars = varSys%nScalars

    res = 0.0_rk
    do iElem = 1, n
      ! if state array is defined level wise then use levelPointer(pos)
      ! to access state array
      pos = elemPos(iElem)
      do iDof = 1, nDofs
        do iComp = 1, fun%nComponents
          res( (iElem-1)*fun%nComponents*nDofs    &
            &  + (iDof-1)*fun%nComponents         &
            &  + iComp ) =                        &
            &  fPtr%state( (pos-1)*nScalars*nDofs &
            &  + (iDof-1)*nScalars                &
            ! position of this variable in the state array
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
    ! -------------------------------------------------------------------------!
    integer :: iElem, iComp, iDof, pos, input_varPos, inVar_nComp
    type(solver_type), pointer :: fPtr
    real(kind=rk) :: dens
    ! -------------------------------------------------------------------------!
    write(*,*) 'derive variable label: ', trim(varSys%varname%val(fun%myPos))
    write(*,*) 'at time: ', time%sim
    write(*,*) 'tree%nElems: ', tree%nElems
    call C_F_POINTER( fun%method_Data, fPtr )

    input_varPos = fun%input_varPos(1)
    inVar_nComp = varSys%method%val(input_varPos)%nComponents
    write(*,*) 'input_varPos ', input_varPos, 'invar_nComp', inVar_nComp
    res = 0.0_rk
    do iElem = 1, n
      ! if state array is defined level wise then use levelPointer(pos)
      ! to access state array
      pos = elemPos(iElem)
      ! use state_varPos(iComp) if state has more than one variable
      do iDof = 1, nDofs
        dens = 0.0_rk
        do iComp = 1, inVar_nComp
          dens = dens + fPtr%state( (pos-1)*varSys%nScalars*nDofs     &
            &                       + (iDof-1)*varSys%nScalars        &
            &                       + varSys%method%val(input_varPos) &
            &                               %state_varPos(iComp)      )
        end do
        res( (iElem-1)*nDofs + iDof )  = dens
      end do
    end do

  end subroutine derive_density
  ! ****************************************************************************!

end program tem_varSys_deriveVar_test


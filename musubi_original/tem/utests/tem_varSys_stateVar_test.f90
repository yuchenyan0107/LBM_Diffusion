! Copyright (c) 2014, 2016 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2018 Robin Weihe <robin.weihe@student.uni-siegen.de>
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
program tem_varSys_stateVar_test
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
    &                                 tem_varSys_proc_point,                   &
    &                                 tem_varSys_proc_element,                 &
    &                                 tem_varSys_proc_setParams,               &
    &                                 tem_varSys_proc_getParams,               &
    &                                 tem_varSys_proc_setupIndices,            &
    &                                 tem_varSys_proc_getValOfIndex
  use tem_dyn_array_module,     only: PositionOfVal

  !mpi!nprocs = 1

  implicit none

  type scheme_type
    character(len=labelLen) :: label
    real(kind=rk), allocatable :: state(:)
    integer :: nDofs
    type(tem_varSys_type) :: varSys
  end type scheme_type


  type solver_type
    type(scheme_type), pointer :: scheme => null()
    type(treelmesh_type), pointer :: tree => null()
    type(tem_bc_prop_type), pointer :: boundary => null()
    type(tem_general_type), pointer :: general => null()
  end type solver_type

  procedure(tem_varSys_proc_point), pointer :: get_point => NULL()
  procedure(tem_varSys_proc_element), pointer :: get_element => NULL()
  procedure(tem_varSys_proc_setParams), pointer :: set_params => null()
  procedure(tem_varSys_proc_getParams), pointer :: get_params => null()
  procedure(tem_varSys_proc_setupIndices), pointer :: &
    &                                      setup_indices => null()
  procedure(tem_varSys_proc_getValOfIndex), pointer :: &
    &                                       get_valOfIndex => null()
  integer :: addedPos, statePos_1, statePos_2
  integer :: iElem, iComp, iDof, iVar, elem_offset, comp_offset, nElems_track
  integer :: nComp, state_varPos
  ! position in global treeID list
  integer, allocatable :: elemPos(:)
  logical :: wasAdded
  real(kind=rk), allocatable :: res(:)
  type(treelmesh_type), target :: tree
  type(tem_bc_prop_type), target:: boundary
  type(tem_general_type), target :: general
  type(scheme_type), allocatable, target :: scheme(:)
  type(solver_type), allocatable, target :: solver(:)
  integer :: iScheme, nScheme

  write(*,*) 'Hello from tem_varSys_test'
  nScheme = 2
  allocate(solver(nScheme))
  do iScheme = 1, nScheme
    solver(iScheme)%tree => tree
    solver(iScheme)%boundary => boundary
    solver(iScheme)%general => general
  end do

  ! load utest mesh
  call load_env( tree     = tree,     &
    &            boundary = boundary, &
    &            general  = general   )

  write(*,*) 'nElems ', solver(1)%tree%nElems
  allocate(scheme(2))
  solver(1)%scheme => scheme(1)
  solver(2)%scheme => scheme(2)
  scheme(1)%label = 'scheme_1'
  scheme(2)%label = 'scheme_2'

  scheme(1)%nDofs = 1
  scheme(2)%nDofs = 4

  do iScheme=1, nScheme
    ! initialize variable system
    write(*,*) 'calling varsys init for ', iScheme, trim(scheme(iScheme)%label)
    write(*,*) ' nDofs ', solver(iScheme)%scheme%nDofs
    call tem_varSys_init(me         = scheme(iScheme)%varSys, &
      &                  systemName = 'utest'//trim(scheme(iScheme)%label) )

    ! add state variable
    write(*,*) 'add state variable '
    get_element => access_state
    call tem_varSys_append_stateVar( me             = scheme(iScheme)%varSys, &
      &                              varName        = 'state1',               &
      &                              nComponents    = 1,                      &
      &                              method_data    = c_loc(solver(iScheme)), &
      &                              get_point      = get_point,              &
      &                              get_element    = get_element,            &
      &                              set_params     = set_params,             &
      &                              get_params     = get_params,             &
      &                              setup_indices  = setup_indices,          &
      &                              get_valOfIndex = get_valOfIndex,         &
      &                              pos            = addedPos,               &
      &                              wasAdded       = wasAdded                )

    call tem_varSys_append_stateVar( me             = scheme(iScheme)%varSys, &
      &                              varName        = 'state2',               &
      &                              nComponents    = 3,                      &
      &                              method_data    = c_loc(solver(iScheme)), &
      &                              get_point      = get_point,              &
      &                              get_element    = get_element,            &
      &                              set_params     = set_params,             &
      &                              get_params     = get_params,             &
      &                              setup_indices  = setup_indices,          &
      &                              get_valOfIndex = get_valOfIndex,         &
      &                              pos            = addedPos,               &
      &                              wasAdded       = wasAdded                )

!KM!  write(*,*) 'addedPos ', addedPos
!KM!  write(*,*) 'added variable ', trim(varSys%varname%val(addedPos))
!KM!  write(*,*) 'state_varPos ', varSys%method%val(addedPos)%state_varPos

    statePos_1 = PositionOfVal( scheme(iScheme)%varSys%varname, 'state1' )
    statePos_2 = PositionOfVal( scheme(iScheme)%varSys%varname, 'state2' )

    allocate( scheme(iScheme)%state( tree%nElems                     &
      &                            * scheme(iScheme)%varSys%nScalars &
      &                            * scheme(iScheme)%nDofs)          )
    ! assign values to state
    do iVar = 1, scheme(iScheme)%varSys%nStateVars
      write(*,*) 'initializing state for ', &
        & trim(scheme(iScheme)%varSys%varname%val(iVar))
      nComp = scheme(iScheme)%varSys%method%val(iVar)%nComponents
      do iElem = 1, tree%nElems
        write(*,*) 'iElem ', iElem
        elem_offset =  (iElem-1)*scheme(iScheme)%varSys%nScalars &
          &         * scheme(iScheme)%nDofs
        do iDof = 1, scheme(iScheme)%nDofs
          comp_offset = elem_offset + (iDof-1)*scheme(iScheme)%varSys%nScalars
          do iComp = 1, scheme(iScheme)%varSys%method%val(iVar)%nComponents
            state_varPos = scheme(iScheme)%varSys%method%val(iVar) &
              &                           %state_varPos(iComp)
            scheme(iScheme)%state( comp_offset + state_varPos ) &
              &   = real( ( (iScheme-1)*nScheme                 &
              &   + comp_offset + state_varPos ), kind=rk )
          end do
          write(*,*) 'iDof ', iDof , ' state ', &
            & scheme(iScheme)%state(comp_offset &
            & + scheme(iScheme)%varSys%method%val(iVar)%state_varPos(1) :  &
            &                       comp_offset &
            & + scheme(iScheme)%varSys%method%val(iVar)%state_varPos(nComp))
        end do !iDof
      end do !iElem
    end do !iVar

  end do !iScheme


!KM!  write(*,*) 'state_pos ', statePos

  nElems_track = 5
  allocate(elemPos(nElems_track))
  elemPos = (/ 1, 3, 5, 7, 8 /)
  !elemPos = (/ 1, 2, 3, 4, 5, 6, 7, 8 /)

  allocate(res(nElems_track &
    & * maxval(scheme(1)%varSys%method%val(:)%nComponents) &
    & * maxval(scheme(:)%nDofs)))

  do iScheme=1, nScheme
    write(*,*) 'Scheme ', trim(scheme(iScheme)%label)
    do iVar = 1, scheme(iScheme)%varSys%nStateVars
      write(*,*) 'Get element for ', trim(scheme(iScheme)%varSys%varname  &
        &                                                       %val(iVar))
      nComp = scheme(iScheme)%varSys%method%val(iVar)%nComponents

      ! access state array for given elemPos
      call scheme(iScheme)%varSys%method%val(iVar)%get_element(            &
        &                                varSys  = scheme(iScheme)%varSys, &
        &                                elemPos = elemPos,                &
        &                                time    = general%simControl%now, &
        &                                tree    = tree,                   &
        &                                nElems  = nElems_track,           &
        &                                nDofs   = scheme(iScheme)%nDofs,  &
        &                                res     = res                     )


      do iElem = 1, nElems_track
        write(*,*) 'iElem ', iElem, 'elemPos ', elemPos(iElem)
        elem_offset =  (iElem-1)*nComp*scheme(iScheme)%nDofs
        do iDof = 1, scheme(iScheme)%nDofs
          comp_offset = elem_offset + (iDof-1)*nComp
          write(*,*) 'iDof ', iDof , ' res ',&
            & res( comp_offset+1: &
            &      comp_offset+scheme(iScheme)%varSys%method%val(iVar) &
            &                                        %nComponents )
        end do
      end do
    end do
  end do

  write(*,*) 'PASSED'
  call fin_env()


contains


  ! ***************************************************************************!
  !> access state variables
  subroutine access_state(fun, varsys, elempos, time, tree, nElems, nDofs, res)
    ! -------------------------------------------------------------------------!
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
    integer, intent(in) :: nElems

    !> Number of degrees of freedom within an element.
    integer, intent(in) :: nDofs

    !> Resulting values for the requested variable.
    !!
    !! The first dimension are the n requested entries and the second
    !! dimension are the components of this variable.
    !! Third dimension are the degrees of freedom.
    real(kind=rk), intent(out) :: res(:)
    ! -------------------------------------------------------------------------!
    integer :: iElem, iComp, iDof, pos, nScalars
    type(solver_type), pointer :: fPtr
    ! -------------------------------------------------------------------------!
    write(*,*) 'access variable label: ', trim(varSys%varname%val(fun%myPos))
    write(*,*) 'at time: ', time%sim
    write(*,*) 'tree%nElems: ', tree%nElems
    call C_F_POINTER( fun%method_Data, fPtr )
    write(*,*) 'scheme label ', trim(fPtr%scheme%label)

    ! number of scalars in state array
    nScalars = varSys%nScalars

    res = 0.0_rk
    do iElem = 1, nElems
      ! if state array is defined level wise then use levelPointer(pos)
      ! to access state array
      pos = elemPos(iElem)
      do iDof = 1, nDofs
        do iComp = 1, fun%nComponents
          res( (iElem-1)*fun%nComponents*nDofs           &
            &  + (iDof-1)*fun%nComponents                &
            &  + iComp ) =                               &
            &  fPtr%scheme%state( (pos-1)*nScalars*nDofs &
            &  + (iDof-1)*nScalars                       &
            ! position of this state variable in the state array
            &  + fun%state_varPos(iComp) )
        end do !iComp
      end do !iDof
    end do !iElem

  end subroutine access_state
  ! ***************************************************************************!

end program tem_varSys_stateVar_test

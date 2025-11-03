! Copyright (c) 2014, 2016 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2014-2015 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2015, 2017-2018 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2016, 2019 Harald Klimach <harald.klimach@uni-siegen.de>
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
program tem_varSys_stfunVar_test
  use, intrinsic :: iso_c_binding, only: C_NEW_LINE, c_ptr, c_loc, c_f_pointer

  use env_module,               only: rk, eps, solSpecLen, labelLen, fin_env
  use tem_bc_prop_module,       only: tem_bc_prop_type
  use tem_general_module,       only: tem_general_type
  use treelmesh_module,         only: treelmesh_type
  use tem_time_module,          only: tem_time_type
  use tem_utestEnv_module,      only: load_env
  use tem_tools_module,         only: upper_to_lower, tem_PositionInSorted
  use tem_varSys_module,        only: tem_varSys_init, tem_varSys_type
  use tem_dyn_array_module,     only: PositionOfVal
  use tem_variable_module,      only: tem_variable_type, tem_variable_load
  use tem_grow_array_module,    only: grw_intArray_type, init, append, destroy
  use tem_dyn_array_module,     only: PositionOfVal
  use tem_spacetime_fun_module, only: tem_st_fun_linkedList_type,              &
    &                                 tem_create_subTree_of_st_funList
  use tem_derived_module,       only: tem_varSys_append_luaVar
  use tem_logging_module,       only: tem_logging_init
  use tem_subTree_module,       only: tem_create_subTree_of
  use tem_geometry_module,      only: tem_BaryOfId

  use aotus_module,          only: open_config_chunk, close_config, flu_state
  use aot_table_module,      only: aot_table_open, aot_table_close,            &
    &                              aot_table_length, aot_get_val

  !mpi!nprocs = 1

  implicit none

  character, parameter :: nl = C_NEW_LINE

  character(len=solSpecLen), parameter :: sysConf = &
    &    'function vel_fun(x,y,z,t)' // nl          &
    & // '  return {x,y,z}' // nl                  &
    & // 'end' // nl                                &
    & // 'function dens_fun(x,y,z,t)' // nl         &
    & // '  return x' // nl                         &
    & // 'end' // nl                                &
    & // 'shape = {' // nl                          &
    & // '  kind = "canoND",' // nl                 &
    & // '  object = {' // nl                       &
    & // '    origin = {0.0,0.0,0.0},' //nl         &
    & // '    vec = {0.0,1.0,0.0},' // nl           &
    & // '    segments = 2' // nl                   &
    & // '  }' // nl                                &
    & // '}' // nl                                  &
    & // 'variable = {' // nl                       &
    & // '  {' // nl                                &
    & // '    name = "density",' // nl              &
    & // '    ncomponents = 1,' // nl               &
    & // '    vartype = "st_fun",' // nl            &
    & // '    st_fun = dens_fun' // nl              &
    & // '  },' // nl                               &
    & // '  {' // nl                                &
    & // '    name = "density_shape",' // nl        &
    & // '    ncomponents = 1,' // nl               &
    & // '    vartype = "st_fun",' // nl            &
    & // '    st_fun = {' // nl                     &
    & // '      fun = dens_fun,' // nl              &
    & // '      shape = shape' // nl                &
    & // '    }' // nl                              &
    & // '  },' // nl                               &
    & // '  {' // nl                                &
    & // '    name = "density_multi_stfun",' // nl  &
    & // '    ncomponents = 1,' // nl               &
    & // '    vartype = "st_fun",' // nl            &
    & // '    st_fun = {' // nl                     &
    & // '      {' // nl                            &
    & // '        fun = dens_fun,' // nl            &
    & // '        shape = shape' // nl              &
    & // '      },' // nl                           &
    & // '      10.0' // nl                         &
    & // '    }' // nl                              &
    & // '  },' // nl                               &
    & // '  {' // nl                                &
    & // '    name = "velocity",' // nl             &
    & // '    ncomponents = 3,' // nl               &
    & // '    vartype = "st_fun",' // nl            &
    & // '    st_fun = vel_fun' // nl               &
    & // '  },' // nl                               &
    & // '  {' // nl                                &
    & // '    name = "velocity_shape",' // nl       &
    & // '    ncomponents = 3,' // nl               &
    & // '    vartype = "st_fun",'// nl             &
    & // '    st_fun = {' // nl                     &
    & // '      const={1.0,2.0,3.0},' // nl         &
    & // '      shape=shape' // nl                  &
    & // '    }' // nl                              &
    & // '  },' // nl                               &
    & // '  {' // nl                                &
    & // '    name = "velocity_multi_stfun",' // nl &
    & // '    ncomponents = 3,' // nl               &
    & // '    vartype = "st_fun",'// nl             &
    & // '    st_fun = {' // nl                     &
    & // '      { 1.0, 2.0, 3.0 },' // nl           &
    & // '      {' // nl                            &
    & // '        fun = vel_fun,' // nl             &
    & // '        shape = shape' // nl              &
    & // '      }' // nl                            &
    & // '    }' // nl                              &
    & // '  }' // nl                                &
    & // '}' // nl                                  &
    & // 'track_variable = {' // nl                 &
    & // '  "density",' // nl                       &
    & // '  "density_shape",' // nl                 &
    & // '  "density_multi_stfun",' // nl           &
    & // '  "velocity_shape",' // nl                &
    & // '  "velocity_multi_stfun",' // nl          &
    & // '  "velocity"' // nl                       &
    & // '}' // nl

  ! write output to screen
  logical, parameter :: dumpRes = .true.
  ! number of elements to track
  integer, parameter :: nElems_track = 5
  ! position in global treeID list to track
  integer, dimension(nElems_track), parameter :: elemPos = (/ 1, 3, 5, 7, 8 /)

  ! reference values
  real(kind=rk), dimension(nElems_track), parameter :: &
    & dens = (/ 0.25_rk,  0.25_rk, 0.25_rk, 0.25_rk, 0.75_rk /)
  real(kind=rk), dimension(nElems_track), parameter :: &
    & dens_shape = (/ 0.25_rk,  0.25_rk, 0.0_rk, 0.0_rk, 0.0_rk /)
  real(kind=rk), dimension(nElems_track), parameter :: &
    & dens_multi_stfun = (/ 10.25_rk,  10.25_rk, 10.0_rk, 10.0_rk, 10.0_rk /)
  real(kind=rk), dimension(nElems_track,3), parameter ::                &
    & vel = reshape((/ 0.25_rk,  0.25_rk, 0.25_rk, 0.25_rk, 0.75_rk,    &
    &                  0.25_rk,  0.75_rk, 0.25_rk, 0.75_rk, 0.75_rk,    &
    &                  0.25_rk,  0.25_rk, 0.75_rk, 0.75_rk, 0.75_rk /), &
    &       (/nElems_track,3/) )
  real(kind=rk), dimension(nElems_track,3), parameter ::                 &
    & vel_shape = reshape((/ 1.0_rk,  1.0_rk, 0.0_rk, 0.0_rk, 0.0_rk,    &
    &                        2.0_rk,  2.0_rk, 0.0_rk, 0.0_rk, 0.0_rk,    &
    &                        3.0_rk,  3.0_rk, 0.0_rk, 0.0_rk, 0.0_rk /), &
    &             (/nElems_track,3/))
  real(kind=rk), dimension(nElems_track,3), parameter ::               &
    & vel_multi_stfun = reshape((/                                     &
    &                   1.25_rk,  1.25_rk, 1.0_rk, 1.0_rk, 1.0_rk,     &
    &                   2.25_rk,  2.75_rk, 2.0_rk, 2.0_rk, 2.0_rk,     &
    &                   3.25_rk,  3.25_rk, 3.0_rk, 3.0_rk, 3.0_rk  /), &
    &                   (/nElems_track,3/))

  type solver_type
    integer :: nDofs
    real(kind=rk), allocatable :: state(:)
    type(treelmesh_type) :: tree
    type(tem_bc_prop_type) :: boundary
    type(tem_general_type) :: general
  end type solver_type

  type tracking_type
    type(grw_intArray_type) :: varPos
    integer :: nRequestedVars
    character(len=labelLen), allocatable :: variable(:)
  end type tracking_type

  type var_index
    integer :: idx(nElems_track)
  end type

  type(tem_varSys_type) :: varSys
  type(tem_st_fun_linkedList_type) :: st_funList
  type(solver_type), target :: solver
  type(tracking_type) :: tracking
  type(tem_variable_type), allocatable :: newVar(:)
  integer :: addedPos
  integer :: iElem
  integer :: iVar
  real(kind=rk), allocatable :: res(:)
  integer, allocatable :: error(:)
  real(kind=rk) :: point(nElems_track,3)
  type(var_index), allocatable :: varIdx(:)

  write(*,*) 'Hello from tem_varSys_test'
  ! load utest mesh
  call load_env( tree     = solver%tree,     &
    &            boundary = solver%boundary, &
    &            general  = solver%general   )

  call tem_logging_init( level = 10, rank = 0 )
  write(*,*) 'nElems ', solver%tree%nElems

  allocate(solver%general%solver%conf(1))
  call load_config( conf             = solver%general%solver%conf(1), &
    &               chunk            = trim(sysConf),                 &
    &               stfun_linkedList = st_funList,                    &
    &               tracking         = tracking,                      &
    &               newVar           = newVar                         )
  solver%nDofs = 1

  ! initialize variable system
  write(*,*) 'calling varsys init'
  call tem_varSys_init(me = varSys, systemName = 'utest')

  write(*,*)
  write(*,*) 'add new variables. nVars:', size(newVar)
  call tem_varSys_append_luaVar( luaVar     = newVar,    &
    &                            varSys     = varSys,    &
    &                            st_funList = st_funList )

  write(*,*)

  write(*,*) 'Create varMap for tracking'
  call init(tracking%varPos)
  do iVar = 1, tracking%nRequestedVars
    addedPos = PositionOfVal( me  = varSys%varname, &
      &                       val = trim(tracking%variable(iVar)) )
    if (addedPos > 0) then
      call append( me = tracking%varPos, val = addedPos )
    else
      write(*,*) 'Variable ', trim(tracking%variable(iVar)), &
        & ' not found in varSys'
    end if
  end do
  write(*,*)

  write(*,*) 'Create subtree for all space time functions stored '
  write(*,*) 'in linked list of spacetime function'
  call tem_create_subTree_of_st_funList( &
    & me      = st_funList,              &
    & tree    = solver%tree,             &
    & bc_prop = solver%boundary          )
  write(*,*) 'Done creating subtree for all spacetime functions'

  write(*,*) 'Setup indices '
  do iElem = 1, nElems_track
    point(iElem,:) = tem_BaryOfId( tree   = solver%tree,            &
      &                            treeID = solver%tree             &
      &                                     %treeID(elemPos(iElem)) )
  end do

  ! idx is different for all variables.
  ! idx>0 only if point is part of subTree shape
  allocate(varIdx(tracking%nRequestedVars))
  do iVar=1,tracking%varPos%nVals
    addedPos = tracking%varPos%val(iVar)
    call varSys%method%val(addedPos)%setup_indices( &
      & varSys = varSys,                            &
      & point  = point,                             &
      & iLevel = 1,                                 &
      & tree   = solver%tree,                       &
      & nPnts  = nElems_track,                      &
      & idx    = varIdx(addedPos)%idx               )
    write(*,*) 'iVar ', iVar, ' name '//trim(newVar(addedPos)%label)
    write(*,*) 'idx ', varIdx(addedPos)%idx
  end do

  allocate(error(tracking%varPos%nVals))
  write(*,*) 'Checking variable extraction '
  do iVar=1,tracking%varPos%nVals
    if (allocated(res))  deallocate(res)
    addedPos = tracking%varPos%val(iVar)
    write(*,*) 'track var: ', trim(varSys%varname%val(addedPos))
    allocate(res(nElems_track*varSys%method%val(addedPos)%nComponents &
      & *solver%nDofs))
    ! access state array for given elemPos
    write(*,*) 'Checking get_element '
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

    write(*,*) 'Checking get_point '
    call varSys%method%val(addedPos)%get_point(                          &
      &                                  varSys = varSys,                &
      &                                  point  = point,                 &
      &                                  time   = solver%general         &
      &                                                 %simControl%now, &
      &                                  tree   = solver%tree,           &
      &                                  nPnts  = nElems_track,          &
      &                                  res    = res                    )

    call check_res( varname = varSys%varname%val(addedPos) , &
      &             res     = res,                           &
      &             error   = error(iVar)                    )

    if (error(iVar) == -1) exit

    write(*,*) 'Checking get_valOfIndex '
    write(*,*) 'idx ', varIdx(iVar)%idx
    call varSys%method%val(addedPos)%get_valOfIndex(                     &
      &                                  varSys = varSys,                &
      &                                  time   = solver%general         &
      &                                                 %simControl%now, &
      &                                  iLevel = 1,                     &
      &                                  idx    = varIdx(addedPos)%idx,  &
      &                                  nVals  = nElems_track,          &
      &                                  res    = res                    )

    call check_res( varname = varSys%varname%val(addedPos) , &
      &             res     = res,                           &
      &             error   = error(iVar)                    )

    if (error(iVar) == -1) exit

    !do iElem = 1, nElems_track
    !  write(*,*) 'iElem ', iElem, 'elemPos ', elemPos(iElem)
    !  write(*,*) 'res ', &
    !    & res( (iElem-1)*varSys%method%val(addedPos)%nComponents + 1 : &
    !    &      iElem*varSys%method%val(addedPos)%nComponents )
    !end do
  end do

  if (any(error/=0)) then
    write(*,*) 'FAILED'
  else
    write(*,*) 'PASSED'
  end if

  call close_config(L=solver%general%solver%conf(1))
  call fin_env()


contains


  ! ***************************************************************************!
  !> load variables defined in sysConf
  subroutine load_config( conf, chunk, stfun_linkedList, tracking, newVar )
    ! -------------------------------------------------------------------------!
    type(flu_state) :: conf
    character(len=*), intent(in) :: chunk
    type(tem_st_fun_linkedList_type), intent(out) :: stfun_linkedList
    type(tracking_type), intent(out) :: tracking
    type(tem_variable_type), allocatable, intent(out) :: newVar(:)
    ! -------------------------------------------------------------------------!
    integer :: nVars, iVar, thandle, iError
    integer, allocatable :: vError(:)
    ! -------------------------------------------------------------------------!
    call open_config_chunk(L=conf, chunk=trim(chunk))

    call tem_variable_load( me     = newVar,     &
      &                     conf   = conf,       &
      &                     key    = 'variable', &
      &                     vError = vError      )
    if ( all(vError /= 0) ) then
      write(*,*) "Loading variable table failed"
      stop
    end if

    call aot_table_open( L       = conf,            &
      &                  thandle = thandle,         &
      &                  key     = 'track_variable' )
    nVars = aot_table_length( L = conf, thandle = thandle )
    tracking%nRequestedVars = nVars
    allocate(tracking%variable(nVars))
    write(*,*) 'track_variables :'
    do iVar = 1, nVars
      ! Get the names of the variable
      call aot_get_val( L       = conf,                    &
        &               thandle = thandle,                 &
        &               val     = tracking%variable(iVar), &
        &               ErrCode = iError,                  &
        &               pos     = iVar                     )
      write(*,*) trim(tracking%variable(iVar))
      if (iError /= 0) then
        write(*,"(A,I2,A)") "Getting value for variable ", ivar, " failed."
        exit
      end if
    end do
    write(*,*)
    call aot_table_close( L = conf, thandle = thandle)

  end subroutine load_config
  ! ***************************************************************************!

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
    real(kind=rk) :: diff_vel(nElems_track,3)
    ! ---------------------------------------------------------------------- !
    error = 0
    select case (trim(varname))
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
    case ('density_shape')
      diff_dens = res(:) - dens_shape
      if ( any(diff_dens > eps) ) then
        error = -1
        write(*,*) 'Refernece value for density_shape does not match'
        write(*,*) 'diff ', diff_dens
      end if
      if (dumpRes) then
        write(*,*) 'reference: ', dens_shape
        write(*,*) '   output: ', res
      end if
    case ('density_multi_stfun')
      diff_dens = res(:) - dens_multi_stfun
      if ( any(diff_dens > eps) ) then
        error = -1
        write(*,*) 'Refernece value for density_multi_stfun does not match'
      end if
      if (dumpRes) then
        write(*,*) 'reference: ', dens_multi_stfun
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
    case ('velocity_shape')
      do iElem = 1, nElems_track
        diff_vel(iElem,:) = res((iElem-1)*3 + 1 : (iElem-1)*3 + 3 ) &
          &               - vel_shape(iElem,:)
      end do
      if ( any(diff_vel > eps) ) then
        error = -1
        write(*,*) 'Refernece value for velocity_shape does not match'
      end if
      if (dumpRes) then
        do iElem = 1, nElems_track
          write(*,*) 'reference: ', vel_shape(iElem,:)
          write(*,*) '   output: ', res( (iElem-1)*3 + 1 : (iElem-1)*3 + 3 )
        end do
      end if
    case ('velocity_multi_stfun')
      do iElem = 1, nElems_track
        diff_vel(iElem,:) = res((iElem-1)*3 + 1 : (iElem-1)*3 + 3 ) &
          &               - vel_multi_stfun(iElem,:)
      end do
      if ( any(diff_vel > eps) ) then
        error = -1
        write(*,*) 'Refernece value for velocity_multi_stfun does not match'
      end if
      if (dumpRes) then
        do iElem = 1, nElems_track
          write(*,*) 'reference: ', vel_multi_stfun(iElem,:)
          write(*,*) '   output: ', res( (iElem-1)*3 + 1 : (iElem-1)*3 + 3 )
        end do
      end if
    case default
      write(*,*) 'Unknown variable name', trim(varname)
      error = -1
    end select

   end subroutine check_res
  ! ***************************************************************************!

end program tem_varSys_stfunVar_test


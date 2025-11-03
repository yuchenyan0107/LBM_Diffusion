! DUMMY Implementation for the PreCICE interface!
module tem_precice_module
  use, intrinsic :: iso_c_binding
  use env_module,                   only: rk, labelLen, globalMaxLevels

  use aotus_module,                 only: flu_State, aot_get_val, &
    &                                     aoterr_NonExistent, aoterr_fatal
  use aot_table_module,             only: aot_table_open, aot_table_close
  use aot_table_ops_module,         only: aot_table_length

  use tem_logging_module,           only: logUnit
  use tem_aux_module,               only: tem_abort, tem_open_distconf
  use tem_tools_module,             only: upper_to_lower
  use tem_time_module,              only: tem_time_type
  use tem_operation_module,         only: tem_indexLvl_type
  use tem_precice_interpolation_module, &
    &                               only: tem_interpolation_type

  implicit none

  private

  logical, parameter :: precice_available = .false.

  type tem_precice_type
    ! general info for precice
    ! name of solver, should be indetical to name in precice_configeFile
    character(len=labelLen) :: accessor
    character(len=labelLen) :: precice_configFile
    !> for initialize data for precice
    integer                 :: action_isrequired
    !> for exchanging data in intermediate timestep, when using ssprk2
    logical                 :: use_RK2_inter
  end type tem_precice_type

  type(tem_precice_type), public :: precice_handle

  !> since points on the bc could be level depended also the positions ID are
  !! levelwise
  type tem_precice_posIDlvl
    !> position IDs for the exchange points, which are unquie per meshID
    !! dimension is npnts on that level
    integer, allocatable :: posIDs(:)
  end type tem_precice_posIDlvl

    !> data type which consists of information for the variables to
  !! couple with precice
  type tem_precice_varInfo
    !> coupling is done via surface, and the ID is the precice given ID of that
    !! surface, which is required for read and write to precice, when using
    !non-equidistant poins, thus we provide for read and write different meshIDs
    integer :: meshID
    !> number of variables to read from precice
    integer :: nVars = 0
    !> name of tha variable which need to be coupled
    character(len=labelLen), allocatable :: Names(:)
    !> name of variable is converted to data ID, which is unquie in precice
    integer, allocatable :: IDs(:)
    !> varPos in the ateles variable system
    integer, allocatable :: varPos(:)
    !> type which contains the posIDs levelwise, this is important to be here,
    !when using equidistant points, since we provide different meshes, thus we
    !also need different posIDs for them
    type(tem_precice_posIDlvl) :: posIDLvl(globalMaxLevels)
  end type tem_precice_varInfo

  !> datatype stored all information required for the coupled variable
  type tem_precice_coupling_type
    ! label for the ateles boundary to write from and read to
    ! required to get the correct exchange points during initialization
    character(len=labelLen) :: boundary
    !> flag for reading from precice
    logical :: read = .false.
    !> flag for writing to  precice
    logical :: write = .false.
    !> data type to store variabels to precice.
    !> data type to stroe variable information to read from precice
    type (tem_precice_varInfo) :: readVar
    !> data type to stroe variable information to read from precice
    type (tem_precice_varInfo) :: writeVar
    !> type which contains indices to access point position
    type(tem_indexLvl_type) :: pntIndex
    type(tem_interpolation_type) :: interpolation
  end type tem_precice_coupling_type

  interface tem_precice_read
    module procedure tem_precice_read_scalar
    module procedure tem_precice_read_vector
  end interface tem_precice_read

  interface tem_precice_write
    module procedure tem_precice_write_scalar
!    module procedure tem_precice_write_vector
  end interface tem_precice_write

  public :: precice_available
  !!public :: precice_handle
  public :: tem_precice_coupling_type
  public :: tem_precice_load, tem_precice_init, tem_precice_create
  !!Ne public :: tem_precice_init, tem_precice_create
  public :: tem_precice_advance
  public :: tem_precice_finalize
  public :: tem_precice_ongoing
  public :: tem_precice_read
  public :: tem_precice_write
  public :: tem_precice_set_vertex_pos
  public :: tem_precice_set_vertices_pos
  public :: tem_precice_set_edge
  public :: tem_precice_set_triangle
  public :: tem_precice_action_required
  public :: tem_precice_action_write_initData
  public :: tem_precice_action_write_iter_checkp
  public :: tem_precice_action_read_iter_checkp
  public :: tem_precice_initialize_data
  public :: tem_precice_fulfilled_Action
  public :: tem_precice_initmeshID
  public :: tem_precice_get_dataIDs
  public :: tem_precice_load_coupling

contains

  ! ************************************************************************ !
   subroutine tem_precice_load( conf )
     ! -------------------------------------------------------------------- !
     !> The filename to read the configuration from.
     type(flu_state), intent(in) :: conf
     ! -------------------------------------------------------------------- !
     integer :: thandle
     ! -------------------------------------------------------------------- !

     ! get the name of the solver and the config file from the lua script
     call aot_table_open( L       = conf,    &
       &                  thandle = thandle, &
       &                  key     = 'precice')

     ! Nothing to do for dummy!
     call aot_table_close( L       = conf,   &
       &                   thandle = thandle )
   end subroutine tem_precice_load
  ! ************************************************************************ !


  ! ************************************************************************ !
  subroutine tem_precice_create(rank, comm_size )
    ! -------------------------------------------------------------------- !
    integer, intent(in) :: rank
    integer, intent(in) :: comm_size
    ! -------------------------------------------------------------------- !
    character(len=labelLen, kind=c_char) :: c_accessor
    character(len=labelLen, kind=c_char) :: c_precice_configFile
    integer(kind=c_int) :: c_my_rank
    integer(kind=c_int) :: c_my_size
    integer(kind=c_int) :: c_accessorNameLength
    integer(kind=c_int) :: c_configFileNameLength
    !> The filename to read the configuration from.
    ! -------------------------------------------------------------------- !
    c_accessor = precice_handle%accessor
    c_precice_configFile = precice_handle%precice_configFile
    c_my_rank = int(rank, kind=c_int)
    c_my_size = int(comm_size, kind=c_int)
    c_accessorNameLength = len_trim(precice_handle%accessor)
    c_configFileNameLength = len_trim(precice_handle%precice_configFile)

  end subroutine tem_precice_create
  ! ************************************************************************ !


  ! ************************************************************************ !
  subroutine tem_precice_load_coupling( me, conf, thandle )
    ! -------------------------------------------------------------------- !
    !> The coupling type which should be filled
    type(tem_precice_coupling_type), intent(inout) :: me
    !> Lua script to obtain the configuration data from.
    type(flu_State), intent(in) :: conf
    !> Boundary condition sub table
    integer, intent(in) :: thandle
    ! -------------------------------------------------------------------- !
    character(len=labelLen) :: preciceMeshRead
    character(len=labelLen) :: preciceMeshWrite
    integer :: iError
    integer, allocatable :: vError(:)
    integer :: iVar
    logical :: varname_failed, mesh_failed
    ! -------------------------------------------------------------------- !
    call aot_get_val( L       = conf,                    &
      &               thandle = thandle,                 &
      &               key     = 'use_NP',                &
      &               val     = me%interpolation%use_NP, &
      &               default = .false.,                 &
      &               ErrCode = iError                   )
    call aot_get_val( L       = conf,                           &
      &               thandle = thandle,                        &
      &               key     = 'use_EQ_points',                &
      &               val     = me%interpolation%use_EQ_points, &
      &               default = .false.,                        &
      &               ErrCode = iError                          )

    ! ---- read from precice ---------------------------------------------------
    ! Check if there are variables to write to precice
    call aot_get_val( val       = me%readVar%names, &
      &               ErrCode   = vError,           &
      &               maxLength = 100,              &
      &               L         = conf,             &
      &               thandle   = thandle,          &
      &               key       = 'read_varname'    )
    if ( any(btest(vError, aoterr_Fatal)) ) then
      call tem_abort( 'Could not load "read_varname" for precice coupling' )
    end if
    write(logUnit(3),*) 'For precice coupling variables to read are ',&
      &                  me%readVar%names
    ! store the amount of variables to read from precice
    me%readVar%nVars = size(me%readVar%names)
    write(logUnit(3),*) 'Number of variables to read from precice ', &
      &                 me%readVar%nVars
    call aot_get_val( L       = conf,               &
      &               thandle = thandle,            &
      &               val     = preciceMeshRead,    &
      &               ErrCode = iError,             &
      &               key     = 'precice_meshRead', &
      &               default = 'none'              )
    if ( iError .ne. 0 ) then
      call tem_abort( 'No coupling mesh for read from precice is defined,' &
        & // ' coupling can not work'                                      )
    end if
    ! ---- end read from precice ----------------------------------------------

    ! ---- write to precice --------------------------------------------------
    ! check if there are variables to write to precice
    call aot_get_val( val       = me%writeVar%names, &
      &               ErrCode   = vError,            &
      &               maxLength = 100,               &
      &               L         = conf,              &
      &               thandle   = thandle,           &
      &               key       = 'write_varname'    )
    ! aot_get_val returns a bitfield containing the errors, if any. But we don't
    ! need to know about a specific error but any. Therefore we just check for
    ! all possible errors.
    !varname_failed = any(btest(vError, aoterr_Fatal))
    varname_failed = size(vError, 1) == 0 .or. sum(vError) /= 0
    call aot_get_val( L       = conf,               &
      &               thandle = thandle,            &
      &               val     = preciceMeshWrite,   &
      &               ErrCode = iError,             &
      &               key     = 'precice_meshWrite',&
      &               default = 'none'              )
    mesh_failed = iError .ne. 0

    ! The user is not forced to specify information for writing to precice.
    ! As we don't know whether he wants to specify it, we can't fail when it is
    ! missing. Thus we take the presence of one of the two needed information
    ! (either mesh or variables) as a hint that he wants to write and only then
    ! fail if the other is missing.
    if ( (varname_failed .and..not. mesh_failed) &
      & .or. (mesh_failed .and..not. varname_failed) ) then
      write(*,*) "mesh_failed ", mesh_failed, " varname_failed ", varname_failed
      call tem_abort( 'Could not load all information necessary to initialize' &
        & // ' writing to precice'                                             )
    end if

    if ( .not. varname_failed .and..not. mesh_failed ) then
      write(logUnit(3),*) 'Variables to write to precice are: ', &
        &                  me%writeVar%names
      ! store the amount of variables to write to precice
      me%writeVar%nVars = size(me%writeVar%names)
      write(logUnit(3),*) 'Number of variables to write: ', me%writeVar%nVars
    else
      me%writeVar%nVars = 0
    end if
    ! ---- end write to precice ------------------------------------------------

    ! ---- initialize external library -----------------------------------------
    ! At the moment the read part is mandatory. Thus we fail when there is
    ! nothing to read.
    if (me%readVar%nVars .eq. 0) then
      call tem_abort("No read variables for coupling specified.")
    end if
    write(logUnit(7),*) 'For precice coupling mesh_name is ', preciceMeshRead
    ! convert name into mesh ID
    call tem_precice_initmeshID( meshName = preciceMeshRead,  &
      &                          meshID   = me%readVar%meshID )
    write(logUnit(7),*) 'For precice coupling mesh_ID is ', me%readVar%meshID
    me%read = .true.
    ! convert the varName to varID
    ! allocate the dataID array before
    allocate( me%readVar%IDs(me%readVar%nVars) )
    allocate( me%readVar%varPos(me%readVar%nVars) )
    write(logUnit(3),*) 'Get dataID for read from precice '
    do iVar = 1, me%readVar%nVars
      call tem_precice_get_dataIDs( dataName = me%readVar%names(iVar), &
        &                           dataID   = me%readVar%IDs(iVar),   &
        &                           meshID   = me%readVar%meshID       )
      write(logUnit(7),*) 'For read from precice: data ID for ',       &
        & trim(me%readVar%names(iVar)), ' is ', me%readVar%IDs(iVar)
    end do

    ! set the write flag if varaibles need to be written
    if  (me%writeVar%nVars .ne. 0) then
      me%write= .true.

      write(logUnit(7),*) 'Precice coupling mesh_name is ', preciceMeshWrite
      ! convert name into mesh ID
      call tem_precice_initmeshID( meshName = preciceMeshWrite,  &
        &                          meshID   = me%writeVar%meshID )
      write(logUnit(7),*) 'Precice coupling mesh_ID is ', me%writeVar%meshID
      ! convert the varName to varID
      ! allocate the dataID array before
      allocate( me%writeVar%IDs(me%writeVar%nVars) )
      allocate( me%writeVar%varPos(me%writeVar%nVars) )
      write(logUnit(3),*) 'Get dataID to write to precice '
      do iVar = 1, me%writeVar%nVars
        call tem_precice_get_dataIDs( dataName = me%writeVar%names(iVar), &
          &                           dataID   = me%writeVar%IDs(iVar),   &
          &                           meshID   = me%writeVar%meshID       )
        write(logUnit(7),*) 'For write to precice: data ID for ', & 
          &  trim(me%writeVar%names(iVar)), ' is ', me%writeVar%IDs(iVar)
      end do
    end if
    ! ---- end initialize external library -------------------------------------
  end subroutine tem_precice_load_coupling
  ! ************************************************************************ !


  ! ************************************************************************ !
  subroutine tem_precice_get_dataIDs( dataName, dataID, meshID )
    ! -------------------------------------------------------------------- !
    character(len=labelLen), intent(in) :: dataName
    integer, intent(inout)              :: dataID
    integer, intent(in)                 :: meshID
    ! -------------------------------------------------------------------- !
    ! -------------------------------------------------------------------- !
    write(*,*) trim(dataName)
    dataID = meshID
  end subroutine tem_precice_get_DataIDs
  ! ************************************************************************ !


  ! ************************************************************************ !
  subroutine tem_precice_initmeshID( meshName, meshID )
    ! -------------------------------------------------------------------- !
    character(len=labelLen), intent(in) :: meshName
    integer, intent(out)                :: meshID
    ! -------------------------------------------------------------------- !
    ! -------------------------------------------------------------------- !
    ! init the meshID
    meshID = len_trim(meshName)

  end subroutine tem_precice_initmeshID
  ! ************************************************************************ !


  ! ************************************************************************ !
  subroutine tem_precice_init( timesteplimit )
    ! -------------------------------------------------------------------- !
    real(kind=rk), intent(inout) :: timestepLimit
    ! -------------------------------------------------------------------- !
    real(kind=c_double) :: c_timestepLimit
    ! -------------------------------------------------------------------- !

    if (precice_handle%use_RK2_inter) then
      c_timestepLimit = timesteplimit*2
    else
      c_timestepLimit = timesteplimit
    end if
    timestepLimit = c_timestepLimit

  end subroutine tem_precice_init
  ! ************************************************************************ !


  ! ************************************************************************ !
  subroutine tem_precice_initialize_data()

    write(logUnit(1),*) 'Exchange initialize data in preCICE'

  end subroutine tem_precice_initialize_data
  ! ************************************************************************ !


  ! ************************************************************************ !
  subroutine tem_precice_advance( timestepLimit )
    ! -------------------------------------------------------------------- !
    real(kind=rk), intent(inout) :: timestepLimit
    ! -------------------------------------------------------------------- !
    real(kind=c_double) :: c_timestepLimit
    ! -------------------------------------------------------------------- !

    c_timestepLimit = real(timestepLimit, kind=c_double)
    ! do i need to convert the timestepLengthlimit from kind=c_double to 
    ! kind=rk ??
    if (precice_handle%use_RK2_inter) then
      timestepLimit=c_timestepLimit/2
    else
      timestepLimit=c_timestepLimit
    end if

  end subroutine tem_precice_advance
  ! ************************************************************************ !


  ! ************************************************************************ !
  subroutine tem_precice_finalize()
    write(*,*) 'tem_precice_finalize'
  end subroutine tem_precice_finalize
  ! ************************************************************************ !


  ! ************************************************************************ !
  subroutine tem_precice_ongoing( isOngoing )
    ! -------------------------------------------------------------------- !
    !---------------------------------------------------------------------------
    integer, intent(inout) :: isOngoing
    ! -------------------------------------------------------------------- !
    integer(kind=c_int) ::c_isOngoing
    ! -------------------------------------------------------------------- !

    c_isOngoing = isOngoing
    isOngoing = c_isOngoing

  end subroutine tem_precice_ongoing
  ! ************************************************************************ !


  ! ************************************************************************ !
  subroutine tem_precice_action_write_initData( nameAction )
    ! -------------------------------------------------------------------- !
    character(len=labelLen), intent(inout) :: nameAction
    ! -------------------------------------------------------------------- !
    character(len=labelLen, kind=c_char) :: c_nameAction
    integer(kind=c_int) :: c_lengthNameAction
    ! -------------------------------------------------------------------- !
    c_lengthNAmeAction = labelLen
    c_nameAction = ''
    NameAction = c_nameAction
  end subroutine tem_precice_action_write_initData
  ! ************************************************************************ !


  ! ************************************************************************ !
  subroutine tem_precice_action_write_iter_checkp( nameAction )
    ! -------------------------------------------------------------------- !
    character(len=labelLen), intent(inout) :: nameAction
    ! -------------------------------------------------------------------- !
    character(len=labelLen, kind=c_char) :: c_nameAction
    integer(kind=c_int) :: c_lengthNameAction
    ! -------------------------------------------------------------------- !
    c_lengthNAmeAction = labelLen
    c_nameAction = ''
    NameAction = c_nameAction
  end subroutine tem_precice_action_write_iter_checkp
  ! ************************************************************************ !


  ! ************************************************************************ !
  subroutine tem_precice_action_read_iter_checkp( nameAction )
    ! -------------------------------------------------------------------- !
    character(len=labelLen), intent(inout) :: nameAction
    ! -------------------------------------------------------------------- !
    character(len=labelLen, kind=c_char) :: c_nameAction
    integer(kind=c_int) :: c_lengthNameAction
    ! -------------------------------------------------------------------- !
    c_lengthNAmeAction = labelLen
    c_nameAction = ''
    NameAction = c_nameAction
  end subroutine tem_precice_action_read_iter_checkp
  ! ************************************************************************ !


  ! ************************************************************************ !
  subroutine tem_precice_fulfilled_Action( nameAction )
    ! -------------------------------------------------------------------- !
    character(len=labelLen), intent(in) :: nameAction
    ! -------------------------------------------------------------------- !
    ! -------------------------------------------------------------------- !
    write(*,*) trim(nameAction)

  end subroutine tem_precice_fulfilled_Action
  ! ************************************************************************ !


  ! ************************************************************************ !
  subroutine tem_precice_action_required( NameAction, isRequired )
    ! -------------------------------------------------------------------- !
    character(len=labelLen), intent(in) :: NameAction
    integer, intent(out) :: isRequired
    ! -------------------------------------------------------------------- !
    ! -------------------------------------------------------------------- !

    isRequired = len_trim(nameAction)

  end subroutine tem_precice_action_required
  ! ************************************************************************ !


  ! ************************************************************************ !
  subroutine tem_precice_set_vertices_pos( meshID, npoints, points, positionID )
    ! -------------------------------------------------------------------- !
    integer, intent(in) :: meshID
    integer, intent(in) :: npoints
    real(kind=rk), intent(in) :: points(3*npoints)
    integer, intent(out) :: positionID(npoints)
    ! -------------------------------------------------------------------- !
    integer(kind=c_int) :: c_meshID
    integer(kind=c_int) :: c_meshsize
    real(kind=c_double) :: c_pos(3*npoints)
    integer(kind=c_int) :: c_posID(npoints)
    ! -------------------------------------------------------------------- !

    c_meshID = int(meshID, kind=c_int)
    c_meshsize = int(npoints, kind=c_int)
    c_pos = real(points, kind=c_double)
    positionID= int(c_posID)

  endsubroutine tem_precice_set_vertices_pos
  ! ************************************************************************ !


  ! ************************************************************************ !
  subroutine tem_precice_set_vertex_pos( meshID, point, positionID )
    ! -------------------------------------------------------------------- !
    integer, intent(in) :: meshID
    real(kind=rk), intent(in) :: point(3)
    integer, intent(out) :: positionID
    ! -------------------------------------------------------------------- !
    integer(kind=c_int) :: c_meshID
    real(kind=c_double) :: c_pos(3)
    ! -------------------------------------------------------------------- !

    c_meshID = int(meshID, kind=c_int)
    c_pos = real(point, kind=c_double)
    positionID= size(c_pos)

  endsubroutine tem_precice_set_vertex_pos
  ! ************************************************************************ !


  ! ************************************************************************ !
  subroutine tem_precice_set_edge( meshID, firstVertexID, secondVertexID, &
    &                              edgeID                                 )
    ! -------------------------------------------------------------------- !
    integer, intent(in) :: meshID
    integer, intent(in) :: firstVertexID
    integer, intent(in) :: secondVertexID
    integer, intent(out) :: edgeID
    ! -------------------------------------------------------------------- !
    integer(kind=c_int) :: c_meshID
    integer(kind=c_int) :: c_firstVertexID
    integer(kind=c_int) :: c_secondVertexID
    ! -------------------------------------------------------------------- !

    c_meshID = int(meshID, kind=c_int)
    c_firstVertexID = int(firstVertexID, kind=c_int)
    c_secondVertexID = int(secondVertexID, kind=c_int)
    edgeID = int(c_firstVertexID + c_secondVertexID, kind=c_int)

  endsubroutine tem_precice_set_edge
  ! ************************************************************************ !


  ! ************************************************************************ !
  subroutine tem_precice_set_triangle( meshID, firstEdgeID, secondEdgeID, &
    &                                  thirdEdgeID                        )
    ! -------------------------------------------------------------------- !
    integer, intent(in) :: meshID
    integer, intent(in) :: firstEdgeID
    integer, intent(in) :: secondEdgeID
    integer, intent(in) :: thirdEdgeID
    ! -------------------------------------------------------------------- !
    integer(kind=c_int) :: c_meshID
    integer(kind=c_int) :: c_firstEdgeID
    integer(kind=c_int) :: c_secondEdgeID
    integer(kind=c_int) :: c_thirdEdgeID
    ! -------------------------------------------------------------------- !

    c_meshID = int(meshID, kind=c_int)
    c_firstEdgeID = int(firstEdgeID, kind=c_int)
    c_secondEdgeID = int(secondEdgeID, kind=c_int)
    c_thirdEdgeID = int(thirdEdgeID, kind=c_int)

  endsubroutine tem_precice_set_triangle
  ! ************************************************************************ !


  ! ************************************************************************ !
  function tem_precice_read_scalar( dataID, npoints, posIDs ) result( res )
    ! -------------------------------------------------------------------- !
    integer, intent(in) :: dataID
    integer, intent(in) :: npoints
    integer, intent(in) :: posIDs(npoints)
    real(kind=rk) :: res(npoints) !< return value
    ! -------------------------------------------------------------------- !
    integer(kind=c_int) :: c_dataID
    integer(kind=c_int) :: c_blocksize
    integer(kind=c_int) :: c_positionIDs(npoints)
    real(kind=c_double) :: c_datavalue(npoints)
    ! -------------------------------------------------------------------- !

    ! get and  convert int PosId to int_c Posid
    c_positionIDs = int(posIDs, kind=c_int)
    c_blocksize = int(npoints, kind=c_int)
    ! get the data id from the dataid array (pressure,velocity)
    c_dataID = int(dataID, kind=c_int)
    res = c_dataValue
  end function tem_precice_read_scalar
  ! ************************************************************************ !


  ! ************************************************************************ !
  function tem_precice_read_vector( dataID, posID, npoints, nscalar ) &
    & result( res )
    ! -------------------------------------------------------------------- !
    integer, intent(in) :: dataID
    integer, intent(in) :: posID(:)
    integer, intent(in) :: nscalar
    integer, intent(in) :: npoints
    real(kind=rk) :: res(npoints,nscalar) !< return value
    ! -------------------------------------------------------------------- !
    integer(kind=c_int) :: c_dataID
    integer(kind=c_int) :: c_positionID
    real(kind=c_double) :: c_datavalue(nscalar)
    integer :: ipoint
    ! -------------------------------------------------------------------- !
    do ipoint= 1, npoints
      ! get and  convert int PosId to int_c Posid
      c_positionID = int(posID(iPoint), kind=c_int)
      ! get the data id from the dataid array (pressure,velocity)
      c_dataID = int(dataID, kind=c_int)

      res(iPoint, :) = c_dataValue(:)
    end do

  end function tem_precice_read_vector
  ! ************************************************************************ !


  ! ************************************************************************ !
  subroutine tem_precice_write_scalar( dataID, posID, npoints, val )
    ! -------------------------------------------------------------------- !
    integer, intent(in)       :: npoints
    integer, intent(in)       :: dataID
    integer, intent(in)       :: posID(npoints)
    real(kind=rk), intent(in) :: val(npoints)
    ! -------------------------------------------------------------------- !
    integer(kind=c_int) :: c_dataID
    integer(kind=c_int) :: c_blocksize
    integer(kind=c_int) :: c_positionID(npoints)
    real(kind=c_double) :: c_datavalue(npoints)
    ! -------------------------------------------------------------------- !

    ! get and  convert int PosId to int_c Posid
    c_positionID(:) = int(posID, kind=c_int)
    ! get the data id from the dataid array (pressure,velocity)
    c_dataID = int(dataID, kind=c_int)
    c_blocksize = int(npoints, kind=c_int)

    ! get and convert the value
    c_datavalue(:) = real(val(:), kind=c_double)

  end subroutine tem_precice_write_scalar
  ! ************************************************************************ !

end module tem_precice_module

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

  !!Ne use PreCICE_solver_if_module,     only: precicef_create,             &
  use precice,     only: precicef_create,                    &
    &                    precicef_initialize,                &
    &                    precicef_initialize_data,           &
    &                    precicef_advance,                   &
    &                    precicef_finalize,                  &
    &                    precicef_set_vertex,                &
    &                    precicef_set_vertices,              &
    &                    precicef_set_edge,                  &
    &                    precicef_set_triangle,              &
    &                    precicef_ongoing,                   &
    &                    precicef_get_data_id,               &
    &                    precicef_read_sdata,                &
    &                    precicef_read_bsdata,               &
    &                    precicef_read_vdata,                &
    &                    precicef_write_sdata,               &
    &                    precicef_write_bsdata,              &
    &                    precicef_write_vdata,               &
    &                    precicef_get_mesh_id,               &
    &                    precicef_action_required,           &
    &                    precicef_mark_action_fulfilled,     &
    &                    precicef_read_data_available,       &
    &                    precicef_action_write_initial_data, &
    &                    precicef_action_write_iter_checkp,  &
    &                    precicef_action_read_iter_checkp

  implicit none

  private

  logical, parameter :: precice_available = .true.

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
     integer :: iError
     ! -------------------------------------------------------------------- !

     ! get the name of the solver and the config file from the lua script
     call aot_table_open( L       = conf,    &
       &                  thandle = thandle, &
       &                  key     = 'precice')

     if (thandle /= 0) then
       call aot_get_val( L       = conf,                    &
         &               thandle = thandle,                 &
         &               key     = 'accessor',              &
         &               val     = precice_handle%accessor, &
         &               default = 'Ateles',                &
         &               ErrCode = iError                   )
       ! use the xml config-file for precice
       call aot_get_val( L       = conf,                              &
         &               thandle = thandle,                           &
         &               key     = 'configFile',                      &
         &               val     = precice_handle%precice_configFile, &
         &               default = 'config.xml',                      &
         &               ErrCode = iError                             )
       write(logUnit(1),*) '* preCICE asseccor = "' &
         & // trim(precice_handle%accessor) // '"'
       write(logUnit(1),*) '* preCICE configFile = "' &
         & // trim(precice_handle%precice_configFile) // '"'
       call aot_get_val( L       = conf,                         &
         &               thandle = thandle,                      &
         &               key     = 'use_RK2_inter',              &
         &               val     = precice_handle%use_RK2_inter, &
         &               default = .false.,                      &
         &               ErrCode = iError                        )

     else ! no precice table
       call tem_abort( 'ERROR while init preCICE: No Precice table is provided' )
     end if
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

    call precicef_create( participantName      = c_accessor,            &
      &                   configFileName       = c_precice_configFile,  &
      &                   solverProcessIndex   = c_my_rank,             &
      &                   solverProcessSize    = c_my_size,             &
      &                   lengthAccessorName   = c_accessorNameLength,  &
      &                   lengthConfigFileName = c_configFileNameLength )

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
    character(len=labelLen) :: preciceMesh
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

    if (me%interpolation%use_EQ_points) then
      ! set how many equidistant points should be provided.
      call aot_get_val( L       = conf,                              &
        &               thandle = thandle,                           &
        &               key     = 'factor_EQ_points',                &
        &               val     = me%interpolation%factor_EQ_points, &
        &               default = 1.0_rk,                                 &
        &               ErrCode = iError                             )
      if (btest(iError, aoterr_Fatal)) then
         write(logUnit(1),*)'FATAL Error occured, while loading factor_EQ_points :'
         if (btest(iError, aoterr_NonExistent))                                   &
           & write(logUnit(1),*)'Variable not existent!'
         write(logUnit(1),*)'STOPPING'
         call tem_abort()
       end if
    end if  

!!VK    !get the boundary label in ateles
!!VK    call aot_get_val( L       = conf,             &
!!VK        &             thandle = thandle,          &
!!VK        &             val     = me%boundary,      &
!!VK        &             ErrCode = iError,           &
!!VK        &             key     = 'boundary_label', &
!!VK        &             default = 'none'          )
!!VK    if ( iError .NE. 0 ) then
!!VK        write(logUnit(1),*) ' No boundary label for precice coupling,' &
!!VK            &                 //' is defined, coupling can not work, abort...'
!!VK        call tem_abort
!!VK    end if

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
    character(len=labelLen, kind=c_char) :: c_dataName
    integer(kind=c_int)                  :: c_dataNameLength
    integer(kind=c_int)                  :: c_meshID
    integer(kind=c_int)                  :: c_dataID
    ! -------------------------------------------------------------------- !
    c_meshID = meshID
    c_dataNameLength = len_trim(dataName)
    c_dataName = dataName
    call precicef_get_data_id( dataName       = c_dataName,      &
      &                        meshID         = c_meshID,        &
      &                        dataID         = c_dataID,        &
      &                        lengthDataName = c_dataNameLength )
    dataID = c_dataID
  end subroutine tem_precice_get_DataIDs
  ! ************************************************************************ !


  ! ************************************************************************ !
  subroutine tem_precice_initmeshID( meshName, meshID )
    ! -------------------------------------------------------------------- !
    character(len=labelLen), intent(in) :: meshName
    integer, intent(out)                :: meshID
    ! -------------------------------------------------------------------- !
    character(len=labelLen, kind=c_char) :: c_MeshName
    integer(kind=c_int)                  :: c_lengthMeshName
    integer(kind=c_int)                  :: c_meshID
    ! -------------------------------------------------------------------- !
    ! init the meshID
    c_lengthMeshName = len_trim(meshName)
    c_meshName = meshName

    call precicef_get_mesh_id( meshID         = c_meshID,        &
      &                        meshName       = c_meshName,      &   
      &                        lengthMeshName = c_lengthMeshName )

    meshID = c_meshID

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
    call precicef_initialize( c_timestepLimit )
    timestepLimit = c_timestepLimit

  end subroutine tem_precice_init
  ! ************************************************************************ !


  ! ************************************************************************ !
  subroutine tem_precice_initialize_data()

    write(logUnit(1),*) 'Exchange initialize data in preCICE'
    call precicef_initialize_data()

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
      call precicef_advance(c_timestepLimit/2)
      timestepLimit=c_timestepLimit/2
    else
      call precicef_advance(c_timestepLimit)
      timestepLimit=c_timestepLimit
    end if

  end subroutine tem_precice_advance
  ! ************************************************************************ !


  ! ************************************************************************ !
  subroutine tem_precice_finalize()
    call precicef_finalize()
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
    call precicef_ongoing(c_isOngoing)
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
    call precicef_action_write_initial_data( c_nameAction, c_lengthNameAction )
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
    call precicef_action_write_iter_checkp( c_nameAction, c_lengthNameAction )
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
    call precicef_action_read_iter_checkp( c_nameAction, c_lengthNameAction )
    NameAction = c_nameAction
  end subroutine tem_precice_action_read_iter_checkp
  ! ************************************************************************ !


  ! ************************************************************************ !
  subroutine tem_precice_fulfilled_Action( nameAction )
    ! -------------------------------------------------------------------- !
    character(len=labelLen), intent(in) :: nameAction
    ! -------------------------------------------------------------------- !
    character(len=labelLen, kind=c_char) :: c_nameAction
    integer(kind=c_int) :: c_lengthNameAction
    ! -------------------------------------------------------------------- !
    c_NameAction = nameAction
    c_LengthNameAction = len_trim(nameAction)

    call precicef_mark_action_fulfilled( c_nameAction, c_lengthNameAction )
  end subroutine tem_precice_fulfilled_Action
  ! ************************************************************************ !


  ! ************************************************************************ !
  subroutine tem_precice_action_required( NameAction, isRequired )
    ! -------------------------------------------------------------------- !
    character(len=labelLen), intent(in) :: NameAction
    integer, intent(out) :: isRequired
    ! -------------------------------------------------------------------- !
    character(len=labelLen, kind=c_char) :: c_NameAction
    integer(kind=c_int) :: c_isRequired
    integer(kind=c_int) :: c_lengthAction
    ! -------------------------------------------------------------------- !

    c_nameaction = nameAction
    c_lengthAction = len_trim(NameAction)

    call precicef_action_required( action       = c_nameaction,  &
      &                            isRequired   = c_isRequired,  &
      &                            lengthAction = c_lengthAction )

    isRequired = c_isRequired

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
    call  precicef_set_vertices( meshID      = c_meshID,   &
      &                          meshsize    = c_meshsize, &
      &                          positions   = c_pos,      &
      &                          positionIDs = c_posID     )
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
    integer(kind=c_int) :: c_posID
    ! -------------------------------------------------------------------- !

    c_meshID = int(meshID, kind=c_int)
    c_pos = real(point, kind=c_double)
    call  precicef_set_vertex( meshID   = c_meshID, &
      &                        position = c_pos,    &
      &                        vertexID = c_posID   )
    positionID= int(c_posID)

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
    integer(kind=c_int) :: c_edgeID
    ! -------------------------------------------------------------------- !

    c_meshID = int(meshID, kind=c_int)
    c_firstVertexID = int(firstVertexID, kind=c_int)
    c_secondVertexID = int(secondVertexID, kind=c_int)
    call  precicef_set_edge( meshID         = c_meshID,         &
      &                      firstVertexID  = c_firstVertexID,  &
      &                      secondVertexID = c_secondVertexID, &
      &                      edgeID         = c_edgeID          )
    edgeID = int(c_edgeID)

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
    call  precicef_set_triangle( meshID         = c_meshID,       &
      &                          firstEdgeID    = c_firstEdgeID,  &
      &                          secondEdgeID   = c_secondEdgeID, &
      &                          thirdEdgeID    = c_thirdEdgeID   )

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
    integer(kind=c_int) :: c_isAvailable
    integer             :: isAvailable
    integer(kind=c_int) :: c_dataID
    integer(kind=c_int) :: c_blocksize
    integer(kind=c_int) :: c_positionIDs(npoints)
    real(kind=c_double) :: c_datavalue(npoints)
    ! -------------------------------------------------------------------- !

!!VK     call precicef_read_data_available(c_isAvailable)
!!VK     isAvailable = c_isAvailable
!!VK     if (isAvailable .eq. 1) then
    !write(*,*) 'posID ', posIDs

    ! get and  convert int PosId to int_c Posid
    c_positionIDs = int(posIDs, kind=c_int)
    c_blocksize = int(npoints, kind=c_int)
    ! get the data id from the dataid array (pressure,velocity)
    c_dataID = int(dataID, kind=c_int)
    call precicef_read_bsdata( dataID       = c_dataID,      &
      &                        blocksize    = c_blocksize,   &
      &                        valueIndices = c_positionIDs, &
      &                        values       = c_datavalue    )
    res = c_dataValue
!!VK     else
!!VK      write(*,*) 'not avaiable'
!!VK    end if
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

      call precicef_read_vdata( dataID     = c_dataID,     &
        &                       valueIndex = c_positionID, &
        &                       dataValue  = c_datavalue   )

      res(iPoint, :) = c_dataValue(:)
    end do

  end function tem_precice_read_vector
  ! ************************************************************************ !


  ! ************************************************************************ !
  subroutine tem_precice_write_scalar( dataID, posID, npoints, val )
    ! -------------------------------------------------------------------- !
    integer, intent(in)       :: dataID
    integer, intent(in)       :: posID(npoints)
    integer, intent(in)       :: npoints
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

    call precicef_write_bsdata( dataID       = c_dataID,     &
      &                         blockSize    = c_blocksize,  &
      &                         valueIndices = c_positionID, &
      &                         values       = c_datavalue   )

  end subroutine tem_precice_write_scalar
  ! ************************************************************************ !


  ! ************************************************************************ !
  subroutine tem_precice_write_vector( dataID, posID, npoints, nscalar, val )
    ! -------------------------------------------------------------------- !
    integer, intent(in)       :: dataID
    integer, intent(in)       :: posID(:)
    integer, intent(in)       :: npoints
    integer, intent(in)       :: nscalar
    real(kind=rk), intent(in) :: val(nscalar,npoints)
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

      ! convert the value
      c_datavalue(:) = real(val(:,ipoint), kind=c_double)

      call precicef_write_vdata( dataID     = c_dataID,      &
        &                        valueIndex = c_positionID,  &
        &                        dataValue  = c_datavalue(:) )
    end do

  end subroutine tem_precice_write_vector
  ! ************************************************************************ !

end module tem_precice_module

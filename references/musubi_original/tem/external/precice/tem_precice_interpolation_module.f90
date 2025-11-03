module tem_precice_interpolation_module
  use env_module,                only: rk, labelLen, globalMaxLevels

  use tem_aux_module,            only: tem_abort
  use tem_logging_module,        only: logUnit, lldebug
  use tem_tools_module,          only: tem_horizontalSpacer

  implicit none

  private

  !> Datatype to represent facewise nodes
  type tem_equadPoints_type
   real(kind=rk), allocatable :: Points(:,:)
   integer :: nQuadPoints
  end type tem_equadPoints_type

  !> Level dependent arguments
  type tem_interpolation_data_type
    !> face is an array of size nDir and iAlign (left/right)
    type(tem_equadPoints_type), allocatable :: faces(:,:)
    !> number of equidistant Points
    integer :: nPoints
    !> number of equidistant Points per direction
    integer :: nPointsPerDir
    !> Number of edges
    integer :: nEdges
    integer, allocatable :: edges(:,:)
    !> Number of triangles
    integer :: nTriangles
    integer, allocatable :: triangles(:,:)
  end type tem_interpolation_data_type

  !> data type, which cotains the information of the interpolation methods, when
  !! using equidistant points for write to precice and using nearest-projection
  !! interpolation
  type tem_interpolation_type
    !> Using the nearest-projection interpolation
    logical :: use_NP
    !> for using equidistant points for the interpolation
    logical :: use_EQ_points
    !> set the number of points, which has to be used for the interpolation with
    !equidistant points
    real(kind=rk) :: factor_EQ_points
    type(tem_interpolation_data_type) :: interpol_data(globalMaxLevels)
  end type tem_interpolation_type

  public :: tem_create_surface_equidistant_points
  public :: tem_equadPoints_type
  public :: tem_interpolation_type
  public :: tem_interpolation_data_type
  public :: build_faceNodes_edges_2D
  public :: build_faceNodes_triangles_3D
  public :: tem_create_edges_triangles

contains

  ! ************************************************************************ !
  !> Routine to generate the face edges, through given nQuadPoints per dircetion
  !! for the 2D testcase. Since we have a 2D case, we have just lines (edges) to
  !!connect the nQuadPoints. d1 is the first coupling domain and d2 the second.
  !!_______
  !!|d1 |d2|
  !!|___|__|
  subroutine build_faceNodes_edges_2D( edges, nEdges, nQuadPointsPerDir )
    ! -------------------------------------------------------------------- !
    integer, allocatable, intent(out) :: edges(:,:)
    integer, intent(out) :: nEdges
    !> quadrature points including oversampling factor
    ! per spatial direction
    integer, intent(in)  :: nQuadPointsPerDir
    ! -------------------------------------------------------------------- !
    integer :: ii ! loop incides
    integer :: iEdge ! initial edge
    ! -------------------------------------------------------------------- !
    nEdges = nQuadPointsPerDir - 1
    allocate(edges(nEdges,2))

    iEdge = 1
    do ii = 1, nEdges
      edges (iEdge,1) = ii
      edges (iEdge,2) = ii+1
      iEdge = iEdge + 1
    end do

  end subroutine build_faceNodes_edges_2D
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Routine to generate the face edges as well as triangles, through given
  !! nQuadPoints per dircetion and the necessary edges. For the coupling in 3D
  !! the coupling face is a plane. From every square two triangles can be
  !! created.
  subroutine build_faceNodes_triangles_3D( triangles, nTriangles, edges, &
    &                                      nEdges, nQuadPointsPerDir     )
    ! -------------------------------------------------------------------- !
    integer, allocatable, intent(out) :: triangles(:,:)
    integer, intent(out) :: nTriangles
    integer, allocatable, intent(out) :: edges(:,:)
    integer, intent(out) :: nEdges
    integer, intent(in) :: nQuadPointsPerDir
    ! -------------------------------------------------------------------- !
    integer :: ii, jj ! loop incides
    integer :: iEdge ! initial edge
    integer :: iTriangle ! initial triangle
    integer :: edgeID_1, edgeID_2, edgeID_3 ! local edgeID created per point
    integer :: nEdgesPerDir ! number of edges per direction
    integer :: offset
    ! -------------------------------------------------------------------- !
    nTriangles = 2*(nQuadPointsPerDir-1)**2
    nEdges = ( (nQuadPointsPerDir-1)*nQuadPointsPerDir ) * 2 &
      &      + (nQuadPointsPerDir-1)**2
    nEdgesPerDir = (nQuadPointsPerDir-1)*3+1
    allocate(edges(nEdges,2))
    allocate(triangles(nTriangles,3))

    iTriangle = 0
    iEdge = 0
    do jj = 1, nQuadPointsPerDir
      do ii = 1, nQuadPointsPerDir
        ! define the vertical edges and give them an ID
        if (ii <= nQuadPointsPerDir .and. jj <= nQuadPointsPerDir-1) then
          iEdge = iEdge + 1
          edges(iEdge,1) = (jj-1)*nQuadPointsPerDir+ii
          edges(iEdge,2) = jj*nQuadPointsPerDir+ii
          edgeID_1 = iEdge
        end if

        ! define the horizantal edges and give them an ID
        if (jj <= nQuadPointsPerDir .and. ii <= nQuadPointsPerDir-1) then
          iEdge = iEdge + 1
          edges(iEdge,1) = (jj-1)*nQuadPointsPerDir+ii
          edges(iEdge,2) = (jj-1)*nQuadPointsPerDir+ii+1
          edgeID_2 = iEdge
        end if

        ! creat diagonal edges, which are needed for the triangels
        if (ii <= nQuadPointsPerDir-1 .and. jj <= nQuadPointsPerDir-1) then
          iEdge = iEdge + 1
          edges(iEdge,1) = (jj-1)*nQuadPointsPerDir+ii
          edges(iEdge,2) = jj*nQuadPointsPerDir+ii+1
          edgeID_3 = iEdge
          ! creat two triangles per quad, each of them needs three edges
          iTriangle = iTriangle + 2
          triangles(iTriangle-1, 1) = edgeID_3
          triangles(iTriangle-1, 2) = edgeID_1
          if (jj == nQuadPointsPerDir-1) then
            if (ii == 1) then
              offset = 1
            else
              offset = offset + 2
            end if
            triangles(iTriangle-1, 3) = edgeID_2 + nEdgesPerDir - offset
          else
            triangles(iTriangle-1, 3) = edgeID_2 + nEdgesPerDir
          end if
          triangles(iTriangle, 1) = edgeID_3
          triangles(iTriangle, 2) = edgeID_1 + 3
          triangles(iTriangle, 3) = edgeID_2
        end if
      end do
    end do

  end subroutine build_faceNodes_triangles_3D
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Routine to for calling the edges and triangles, when using
  !! Nearest-Projection for the interpolation, in 2D and 3D.
  subroutine tem_create_edges_triangles( me, nQuadPointsPerDir, scheme_dim )
    ! -------------------------------------------------------------------- !
    type(tem_interpolation_data_type), intent(inout) :: me
    integer, intent(in) :: nQuadPointsPerDir
    integer, intent(in) :: scheme_dim
    ! -------------------------------------------------------------------- !

    ! Call for the edges and triangles, which depends on the dimension
    if (scheme_dim == 2) then
      call build_faceNodes_edges_2D(            &
        & edges             = me%edges,         &
        & nEdges            = me%nEdges,        &
        & nQuadPointsPerDir = nQuadPointsPerDir )
      write(logunit(lldebug),*) 'Done initializing projection: Edges'
    end if
    if (scheme_dim == 3) then
    ! Routine to build up list of edges and triangles
      call build_faceNodes_triangles_3D(        &
        & triangles         = me%triangles,     &
        & nTriangles        = me%nTriangles,    &
        & edges             = me%edges,         &
        & nEdges            = me%nEdges,        &
        & nQuadPointsPerDir = nQuadPointsPerDir )
      write(logunit(lldebug),*) 'Done initializing projection: ' &
        & // 'Edges and Triangles'
    end if

  end subroutine tem_create_edges_triangles
  ! ************************************************************************ !


  ! ************************************************************************ !
  !>create points on reference element for 2D testcases
  function get2DFacePnts( iDir, iAlign, nPoints ) result( Points )
    ! -------------------------------------------------------------------- !
    integer, intent(in) :: iAlign
    integer, intent(in) :: idir
    integer, intent(in) :: nPoints
    real:: Points(nPoints,3)
    ! -------------------------------------------------------------------- !
    integer :: Iinterval
    ! -------------------------------------------------------------------- !
    select case(iDir)
    !> face in x direction, x coord is fixed
    case(1)
      do iInterval = 1, nPoints
        Points(iInterval, 1) = (-1.0_rk)**iAlign
        Points(iInterval, 2) = (real((iInterval - 1) * 2 + 1) / nPoints) &
          &                    - 1.0_rk
        Points(iInterval, 3) = 0.0_rk
      end do !iInterval
    !> face in y direction, y coord is fixed
    case(2)
      do iInterval = 1, nPoints
        Points(iInterval, 1) = (real((iInterval - 1) * 2 + 1) / nPoints) &
          &                    - 1.0_rk
        Points(iInterval, 2) = (-1.0_rk)**iAlign
        Points(iInterval, 3) = 0.0_rk
      end do
    end select
  end function get2DFacePnts
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> create the points on the reference element, for 3D testcase
  function get3DFacePnts( idir, iAlign, nPoints, nPointsPerDir ) result( Points )
    ! -------------------------------------------------------------------- !
    integer, intent(in) :: iAlign
    integer, intent(in) :: idir
    integer, intent(in) :: nPoints
    integer, intent(in) :: nPointsPerDir
    real :: Points(nPoints,3)
    ! -------------------------------------------------------------------- !
    integer :: iInterval
    integer :: jInterval
    integer :: iPnt
    ! -------------------------------------------------------------------- !
    iPnt = 0
    select case(idir)
    !> face in x direction, x coord is fixed
    case(1)
      do iInterval = 1, nPointsPerDir
        do jInterval = 1, nPointsPerDir
          iPnt = iPnt + 1
          Points(iPnt, 1) = (-1.0_rk)**iAlign
          Points(iPnt, 2) = (real((iInterval - 1) * 2 + 1) / nPointsPerDir) &
            &                    - 1.0_rk
          Points(iPnt, 3) = (real((jInterval - 1) * 2 + 1) / nPointsPerDir) &
            &                    - 1.0_rk
        end do
      end do
    !> face in y direction, y coord is fixed
    case(2)
      do iInterval = 1, nPointsPerDir
        do jInterval = 1, nPointsPerDir
          iPnt = iPnt + 1
          Points(iPnt, 1) = (real((iInterval -1) * 2 + 1) / nPointsPerDir) &
            &                    - 1.0_rk
          Points(iPnt, 2) = (-1.0_rk)**iAlign
          Points(iPnt, 3) = (real((jInterval - 1) * 2 + 1) / nPointsPerDir) &
            &                    - 1.0_rk
        end do
      end do
    !> face in z direction, z coord is fixed
    case(3)
      do iInterval = 1, nPointsPerDir
        do jInterval = 1, nPointsPerDir
          iPnt = iPnt + 1
          Points(iPnt, 1) = (real((iInterval - 1) * 2 + 1) / nPointsPerDir) &
            &                    - 1.0_rk
          Points(iPnt, 2) = (real((jInterval - 1) * 2 + 1) / nPointsPerDir) &
            &                    - 1.0_rk
          Points(iPnt, 3) = (-1.0_rk)**iAlign
        end do
      end do
    end select
  end function get3DFacePnts
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Routine to generate a list of equidistant points
  subroutine tem_create_surface_equidistant_points( me, nDir, nPoly, nfactor )
    ! -------------------------------------------------------------------- !
    type(tem_interpolation_data_type), intent(out) :: me
    !> store points to the corresponding spatial directions
    integer, intent(in) :: nDir
    !> polynomial degree
    integer, intent(in) :: nPoly
    type(tem_interpolation_type), intent(in) :: nfactor
    ! -------------------------------------------------------------------- !
    integer :: iAlign
    integer :: idir
    !> loop variable for the number of intervals
    integer :: iInterval
    !> loop variable for the number of intervals
    integer :: jInterval
    !>factor to increase number of points, default 1.0
    real(kind=rk) :: factor
    integer :: iPnt
    ! -------------------------------------------------------------------- !
    factor= nfactor%factor_EQ_points
    write(logunit(lldebug),*) "Chosen factor for the equidistant points :", &
      & factor
    if (nDir>1) then
      me%nPointsPerDir = nint((nPoly + 1.0) * factor) 
      me%nPoints = (nint((nPoly + 1.0) * factor))**(nDir - 1)
    else
      me%nPoints = 1
    end if
    write(logunit(lldebug),*) ' Number of equidistant Points per Face', &
      & me%nPoints
    allocate( me%faces(nDir,2) )
    do idir = 1, nDir
      do iAlign = 1,2
        allocate( me%faces(idir,iAlign)%Points(me%nPoints, 3) )
        select case (nDir)
        case(1)
          me%faces(idir, iAlign)%Points(1,1) = (-1.0_rk)**iAlign
          me%faces(idir, iAlign)%Points(1,2) = 0.0_rk
          me%faces(idir, iAlign)%Points(1,3) = 0.0_rk
        case(2)
          me%faces(idir,iAlign)%Points = get2DFacePnts( &
            & iDir    = idir,                           &
            & iAlign  = iAlign,                         &
            & nPoints = me%nPointsPerDir                )
        case(3)
          me%faces(idir,iAlign)%Points = get3DFacePnts( & 
            & iDir    = idir,                           &
            & iAlign  = iAlign,                         &
            & nPoints = me%nPoints,                     &
            & nPointsPerDir = me%nPointsPerDir          )
        do iPnt = 1, me%nPoints 
        end do   

        end select
      end do
    end do
  end subroutine tem_create_surface_equidistant_points
  ! ************************************************************************ !

end module tem_precice_interpolation_module

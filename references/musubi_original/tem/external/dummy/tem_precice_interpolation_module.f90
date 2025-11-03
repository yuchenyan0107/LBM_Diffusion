module tem_precice_interpolation_module
  use env_module,                only: rk, labelLen, globalMaxLevels

  use tem_aux_module,            only: tem_abort
  use tem_logging_module,        only: logUnit
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
  !using equidistant points for write to precice and using nearest-projection
  !interpolation
  type tem_interpolation_type
    !> Using the nearest-projection interpolation
    logical :: use_NP
    !> for using equidistant points for the interpolation
    logical :: use_EQ_points
    !> set the number of points, which has to be used for the interpolation with
    !equidistant points
    integer :: factor_EQ_points
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
    nEdges = size(edges)*nQuadPointsPerDir ! Avoid warnings about unused
  end subroutine build_faceNodes_edges_2D
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Routine to generate the face edges as well as triangles, through given
  !! nQuadPoints per dircetion and the necessary edges. For the coupling in 3D
  !! the coupling face is a plane. From every square two triangles can be
  !! created.
  subroutine build_faceNodes_triangles_3D( triangles, nTriangles, edges, &
    &                                      nEdges, nQuadPointsPerDir )
    ! -------------------------------------------------------------------- !
    integer, allocatable, intent(out) :: triangles(:,:)
    integer, intent(out) :: nTriangles
    integer, allocatable, intent(out) :: edges(:,:)
    integer, intent(out) :: nEdges
    integer, intent(in) :: nQuadPointsPerDir
    ! -------------------------------------------------------------------- !
    ! Avoid warnings about unused
    allocate(triangles(0,nQuadPointsPerDir))
    allocate(edges(0,nQuadPointsPerDir))
    nTriangles = size(triangles)
    nEdges = size(edges)
  end subroutine build_faceNodes_triangles_3D
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Routine to for calling the edges and triangles, when using
  !! Nearest-Projection for the interpolation, in 2D and 3D.
  subroutine tem_create_edges_triangles (me, nQuadPointsPerDir, scheme_dim)
    ! -------------------------------------------------------------------- !
    type(tem_interpolation_data_type), intent(inout) :: me
    integer, intent(in) :: nQuadPointsPerDir
    integer, intent(in) :: scheme_dim
    ! -------------------------------------------------------------------- !
    me%nPoints = scheme_dim*nQuadPointsPerDir
  end subroutine tem_create_edges_triangles
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Routine to generate a list of equidistant points
  subroutine tem_create_surface_equidistant_points( me, nDir, nPoly, nfactor )
    ! -------------------------------------------------------------------- !
    type(tem_interpolation_data_type), intent(out) :: me
    !> store points to the corresponding spatial directions
    integer, intent(in)  :: nDir
    !> polynomial degree
    integer, intent(in) :: nPoly
    type(tem_interpolation_type), intent(in) :: nfactor
    ! -------------------------------------------------------------------- !
    ! Just to avoid a compiler warning for unused arguments
    me%nPoints = nDir*nPoly*nFactor%factor_EQ_points
  end subroutine tem_create_surface_equidistant_points
  ! ************************************************************************ !

end module tem_precice_interpolation_module

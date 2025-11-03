! Copyright (c) 2019 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
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
! *******************************************************************************!
!> summary: This module contains \ref points description and growing array
!! author: Kannan Masilamani
!! description for points.

?? include 'arrayMacros.inc'

module tem_point_module
  use env_module,         only: rk, minLength, zeroLength
  use tem_cube_module,    only: tem_cube_type
  use tem_logging_module, only: tem_toStr, logUnit

  implicit none

  private

  public :: grw_pointArray_type
  public :: tem_point_type
  public :: init, append, truncate, destroy, empty, placeAt
  public :: tem_pointCubeOverlap

  !> This type contains coordinate of a point
  type tem_point_type
    real(kind=rk) :: coord(3) !< real world coordinate of a point
  end type tem_point_type

?? copy :: GA_decltxt(point, type(tem_point_type))

contains

! **************************************************************************** !
  !> This function checks whether the given point is found inside given cube.
  !!
  !! Point is inside the cube only if the point is >= cube origin and
  !! < cube max. Point lying on the cube max is not part of the cube
  function tem_pointCubeOverlap(point, cube) result(overlap)
    ! --------------------------------------------------------------------------!
    !> Coordinate of the point to check for intersection.
    type(tem_point_type) :: point
    !> Cube to intersect with.
    type(tem_cube_type) :: cube
    logical :: overlap !< true if point lies inside else false
    ! --------------------------------------------------------------------------!
    logical :: dirrange(3)

    ! Check interval in all 3 directions.
    dirrange = (point%coord >= cube%origin) &
      &        .and. (point%coord < cube%endPnt)

    ! Overlap depends on all 3 intervals.
    overlap = all(dirrange)

  end function tem_pointCubeOverlap
! **************************************************************************** !

?? copy :: GA_impltxt(point, type(tem_point_type), type(tem_point_type))

end module tem_point_module

!> \page point Points
!! Points are defined in the configuration file through canonoical
!! geometry kind with an origin.\n
!! Valid definition:
!! \li Single point
!! \verbatim
!! geometry = {
!!   kind = 'canoND',
!!   object = {
!!     origin = { 0.0,0.0,0.0 }
!!   }
!! }
!! \endverbatim
!! \li Multiple point
!! \verbatim
!! geometry = {
!!   kind = 'canoND',
!!   object = {
!!     {
!!     origin = { 0.0,0.0,0.0 }
!!     },
!!     {
!!     origin = { 1.0,0.0,0.0 }
!!     },
!!   }
!! }
!! \endverbatim
!!\n\n
!! Seeder file to generate the mesh with multiple point obstacle is given below:
!! \include testsuite/point/seeder.lua
!! \n\n
!! The image generated with multiple point obstacles from the above code:
!! \image html tem_point.png
!! Example lua file is available at \link testsuite/point/seeder.lua
!! \example testsuite/point/seeder.lua

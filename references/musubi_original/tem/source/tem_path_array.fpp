! Copyright (c) 2013 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2013, 2019 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2013 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
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
! ****************************************************************************** !
!> Module to provide dynamic structures for treeID elements
!!
!! The dynamic arrays provided by this module are
!! capable of handling lists of values, which might
!! need to grow over time.
!! Removal of entries is not possible directly.
!! A ranking list for the indices is maintained to
!! fast lookups of given values by a binary search.
!! This allows for the efficient handling of lists
!! with unique entries.
!! The complete module might be put into a CoCo Text
!! template, to create new modules of this object
!! for different types. For now, two different
!! templates are used for the declaration part and
!! the implementation part.
!!
?? include 'arrayMacros.inc'
module tem_path_array_module

  ! include treelm modules
  use env_module,          only: minLength, zeroLength
  use tem_topology_module, only: tem_PathComparison, tem_path_type

  implicit none

! -----------------------------------------------------------------
! In order to define this data structure easily for several data
! type, we use the Coco text copying feature here, to duplicate
! the necessary declarations.
! tname: indicates type of dynamic array (long, int, real, ...)
! tstring: is the actual string describing the type specification

  public :: dyn_pathArray_type
  public :: init, append, destroy, empty, sortTruncate, expand
  public :: positionOfVal, sortedPosOfVal
  public :: operator(==), operator(/=)
  public :: operator(<), operator(<=), operator(>), operator(>=)

  interface operator(==)
    module procedure isEqual
  end interface

  interface operator(/=)
    module procedure isUnequal
  end interface

  interface operator(<)
    module procedure isSmaller
  end interface

  interface operator(<=)
    module procedure isSmallerOrEqual
  end interface

  interface operator(>)
    module procedure isGreater
  end interface

  interface operator(>=)
    module procedure isGreaterOrEqual
  end interface

?? copy :: DA_decltxt(path, type(tem_path_type))

  contains
! ****************************************************************************** !
  !> Include the subroutines for the dynamic array.
?? copy :: DA_impltxt(path, type(tem_path_type), type(tem_path_type))

! ****************************************************************************** !
  !> This function provides the test for equality of two treeIDs.
  !!
  !! Two treeIDs are equal if their integer is equal
  function isEqual(left, right) result(equality)
    ! ---------------------------------------------------------------------------
    !> path to compare
    type(tem_path_type), intent(in) :: left
    !> path to compare against
    type(tem_path_type), intent(in) :: right
    !> is equal??
    logical :: equality
    ! ---------------------------------------------------------------------------

    equality = ( tem_pathComparison( left, right ) .eq. 0 )

  end function isEqual
! ****************************************************************************** !


! ****************************************************************************** !
  !> This function provides the test for unequality of two path.
  function isUnequal(left, right) result(unequality)
    ! ---------------------------------------------------------------------------
    !> path to compare
    type(tem_path_type), intent(in) :: left
    !> path to compare against
    type(tem_path_type), intent(in) :: right
    !> is unequal??
    logical :: unequality
    ! ---------------------------------------------------------------------------

    unequality = ( tem_pathComparison( left, right ) .ne. 0 )

  end function isUnequal
! ****************************************************************************** !


! ****************************************************************************** !
  !> This function provides a comparison of two paths.
  function isSmaller(left, right) result(small)
    ! ---------------------------------------------------------------------------
    !> path to compare
    type(tem_path_type), intent(in) :: left
    !> path to compare against
    type(tem_path_type), intent(in) :: right
    !> is smaller??
    logical :: small
    ! ---------------------------------------------------------------------------

    small = ( tem_pathComparison( left, right ) .lt. 0 )

  end function isSmaller
! ****************************************************************************** !


! ****************************************************************************** !
  !> This function provides a comparison of two paths.
  function isSmallerOrEqual(left, right) result(small)
    ! ---------------------------------------------------------------------------
    !> path to compare
    type(tem_path_type), intent(in) :: left
    !> path to compare against
    type(tem_path_type), intent(in) :: right
    !> is smaller??
    logical :: small
    ! ---------------------------------------------------------------------------

    small = ( tem_pathComparison( left, right ) .le. 0 )

  end function isSmallerOrEqual
! ****************************************************************************** !


! ****************************************************************************** !
  !> This function provides a comparison of two paths.
  function isGreater(left, right) result(great)
    ! ---------------------------------------------------------------------------
    !> path to compare
    type(tem_path_type), intent(in) :: left
    !> path to compare against
    type(tem_path_type), intent(in) :: right
    !> is greater??
    logical :: great
    ! ---------------------------------------------------------------------------
    great = ( tem_pathComparison( left, right ) .gt. 0 )

  end function isGreater
! ****************************************************************************** !


! ****************************************************************************** !
  !> This function provides a comparison of two paths.
  function isGreaterOrEqual(left, right) result(great)
    ! ---------------------------------------------------------------------------
    !> path to compare
    type(tem_path_type), intent(in) :: left
    !> path to compare against
    type(tem_path_type), intent(in) :: right
    !> is greater??
    logical :: great
    ! ---------------------------------------------------------------------------

    great = ( tem_pathComparison( left, right ) .ge. 0 )

  end function isGreaterOrEqual
! ****************************************************************************** !


end module tem_path_array_module
! ****************************************************************************** !

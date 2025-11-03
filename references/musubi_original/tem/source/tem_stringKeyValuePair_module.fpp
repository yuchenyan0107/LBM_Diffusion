! Copyright (c) 2015-2016 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2016 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2016 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
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
?? include 'arrayMacros.inc'
module tem_stringKeyValuePair_module

  use env_module, only: labellen, minlength, zerolength

  implicit none

  private

  ! **************************************************************************** !
  !> Defines a key/value pair of strings that can be set or retrieved.
  type tem_stringKeyValuePair_type
    !> The key in the key/value pair. It's length is limited to /ref labellen.
    character(len=labellen) :: key
    !> The value in the key/value pair. It's length is limited to /ref labellen.
    character(len=labellen) :: value
  end type tem_stringKeyValuePair_type
  ! **************************************************************************** !

  interface operator (==)
    module procedure tem_stringKVP_equals
  end interface

  interface operator (/=)
    module procedure tem_stringKVP_notEquals
  end interface

  public :: tem_stringKeyValuePair_type
  public :: grw_stringKeyValuePairArray_type
  public :: append, truncate, init, destroy, empty, placeAt
  public :: operator(==)
  public :: operator(/=)

?? copy :: GA_decltxt(stringKeyValuePair, type(tem_stringKeyValuePair_type))


contains


  ! *****************************************************************************
  !> Indicates whether this instance and a specified object are equal.
  pure function tem_stringKVP_equals(me, other) result (res)
    ! ---------------------------------------------------------------------------
    !> THe current instance.
    type(tem_stringKeyValuePair_type), intent(in) :: me
    !> The instance to compare with the current instance.
    type(tem_stringKeyValuePair_type), intent(in) :: other
    !> .true. when both instances are equal, otherwise .false.
    logical :: res
    ! ---------------------------------------------------------------------------

    res = (me%key == other%key) .and. (me%value == other%value)

  end function tem_stringKVP_equals
  ! *****************************************************************************

  ! *****************************************************************************
  !> Indicates whether this instance and a specified object are equal.
  pure function tem_stringKVP_notEquals(me, other) result (res)
    ! ---------------------------------------------------------------------------
    !> THe current instance.
    type(tem_stringKeyValuePair_type), intent(in) :: me
    !> The instance to compare with the current instance.
    type(tem_stringKeyValuePair_type), intent(in) :: other
    !> .true. when both instances are equal, otherwise .false.
    logical :: res
    ! ---------------------------------------------------------------------------

    res = .not. tem_stringKVP_equals(me, other)

  end function tem_stringKVP_notEquals
  ! *****************************************************************************

?? copy :: GA_impltxt(stringKeyValuePair, type(tem_stringKeyValuePair_type))

end module tem_stringKeyValuePair_module

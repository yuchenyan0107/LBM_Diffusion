! Copyright (c) 2015 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
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
! ******************************************************************************!
module hvs_sizeof_module

  ! include treelm modules
  use env_module, only: rk, long_k

  use, intrinsic :: iso_c_binding, only: c_double, c_float, c_int, c_char,     &
    &                                    c_int_least8_t, c_int_least64_t

  implicit none

  private

  public :: c_sizeof

  interface c_sizeof
    module procedure sizeof_float
    module procedure sizeof_double
    module procedure sizeof_int
    module procedure sizeof_char
    module procedure sizeof_int_least8
    module procedure sizeof_int_least64
  end interface c_sizeof

contains

! ******************************************************************************!
  function sizeof_float(arg) result(argsize)
    ! ---------------------------------------------------------------------------
    real(kind=c_float), intent(in) :: arg
    integer :: argsize
    ! ---------------------------------------------------------------------------
    argsize = 4
  end function sizeof_float
! ******************************************************************************!


! ******************************************************************************!
  function sizeof_double(arg) result(argsize)
    ! ---------------------------------------------------------------------------
    real(kind=c_double), intent(in) :: arg
    integer :: argsize
    ! ---------------------------------------------------------------------------
    argsize = 8
  end function sizeof_double
! ******************************************************************************!


! ******************************************************************************!
  function sizeof_int(arg) result(argsize)
    ! ---------------------------------------------------------------------------
    integer(kind=c_int), intent(in) :: arg
    integer :: argsize
    ! ---------------------------------------------------------------------------
    argsize = 4
  end function sizeof_int
! ******************************************************************************!


! ******************************************************************************!
  function sizeof_int_least8(arg) result(argsize)
    ! ---------------------------------------------------------------------------
    integer(kind=c_int_least8_t), intent(in) :: arg
    integer :: argsize
    ! ---------------------------------------------------------------------------
    argsize = 1
  end function sizeof_int_least8
! ******************************************************************************!


! ******************************************************************************!
  function sizeof_int_least64(arg) result(argsize)
    ! ---------------------------------------------------------------------------
    integer(kind=c_int_least64_t), intent(in) :: arg
    integer :: argsize
    ! ---------------------------------------------------------------------------
    argsize = 8
  end function sizeof_int_least64
! ******************************************************************************!


! ******************************************************************************!
  function sizeof_char(arg) result(argsize)
    ! ---------------------------------------------------------------------------
    character(kind=c_char), intent(in) :: arg
    integer :: argsize
    ! ---------------------------------------------------------------------------
    argsize = 1
  end function sizeof_char
! ******************************************************************************!


end module hvs_sizeof_module
! ******************************************************************************!

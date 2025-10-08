! Copyright (c) 2013 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2013 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2013 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
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
!> Module to provide access to the gamma function.
!!
!! This dummy module is used, if the compiler does not
!! support this F2008 feature, and the NAG implementation
!! is not found.
!! We fallback to the gamma function provided by SLATEC.
!!
module tem_gamma_module

  ! include treelm modules
  use env_module, only: rk

  implicit none

  private

  logical, parameter, public :: elemental_gamma = .false.

  public :: gamma

contains


! ****************************************************************************** !
  !> Wrapper around the SLATEC gamma implementation.
  !!
  elemental function gamma(xx)
    ! ---------------------------------------------------------------------------
    real(kind=rk),intent(in) :: xx
    real(kind=rk) :: gamma
    ! ---------------------------------------------------------------------------
    double precision, external :: dgamma

    gamma = real(dgamma(xx), kind=rk)

  end function gamma
! ****************************************************************************** !


end module tem_gamma_module
! ****************************************************************************** !

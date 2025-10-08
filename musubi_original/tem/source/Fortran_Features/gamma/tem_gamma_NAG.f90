! Copyright (c) 2013 Harald Klimach <harald.klimach@uni-siegen.de>
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
!> Module to provide access to the gamma function.
!!
!! This module is used, to provide an adapter to the
!! Gamma function implementation by NAG.
!! It is only used if the NAG implementation is supported,
!! and the compiler itself lacks this F2008 feature.
!!
module tem_gamma_module

  ! include treelm modules
  use nag_library, only: s14aaf, nag_wp

  implicit none

  logical,parameter :: elemental_gamma = .false.

  contains

! ****************************************************************************** !
  !> @todo: Add description
  !!
  function gamma(xx)
    ! ---------------------------------------------------------------------------
    real(kind=nag_wp), intent(in) :: xx
    real(kind=nag_wp) :: gamma
    ! ---------------------------------------------------------------------------
    integer :: iFail
    ! ---------------------------------------------------------------------------

    gamma = s14aaf(xx, iFail)

  end function gamma
! ****************************************************************************** !


end module tem_gamma_module
! ****************************************************************************** !

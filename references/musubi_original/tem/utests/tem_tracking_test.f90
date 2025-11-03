! Copyright (c) 2013-2016 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2013, 2019 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2013 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2017 Peter Vitt <peter.vitt2@uni-siegen.de>
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
program tem_tracking_test
  use, intrinsic :: iso_c_binding, only: C_NEW_LINE
  use mpi
  ! incude treelm modules
  use env_module,             only: solSpecLen
  use tem_tracking_module,    only: tem_tracking_type, &
    &                               tem_load_tracking

  ! include aotus modules
  use aotus_module,           only: flu_State, close_config,      &
    &                               open_config_chunk

  !mpi!nprocs = 1

  implicit none

  character, parameter :: nl = C_NEW_LINE
  type(tem_tracking_type) :: track
  type(flu_State) :: conf
  integer :: iError

  character(len=solSpecLen), parameter :: sysConf =            &
    &    'tracking = {'                                  // nl &
    & // '  label = "point",'                            // nl &
    & // '  variable = { "density", "velocity" },'       // nl &
    & // '  folder = "tracking_",'                       // nl &
    & // '  shapce = { kind = "all" },'                  // nl &
    & // '  output = { format = "ascii" },'              // nl &
    & // '  time = { min = 1, max = 10, interval = 1 },' // nl &
    & // '}'

  call MPI_Init(iError)

  call open_config_chunk( L     = conf,         &
    &                     chunk = trim(sysConf) )

  ! load tracking
  call tem_load_tracking( me   = track, &
    &                     conf = conf   )

  call close_config( L = conf )
  write(*,*) 'PASSED'

  call MPI_Finalize(iError)

end program tem_tracking_test

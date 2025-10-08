! Copyright (c) 2012 Sathish Krishnan P S <s.krishnan@grs-sim.de>
! Copyright (c) 2012-2013 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2013-2014, 2016 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2016 Peter Vitt <peter.vitt2@uni-siegen.de>
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
!> Module to create Treelm environment for unit testcases.
module tem_utestEnv_module
  use, intrinsic :: iso_c_binding, only: C_NEW_LINE
  use tem_bc_prop_module,    only: tem_bc_prop_type, load_tem_bc_prop,         &
    &                              tem_empty_bc_prop
  use tem_property_module,   only: prp_hasbnd
  use treelmesh_module,      only: treelmesh_type, load_tem
  use tem_general_module,    only: tem_general_type, tem_load_general, tem_start
  use tem_logging_module,    only: tem_logging_load_primary

  use aotus_module,          only: open_config_chunk, close_config
  use flu_binding,           only: flu_state
  implicit none

  character, parameter :: nl = C_NEW_LINE

  character(len=400), parameter :: cubeconf =                                  &
    &   'mesh = {' // nl                                                       &
    & //'         predefined = "cube",' // nl                                  &
    & //'         origin = {0.0, 0.0, 0.0},' // nl                             &
    & //'         length=1.0,' // nl                                           &
    & //'         refinementLevel = 1' // nl                                   &
    & //'       }' // nl                                                       &
    & //'NOdebug = {logging = {' // nl                                         &
    & //'            level=10, '// nl                                          &
    & //'            filename="dbgutest",' // nl                               &
    & //'            root_only = false'// nl                                   &
    & //'          }}' // nl

contains

  !> Create the treelm environment for unit tests.
  subroutine load_env(tree, boundary, general )
    !---------------------------------------------------------------------------
    type(treelmesh_type), intent(out) :: tree
    type(tem_bc_prop_type), intent(out) :: boundary
    type(tem_general_type), intent(out) :: general
    !---------------------------------------------------------------------------
    type(flu_state) :: conf
    integer :: iProp
    !---------------------------------------------------------------------------

    ! Init the Treelm environment
    call tem_start('TREELM unit test', general)

    print *, "Hello from loadenv"
    ! Open the configuration file
    call open_config_chunk(L = conf, chunk = trim(cubeconf))

    ! load and initialize logUnit
    call tem_logging_load_primary(conf = conf,              &
      &                           rank = general%proc%rank  )

    call tem_load_general( me = general, conf = conf )

    ! Load the mesh first.
    call load_tem(me = tree, conf = conf, myPart = general%proc%rank, &
      &           nParts = general%proc%comm_size, comm = general%proc%comm )

    call close_config(conf)

    ! Load the properties now.
    do iprop = 1, size(tree%Property)
      ! Is the current property setting the prp_hasBnd bit?
      if (tree%global%property(iprop)%bitpos == prp_hasBnd) then
        ! Found the property setting the prp_hasBnd, load the
        ! Boundaries accordingly and leave the loop.
        call load_tem_bc_prop(me = boundary,                &
          &    offset = tree%Property(iprop)%Offset,        &
          &    nElems = tree%Property(iprop)%nElems,        &
          &    basename = trim(tree%global%dirname)//'bnd', &
          &    comm = general%proc%comm,                    &
          &    mypart = general%proc%rank )
        exit
      end if
    end do

  end subroutine load_env


end module tem_utestEnv_module

! Copyright (c) 2012 Sathish Krishnan P S <s.krishnan@grs-sim.de>
! Copyright (c) 2013-2014, 2019 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2013 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2013-2014, 2016 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2014 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
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
program tem_face_test

  ! include treelm modules
  use env_module,           only: long_k, fin_env
  use tem_utestEnv_module,  only: cubeconf
  use tem_comm_module,      only: tem_commPattern_type, tem_load_commPattern
  use tem_faceData_module,  only: tem_face_type
  use tem_face_module,      only: tem_build_face_info
  use treelmesh_module,     only: treelmesh_type
  use tem_bc_prop_module,   only: tem_bc_prop_type
  use tem_dyn_array_module, only: PositionOfVal
  use tem_aux_module,       only: tem_abort
  use tem_logging_module,   only: logUnit, tem_logging_load_primary
  use tem_general_module,   only: tem_general_type, tem_load_general, &
    &                             tem_start

  use aotus_module,          only: open_config_chunk
  use flu_binding,           only: flu_state

  !mpi!nprocs = 1

  implicit none
 
  write(*,*) 'Running tem_construction_test...'
  call check_serial_multilevel_faceDesc()
  call fin_env()


contains


  subroutine check_serial_multilevel_faceDesc()
    integer :: pos
    type(tem_general_type) :: general
    type(treelmesh_type) :: tree
    type(tem_bc_prop_type) :: boundary 
    type(tem_face_type),allocatable :: faces(:)
    type(flu_state) :: conf

    call tem_start('TREELM unit test', general)

    ! Open the configuration file 
    call open_config_chunk(L = conf, chunk = trim(cubeconf))
    ! load and initialize logUnit
    call tem_logging_load_primary(conf = conf,              &
      &                           rank = general%proc%rank  )
    ! Load Commpattern
    call tem_load_general( me = general, conf = conf )

    allocate(tree%treeID(15), tree%part_First(1), tree%part_Last(1), & 
             & tree%global%property(0), tree%property(0), tree%elempropertybits(15))
    tree%nElems = 15
    tree%treeID = (/1,2,3,33,34,35,36,37,38,39,40,5,6,7,8/)
    tree%part_First(1) = 1
    tree%part_Last(1) = 8 
    tree%global%minlevel = 1
    tree%global%maxlevel = 2
    tree%global%nElems = 15
    tree%global%nParts = 1
    tree%global%myPart = 0
    tree%elempropertybits(:) = 0
    tree%global%nProperties = 0

    tree%global%comm = general%proc%comm

    boundary%nBCTypes = 0

    write(logUnit(1),*) 'Building face info'

    ! Create face description
    call tem_build_face_info( tree = tree, boundary = boundary, &
                            & commPattern = general%commPattern, &
                            & proc = general%proc, &
                            & faces = faces, nEligibleChildren = 4 )

    ! Start our checks:
    
    !!!!!!!!!!!!!! Bottom Plane z=0 !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! We only have single level, so iterate over all these faces.
    ! First, we start with the faces aligned in x direction.

    print *, size(faces(1)%faces(1)%computeFace%leftPos)

    
    ! Let us have a look if the number of faces is 15.
    write(logUnit(1),*)'X-direction'

!    print *, "Level 1"
    if( size(faces(1)%faces(1)%computeFace%leftPos).ne.6 ) then
      write(logUnit(1),*)'ERROR in check_serial_multilevel_faceDesc : ' // &
                  & 'expected 6 faces, but face description is wrong ...'
      call tem_abort()
    else 
      write(logUnit(1),*)'SUCCESS in check_serial_multilevel_faceDesc : ' // &
                  & 'number of compute faces is 6.'
    end if
    print *, "Level 2"
    if( size(faces(2)%faces(1)%computeFace%leftPos).ne.12 ) then
      write(logUnit(1),*)'ERROR in check_serial_multilevel_faceDesc : ' // &
                  & 'expected 12 faces, but face description is wrong ...'
      call tem_abort()
    else 
      write(logUnit(1),*)'SUCCESS in check_serial_multilevel_faceDesc : ' // &
                  & 'number of compute faces is 12.'
    end if
     
    ! Let us have a look if the face 2-1 (left elem id - right elem id) is
    ! in the list of faces.
    pos = PositionOfVal( me = faces(1)%faces(1)%faceList%faceId, val = 2_long_k )
    if(pos .le. 0) then
      write(logUnit(1),*)'ERROR in check_serial_multilevel_faceDesc : ' // &
                  & 'not able to locate face 2-1 in face list ...'
      call tem_abort()
    else 
      write(logUnit(1),*)'SUCCESS in check_serial_multilevel_faceDesc : ' // &
                  & 'found face 2-1 in face list.'
      if(faces(1)%faces(1)%faceList%rightElemId%val(pos).ne.1_long_k) then
        write(logUnit(1),*)'ERROR in check_serial_multilevel_faceDesc : ' // &
                  & 'not able to locate face 2-1 in face list (wrong ' // &
                  & 'right elem id) ...'
        call tem_abort()
      else 
        write(logUnit(1),*)'SUCCESS in check_serial_multilevel_faceDesc : ' // &
                  & 'found face 2-1 in face list, correct right elem id.'
      end if
    end if

    ! Let us have a look if the face 5-6 (left elem id - right elem id) is
    ! in the list of faces.
    pos = PositionOfVal( me = faces(1)%faces(1)%faceList%faceId, val = 5_long_k )
    if(pos .le. 0) then
      write(logUnit(1),*)'ERROR in check_serial_multilevel_faceDesc : ' // &
                  & 'not able to locate face 5-6 in face list ...'
      call tem_abort()
    else 
      write(logUnit(1),*)'SUCCESS in check_serial_multilevel_faceDesc : ' // &
                  & 'found face 5-6 in face list.'
      if(faces(1)%faces(1)%faceList%rightElemId%val(pos).ne.6_long_k) then
        write(logUnit(1),*)'ERROR in check_serial_multilevel_faceDesc : ' // &
                  & 'not able to locate face 5-6 in face list (wrong ' // &
                  & 'right elem id) ...'
        call tem_abort()
      else 
        write(logUnit(1),*)'SUCCESS in check_serial_multilevel_faceDesc : ' // &
                  & 'found face 5-6 in face list, correct right elem id.'
      end if
    end if

    ! Let us have a look if the face 6-5 (left elem id - right elem id) is
    ! in the list of faces.
    pos = PositionOfVal( me = faces(1)%faces(1)%faceList%faceId, val = 6_long_k )
    if(pos .le. 0) then
      write(logUnit(1),*)'ERROR in check_serial_multilevel_faceDesc : ' // &
                  & 'not able to locate face 6-5 in face list ...'
      call tem_abort()
    else 
      write(logUnit(1),*)'SUCCESS in check_serial_multilevel_faceDesc : ' // &
                  & 'found face 6-5 in face list.'
      if(faces(1)%faces(1)%faceList%rightElemId%val(pos).ne.5_long_k) then
        write(logUnit(1),*)'ERROR in check_serial_multilevel_faceDesc : ' // &
                  & 'not able to locate face 6-5 in face list (wrong ' // &
                  & 'right elem id) ...'
        call tem_abort()
      else 
        write(logUnit(1),*)'SUCCESS in check_serial_multilevel_faceDesc : ' // &
                  & 'found face 6-5 in face list, correct right elem id.'
      end if
    end if

    ! Let us have a look if the face 7-8 (left elem id - right elem id) is
    ! in the list of faces.
    pos = PositionOfVal( me = faces(1)%faces(1)%faceList%faceId, val = 7_long_k )
    if(pos .le. 0) then
      write(logUnit(1),*)'ERROR in check_serial_multilevel_faceDesc : ' // &
                  & 'not able to locate face 7-8 in face list ...'
      call tem_abort()
    else 
      write(logUnit(1),*)'SUCCESS in check_serial_multilevel_faceDesc : ' // &
                  & 'found face 7-8 in face list.'
      if(faces(1)%faces(1)%faceList%rightElemId%val(pos).ne.8_long_k) then
        write(logUnit(1),*)'ERROR in check_serial_multilevel_faceDesc : ' // &
                  & 'not able to locate face 7-8 in face list (wrong ' // &
                  & 'right elem id) ...'
        call tem_abort()
      else 
        write(logUnit(1),*)'SUCCESS in check_serial_multilevel_faceDesc : ' // &
                  & 'found face 7-8 in face list, correct right elem id.'
      end if
    end if

    ! Let us have a look if the face 26-33 (left elem id - right elem id) is
    ! in the list of faces. Z=0. Level 2
    pos = PositionOfVal( me = faces(2)%faces(1)%faceList%faceId, val = 26_long_k )
    if(pos .le. 0) then
      write(logUnit(1),*)'ERROR in check_serial_multilevel_faceDesc : ' // &
                  & 'not able to locate face 26-33 in face list ...'
      call tem_abort()
    else 
      write(logUnit(1),*)'SUCCESS in check_serial_multilevel_faceDesc : ' // &
                  & 'found face 26-33 in face list.'
      if(faces(2)%faces(1)%faceList%rightElemId%val(pos).ne.33_long_k) then
        write(logUnit(1),*)'ERROR in check_serial_multilevel_faceDesc : ' // &
                  & 'not able to locate face 2-1 in face list (wrong ' // &
                  & 'right elem id) ...'
        call tem_abort()
      else 
        write(logUnit(1),*)'SUCCESS in check_serial_multilevel_faceDesc : ' // &
                  & 'found face 26-33 in face list, correct right elem id.'
      end if
    end if

    ! Let us have a look if the face 28-35 (left elem id - right elem id) is
    ! in the list of faces. Z=0. Level 2
    pos = PositionOfVal( me = faces(2)%faces(1)%faceList%faceId, val = 28_long_k )
    if(pos .le. 0) then
      write(logUnit(1),*)'ERROR in check_serial_multilevel_faceDesc : ' // &
                  & 'not able to locate face 28-35 in face list ...'
      call tem_abort()
    else 
      write(logUnit(1),*)'SUCCESS in check_serial_multilevel_faceDesc : ' // &
                  & 'found face 28-35 in face list.'
      if(faces(2)%faces(1)%faceList%rightElemId%val(pos).ne.35_long_k) then
        write(logUnit(1),*)'ERROR in check_serial_multilevel_faceDesc : ' // &
                  & 'not able to locate face 2-1 in face list (wrong ' // &
                  & 'right elem id) ...'
        call tem_abort()
      else 
        write(logUnit(1),*)'SUCCESS in check_serial_multilevel_faceDesc : ' // &
                  & 'found face 28-35 in face list, correct right elem id.'
      end if
    end if

    ! Let us have a look if the face 30-37 (left elem id - right elem id) is
    ! in the list of faces. Z=1/2
    pos = PositionOfVal( me = faces(2)%faces(1)%faceList%faceId, val = 30_long_k )
    if(pos .le. 0) then
      write(logUnit(1),*)'ERROR in check_serial_multilevel_faceDesc : ' // &
                  & 'not able to locate face 30-37 in face list ...'
      call tem_abort()
    else 
      write(logUnit(1),*)'SUCCESS in check_serial_multilevel_faceDesc : ' // &
                  & 'found face 30-37 in face list.'
      if(faces(2)%faces(1)%faceList%rightElemId%val(pos).ne.37_long_k) then
        write(logUnit(1),*)'ERROR in check_serial_multilevel_faceDesc : ' // &
                  & 'not able to locate face 2-1 in face list (wrong ' // &
                  & 'right elem id) ...'
        call tem_abort()
      else 
        write(logUnit(1),*)'SUCCESS in check_serial_multilevel_faceDesc : ' // &
                  & 'found face 30-37 in face list, correct right elem id.'
      end if
    end if


    ! Let us have a look if the face 32-39 (left elem id - right elem id) is
    ! in the list of faces. Z=1/2
    pos = PositionOfVal( me = faces(2)%faces(1)%faceList%faceId, val = 32_long_k )
    if(pos .le. 0) then
      write(logUnit(1),*)'ERROR in check_serial_multilevel_faceDesc : ' // &
                  & 'not able to locate face 32-39 in face list ...'
      call tem_abort()
    else 
      write(logUnit(1),*)'SUCCESS in check_serial_multilevel_faceDesc : ' // &
                  & 'found face 32-39 in face list.'
      if(faces(2)%faces(1)%faceList%rightElemId%val(pos).ne.39_long_k) then
        write(logUnit(1),*)'ERROR in check_serial_multilevel_faceDesc : ' // &
                  & 'not able to locate face 2-1 in face list (wrong ' // &
                  & 'right elem id) ...'
        call tem_abort()
      else 
        write(logUnit(1),*)'SUCCESS in check_serial_multilevel_faceDesc : ' // &
                  & 'found face 32-39 in face list, correct right elem id.'
      end if
    end if

    ! Let us have a look if the face 34-25 (left elem id - right elem id) is
    ! in the list of faces. Z=1/2
    pos = PositionOfVal( me = faces(2)%faces(1)%faceList%faceId, val = 34_long_k )
    if(pos .le. 0) then
      write(logUnit(1),*)'ERROR in check_serial_multilevel_faceDesc : ' // &
                  & 'not able to locate face 34-25 in face list ...'
      call tem_abort()
    else 
      write(logUnit(1),*)'SUCCESS in check_serial_multilevel_faceDesc : ' // &
                  & 'found face 34-25 in face list.'
      if(faces(2)%faces(1)%faceList%rightElemId%val(pos).ne.25_long_k) then
        write(logUnit(1),*)'ERROR in check_serial_multilevel_faceDesc : ' // &
                  & 'not able to locate face 2-1 in face list (wrong ' // &
                  & 'right elem id) ...'
        call tem_abort()
      else 
        write(logUnit(1),*)'SUCCESS in check_serial_multilevel_faceDesc : ' // &
                  & 'found face 34-25 in face list, correct right elem id.'
      end if
    end if


    ! Let us have a look if the face 36-27 (left elem id - right elem id) is
    ! in the list of faces. Z=1/2
    pos = PositionOfVal( me = faces(2)%faces(1)%faceList%faceId, val = 36_long_k )
    if(pos .le. 0) then
      write(logUnit(1),*)'ERROR in check_serial_multilevel_faceDesc : ' // &
                  & 'not able to locate face 36-27 in face list ...'
      call tem_abort()
    else 
      write(logUnit(1),*)'SUCCESS in check_serial_multilevel_faceDesc : ' // &
                  & 'found face 36-27 in face list.'
      if(faces(2)%faces(1)%faceList%rightElemId%val(pos).ne.27_long_k) then
        write(logUnit(1),*)'ERROR in check_serial_multilevel_faceDesc : ' // &
                  & 'not able to locate face 2-1 in face list (wrong ' // &
                  & 'right elem id) ...'
        call tem_abort()
      else 
        write(logUnit(1),*)'SUCCESS in check_serial_multilevel_faceDesc : ' // &
                  & 'found face 36-27 in face list, correct right elem id.'
      end if
    end if


    ! Let us have a look if the face 38-29 (left elem id - right elem id) is
    ! in the list of faces. Z=1/2
    pos = PositionOfVal( me = faces(2)%faces(1)%faceList%faceId, val = 38_long_k )
    if(pos .le. 0) then
      write(logUnit(1),*)'ERROR in check_serial_multilevel_faceDesc : ' // &
                  & 'not able to locate face 38-29 in face list ...'
      call tem_abort()
    else 
      write(logUnit(1),*)'SUCCESS in check_serial_multilevel_faceDesc : ' // &
                  & 'found face 38-29 in face list.'
      if(faces(2)%faces(1)%faceList%rightElemId%val(pos).ne.29_long_k) then
        write(logUnit(1),*)'ERROR in check_serial_multilevel_faceDesc : ' // &
                  & 'not able to locate face 2-1 in face list (wrong ' // &
                  & 'right elem id) ...'
        call tem_abort()
      else 
        write(logUnit(1),*)'SUCCESS in check_serial_multilevel_faceDesc : ' // &
                  & 'found face 38-29 in face list, correct right elem id.'
      end if
    end if


    ! Let us have a look if the face 40-31 (left elem id - right elem id) is
    ! in the list of faces. Z=1/2
    pos = PositionOfVal( me = faces(2)%faces(1)%faceList%faceId, val = 40_long_k )
    if(pos .le. 0) then
      write(logUnit(1),*)'ERROR in check_serial_multilevel_faceDesc : ' // &
                  & 'not able to locate face 40-31 in face list ...'
      call tem_abort()
    else 
      write(logUnit(1),*)'SUCCESS in check_serial_multilevel_faceDesc : ' // &
                  & 'found face 40-31 in face list.'
      if(faces(2)%faces(1)%faceList%rightElemId%val(pos).ne.31_long_k) then
        write(logUnit(1),*)'ERROR in check_serial_multilevel_faceDesc : ' // &
                  & 'not able to locate face 2-1 in face list (wrong ' // &
                  & 'right elem id) ...'
        call tem_abort()
      else 
        write(logUnit(1),*)'SUCCESS in check_serial_multilevel_faceDesc : ' // &
                  & 'found face 40-31 in face list, correct right elem id.'
      end if
    end if


    !! Y-direction
    print *, "Level 1"
    if( size(faces(1)%faces(2)%computeFace%leftPos).ne.6 ) then
      write(logUnit(1),*)'ERROR in check_serial_multilevel_faceDesc : ' // &
                  & 'expected 6 faces, but face description is wrong ...'
      call tem_abort()
    else 
      write(logUnit(1),*)'SUCCESS in check_serial_multilevel_faceDesc : ' // &
                  & 'number of compute faces is 6.'
    end if
    print *, "Level 2"
    if( size(faces(2)%faces(2)%computeFace%leftPos).ne.12 ) then
      write(logUnit(1),*)'ERROR in check_serial_multilevel_faceDesc : ' // &
                  & 'expected 12 faces, but face description is wrong ...'
      call tem_abort()
    else 
      write(logUnit(1),*)'SUCCESS in check_serial_multilevel_faceDesc : ' // &
                  & 'number of compute faces is 12.'
    end if



     
    ! Let us have a look if the face 3-1 (left elem id - right elem id) is
    ! in the list of faces.
    pos = PositionOfVal( me = faces(1)%faces(2)%faceList%faceId, val = 3_long_k )
    if(pos .le. 0) then
      write(logUnit(1),*)'ERROR in check_serial_multilevel_faceDesc : ' // &
                  & 'not able to locate face 3-1 in face list ...'
      call tem_abort()
    else 
      write(logUnit(1),*)'SUCCESS in check_serial_multilevel_faceDesc : ' // &
                  & 'found face 3-1 in face list.'

      if(faces(1)%faces(2)%faceList%rightElemId%val(pos).ne.1_long_k) then
        write(logUnit(1),*)'ERROR in check_serial_multilevel_faceDesc : ' // &
                  & 'not able to locate face 3-1 in face list (wrong ' // &
                  & 'right elem id) ...'
        call tem_abort()
      else 
        write(logUnit(1),*)'SUCCESS in check_serial_multilevel_faceDesc : ' // &
                  & 'found face 3-1 in face list, correct right elem id.'
      end if
    end if


    ! Let us have a look if the face 5-7 (left elem id - right elem id) is
    ! in the list of faces.
    pos = PositionOfVal( me = faces(1)%faces(2)%faceList%faceId, val = 5_long_k )
    if(pos .le. 0) then
      write(logUnit(1),*)'ERROR in check_serial_multilevel_faceDesc : ' // &
                  & 'not able to locate face 5-7 in face list ...'
      call tem_abort()
    else 
      write(logUnit(1),*)'SUCCESS in check_serial_multilevel_faceDesc : ' // &
                  & 'found face 5-7 in face list.'

      if(faces(1)%faces(2)%faceList%rightElemId%val(pos).ne.7_long_k) then
        write(logUnit(1),*)'ERROR in check_serial_multilevel_faceDesc : ' // &
                  & 'not able to locate face 5-7 in face list (wrong ' // &
                  & 'right elem id) ...'
        call tem_abort()
      else 
        write(logUnit(1),*)'SUCCESS in check_serial_multilevel_faceDesc : ' // &
                  & 'found face 5-7 in face list, correct right elem id.'
      end if
    end if
 
    ! Let us have a look if the face 7-5 (left elem id - right elem id) is
    ! in the list of faces.
    pos = PositionOfVal( me = faces(1)%faces(2)%faceList%faceId, val = 7_long_k )
    if(pos .le. 0) then
      write(logUnit(1),*)'ERROR in check_serial_multilevel_faceDesc : ' // &
                  & 'not able to locate face 7-5 in face list ...'
      call tem_abort()
    else 
      write(logUnit(1),*)'SUCCESS in check_serial_multilevel_faceDesc : ' // &
                  & 'found face 7-5 in face list.'

      if(faces(1)%faces(2)%faceList%rightElemId%val(pos).ne.5_long_k) then
        write(logUnit(1),*)'ERROR in check_serial_multilevel_faceDesc : ' // &
                  & 'not able to locate face 7-5 in face list (wrong ' // &
                  & 'right elem id) ...'
        call tem_abort()
      else 
        write(logUnit(1),*)'SUCCESS in check_serial_multilevel_faceDesc : ' // &
                  & 'found face 7-5 in face list, correct right elem id.'
      end if
    end if

    ! Let us have a look if the face 6-8 (left elem id - right elem id) is
    ! in the list of faces.
    pos = PositionOfVal( me = faces(1)%faces(2)%faceList%faceId, val = 6_long_k )
    if(pos .le. 0) then
      write(logUnit(1),*)'ERROR in check_serial_multilevel_faceDesc : ' // &
                  & 'not able to locate face 6-8 in face list ...'
      call tem_abort()
    else 
      write(logUnit(1),*)'SUCCESS in check_serial_multilevel_faceDesc : ' // &
                  & 'found face 6-8 in face list.'

      if(faces(1)%faces(2)%faceList%rightElemId%val(pos).ne.8_long_k) then
        write(logUnit(1),*)'ERROR in check_serial_multilevel_faceDesc : ' // &
                  & 'not able to locate face 6-8 in face list (wrong ' // &
                  & 'right elem id) ...'
        call tem_abort()
      else 
        write(logUnit(1),*)'SUCCESS in check_serial_multilevel_faceDesc : ' // &
                  & 'found face 6-8 in face list, correct right elem id.'
      end if
    end if

    ! Let us have a look if the face 19-33 (left elem id - right elem id) is
    ! in the list of faces.
    pos = PositionOfVal( me = faces(2)%faces(2)%faceList%faceId, val = 19_long_k )
    if(pos .le. 0) then
      write(logUnit(1),*)'ERROR in check_serial_multilevel_faceDesc : ' // &
                  & 'not able to locate face 19-33 in face list ...'
      call tem_abort()
    else 
      write(logUnit(1),*)'SUCCESS in check_serial_multilevel_faceDesc : ' // &
                  & 'found face 19-33 in face list.'

      if(faces(2)%faces(2)%faceList%rightElemId%val(pos).ne.33_long_k) then
        write(logUnit(1),*)'ERROR in check_serial_multilevel_faceDesc : ' // &
                  & 'not able to locate face 19-33 in face list (wrong ' // &
                  & 'right elem id) ...'
        call tem_abort()
      else 
        write(logUnit(1),*)'SUCCESS in check_serial_multilevel_faceDesc : ' // &
                  & 'found face 19-33 in face list, correct right elem id.'
      end if
     end if

    ! Let us have a look if the face 20-34 (left elem id - right elem id) is
    ! in the list of faces.
    pos = PositionOfVal( me = faces(2)%faces(2)%faceList%faceId, val = 20_long_k )
    if(pos .le. 0) then
      write(logUnit(1),*)'ERROR in check_serial_multilevel_faceDesc : ' // &
                  & 'not able to locate face 20-34 in face list ...'
      call tem_abort()
    else 
      write(logUnit(1),*)'SUCCESS in check_serial_multilevel_faceDesc : ' // &
                  & 'found face 20-34 in face list.'

      if(faces(2)%faces(2)%faceList%rightElemId%val(pos).ne.34_long_k) then
        write(logUnit(1),*)'ERROR in check_serial_multilevel_faceDesc : ' // &
                  & 'not able to locate face 20-34 in face list (wrong ' // &
                  & 'right elem id) ...'
        call tem_abort()
      else 
        write(logUnit(1),*)'SUCCESS in check_serial_multilevel_faceDesc : ' // &
                  & 'found face 20-34 in face list, correct right elem id.'
      end if
     end if

    ! Let us have a look if the face 35-17 (left elem id - right elem id) is
    ! in the list of faces.
    pos = PositionOfVal( me = faces(2)%faces(2)%faceList%faceId, val = 35_long_k )
    if(pos .le. 0) then
      write(logUnit(1),*)'ERROR in check_serial_multilevel_faceDesc : ' // &
                  & 'not able to locate face 35-17 in face list ...'
      call tem_abort()
    else 
      write(logUnit(1),*)'SUCCESS in check_serial_multilevel_faceDesc : ' // &
                  & 'found face 35-17 in face list.'

      if(faces(2)%faces(2)%faceList%rightElemId%val(pos).ne.17_long_k) then
        write(logUnit(1),*)'ERROR in check_serial_multilevel_faceDesc : ' // &
                  & 'not able to locate face 35-17 in face list (wrong ' // &
                  & 'right elem id) ...'
        call tem_abort()
      else 
        write(logUnit(1),*)'SUCCESS in check_serial_multilevel_faceDesc : ' // &
                  & 'found face 35-17 in face list, correct right elem id.'
      end if
     end if



    !! Z-direction
    print *, "Level 1"
    if( size(faces(1)%faces(3)%computeFace%leftPos).ne.6 ) then
      write(logUnit(1),*)'ERROR in check_serial_multilevel_faceDesc : ' // &
                  & 'expected 6 faces, but face description is wrong ...'
      call tem_abort()
    else 
      write(logUnit(1),*)'SUCCESS in check_serial_multilevel_faceDesc : ' // &
                  & 'number of compute faces is 6.'
    end if
    print *, "Level 2"
    if( size(faces(2)%faces(3)%computeFace%leftPos).ne.12 ) then
      write(logUnit(1),*)'ERROR in check_serial_multilevel_faceDesc : ' // &
                  & 'expected 12 faces, but face description is wrong ...'
      call tem_abort()
    else 
      write(logUnit(1),*)'SUCCESS in check_serial_multilevel_faceDesc : ' // &
                  & 'number of compute faces is 12.'
    end if



 
    ! Let us have a look if the face 1-5 (left elem id - right elem id) is
    ! in the list of faces.
    pos = PositionOfVal( me = faces(1)%faces(3)%faceList%faceId, val = 1_long_k )
    if(pos .le. 0) then
      write(logUnit(1),*)'ERROR in check_serial_multilevel_faceDesc : ' // &
                  & 'not able to locate face 1-5 in face list ...'
      call tem_abort()
    else 
      write(logUnit(1),*)'SUCCESS in check_serial_multilevel_faceDesc : ' // &
                  & 'found face 1-5 in face list.'

      if(faces(1)%faces(3)%faceList%rightElemId%val(pos).ne.5_long_k) then
        write(logUnit(1),*)'ERROR in check_serial_multilevel_faceDesc : ' // &
                  & 'not able to locate face 1-5 in face list (wrong ' // &
                  & 'right elem id) ...'
        call tem_abort()
      else 
        write(logUnit(1),*)'SUCCESS in check_serial_multilevel_faceDesc : ' // &
                  & 'found face 1-5 in face list, correct right elem id.'
      end if
    end if


    ! Let us have a look if the face 2-6 (left elem id - right elem id) is
    ! in the list of faces.
    pos = PositionOfVal( me = faces(1)%faces(3)%faceList%faceId, val = 2_long_k )
    if(pos .le. 0) then
      write(logUnit(1),*)'ERROR in check_serial_multilevel_faceDesc : ' // &
                  & 'not able to locate face 2-6 in face list ...'
      call tem_abort()
    else 
      write(logUnit(1),*)'SUCCESS in check_serial_multilevel_faceDesc : ' // &
                  & 'found face 2-6 in face list.'

      if(faces(1)%faces(3)%faceList%rightElemId%val(pos).ne.6_long_k) then
        write(logUnit(1),*)'ERROR in check_serial_multilevel_faceDesc : ' // &
                  & 'not able to locate face 2-6 in face list (wrong ' // &
                  & 'right elem id) ...'
        call tem_abort()
      else 
        write(logUnit(1),*)'SUCCESS in check_serial_multilevel_faceDesc : ' // &
                  & 'found face 2-6 in face list, correct right elem id.'
      end if
    end if

    ! Let us have a look if the face 3-7 (left elem id - right elem id) is
    ! in the list of faces.
    pos = PositionOfVal( me = faces(1)%faces(3)%faceList%faceId, val = 3_long_k )
    if(pos .le. 0) then
      write(logUnit(1),*)'ERROR in check_serial_multilevel_faceDesc : ' // &
                  & 'not able to locate face 3-7 in face list ...'
      call tem_abort()
    else 
      write(logUnit(1),*)'SUCCESS in check_serial_multilevel_faceDesc : ' // &
                  & 'found face 3-7 in face list.'

      if(faces(1)%faces(3)%faceList%rightElemId%val(pos).ne.7_long_k) then
        write(logUnit(1),*)'ERROR in check_serial_multilevel_faceDesc : ' // &
                  & 'not able to locate face 3-7 in face list (wrong ' // &
                  & 'right elem id) ...'
        call tem_abort()
      else 
        write(logUnit(1),*)'SUCCESS in check_serial_multilevel_faceDesc : ' // &
                  & 'found face 3-7 in face list, correct right elem id.'
      end if
    end if


    ! Let us have a look if the face 69-33 (left elem id - right elem id) is
    ! in the list of faces.
    pos = PositionOfVal( me = faces(2)%faces(3)%faceList%faceId, val = 69_long_k )
    if(pos .le. 0) then
      write(logUnit(1),*)'ERROR in check_serial_multilevel_faceDesc : ' // &
                  & 'not able to locate face 69-33 in face list ...'
      call tem_abort()
    else 
      write(logUnit(1),*)'SUCCESS in check_serial_multilevel_faceDesc : ' // &
                  & 'found face 69-33 in face list.'

      if(faces(2)%faces(3)%faceList%rightElemId%val(pos).ne.33_long_k) then
        write(logUnit(1),*)'ERROR in check_serial_multilevel_faceDesc : ' // &
                  & 'not able to locate face 69-33 in face list (wrong ' // &
                  & 'right elem id) ...'
        call tem_abort()
      else 
        write(logUnit(1),*)'SUCCESS in check_serial_multilevel_faceDesc : ' // &
                  & 'found face 69-33 in face list, correct right elem id.'
      end if
    end if


    ! Let us have a look if the face 70-34 (left elem id - right elem id) is
    ! in the list of faces.
    pos = PositionOfVal( me = faces(2)%faces(3)%faceList%faceId, val = 70_long_k )
    if(pos .le. 0) then
      write(logUnit(1),*)'ERROR in check_serial_multilevel_faceDesc : ' // &
                  & 'not able to locate face 70-34 in face list ...'
      call tem_abort()
    else 
      write(logUnit(1),*)'SUCCESS in check_serial_multilevel_faceDesc : ' // &
                  & 'found face 70-34 in face list.'

      if(faces(2)%faces(3)%faceList%rightElemId%val(pos).ne.34_long_k) then
        write(logUnit(1),*)'ERROR in check_serial_multilevel_faceDesc : ' // &
                  & 'not able to locate face 70-34 in face list (wrong ' // &
                  & 'right elem id) ...'
        call tem_abort()
      else 
        write(logUnit(1),*)'SUCCESS in check_serial_multilevel_faceDesc : ' // &
                  & 'found face 70-34 in face list, correct right elem id.'
      end if
    end if


    write(logUnit(1),*)'PASSED'

  end subroutine check_serial_multilevel_faceDesc

end program tem_face_test

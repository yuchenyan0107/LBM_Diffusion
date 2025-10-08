! Copyright (c) 2012 Sathish Krishnan P S <s.krishnan@grs-sim.de>
! Copyright (c) 2012-2014, 2019 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2013 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2013-2014 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2014 Jens Zudrop <j.zudrop@grs-sim.de>
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
  use tem_utestEnv_module,  only: load_env
  use tem_faceData_module,  only: tem_face_type
  use tem_face_module,      only: tem_build_face_info
  use treelmesh_module,     only: treelmesh_type 
  use tem_bc_prop_module,   only: tem_bc_prop_type
  use env_module,           only: long_k
  use tem_dyn_array_module, only: PositionOfVal 
  use tem_aux_module,       only: tem_abort
  use tem_general_module,   only: tem_general_type, tem_finalize
  use tem_logging_module,   only: logUnit

!use tem_element_module, only: PositionOfVal

  !mpi!nprocs = 1

  implicit none
 
  write(*,*) 'Running tem_construction_test...'
  call check_serial_singlelevel_faceDesc()

contains

  subroutine check_serial_singlelevel_faceDesc()
!    character(len=*) :: filename
    integer i
    integer :: pos
    type(treelmesh_type) :: tree
    type(tem_bc_prop_type) :: boundary 
    type(tem_general_type) :: general
    type(tem_face_type), allocatable :: faces(:)
    integer :: level
 !   type(tem_face_type) :: faces(tree%global%minLevel:tree%global%maxLevel)
      

    ! Load the environment.
    call load_env(tree = tree, boundary = boundary, general = general ) 
    level = tree%global%minLevel

    write(*,*) "Building face info"
    ! Create face description
    call tem_build_face_info( tree = tree, boundary = boundary, &
                            & commPattern = general%commPattern, &
                            & proc = general%proc, &
                            & faces = faces, nEligibleChildren = 4 )

    write(*,*) tree%nElems

    do i=1, tree%nElems
       write(*,*) tree%treeID(i)
    end do

    write(*,*) tree%part_Last(1)
    write(*,*) tree%elemOffset


    ! Start our checks:
    
    !!!!!!!!!!!!!! Bottom Plane z=0 !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! We only have single level, so iterate over all these faces.
    ! First, we start with the faces aligned in x direction.
    
    ! Let us have a look if the number of faces is 8.
    write(logUnit(1),*)'X-direction!'
    if( size(faces(level)%faces(1)%computeFace%leftPos).ne.8 ) then
      write(logUnit(1),*)'ERROR in check_serial_singlelevel_faceDesc : ' // &
                  & 'expected 8 faces, but face description is wrong ...'
      call tem_abort()
    else 
      write(logUnit(1),*)'SUCCESS in check_serial_singlelevel_faceDesc : ' // &
                  & 'number of compute faces is 8.'
    end if
     
    ! Let us have a look if the face 2-1 (left elem id - right elem id) is
    ! in the list of faces.
    pos = PositionOfVal( me = faces(level)%faces(1)%faceList%faceId, val = 2_long_k )
    if(pos .le. 0) then
      write(logUnit(1),*)'ERROR in check_serial_singlelevel_faceDesc : ' // &
                  & 'not able to locate face 2-1 in face list ...'
      call tem_abort()
    else 
      write(logUnit(1),*)'SUCCESS in check_serial_singlelevel_faceDesc : ' // &
                  & 'found face 2-1 in face list.'
      if(faces(level)%faces(1)%faceList%rightElemId%val(pos).ne.1_long_k) then
        write(logUnit(1),*)'ERROR in check_serial_singlelevel_faceDesc : ' // &
                  & 'not able to locate face 2-1 in face list (wrong ' // &
                  & 'right elem id) ...'
        call tem_abort()
      else 
        write(logUnit(1),*)'SUCCESS in check_serial_singlelevel_faceDesc : ' // &
                  & 'found face 2-1 in face list, correct right elem id.'
      end if
    end if


    ! Checking in y-direction    
    ! Let us have a look if the number of faces is 8.
    write(logUnit(1),*)'Y-direction'
    if( size(faces(level)%faces(2)%computeFace%leftPos).ne.8 ) then
      write(logUnit(1),*)'ERROR in check_serial_singlelevel_faceDesc : ' // &
                  & 'expected 8 faces, but face description is wrong ...'
      call tem_abort()
    else 
      write(logUnit(1),*)'SUCCESS in check_serial_singlelevel_faceDesc : ' // &
                  & 'number of compute faces is 8.'
    end if
     
    ! Let us have a look if the face 3-1 (left elem id - right elem id) is
    ! in the list of faces.
    pos = PositionOfVal( me = faces(level)%faces(2)%faceList%faceId, val = 3_long_k )
    if(pos .le. 0) then
      write(logUnit(1),*)'ERROR in check_serial_singlelevel_faceDesc : ' // &
                  & 'not able to locate face 3-1 in face list ...'
      call tem_abort()
    else 
      write(logUnit(1),*)'SUCCESS in check_serial_singlelevel_faceDesc : ' // &
                  & 'found face 3-1 in face list.'

      if(faces(level)%faces(2)%faceList%rightElemId%val(pos).ne.1_long_k) then
        write(logUnit(1),*)'ERROR in check_serial_singlelevel_faceDesc : ' // &
                  & 'not able to locate face 3-1 in face list (wrong ' // &
                  & 'right elem id) ...'
        call tem_abort()
      else 
        write(logUnit(1),*)'SUCCESS in check_serial_singlelevel_faceDesc : ' // &
                  & 'found face 3-1 in face list, correct right elem id.'
      end if
    end if




    ! Checking in z-direction    
    ! Let us have a look if the number of faces is 8.
    write(logUnit(1),*)'Z-direction'
    if( size(faces(level)%faces(3)%computeFace%leftPos).ne.8 ) then
      write(logUnit(1),*)'ERROR in check_serial_singlelevel_faceDesc : ' // &
                  & 'expected 8 faces, but face description is wrong ...'
      call tem_abort()
    else 
      write(logUnit(1),*)'SUCCESS in check_serial_singlelevel_faceDesc : ' // &
                  & 'number of compute faces is 8.'
    end if
     
    ! Let us have a look if the face 1-5 (left elem id - right elem id) is
    ! in the list of faces.
    pos = PositionOfVal( me = faces(level)%faces(3)%faceList%faceId, val = 1_long_k )
    if(pos .le. 0) then
      write(logUnit(1),*)'ERROR in check_serial_singlelevel_faceDesc : ' // &
                  & 'not able to locate face 1-5 in face list ...'
      call tem_abort()
    else 
      write(logUnit(1),*)'SUCCESS in check_serial_singlelevel_faceDesc : ' // &
                  & 'found face 1-5 in face list.'

      if(faces(level)%faces(3)%faceList%rightElemId%val(pos).ne.5_long_k) then
        write(logUnit(1),*)'ERROR in check_serial_singlelevel_faceDesc : ' // &
                  & 'not able to locate face 1-5 in face list (wrong ' // &
                  & 'right elem id) ...'
        call tem_abort()
      else 
        write(logUnit(1),*)'SUCCESS in check_serial_singlelevel_faceDesc : ' // &
                  & 'found face 1-5 in face list, correct right elem id.'
      end if
    end if



    ! Checking in z-direction    
    ! Let us have a look if the face 4-8 (left elem id - right elem id) is
    ! in the list of faces.
    pos = PositionOfVal( me = faces(level)%faces(3)%faceList%faceId, val = 4_long_k )
    if(pos .le. 0) then
      write(logUnit(1),*)'ERROR in check_serial_singlelevel_faceDesc : ' // &
                  & 'not able to locate face 4-8 in face list ...'
      call tem_abort()
    else 
      write(logUnit(1),*)'SUCCESS in check_serial_singlelevel_faceDesc : ' // &
                  & 'found face 4-8 in face list.'

      if(faces(level)%faces(3)%faceList%rightElemId%val(pos).ne.8_long_k) then
        write(logUnit(1),*)'ERROR in check_serial_singlelevel_faceDesc : ' // &
                  & 'not able to locate face 4-8 in face list (wrong ' // &
                  & 'right elem id) ...'
        call tem_abort()
      else 
        write(logUnit(1),*)'SUCCESS in check_serial_singlelevel_faceDesc : ' // &
                  & 'found face 4-8 in face list, correct right elem id.'
      end if
    end if


    !!!!!!!!!!!!!!!!! Top plane Z=1 !!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! We only have single level, so iterate over all these faces.
    ! First, we start with the faces aligned in x direction.
    write(logUnit(1),*)'X-direction'
     
    ! Let us have a look if the face 6-5 (left elem id - right elem id) is
    ! in the list of faces.
    pos = PositionOfVal( me = faces(level)%faces(1)%faceList%faceId, val = 6_long_k )
    if(pos .le. 0) then
      write(logUnit(1),*)'ERROR in check_serial_singlelevel_faceDesc : ' // &
                  & 'not able to locate face 6-5 in face list ...'
      call tem_abort()
    else 
      write(logUnit(1),*)'SUCCESS in check_serial_singlelevel_faceDesc : ' // &
                  & 'found face 6-5 in face list.'
      if(faces(level)%faces(1)%faceList%rightElemId%val(pos).ne.5_long_k) then
        write(logUnit(1),*)'ERROR in check_serial_singlelevel_faceDesc : ' // &
                  & 'not able to locate face 6-5 in face list (wrong ' // &
                  & 'right elem id) ...'
        call tem_abort()
      else 
        write(logUnit(1),*)'SUCCESS in check_serial_singlelevel_faceDesc : ' // &
                  & 'found face 6-5 in face list, correct right elem id.'
      end if
    end if


    ! Checking in y-direction    
    ! Let us have a look if the number of faces is 8.
    write(logUnit(1),*)'Y-direction'
    if( size(faces(level)%faces(2)%computeFace%leftPos).ne.8 ) then
      write(logUnit(1),*)'ERROR in check_serial_singlelevel_faceDesc : ' // &
                  & 'expected 8 faces, but face description is wrong ...'
      call tem_abort()
    else 
      write(logUnit(1),*)'SUCCESS in check_serial_singlelevel_faceDesc : ' // &
                  & 'number of compute faces is 8.'
    end if
     
    ! Let us have a look if the face 7-5 (left elem id - right elem id) is
    ! in the list of faces.
    pos = PositionOfVal( me = faces(level)%faces(2)%faceList%faceId, val = 7_long_k )
    if(pos .le. 0) then
      write(logUnit(1),*)'ERROR in check_serial_singlelevel_faceDesc : ' // &
                  & 'not able to locate face 7-5 in face list ...'
      call tem_abort()
    else 
      write(logUnit(1),*)'SUCCESS in check_serial_singlelevel_faceDesc : ' // &
                  & 'found face 7-5 in face list.'

      if(faces(level)%faces(2)%faceList%rightElemId%val(pos).ne.5_long_k) then
        write(logUnit(1),*)'ERROR in check_serial_singlelevel_faceDesc : ' // &
                  & 'not able to locate face 7-5 in face list (wrong ' // &
                  & 'right elem id) ...'
        call tem_abort()
      else 
        write(logUnit(1),*)'SUCCESS in check_serial_singlelevel_faceDesc : ' // &
                  & 'found face 7-5 in face list, correct right elem id.'
      end if
    end if

    call tem_finalize(general)

    ! When we reach this point, all tests have passed...
    write(logUnit(1),*)'PASSED'

  end subroutine check_serial_singlelevel_faceDesc

end program tem_face_test

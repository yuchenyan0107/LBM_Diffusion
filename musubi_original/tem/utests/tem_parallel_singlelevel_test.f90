! Copyright (c) 2013 Sathish Krishnan P S <s.krishnan@grs-sim.de>
! Copyright (c) 2013-2014 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2013 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2013 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2014 Jens Zudrop <j.zudrop@grs-sim.de>
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
program tem_face_test
  use tem_utestEnv_module,  only: load_env
  use tem_faceData_module,  only: tem_face_type
  use tem_face_module,      only: tem_build_face_info
  use tem_comm_env_module,  only: tem_comm_env_type
  use treelmesh_module,     only: treelmesh_type 
  use tem_bc_prop_module,   only: tem_bc_prop_type
  use env_module,           only: long_k
  use tem_dyn_array_module, only: PositionOfVal
  use tem_aux_module,       only: tem_abort
  use env_module,           only: fin_env
  use tem_general_module,   only: tem_general_type
  use tem_logging_module,   only: logUnit

  !mpi!nprocs = 2 


  implicit none
 
  write(*,*) 'Running tem_construction_test...'
  call check_parallel_singlelevel_faceDesc()

  call fin_env()

  write(*,*) 'PASSED'

contains

  subroutine check_parallel_singlelevel_faceDesc()
    integer :: pos
    type(tem_general_type) :: general
    type(treelmesh_type) :: tree
    type(tem_bc_prop_type) :: boundary 
    type(tem_face_type), allocatable :: faces(:)
    integer :: level
      
    ! Load the environment.
    call load_env(tree = tree, boundary = boundary, general = general ) 
    level = tree%global%minLevel

    ! Create face description
    write(*,*) 'Building face info'
    call tem_build_face_info( tree = tree, boundary = boundary, &
                            & commPattern = general%commPattern, &
                            & proc = general%proc, &
                            & faces = faces, nEligibleChildren = 4 )

    write(*,*) 'Rank',  general%proc%rank

   ! Start our checks:
  
   !!!!!!!!!!!!!! Bottom Plane z=0 !!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! process_0 contains plane z=0  for nprocs=2
   ! Check faces in x and y-direction 
   if(general%proc%rank.eq.0)  then
      ! Let us have a look if the number of faces is 4.
      write(*,*) 'Process : ', general%proc%rank
      write(*,*) "X-direction!"
      if( faces(level)%faces(1)%faceList%faceId%nvals.ne.4 ) then
        write(*,*) 'ERROR in check_parallel_singlelevel_faceDesc : //&
                    &expected 4 faces, but face description is wrong ...'
        call tem_abort()
      else 
        write(*,*) 'SUCCESS in check_parallel_singlelevel_faceDesc : ' // &
                    & 'number of nvals is 4.'
      end if
       
      ! Let us have a look if the face 2-1 (left elem id - right elem id) is
      ! in the list of faces.
      pos = PositionOfVal( me = faces(level)%faces(1)%faceList%faceId, val = 2_long_k )
      if(pos .le. 0) then
        write(*,*) 'ERROR in check_parallel_singlelevel_faceDesc : ' // &
                    & 'not able to locate face 2-1 in face list ...'
        call tem_abort()
      else 
        write(*,*) 'SUCCESS in check_parallel_singlelevel_faceDesc : ' // &
                    & 'found face 2-1 in face list.'
        if(faces(level)%faces(1)%faceList%rightElemId%val(pos).ne.1_long_k) then
          write(*,*) 'ERROR in check_parallel_singlelevel_faceDesc : ' // &
                    & 'not able to locate face 2-1 in face list (wrong ' // &
                    & 'right elem id) ...'
          call tem_abort()
        else 
          write(*,*) 'SUCCESS in check_parallel_singlelevel_faceDesc : ' // &
                    & 'found face 2-1 in face list, correct right elem id.'
        end if
      end if

      ! Checking in y-direction    
      ! Let us have a look if the number of faces is 4.
      write(*,*) 'Y-direction'
      if(  faces(level)%faces(2)%faceList%faceId%nvals.ne.4 ) then
        write(*,*) 'ERROR in check_parallel_singlelevel_faceDesc : ' // &
                    & 'expected 4 faces, but face description is wrong ...'
        call tem_abort()
      else 
        write(*,*) 'SUCCESS in check_parallel_singlelevel_faceDesc : ' // &
                    & 'number of nvals is 4.'
      end if
       
      ! Let us have a look if the face 3-1 (left elem id - right elem id) is
      ! in the list of faces.
      pos = PositionOfVal( me = faces(level)%faces(2)%faceList%faceId, val = 3_long_k )
      if(pos .le. 0) then
        write(*,*) 'ERROR in check_parallel_singlelevel_faceDesc : ' // &
                    & 'not able to locate face 3-1 in face list ...'
        call tem_abort()
      else 
        write(*,*) 'SUCCESS in check_parallel_singlelevel_faceDesc : ' // &
                    & 'found face 3-1 in face list.'
  
        if(faces(level)%faces(2)%faceList%rightElemId%val(pos).ne.1_long_k) then
          write(*,*) 'ERROR in check_parallel_singlelevel_faceDesc : ' // &
                    & 'not able to locate face 3-1 in face list (wrong ' // &
                    & 'right elem id) ...'
          call tem_abort()
        else 
          write(*,*) 'SUCCESS in check_parallel_singlelevel_faceDesc : ' // &
                    & 'found face 3-1 in face list, correct right elem id.'
        end if
      end if
    end if



    ! Checking in z-direction    
    ! Let us have a look if the number of faces is 8.
    write(logUnit(1),*)'Z-direction'
    if( faces(level)%faces(3)%faceList%faceId%nvals.ne.8 ) then
      write(logUnit(1),*)'ERROR in check_parallel_singlelevel_faceDesc : ' // &
                  & 'expected 8 faces, but face description is wrong ...'
      call tem_abort()
    else 
      write(logUnit(1),*)'SUCCESS in check_parallel_singlelevel_faceDesc : ' // &
                  & 'number of nvals is 8.'
    end if
     
    ! Let us have a look if the face 1-5 (left elem id - right elem id) is
    ! in the list of faces.
    pos = PositionOfVal( me = faces(level)%faces(3)%faceList%faceId, val = 1_long_k )
    if(pos .le. 0) then
      write(logUnit(1),*)'ERROR in check_parallel_singlelevel_faceDesc : ' // &
                  & 'not able to locate face 1-5 in face list ...'
      call tem_abort()
    else 
      write(logUnit(1),*)'SUCCESS in check_parallel_singlelevel_faceDesc : ' // &
                  & 'found face 1-5 in face list.'

      if(faces(level)%faces(3)%faceList%rightElemId%val(pos).ne.5_long_k) then
        write(logUnit(1),*)'ERROR in check_parallel_singlelevel_faceDesc : ' // &
                  & 'not able to locate face 1-5 in face list (wrong ' // &
                  & 'right elem id) ...'
        call tem_abort()
      else 
        write(logUnit(1),*)'SUCCESS in check_parallel_singlelevel_faceDesc : ' // &
                  & 'found face 1-5 in face list, correct right elem id.'
      end if
    end if

    ! Let us have a look if the face 4-8 (left elem id - right elem id) is
    ! in the list of faces.
    pos = PositionOfVal( me = faces(level)%faces(3)%faceList%faceId, val = 4_long_k )
    if(pos .le. 0) then
      write(logUnit(1),*)'ERROR in check_parallel_singlelevel_faceDesc : ' // &
                  & 'not able to locate face 4-8 in face list ...'
      call tem_abort()
    else 
      write(logUnit(1),*)'SUCCESS in check_parallel_singlelevel_faceDesc : ' // &
                  & 'found face 4-8 in face list.'

      if(faces(level)%faces(3)%faceList%rightElemId%val(pos).ne.8_long_k) then
        write(logUnit(1),*)'ERROR in check_parallel_singlelevel_faceDesc : ' // &
                  & 'not able to locate face 4-8 in face list (wrong ' // &
                  & 'right elem id) ...'
        call tem_abort()
      else 
        write(logUnit(1),*)'SUCCESS in check_parallel_singlelevel_faceDesc : ' // &
                  & 'found face 4-8 in face list, correct right elem id.'
      end if
    end if


   !!!!!!!!!!!!!!!!! Top plane Z=1 !!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! process_1 contains the plane Z=1 for nprocs=2
   ! Check faces in x and y-direction
    if(general%proc%rank.eq.1) then
      ! x-direction
      write(*,*) 'Process :', general%proc%rank
      write(*,*) 'X-direction'
      if( faces(level)%faces(1)%faceList%faceId%nvals.ne.4 ) then
        write(*,*) 'ERROR in check_parallel_singlelevel_faceDesc : ' // &
                    & 'expected 4 faces, but face description is wrong ...'
        call tem_abort()
      else 
        write(*,*) 'SUCCESS in check_parallel_singlelevel_faceDesc : ' // &
                    & 'number of nvals is 4.'
      end if

      ! Let us have a look if the face 6-5 (left elem id - right elem id) is
      ! in the list of faces.
      pos = PositionOfVal( me = faces(level)%faces(1)%faceList%faceId, val = 6_long_k )
      if(pos .le. 0) then
        write(*,*) 'ERROR in check_parallel_singlelevel_faceDesc : ' // &
                    & 'not able to locate face 6-5 in face list ...'
        call tem_abort()
      else 
        write(*,*) 'SUCCESS in check_parallel_singlelevel_faceDesc : ' // &
                    & 'found face 6-5 in face list.'
        if(faces(level)%faces(1)%faceList%rightElemId%val(pos).ne.5_long_k) then
          write(*,*) 'ERROR in check_parallel_singlelevel_faceDesc : ' // &
                     & 'not able to locate face 6-5 in face list (wrong ' // &
                     & 'right elem id) ...'
          call tem_abort()
        else 
          write(*,*) 'SUCCESS in check_parallel_singlelevel_faceDesc : ' // &
                     & 'found face 6-5 in face list, correct right elem id.'
        end if
      end if
  
  
      ! Checking in y-direction    
      write(*,*) 'y-direction'
      if( faces(level)%faces(2)%faceList%faceId%nvals.ne.4 ) then
        write(*,*) 'ERROR in check_parallel_singlelevel_faceDesc : ' // &
                    & 'expected 4 faces, but face description is wrong ...'
        call tem_abort()
      else 
        write(*,*) 'SUCCESS in check_parallel_singlelevel_faceDesc : ' // &
                    & 'number of nvals is 4.'
      end if
      ! Let us have a look if the face 7-5 (left elem id - right elem id) is
      ! in the list of faces.
      pos = PositionOfVal( me = faces(level)%faces(2)%faceList%faceId, val = 7_long_k )
      if(pos .le. 0) then
        write(*,*) 'ERROR in check_parallel_singlelevel_faceDesc : ' // &
                    & 'not able to locate face 7-5 in face list ...'
        call tem_abort()
      else 
        write(*,*) 'SUCCESS in check_parallel_singlelevel_faceDesc : ' // &
                    & 'found face 7-5 in face list.'
  
        if(faces(level)%faces(2)%faceList%rightElemId%val(pos).ne.5_long_k) then
          write(*,*) 'ERROR in check_parallel_singlelevel_faceDesc : ' // &
                     & 'not able to locate face 7-5 in face list (wrong ' // &
                     & 'right elem id) ...'
          call tem_abort()
        else 
          write(*,*) 'SUCCESS in check_parallel_singlelevel_faceDesc : ' // &
                     & 'found face 7-5 in face list, correct right elem id.'
        end if
      end if
    end if



    !! ComputeFace checks for
    ! process 0
    write(logUnit(1),*)"ComputeFace Checks in z-direction"
    call ComputeFace_to_TreeIDcheck(1,5,0, general%proc, faces=faces)
    call ComputeFace_to_TreeIDcheck(2,6,0, general%proc, faces=faces)
    call ComputeFace_to_TreeIDcheck(3,7,0, general%proc, faces=faces)
    call ComputeFace_to_TreeIDcheck(4,8,0, general%proc, faces=faces)
    ! process 1
    call ComputeFace_to_TreeIDcheck(5,1,1, general%proc, faces=faces)
    call ComputeFace_to_TreeIDcheck(6,2,1, general%proc, faces=faces)
    call ComputeFace_to_TreeIDcheck(7,3,1, general%proc, faces=faces)
    call ComputeFace_to_TreeIDcheck(8,4,1, general%proc, faces=faces)

  end subroutine check_parallel_singlelevel_faceDesc



   !Compute Face checks
   subroutine ComputeFace_to_TreeIDcheck(leftID, rightID, rank, proc, faces)
   !--------------------------------------------------------------------------
   integer, intent(in) :: rank
   integer, intent(in) :: leftID
   integer, intent(in) :: rightID
   type(tem_comm_env_type), intent(in) :: proc
   type(tem_face_type) :: faces(:)
   !--------------------------------------------------------------------------
   integer :: i, nfaces
   integer :: temp, found
   !--------------------------------------------------------------------------

   found = 0
   if(proc%rank.eq.rank) then
   nfaces=size(faces(1)%faces(3)%computeFace%facePos)
   do i=1,nfaces
      temp=faces(1)%faces(3)%computeFace%facePos(i)
      if(faces(1)%faces(3)%faceList%faceId%val(temp).eq.leftID) then
         found=1
         if(faces(1)%faces(3)%faceList%rightElemId%val(temp).eq.rightID) then
            write(*,*) 'SUCCESS in ComputeFace_to_TreeIDcheck : '//&
                        & 'found face',leftID,rightID,'in face list'
         else
            write(*,*) 'ERROR in ComputeFace_to_TreeIDcheck  : '//&
                        & 'Right ElemID', rightID, 'not found in FaceList'
            call tem_abort()
         end if
      end if
   end do
 
   if(found.ne.1) then
      write(*,*) 'ERROR in ComputeFace_to_TreeIDcheck : '//&
                  & 'Left ElemID', leftID, 'not found in Process', rank
      call tem_abort()
   end if

   end if


   end subroutine ComputeFace_to_TreeIDcheck


end program tem_face_test

! Copyright (c) 2016 Harald Klimach <harald.klimach@uni-siegen.de>
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
module mem_for_mpi_module
  use mpi_f08

  implicit none


contains


  subroutine alloc_mpif_mem(asize, baseptr, ierr)
    integer(kind=MPI_ADRESS_KIND), intent(in) :: asize
    type(c_ptr), intent(out) :: baseptr
    integer, optional, intent(out) :: ierr

    call MPI_Alloc_mem( size    = asize,         &
      &                 info    = MPI_INFO_NULL, &
      &                 baseptr = baseptr,       &
      &                 ierr    = ierr           )

  end subroutine alloc_mpif_mem


  subroutine free_mpif_mem(baseptr, ierr)
    type(c_ptr) :: baseptr
    integer, optional, intent(out) :: ierr

    call MPI_Free_mem( baseptr = baseptr, &
      &                ierr    = ierr     )

  end subroutine free_mpif_mem

end module mem_for_mpi_module

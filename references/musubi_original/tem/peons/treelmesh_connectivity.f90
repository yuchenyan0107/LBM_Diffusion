!> Helping module to find neighbors in the mesh and identify
!! the connectivity.
module connectivity_module
  use env_module,          only: long_k, newunit
  use treelmesh_module,    only: treelmesh_type
  use tem_bc_prop_module,  only: tem_BC_prop_type
  use tem_geometry_module, only: tem_posofid
  use tem_param_module,    only: qOffset
  use tem_property_module, only: prp_hasBnd
  use tem_topology_module, only: tem_coordofid, tem_idofcoord

  implicit none

  private

  public :: search_neighbors

contains

  !> Searching and printing neighbors for all elements.
  !!
  !! For each elements the neighbors in all directions are
  !! looked up, and if there is not a boundary condition
  !! in the given direction, the element in the neighboring
  !! direction is looked up.
  !! If the neighboring element is a ghost from finer,
  !! all children are recursively looked up and printed
  !! instead of the ghost element.
  !! The idea here is, that the ghost element depends on all
  !! these children and thus, they kind of belong to the
  !! dependency of the neighbor.
  !!
  !! Neighbors are defined according the the qOffset in
  !! the tem_params_module.
  subroutine search_neighbors(tree, BC, nNeighbors, outfile)
    type(treelmesh_type), intent(in)   :: tree
    type(tem_bc_prop_type), intent(in) :: BC
    integer, intent(in)                :: nNeighbors
    character(len=*), intent(in)       :: outfile

    integer :: iElem
    integer :: iBnd
    integer :: iNeighbor
    integer(kind=long_k) :: nTID
    integer :: maxneighbors
    integer :: elemcoord(4)
    integer :: neighoffset(4)
    integer :: neighCoord(4)
    integer :: fu

    fu = newunit()
    open(unit=fu, file=trim(outfile))

    maxneighbors = min(nNeighbors, BC%nSides)
    neighOffset = 0

    iBnd = 0
    do iElem=1,tree%nElems
      elemcoord = tem_coordofid(tree%treeID(iElem))

      if (btest(tree%ElemPropertyBits(iElem), prp_hasBnd)) then
        ! This element has boundaries, increase the boundary
        ! element counter and ensure, that directions with
        ! boundaries are ignored in the neighbor search.

        iBnd = iBnd + 1
        do iNeighbor=1,maxneighbors
          ! Only treat sides without boundaries:
          if (BC%boundary_ID(iNeighbor, iBnd) <= 0_long_k) then
            ! Negative values indicate explicit dependencies
            nTID = - BC%boundary_ID(iNeighbor, iBnd)
            if (nTID == 0_long_k) then
              neighOffset(1:3) = qOffset(iNeighbor,:)
              neighCoord = elemCoord + neighOffset
              nTID = tem_idofcoord(neighCoord)
            end if
            call find_and_write_tid(iElem, nTID, tree%treeID, fu)
          end if
        end do

      else

        ! No boundaries for this element, look up all neighbors
        do iNeighbor=1,nNeighbors
          neighOffset(1:3) = qOffset(iNeighbor,:)
          neighCoord = elemCoord + neighOffset
          nTID = tem_idofcoord(neighCoord)
          call find_and_write_tid(iElem, nTID, tree%treeID, fu)
        end do

      end if

    end do

    close(fu)

  end subroutine search_neighbors


  !> Recursive helping routine to find possibly refined neighbor
  !! elements.
  !!
  !! Found neighbors are immediatly written to the file provided
  !! in fu, to minimize overheads in memory.
  recursive subroutine find_and_write_tid(iElem, tid, treeid_list, fu)
    integer, intent(in)              :: iElem
    integer(kind=long_k), intent(in) :: tid
    integer(kind=long_k), intent(in) :: treeid_list(:)
    integer, intent(in)              :: fu

    integer :: iPos
    integer :: iChild
    integer(kind=long_k) :: childid
    integer(kind=long_k) :: off 
  
    iPos = tem_posofid(tid, treeid_list)
    if (iPos < 0) then
      off = tID * 8_long_k
      do iChild = 1,8
        childid = off + int(iChild, long_k)
        call find_and_write_tid(iElem, childid, treeid_list, fu)
      end do
    else if (iPos > 0) then
      write(fu, *) iElem, iPos
    end if
  end subroutine find_and_write_tid

end module connectivity_module


!> Short program to print the connectivity of a treelmesh.
!!
!! This program takes a treelmesh and generates the list
!! of neighbors for each element.
!! For each element-neighbor pairing, there is a seperate
!! line:
!! Element Neighbor
!! Each number in that list refers to the element position
!! in the list of all elements.
!!
!! ATTENTION: This program is not working in parallel, only
!!            start it with a single rank!
program treelmesh_connectivity
  use env_module,          only: fin_env, pathLen
  use aotus_module,        only: flu_State, aoterr_Fatal, close_config,        &
    &                             aot_get_val

  use treelmesh_module,    only: treelmesh_type, load_tem
  use tem_aux_module,      only: tem_abort, tem_open_distconf
  use tem_bc_prop_module,  only: tem_BC_prop_type, init_tem_bc_prop
  use tem_general_module,  only: tem_general_type, tem_load_general,           &
    &                            tem_start, tem_finalize
  use tem_logging_module,  only: logUnit

  use connectivity_module,     only: search_neighbors

  implicit none

  character(len=pathlen) :: filename
  type(flu_State) :: conf
  type(tem_general_type) :: general
  integer :: iError

  character(len=pathlen) :: outfile
  type(treelmesh_type) :: tree
  type(tem_BC_prop_type) :: BC_property
  integer :: nNeighbors

  call tem_start('connectivity', general)

  write(logUnit(1),*)'ATTENTION: Does not work in parallel!'

  call get_command_argument( 1,filename )
  if( trim(filename) .eq. '' ) then
    filename = 'connectivity.lua'
  end if
  write(logUnit(1),*) "Loading configuration file "//trim( filename )

  call tem_open_distconf(L = conf, filename = trim(filename),                  &
    &                    proc = general%proc)

  call tem_load_general( me = general, conf = conf )

  call load_tem( me     = tree,              &
    &            conf   = conf,              &
    &            myPart = general%proc%rank,      &
    &            nParts = general%proc%comm_size, &
    &            comm   = general%proc%comm       )

  call init_tem_bc_prop(tree   = tree,         &
    &                   mypart = general%proc%rank, &
    &                   comm   = general%proc%comm, &
    &                   bc     = BC_property   )

  call aot_get_val( L       = conf,        &
    &               key     = 'neighbors', &
    &               val     = nNeighbors,  &
    &               ErrCode = iError,      &
    &               default = 6            )

  write(logUnit(1),*) 'Looking for ', nNeighbors, ' nNeighbors'

  call aot_get_val( L       = conf,      &
    &               key     = 'outfile', &
    &               val     = outfile,   &
    &               ErrCode = iError     )

  if (btest(iError, aoterr_Fatal) .or. trim(outfile)=='') then
    write(logUnit(0),*) 'ERROR: you have to provide an outfile to which the'
    write(logUnit(0),*) '       results are to be written!'
    call tem_abort()
  end if

  write(logUnit(1),*) 'Writing results to '//trim(outfile)

  call close_config(conf)

  write(logUnit(2),*) 'Successfully read configuration'
  write(logUnit(2),*) 'Starting search for neighbors...'

  call search_neighbors(tree, BC_property, nNeighbors, trim(outfile))

  write(logUnit(1),*) 'Finished search for neighbors SUCCESSFULLY!'

  call tem_finalize(general)
  call fin_env()

end program treelmesh_connectivity

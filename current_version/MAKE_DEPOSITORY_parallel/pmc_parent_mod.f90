 MODULE pmc_parent

!------------------------------------------------------------------------------!
! This file is part of PALM.
!
! PALM is free software: you can redistribute it and/or modify it under the
! terms of the GNU General Public License as published by the Free Software
! Foundation, either version 3 of the License, or (at your option) any later
! version.
!
! PALM is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
! A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along with
! PALM. If not, see <http://www.gnu.org/licenses/>.
!
! Copyright 1997-2016 Leibniz Universitaet Hannover
!------------------------------------------------------------------------------!
!
! Current revisions:
! ------------------
! 
! 
! Former revisions:
! -----------------
! $Id: pmc_parent_mod.f90 2001 2016-08-20 18:41:22Z knoop $
!
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
! 
! 1938 2016-06-13 15:26:05Z hellstea
! Minor clean up.
! 
! 1901 2016-05-04 15:39:38Z raasch
! Module renamed. Code clean up. The words server/client changed to parent/child. 
!
! 1900 2016-05-04 15:27:53Z raasch
! re-formatted to match PALM style
!
! 1850 2016-04-08 13:29:27Z maronga
! Module renamed
! 
! 
! 1833 2016-04-07 14:23:03Z raasch
! gfortran requires pointer attributes for some array declarations,
! long line wrapped
!
! 1808 2016-04-05 19:44:00Z raasch
! MPI module used by default on all machines
!
! 1797 2016-03-21 16:50:28Z raasch
! introduction of different datatransfer modes
!
! 1791 2016-03-11 10:41:25Z raasch
! Debug write-statements commented out
!
! 1786 2016-03-08 05:49:27Z raasch
! change in child-parent data transfer: parent now gets data from child
! instead that child put's it to the parent
!
! 1779 2016-03-03 08:01:28Z raasch
! kind=dp replaced by wp,
! error messages removed or changed to PALM style, dim_order removed
! array management changed from linked list to sequential loop
!
! 1766 2016-02-29 08:37:15Z raasch
! modifications to allow for using PALM's pointer version
! +new routine pmc_s_set_active_data_array
!
! 1764 2016-02-28 12:45:19Z raasch
! cpp-statement added (nesting can only be used in parallel mode)
!
! 1762 2016-02-25 12:31:13Z hellstea
! Initial revision by K. Ketelsen
!
! Description:
! ------------
!
! Parent part of Palm Model Coupler
!-------------------------------------------------------------------------------!

#if defined( __parallel )
    USE, INTRINSIC ::  ISO_C_BINDING

#if defined( __mpifh )
    INCLUDE "mpif.h"
#else
    USE MPI
#endif
    USE kinds
    USE pmc_general,                                                            &
        ONLY: arraydef, childdef, da_namedef, da_namelen, pedef,                &
              pmc_g_setname, pmc_max_array, pmc_max_models, pmc_sort

    USE pmc_handle_communicator,                                                &
        ONLY: m_model_comm,m_model_rank,m_model_npes, m_to_child_comm,          &
              m_world_rank, pmc_parent_for_child

    USE pmc_mpi_wrapper,                                                        &
        ONLY: pmc_alloc_mem, pmc_bcast, pmc_time

   IMPLICIT NONE

   PRIVATE
   SAVE

   TYPE childindexdef
      INTEGER                              ::  nrpoints       !<
      INTEGER, DIMENSION(:,:), ALLOCATABLE ::  index_list_2d  !<
   END TYPE childindexdef

   TYPE(childdef), DIMENSION(pmc_max_models)       ::  children     !<
   TYPE(childindexdef), DIMENSION(pmc_max_models)  ::  indchildren  !<

   INTEGER ::  next_array_in_list = 0  !<


   PUBLIC pmc_parent_for_child


   INTERFACE pmc_parentinit
      MODULE PROCEDURE  pmc_parentinit
   END INTERFACE pmc_parentinit

    INTERFACE pmc_s_set_2d_index_list
        MODULE PROCEDURE pmc_s_set_2d_index_list
    END INTERFACE pmc_s_set_2d_index_list

    INTERFACE pmc_s_clear_next_array_list
        MODULE PROCEDURE pmc_s_clear_next_array_list
    END INTERFACE pmc_s_clear_next_array_list

    INTERFACE pmc_s_getnextarray
        MODULE PROCEDURE pmc_s_getnextarray
    END INTERFACE pmc_s_getnextarray

    INTERFACE pmc_s_set_dataarray
        MODULE PROCEDURE pmc_s_set_dataarray_2d
        MODULE PROCEDURE pmc_s_set_dataarray_3d
    END INTERFACE pmc_s_set_dataarray

    INTERFACE pmc_s_setind_and_allocmem
        MODULE PROCEDURE pmc_s_setind_and_allocmem
    END INTERFACE pmc_s_setind_and_allocmem

    INTERFACE pmc_s_fillbuffer
        MODULE PROCEDURE pmc_s_fillbuffer
    END INTERFACE pmc_s_fillbuffer

    INTERFACE pmc_s_getdata_from_buffer
        MODULE PROCEDURE pmc_s_getdata_from_buffer
    END INTERFACE pmc_s_getdata_from_buffer

    INTERFACE pmc_s_set_active_data_array
        MODULE PROCEDURE pmc_s_set_active_data_array
    END INTERFACE pmc_s_set_active_data_array

    PUBLIC pmc_parentinit, pmc_s_clear_next_array_list, pmc_s_fillbuffer,       &
           pmc_s_getdata_from_buffer, pmc_s_getnextarray,                       &
           pmc_s_setind_and_allocmem, pmc_s_set_active_data_array,              &
           pmc_s_set_dataarray, pmc_s_set_2d_index_list

 CONTAINS


 SUBROUTINE pmc_parentinit

    IMPLICIT NONE

    INTEGER ::  childid   !<
    INTEGER ::  i         !<
    INTEGER ::  j         !<
    INTEGER ::  istat     !<


    DO  i = 1, SIZE( pmc_parent_for_child )-1

       childid = pmc_parent_for_child( i )

       children(childid)%model_comm = m_model_comm
       children(childid)%inter_comm = m_to_child_comm(childid)

!
!--    Get rank and size
       CALL MPI_COMM_RANK( children(childid)%model_comm,                        &
                           children(childid)%model_rank, istat )
       CALL MPI_COMM_SIZE( children(childid)%model_comm,                        &
                           children(childid)%model_npes, istat )
       CALL MPI_COMM_REMOTE_SIZE( children(childid)%inter_comm,                 &
                                  children(childid)%inter_npes, istat )

!
!--    Intra communicater is used for MPI_GET
       CALL MPI_INTERCOMM_MERGE( children(childid)%inter_comm, .FALSE.,         &
                                 children(childid)%intra_comm, istat )
       CALL MPI_COMM_RANK( children(childid)%intra_comm,                        &
                           children(childid)%intra_rank, istat )

       ALLOCATE( children(childid)%pes(children(childid)%inter_npes))

!
!--    Allocate array of TYPE arraydef for all child PEs to store information
!--    of the transfer array
       DO  j = 1, children(childid)%inter_npes
         ALLOCATE( children(childid)%pes(j)%array_list(pmc_max_array) )
       ENDDO

       CALL get_da_names_from_child (childid)

    ENDDO

 END SUBROUTINE pmc_parentinit



 SUBROUTINE pmc_s_set_2d_index_list( childid, index_list )

     IMPLICIT NONE

     INTEGER, INTENT(IN)                    :: childid     !<
     INTEGER, DIMENSION(:,:), INTENT(INOUT) :: index_list  !<

     INTEGER ::  ian    !<
     INTEGER ::  ic     !<
     INTEGER ::  ie     !<
     INTEGER ::  ip     !<
     INTEGER ::  is     !<
     INTEGER ::  istat  !<
     INTEGER ::  n      !<


     IF ( m_model_rank == 0 )  THEN

!
!--     Sort to ascending parent PE order
        CALL pmc_sort( index_list, 6 )

        is = 1
        DO  ip = 0, m_model_npes-1

!
!--        Split into parent PEs
           ie = is - 1

!
!--        There may be no entry for this PE
           IF ( is <= SIZE( index_list,2 )  .AND.  ie >= 0 )  THEN

              DO WHILE ( index_list(6,ie+1 ) == ip )
                 ie = ie + 1
                 IF ( ie == SIZE( index_list,2 ) )  EXIT
              ENDDO

              ian = ie - is + 1

           ELSE
              is  = -1
              ie  = -2
              ian =  0
           ENDIF

!
!--        Send data to other parent PEs
           IF ( ip == 0 )  THEN
              indchildren(childid)%nrpoints = ian
              IF ( ian > 0)  THEN
                  ALLOCATE( indchildren(childid)%index_list_2d(6,ian) )
                  indchildren(childid)%index_list_2d(:,1:ian) =                 &
                                                             index_list(:,is:ie)
              ENDIF
           ELSE
              CALL MPI_SEND( ian, 1, MPI_INTEGER, ip, 1000, m_model_comm,       &
                             istat )
              IF ( ian > 0)  THEN
                  CALL MPI_SEND( index_list(1,is), 6*ian, MPI_INTEGER, ip,      &
                                 1001, m_model_comm, istat )
              ENDIF
           ENDIF
           is = ie + 1

        ENDDO

     ELSE

        CALL MPI_RECV( indchildren(childid)%nrpoints, 1, MPI_INTEGER, 0, 1000,  &
                       m_model_comm, MPI_STATUS_IGNORE, istat )
        ian = indchildren(childid)%nrpoints

        IF ( ian > 0 )  THEN
           ALLOCATE( indchildren(childid)%index_list_2d(6,ian) )
           CALL MPI_RECV( indchildren(childid)%index_list_2d, 6*ian,            &
                          MPI_INTEGER, 0, 1001, m_model_comm,                   &
                          MPI_STATUS_IGNORE, istat)
        ENDIF

     ENDIF

     CALL set_pe_index_list( childid, children(childid),                        &
                             indchildren(childid)%index_list_2d,                &
                             indchildren(childid)%nrpoints )

 END SUBROUTINE pmc_s_set_2d_index_list



 SUBROUTINE pmc_s_clear_next_array_list

    IMPLICIT NONE

    next_array_in_list = 0

 END SUBROUTINE pmc_s_clear_next_array_list



 LOGICAL FUNCTION pmc_s_getnextarray( childid, myname )

!
!-- List handling is still required to get minimal interaction with
!-- pmc_interface
!-- TODO: what does "still" mean? Is there a chance to change this!
    CHARACTER(LEN=*), INTENT(OUT) ::  myname    !<
    INTEGER(iwp), INTENT(IN)      ::  childid   !<

    TYPE(arraydef), POINTER :: ar
    TYPE(pedef), POINTER    :: ape

    next_array_in_list = next_array_in_list + 1

!
!-- Array names are the same on all children PEs, so take first PE to get the name
    ape => children(childid)%pes(1)

    IF ( next_array_in_list > ape%nr_arrays )  THEN

!
!--    All arrays are done
       pmc_s_getnextarray = .FALSE.
       RETURN
    ENDIF

    ar => ape%array_list(next_array_in_list)
    myname = ar%name

!
!-- Return true if legal array
!-- TODO: what does this comment mean? Can there be non-legal arrays??
    pmc_s_getnextarray = .TRUE.

 END FUNCTION pmc_s_getnextarray



 SUBROUTINE pmc_s_set_dataarray_2d( childid, array, array_2 )

    IMPLICIT NONE

    INTEGER,INTENT(IN) ::  childid   !<

    REAL(wp), INTENT(IN), DIMENSION(:,:), POINTER           ::  array    !<
    REAL(wp), INTENT(IN), DIMENSION(:,:), POINTER, OPTIONAL ::  array_2  !<

    INTEGER               ::  nrdims      !<
    INTEGER, DIMENSION(4) ::  dims        !<
    TYPE(C_PTR)           ::  array_adr   !<
    TYPE(C_PTR)           ::  second_adr  !<


    dims      = 1
    nrdims    = 2
    dims(1)   = SIZE( array,1 )
    dims(2)   = SIZE( array,2 )
    array_adr = C_LOC( array )

    IF ( PRESENT( array_2 ) )  THEN
       second_adr = C_LOC(array_2)
       CALL pmc_s_setarray( childid, nrdims, dims, array_adr,                   &
                            second_adr = second_adr)
    ELSE
       CALL pmc_s_setarray( childid, nrdims, dims, array_adr )
    ENDIF

 END SUBROUTINE pmc_s_set_dataarray_2d



 SUBROUTINE pmc_s_set_dataarray_3d( childid, array, nz_cl, nz, array_2 )

    IMPLICIT NONE

    INTEGER, INTENT(IN) ::  childid   !<
    INTEGER, INTENT(IN) ::  nz        !<
    INTEGER, INTENT(IN) ::  nz_cl     !<

    REAL(wp), INTENT(IN), DIMENSION(:,:,:), POINTER           ::  array    !<
    REAL(wp), INTENT(IN), DIMENSION(:,:,:), POINTER, OPTIONAL ::  array_2  !<

    INTEGER               ::  nrdims      !<
    INTEGER, DIMENSION(4) ::  dims        !<
    TYPE(C_PTR)           ::  array_adr   !<
    TYPE(C_PTR)           ::  second_adr  !<

!
!-- TODO: the next assignment seems to be obsolete. Please check!
    dims      = 1
    dims      = 0
    nrdims    = 3
    dims(1)   = SIZE( array,1 )
    dims(2)   = SIZE( array,2 )
    dims(3)   = SIZE( array,3 )
    dims(4)   = nz_cl+dims(1)-nz  ! works for first dimension 1:nz and 0:nz+1

    array_adr = C_LOC(array)

!
!-- In PALM's pointer version, two indices have to be stored internally.
!-- The active address of the data array is set in swap_timelevel.
    IF ( PRESENT( array_2 ) )  THEN
      second_adr = C_LOC( array_2 )
      CALL pmc_s_setarray( childid, nrdims, dims, array_adr,                    &
                           second_adr = second_adr)
    ELSE
       CALL pmc_s_setarray( childid, nrdims, dims, array_adr )
    ENDIF

 END SUBROUTINE pmc_s_set_dataarray_3d



 SUBROUTINE pmc_s_setind_and_allocmem( childid )

    USE control_parameters,                                                     &
        ONLY:  message_string

    IMPLICIT NONE

!
!-- Naming convention for appendices:   _pc  -> parent to child transfer
!--                                     _cp  -> child to parent transfer
!--                                     send -> parent to child transfer
!--                                     recv -> child to parent transfer
    INTEGER, INTENT(IN) ::  childid   !<

    INTEGER                        ::  arlen    !<
    INTEGER                        ::  i        !<
    INTEGER                        ::  ierr     !<
    INTEGER                        ::  istat    !<
    INTEGER                        ::  j        !<
    INTEGER                        ::  myindex  !<
    INTEGER                        ::  rcount   !< count MPI requests
    INTEGER                        ::  tag      !<

    INTEGER(idp)                   ::  bufsize  !< size of MPI data window
    INTEGER(KIND=MPI_ADDRESS_KIND) ::  winsize  !<

    INTEGER, DIMENSION(1024)       ::  req      !<

    TYPE(C_PTR)             ::  base_ptr  !<
    TYPE(pedef), POINTER    ::  ape       !<
    TYPE(arraydef), POINTER ::  ar        !<

    REAL(wp),DIMENSION(:), POINTER, SAVE ::  base_array_pc  !< base array for parent to child transfer
    REAL(wp),DIMENSION(:), POINTER, SAVE ::  base_array_cp  !< base array for child to parent transfer

!
!-- Parent to child direction
    myindex = 1
    rcount  = 0
    bufsize = 8

!
!-- First stride: compute size and set index
    DO  i = 1, children(childid)%inter_npes

       ape => children(childid)%pes(i)
       tag = 200

       DO  j = 1, ape%nr_arrays

          ar  => ape%array_list(j)
          IF ( ar%nrdims == 2 )  THEN
             arlen = ape%nrele
          ELSEIF ( ar%nrdims == 3 )  THEN
             arlen = ape%nrele * ar%a_dim(4)
          ELSE
             arlen = -1
          ENDIF
          ar%sendindex = myindex

          tag    = tag + 1
          rcount = rcount + 1
          CALL MPI_ISEND( myindex, 1, MPI_INTEGER, i-1, tag,                    &
                          children(childid)%inter_comm, req(rcount), ierr )

!
!--       Maximum of 1024 outstanding requests
!--       TODO: what does this limit mean?
          IF ( rcount == 1024 )  THEN
             CALL MPI_WAITALL( rcount, req, MPI_STATUSES_IGNORE, ierr )
             rcount = 0
          ENDIF

          myindex = myindex + arlen
          bufsize = bufsize + arlen
          ar%sendsize = arlen

       ENDDO

       IF ( rcount > 0 )  THEN
          CALL MPI_WAITALL( rcount, req, MPI_STATUSES_IGNORE, ierr )
       ENDIF

    ENDDO

!
!-- Create RMA (One Sided Communication) window for data buffer parent to
!-- child transfer.
!-- The buffer of MPI_GET (counterpart of transfer) can be PE-local, i.e.
!-- it can but must not be part of the MPI RMA window. Only one RMA window is
!-- required to prepare the data for
!--                       parent -> child transfer on the parent side
!-- and for
!--                       child -> parent transfer on the child side
    CALL pmc_alloc_mem( base_array_pc, bufsize )
    children(childid)%totalbuffersize = bufsize * wp

    winsize = bufsize * wp
    CALL MPI_WIN_CREATE( base_array_pc, winsize, wp, MPI_INFO_NULL,             &
                         children(childid)%intra_comm,                          &
                         children(childid)%win_parent_child, ierr )

!
!-- Open window to set data
    CALL MPI_WIN_FENCE( 0, children(childid)%win_parent_child, ierr )

!
!-- Second stride: set buffer pointer
    DO  i = 1, children(childid)%inter_npes

       ape => children(childid)%pes(i)

       DO  j = 1, ape%nr_arrays

          ar => ape%array_list(j)
          ar%sendbuf = C_LOC( base_array_pc(ar%sendindex) )

          IF ( ar%sendindex + ar%sendsize > bufsize )  THEN             
             WRITE( message_string, '(a,i4,4i7,1x,a)' )                         &
                    'parent buffer too small ',i,                               &
                    ar%sendindex,ar%sendsize,ar%sendindex+ar%sendsize,          &
                    bufsize,trim(ar%name)
             CALL message( 'pmc_s_setind_and_allocmem', 'PA0429', 3, 2, 0, 6, 0 )
          ENDIF
       ENDDO
    ENDDO

!
!-- Child to parent direction
    bufsize = 8

!
!-- First stride: compute size and set index
    DO  i = 1, children(childid)%inter_npes

       ape => children(childid)%pes(i)
       tag = 300

       DO  j = 1, ape%nr_arrays

          ar => ape%array_list(j)

!
!--       Receive index from child
          tag = tag + 1
          CALL MPI_RECV( myindex, 1, MPI_INTEGER, i-1, tag,                     &
                         children(childid)%inter_comm, MPI_STATUS_IGNORE, ierr )

          IF ( ar%nrdims == 3 )  THEN
             bufsize = MAX( bufsize, ape%nrele * ar%a_dim(4) )
          ELSE
             bufsize = MAX( bufsize, ape%nrele )
          ENDIF
          ar%recvindex = myindex

        ENDDO

    ENDDO

!
!-- Create RMA (one sided communication) data buffer.
!-- The buffer for MPI_GET can be PE local, i.e. it can but must not be part of
!-- the MPI RMA window
    CALL pmc_alloc_mem( base_array_cp, bufsize, base_ptr )
    children(childid)%totalbuffersize = bufsize * wp

    CALL MPI_BARRIER( children(childid)%intra_comm, ierr )

!
!-- Second stride: set buffer pointer
    DO  i = 1, children(childid)%inter_npes

       ape => children(childid)%pes(i)

       DO  j = 1, ape%nr_arrays
          ar => ape%array_list(j)
          ar%recvbuf = base_ptr
       ENDDO

    ENDDO

 END SUBROUTINE pmc_s_setind_and_allocmem



 SUBROUTINE pmc_s_fillbuffer( childid, waittime )

    IMPLICIT NONE

    INTEGER, INTENT(IN)             ::  childid   !<

    REAL(wp), INTENT(OUT), OPTIONAL ::  waittime  !<

    INTEGER               ::  ierr     !<
    INTEGER               ::  ij       !<
    INTEGER               ::  ip       !<
    INTEGER               ::  istat    !<
    INTEGER               ::  j        !<
    INTEGER               ::  myindex  !<

    INTEGER, DIMENSION(1) ::  buf_shape

    REAL(wp)                            ::  t1       !<
    REAL(wp)                            ::  t2       !<
    REAL(wp), POINTER, DIMENSION(:)     ::  buf      !<
    REAL(wp), POINTER, DIMENSION(:,:)   ::  data_2d  !<
    REAL(wp), POINTER, DIMENSION(:,:,:) ::  data_3d  !<

    TYPE(pedef), POINTER    ::  ape  !<
    TYPE(arraydef), POINTER ::  ar   !<

!
!-- Synchronization of the model is done in pmci_synchronize.
!-- Therefor the RMA window can be filled without
!-- sychronization at this point and a barrier is not necessary.
!-- Please note that waittime has to be set in pmc_s_fillbuffer AND
!-- pmc_c_getbuffer
    IF ( PRESENT( waittime) )  THEN
      t1 = pmc_time()
      CALL MPI_BARRIER( children(childid)%intra_comm, ierr )
      t2 = pmc_time()
      waittime = t2- t1
    ENDIF

    DO  ip = 1, children(childid)%inter_npes

       ape => children(childid)%pes(ip)

       DO  j = 1, ape%nr_arrays

          ar => ape%array_list(j)
          myindex = 1

          IF ( ar%nrdims == 2 )  THEN

             buf_shape(1) = ape%nrele
             CALL C_F_POINTER( ar%sendbuf, buf, buf_shape )
             CALL C_F_POINTER( ar%data, data_2d, ar%a_dim(1:2) )
             DO  ij = 1, ape%nrele
                buf(myindex) = data_2d(ape%locind(ij)%j,ape%locind(ij)%i)
                myindex = myindex + 1
             ENDDO

          ELSEIF ( ar%nrdims == 3 )  THEN

             buf_shape(1) = ape%nrele*ar%a_dim(4)
             CALL C_F_POINTER( ar%sendbuf, buf, buf_shape )
             CALL C_F_POINTER( ar%data, data_3d, ar%a_dim(1:3) )
             DO  ij = 1, ape%nrele
                buf(myindex:myindex+ar%a_dim(4)-1) =                            &
                        data_3d(1:ar%a_dim(4),ape%locind(ij)%j,ape%locind(ij)%i)
                myindex = myindex + ar%a_dim(4)
             ENDDO

          ENDIF

        ENDDO

    ENDDO

!
!-- Buffer is filled
    CALL MPI_BARRIER( children(childid)%intra_comm, ierr )

 END SUBROUTINE pmc_s_fillbuffer



 SUBROUTINE pmc_s_getdata_from_buffer( childid, waittime )

    IMPLICIT NONE

    INTEGER, INTENT(IN)             ::  childid      !<
    REAL(wp), INTENT(OUT), OPTIONAL ::  waittime     !<

    INTEGER                        ::  ierr          !<
    INTEGER                        ::  ij            !<
    INTEGER                        ::  ip            !<
    INTEGER                        ::  istat         !<
    INTEGER                        ::  j             !<
    INTEGER                        ::  myindex       !<
    INTEGER                        ::  nr            !<
    INTEGER                        ::  target_pe     !<
    INTEGER(kind=MPI_ADDRESS_KIND) ::  target_disp   !<

    INTEGER, DIMENSION(1)          ::  buf_shape     !<

    REAL(wp)                            ::  t1       !<
    REAL(wp)                            ::  t2       !<
    REAL(wp), POINTER, DIMENSION(:)     ::  buf      !<
    REAL(wp), POINTER, DIMENSION(:,:)   ::  data_2d  !<
    REAL(wp), POINTER, DIMENSION(:,:,:) ::  data_3d  !<

    TYPE(pedef), POINTER    ::  ape  !<
    TYPE(arraydef), POINTER ::  ar   !<


    t1 = pmc_time()

!
!-- Wait for child to fill buffer
    CALL MPI_BARRIER( children(childid)%intra_comm, ierr )
    t2 = pmc_time() - t1
    IF ( PRESENT( waittime ) )  waittime = t2

!
!-- TODO: check next statement
!-- Fence might do it, test later
!-- CALL MPI_WIN_FENCE( 0, children(childid)%win_parent_child, ierr)
    CALL MPI_BARRIER( children(childid)%intra_comm, ierr )

    DO  ip = 1, children(childid)%inter_npes

       ape => children(childid)%pes(ip)

       DO  j = 1, ape%nr_arrays

          ar => ape%array_list(j)

          IF ( ar%recvindex < 0 )  CYCLE

          IF ( ar%nrdims == 2 )  THEN
             nr = ape%nrele
          ELSEIF ( ar%nrdims == 3 )  THEN
             nr = ape%nrele * ar%a_dim(4)
          ENDIF

          buf_shape(1) = nr
          CALL C_F_POINTER( ar%recvbuf, buf, buf_shape )

!
!--       MPI passive target RMA
          IF ( nr > 0 )  THEN
             target_disp = ar%recvindex - 1

!
!--          Child PEs are located behind parent PEs
             target_pe = ip - 1 + m_model_npes
             CALL MPI_WIN_LOCK( MPI_LOCK_SHARED, target_pe, 0,                  &
                                children(childid)%win_parent_child, ierr )
             CALL MPI_GET( buf, nr, MPI_REAL, target_pe, target_disp, nr,       &
                           MPI_REAL, children(childid)%win_parent_child, ierr )
             CALL MPI_WIN_UNLOCK( target_pe,                                    &
                                  children(childid)%win_parent_child, ierr )
          ENDIF

          myindex = 1
          IF ( ar%nrdims == 2 )  THEN

             CALL C_F_POINTER( ar%data, data_2d, ar%a_dim(1:2) )
             DO  ij = 1, ape%nrele
                data_2d(ape%locind(ij)%j,ape%locind(ij)%i) = buf(myindex)
                myindex = myindex + 1
             ENDDO

          ELSEIF ( ar%nrdims == 3 )  THEN

             CALL C_F_POINTER( ar%data, data_3d, ar%a_dim(1:3))
             DO  ij = 1, ape%nrele
                data_3d(1:ar%a_dim(4),ape%locind(ij)%j,ape%locind(ij)%i) =      &
                                              buf(myindex:myindex+ar%a_dim(4)-1)
                myindex = myindex + ar%a_dim(4)
             ENDDO

          ENDIF

       ENDDO

    ENDDO

 END SUBROUTINE pmc_s_getdata_from_buffer



 SUBROUTINE get_da_names_from_child( childid )

!
!-- Get data array description and name from child
    IMPLICIT NONE

    INTEGER, INTENT(IN) ::  childid  !<

    TYPE(da_namedef) ::  myname  !<

    DO
       CALL pmc_bcast( myname%couple_index, 0, comm=m_to_child_comm(childid) )
       IF ( myname%couple_index == -1 )  EXIT
       CALL pmc_bcast( myname%parentdesc,   0, comm=m_to_child_comm(childid) )
       CALL pmc_bcast( myname%nameonparent, 0, comm=m_to_child_comm(childid) )
       CALL pmc_bcast( myname%childdesc,    0, comm=m_to_child_comm(childid) )
       CALL pmc_bcast( myname%nameonchild,  0, comm=m_to_child_comm(childid) )

       CALL pmc_g_setname( children(childid), myname%couple_index,              &
                           myname%nameonparent )
   ENDDO

 END SUBROUTINE get_da_names_from_child



 SUBROUTINE pmc_s_setarray(childid, nrdims, dims, array_adr, second_adr )

!
!-- Set array for child inter PE 0
    IMPLICIT NONE

    INTEGER, INTENT(IN)               ::  childid    !<
    INTEGER, INTENT(IN)               ::  nrdims     !<
    INTEGER, INTENT(IN), DIMENSION(:) ::  dims       !<

    TYPE(C_PTR), INTENT(IN)           :: array_adr   !<
    TYPE(C_PTR), INTENT(IN), OPTIONAL :: second_adr  !<

    INTEGER ::  i  !< local counter

    TYPE(pedef), POINTER    ::  ape  !<
    TYPE(arraydef), POINTER ::  ar   !<


    DO  i = 1, children(childid)%inter_npes

       ape => children(childid)%pes(i)
       ar  => ape%array_list(next_array_in_list)
       ar%nrdims = nrdims
       ar%a_dim  = dims
       ar%data   = array_adr

       IF ( PRESENT( second_adr ) )  THEN
          ar%po_data(1) = array_adr
          ar%po_data(2) = second_adr
       ELSE
          ar%po_data(1) = C_NULL_PTR
          ar%po_data(2) = C_NULL_PTR
       ENDIF

    ENDDO

 END SUBROUTINE pmc_s_setarray



 SUBROUTINE pmc_s_set_active_data_array( childid, iactive )

    IMPLICIT NONE

    INTEGER, INTENT(IN) ::  childid   !<
    INTEGER, INTENT(IN) ::  iactive   !<

    INTEGER :: i   !<
    INTEGER :: ip  !<
    INTEGER :: j   !<

    TYPE(pedef), POINTER    ::  ape  !<
    TYPE(arraydef), POINTER ::  ar   !<

    DO  ip = 1, children(childid)%inter_npes

       ape => children(childid)%pes(ip)

       DO  j = 1, ape%nr_arrays

          ar => ape%array_list(j)
          IF ( iactive == 1  .OR.  iactive == 2 )  THEN
             ar%data = ar%po_data(iactive)
          ENDIF

       ENDDO

    ENDDO

 END SUBROUTINE pmc_s_set_active_data_array



 SUBROUTINE set_pe_index_list( childid, mychild, index_list, nrp )

    IMPLICIT NONE

    INTEGER, INTENT(IN)                 ::  childid     !<
    INTEGER, INTENT(IN), DIMENSION(:,:) ::  index_list  !<
    INTEGER, INTENT(IN)                 ::  nrp         !<

    TYPE(childdef), INTENT(INOUT)       ::  mychild     !<

    INTEGER                                 :: i        !<
    INTEGER                                 :: ierr     !<
    INTEGER                                 :: ind      !<
    INTEGER                                 :: indwin   !<
    INTEGER                                 :: indwin2  !<
    INTEGER                                 :: i2       !<
    INTEGER                                 :: j        !<
    INTEGER                                 :: rempe    !<
    INTEGER(KIND=MPI_ADDRESS_KIND)          :: winsize  !<

    INTEGER, DIMENSION(mychild%inter_npes)  :: remind   !<

    INTEGER, DIMENSION(:), POINTER          :: remindw  !<
    INTEGER, DIMENSION(:), POINTER          :: rldef    !<

    TYPE(pedef), POINTER                    :: ape      !<

!
!-- First, count entries for every remote child PE
    DO  i = 1, mychild%inter_npes
       ape => mychild%pes(i)
       ape%nrele = 0
    ENDDO

!
!-- Loop over number of coarse grid cells
    DO  j = 1, nrp
       rempe = index_list(5,j) + 1   ! PE number on remote PE
       ape => mychild%pes(rempe)
       ape%nrele = ape%nrele + 1     ! Increment number of elements for this child PE
    ENDDO

    DO  i = 1, mychild%inter_npes
       ape => mychild%pes(i)
       ALLOCATE( ape%locind(ape%nrele) )
    ENDDO

    remind = 0

!
!-- Second, create lists
!-- Loop over number of coarse grid cells
    DO  j = 1, nrp
       rempe = index_list(5,j) + 1
       ape => mychild%pes(rempe)
       remind(rempe)     = remind(rempe)+1
       ind               = remind(rempe)
       ape%locind(ind)%i = index_list(1,j)
       ape%locind(ind)%j = index_list(2,j)
    ENDDO

!
!-- Prepare number of elements for children PEs
    CALL pmc_alloc_mem( rldef, mychild%inter_npes*2 )

!
!-- Number of child PEs * size of INTEGER (i just arbitrary INTEGER)
    winsize = mychild%inter_npes*c_sizeof(i)*2

    CALL MPI_WIN_CREATE( rldef, winsize, iwp, MPI_INFO_NULL,                    &
                         mychild%intra_comm, indwin, ierr )

!
!-- Open window to set data
    CALL MPI_WIN_FENCE( 0, indwin, ierr )

    rldef(1) = 0            ! index on remote PE 0
    rldef(2) = remind(1)    ! number of elements on remote PE 0

!
!-- Reserve buffer for index array
    DO  i = 2, mychild%inter_npes
       i2          = (i-1) * 2 + 1
       rldef(i2)   = rldef(i2-2) + rldef(i2-1) * 2  ! index on remote PE
       rldef(i2+1) = remind(i)                      ! number of elements on remote PE
    ENDDO

!
!-- Close window to allow child to access data
    CALL MPI_WIN_FENCE( 0, indwin, ierr )

!
!-- Child has retrieved data
    CALL MPI_WIN_FENCE( 0, indwin, ierr )

    i2 = 2 * mychild%inter_npes - 1
    winsize = ( rldef(i2) + rldef(i2+1) ) * 2

!
!-- Make sure, MPI_ALLOC_MEM works
    winsize = MAX( winsize, 1 )

    CALL pmc_alloc_mem( remindw, INT( winsize ) )

    CALL MPI_BARRIER( m_model_comm, ierr )
    CALL MPI_WIN_CREATE( remindw, winsize*c_sizeof(i), iwp, MPI_INFO_NULL,      &
                         mychild%intra_comm, indwin2, ierr )
!
!-- Open window to set data
    CALL MPI_WIN_FENCE( 0, indwin2, ierr )

!
!-- Create the 2D index list
    DO  j = 1, nrp
       rempe = index_list(5,j) + 1    ! PE number on remote PE
       ape => mychild%pes(rempe)
       i2    = rempe * 2 - 1
       ind   = rldef(i2) + 1
       remindw(ind)   = index_list(3,j)
       remindw(ind+1) = index_list(4,j)
       rldef(i2)      = rldef(i2)+2
    ENDDO

!
!-- All data are set
    CALL MPI_WIN_FENCE( 0, indwin2, ierr )

!
!-- Don't know why, but this barrier is necessary before windows can be freed
!-- TODO: find out why this is required
    CALL MPI_BARRIER( mychild%intra_comm, ierr )

    CALL MPI_WIN_FREE( indwin, ierr )
    CALL MPI_WIN_FREE( indwin2, ierr )

!
!-- TODO: check if the following idea needs to be done
!-- Sollte funktionieren, Problem mit MPI implementation
!-- https://www.lrz.de/services/software/parallel/mpi/onesided
!-- CALL MPI_Free_mem (remindw, ierr)

 END SUBROUTINE set_pe_index_list

#endif
 END MODULE pmc_parent

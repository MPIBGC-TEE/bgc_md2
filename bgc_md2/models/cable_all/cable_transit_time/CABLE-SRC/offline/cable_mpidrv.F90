!==============================================================================
! This source code is part of the 
! Australian Community Atmosphere Biosphere Land Exchange (CABLE) model.
! This work is licensed under the CABLE Academic User Licence Agreement 
! (the "Licence").
! You may not use this file except in compliance with the Licence.
! A copy of the Licence and registration form can be obtained from 
! http://www.accessimulator.org.au/cable
! You need to register and read the Licence agreement before use.
! Please contact cable_help@nf.nci.org.au for any questions on 
! registration and the Licence.
!
! Unless required by applicable law or agreed to in writing, 
! software distributed under the Licence is distributed on an "AS IS" BASIS, 
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the Licence for the specific language governing permissions and 
! limitations under the Licence.
! ==============================================================================
!
! Purpose: Bare bones MPI driver for CABLE
!
! Contact: Bernard.Pak@csiro.au
!
! History: MPI wrapper developed by Maciej Golebiewski (2012)
!
! ==============================================================================
!
PROGRAM mpi_driver

  USE mpi

  USE cable_mpicommon
  USE cable_mpimaster
  USE cable_mpiworker

  IMPLICIT NONE

  INTEGER:: comm, np, rank, ierr
  LOGICAL:: flag

  CALL MPI_Init (ierr)
  CALL MPI_Comm_dup (MPI_COMM_WORLD, comm, ierr)
  CALL MPI_Comm_size (comm, np, ierr)

  IF (np < 2) THEN
    WRITE (*,*) 'This program needs at least 2 processes to run!'
    CALL MPI_Abort (comm, 0, ierr)
  END IF

  CALL MPI_Comm_rank (comm, rank, ierr)
  !  print*,'after MPICommrank, ierr',ierr, rank

  IF (rank == 0) THEN
    WRITE(*,*) 'start to run master'
    CALL mpidrv_master (comm)
    !        print*,'finish run master'
  ELSE
    WRITE(*,*) 'start to run worker'
    CALL mpidrv_worker (comm, rank)
    !        print*,'finish run worker',rank
  END IF
  print*,'error index, rank',ierr, rank
  CALL MPI_Finalize (ierr)
  print*,'before finalized',flag, rank
  CALL MPI_Finalized (flag, ierr)
  print*,'after finalized',flag, rank

END PROGRAM mpi_driver


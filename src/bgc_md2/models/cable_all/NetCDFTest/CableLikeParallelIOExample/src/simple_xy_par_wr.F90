!     http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-tutorial
module constants
  real(kind=8), parameter    :: pi = 3.1415927                
end module

program simple_xy_par_wr

  use netcdf
  use mpi
  use constants

  implicit none
 
  INTEGER :: ierr,wnp,rank_in_world,nproc
  integer :: world_group_id,worker_group_id,rank_in_worker_comm,wnptest,worker_comm_id
  character (len = *), parameter :: SERIAL_FILE_NAME = "simple_xy.nc"
  character (len = *), parameter :: PARALLEL_FILE_NAME = "simple_xy_par.nc"
  ! We assume the points to be threaded on a single string
  ! We are writing a temperature timeline 1D
  ! for every point. So the dimension for the temperature data is 2
  integer,      parameter :: NDIMS_NPP= 2 ,NDIMS_temperature= 2 

  integer                 ::t,i, j,k,l,mland,mp

  ! When we create netCDF files, variables and dimensions, we get back
  ! an ID for each one.
  integer :: ncid_par,ncid_ser, varid_NPP, varid_temperature&
            ,dimids_temperature(NDIMS_temperature)&
            ,dimids_NPP(NDIMS_NPP)&
            ,landpoint_dimid, patch_dimid,time_dimid&
            ,my_lp_index_offset ,my_lp_number&
            ,my_patch_ind_offset ,my_patch_number&
            ,patch_dimid_NPP, time_dimid_NPP
            
  ! add chunk size for unlimited variables
  integer :: chunk_size_temperature(NDIMS_temperature)
  integer :: chunk_size_NPP(NDIMS_NPP)

  ! These will tell where in the data file this processor should write.
  integer :: start_NPP(NDIMS_NPP)                 , count_NPP(NDIMS_NPP)
  integer :: start_temperature(NDIMS_temperature) , count_temperature(NDIMS_temperature)
  
  ! This is the data array we will write. It will just be filled with
  ! the rank_in_world of this processor.
  integer, allocatable :: data_out(:)
  real(kind=8), allocatable :: temperature_out(:) ,NPP_out(:)
  integer, allocatable :: worker_ranks(:) !array of the ranks of the workers. Note that the entries are w.r.t. COMM_WORLD and not
  real, allocatable :: temp_out(:,:)
  
  real(kind=8), parameter                       :: min_x=0,max_x=1,span_x=max_x-min_x
  real(kind=8)                                  :: d_t,d_x ,my_offset ,patch_extend

  real(kind=8), allocatable, dimension(:)       :: my_xs
  real(kind=8), allocatable, dimension(:,:)     :: my_temperatures
  logical :: i_am_a_worker,i_am_master


  ! Initialize MPI, learn local rank and total number of processors.
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, rank_in_world, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)
  call MPI_Comm_group ( MPI_COMM_WORLD, world_group_id, ierr )
  wnp=nproc-1
  allocate ( worker_ranks(1:wnp)) 
  do i = 1, wnp ! we start with 1 since the master will be excluded 
    worker_ranks(i) = i
  end do
  
  !if (rank_in_world ==1) print *,'worker_ranks: ',worker_ranks
  call MPI_Comm_group ( MPI_COMM_WORLD, world_group_id, ierr ) !world group is now the group of all processes in MPI_COMM_WORLD
  call MPI_Group_incl ( world_group_id, wnp,worker_ranks, worker_group_id, ierr )
  ! create the new communicator, the variable worker_comm_id
  ! will be equal to MPI_Comm_null if our process does not belong to the
  ! group and some integer otherwise
  call MPI_Comm_create ( MPI_COMM_WORLD, worker_group_id, worker_comm_id, ierr ) 
  i_am_a_worker=worker_comm_id /= MPI_Comm_null
  i_am_master =rank_in_world==0
  if (i_am_master) then
    ! open a different netcdf file in serial mode and write some 
    ! information to it that does not depend on the workers
    ! we could also use the same file as the workers but would have
    ! to close it first in the master and AFTERWARDS open it in parallel
    ! mode in the workers which would force us to syncronize.
    ! If we use an extra file for the master we can access it 
    ! asyncronously.
    ! After the program is finished we can add the serial variables to the
    ! file that was previously opened for parallel access.
    ! 
    ! I t might be possibe to achive the same effect only using 
    ! the file in parallel mode and let the master write a different variable
    ! but this is not clear yet

    call check(nf90_create(SERIAL_FILE_NAME, NF90_NETCDF4, ncid_ser))
    call check( nf90_close(ncid_ser) )

  endif

  if (i_am_a_worker) then
    ! only workers can get their rank and the size of their communicator
    !call MPI_Comm_rank(worker_comm_id, rank_in_worker_comm, ierr)
    !call MPI_Comm_size(worker_comm_id, wnptest, ierr)
    
    ! get cable like information about which worker owns which landpoints and
    ! patches
    call get_offsets_and_sizes(&
      my_lp_index_offset&
      ,my_lp_number&
      ,my_patch_ind_offset&
      ,my_patch_number&
      ,mland&
      ,mp&
      ,rank_in_worker_comm&
      ,wnp&
      ,worker_comm_id&
    )
    
    ! HDF5/netCDF-4 file to be created. The comm and info parameters
    ! cause parallel I/O to be enabled. Use either NF90_MPIIO or
    ! NF90_MPIPOSIX to select between MPI/IO and MPI/POSIX.
    call check(nf90_create(PARALLEL_FILE_NAME, IOR(NF90_NETCDF4, NF90_MPIIO), ncid_par, &
         comm = worker_comm_id, info = MPI_INFO_NULL))

    ! Define the dimensions. NetCDF will hand back an ID for
    ! each. Metadata operations must take place on all processors.
    call check(nf90_def_dim(ncid_par, "lp_temperature", mland, landpoint_dimid))
    call check(nf90_def_dim(ncid_par, "time_temperature", NF90_UNLIMITED, time_dimid)) 

    call check(nf90_def_dim(ncid_par, "patch_NPP", mp, patch_dimid))

    ! Note that in fortran arrays are stored in
    ! column-major format.
    dimids_temperature= (/ landpoint_dimid,time_dimid/)
    dimids_NPP= (/ patch_dimid,time_dimid/)

    ! define the chunk size (1 ayg unlimited time dimension)
    ! chunk_size_temperature= (/ my_lp_number, 1/) !This seems to be difficult since here
    ! the chunksize depends on the rank which leads to value displacement 
    !chunk_size_temperature= (/ 1, 1/) !This seems to be a safe bet if this is really a global variable
    chunk_size_temperature= (/ mland, 1/) !This seems to be a safe bet too and is the default setting
    chunk_size_NPP= (/ mp, 1/) !This seems to be a safe bet too and is the default setting

    call check(&
      nf90_def_var(&
        ncid_par&
        ,"temperature"&
        ,NF90_DOUBLE&
        ,dimids_temperature&
        ,varid_temperature&
        ,chunksizes=chunk_size_temperature& 
      )&
    )
    call check(&
      nf90_def_var(&
        ncid_par&
        ,"NPP"&
        ,NF90_DOUBLE&
        ,dimids_NPP&
        ,varid_NPP&
        ,chunksizes=chunk_size_NPP& 
      )&
    )

    ! End define mode. This tells netCDF we are done defining
    ! metadata. This operation is collective and all processors will
    ! write their metadata to disk.
    call check(nf90_enddef(ncid_par))
    
    ! Unlimited dimensions require collective writes
    call check(nf90_var_par_access(ncid_par, varid_temperature, nf90_collective))
    call check(nf90_var_par_access(ncid_par, varid_NPP, nf90_collective))

    !)
    allocate(temperature_out(1:my_lp_number), stat = ierr)
    if (ierr .ne. 0) stop 3
    allocate(NPP_out(1:my_patch_number), stat = ierr)
    if (ierr .ne. 0) stop 3
    !count_temperature=(/ my_lp_number,1 /)
    !count_NPP=(/ my_patch_number,1 /)
    ! we now simulate a computation with timesteps
    ! this is the main loop
    do t=1,2
      ! in cable we would call some subroutines in the main loop
      ! Create some pretend data. We just need one row.
      ! we can either write the data here in the main loop
      ! or directly in the subroutines in which case we have to
      ! hand them the file id. 
      ! here we pretend that we can do the io in this main loop
      !start_temperature=(/ my_lp_index_offset+1,t /)
      do k = 1, my_lp_number
        temperature_out(k) = t*(k+my_lp_index_offset)
      end do
      
      ! Write the pretend data to the file. Each processor writes one row.
      call write_landpoint_var_netcdf(temperature_out,ncid_par,varid_temperature,t,worker_comm_id)
      
      !start_temperature=(/ my_lp_index_offset+1,t /)
      do k = 1, my_patch_number
        NPP_out(k) = t*(k+my_patch_ind_offset)
      end do
      !call check(&
      !  nf90_put_var(&
      !    ncid_par&
      !    ,varid_NPP&
      !    ,NPP_out& 
      !    ,start = (/ my_lp_index_offset+1,t /)&
      !    ,count = (/ my_patch_number,1 /)&
      !  )&
      !)
      call write_patch_var_netcdf(NPP_out,ncid_par,varid_NPP,t,worker_comm_id)
    enddo
    ! Close the file. This frees up any internal netCDF resources
    ! associated with the file, and flushes any buffers.
    call check( nf90_close(ncid_par) )

    ! Free my local memory.
    deallocate(temperature_out)
    deallocate(NPP_out)
  endif

  ! MPI library must be shut down.
  call MPI_Finalize(ierr)

  if (rank_in_world .eq. 0) then
    print *, "*** SUCCESS writing example files ", SERIAL_FILE_NAME,&
      " and ",PARALLEL_FILE_NAME,"! "
  endif 

contains
  subroutine check(status)
    integer, intent ( in) :: status
    
    if(status /= nf90_noerr) then 
      print *, trim(nf90_strerror(status))
      stop 2
    end if
  end subroutine check  
  
  subroutine write_landpoint_var_netcdf(lpvar,file_id,varid,t,worker_comm_id)
    integer,intent(in)            ::file_id,varid,t,worker_comm_id
    real(kind=8),intent(in),dimension(:) ::lpvar
    integer                       ::my_lp_index_offset ,my_lp_number,mland,mp,rank_in_worker_comm,wnp
    ! get cable like information about which worker owns which landpoints and
    ! patches This is not necessary in CABLE since the workers are told about this by the master
    call get_offsets_and_sizes(my_lp_index_offset ,my_lp_number&
      ,my_patch_ind_offset ,my_patch_number ,mland,mp,rank_in_worker_comm,wnp,worker_comm_id)
    call check(&
        nf90_put_var(&
          file_id&
          ,varid&
          ,lpvar& 
          ,start=(/ my_lp_index_offset+1,t /)&
          ,count = (/ my_lp_number,1 /)&
        )&
    )
  end subroutine 

  subroutine write_patch_var_netcdf(patch_var,file_id,varid,t,worker_comm_id)
    integer,intent(in)            ::file_id,varid,t,worker_comm_id
    real(kind=8),intent(in),dimension(:) ::patch_var
    integer                       ::my_lp_index_offset ,my_lp_number,mland,mp,rank_in_worker_comm,wnp
    ! get cable like information about which worker owns which landpoints and
    ! patches This is not necessary in CABLE since the workers are told about this by the master
    call get_offsets_and_sizes(my_lp_index_offset ,my_lp_number&
      ,my_patch_ind_offset ,my_patch_number ,mland,mp,rank_in_worker_comm,wnp,worker_comm_id)
    call sleep(rank_in_worker_comm)
    print*,'#################################'
    print*,'rank_in_worker_comm',rank_in_worker_comm
    print*,'count_temperature=',count_temperature
    print*,'chunk_size_temperature=',chunk_size_temperature
    print*,'my_patch_ind_offset',my_patch_ind_offset
    print*,'my_patch_number',my_patch_number
    print*,'mland',mland
    print*,'mp',mp
    print*,'patch_var',patch_var
    call check(&
        nf90_put_var(&
          file_id&
          ,varid&
          ,patch_var& 
          ,start=(/ my_patch_ind_offset+1,t /)&
          ,count = (/ my_patch_number,1 /)&
        )&
    )
  end subroutine 

  subroutine get_offsets_and_sizes(&
    my_lp_index_offset&
    ,my_lp_number&
    ,my_patch_ind_offset&
    ,my_patch_number&
    ,mland&
    ,mp&
    ,rank_in_worker_comm&
    ,wnp&
    ,worker_comm_id&
  )
    ! this is a pretend subroutine that mimics the behaviour in cable
    ! where the number of patches/per landpoint is input from a file
    ! which is eveluated by the master process and send to the workers
    ! So every worker gets a message telling it how many patches it 
    ! is responsible for and where the first of these is (offset).
    ! and the same for the landpoints. 
    ! While most subroutines will write out data patchwise some will write
    ! landpoint wise
    ! all this information can be inferred from the patches array which 
    ! represents the cable input file the number of workers determined by
    ! the n parameter with which mpirun was startded and the
    ! rank_in_worker_comm of the worker as returned by MPI_COMM_RANK at runtime 
    ! possibly reduced by one if rank 0 is the master.
    use mpi
    INTEGER,intent(out):: mland,mp
    integer,intent(out):: my_lp_index_offset ,my_lp_number,my_patch_ind_offset ,my_patch_number 
    integer,intent(out):: rank_in_worker_comm,wnp 
    integer,intent(in) :: worker_comm_id
    integer rest,lpw,wk,ierr
    integer,parameter,dimension(10)  :: ps_per_lp= (/ 1,1,1,3,2,1,1,3,4,1/)
    integer,allocatable,dimension(:)::lps&
                                      ,lp_index_offsets&
                                      ,ps_per_process&
                                      ,ps_index_offsets
    ! only workers can get their rank and the size of their communicator
    call MPI_Comm_rank(worker_comm_id, rank_in_worker_comm, ierr)
    call MPI_Comm_size(worker_comm_id, wnp, ierr)
    allocate(lps(0:wnp-1))
    allocate(lp_index_offsets(0:wnp-1))
    allocate(ps_per_process(0:wnp-1))
    allocate(ps_index_offsets(0:wnp-1))
    mland=size(ps_per_lp)
    mp=sum(ps_per_lp)
    
    rest=mod(mland,wnp)
    lpw=mland/wnp
    if (rank_in_worker_comm .eq. 0) then
      print *,"wnp: ", wnp
      print *,"rest: ", rest
      print *,"lpw: ", lpw
      print *,"ps_per_lp: ", ps_per_lp
    endif 
    ! this routine is called by every worker but it always has to consider
    ! what would be assinged to the other workers with lesser index since
    ! the offsets are cumulative
    ! first compute the landpoint load balancing
    do wk=0,(wnp-1)
      if (rest==0) then
        lps(wk)=lpw
      else
        lps(wk)=lpw+1
        rest=rest-1
      endif
    end do
    lp_index_offsets(0)=0
    do i=1,wnp-1
      lp_index_offsets(i)=sum(lps(0:(i-1)))
    enddo
    ! now we do the same for the patches
    do wk=0,(wnp-1)
      ps_per_process(wk)=sum(&
        ps_per_lp(&
          (lp_index_offsets(wk)+1):(lp_index_offsets(wk)+lps(wk))&
        )&
      )
    enddo
    ps_index_offsets(0)=0
    do wk=1,wnp-1
      ps_index_offsets(wk)=sum(ps_per_process(0:wk-1))
    enddo
    
    if (rank_in_worker_comm .eq. 0) then
      do wk=0,(wnp-1)
        print *, '##################################'
        print *, 'wk',wk
        print *, 'lps(wk)',lps(wk)
        print *, 'lp_index_offsets(wk)',lp_index_offsets(wk)
        print *, 'ps_per_process(wk)',ps_per_process(wk)
        print *, 'ps_index_offsets(wk)',ps_index_offsets(wk)
      end do
      endif
      my_lp_index_offset=lp_index_offsets(rank_in_worker_comm)
      my_lp_number=lps(rank_in_worker_comm)
      my_patch_ind_offset=ps_index_offsets(rank_in_worker_comm)
      my_patch_number=ps_per_process(rank_in_worker_comm)
     
  end subroutine get_offsets_and_sizes
    
  !function i_am_a_worker (world_rank)
  !  integer :: world_rank 
  !  logical :: i_am_a_worker
  !  i_am_a_worker = world_rank .ne. 0
  !  RETURN 
  !end 
    
end program simple_xy_par_wr

program INC2FV3
    
    use netcdf
    use machine,                  only: kind_phys, kind_dyn

    implicit none

    include 'mpif.h'

    Integer                 :: n_S, n_b
    INTEGER, ALLOCATABLE    :: row_esmf(:), col_esmf(:), mask_b(:) 
    REAL, ALLOCATABLE       :: S_esmf(:), frac_b(:)
    real, allocatable    :: stc_inc_gauss(:, :, :, :), slc_inc_gauss(:, :, :, :)  !t,k,y,x
    Real, ALLOCATABLE    :: stc_inc_reg(:,:, :), slc_inc_reg(:, :, :), obs_in(:), obs_ref(:)

    character(len=32), dimension(4) :: stc_vars = [character(len=32) :: 'soilt1_inc', 'soilt2_inc', 'soilt3_inc', 'soilt4_inc']
    character(len=32), dimension(4) :: slc_vars = [character(len=32) :: 'slc1_inc', 'slc2_inc', 'slc3_inc', 'slc4_inc']
    character(len=32) :: slsn_mask = "soilsnow_mask"
    ! stc_inc_vars = ["soilt1_inc", "soilt2_inc", "soilt3_inc", "soilt4_inc"]
    ! slc_inc_vars = ["slc1_inc", "slc2_inc", "slc3_inc", "slc4_inc"]

    character(len=32), dimension(3) :: gaussian_inc_files = [character(len=32) :: "fv_sfc_increment3.nc", "fv_sfc_increment6.nc", "fv_sfc_increment9.nc"]
    character(len=32) :: fv3_inc_prefix = "sfc_inc"
    real, dimension(3) :: time_list = [3.0, 6.0, 9.0]
    integer, dimension(3) :: fhr_int = [3, 6, 9]
    integer  :: nt = 3
    integer  :: nk = 4
    integer  :: ny = 48
    integer  :: nx = 48
    integer  :: xg, yg
    integer  :: it, k, j, y, x
    integer :: ierr, nprocs, myrank
    integer              :: nlunit

    character(len=32)    :: weight_file = "wgf48.nc"   != inc_dir + "wgf48.nc"       
    character(len=128)   :: gaussian_sfc_inc_prefix = "sfcincr_2021122100_fhr0"  !sfcincr_2021122100_fhr06_mem006
    character(len=128)   :: fv3_sfc_inc_prefix = "sfc_inc_2021122100_mem"           !sfc_inc_2021122100_mem006_tile6.nc
    integer              :: ens_size = 20

    NAMELIST/NAMSNO/ weight_file, gaussian_sfc_inc_prefix, fv3_sfc_inc_prefix, ens_size 

    CALL MPI_INIT(IERR)
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, NPROCS, IERR)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD, MYRANK, IERR)

    PRINT*,"starting Convert Gaussian Increment to FV3 ON RANK ", MYRANK, " RUNNING WITH ", NPROCS, "TASKS"
    
    IF (MYRANK==0) PRINT*,"READING NAMSNO NAMELIST."

    nlunit=23
    open (unit=nlunit, file='fort.3600', READONLY, status='OLD')
    read(nlunit, NAMSNO)
    close(nlunit)
    IF (MYRANK==0) WRITE(6, NAMSNO)

    if (myrank+1 > ens_size) then
        print*, "Not running on proc ", myrank
        goto 999
    endif

    IF (MYRANK==0) print*, " Start reading esmf weights"
    
    call read_back_esmf_weights(weight_file, n_S, row_esmf, col_esmf, mask_b, S_esmf, frac_b)
    print*, "Finished reading ESMF weights"

    n_b = size(frac_b)
    if (n_b .ne. (6 * ny * nx)) then
        print*, "proc ", myrank, " length of target array nb ", n_b, " not equal to num tiles * ny * nx ", 6*ny*nx
        STOP
    endif

    allocate(stc_inc_reg(nt, nk, n_b))
    allocate(slc_inc_reg(nt, nk, n_b))
    allocate(obs_in(n_S))
    allocate(obs_ref(n_b))

    do irank = myrank+1, ens_size, nprocs 

        print*, "proc ", myrank, " loop ", irank

        call read_increments_atrank(irank, nt, nk, fhr_int, gaussian_sfc_inc_prefix, stc_inc_gauss, slc_inc_gauss, xg, yg)
        ! call read_increments(nt, nk, gaussian_inc_files, stc_inc_gauss, slc_inc_gauss, xg, yg)
        print*, "proc ", myrank, " loop ", irank, " Finished reading increments starting interpolations"
        do it = 1, nt
            do k = 1, 4
                ! print*, "it = ", it, " k = ", k

                obs_in = reshape(stc_inc_gauss(it, k, :, :), [xg * yg])  !(/xg*yg/)
                obs_ref = 0.0
                do j = 1, n_S
                    obs_ref(row_esmf(j)) = obs_ref(row_esmf(j)) + S_esmf(j) * obs_in(col_esmf(j))
                enddo
                do j = 1, n_b
                    if ((frac_b(j) .ne. 0) .and. (frac_b(j) .ne. 1)) then 
                        ! print*, stc_inc_vars[k], " frac adj applied"
                        obs_ref(j) = obs_ref(j) / frac_b(j)
                    endif
                enddo
                stc_inc_reg(it, k, :) = obs_ref(:)

                obs_in = reshape(slc_inc_gauss(it, k, :, :), [xg * yg])
                obs_ref = 0.0
                do j = 1, n_S
                    obs_ref(row_esmf(j)) = obs_ref(row_esmf(j)) + S_esmf(j) * obs_in(col_esmf(j))
                enddo
                do j = 1, n_b
                    if ((frac_b(j) .ne. 0) .and. (frac_b(j) .ne. 1)) then 
                        ! print*, stc_inc_vars[k], " frac adj applied"
                        obs_ref(j) = obs_ref(j) / frac_b(j)
                    endif
                enddo
                slc_inc_reg(it, k, :) = obs_ref(:)
            enddo
        enddo

        print*, "proc ", myrank, " loop ", irank, " Finished interpolation, starting writing regridded data"

        call write_regridded_inc_atrank(irank, fv3_sfc_inc_prefix, nt, nk, ny, nx, n_b, stc_inc_reg, slc_inc_reg)
        ! call write_regridded_inc(fv3_sfc_inc_prefix, nt, nk, ny, nx, n_b, stc_inc_reg, slc_inc_reg)

        print*, "proc ", myrank, " loop ", irank, " deallocating memory"  
        
        if (allocated(stc_inc_gauss)) deallocate(stc_inc_gauss)
        if (allocated(slc_inc_gauss)) deallocate(slc_inc_gauss)

    enddo
    
    deallocate(stc_inc_reg, slc_inc_reg, obs_in, obs_ref)
    if (allocated(row_esmf)) deallocate(row_esmf)
    if (allocated(col_esmf)) deallocate(col_esmf)
    if (allocated(mask_b)) deallocate(mask_b)
    if (allocated(S_esmf)) deallocate(S_esmf)
    if (allocated(frac_b)) deallocate(frac_b)

    print*, "Done on proc ", myrank

999 continue

    CALL MPI_FINALIZE(IERR)

    STOP

  contains
    
    subroutine read_increments_atrank(irank, nt, nk, gaussian_inc_prefix, stc_inc_gauss, slc_inc_gauss, xg, yg)

        implicit none
        
        integer :: irank, nt, nk
        integer, dimension(nt) :: fhr_int
        character(len=*)    :: gaussian_inc_prefix  ! = [character(len=32) :: "fv_sfc_increment3.nc", "fv_sfc_increment6.nc", "fv_sfc_increment9.nc"]
        real, allocatable    :: stc_inc_gauss(:, :, :, :), slc_inc_gauss(:, :, :, :)
        integer  :: xg, yg

        character(len=500)    :: gaussian_inc_file
        character(len=1)    :: tile_str
        character(len=3)    :: mem_str
        character(len=2)    :: fhr_str
        character(len=120)  :: file_out
        logical  :: exists
        integer  :: ncid, status, varid
        integer  :: ierr
        integer  :: i, it, k, tl

        character(len=32), dimension(4) :: stc_vars = [character(len=32) :: 'soilt1_inc', 'soilt2_inc', 'soilt3_inc', 'soilt4_inc']
        character(len=32), dimension(4) :: slc_vars = [character(len=32) :: 'slc1_inc', 'slc2_inc', 'slc3_inc', 'slc4_inc']

        character(len=500) :: errmsg
        integer          :: errflg
        
        write(fhr_str, "(I2.2)") fhr_int(1)
        write(mem_str, '(I3.3)') irank
        gaussian_inc_file = trim(gaussian_inc_prefix)//fhr_str//"_mem"//mem_str       !sfcincr_2021122100_fhr09_mem008

        inquire (file=trim(gaussian_inc_file), exist=exists)    
        if (exists) then
            status = nf90_open(trim(gaussian_inc_file), NF90_NOWRITE, ncid)  ! open the file
            call netcdf_err(status, ' opening file '//trim(gaussian_inc_file), errflg, errmsg) 
            if (errflg .ne. 0) then 
                print*, trim(errmsg)
                stop
            endif
        else
            errmsg = 'FATAL Error Expected file '//trim(gaussian_inc_file)//' for DA increment does not exist'
            errflg = 1
            stop
        endif

        call get_nc_dimlen(ncid, "longitude", xg, errflg, errmsg) 
        if (errflg .ne. 0) then 
            print*, trim(errmsg)
            stop
        endif
        call get_nc_dimlen(ncid, "latitude", yg, errflg, errmsg) 
        if (errflg .ne. 0) then 
            print*, trim(errmsg)
            stop
        endif

        status =nf90_close(ncid) 
        call netcdf_err(status, 'closing file '//trim(gaussian_inc_file), errflg, errmsg)
        if (errflg .ne. 0) then 
            print*, trim(errmsg)
            stop
        endif

        allocate(stc_inc_gauss(nt, nk, xg, yg))
        allocate(slc_inc_gauss(nt, nk, xg, yg))

        do it = 1, nt

            write(fhr_str, "(I2.2)") fhr_int(it)
            gaussian_inc_file = trim(gaussian_inc_prefix)//fhr_str//"_mem"//mem_str

            inquire (file=trim(gaussian_inc_file), exist=exists)    
            if (exists) then
                status = nf90_open(trim(gaussian_inc_file), NF90_NOWRITE, ncid)  ! open the file
                call netcdf_err(status, ' opening file '//trim(gaussian_inc_file), errflg, errmsg) 
                if (errflg .ne. 0) then 
                    print*, trim(errmsg)
                    stop
                endif
            else
                errmsg = 'FATAL Error Expected file '//trim(gaussian_inc_file)//' for DA increment does not exist'
                errflg = 1
                stop
            endif

            do k = 1, nk
                status = nf90_inq_varid(ncid, stc_vars(k), varid)
                call netcdf_err(status, ' getting varid for '//trim(stc_vars(k))//' from file '//trim(gaussian_inc_file), errflg, errmsg) 
                if (errflg .ne. 0) then 
                    print*, trim(errmsg)
                    stop
                endif
                ! soilt1_inc(latitude, longitude)
                status = nf90_get_var(ncid, varid, stc_inc_gauss(it, k, :, :),  start = (/1, 1/), count = (/xg, yg/))
                call netcdf_err(status, ' reading data for '//trim(stc_vars(k))//' from file '//trim(gaussian_inc_file)), errflg, errmsg) 
                if (errflg .ne. 0) then 
                    print*, trim(errmsg)
                    stop
                endif

                status = nf90_inq_varid(ncid, slc_vars(k), varid)
                call netcdf_err(status, ' getting varid for '//trim(slc_vars(k))//' from file '//trim(gaussian_inc_file), errflg, errmsg) 
                if (errflg .ne. 0) then 
                    print*, trim(errmsg)
                    stop
                endif
                ! soilt1_inc(latitude, longitude)
                status = nf90_get_var(ncid, varid, slc_inc_gauss(it, k, :, :),  start = (/1, 1/), count = (/xg, yg/))
                call netcdf_err(status, ' reading data for '//trim(slc_vars(k))//' from file '//trim(gaussian_inc_file), errflg, errmsg) 
                if (errflg .ne. 0) then 
                    print*, trim(errmsg)
                    stop
                endif
            enddo

            status =nf90_close(ncid) 
            call netcdf_err(status, 'closing file '//trim(gaussian_inc_file), errflg, errmsg)
            if (errflg .ne. 0) then 
                print*, trim(errmsg)
                stop
            endif

        enddo

    end subroutine read_increments_atrank

    subroutine write_regridded_inc_atrank(irank, fv3_sfc_inc_prefix, nt, nk, ny, nx, nb, stc_inc_reg, slc_inc_reg)

        implicit none
        
        integer           :: irank
        character(len=*)  :: fv3_sfc_inc_prefix
        integer :: nt, nk, ny, nx, nb
        real    :: stc_inc_reg(nt, nk, nb), slc_inc_reg(nt, nk, nb)
        
        character(len=1)    :: tile_str
        character(len=3)    :: mem_str
        character(len=500)  :: file_out
        logical  :: exists
        integer  :: ncid, status, varid
        integer  :: ierr
        integer  :: it, k, tl
        real     :: val_tile (ny*nx)

        character(len=500) :: errmsg
        integer          :: errflg

        character(len=32), dimension(4) :: stc_vars = [character(len=32) :: 'soilt1_inc', 'soilt2_inc', 'soilt3_inc', 'soilt4_inc']
        character(len=32), dimension(4) :: slc_vars = [character(len=32) :: 'slc1_inc', 'slc2_inc', 'slc3_inc', 'slc4_inc']
 
        write(mem_str, '(I3.3)') irank

        do tl = 1, 6
            print*, "tile ", tl

            write(tile_str, '(I0)') tl
            file_out = trim(fv3_sfc_inc_prefix)//mem_str//"_tile"//tile_str     !//".nc"    !sfc_inc_2021122100_mem006_tile6.nc

            inquire (file=trim(file_out), exist=exists)    
            if (exists) then
                status = nf90_open(trim(file_out), NF90_WRITE, ncid)  ! open the file
                call netcdf_err(status, ' opening file '//trim(file_out), errflg, errmsg) 
                if (errflg .ne. 0) then 
                    print*, trim(errmsg)
                    stop
                endif
            else
                print*, 'FATAL Error in INC2FV3 Expected file '//trim(file_out)//' for DA increment does not exist'
                return
            endif
            
            do k = 1, nk
                print*, "k = ", k

                ! ncOut.createVariable('soilt1_inc', 'f4', ('Time', 'yaxis_1', 'xaxis_1',), fill_value=9.96921e+36)
                status = nf90_inq_varid(ncid, stc_vars(k), varid)
                call netcdf_err(status, ' getting varid for '//trim(stc_vars(k)), errflg, errmsg) 
                if (errflg .ne. 0) then 
                    print*, trim(errmsg)
                    stop
                endif
                
                do it = 1, nt
                    val_tile = stc_inc_reg(it, k, (tl-1)*ny*nx+1:tl*ny*nx)
                    ! var stored as soilt1_inc(Time, yaxis_1, xaxis_1)
                    status = nf90_put_var(ncid, varid , reshape(val_tile, [nx, ny]), start = (/1, 1, it/), count = (/nx, ny, 1/))
                    call netcdf_err(status, ' writing array for '//trim(stc_vars(k)), errflg, errmsg) 
                    if (errflg .ne. 0) then 
                        print*, trim(errmsg), 'time step ', it
                        stop
                    endif
                enddo

                ! ncOut.createVariable('soilt1_inc', 'f4', ('Time', 'yaxis_1', 'xaxis_1',), fill_value=9.96921e+36)
                status = nf90_inq_varid(ncid, slc_vars(k), varid)
                call netcdf_err(status, ' getting varid for '//trim(slc_vars(k)), errflg, errmsg) 
                if (errflg .ne. 0) then 
                    print*, trim(errmsg)
                    stop
                endif
                
                do it = 1, nt
                    print*, "it = ", it

                    val_tile = slc_inc_reg(it, k, (tl-1)*ny*nx+1:tl*ny*nx)
                    ! var stored as soilt1_inc(Time, yaxis_1, xaxis_1)
                    status = nf90_put_var(ncid, varid , reshape(val_tile, [nx, ny]), start = (/1, 1, it/), count = (/nx, ny, 1/))
                    call netcdf_err(status, ' writing array for '//trim(slc_vars(k)), errflg, errmsg) 
                    if (errflg .ne. 0) then 
                        print*, trim(errmsg), 'time step ', it
                        stop
                    endif
                enddo

            enddo

            status =nf90_close(ncid) 
            call netcdf_err(status, 'closing file '//trim(file_out), errflg, errmsg)
            if (errflg .ne. 0) then 
                print*, trim(errmsg)
                stop
            endif 

        enddo
                

    end subroutine write_regridded_inc_atrank

    subroutine read_increments(nt, nk, gaussian_inc_files, stc_inc_gauss, slc_inc_gauss, xg, yg)

        implicit none
        
        integer :: nt, nk
        character(len=32), dimension(nt) :: gaussian_inc_files  ! = [character(len=32) :: "fv_sfc_increment3.nc", "fv_sfc_increment6.nc", "fv_sfc_increment9.nc"]
        real, allocatable    :: stc_inc_gauss(:, :, :, :), slc_inc_gauss(:, :, :, :)
        integer  :: xg, yg

        character(len=1)    :: tile_str
        character(len=120)  :: file_out
        logical  :: exists
        integer  :: ncid, status, varid
        integer  :: ierr
        integer  :: i, it, k, tl

        character(len=32), dimension(4) :: stc_vars = [character(len=32) :: 'soilt1_inc', 'soilt2_inc', 'soilt3_inc', 'soilt4_inc']
        character(len=32), dimension(4) :: slc_vars = [character(len=32) :: 'slc1_inc', 'slc2_inc', 'slc3_inc', 'slc4_inc']

        character(len=500) :: errmsg
        integer          :: errflg

        inquire (file=trim(gaussian_inc_files(1)), exist=exists)    
        if (exists) then
            status = nf90_open(trim(gaussian_inc_files(1)), NF90_NOWRITE, ncid)  ! open the file
            call netcdf_err(status, ' opening file '//trim(gaussian_inc_files(1)), errflg, errmsg) 
            if (errflg .ne. 0) then 
                print*, trim(errmsg)
                stop
            endif
        else
            errmsg = 'FATAL Error Expected file '//trim(gaussian_inc_files(1))//' for DA increment does not exist'
            errflg = 1
            stop
        endif

        call get_nc_dimlen(ncid, "longitude", xg, errflg, errmsg) 
        if (errflg .ne. 0) then 
            print*, trim(errmsg)
            stop
        endif
        call get_nc_dimlen(ncid, "latitude", yg, errflg, errmsg) 
        if (errflg .ne. 0) then 
            print*, trim(errmsg)
            stop
        endif

        status =nf90_close(ncid) 
        call netcdf_err(status, 'closing file '//trim(gaussian_inc_files(1)), errflg, errmsg)
        if (errflg .ne. 0) then 
            print*, trim(errmsg)
            stop
        endif

        allocate(stc_inc_gauss(nt, nk, xg, yg))
        allocate(slc_inc_gauss(nt, nk, xg, yg))

        do it = 1, nt
            inquire (file=trim(gaussian_inc_files(it)), exist=exists)    
            if (exists) then
                status = nf90_open(trim(gaussian_inc_files(it)), NF90_NOWRITE, ncid)  ! open the file
                call netcdf_err(status, ' opening file '//trim(gaussian_inc_files(it)), errflg, errmsg) 
                if (errflg .ne. 0) then 
                    print*, trim(errmsg)
                    stop
                endif
            else
                errmsg = 'FATAL Error Expected file '//trim(gaussian_inc_files(it))//' for DA increment does not exist'
                errflg = 1
                stop
            endif

            do k = 1, nk
                status = nf90_inq_varid(ncid, stc_vars(k), varid)
                call netcdf_err(status, ' getting varid for '//trim(stc_vars(k))//' from file '//trim(gaussian_inc_files(it)), errflg, errmsg) 
                if (errflg .ne. 0) then 
                    print*, trim(errmsg)
                    stop
                endif
                ! soilt1_inc(latitude, longitude)
                status = nf90_get_var(ncid, varid, stc_inc_gauss(it, k, :, :),  start = (/1, 1/), count = (/xg, yg/))
                call netcdf_err(status, ' reading data for '//trim(stc_vars(k))//' from file '//trim(gaussian_inc_files(it)), errflg, errmsg) 
                if (errflg .ne. 0) then 
                    print*, trim(errmsg)
                    stop
                endif

                status = nf90_inq_varid(ncid, slc_vars(k), varid)
                call netcdf_err(status, ' getting varid for '//trim(slc_vars(k))//' from file '//trim(gaussian_inc_files(it)), errflg, errmsg) 
                if (errflg .ne. 0) then 
                    print*, trim(errmsg)
                    stop
                endif
                ! soilt1_inc(latitude, longitude)
                status = nf90_get_var(ncid, varid, slc_inc_gauss(it, k, :, :),  start = (/1, 1/), count = (/xg, yg/))
                call netcdf_err(status, ' reading data for '//trim(slc_vars(k))//' from file '//trim(gaussian_inc_files(it)), errflg, errmsg) 
                if (errflg .ne. 0) then 
                    print*, trim(errmsg)
                    stop
                endif
            enddo

            status =nf90_close(ncid) 
            call netcdf_err(status, 'closing file '//trim(gaussian_inc_files(it)), errflg, errmsg)
            if (errflg .ne. 0) then 
                print*, trim(errmsg)
                stop
            endif

        enddo

    end subroutine read_increments

    subroutine write_regridded_inc(inc_out, nt, nk, ny, nx, nb, stc_inc_reg, slc_inc_reg)

        implicit none
        
        character(len=*)  :: inc_out
        integer :: nt, nk, ny, nx, nb
        real    :: stc_inc_reg(nt, nk, nb), slc_inc_reg(nt, nk, nb)
        
        character(len=1)    :: tile_str
        character(len=120)  :: file_out
        logical  :: exists
        integer  :: ncid, status, varid
        integer  :: ierr
        integer  :: it, k, tl
        real     :: val_tile (ny*nx)

        character(len=500) :: errmsg
        integer          :: errflg

        character(len=32), dimension(4) :: stc_vars = [character(len=32) :: 'soilt1_inc', 'soilt2_inc', 'soilt3_inc', 'soilt4_inc']
        character(len=32), dimension(4) :: slc_vars = [character(len=32) :: 'slc1_inc', 'slc2_inc', 'slc3_inc', 'slc4_inc']


        do tl = 1, 6
            print*, "tile ", tl

            write(tile_str, '(I0)') tl
            file_out = trim(inc_out) // ".tile" // tile_str // ".nc"

            inquire (file=trim(file_out), exist=exists)    
            if (exists) then
                status = nf90_open(trim(file_out), NF90_WRITE, ncid)  ! open the file
                call netcdf_err(status, ' opening file '//trim(file_out), errflg, errmsg) 
                if (errflg .ne. 0) then 
                    print*, trim(errmsg)
                    stop
                endif
            else
                print*, 'FATAL Error in INC2FV3 Expected file '//trim(file_out)//' for DA increment does not exist'
                return
            endif
            
            do k = 1, nk
                print*, "k = ", k

                ! ncOut.createVariable('soilt1_inc', 'f4', ('Time', 'yaxis_1', 'xaxis_1',), fill_value=9.96921e+36)
                status = nf90_inq_varid(ncid, stc_vars(k), varid)
                call netcdf_err(status, ' getting varid for '//trim(stc_vars(k)), errflg, errmsg) 
                if (errflg .ne. 0) then 
                    print*, trim(errmsg)
                    stop
                endif
                
                do it = 1, nt
                    val_tile = stc_inc_reg(it, k, (tl-1)*ny*nx+1:tl*ny*nx)
                    ! var stored as soilt1_inc(Time, yaxis_1, xaxis_1)
                    status = nf90_put_var(ncid, varid , reshape(val_tile, [nx, ny]), start = (/1, 1, it/), count = (/nx, ny, 1/))
                    call netcdf_err(status, ' writing array for '//trim(stc_vars(k)), errflg, errmsg) 
                    if (errflg .ne. 0) then 
                        print*, trim(errmsg), 'time step ', it
                        stop
                    endif
                enddo

                ! ncOut.createVariable('soilt1_inc', 'f4', ('Time', 'yaxis_1', 'xaxis_1',), fill_value=9.96921e+36)
                status = nf90_inq_varid(ncid, slc_vars(k), varid)
                call netcdf_err(status, ' getting varid for '//trim(slc_vars(k)), errflg, errmsg) 
                if (errflg .ne. 0) then 
                    print*, trim(errmsg)
                    stop
                endif
                
                do it = 1, nt
                    print*, "it = ", it

                    val_tile = slc_inc_reg(it, k, (tl-1)*ny*nx+1:tl*ny*nx)
                    ! var stored as soilt1_inc(Time, yaxis_1, xaxis_1)
                    status = nf90_put_var(ncid, varid , reshape(val_tile, [nx, ny]), start = (/1, 1, it/), count = (/nx, ny, 1/))
                    call netcdf_err(status, ' writing array for '//trim(slc_vars(k)), errflg, errmsg) 
                    if (errflg .ne. 0) then 
                        print*, trim(errmsg), 'time step ', it
                        stop
                    endif
                enddo

            enddo

            status =nf90_close(ncid) 
            call netcdf_err(status, 'closing file '//trim(file_out), errflg, errmsg)
            if (errflg .ne. 0) then 
                print*, trim(errmsg)
                stop
            endif 

        enddo
                

    end subroutine write_regridded_inc

            ! assumes standard esmf names col(n_s) ;
    SUBROUTINE read_back_esmf_weights(inp_file, dim_size, row_esmf, col_esmf, mask_b, S_esmf, frac_b)
        
        IMPLICIT NONE    

        CHARACTER(LEN=*), Intent(In)      :: inp_file
        INTEGER, Intent(Out)              :: dim_size
        INTEGER, ALLOCATABLE, Intent(Out)    :: row_esmf(:), col_esmf(:), mask_b(:)
        REAL, ALLOCATABLE, Intent(Out)    :: S_esmf(:), frac_b(:)
    
        INTEGER                :: ERROR, NCID, grp_ncid, ID_DIM, ID_VAR, n_b
        LOGICAL                :: file_exists
        character(len=32)      :: dim_name = "n_s"
        character(len=500) :: errmsg
        integer          :: errflg

        INQUIRE(FILE=trim(inp_file), EXIST=file_exists)
        if (.not. file_exists) then
                print *, 'read_esmf_weights error, file does not exist', trim(inp_file) 
                stop
        endif
    
        ERROR=NF90_OPEN(TRIM(inp_file),NF90_NOWRITE,NCID)
        CALL NETCDF_ERR(ERROR, 'OPENING FILE: '//TRIM(inp_file), errflg, errmsg )
        if (errflg .ne. 0) then 
            print*, trim(errmsg)
            stop
        endif 
      
        ERROR=NF90_INQ_DIMID(NCID, TRIM(dim_name), ID_DIM)
        CALL NETCDF_ERR(ERROR, 'ERROR READING Dimension'//trim(dim_name), errflg, errmsg )  
        if (errflg .ne. 0) then 
            print*, trim(errmsg)
            stop
        endif   
        ERROR=NF90_INQUIRE_DIMENSION(NCID,ID_DIM,LEN=dim_size)
        CALL NETCDF_ERR(ERROR, 'ERROR READING Size of Dimension', errflg, errmsg )
        if (errflg .ne. 0) then 
            print*, trim(errmsg)
            stop
        endif 

        ERROR=NF90_INQ_DIMID(NCID, 'n_b', ID_DIM)
        CALL NETCDF_ERR(ERROR, 'ERROR READING Dimension n_b', errflg, errmsg )   
        if (errflg .ne. 0) then 
            print*, trim(errmsg)
            stop
        endif  
        ERROR=NF90_INQUIRE_DIMENSION(NCID,ID_DIM,LEN=n_b)
        CALL NETCDF_ERR(ERROR, 'ERROR READING Size of Dimension n_b', errflg, errmsg ) 
        if (errflg .ne. 0) then 
            print*, trim(errmsg)
            stop
        endif    

        ALLOCATE(row_esmf(dim_size))
        ALLOCATE(col_esmf(dim_size))
        ALLOCATE(S_esmf(dim_size))
        
        ALLOCATE(mask_b(n_b))
        ALLOCATE(frac_b(n_b))

        ERROR=NF90_INQ_VARID(ncid, 'col', ID_VAR)
        CALL NETCDF_ERR(ERROR, 'ERROR READING col ID', errflg, errmsg )
        if (errflg .ne. 0) then 
            print*, trim(errmsg)
            stop
        endif 
        ERROR=NF90_GET_VAR(ncid, ID_VAR, col_esmf)
        CALL NETCDF_ERR(ERROR, 'ERROR READING col RECORD', errflg, errmsg )
        if (errflg .ne. 0) then 
            print*, trim(errmsg)
            stop
        endif 

        ERROR=NF90_INQ_VARID(ncid, 'row', ID_VAR)
        CALL NETCDF_ERR(ERROR, 'ERROR READING row ID', errflg, errmsg )
        if (errflg .ne. 0) then 
            print*, trim(errmsg)
            stop
        endif 
        ERROR=NF90_GET_VAR(ncid, ID_VAR, row_esmf)
        CALL NETCDF_ERR(ERROR, 'ERROR READING row RECORD', errflg, errmsg )
        if (errflg .ne. 0) then 
            print*, trim(errmsg)
            stop
        endif 

        ERROR=NF90_INQ_VARID(ncid, 'mask_b', ID_VAR)
        CALL NETCDF_ERR(ERROR, 'ERROR READING mask_b ID', errflg, errmsg )
        if (errflg .ne. 0) then 
            print*, trim(errmsg)
            stop
        endif 
        ERROR=NF90_GET_VAR(ncid, ID_VAR, mask_b)
        CALL NETCDF_ERR(ERROR, 'ERROR READING mask_b RECORD', errflg, errmsg )
        if (errflg .ne. 0) then 
            print*, trim(errmsg)
            stop
        endif 

        ERROR=NF90_INQ_VARID(ncid, 'S', ID_VAR)
        CALL NETCDF_ERR(ERROR, 'ERROR READING S ID', errflg, errmsg )
        if (errflg .ne. 0) then 
            print*, trim(errmsg)
            stop
        endif 
        ERROR=NF90_GET_VAR(ncid, ID_VAR, S_esmf)
        CALL NETCDF_ERR(ERROR, 'ERROR READING S RECORD', errflg, errmsg )
        if (errflg .ne. 0) then 
            print*, trim(errmsg)
            stop
        endif 

        ERROR=NF90_INQ_VARID(ncid, 'frac_b', ID_VAR)
        CALL NETCDF_ERR(ERROR, 'ERROR READING frac_b ID', errflg, errmsg )
        if (errflg .ne. 0) then 
            print*, trim(errmsg)
            stop
        endif 
        ERROR=NF90_GET_VAR(ncid, ID_VAR, frac_b)
        CALL NETCDF_ERR(ERROR, 'ERROR READING frac_b RECORD', errflg, errmsg )
        if (errflg .ne. 0) then 
            print*, trim(errmsg)
            stop
        endif 
    
        ERROR = NF90_CLOSE(NCID)
        CALL NETCDF_ERR(ERROR, 'ERROR closing file'//TRIM(inp_file), errflg, errmsg )
        if (errflg .ne. 0) then 
            print*, trim(errmsg)
            stop
        endif 
                  
        RETURN
        
     End SUBROUTINE read_back_esmf_weights 

     SUBROUTINE NETCDF_ERR(ERR, STRING, errflg, errmsg_out)

    !--------------------------------------------------------------
    ! IF AT NETCDF CALL RETURNS AN ERROR, PRINT OUT A MESSAGE
    ! AND STOP PROCESSING.
    !--------------------------------------------------------------
        IMPLICIT NONE
    
        INTEGER, INTENT(IN) :: ERR
        CHARACTER(LEN=*), INTENT(IN) :: STRING
        CHARACTER(LEN=80) :: ERRMSG
        integer :: errflg
        character(len=*) :: errmsg_out
    
        !Errors messages handled through CCPP error handling variables
        errmsg_out = ''
        errflg = 0
    
        IF (ERR == NF90_NOERR) RETURN
        ERRMSG = NF90_STRERROR(ERR)
        ! PRINT*,'FATAL ERROR in INC2FV3 ', TRIM(STRING), ': ', TRIM(ERRMSG)
        errmsg_out = 'FATAL ERROR in INC2FV3 '//TRIM(STRING)//': '//TRIM(ERRMSG)
        errflg = 1
        return
    
    END SUBROUTINE NETCDF_ERR
    
    subroutine get_nc_dimlen(ncid, dim_name, dim_len, errflg, errmsg_out )
        integer, intent(in):: ncid
        character(len=*), intent(in)::  dim_name
        integer, intent(out):: dim_len
        integer :: dimid
        integer :: errflg
        character(len=*) :: errmsg_out
        integer :: status
    
        !Errors messages handled through CCPP error handling variables
        errmsg_out = ''
        errflg = 0
    
        status = nf90_inq_dimid(ncid, dim_name, dimid)
        CALL netcdf_err(status, 'reading dim id '//trim(dim_name), errflg, errmsg_out)
        if (errflg .ne. 0) return
        status = nf90_inquire_dimension(ncid, dimid, len = dim_len)
        CALL netcdf_err(status, 'reading dim length '//trim(dim_name), errflg, errmsg_out)
    
    end subroutine get_nc_dimlen
        ! status = nf90_inq_dimid(ncid, "longitude", dimid)
        ! CALL netcdf_err(status, 'reading longitude dim id')
        ! status = nf90_inquire_dimension(ncid, dimid, len = im)
        ! CALL netcdf_err(status, 'reading dim longitude')
        ! status = nf90_inq_dimid(ncid, "latitude", dimid)
        ! CALL netcdf_err(status, 'reading latitude dim id')
        ! status = nf90_inquire_dimension(ncid, dimid, len = jm)
        ! CALL netcdf_err(status, 'reading dim latitude')
    subroutine get_var1d(ncid, dim_len, var_name, var_arr, errflg, errmsg_out)
        integer, intent(in):: ncid, dim_len
        character(len=*), intent(in)::  var_name
        real(kind=kind_phys), intent(out):: var_arr(dim_len)
        integer :: errflg
        character(len=*) :: errmsg_out
        integer :: varid, status
    
        !Errors messages handled through CCPP error handling variables
        errmsg_out = ''
        errflg = 0
    
        status = nf90_inq_varid(ncid, trim(var_name), varid)
        CALL NETCDF_ERR(status, 'getting varid: '//trim(var_name), errflg, errmsg_out)
        if (errflg .ne. 0) return
        status = nf90_get_var(ncid, varid, var_arr)
                    ! start = (/1/), count = (/dim_len/))
        CALL NETCDF_ERR(status, 'reading var: '//trim(var_name), errflg, errmsg_out)
    
    end subroutine get_var1d
    
    subroutine get_var3d_values(ncid, varid, is,ix, js,jy, ks,kz, var3d, status)
        integer, intent(in):: ncid, varid
        integer, intent(in):: is, ix, js, jy, ks,kz
        real(kind=kind_phys), intent(out):: var3d(ix, jy, kz)   !var3d(is:ie,js:je,ks:ke)
        integer, intent(out):: status 
        ! integer, dimension(3):: start, nreco
        ! start(1) = is; start(2) = js; start(3) = ks
        ! nreco(1) = ie - is + 1
        ! nreco(2) = je - js + 1
        ! nreco(3) = ke - ks + 1
    
        status = nf90_get_var(ncid, varid, var3d, &  !start = start, count = nreco)
                start = (/is, js, ks/), count = (/ix, jy, kz/))
                ! start = (/is, js, ks/), count = (/ie - is + 1, je - js + 1, ke - ks + 1/))
    
    end subroutine get_var3d_values
    
    subroutine get_var3d_values_int(ncid, varid, is,ix, js,jy, ks,kz, var3d, status)
        integer, intent(in):: ncid, varid
        integer, intent(in):: is, ix, js, jy, ks,kz
        integer, intent(out):: var3d(ix, jy, kz)   !var3d(is:ie,js:je,ks:ke)
        integer, intent(out):: status 
        ! integer, dimension(3):: start, nreco
        ! start(1) = is; start(2) = js; start(3) = ks
        ! nreco(1) = ie - is + 1
        ! nreco(2) = je - js + 1
        ! nreco(3) = ke - ks + 1
    
        status = nf90_get_var(ncid, varid, var3d, &  !start = start, count = nreco)
                start = (/is, js, ks/), count = (/ix, jy, kz/))
                ! start = (/is, js, ks/), count = (/ie - is + 1, je - js + 1, ke - ks + 1/))
    
    end subroutine get_var3d_values_int
    

end program INC2FV3

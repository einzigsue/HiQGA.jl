!-----------------------------------------------------------------------------------------------------------------------------------
    module mare2dem_lib_interface
!
! Interface for calling MARE2DEM from external codes for forward computations
!
! General sequence of calling from external routine:
!    (1) mare2dem_load_files                -  loads mare2dem format data and model files
!    (2) mare2dem_forward(freeResistivity)  - computes forward response and misfit for freeResistivity array of free parameter vals
!
! Only need to do step (1) once at the start. Then step (2) can be done any number of times for varying freeResistivity arrays.
!
! Optional:
!    mare2dem_get_params() can be called after (1) to get y,z of free parameter centroids.
!
!-----------------------------------------------------------------------------------------------------------------------------------

    use occam
    use mare2d_mpi_definitions    ! for MPI communicator tags etc
    use mare2dem_global           ! for setting defaults and storing worker status array
    use mare2dem_io
    !use mare2dem_input_data_params

    use, intrinsic :: iso_fortran_env, only: int64, real64

    implicit none

    public :: mare2dem_load_files, mare2dem_forward

    contains

!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------

   subroutine say_hello_pids(pid) bind(C, name="say_hello_pids")

   integer(int64), intent(in)   :: pid

   write(*,*) 'hello, I am pid #',pid

   end subroutine say_hello_pids

!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------

    subroutine mare2dem_load_files_Test(nFreeRegions,nParamsPerRegion) bind(C, name="mare2dem_load_files_Test")

    ! test version of routine below

    integer(int64), intent(out)    :: nFreeRegions,nParamsPerRegion

!
! Copy to model shared variable:
!
    resistivityFile = 'Demo.2'

!
! Load MARE2DEM .resistivity, .poly, .settings and data files:
!
    call displayBanner
    call readModel
    call readData
    if ( nFree < 1000 ) lUseInversionMeshCoarsening = .false. ! don't use coarsening if small model

    nParamsPerRegion = nRhoPerRegion
    nFreeRegions     = nFree/nRhoPerRegion

    write(*,*) 'nFreeRegions: ',nFreeRegions
    write(*,*) 'nParamsPerRegion: ',nParamsPerRegion

    end subroutine mare2dem_load_files_Test
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
    subroutine mare2dem_load_files(inputResistivityFile,fileLength,nFreeRegions,nParamsPerRegion) bind(C, name="mare2dem_load_files")
!
! Loads the .resistivity, .poly, data and settings files into mare2dem memory
! Returns number of free regions and the number of resistivity params per region (so total free params = nFreeRegions*nParamsPerRegion)
!

    character(1), dimension(256)   :: inputResistivityFile            ! e.g. Demo.0 . note mare2dem wants this without the .resistivity extension
    integer(int64), intent(out)    :: nFreeRegions,nParamsPerRegion
    integer(int64), intent(in)  :: fileLength

!
! Copy to model shared variable:
!
    resistivityFile = inputResistivityFile(1)(1:256) ! clunky way to deal with bind(C) only allowing character(1) string,
    ! so we use array of 256 of them, and this inserts those into a character(256) string
!    resistivityFile = inputResistivityFile(1)(1:fileLength)

!    write(*,*) 'file name length: ',fileLength
    write(*,*) 'resistivityFile: ',resistivityFile !kwk debug
!    write(*,*) 'resistivityFile: ',inputResistivityFile(1)(1:fileLength) !kwk debug

!
! Load MARE2DEM .resistivity, .poly, .settings and data files:
!
    call displayBanner
    call readModel
    call readData
    if ( nFree < 1000 ) lUseInversionMeshCoarsening = .false. ! don't use coarsening if small model

    nParamsPerRegion = nRhoPerRegion
    nFreeRegions     = nFree/nRhoPerRegion

    write(*,*) 'nFreeRegions: ',nFreeRegions
    write(*,*) 'nParamsPerRegion: ',nParamsPerRegion

    end subroutine mare2dem_load_files

!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
    subroutine mare2dem_get_params(nLog10RhoFree,freeCentroids)  bind(C, name="mare2dem_get_params")

    use kx_io

    real(real64), dimension(:),allocatable,intent(out)   :: nLog10RhoFree
    real(real64), dimension(:,:),allocatable,intent(out) :: freeCentroids

    integer :: i,j,ict
!
! Copy free parameters to output:
!
    allocate(nLog10RhoFree(nFree))
    ict = 0
    do i = 1,nRegions
        do j = 1,nRhoPerRegion
            if (iFreeParam((i-1)*nRhoPerRegion+j)>0) then ! note this is NOT setup for cole-cole ip
                ict = ict + 1
                nLog10RhoFree(ict) = log10(rhoParams((i-1)*nRhoPerRegion+j))
            endif
        enddo
    enddo

!
! Get free region centroids:
!
    call getFreeRegionCentroids()
    allocate(freeCentroids(size(freeRegionCentroids,1),2))
    freeCentroids = freeRegionCentroids

    end subroutine mare2dem_get_params
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
    subroutine mare2dem_forward(nproc_per_team,nLog10RhoFreeIn,chiSquaredMisfit) bind(C, name="mare2dem_forward")
!
! Routine for computing MARE2DEM forward response that can be called from external libraries.
! You must first call mare2dem_load_files(inputResistivityFile) on the rank 0 processor
!
! Inputs:   nproc_per_team   variable to use to create local mpi communicator in Fortran in same way as done in
!                            external calling code, which may be running multiple forward models in parallel
!                            on separate teams of processors. e.g.,Julia  PT chains running using nproc_per_team cores
!                            for each PT chain. The calling routine needs to have teams set up identically
!                            (we recreate the MPI communicator here so its available to mare2dem's fortran routines)
!           nLog10RhoFreeIn  1D array (size nFree) of free parameter log10 resistivities (log10 ohm-m)
!
! Outputs:  chiSquaredMisfit data fit to input model using nRhoFreeIn
!

    integer(int64),                 intent(in)   :: nproc_per_team
    real(real64), dimension(nFree), intent(in)   :: nLog10RhoFreeIn
    real(real64),                   intent(out)  :: chiSquaredMisfit

    integer             :: rank,nproc,nteams,team, ierr

!
! Intel MKL specific: make sure one thread per MPI process so mkl solver calls don't oversubscribe
!
    call mkl_set_num_threads ( 1 )

!
! Get rank and setup communicator team(s)
!
    call mpi_comm_rank( MPI_COMM_WORLD, rank, ierr )
    call mpi_comm_size( MPI_COMM_WORLD, nproc, ierr )
    if (nproc == 1) then
        write(*,*) ' '
        write(*,*) '!!!!! Error: MARE2DEM needs to be run using at least 2 processors !!!'
        write(*,*) '             You specified only one, sorry, stopping!'
        write(*,*) ' '
        call exitMARE2DEM()
    endif
    nworkers = nproc - 1

    nteams = ceiling(dble(nworkers)/dble(nproc_per_team))
    team   = (rank-1)/nproc_per_team

    call mpi_comm_split(MPI_COMM_WORLD, team, rank, mcomm ,ierr) ! assigns this rank to be part of mcomm for team.
                                                                 ! so each team has different mcomm

    write(*,*) 'rank, team: ', rank,team

    call mpi_barrier(MPI_COMM_WORLD,ierr)

    !
    ! Launch manager and worker controllers:
    !
    if (rank /= manager) then ! worker

        call mpi_worker_controller

    else ! manager

        !
        ! Initialize the worker status array:
        !
        allocate (lWorker_status(nworkers))
        lWorker_status  = .true.

        !
        ! Compute forward response:
        !
        linversion = .false. ! tell mare2dem to NOT compute Jacobians otherwise it'll run way more slowly

        call computeFwd( linversion, nLog10RhoFreeIn )

        !
        ! Compute chi-squared misfit
        !
        chiSquaredMisfit =  sum( ((d - dm)/sd)**2 )

        ! debug check:
        write(*,*) 'rms: ',sqrt(chiSquaredMisfit/size(d))

        !
        ! Tell workers to exit mpi_worker_controller:
        !
        call mpi_shutDownWorkers()

    endif

    !
    ! Broadcast chiSquaredMisfit to all processors on local communicator
    !
    call mpi_bcast(chiSquaredMisfit, 1, MPI_DOUBLE_PRECISION, manager, mcomm, ierr)

    !debug: write(*,*) 'rank, chiSquaredMisfit: ',rank, chiSquaredMisfit

    end subroutine mare2dem_forward


!==================================================================================================================================!
!===================================================================================================================== displayBanner
!==================================================================================================================================!
    subroutine displayBanner

    implicit none


    write(*,*) ' '
    write(*,*) '============================= MARE2DEM ==================================='
    write(*,*) ' '
    write(*,*) ' MARE2DEM: Modeling with Adaptively Refined Elements for 2.5D EM'
    write(*,*) ''
    write(*,*) ' Version: 4.0 September 9, 2018'
    write(*,*) ' '
    write(*,*) ' A parallel goal-oriented adaptive finite element forward and inverse'
    write(*,*) ' modeling code for electromagnetic fields from electric dipoles, magnetic'
    write(*,*) ' dipoles and magnetotelluric sources in triaxially anisotropic conducting'
    write(*,*) ' media. Iterative adaptive mesh refinement is accomplished using the'
    write(*,*) ' goal-oriented error estimation method described in Key and Ovall (2011) '
    write(*,*) ' Inversion is accomplished with Occam''s method (Constable et al., 1987).'
    write(*,*) ' Key (2016) describes most of the features in the current 2018 version '
    write(*,*) ' of the code.'
    write(*,*) ' '
    write(*,*) ' When citing the code, please use the most recent reference:'
    write(*,*) ' '
    write(*,*) ' Key, K. MARE2DEM: a 2-D inversion code for controlled-source electromagnetic '
    write(*,*) '     and magnetotelluric data. Geophysical Journal International 207, '
    write(*,*) '     571–588 (2016).  '
    write(*,*) ''
    write(*,*) ' This work is currently supported by: '
    write(*,*) ''
    write(*,*) ' Electromagnetic Methods Research Consortium'
    write(*,*) ' Lamont-Doherty Earth Observatory'
    write(*,*) ' Columbia University'
    write(*,*) ' http://emrc.ldeo.columbia.edu'
    write(*,*) ' '
    write(*,*) ' Originally funded by:'
    write(*,*) ''
    write(*,*) ' Seafloor Electromagnetic Methods Consortium '
    write(*,*) ' Scripps Institution of Oceanography '
    write(*,*) ' University of California San Diego'
    write(*,*) ''
    write(*,*) ' Copyright (C) 2017-2018'
    write(*,*) ' Kerry Key'
    write(*,*) ' Lamont-Doherty Earth Observatory'
    write(*,*) ' Columbia University'
    write(*,*) ' http://emlab.ldeo.columbia.edu'
    write(*,*) ' '
    write(*,*) ' Copyright (C) 2008-2016'
    write(*,*) ' Kerry Key'
    write(*,*) ' Scripps Institution of Oceanography'
    write(*,*) ' University of California, San Diego'
    write(*,*) ''
    write(*,*) ''
    write(*,*) ' This file is part of MARE2DEM.'
    write(*,*) ''
    write(*,*) ' MARE2DEM is free software: you can redistribute it and/or modify'
    write(*,*) ' it under the terms of the GNU General Public License as published by'
    write(*,*) ' the Free Software Foundation, either version 3 of the License, or'
    write(*,*) ' (at your option) any later version.'
    write(*,*) ''
    write(*,*) ' MARE2DEM is distributed in the hope that it will be useful,'
    write(*,*) ' but WITHOUT ANY WARRANTY; without even the implied warranty of'
    write(*,*) ' MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the'
    write(*,*) ' GNU General Public License for more details.'
    write(*,*) ''
    write(*,*) ' You should have received a copy of the GNU General Public License'
    write(*,*) ' along with MARE2DEM.  If not, see <http://www.gnu.org/licenses/>. '
    write(*,*) ' '
    write(*,*) '=========================================================================='

    end subroutine displayBanner

end module mare2dem_lib_interface

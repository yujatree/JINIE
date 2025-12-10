module READER

   ! Variables Declaration----------------------------------------------------------------------	
	
   implicit none
        
   ! For loop
   integer :: i, j, k, l, nfile
   
   ! File I/O system
   real                                           :: dummyr
   integer, dimension(20)                         :: dummyi(20)
   character(4)                                   :: dummyc
   character(80), allocatable, dimension(:)       :: filename
   integer                                        :: num_of_args

   ! Checking Trajectory variables
   integer 	    			          :: init_check, dt_check, end_check, num_check

   ! Data variables

   real, allocatable, dimension(:)                :: xr, yr, zr
   double precision, allocatable, dimension(:,:)  :: x, y, z
   double precision, allocatable, dimension(:,:)  :: box
   double precision, allocatable, dimension(:)    :: volume
   
   integer                                        :: nmedia, total_traj, atom_idx, traj_idx
   integer, allocatable, dimension(:)             :: ntraj
					           ! num_of_args : number of trajectories to link
   				                   ! nmedia : number of particles
   				                   ! total_traj : number of trajectory

   ! For LPS system
   integer				          :: num_Li, num_P, num_S	

   contains

! ----------------------------------------------------------------------------------------------
	
   subroutine dcd_reader

   !------------------------------------------------!
   !                                                !
   !                Main Link Program               !
   !              made by Park HyungShick           !
   !                modified by EunHyeok            !
   !                modified by Jina    	    !
   !				                    !
   !               Ver.1 2020.06.26                 !
   !               Ver.2 2022.01.03 		    !
   !   		   Ver.3 2025.04.27/12.08           !
   !				                    !
   !------------------------------------------------!

   ! Refer to JENOGMIN KIM READ_DCD program
   !          Mar. 18, 2013

   
   ! Argument-----------------------------------------------------------------------------------	

        num_of_args=iargc() ! # of files

        if (num_of_args .eq. 0) then
           print*, "------------ERROR------------"
           print*, "USAGE : ./run.x dcd1 dcd2 ..."
           print*, "-----------------------------"
           stop
        endif

	allocate(filename(num_of_args), ntraj(num_of_args))

	print*
        print*, "============================="
	print*, " DCD trajectory file READING "
        print*, "============================="
	print*
	print*, num_of_args, " files are loaded"
	print*
	do i=1,num_of_args
           call getarg(i, filename(i)) ! ith arg is saved at filename(i)
	enddo
	
   ! Read DCD file------------------------------------------------------------------------------	
    
	! Checking ensembles to have same info--------------------------------------------------

	do nfile=1,num_of_args
	
           open(10+nfile,file=trim(filename(nfile)),status='old',form='unformatted')

           read(10+nfile) dummyc, ntraj(nfile), (dummyi(i),i=1,8), dummyr, (dummyi(i),i=9,17)
           read(10+nfile) dummyi(18), dummyr
           read(10+nfile) nmedia
           print*, "num of frames: ", ntraj(nfile)
	   print*, "num of atoms : ", nmedia
           ! dummyi(1); trajectory starting point
           ! dummyi(2); dt
           ! dummyi(3); trajectory end point

	   if (nfile .eq. 1) then
	      init_check = dummyi(1)
	      dt_check   = dummyi(2)
	      end_check  = dummyi(3)
	      num_check  = nmedia
	      print*, 'Trajectory end : ', end_check
	   else
	      ! Trajectory contininueous check
	      if (end_check .ne. dummyi(1)) then
                 print*, "-------------------ERROR-------------------"
                 print*, "Trajectories should be linked continueoulsy"
                 print*, "-------------------------------------------"
                 stop
              endif

	      ! Trajectory interval check
	      if (dt_check .ne. dummyi(2)) then
                 print*, "---------------ERROR---------------"
                 print*, "Trajectory interval should be equal"
                 print*, "-----------------------------------"
                 stop
              endif

	      ! Nmedia check
	      if (num_check .ne. nmedia) then
                 print*, "-------------ERROR-------------"
                 print*, "Number of atoms should be equal"
                 print*, "-------------------------------"
                 stop
              endif

	      end_check = dummyi(3)
	      print*, 'Trajectory end : ', end_check
	   endif
	enddo

	print*
	!---------------------------------------------------------------------------------------

	total_traj=int((end_check-init_check)/dt_check)+1
	
	print*, "------------INFO-------------"
	print*, "Initial time    |", init_check
	print*, "End time        |", end_check
	print*, "Time interval   |", dt_check
	print*, "Number of trajs |", total_traj
	print*, "Number of atoms |", nmedia
        print*, "-----------------------------"
	print*
	
	! X,Y,Z(Traj,atom), Box(traj,size)   
	allocate(xr(nmedia), yr(nmedia), zr(nmedia))        
        allocate(x(total_traj,nmedia), y(total_traj,nmedia), z(total_traj,nmedia),&
		 box(total_traj,6), volume(total_traj))

	if (mod(nmedia,8)==0) then
	   num_Li=nmedia*3/8
	   num_P=nmedia/8
	   num_S=nmedia/2
	else
	   print*, "There might be some defects ..."
	   print*
	   stop
	endif

	print*
	print*, "---------LPS System----------"
	print*, " Number of Li+ | ", num_Li
	print*, " Number of P5+ | ", num_P
	print*, " Number of S2- | ", num_S
        print*, "-----------------------------"

	! Reading Trajectory Data---------------------------------------------------------------

	traj_idx=0

	do nfile=1, num_of_args
	   
	   print*, nfile, "th file READING.."
	
	   do i=1, ntraj(nfile)
	      read(10+nfile) (box(i+traj_idx,j), j=1,6) ! xx xy yy yz zx zz
	      read(10+nfile) (xr(j), j=1,nmedia)
	      read(10+nfile) (yr(j), j=1,nmedia)
	      read(10+nfile) (zr(j), j=1,nmedia)
	      ! Improving precision : single > double
              do j = 1, nmedia
                  x(i+traj_idx,j) = real(xr(j), kind=8)
                  y(i+traj_idx,j) = real(yr(j), kind=8)
                  z(i+traj_idx,j) = real(zr(j), kind=8)
              end do
	      volume(i+traj_idx)=box(i+traj_idx,1)*box(i+traj_idx,3)*box(i+traj_idx,6)
	   enddo

	   print*, nfile, "th file DONE.."
	   traj_idx = traj_idx + ntraj(nfile) - 1

	enddo

	print*
        print*, "-----------------------------"
	print*, "          SUCCESS !"
        print*, "-----------------------------"
	print*

	do nfile=1, num_of_args
	   close(10+nfile)
	enddo

	!---------------------------------------------------------------------------------------

   end subroutine dcd_reader

! ----------------------------------------------------------------------------------------------
	
   subroutine lammpstraj_reader

   end subroutine lammpstraj_reader

! ----------------------------------------------------------------------------------------------

end module READER

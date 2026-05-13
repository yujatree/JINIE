program conditional_distance_distribution_function 

! Variables Declaration-----------------------------------------------------------------
        
        use READER

        implicit none

        logical					      :: file_exists

        integer 			              :: bin, bin_num, hop_traj, &
							 count_h, count_coord, idx1,&
							 nsurr
        double precision	                      :: dist, bin_size, &
							 delta_t, h_c, r
        double precision, parameter 	              :: pi=datan(1.d0)*4.d0, &
							 coord_rcut=8.d0
        double precision, dimension(100)	      :: dist_eq, dist_h
        double precision, dimension(:,:), allocatable :: hop, D_eq, D_h, &
							 D_eq_tot, D_h_tot

        character(10)				      :: filenum

! Setting-------------------------------------------------------------------------------

        call dcd_reader

999	format(F7.3," ",*(G0,1X),/)
998	format(A,I5,A)

	open(50, file='/mnt/CALC_SSD/jina/CODES/JINIE/inp/bin.inp', status='old')
	open(51, file='/mnt/CALC_SSD/jina/CODES/JINIE/inp/hop.inp', status='old')

	do nfile=1,num_of_args
	   write(filenum,'(I1)') nfile 

	   inquire(file=trim(adjustl(filenum))//"/hop.out", exist=file_exists)

	   if (file_exists) then
	      open(52+nfile, file=trim(adjustl(filenum))//'/hop.out', status='old')
	   else
	      open(52+nfile, file=trim(adjustl(filenum))//'/OUTPUT/hop.out', status='old')
	   endif
	enddo 
	open(100, file='CDDF_Li.out', status='unknown')

	write(100,'(A/)',advance='no') &
           "# r(A) 1st_eq 1st_h* 2nd_eq 2nd_h* 3rd_eq 3rd_h* 4th_eq 4th_h* 5th_eq 5th_h* ... "
	
	print*
        print*, "-----------------------------"
        print*, "        CDDF START !"
        print*, "-----------------------------"
	print*
	read(50,*) bin_size
	read(50,*) bin
	print*, " Current size/num of bins :"
	write(*,'(22X,F7.3)'), bin_size
	write(*,'(22X,I7)'), bin
	print*
	print*, " CAUTION : valid only w/ NVT"
	read(51,*) delta_t
	read(51,*) h_c
	read(51,*) nsurr
	print*, " Delta t : "
	write(*,'(27X,I2)'), int(delta_t)
	print*, " h_c : "
	write(*,'(27X,F3.1)'), h_c
	print*, " nsurr : "
	write(*,'(27X,I2)'), nsurr

	hop_traj=total_traj-4*int(delta_t/2.d0)

	print*, " Hop traj :", hop_traj

        allocate(D_eq(nsurr,bin), D_h(nsurr,bin)) ! 1st, 2nd, 3rd, 4th, 5th coordinated Li
        allocate(D_eq_tot(nsurr,bin), D_h_tot(nsurr,bin)) ! For all ensembles
	allocate(hop(hop_traj,num_Li))

!---------------------------------------------------------------------------------------

	D_eq_tot=0.d0 ; D_h_tot=0.d0

	do nfile=1, num_of_args

      	   D_eq=0.d0 ; D_h=0.d0

	   count_h=0

	   read(52+nfile,*) 

	   do i=1, hop_traj

	      if (mod(i,10)==0) write(*,998) " Processing ", i, " th trajectory"

	      read(52+nfile,*) dummyi(1), (hop(i,j), j=1,num_Li)

	      do j=1, num_Li!==============================================

                 count_coord=0

	         if (hop(i,j)>h_c) then ! In case of hopping--------------

	   	 count_h=count_h+1

	            do k=1, num_Li
	               if (j==k) cycle 

	               call pbc_dist(nfile, i, j, k, dist)

	               if (dist<coord_rcut) then
	                  count_coord=count_coord+1
	                  dist_eq(count_coord)=dist
	                  dist_h(count_coord)=dist
	               endif

	            enddo
	   
	            call sort_pick(dist_eq(1:count_coord))
	            call sort_pick(dist_h(1:count_coord))

	            do l=1,nsurr  ! only coordinated Li
	               bin_num=int(dist_eq(l)/bin_size)+1
	               D_eq(l,bin_num)=D_eq(l,bin_num)+1.d0
	   	    D_h(l,bin_num)=D_h(l,bin_num)+1.d0
	            enddo

	         else ! Everage distribution---------------------------------

	            do k=1, num_Li
	   	    if (j==k) cycle

	               call pbc_dist(nfile, i, j, k, dist)

	               if (dist<coord_rcut) then
	                  count_coord=count_coord+1
	                  dist_eq(count_coord)=dist
	               endif
	            enddo
	   
	            call sort_pick(dist_eq(1:count_coord))

	            do l=1,nsurr  ! only coordinated Li
	               bin_num=int(dist_eq(l)/bin_size)+1
	               D_eq(l,bin_num)=D_eq(l,bin_num)+1.d0
	            enddo

	         endif ! ----------------------------------------------------

	      enddo!==========================================================

	   enddo

	   D_eq_tot=D_eq_tot+D_eq/dble(num_Li)/dble(hop_traj)
           D_h_tot=D_h_tot+D_h/dble(count_h)

        enddo	

	D_eq_tot=D_eq_tot/dble(num_of_args)
	D_h_tot=D_h_tot/dble(num_of_args)

	do i=1, bin
	   r=(dble(i)-0.5d0)*bin_size
	   write(100,999) r, (D_eq_tot(j,i), D_h_tot(j,i), j=1,nsurr)
	enddo

	print*
	print*, "   Saved as 'CDDF_Li.out'    "
        print*, "-----------------------------"
	print*, "          SUCCESS !"
        print*, "-----------------------------"
	print*

	close(100)

!---------------------------------------------------------------------------------------

contains
  
!--------------------------------------------------------------------------------------------

   subroutine pbc_dist(n, time, atom1, atom2, dist)
        
        use READER

        implicit none

        integer          :: n, time, atom1, atom2
        double precision :: dx, dy, dz, dist	

        dx=ens(n)%x(time,atom1)-ens(n)%x(time,atom2)
        dy=ens(n)%y(time,atom1)-ens(n)%y(time,atom2)
        dz=ens(n)%z(time,atom1)-ens(n)%z(time,atom2)
	dx=dx-nint(dx/ens(n)%box(time,1))*ens(n)%box(time,1)
	dy=dy-nint(dy/ens(n)%box(time,3))*ens(n)%box(time,3)
	dz=dz-nint(dz/ens(n)%box(time,6))*ens(n)%box(time,6)
	dist=dsqrt(dx*dx+dy*dy+dz*dz)

   end subroutine pbc_dist

!--------------------------------------------------------------------------------------------

   subroutine sort_pick(arr) 

        implicit none

        integer :: i, j, n
        double precision :: a
        double precision, dimension(:) :: arr

        n = size(arr) 

        do i=2,n 
           a=arr(i)
           do j=i-1,1,-1 
              if (arr(j)<=a) exit
              arr(j+1)=arr(j)
           enddo
           arr(j+1)=a 
        enddo

   end subroutine sort_pick

!--------------------------------------------------------------------------------------------

end program conditional_distance_distribution_function

program radial_distribution_function 

! Variables Declaration-----------------------------------------------------------------
        
        use READER

        implicit none
 
        integer 			              :: bin, bin_num, idx1, idx2, &
						         time_window
        double precision	                      :: dist, bin_size, r
        double precision, parameter 	              :: pi=datan(1.d0)*4.d0
        double precision, dimension(:,:), allocatable :: g

! Setting-------------------------------------------------------------------------------

        call dcd_reader

999	format(F7.3," ",7(G0,1X)/)
998	format(A,I4,A)

	open(50, file='/mnt/CALC_SSD/jina/CODES/JINIE/inp/bin.inp', status='old')
	open(100, file='RDF.out', status='unknown')
	write(100,*) "# r(A) Li_Li P_P S_S Li_P Li_S P_S Total"	

	print*
        print*, "-----------------------------"
        print*, "         RDF START !"
        print*, "-----------------------------"
	print*
	read(50,*) bin_size
	read(50,*) bin
	read(50,*) time_window
	print*, " Current size/num of bins :"
	write(*,'(22X,F7.3)'), bin_size
	write(*,'(22X,I7)'), bin
	print*, " Trajectories to average  : "
	write(*,'(20X,I7)'), time_window
time_window=1
	print*
	print*, " CAUTION : valid only w/ NVT"

        allocate(g(7,bin))

!---------------------------------------------------------------------------------------

	g=0.d0

	do i=total_traj-time_window+1, total_traj

	   if (mod(i,10)==0) write(*,998) " Processing ", i, " th trajectory"

	   ! Li-Li RDF
	   do j=1, num_Li-1
	      do k=j+1, num_Li
	         call pbc_dist(i, j, k, dist)
	         bin_num=int(dist/bin_size)+1
	         g(1,bin_num)=g(1,bin_num)+2.d0
	      enddo
	   enddo
	   ! P-P RDF
	   do j=1, num_P-1
	      idx1=num_Li+j    ! To match P idx
	      do k=j+1, num_P
	         idx2=num_Li+k
	         call pbc_dist(i, idx1, idx2, dist)
	         bin_num=int(dist/bin_size)+1
	         g(2,bin_num)=g(2,bin_num)+2.d0
	      enddo
	   enddo
	   ! S-S RDF
	   do j=1, num_S-1
	      idx1=num_Li+num_P+j    ! To match P idx
	      do k=j+1, num_S
	         idx2=num_Li+num_P+k
	         call pbc_dist(i, idx1, idx2, dist)
	         bin_num=int(dist/bin_size)+1
	         g(3,bin_num)=g(3,bin_num)+2.d0
	      enddo
	   enddo
	   ! Li-P RDF
	   do j=1, num_Li
	      do k=1, num_P
	         idx1=num_Li+k ! To match P idx
	         call pbc_dist(i, j, idx1, dist)
	         bin_num=int(dist/bin_size)+1
	         g(4,bin_num)=g(4,bin_num)+1.d0
	      enddo
	   enddo
	   ! Li-S RDF
	   do j=1, num_Li
	      do k=1, num_S
	         idx1=num_Li+num_P+k ! To match S idx
	         call pbc_dist(i, j, idx1, dist)
	         bin_num=int(dist/bin_size)+1
	         g(5,bin_num)=g(5,bin_num)+1.d0
	      enddo
	   enddo
	   ! P-S RDF
	   do j=1, num_P
	      idx1=num_Li+j
	      do k=1, num_S
	         idx2=num_Li+num_P+k ! To match P idx
	         call pbc_dist(i, idx1, idx2, dist)
	         bin_num=int(dist/bin_size)+1
	         g(6,bin_num)=g(6,bin_num)+1.d0
	      enddo
	   enddo
	
	enddo

	g=g*volume(total_traj)/dble(time_window)

	do i=1, bin

	   g(7,i)=(g(1,i)+g(2,i)+g(3,i)+g(4,i)+g(5,i)+g(6,i))/dble(nmedia*nmedia)
	   g(1,i)=g(1,i)/dble(num_Li*(num_Li-1))
	   g(2,i)=g(2,i)/dble(num_P*(num_P-1))
	   g(3,i)=g(3,i)/dble(num_S*(num_S-1))
	   g(4,i)=g(4,i)/dble(num_Li*num_P)
	   g(5,i)=g(5,i)/dble(num_Li*num_S)
	   g(6,i)=g(6,i)/dble(num_P*num_S)

	   r=(dble(i)-0.5)*bin_size

	   do j=1,7
	      if (i==1) then
	         g(j,i)=0.d0
	      else
	         g(j,i)=g(j,i)/(4.d0*pi*r*r*bin_size)
	      endif
	   enddo

	   write(100,999,advance='no') r, (g(j,i), j=1,7)

	enddo

	print*
	print*, "   Saved as 'RDF.out' file   "
        print*, "-----------------------------"
	print*, "          SUCCESS !"
        print*, "-----------------------------"
	print*

	close(100)

!---------------------------------------------------------------------------------------

contains
  
!--------------------------------------------------------------------------------------------

   subroutine pbc_dist(time, atom1, atom2, dist)
        
        use READER

        implicit none

        integer          :: time, atom1, atom2
        double precision :: dx, dy, dz, dist	

        dx=x(time,atom1)-x(time,atom2)
        dy=y(time,atom1)-y(time,atom2)
        dz=z(time,atom1)-z(time,atom2)
	dx=dx-nint(dx/box(time,1))*box(time,1)
	dy=dy-nint(dy/box(time,3))*box(time,3)
	dz=dz-nint(dz/box(time,6))*box(time,6)
	dist=dsqrt(dx*dx+dy*dy+dz*dz)

   end subroutine pbc_dist

!--------------------------------------------------------------------------------------------

end program radial_distribution_function

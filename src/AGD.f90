program angle_distribution

! Variables Declaration-----------------------------------------------------------------
        
        use READER

        implicit none
 
        integer				    	      :: bin, bin_num, time_window, &
							 count, pidx, sidx, total_count
        double precision	           	      :: bin_size, dist, &
							 tmp_dx, tmp_dy, tmp_dz, &
							 ori_dx, ori_dy, ori_dz, theta,&
							 mean, total_mean
	double precision, parameter	   	      :: pi=datan(1.d0)*4.d0, &
							 rad2deg=57.2958d0 , &
							 P_S=2.5d0
	double precision, dimension(:), allocatable   :: angle ! histogram
	double precision, dimension(:,:), allocatable :: dx, dy, dz 

! Setting-------------------------------------------------------------------------------

        call dcd_reader

999	format(F7.3," ",G0,1X/)
998	format(A,I5,A)

	open(50, file='/mnt/CALC_SSD/jina/CODES/JINIE/inp/bin.inp', status='old')
	open(100, file='AGD.out', status='unknown')
	write(100,*) "# angle(deg) probability"	

	print*
        print*, "-----------------------------"
        print*, "         AGD START !"
        print*, "-----------------------------"
	print*

	read(50,*) bin_size
	read(50,*)
	read(50,*) time_window
	bin_size=bin_size/5.d0
	bin=int(pi/bin_size)+1
	print*, " Current size/num of bins :"
	write(*,'(22X,F7.3)') bin_size
	write(*,'(22X,I7)'), bin
	print*, " Trajectories to average  : "
	write(*,'(20X,I7)'), time_window

	allocate(angle(bin))
	allocate(dx(num_P,num_S),dy(num_P,num_S),dz(num_P,num_S))

!---------------------------------------------------------------------------------------

        angle=0.d0
	total_mean=0.d0 ; total_count=0

	do i=total_traj-time_window+1, total_traj

	   if (mod(i,10)==0) write(*,998) " Processing ", i, " th trajectory"

	   dx=0.d0 ; dy=0.d0 ; dz=0.d0

	   ! initialize (collect all the P-S vectors in polyanion)
	   do j=1, num_P
	      pidx=num_Li+j
	      do k=1, num_S
		 sidx=num_Li+num_P+k
	         call pbc_dist(i, pidx, sidx, tmp_dx, tmp_dy, tmp_dz, dist)
	         if (dist < P_S) then
		    dx(j,k)=tmp_dx
	            dy(j,k)=tmp_dy
		    dz(j,k)=tmp_dz
		 endif
	      enddo
	   enddo	
	     
	   do j=1, num_P
	      ori_dx=0.d0
	      mean=0.d0
	      count=0
	      do k=1, num_S
	         if (dx(j,k)==0.d0) goto 1000 ! No bonding

		 if (ori_dx==0.d0) then       ! Setting an Original Vector
		    ori_dx=dx(j,k) ; ori_dy=dy(j,k) ; ori_dz=dz(j,k)
		 else			   ! Calculating angles w.r.t Original vector
		    call get_angle(ori_dx, ori_dy, ori_dz, dx(j,k), dy(j,k), dz(j,k), theta)
		    mean=mean+theta
		    count=count+1
		 endif
1000             continue
	      enddo
	      if (count==0) cycle
	      bin_num=int(mean/dble(count)/bin_size)+1
	      total_mean=total_mean+mean/dble(count)
	      total_count=total_count+1
	      if (bin_num < 1 .or. bin_num > bin) then
	         write(*,*) " ERROR: bin_num out of bounds:", bin_num, " (max=", bin, ")"
	      else
	         angle(bin_num)=angle(bin_num)+1.d0
	      endif
	   enddo
	enddo      

	angle=angle/dble(num_P)/dble(time_window)

	do i=1, bin
	   theta=(dble(i)-0.5d0)*bin_size*rad2deg
	   write(100,999,advance='no') theta, angle(i)
	enddo
	if (total_count > 0) then
	   write(100,'(A,1X,F7.3)') "# AVERAGE OF THETA : ", total_mean/dble(total_count)*rad2deg
	else
	   write(100,'(A)') "# AVERAGE OF THETA : N/A"
	endif

	print*
	print*, "   Saved as 'AGD.out' file   "
        print*, "-----------------------------"
	print*, "          SUCCESS !"
        print*, "-----------------------------"
	print*

	close(50)
	close(100)
contains
  
!--------------------------------------------------------------------------------------------

   subroutine pbc_dist(time, atom1, atom2, dx, dy, dz, dist)
        
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

   subroutine get_angle(r0_dx, r0_dy, r0_dz, r1_dx, r1_dy, r1_dz, theta)
        
        use READER

        implicit none

        double precision :: r0_dx, r0_dy, r0_dz, r1_dx, r1_dy, r1_dz, &
			    inner, r0, r1, theta

	inner=r0_dx*r1_dx+r0_dy*r1_dy+r0_dz*r1_dz
	r0=dsqrt(r0_dx*r0_dx+r0_dy*r0_dy+r0_dz*r0_dz)
	r1=dsqrt(r1_dx*r1_dx+r1_dy*r1_dy+r1_dz*r1_dz)
	theta=dacos(inner/(r0*r1))

   end subroutine get_angle

!--------------------------------------------------------------------------------------------

end program angle_distribution

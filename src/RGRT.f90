program self_part_van_hove_correlation_function_of_rotation

! Variables Declaration-----------------------------------------------------------------
        
        use READER

        implicit none
 
        integer 			              :: bin, bin_num, pidx, sidx
        double precision	                      :: dist, bin_size, &
							 tmp_dx, tmp_dy, tmp_dz, &
							 ori_dx, ori_dy, ori_dz, theta
	double precision, parameter	   	      :: pi=datan(1.d0)*4.d0, &
							 rad2deg=57.2958d0 , &
							 P_S=2.5d0

        integer, dimension(:), allocatable 	      :: counter, pcounter
	type					      :: ptr
           double precision, pointer, dimension(:,:)  :: dx, dy, dz
	endtype
	type(ptr), allocatable, dimension(:)	      :: t

        double precision, dimension(:,:), allocatable :: angle

        character(len=100) :: num_str, fmt_str, head_num, head_str

! Setting-------------------------------------------------------------------------------

        call dcd_reader

998	format(A,I5,A)

	open(50, file='/mnt/CALC_SSD/jina/CODES/JINIE/inp/bin.inp', status='old')
	open(100, file='RGRT.out', status='unknown')

        write(num_str, '(I0)') int(total_traj/3)
	write(head_num, '(I0)') int(total_traj/3)
	fmt_str = '(F7.3," ",'//trim(num_str)//'(G0,1X))'
	head_str = '('//trim(head_num)//'(I6))'
	write(100,fmt=head_str,advance='no') (i-1, i=1,int(total_traj/3))
	write(100,*)

	print*
        print*, "-----------------------------"
        print*, "    Rotatinal GSRT START !"
        print*, "-----------------------------"
	print*

	read(50,*) bin_size
	read(50,*) bin

	bin_size=bin_size/5.d0
	bin=int(pi/bin_size)+1

	print*, " Current size/num of bins :"
	write(*,'(22X,F7.3)'), bin_size
	write(*,'(22X,I7)'), bin
	print*

	allocate(t(total_traj))
	do i=1, total_traj
	   allocate(t(i)%dx(num_P,4),&
		    t(i)%dy(num_P,4),&
		    t(i)%dz(num_P,4))
	   t(i)%dx=0.d0 ; t(i)%dy=0.d0 ; t(i)%dz=0.d0
	enddo

	allocate(counter(total_traj))
	allocate(pcounter(total_traj))
	allocate(angle(total_traj,bin))

!---------------------------------------------------------------------------------------

	! Initialization : dx(t,nth bond) vector
	do i=1, total_traj
	   if (mod(i,10)==0) write(*,998) " Initializing ", i, " th traj"
	   pcounter=0
	   do j=1, num_P
	      pidx=num_Li+j
	      do k=1, num_S
	         sidx=num_Li+num_P+k
	         call pbc_dist(i, pidx, sidx, tmp_dx, tmp_dy, tmp_dz, dist)
	         if (dist < P_S) then
	  	    pcounter(j)=pcounter(j)+1
		    if (pcounter(j)>4) then
		       print*, " There're more than four S atoms around one P atom...! "
		       print*, "Exit..."
		       stop
		    endif
		    t(i)%dx(j,pcounter(j))=tmp_dx
		    t(i)%dy(j,pcounter(j))=tmp_dy
		    t(i)%dz(j,pcounter(j))=tmp_dz
	         endif
	      enddo
	   enddo
	enddo

	print*, " Initialization done!"

	! ROTATIONAL GSRT
	
	counter=0
	angle=0.d0

	do i=1, int(total_traj/3) ! lag t, delta t

	   if (mod(i,10)==0) write(*,998) " Processing ", i, " th trajectory"

	   do j=1, total_traj-i   ! origin t

	      do k=1, num_P
	         pidx=num_Li+k
	         do l=1, 4
		    ! in case of no bonding
		    if ((t(j)%dx(k,l)==0.d0).or.(t(j+i)%dx(k,l)==0.d0)) cycle
 
	            call get_angle(t(j)%dx(k,l),t(j)%dy(k,l),t(j)%dz(k,l),&
				   t(j+i)%dx(k,l),t(j+i)%dy(k,l),t(j+i)%dz(k,l),&
				   theta)

	            bin_num=int(theta/bin_size)+1
	            angle(i,bin_num)=angle(i,bin_num)+1.d0
		    counter(i)=counter(i)+1
		 enddo
	      enddo

	   enddo

	   do j=1, bin
	      if (counter(i)>0) then
	         angle(i,j)=angle(i,j)/dble(counter(i))/bin_size
	      else
	         angle(i,j)=0.d0
	      endif
	   enddo

	enddo

	do i=1, bin
	   theta=(dble(i)-0.d0)*bin_size*rad2deg
   	   write(100,fmt=fmt_str,advance='no') theta, (angle(j,i), j=1,int(total_traj/3))	
  	   write(100,*)	   
	enddo
	
	print*
	print*, "   Saved as 'RGRT.out' file   "
        print*, "-----------------------------"
	print*, "          SUCCESS !"
        print*, "-----------------------------"
	print*

	close(100)

!---------------------------------------------------------------------------------------

contains
  
!--------------------------------------------------------------------------------------------

   subroutine pbc_dist(t, atom1, atom2, dx, dy, dz, dist)
        
        use READER

        implicit none

        integer          :: t, atom1, atom2
        double precision :: dx, dy, dz, dist	

        dx=x(t,atom1)-x(t,atom2)
        dy=y(t,atom1)-y(t,atom2)
        dz=z(t,atom1)-z(t,atom2)
	dx=dx-nint(dx/box(t,1))*box(t,1)
	dy=dy-nint(dy/box(t,3))*box(t,3)
	dz=dz-nint(dz/box(t,6))*box(t,6)
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

end program self_part_van_hove_correlation_function_of_rotation

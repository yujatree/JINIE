program self_part_van_hove_correlation_function

! Variables Declaration-----------------------------------------------------------------
        
        use READER

        implicit none
 
        integer 			              :: bin, bin_num
        double precision	                      :: dist, bin_size
        double precision, parameter 	              :: pi=datan(1.d0)*4.d0
        double precision, dimension(:), allocatable   :: r, norm_factor
        double precision, dimension(:,:), allocatable :: g

	type :: ptr
	   double precision, pointer, dimension(:,:)  :: p
	end type

	type(ptr), dimension(3)          :: gsrt


        character(len=100) :: num_str, fmt_str, head_num, head_str

! Setting-------------------------------------------------------------------------------

        call dcd_reader

999	format(F7.3," ",7(G0,1X)/)
998	format(A,I4,A)

	open(50, file='/mnt/CALC_SSD/jina/CODES/JINIE/inp/bin.inp', status='old')
	open(100, file='GSRT_Li.out', status='unknown')
	open(101, file='GSRT_P.out', status='unknown')
	open(102, file='GSRT_S.out', status='unknown')

        write(num_str, '(I0)') total_traj
	write(head_num, '(I0)') total_traj
	fmt_str = '(F7.3," ",'//trim(num_str)//'(G0,1X))'
	head_str = '('//trim(head_num)//'(I6))'
	write(100,fmt=head_str,advance='no') (i-1, i=1,total_traj)
	write(100,*)
	write(101,fmt=head_str,advance='no') (i-1, i=1,total_traj)
	write(101,*)
	write(102,fmt=head_str,advance='no') (i-1, i=1,total_traj)
	write(102,*)

	print*
        print*, "-----------------------------"
        print*, "         GSRT START !"
        print*, "-----------------------------"
	print*

	read(50,*) bin_size
	read(50,*) bin
	print*, " Current size/num of bins :"
	write(*,'(22X,F7.3)'), bin_size
	write(*,'(22X,I7)'), bin
	print*

	allocate(r(bin),norm_factor(bin))
	allocate(g(3,bin))
	do i=1,3
	   allocate(gsrt(i)%p(total_traj,bin))
	enddo

!---------------------------------------------------------------------------------------


	do i=1, bin
	   r(i)=(dble(i)-0.5d0)*bin_size
	   norm_factor(i)=4.d0*pi*r(i)*r(i)*bin_size
	enddo

	do i=1, total_traj

	   if (mod(i,10)==0) write(*,998) " Processing ", i, " th trajectory"

	   g=0.d0

	   do j=1, num_Li
	      call get_dist(i, 1, j, dist)
	      bin_num=int(dist/bin_size)+1
	      g(1,bin_num)=g(1,bin_num)+1.d0
	   enddo
	   do j=1, num_P
	      call get_dist(i, 1, num_Li+j, dist)
	      bin_num=int(dist/bin_size)+1
	      g(2,bin_num)=g(2,bin_num)+1.d0
	   enddo
	   do j=1, num_S
	      call get_dist(i, 1, num_Li+num_P+j, dist)
	      bin_num=int(dist/bin_size)+1
	      g(3,bin_num)=g(3,bin_num)+1.d0
	   enddo
	
	   do j=1, bin
	      gsrt(1)%p(i,j)=g(1,j)/dble(num_Li)/norm_factor(j)
	      gsrt(2)%p(i,j)=g(2,j)/dble(num_P)/norm_factor(j)
	      gsrt(3)%p(i,j)=g(3,j)/dble(num_S)/norm_factor(j)
	   enddo

	enddo

	do i=1, bin
 	   write(100,fmt=fmt_str,advance='no') r(i), (gsrt(1)%p(j,i), j=1,total_traj)	
	   write(100,*)	   
 	   write(101,fmt=fmt_str,advance='no') r(i), (gsrt(2)%p(j,i), j=1,total_traj)		   
	   write(101,*)	   
 	   write(102,fmt=fmt_str,advance='no') r(i), (gsrt(3)%p(j,i), j=1,total_traj)		   
	   write(102,*)	   
	enddo

	print*
	print*, "   Saved as 'GSRT.out' file   "
        print*, "-----------------------------"
	print*, "          SUCCESS !"
        print*, "-----------------------------"
	print*

	close(100)
	close(101)
	close(102)

!---------------------------------------------------------------------------------------

contains
  
!--------------------------------------------------------------------------------------------

   subroutine get_dist(time1, time2, atom, dist)
        
        use READER

        implicit none

        integer          :: time1, time2, atom
        double precision :: dx, dy, dz, dist	

        dx=x(time1,atom)-x(time2,atom)
        dy=y(time1,atom)-y(time2,atom)
        dz=z(time1,atom)-z(time2,atom)
	dist=dsqrt(dx*dx+dy*dy+dz*dz)

   end subroutine get_dist

!--------------------------------------------------------------------------------------------

end program self_part_van_hove_correlation_function

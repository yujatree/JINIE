program hop_function

        
! Variables Declaration-----------------------------------------------------------------

        use READER

        implicit none

        integer 				      :: deltaT, GAP
	double precision 			      :: h_c, &
							 norm, dA, dB, dx, dy, dz

	type :: ptr
	   double precision, pointer, dimension(:,:)  :: p
	end type

	type(ptr), dimension(3)    		      :: posA, posB, avgPosA, avgPosB
	double precision, allocatable, dimension(:,:) :: hop
        integer, allocatable, dimension(:,:)          :: event

	character(len=100)		              :: num_str, fmt_F, fmt_I

! Setting-------------------------------------------------------------------------------

	call dcd_reader

	open(50, file='/mnt/CALC_SSD/jina/CODES/JINIE/inp/hop.inp', status='old')
	open(101, file='HOP.out', status='unknown')
	open(102, file='EVENT.out', status='unknown')

        write(num_str, '(I0)') num_Li
	fmt_F = '(I0,X,'//trim(num_str)//'(F12.6,X))'
	fmt_I = '(I0,X,'//trim(num_str)//'(I0,X))'

998	format(A,I5,A)

	read(50,*) deltaT
	read(50,*) h_c
	print*, 'deltaT is', deltaT
	print*, 'h_c is', h_c

	do i=1,3
	   allocate(posA(i)%p(total_traj,num_Li),posB(i)%p(total_traj,num_Li),&
		    avgPosA(i)%p(total_traj,num_Li),avgPosB(i)%p(total_traj,num_Li))
	enddo
	allocate(hop(total_traj,num_Li))
	allocate(event(total_traj,num_Li))
	
	print*, "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
	print*, ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
	print*
        print*, "-----------------------------"
	print*, "(  HOP ANALYSIS  )> \(o_o ;;)"
        print*, "-----------------------------"
	print*

	!---------------------------------------------------------------------------------------
	
	!Initialize Averaged Position : <r_i(t)>_A,B

	GAP=deltaT/2 ! hope 'deltaT' value is even...

        norm=dble(GAP+1)

        do i=1, 3
           posA(i)%p=0.d0
           posB(i)%p=0.d0
        enddo

	do i=1+GAP+GAP, total_traj-GAP-GAP

	   do j=1, num_Li
	     
	      do k=0, GAP  ! Or GAP and divide with GAP+1
		 posA(1)%p(i,j)=posA(1)%p(i,j)+x(i-k,j)
		 posA(2)%p(i,j)=posA(2)%p(i,j)+y(i-k,j)
		 posA(3)%p(i,j)=posA(3)%p(i,j)+z(i-k,j)
		 posB(1)%p(i,j)=posB(1)%p(i,j)+x(i+k,j)
		 posB(2)%p(i,j)=posB(2)%p(i,j)+y(i+k,j)
		 posB(3)%p(i,j)=posB(3)%p(i,j)+z(i+k,j)
	      enddo
	      do k=1, 3
		 avgPosA(k)%p(i,j)=posA(k)%p(i,j)/norm
		 avgPosB(k)%p(i,j)=posB(k)%p(i,j)/norm
	      enddo
	   enddo 
	enddo

	do i=1, 3
	   deallocate(posA(i)%p, posB(i)%p)
	enddo

	! Hop function
	
	do i=1+GAP+GAP, total_traj-GAP-GAP

	   if (mod(i,100)==0) write(*,998) " Processing ", i, " th trajectory"

	   do j=1, num_Li

	      dA=0.d0 ; dB=0.d0

	      do k=0, GAP
		 dx=x(i+k,j)-avgPosA(1)%p(i+k,j)
		 dy=y(i+k,j)-avgPosA(2)%p(i+k,j)
		 dz=z(i+k,j)-avgPosA(3)%p(i+k,j)
		 dB=dB+dx*dx+dy*dy+dz*dz
		
		 dx=x(i-k,j)-avgPosB(1)%p(i-k,j)
		 dy=y(i-k,j)-avgPosB(2)%p(i-k,j)
		 dz=z(i-k,j)-avgPosB(3)%p(i-k,j)
		 dA=dA+dx*dx+dy*dy+dz*dz
	      enddo
	
	      hop(i,j)=sqrt(dA*dB/norm/norm)
	      event(i,j)=-1*(hop(i,j).gt.h_c)

	   enddo

	   write(101,fmt=fmt_F,advance='no') i, (real(hop(i,j)), j=1,num_Li)
	   write(101,*)
	   write(102,fmt=fmt_I,advance='no') i, (event(i,j), j=1,num_Li)
	   write(102,*)

	enddo

	!---------------------------------------------------------------------------------------

	print*
	print*, "   Saved as 'HOP.out' file  "
	print*
	close(50) ; close(101) ; close(102)

!---------------------------------------------------------------------------------------

end program hop_function

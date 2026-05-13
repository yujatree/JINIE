program conditional_radial_distribution_function 

! Variables Declaration-----------------------------------------------------------------
        
        use READER

        implicit none

        logical					      :: file_exists

        integer 			              :: bin, bin_num, hop_traj, &
							 htype, idx1, idx2
        integer, dimension(:,:), allocatable	      :: hopidx
        double precision	                      :: dist, bin_size, delta_t, r
        double precision, parameter 	              :: pi=datan(1.d0)*4.d0
        double precision, dimension(6)		      :: count_hidx
        double precision, dimension(:,:), allocatable :: hop, g_Li, g_P, g_S, &
						         g_Li_tot, g_P_tot, g_S_tot
        character(10)				      :: filenum

! Setting-------------------------------------------------------------------------------

        call dcd_reader

999	format(F7.3," ",6(G0,1X)/)
998	format(A,I5,A)

	open(50, file='/mnt/CALC_SSD/jina/CODES/JINIE/inp/bin.inp', status='old')
	open(51, file='/mnt/CALC_SSD/jina/CODES/JINIE/inp/hop.inp', status='old')

	do nfile=1,num_of_args
	   write(filenum,'(I1)') nfile 

	   inquire(file=trim(adjustl(filenum))//"/HOP.out", exist=file_exists)

	   if (file_exists) then
	      open(52+nfile, file=trim(adjustl(filenum))//'/HOP.out', status='old')
	   else
	      open(52+nfile, file=trim(adjustl(filenum))//'/OUTPUT/HOP.out', status='old')
	   endif
	enddo 
   
	open(100, file='CRDF_Li.out', status='unknown')
	open(101, file='CRDF_P.out', status='unknown')
	open(102, file='CRDF_S.out', status='unknown')

	do i=100,102
	   write(i,*) "# r(A) 0.0-1.0 1.0-2.0 2.0-3.0 3.0-4.0 4.0-5.0 total"
	enddo

	print*
        print*, "-----------------------------"
        print*, "        CRDF START !"
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
	print*, " Delta t : "
	write(*,'(27X,I2)'), int(delta_t)

	hop_traj=total_traj-4*int(delta_t/2.d0)

	print*, " Hop traj :", hop_traj

        allocate(g_Li(6,bin),g_P(6,bin),g_S(6,bin)) ! total, 0th, 1st, 2nd, 3rd, 4th > 6
	allocate(g_Li_tot(6,bin),g_P_tot(6,bin),g_S_tot(6,bin)) 
	allocate(hop(hop_traj,num_Li),hopidx(hop_traj,num_Li))

!---------------------------------------------------------------------------------------

	g_Li_tot=0.d0 ; g_P_tot=0.d0 ; g_S_tot=0.d0

	do nfile=1, num_of_args

	   g_Li=0.d0 ; g_P=0.d0 ; g_S=0.d0

	   count_hidx=0.d0

	   read(52+nfile,*)

	   do i=1, hop_traj

	      if (mod(i,10)==0) write(*,998) " Processing ", i, " th trajectory"

	      read(52+nfile,*) dummyi(1), (hop(i,j), j=1,num_Li)

	      do j=1, num_Li
	         if (hop(i,j)>=5) hop(i,j)=4.d0
	         htype=int(hop(i,j))+1 ! if 0.0-1.0:1, 1.0-2.0:2, ... g_Li(6,:) is total!
	         count_hidx(htype)=count_hidx(htype)+1
	         ! Li-Li
	         do k=1, num_Li
	   	 if (j==k) goto 100
	            call pbc_dist(nfile, i, j, k, dist)
	            bin_num=int(dist/bin_size)+1
	            g_Li(htype,bin_num)=g_Li(htype,bin_num)+1.d0
	            g_Li(6,bin_num)=g_Li(6,bin_num)+1.d0
100	         enddo
	         ! Li-P
	         do k=1, num_P
	            idx1=num_Li+k ! To match P idx
	            call pbc_dist(nfile, i, j, idx1, dist)
	            bin_num=int(dist/bin_size)+1
	            g_P(htype,bin_num)=g_P(htype,bin_num)+1.d0
	            g_P(6,bin_num)=g_P(6,bin_num)+1.d0
	         enddo
	         ! Li-S
	         do k=1, num_S
	            idx1=num_Li+num_P+k ! To match S idx
	            call pbc_dist(nfile, i, j, idx1, dist)
	            bin_num=int(dist/bin_size)+1
	            g_S(htype,bin_num)=g_S(htype,bin_num)+1.d0
	            g_S(6,bin_num)=g_S(6,bin_num)+1.d0
	         enddo
	      enddo
	   enddo

	   g_Li=g_Li*ens(nfile)%vol(total_traj)
	   g_P=g_P*ens(nfile)%vol(total_traj)
	   g_S=g_S*ens(nfile)%vol(total_traj)

	   count_hidx(6)=dble(hop_traj*num_Li)

	   do i=1, bin
	         
	      r=(dble(i)-0.5d0)*bin_size

	      do j=1,6
	         if (count_hidx(j).ne.0.d0) then
	            g_Li(j,i)=g_Li(j,i)/count_hidx(j)/dble(num_Li)
	            g_P(j,i)=g_P(j,i)/count_hidx(j)/dble(num_P)
	            g_S(j,i)=g_S(j,i)/count_hidx(j)/dble(num_S)
	         else
	   	 g_Li(j,i)=0.d0 ; g_P(j,i)=0.d0 ; g_S(j,i)=0.d0
	         endif
	      enddo
	   enddo
	   
	   g_Li_tot=g_Li_tot+g_Li
	   g_P_tot =g_P_tot +g_P
	   g_S_tot =g_S_tot +g_S

	enddo

	g_Li_tot=g_Li_tot/dble(num_of_args)
	g_P_tot=g_P_tot/dble(num_of_args)
	g_S_tot=g_S_tot/dble(num_of_args)

	do i=1, bin

	   r=(dble(i)-0.5d0)*bin_size

	   do j=1,6
	      if (i==1) then
	         g_Li_tot(j,i)=0.d0
	         g_P_tot(j,i)=0.d0
	         g_S_tot(j,i)=0.d0
	      else
	         g_Li_tot(j,i)=g_Li_tot(j,i)/(4.d0*pi*r*r*bin_size)
	         g_P_tot(j,i)=g_P_tot(j,i)/(4.d0*pi*r*r*bin_size)
	         g_S_tot(j,i)=g_S_tot(j,i)/(4.d0*pi*r*r*bin_size)
	      endif
	   enddo

	   write(100,999,advance='no') r, (g_Li_tot(j,i), j=1,6)
	   write(101,999,advance='no') r, (g_P_tot(j,i), j=1,6)
	   write(102,999,advance='no') r, (g_S_tot(j,i), j=1,6)

	enddo


	print*
	print*, "   Saved as 'CRDF.out' file   "
        print*, "-----------------------------"
	print*, "          SUCCESS !"
        print*, "-----------------------------"
	print*

	do nfile=1, num_of_args
	   close(52+nfile)
	enddo
	close(100) ; close(101) ; close(102)

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

end program conditional_radial_distribution_function

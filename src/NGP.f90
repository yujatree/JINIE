program non_gaussian_parameter

! Variables Declaration-----------------------------------------------------------------
        
        use READER

        implicit none
 
        double precision	            :: displ
        double precision, dimension(4)      :: msd, mqd, ngp

! Setting-------------------------------------------------------------------------------

        call dcd_reader

999	format(I6," ",4(G0,1X)/)
998	format(A,I4,A)

	open(100, file='NGP.out', status='unknown')
	write(100,*) "# time(ps) Li P S Total"	

	print*
        print*, "-----------------------------"
        print*, "         NGP START !"
        print*, "-----------------------------"
	print*

!---------------------------------------------------------------------------------------

	do i=1, int(total_traj/3)
	  
	   if (mod(i,10)==0) write(*,998) " Processing ", i, " th trajectory"

	   msd=0.d0 ; mqd=0.d0

	   do j=1,total_traj-i

	      do k=1, num_Li
	         call sq_displ(i+j, k, j, k, displ)
	         msd(1)=msd(1)+displ
	         mqd(1)=mqd(1)+displ*displ
	      enddo
	      do k=1, num_P
	         atom_idx=num_Li+k
	         call sq_displ(i+j, atom_idx, j, atom_idx, displ)
	         msd(2)=msd(2)+displ
	         mqd(2)=mqd(2)+displ*displ
	      enddo
	      do k=1, num_S
	         atom_idx=num_Li+num_P+k
	         call sq_displ(i+j, atom_idx, j, atom_idx, displ)
	         msd(3)=msd(3)+displ
	         mqd(3)=mqd(3)+displ*displ
	      enddo
	      msd(4)=msd(1)+msd(2)+msd(3)
	      mqd(4)=mqd(1)+mqd(2)+mqd(3)

	   enddo

	   msd(1)=msd(1)/dble(total_traj-i)/dble(num_Li)
	   msd(2)=msd(2)/dble(total_traj-i)/dble(num_P)
	   msd(3)=msd(3)/dble(total_traj-i)/dble(num_S)
	   msd(4)=msd(4)/dble(total_traj-i)/dble(nmedia)
	   mqd(1)=mqd(1)/dble(total_traj-i)/dble(num_Li)
	   mqd(2)=mqd(2)/dble(total_traj-i)/dble(num_P)
	   mqd(3)=mqd(3)/dble(total_traj-i)/dble(num_S)
	   mqd(4)=mqd(4)/dble(total_traj-i)/dble(nmedia)

	   do j=1,4
	      ngp(j)=0.6d0*mqd(j)/msd(j)/msd(j)-1.d0
	   enddo

	   write(100,999,advance='no') i, ngp(1), ngp(2), ngp(3), ngp(4)

	enddo      

	print*
	print*, "   Saved as 'NGP.out' file   "
        print*, "-----------------------------"
	print*, "          SUCCESS !"
        print*, "-----------------------------"
	print*

	close(100)

!---------------------------------------------------------------------------------------

contains
  
!--------------------------------------------------------------------------------------------

   subroutine sq_displ(time1, atom1, time2, atom2, displ)
        
        use READER

        implicit none

        integer          :: time1, atom1, time2, atom2
        double precision :: dx, dy, dz, displ	

        dx=x(time1,atom1)-x(time2,atom2)
        dy=y(time1,atom1)-y(time2,atom2)
        dz=z(time1,atom1)-z(time2,atom2)
	displ=dx*dx+dy*dy+dz*dz

   end subroutine sq_displ

!--------------------------------------------------------------------------------------------

end program non_gaussian_parameter

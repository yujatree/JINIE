program mean_squared_displacement

! Variables Declaration-----------------------------------------------------------------
        
        use READER

        implicit none
 
        double precision	            :: displ
        double precision, dimension(4)      :: msd

! Setting-------------------------------------------------------------------------------

        call dcd_reader

999	format(I6," ",4(G0,1X)/)
998	format(A,I4,A)

	open(100, file='MSD.out', status='unknown')
	write(100,*) "# time(ps) Li P S Total"	

	print*
        print*, "-----------------------------"
        print*, "         MSD START !"
        print*, "-----------------------------"
	print*

!---------------------------------------------------------------------------------------

	do i=1, int(total_traj/3)
	  
	   if (mod(i,10)==0) write(*,998) " Processing ", i, " th trajectory"

	   msd=0.d0

	   do j=1,total_traj-i

	      do k=1, num_Li
	         call sq_displ(i+j, k, j, k, displ)
	         msd(1)=msd(1)+displ
	      enddo
	      do k=1, num_P
	         atom_idx=num_Li+k
	         call sq_displ(i+j, atom_idx, j, atom_idx, displ)
	         msd(2)=msd(2)+displ
	      enddo
	      do k=1, num_S
	         atom_idx=num_Li+num_P+k
	         call sq_displ(i+j, atom_idx, j, atom_idx, displ)
	         msd(3)=msd(3)+displ
	      enddo
	      msd(4)=msd(1)+msd(2)+msd(3)

	   enddo

	   msd(1)=msd(1)/dble(total_traj-i)/dble(num_Li)
	   msd(2)=msd(2)/dble(total_traj-i)/dble(num_P)
	   msd(3)=msd(3)/dble(total_traj-i)/dble(num_S)
	   msd(4)=msd(4)/dble(total_traj-i)/dble(nmedia)
				      
	   write(100,999,advance='no') i, msd(1), msd(2), msd(3), msd(4)

	enddo      

	print*
	print*, "   Saved as 'MSD.out' file   "
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

end program mean_squared_displacement

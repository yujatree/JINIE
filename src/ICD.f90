program ionic_conductivity

! Variables Declaration-----------------------------------------------------------------
        
        use READER

        implicit none
 
        double precision	    :: displ, msd, const, D, icd
        double precision, parameter :: q=1.60217663d-19, &
				       kb=1.380649d-23,&             ! Constant N-E eq
				       temp=300.d0, &	             ! Default : 300 K
				       A_to_cm=1.d-8, ps_to_s=1.d-12 ! Optimized to ps

! Setting-------------------------------------------------------------------------------

        call dcd_reader

999	format(I6," ",G0/)
998	format(A,I4,A)

	open(100, file='ICD.out', status='unknown')

	print*
        print*, "-----------------------------"
        print*, "         ICD START !"
        print*, "-----------------------------"
	print*

!---------------------------------------------------------------------------------------

	! READ MSD.out
        do i=1, int(total_traj/3)
           if (mod(i,10)==0) write(*,998) " Processing ", i, " th trajectory"
           msd=0.d0
           do j=1,total_traj-i
              do k=1, num_Li
                 call sq_displ(i+j, k, j, k, displ)
                 msd=msd+displ
              enddo
           enddo
           msd=msd/dble(total_traj-i)/dble(num_Li)
           write(100,999,advance='no') i, msd
        enddo

	volume(i)=volume(i)*A_to_cm*A_to_cm*A_to_cm

	const=(dble(num_Li)/volume(i))*q*q/kb/temp ! (C^2/J.cm^3)
	D=msd/6.d0/dble(i)			   ! (A^2/ps)
	D=D*A_to_cm*A_to_cm/ps_to_s		   ! (cm^2/s)
        icd=const*D				   ! (A/J.cm=S/cm)
       print*,const 
	write(100,'(A," ",G0," ",A)') 'Diffusion coefficient :', D, '(cm^2/s)'
	write(100,'(A," ",G0," ",A)') 'Ionic conductivity    :', icd, '(S/cm)'

	print*
	print*, "   Saved as 'ICD.out' file   "
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

end program ionic_conductivity

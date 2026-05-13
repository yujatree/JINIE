program ionic_conductivity

! Variables Declaration-----------------------------------------------------------------

        use READER

        implicit none

        integer                         :: t_lag
        double precision	            :: msd, dummy, const, D, icd
        double precision, parameter :: q=1.60217663d-19, &
				       kb=1.380649d-23,             &
				       temp=300.d0, 	             &
				       A_to_cm=1.d-8, ps_to_s=1.d-12
        logical :: ex

! Setting-------------------------------------------------------------------------------

        call dcd_reader

998	format(A,I5,A)

        inquire(file='MSD.out', exist=ex)

        if (ex) then
           open(50, file='MSD.out', status='old', action='read')

        else
           inquire(file='OUTPUT/MSD.out', exist=ex)

           if (ex) then
              open(50, file='OUTPUT/MSD.out', status='old', action='read')

           else
              print*, ' ERROR: MSD.out not found.'
	      print*
              stop
           end if
        end if

	open(100, file='ICD.out', status='unknown')

	print*
        print*, "-----------------------------"
        print*, "         ICD START !"
        print*, "-----------------------------"
	print*

!---------------------------------------------------------------------------------------

        ! Read Li MSD from MSD.out
        read(50,*)
        do i=1, int(total_traj/3)
           if (mod(i,10)==0) write(*,998) " Reading  ", i, " th trajectory"
           read(50,*) t_lag, msd, dummy, dummy, dummy
        enddo

        close(50)

	volume(total_traj)=volume(total_traj)*A_to_cm*A_to_cm*A_to_cm
	const=(dble(num_Li)/volume(total_traj))*q*q/kb/temp
	D=msd/6.d0/dble(t_lag)
	D=D*A_to_cm*A_to_cm/ps_to_s
        icd=const*D

	write(100,'(A," ",G0," ",A)') 'Diffusion coefficient :', D,   '(cm^2/s)'
	write(100,'(A," ",G0," ",A)') 'Ionic conductivity    :', icd, '(S/cm)'

	print*
	print*, "   Saved as 'ICD.out' file   "
        print*, "-----------------------------"
	print*, "          SUCCESS !"
        print*, "-----------------------------"
	print*

	close(100)

!---------------------------------------------------------------------------------------

end program ionic_conductivity

program survival_probability

! Variables Declaration-----------------------------------------------------------------

        use READER

        implicit none

        logical                                       :: ex

        integer                                       :: deltaT, count_h, dt
        double precision                              :: mean, sqmean, kfast, t, &
                                                         cumm, sum_crt

        integer, allocatable, dimension(:,:)          :: hop
        integer, allocatable, dimension(:)            :: prev, Life
        double precision, allocatable, dimension(:)   :: Crt

        character(10)                                 :: filenum

! Setting-------------------------------------------------------------------------------

        call dcd_reader

999	format(F7.3," ",*(G0,1X),/)
998	format(A,I5,A)

	open(50, file='/mnt/CALC_SSD/jina/CODES/JINIE/inp/hop.inp', status='old')

	do nfile=1,num_of_args
	   write(filenum,'(I1)') nfile

	   inquire(file=trim(adjustl(filenum))//"/EVT.out", exist=ex)

	   if (ex) then
	      open(52+nfile, file=trim(adjustl(filenum))//'/EVT.out', status='old')
	   else
	      open(52+nfile, file=trim(adjustl(filenum))//'/OUTPUT/EVT.out', status='old')
	   endif
	enddo
	open(100, file='CST_Lifetime.out', status='unknown')
	open(101, file='CST_Survival.out', status='unknown')
	open(102, file='CST_Randomness.out', status='unknown')

	print*
        print*, "-----------------------------"
        print*, "       C_S(t) START !"
        print*, "-----------------------------"
	print*
	read(50,*) deltaT
	write(*,'(14X,"deltaT :",I5)') deltaT
	print*
	close(50)

	allocate(hop(total_traj, num_Li))
	allocate(prev(num_Li))
	allocate(Life(total_traj), Crt(total_traj))

!-Lifetime Distribution-----------------------------------------------------------------

	Life=0 ; count_h=0

	do nfile=1,num_of_args

 	   hop=0

	   read(52+nfile,*)  ! skip header

           do i=1+deltaT, total_traj-deltaT
	      read(52+nfile,*) dummyi(1), (hop(i,j), j=1,num_Li)
	   enddo

	   prev(:)=1+deltaT

	   do i=1+deltaT, total_traj-deltaT
	      do j=1,num_Li
	         if (hop(i,j)==1) then
	   	    count_h=count_h+1
	   	    dt=i-prev(j)
	   	    prev(j)=i
	   	    if (dt >= 1) Life(dt)=Life(dt)+1
	         endif
	      enddo
	   enddo

	   close(52+nfile)

	enddo


!-Randomness parameter------------------------------------------------------------------

	mean=0.d0
	sqmean=0.d0

	do i=1,total_traj

	   t=dble(i)-0.5d0
	   write(100,*) t, dble(Life(i))/dble(count_h)

	   mean=mean+t*dble(Life(i))
	   sqmean=sqmean+t*t*dble(Life(i))

	enddo

	kfast=dble(count_h)/dble(total_traj-2*deltaT)/dble(num_Li)
	mean=mean/dble(count_h)
	sqmean=sqmean/dble(count_h)

	write(102,*) "Number of sampled lifetimes: ", count_h
	write(102,*) "K_fast: ", kfast
	write(102,*) "<t>: ", mean
	write(102,*) "<t^2>: ", sqmean
	write(102,*) "stdev: ", sqrt(sqmean-mean*mean)
	write(102,*) "R: ", sqmean/mean/mean-1.d0

!-Survival probability------------------------------------------------------------------

	cumm=0.d0

	do i=1,total_traj
	   cumm=cumm+dble(Life(i))/dble(count_h)
	   Crt(i)=1.d0-cumm
	enddo

	cumm=0.d0
	sum_crt=sum(Crt)
	do i=1,total_traj
	   cumm=cumm+Crt(i)
	   write(101,*) dble(i), Crt(i), (sum_crt-cumm)/mean
	enddo

	print*
        print*, "   Saved as 'CST_life/surv/rand.out' files  "
        print*, "-----------------------------"
	print*, "          SUCCESS !"
        print*, "-----------------------------"
	print*

	close(100) ; close(101) ; close(102)

!--------------------------------------------------------------------------------------------

end program survival_probability

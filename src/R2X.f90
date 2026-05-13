program rdf_to_structure_factor_xray

! Variables Declaration-----------------------------------------------------------------
  
        use READER
        
        implicit none

        integer :: bin_r, bin_q, a, bidx
        double precision :: bin_size_r, bin_size_q, &
                            rho, kernel, sumr, norm, sa, factor
        double precision, parameter :: pi=datan(1.d0)*4.d0
        double precision :: fx(3), c(3)

        double precision, allocatable, dimension(:)   :: r, q, S
        double precision, allocatable, dimension(:,:) :: gij, Sij

        logical :: ex

! Setting-------------------------------------------------------------------------------

        call dcd_reader

        inquire(file='RDF.out', exist=ex)
      
        if (ex) then
           open(100, file='RDF.out', status='old', action='read')
      
        else
           inquire(file='OUTPUT/RDF.out', exist=ex)
      
           if (ex) then
              open(100, file='OUTPUT/RDF.out', status='old', action='read')
      
           else
              print*, ' ERROR: RDF.out not found.'
	      print*
              stop
           end if
        end if

	open(50, file='/mnt/CALC_SSD/jina/CODES/JINIE/inp/bin.inp', status='old')
	read(50,*) bin_size_r
	read(50,*) bin_r
        read(50,*)
	read(50,*) bin_size_q
	read(50,*) bin_q
	open(101, file='XSF.out', status='unknown')

	print*
        print*, "-----------------------------"
        print*, "         R2X START !"
        print*, "-----------------------------"
	print*, " HOPE THAT U CALCULATE THIS"
	print*, "   WITH STATIC CELL (NO NPT)!"
	print*
	print*, " Size/num of bins of r :"
	write(*,'(22X,F7.3)'), bin_size_r
	write(*,'(22X,I7)'), bin_r
	print*, " Size/num of bins of q :"
	write(*,'(22X,F7.3)'), bin_size_q
	write(*,'(22X,I7)'), bin_q
	print*

	allocate(r(bin_r), q(bin_q), gij(6,bin_r), Sij(6,bin_q), S(bin_q))
	gij = 1.d0
	r   = 0.d0

	read(100,*)
	do i=1, bin_r
	   read(100,*) r(i), (gij(j,i), j=1,6), dummyr
	   if (r(i)>19.5) goto 1000
	enddo
1000	continue

!---------------------------------------------------------------------------------------

        rho = dble(nmedia)/volume(1)
        
        c(1)=dble(num_Li)/dble(nmedia)
        c(2)=dble(num_P)/dble(nmedia)
        c(3)=dble(num_S)/dble(nmedia)

        do i=1,bin_q
           q(i) = dble(i) * bin_size_q
        enddo

        do i=1,6
           do j=1,bin_q

              sumr = 0.d0
              do k=1,bin_r
                 if (q(j)*r(k)>1.d-8) then
                    kernel = sin(q(j)*r(k)) / (q(j)*r(k))
                 else
                    kernel = 1.d0
                 endif
                 sumr = sumr + 4.d0*pi*r(k)**2.d0 &
                        * (gij(i,k)-1.d0) * kernel
              enddo

              Sij(i,j) = 1.d0 + rho * sumr * bin_size_r

           enddo
        enddo
        
!---------------------------------------------------------------------------------------

        S = 0.d0

        do j=1,bin_q

           sa = q(j) / (4.d0*pi)

           ! Li (Z=3)
           fx(1) = 0.413048d0*exp(-15.801d0*sa*sa) &
                 + 0.294953d0*exp(-6.647d0*sa*sa)  &
                 + 0.187491d0*exp(-2.001d0*sa*sa)  &
                 + 0.080701d0*exp(-0.312d0*sa*sa)  &
                 + 0.023736d0

           ! P (Z=15)
           fx(2) = 1.11427d0*exp(-15.034d0*sa*sa) &
                 + 0.83547d0*exp(-4.499d0*sa*sa)  &
                 + 0.36278d0*exp(-1.328d0*sa*sa)  &
                 + 0.18530d0*exp(-0.240d0*sa*sa)  &
                 + 0.03874d0

           ! S (Z=16)
           fx(3) = 0.96714d0*exp(-13.068d0*sa*sa) &
                 + 0.77253d0*exp(-3.753d0*sa*sa)  &
                 + 0.29792d0*exp(-1.143d0*sa*sa)  &
                 + 0.16535d0*exp(-0.217d0*sa*sa)  &
                 + 0.03016d0

           norm = ( c(1)*fx(1) + c(2)*fx(2) + c(3)*fx(3) )**2
           S(j) = 1.d0

           do i=1,6
              select case(i)
              case(1); a=1; bidx=1; factor=1.d0   ! Li-Li
              case(2); a=2; bidx=2; factor=1.d0   ! P-P
              case(3); a=3; bidx=3; factor=1.d0   ! S-S
              case(4); a=1; bidx=2; factor=2.d0   ! Li-P  (cross: x2)
              case(5); a=1; bidx=3; factor=2.d0   ! Li-S  (cross: x2)
              case(6); a=2; bidx=3; factor=2.d0   ! P-S   (cross: x2)
              end select

              S(j) = S(j) &
                   + factor * ( c(a)*c(bidx)*fx(a)*fx(bidx) / norm ) &
                   * ( Sij(i,j) - 1.d0 )
           enddo

        enddo

! ----------------------------------------------------------------------
! Output
! ----------------------------------------------------------------------

        do i=1,bin_q
           write(101,*) q(i), S(i)
        enddo

        close(100)
        close(101)
        close(50)

	print*
	print*, "   Saved as 'XSF.out' file   "
        print*, "-----------------------------"
	print*, "          SUCCESS !"
        print*, "-----------------------------"
	print*

end program rdf_to_structure_factor_xray

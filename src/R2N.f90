program rdf_to_neutron_structure_factor

! Variables Declaration-----------------------------------------------------------------

        use READER

        implicit none

        integer 				      :: bin_r, bin_q

        double precision 			      :: bin_size_r, bin_size_q, &
							 rho, kernel, sumr
        double precision, parameter		      :: pi=datan(1.d0)*4.d0

        double precision           		      :: c(3), b(3), cc(6), bb(6)
				  			       ! c: atomic fraction,
							       ! b: neutron scattering lengths
	double precision, allocatable, dimension(:)   :: r, q, S
	double precision, allocatable, dimension(:,:) :: gij, Sij

        logical 				      :: ex

! Setting-------------------------------------------------------------------------------

        call dcd_reader ! to get volume ...

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
	open(101, file='NSF.out', status='unknown')

	print*
        print*, "-----------------------------"
        print*, "         R2N START !"
        print*, "-----------------------------"
	print*, " HOPE THAT U CALCULATE THIS"
	print*, "   WITH STATIC CELL (NO NPT)!"
	print*
	print*, " Size/num of bins of r :"
	write(*,'(22X,F7.3)') bin_size_r
	write(*,'(22X,I7)') bin_r
	print*, " Size/num of bins of q :"
	write(*,'(22X,F7.3)') bin_size_q
	write(*,'(22X,I7)') bin_q
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

	rho=dble(nmedia)/volume(1)

	c(1)=dble(num_Li)/dble(nmedia)
	c(2)=dble(num_P)/dble(nmedia)
	c(3)=dble(num_S)/dble(nmedia)
	cc(1)=c(1)*c(1)
	cc(2)=c(2)*c(2)
	cc(3)=c(3)*c(3)
	cc(4)=2.d0*c(1)*c(2)
	cc(5)=2.d0*c(1)*c(3)
	cc(6)=2.d0*c(2)*c(3)

        b(1)=-1.90d0  ! Li (fm)
        b(2)= 5.13d0  ! P  (fm)
        b(3)= 2.847d0 ! S  (fm)
	bb(1)=b(1)*b(1)
	bb(2)=b(2)*b(2)
	bb(3)=b(3)*b(3)
	bb(4)=b(1)*b(2)
	bb(5)=b(1)*b(3)
	bb(6)=b(2)*b(3)

	do i=1, bin_q
	   q(i)=i*bin_size_q
	enddo

!---------------------------------------------------------------------------------------

        S = 0.d0

        do i=1, 6
           do j=1, bin_q

              sumr=0.d0

              do k=1, bin_r

                 if (q(j)*r(k)>1.d-8) then
                    kernel=sin(q(j)*r(k))/(q(j)*r(k))
                 else
                    kernel=1.d0
                 endif

                 sumr=sumr+4.d0*pi*r(k)**2.d0*(gij(i,k)-1.d0)*kernel

              enddo

              Sij(i,j)=1.d0+rho*sumr*bin_size_r

              S(j)=S(j)+cc(i)*bb(i)*(Sij(i,j)-1.d0)

           enddo
        enddo

        S=S+1.d0

	do i=1, bin_q
	   write(101,*) q(i), S(i)
	enddo

	print*
	print*, "   Saved as 'NSF.out' file   "
        print*, "-----------------------------"
	print*, "          SUCCESS !"
        print*, "-----------------------------"
	print*

	close(50)
	close(100)
	close(101)

!---------------------------------------------------------------------------------------

end program rdf_to_neutron_structure_factor

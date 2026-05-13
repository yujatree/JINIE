program rotational_non_gaussian_parameter

! Variables Declaration-----------------------------------------------------------------

        use READER

        implicit none

        double precision                              :: theta, msad, mqad, rngp
        double precision                              :: tmp_dx, tmp_dy, tmp_dz, dist
        double precision, parameter                   :: pi=datan(1.d0)*4.d0, &
                                                         P_S=2.5d0

        integer                                       :: pidx, sidx, bond_count

        integer, allocatable, dimension(:)            :: pcounter
        type :: ptr
           double precision, pointer, dimension(:,:)  :: dx, dy, dz
        end type
        type(ptr), allocatable, dimension(:)          :: t

! Setting-------------------------------------------------------------------------------

        call dcd_reader

999	format(I6," ",G0/)
998	format(A,I5,A)

	open(100, file='RNGP.out', status='unknown')
	write(100,*) "# time(steps) RNGP_P"

	print*
        print*, "-----------------------------"
        print*, "         RNGP START !"
        print*, "-----------------------------"
	print*

!-P-S bond vector initialization--------------------------------------------------------

	allocate(t(total_traj))
	do i=1, total_traj
	   allocate(t(i)%dx(num_P,4), t(i)%dy(num_P,4), t(i)%dz(num_P,4))
	   t(i)%dx=0.d0 ; t(i)%dy=0.d0 ; t(i)%dz=0.d0
	enddo
	allocate(pcounter(num_P))

	do i=1, total_traj
	   if (mod(i,100)==0) write(*,998) " Initializing ", i, " th traj"
	   pcounter=0
	   do j=1, num_P
	      pidx=num_Li+j
	      do k=1, num_S
	         sidx=num_Li+num_P+k
	         call pbc_dist(i, pidx, sidx, tmp_dx, tmp_dy, tmp_dz, dist)
	         if (dist < P_S) then
	            pcounter(j)=pcounter(j)+1
	            if (pcounter(j)>4) then
	               print*, " More than four S atoms around one P atom!"
	               stop
	            endif
	            t(i)%dx(j,pcounter(j))=tmp_dx
	            t(i)%dy(j,pcounter(j))=tmp_dy
	            t(i)%dz(j,pcounter(j))=tmp_dz
	         endif
	      enddo
	   enddo
	enddo

	print*, " Initialization done!"
	print*

!-Rotational NGP------------------------------------------------------------------------

	do i=1, int(total_traj/3)

	   if (mod(i,10)==0) write(*,998) " Processing ", i, " th trajectory"

	   msad=0.d0 ; mqad=0.d0 ; bond_count=0

	   do j=1, total_traj-i
	      do k=1, num_P
	         do l=1, 4
	            if ((t(j)%dx(k,l)==0.d0).or.(t(j+i)%dx(k,l)==0.d0)) cycle

	            call get_angle(t(j)%dx(k,l), t(j)%dy(k,l), t(j)%dz(k,l), &
	                           t(j+i)%dx(k,l), t(j+i)%dy(k,l), t(j+i)%dz(k,l), &
	                           theta)

	            msad=msad+theta*theta
	            mqad=mqad+theta*theta*theta*theta
	            bond_count=bond_count+1
	         enddo
	      enddo
	   enddo

	   if (bond_count > 0) then
	      msad=msad/dble(bond_count)
	      mqad=mqad/dble(bond_count)
	      rngp=mqad/3.d0/msad/msad-1.d0
	   else
	      rngp=0.d0
	   endif

	   write(100,999) i, rngp

	enddo

	print*
	print*, "   Saved as 'RNGP.out' file   "
        print*, "-----------------------------"
	print*, "          SUCCESS !"
        print*, "-----------------------------"
	print*

	close(100)

!---------------------------------------------------------------------------------------

contains

!--------------------------------------------------------------------------------------------

   subroutine pbc_dist(t, atom1, atom2, dx, dy, dz, dist)

        use READER

        implicit none

        integer          :: t, atom1, atom2
        double precision :: dx, dy, dz, dist

        dx=x(t,atom1)-x(t,atom2)
        dy=y(t,atom1)-y(t,atom2)
        dz=z(t,atom1)-z(t,atom2)
	dx=dx-nint(dx/box(t,1))*box(t,1)
	dy=dy-nint(dy/box(t,3))*box(t,3)
	dz=dz-nint(dz/box(t,6))*box(t,6)
	dist=dsqrt(dx*dx+dy*dy+dz*dz)

   end subroutine pbc_dist

!--------------------------------------------------------------------------------------------

   subroutine get_angle(r0_dx, r0_dy, r0_dz, r1_dx, r1_dy, r1_dz, theta)

        implicit none

        double precision :: r0_dx, r0_dy, r0_dz, r1_dx, r1_dy, r1_dz, &
                            inner, r0, r1, theta

        inner=r0_dx*r1_dx+r0_dy*r1_dy+r0_dz*r1_dz
        r0=dsqrt(r0_dx*r0_dx+r0_dy*r0_dy+r0_dz*r0_dz)
        r1=dsqrt(r1_dx*r1_dx+r1_dy*r1_dy+r1_dz*r1_dz)
        theta=dacos(inner/(r0*r1))

   end subroutine get_angle

!--------------------------------------------------------------------------------------------

end program rotational_non_gaussian_parameter

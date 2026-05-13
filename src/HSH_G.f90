program hop_shell_glass

! Variables Declaration-----------------------------------------------------------------

        use READER

        implicit none

        integer                         :: deltaT, GAP, nsurr
        integer                         :: bin_num, shl, t_lag
        double precision                :: h_c, binsize, dist, h_val
        double precision, parameter     :: max_h=5.d0
        integer, parameter              :: num_bins=1000, nshell=6

        double precision, allocatable, dimension(:,:) :: hop
        double precision, allocatable, dimension(:,:) :: hs
        integer, allocatable, dimension(:,:)          :: npts

        double precision                :: dx, dy, dz
        logical                         :: ex

! Setting-------------------------------------------------------------------------------

        call dcd_reader

998	format(A,I5,A)

        inquire(file='HOP.out', exist=ex)

        if (ex) then
           open(50, file='HOP.out', status='old', action='read')

        else
           inquire(file='OUTPUT/HOP.out', exist=ex)

           if (ex) then
              open(50, file='OUTPUT/HOP.out', status='old', action='read')

           else
              print*, ' ERROR: HOP.out not found.'
              print*
              stop
           end if
        end if

        open(51, file='/mnt/CALC_SSD/jina/CODES/JINIE/inp/hop.inp', status='old')
        open(100, file='HSH_G.out', status='unknown')

        read(51,*) deltaT
        read(51,*) h_c
        read(51,*) nsurr

        GAP     = deltaT / 2
        binsize = max_h / dble(num_bins)

        print*
        print*, "-----------------------------"
        print*, "        HSH_G START !"
        print*, "-----------------------------"
        print*
        write(*,'(14X,"deltaT :",I5)') deltaT
        write(*,'(14X,"h_c    :",F6.3)') h_c
        print*

!---------------------------------------------------------------------------------------

        allocate(hop(total_traj, num_Li))
        allocate(hs(nshell, num_bins), npts(nshell, num_bins))
        hop = 0.d0 ; hs = 0.d0 ; npts = 0

        do i = 1+GAP+GAP, total_traj-GAP-GAP
           if (mod(i,100)==0) write(*,998) " Reading  ", i, " th trajectory"
           read(50,*) t_lag, (hop(i,j), j=1,num_Li)
        enddo

        close(50)

        write(100,'(A)') "# h  total  1st  2nd  3rd  4th  5+"

        do i = 1+GAP+GAP, total_traj-GAP-GAP

           if (mod(i,100)==0) write(*,998) " Processing ", i, " th trajectory"

           do j = 1, num_Li

              if (hop(i,j) > max_h) cycle

              bin_num = int(hop(i,j)/binsize) + 1
              if (bin_num > num_bins) cycle

              do k = 1, num_Li
                 if (k == j) cycle
                 call pbc_dist(i, j, k, dx, dy, dz, dist)
                 shl = shell_g(dist)
                 hs(1, bin_num)   = hs(1, bin_num)   + hop(i,k)
                 npts(1, bin_num) = npts(1, bin_num)  + 1
                 if (shl > 0) then
                    hs(shl+1, bin_num)   = hs(shl+1, bin_num)   + hop(i,k)
                    npts(shl+1, bin_num) = npts(shl+1, bin_num) + 1
                 endif
              enddo

           enddo
        enddo

        do i = 1, num_bins
           h_val = (dble(i)-0.5d0) * binsize
           write(100,'(F7.4)',advance='no') h_val
           do j = 1, nshell
              if (npts(j,i) > 0) then
                 write(100,'(" ",G0)',advance='no') hs(j,i)/dble(npts(j,i))
              else
                 write(100,'(" 0.000000")',advance='no')
              endif
           enddo
           write(100,*)
        enddo

        print*
        print*, "   Saved as 'HSH_G.out' file  "
        print*, "-----------------------------"
        print*, "          SUCCESS !"
        print*, "-----------------------------"
        print*

        close(51) ; close(100)

!---------------------------------------------------------------------------------------

contains

!--------------------------------------------------------------------------------------------

   function shell_g(dist) result(shl)

        double precision, intent(in) :: dist
        integer :: shl

        if      (dist > 17.375d0) then ; shl = 0
        else if (dist > 14.375d0) then ; shl = 5
        else if (dist > 11.375d0) then ; shl = 4
        else if (dist >  8.225d0) then ; shl = 3
        else if (dist >  5.425d0) then ; shl = 2
        else                           ; shl = 1
        end if

   end function shell_g

!--------------------------------------------------------------------------------------------

   subroutine pbc_dist(time, atom1, atom2, dx, dy, dz, dist)

        use READER

        implicit none

        integer          :: time, atom1, atom2
        double precision :: dx, dy, dz, dist

        dx=x(time,atom1)-x(time,atom2)
        dy=y(time,atom1)-y(time,atom2)
        dz=z(time,atom1)-z(time,atom2)
        dx=dx-nint(dx/box(time,1))*box(time,1)
        dy=dy-nint(dy/box(time,3))*box(time,3)
        dz=dz-nint(dz/box(time,6))*box(time,6)
        dist=dsqrt(dx*dx+dy*dy+dz*dz)

   end subroutine pbc_dist

!--------------------------------------------------------------------------------------------

end program hop_shell_glass

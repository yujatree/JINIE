program neighbor_visualization

! Variables Declaration-----------------------------------------------------------------

        use READER

        implicit none

        logical                         :: ex

        integer                         :: deltaT, GAP, nsurr, id, idx1, m
        integer                         :: count_Li, count_P, count_S
        integer                         :: idx_tmp
        double precision                :: h_c, dist, dist_tmp

        integer, allocatable            :: sieve(:)
        integer, allocatable            :: event(:,:)
        double precision, allocatable   :: hop_row(:)
        integer, allocatable            :: assigned_Li(:,:), assigned_P(:,:), assigned_S(:,:)

        double precision, allocatable   :: dist_Li(:), dist_P(:), dist_S(:)
        integer, allocatable            :: idx_Li(:), idx_P(:), idx_S(:)

! Setting-------------------------------------------------------------------------------

        call dcd_reader

998	format(A,I5,A)

        open(50, file='/mnt/CALC_SSD/jina/CODES/JINIE/inp/hop.inp', status='old')

        inquire(file='HOP.out', exist=ex)

        if (ex) then
           open(51, file='HOP.out', status='old', action='read')

        else
           inquire(file='OUTPUT/HOP.out', exist=ex)

           if (ex) then
              open(51, file='OUTPUT/HOP.out', status='old', action='read')

           else
              print*, ' ERROR: HOP.out not found.'
              print*
              stop
           end if
        end if

        print*
        print*, "-----------------------------"
        print*, "        NVIS START !"
        print*, "-----------------------------"
        print*

        read(50,*) deltaT
        read(50,*) h_c
        read(50,*) nsurr

        close(50)

        GAP = deltaT / 2

        write(*,'(14X,"deltaT :",I5)') deltaT
        write(*,'(14X,"h_c    :",F6.3)') h_c
        write(*,'(14X,"nsurr  :",I5)') nsurr
        print*

!---------------------------------------------------------------------------------------

        allocate(sieve(num_Li))
        allocate(event(total_traj, num_Li))
        allocate(hop_row(num_Li))
        allocate(assigned_Li(total_traj, num_Li))
        allocate(assigned_P(total_traj, num_P))
        allocate(assigned_S(total_traj, num_S))

        allocate(dist_Li(num_Li), idx_Li(num_Li))
        allocate(dist_P(num_P),   idx_P(num_P))
        allocate(dist_S(num_S),   idx_S(num_S))

        sieve = 0 ; event = 0
        assigned_Li = 0 ; assigned_P = 0 ; assigned_S = 0

        read(51,*)  ! skip header

        do i = 1+GAP+GAP, total_traj-GAP-GAP
           read(51,*) dummyi(1), (hop_row(j), j=1,num_Li)
           do j = 1, num_Li
              if (hop_row(j) >= h_c) then
                 event(i,j) = 1
                 sieve(j)   = 1
              endif
           enddo
        enddo

        close(51)

!-Assign nearest neighbors for each hopping event--------------------------------------

        do j = 1, num_Li

           if (mod(j,10)==0) write(*,998) " Assigning Li", j, " ..."

           do i = 1, total_traj

              if (event(i,j) == 1) then

                 ! Li neighbors
                 count_Li = 0
                 do k = 1, num_Li
                    if (k==j) cycle
                    call pbc_dist(i, j, k, dist)
                    count_Li = count_Li + 1
                    dist_Li(count_Li) = dist
                    idx_Li(count_Li)  = k
                 enddo
                 call sort_insert(dist_Li, idx_Li, count_Li)
                 do l = 1, min(nsurr, count_Li)
                    k = idx_Li(l)
                    if (assigned_Li(i,k) /= 0) then
                       write(*,*) " OVERLAP: traj", i, " Li", k
                    else
                       assigned_Li(i,k) = l
                    endif
                 enddo

                 ! P neighbors
                 count_P = 0
                 do k = 1, num_P
                    idx1 = num_Li + k
                    call pbc_dist(i, j, idx1, dist)
                    count_P = count_P + 1
                    dist_P(count_P) = dist
                    idx_P(count_P)  = k
                 enddo
                 call sort_insert(dist_P, idx_P, count_P)
                 do l = 1, min(nsurr, count_P)
                    k = idx_P(l)
                    if (assigned_P(i,k) /= 0) then
                       write(*,*) " OVERLAP: traj", i, " P", k
                    else
                       assigned_P(i,k) = l
                    endif
                 enddo

                 ! S neighbors
                 count_S = 0
                 do k = 1, num_S
                    idx1 = num_Li + num_P + k
                    call pbc_dist(i, j, idx1, dist)
                    count_S = count_S + 1
                    dist_S(count_S) = dist
                    idx_S(count_S)  = k
                 enddo
                 call sort_insert(dist_S, idx_S, count_S)
                 do l = 1, min(nsurr, count_S)
                    k = idx_S(l)
                    if (assigned_S(i,k) /= 0) then
                       write(*,*) " OVERLAP: traj", i, " S", k
                    else
                       assigned_S(i,k) = l
                    endif
                 enddo

              endif
           enddo
        enddo

! writing Trajectory Data---------------------------------------------------------------

        open(100, file='neighbor.lammpstraj', status='unknown')

        print*
        write(*,'(14X,"Writing neighbor.lammpstraj ...")')
        print*

        do i = 1, total_traj

           if (mod(i,100)==0) write(*,998) " Writing   ", i, " th trajectory"

           write(100,'(A)') "ITEM: TIMESTEP"
           write(100,'(I8)') init_check + (i-1)*dt_check
           write(100,'(A)') "ITEM: NUMBER OF ATOMS"
           write(100,'(I8)') nmedia
           write(100,'(A)') "ITEM: BOX BOUNDS pp pp pp"
           write(100,'(E25.16," ",G0)') 0.d0, box(i,1)
           write(100,'(E25.16," ",G0)') 0.d0, box(i,3)
           write(100,'(E25.16," ",G0)') 0.d0, box(i,6)
           write(100,'(A)') "ITEM: ATOMS id type x y z sieve event assigned"

           do j = 1, nmedia

              if (j <= num_Li) then
                 id = 1
                 write(100,'(I4," ",I1," ",3(G0,1X),I1,1X,I1,1X,I0)') &
                       j, id, x(i,j), y(i,j), z(i,j), sieve(j), event(i,j), assigned_Li(i,j)

              else if (j <= num_Li+num_P) then
                 id = 2
                 write(100,'(I4," ",I1," ",3(G0,1X),I1,1X,I1,1X,I0)') &
                       j, id, x(i,j), y(i,j), z(i,j), 0, 0, assigned_P(i, j-num_Li)

              else
                 id = 3
                 write(100,'(I4," ",I1," ",3(G0,1X),I1,1X,I1,1X,I0)') &
                       j, id, x(i,j), y(i,j), z(i,j), 0, 0, assigned_S(i, j-num_Li-num_P)
              endif

           enddo

        enddo

        close(100)

        print*
        print*, "   Saved as 'neighbor.lammpstraj' file  "
        print*, "-----------------------------"
        print*, "          SUCCESS !"
        print*, "-----------------------------"
        print*

!---------------------------------------------------------------------------------------

contains

!--------------------------------------------------------------------------------------------

   subroutine sort_insert(darr, iarr, n)

        implicit none

        integer, intent(in)                            :: n
        double precision, dimension(:), intent(inout)  :: darr
        integer, dimension(:), intent(inout)           :: iarr

        integer          :: m, l, idx_tmp
        double precision :: dist_tmp

        do m = 2, n
           dist_tmp = darr(m)
           idx_tmp  = iarr(m)
           l = m-1
           do while (l >= 1 .and. darr(l) > dist_tmp)
              darr(l+1) = darr(l)
              iarr(l+1) = iarr(l)
              l = l-1
           enddo
           darr(l+1) = dist_tmp
           iarr(l+1) = idx_tmp
        enddo

   end subroutine sort_insert

!--------------------------------------------------------------------------------------------

   subroutine pbc_dist(time, atom1, atom2, dist)

        use READER

        implicit none

        integer          :: time, atom1, atom2
        double precision :: dx, dy, dz, dist

        dx = x(time,atom1) - x(time,atom2)
        dy = y(time,atom1) - y(time,atom2)
        dz = z(time,atom1) - z(time,atom2)
        dx = dx - nint(dx/box(time,1))*box(time,1)
        dy = dy - nint(dy/box(time,3))*box(time,3)
        dz = dz - nint(dz/box(time,6))*box(time,6)
        dist = dsqrt(dx*dx + dy*dy + dz*dz)

   end subroutine pbc_dist

!--------------------------------------------------------------------------------------------

end program neighbor_visualization

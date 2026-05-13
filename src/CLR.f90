program clear_trajectory

! Variables Declaration-----------------------------------------------------------------

        use READER

        implicit none

        logical                                       :: ex

        integer                                       :: deltaT, id, count_cage, cage_start
        double precision                              :: meanx, meany, meanz
        double precision, dimension(3)                :: tmp_pos

        integer, allocatable, dimension(:)            :: sieved
        integer, allocatable, dimension(:,:)          :: hop
        double precision, allocatable, dimension(:,:) :: cx, cy, cz

! Setting-------------------------------------------------------------------------------

        call dcd_reader

998	format(A,I5,A)

        open(50, file='/mnt/CALC_SSD/jina/CODES/JINIE/inp/hop.inp', status='old')

        inquire(file='EVT.out', exist=ex)

        if (ex) then
           open(51, file='EVT.out', status='old', action='read')

        else
           inquire(file='OUTPUT/EVT.out', exist=ex)

           if (ex) then
              open(51, file='OUTPUT/EVT.out', status='old', action='read')

           else
              print*, ' ERROR: EVT.out not found.'
              print*
              stop
           end if
        end if

        print*
        print*, "-----------------------------"
        print*, "        CLR START !"
        print*, "-----------------------------"
        print*

        read(50,*) deltaT
        write(*,'(14X,"deltaT :",I5)') deltaT
        print*

        close(50)

        allocate(sieved(num_Li))
        allocate(hop(total_traj, num_Li))
        allocate(cx(total_traj, num_Li), cy(total_traj, num_Li), cz(total_traj, num_Li))

!---------------------------------------------------------------------------------------

        hop = 0 ; sieved = 0

        read(51,*)  ! skip header

        do i = 1+deltaT, total_traj-deltaT
           if (mod(i,100)==0) write(*,998) " Reading   ", i, " th trajectory"
           read(51,*) dummyi(1), (hop(i,j), j=1,num_Li)
           do j = 1, num_Li
              if (hop(i,j)==1) sieved(j) = 1
           enddo
        enddo

        close(51)

!-Averaged position within the cage residence time--------------------------------------

        do j = 1, num_Li

           cage_start = 1

           do i = 1, total_traj

              if (hop(i,j)==1) then

                 tmp_pos = 0.d0
                 count_cage = 0

                 do k = cage_start, i-1
                    tmp_pos(1) = tmp_pos(1) + x(k,j)
                    tmp_pos(2) = tmp_pos(2) + y(k,j)
                    tmp_pos(3) = tmp_pos(3) + z(k,j)
                    count_cage = count_cage + 1
                 enddo

                 if (count_cage > 0) then
                    meanx = tmp_pos(1) / count_cage
                    meany = tmp_pos(2) / count_cage
                    meanz = tmp_pos(3) / count_cage
                    do k = cage_start, i-1
                       cx(k,j) = meanx
                       cy(k,j) = meany
                       cz(k,j) = meanz
                    enddo
                 endif

                 cx(i,j) = x(i,j)
                 cy(i,j) = y(i,j)
                 cz(i,j) = z(i,j)
                 cage_start = i+1

              endif
           enddo

           ! Last cage segment
           tmp_pos = 0.d0
           count_cage = 0

           do k = cage_start, total_traj
              tmp_pos(1) = tmp_pos(1) + x(k,j)
              tmp_pos(2) = tmp_pos(2) + y(k,j)
              tmp_pos(3) = tmp_pos(3) + z(k,j)
              count_cage = count_cage + 1
           enddo

           if (count_cage > 0) then
              meanx = tmp_pos(1) / count_cage
              meany = tmp_pos(2) / count_cage
              meanz = tmp_pos(3) / count_cage
              do k = cage_start, total_traj
                 cx(k,j) = meanx
                 cy(k,j) = meany
                 cz(k,j) = meanz
              enddo
           endif

        enddo

! writing Trajectory Data---------------------------------------------------------------

        open(100, file='clear.lammpstraj', status='unknown')

        print*
        write(*,'(14X,"Writing clear.lammpstraj ...")')
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
           write(100,'(A)') "ITEM: ATOMS id type x y z sieve"

           do j = 1, nmedia

              if (j <= num_Li) then
                 id = 1
                 write(100,'(I4," ",I1," ",3(G0,1X),I1)') j, id, cx(i,j), cy(i,j), cz(i,j), sieved(j)
              else if (j <= num_Li+num_P) then
                 id = 2
                 write(100,'(I4," ",I1," ",3(G0,1X),I1)') j, id, x(i,j), y(i,j), z(i,j), 0
              else
                 id = 3
                 write(100,'(I4," ",I1," ",3(G0,1X),I1)') j, id, x(i,j), y(i,j), z(i,j), 0
              endif

           enddo

        enddo

        close(100)

        print*
        print*, "   Saved as 'clear.lammpstraj' file  "
        print*, "-----------------------------"
        print*, "          SUCCESS !"
        print*, "-----------------------------"
        print*

!---------------------------------------------------------------------------------------

end program clear_trajectory

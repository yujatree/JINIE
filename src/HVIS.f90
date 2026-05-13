program hop_visualization

! Variables Declaration-----------------------------------------------------------------

        use READER

        implicit none

        integer                         :: deltaT, GAP, nsurr, id
        double precision                :: h_c, norm, dA, dB, dx, dy, dz

        type :: ptr
           double precision, pointer, dimension(:,:) :: p
        end type

        type(ptr), dimension(3)                       :: posA, posB, avgPosA, avgPosB
        double precision, allocatable, dimension(:,:) :: hop
        integer, allocatable, dimension(:,:)          :: event
        integer, allocatable, dimension(:)            :: sieve

! Setting-------------------------------------------------------------------------------

        call dcd_reader

        open(51, file='/mnt/CALC_SSD/jina/CODES/JINIE/inp/hop.inp', status='old')

998	format(A,I5,A)

        read(51,*) deltaT
        read(51,*) h_c
        read(51,*) nsurr

        GAP  = deltaT / 2
        norm = dble(GAP+1)

        print*
        print*, "-----------------------------"
        print*, "        HVIS START !"
        print*, "-----------------------------"
        print*
        write(*,'(14X,"deltaT :",I5)') deltaT
        write(*,'(14X,"h_c    :",F6.3)') h_c
        print*

!---------------------------------------------------------------------------------------

        do i=1,3
           allocate(posA(i)%p(total_traj,num_Li), posB(i)%p(total_traj,num_Li), &
                    avgPosA(i)%p(total_traj,num_Li), avgPosB(i)%p(total_traj,num_Li))
           posA(i)%p = 0.d0 ; posB(i)%p = 0.d0
           avgPosA(i)%p = 0.d0 ; avgPosB(i)%p = 0.d0
        enddo

        allocate(hop(total_traj, num_Li))
        allocate(event(total_traj, num_Li))
        allocate(sieve(num_Li))
        hop = 0.d0 ; event = 0 ; sieve = 0

        ! Averaged position : <r_i(t)>_A,B
        do i = GAP+1, total_traj-GAP
           do j = 1, num_Li
              do k = 0, GAP
                 posA(1)%p(i,j) = posA(1)%p(i,j) + x(i-k,j)
                 posA(2)%p(i,j) = posA(2)%p(i,j) + y(i-k,j)
                 posA(3)%p(i,j) = posA(3)%p(i,j) + z(i-k,j)
                 posB(1)%p(i,j) = posB(1)%p(i,j) + x(i+k,j)
                 posB(2)%p(i,j) = posB(2)%p(i,j) + y(i+k,j)
                 posB(3)%p(i,j) = posB(3)%p(i,j) + z(i+k,j)
              enddo
              do k = 1, 3
                 avgPosA(k)%p(i,j) = posA(k)%p(i,j) / norm
                 avgPosB(k)%p(i,j) = posB(k)%p(i,j) / norm
              enddo
           enddo
        enddo

        ! Hop function
        do i = 1+GAP+GAP, total_traj-GAP-GAP

           if (mod(i,100)==0) write(*,998) " Processing ", i, " th trajectory"

           do j = 1, num_Li

              dA = 0.d0 ; dB = 0.d0

              do k = 0, GAP
                 dx = x(i+k,j) - avgPosA(1)%p(i+k,j)
                 dy = y(i+k,j) - avgPosA(2)%p(i+k,j)
                 dz = z(i+k,j) - avgPosA(3)%p(i+k,j)
                 dB = dB + dx*dx + dy*dy + dz*dz

                 dx = x(i-k,j) - avgPosB(1)%p(i-k,j)
                 dy = y(i-k,j) - avgPosB(2)%p(i-k,j)
                 dz = z(i-k,j) - avgPosB(3)%p(i-k,j)
                 dA = dA + dx*dx + dy*dy + dz*dz
              enddo

              hop(i,j) = sqrt(dA*dB/norm/norm)
              if (hop(i,j) > h_c) event(i,j) = 1

           enddo
        enddo

        ! Sieve : 1 if Li ever hopped during the trajectory
        do j = 1, num_Li
           do i = 1+GAP+GAP, total_traj-GAP-GAP
              if (event(i,j) == 1) then
                 sieve(j) = 1
                 exit
              endif
           enddo
        enddo

!---------------------------------------------------------------------------------------

        open(100, file='hop.lammpstraj', status='unknown')

        print*
        write(*,'(14X,"Writing hop.lammpstraj ...")')
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
           write(100,'(A)') "ITEM: ATOMS id type x y z hop event sieve"

           id = 1
           do j = 1, num_Li
              write(100,'(I4," ",I1," ",3(G0,1X),F12.6,1X,I1,1X,I1)') &
                    j, id, x(i,j), y(i,j), z(i,j), hop(i,j), event(i,j), sieve(j)
           enddo

           do j = num_Li+1, nmedia
              if (j <= num_Li+num_P) then
                 id = 2
              else
                 id = 3
              endif
              write(100,'(I4," ",I1," ",3(G0,1X),F12.6,1X,I1,1X,I1)') &
                    j, id, x(i,j), y(i,j), z(i,j), 0.d0, 0, 0
           enddo

        enddo

        close(51) ; close(100)

        print*
        print*, "   Saved as 'hop.lammpstraj' file  "
        print*, "-----------------------------"
        print*, "          SUCCESS !"
        print*, "-----------------------------"
        print*

!---------------------------------------------------------------------------------------

end program hop_visualization

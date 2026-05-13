program effective_hopping

! Variables Declaration-----------------------------------------------------------------

        use READER

        implicit none

        logical                         :: ex

        integer                         :: deltaT, nsurr, num, typ
        integer                         :: count_Li, count_back, count_forw, count_cage
        double precision                :: h_c, dx, dy, dz, displ_sq
        double precision, parameter     :: cutoff_sq = 1.5d0   ! squared cutoff for same-site (~1.22 A)

        double precision, allocatable   :: cx(:,:), cy(:,:), cz(:,:)
        integer, allocatable            :: sieved(:)
        double precision, allocatable   :: cage_x(:), cage_y(:), cage_z(:)

! Setting-------------------------------------------------------------------------------

        call dcd_reader

998	format(A,I5,A)

        open(50, file='/mnt/CALC_SSD/jina/CODES/JINIE/inp/hop.inp', status='old')

        inquire(file='clear.lammpstraj', exist=ex)

        if (ex) then
           open(52, file='clear.lammpstraj', status='old', action='read')

        else
           inquire(file='OUTPUT/clear.lammpstraj', exist=ex)

           if (ex) then
              open(52, file='OUTPUT/clear.lammpstraj', status='old', action='read')

           else
              print*, ' ERROR: clear.lammpstraj not found.'
              print*
              stop
           end if
        end if

        open(100, file='EFH.out', status='unknown')

        print*
        print*, "-----------------------------"
        print*, "        EFH START !"
        print*, "-----------------------------"
        print*

        read(50,*) deltaT
        read(50,*) h_c
        read(50,*) nsurr
        close(50)

        write(*,'(14X,"deltaT :",I5)') deltaT
        write(*,'(14X,"h_c    :",F6.3)') h_c
        print*

!---------------------------------------------------------------------------------------

        allocate(cx(total_traj, num_Li), cy(total_traj, num_Li), cz(total_traj, num_Li))
        allocate(sieved(num_Li))
        allocate(cage_x(total_traj), cage_y(total_traj), cage_z(total_traj))

        sieved = 0

        do i = 1, total_traj
           if (mod(i,100)==0) write(*,998) " Reading   ", i, " th trajectory"
           do j = 1, 9
              read(52,*)
           enddo
           do j = 1, num_Li
              read(52,*) num, typ, cx(i,j), cy(i,j), cz(i,j), sieved(j)
           enddo
           do j = 1, num_P+num_S
              read(52,*)
           enddo
        enddo

        close(52)

!-Backward hop analysis-----------------------------------------------------------------

        count_back = 0
        count_forw = 0
        count_Li   = 0

        do j = 1, num_Li

           if (sieved(j)==0) cycle

           count_Li = count_Li + 1
           if (mod(count_Li,5)==0) write(*,998) " Analyzing Li", count_Li, " ..."

           ! Build cage position sequence from stable runs in the cleared trajectory
           cage_x(1) = cx(1,j)
           cage_y(1) = cy(1,j)
           cage_z(1) = cz(1,j)
           count_cage = 1

           do i = 2, total_traj

              if (abs(cx(i,j)-cage_x(count_cage)) > 1d-6 .or. &
                  abs(cy(i,j)-cage_y(count_cage)) > 1d-6 .or. &
                  abs(cz(i,j)-cage_z(count_cage)) > 1d-6) then

                 if (i == total_traj) then
                    count_cage = count_cage + 1
                    cage_x(count_cage) = cx(i,j)
                    cage_y(count_cage) = cy(i,j)
                    cage_z(count_cage) = cz(i,j)

                 else if (abs(cx(i+1,j)-cx(i,j)) < 1d-6 .and. &
                          abs(cy(i+1,j)-cy(i,j)) < 1d-6 .and. &
                          abs(cz(i+1,j)-cz(i,j)) < 1d-6) then
                    count_cage = count_cage + 1
                    cage_x(count_cage) = cx(i,j)
                    cage_y(count_cage) = cy(i,j)
                    cage_z(count_cage) = cz(i,j)
                 endif

              endif
           enddo

           ! First hop (cage[1] -> cage[2]) is always forward (no cage[0] to compare)
           if (count_cage >= 2) count_forw = count_forw + 1

           ! Hop n (cage[n-1] -> cage[n]): backward if cage[n] ~ cage[n-2]
           do i = 3, count_cage
              dx = cage_x(i) - cage_x(i-2)
              dy = cage_y(i) - cage_y(i-2)
              dz = cage_z(i) - cage_z(i-2)
              displ_sq = dx*dx + dy*dy + dz*dz
              if (displ_sq < cutoff_sq) then
                 count_back = count_back + 1
              else
                 count_forw = count_forw + 1
              endif
           enddo

        enddo

!-Output--------------------------------------------------------------------------------

        print*
        write(*,'(A,I8)')  "  Forward hops  :", count_forw
        write(*,'(A,I8)')  "  Backward hops :", count_back
        write(*,'(A,I8)')  "  Total hops    :", count_forw+count_back
        print*
        if (count_forw+count_back > 0) then
           write(*,'(A,F8.4)') "  Forward ratio :", dble(count_forw)/dble(count_forw+count_back)
           write(*,'(A,F8.4)') "  Backward ratio:", dble(count_back)/dble(count_forw+count_back)
        endif
        print*

        write(100,'(A,I8)')  " Forward hops  :", count_forw
        write(100,'(A,I8)')  " Backward hops :", count_back
        write(100,'(A,I8)')  " Total hops    :", count_forw+count_back
        if (count_forw+count_back > 0) then
           write(100,'(A,F8.4)') " Forward ratio :", dble(count_forw)/dble(count_forw+count_back)
           write(100,'(A,F8.4)') " Backward ratio:", dble(count_back)/dble(count_forw+count_back)
        endif

        close(100)

        print*
        print*, "   Saved as 'EFH.out' file  "
        print*, "-----------------------------"
        print*, "          SUCCESS !"
        print*, "-----------------------------"
        print*

!---------------------------------------------------------------------------------------

end program effective_hopping

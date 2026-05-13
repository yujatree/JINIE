program event_classification

! Variables Declaration-----------------------------------------------------------------

        use READER

        implicit none

        integer                         :: deltaT, GAP, nsurr, t_lag
        double precision                :: h_c
        double precision, allocatable   :: hop_row(:)
        integer, allocatable            :: evt(:)
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
        open(100, file='EVT.out', status='unknown')

        read(51,*) deltaT
        read(51,*) h_c
        read(51,*) nsurr

        GAP = deltaT / 2

        print*
        print*, "-----------------------------"
        print*, "        EVT START !"
        print*, "-----------------------------"
        print*
        write(*,'(14X,"deltaT :",I5)') deltaT
        write(*,'(14X,"h_c    :",F6.3)') h_c
        print*

!---------------------------------------------------------------------------------------

        allocate(hop_row(num_Li))
        allocate(evt(num_Li))

        write(100,'(A)') "# t_step  event per Li  (0=equil, 1=hop)"

        do i = 1+GAP+GAP, total_traj-GAP-GAP

           if (mod(i,100)==0) write(*,998) " Processing ", i, " th trajectory"

           read(50,*) t_lag, (hop_row(j), j=1,num_Li)

           do j = 1, num_Li
              if (hop_row(j) > h_c) then
                 evt(j) = 1
              else
                 evt(j) = 0
              end if
           enddo

           write(100,'(I8,*(1X,I1))') t_lag, (evt(j), j=1,num_Li)

        enddo

        close(50)

        print*
        print*, "   Saved as 'EVT.out' file  "
        print*, "-----------------------------"
        print*, "          SUCCESS !"
        print*, "-----------------------------"
        print*

        close(51) ; close(100)

!---------------------------------------------------------------------------------------

end program event_classification

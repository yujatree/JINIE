program cluster_fraction

! Variables Declaration-----------------------------------------------------------------

        use READER

        implicit none

        integer            :: clst_idx, S_idx, time_window, max_siz, &
                              total_clst, siz, root_j, root_k
        double precision   :: PS_dist, PP_dist, SS_dist, criterion, &
                              dx, dy, dz, dist

        integer, allocatable :: parent(:), clst_siz(:), counter(:)

! Setting-------------------------------------------------------------------------------

        call dcd_reader

999     format(I6," ",I8," ",F10.6/)
998     format(A,I5,A)

        open(50, file='/mnt/CALC_SSD/jina/CODES/JINIE/inp/cluster.inp', status='old')
        open(100, file='CLF.out', status='unknown')
        write(100,*) "# size count probability"

        print*
        print*, "-----------------------------"
        print*, "         CLF START !"
        print*, "-----------------------------"
        print*

        read(50,*) PS_dist
        read(50,*) PP_dist
        read(50,*) SS_dist
        read(50,*)
        read(50,*) time_window

        print*, " Cutoff distances :"
        write(*,'(14X,"P-S :",F6.3)') PS_dist
        write(*,'(14X,"P-P :",F6.3)') PP_dist
        write(*,'(14X,"S-S :",F6.3)') SS_dist
        print*, " Trajectories to average  :"
        write(*,'(20X,I7)') time_window

        S_idx = num_Li + num_P + 1
        max_siz = num_P + num_S

        allocate(parent(nmedia), clst_siz(nmedia), counter(max_siz))
        counter = 0
        total_clst = 0

!---------------------------------------------------------------------------------------

        do i = total_traj-time_window+1, total_traj

           if (mod(i,10)==0) write(*,998) " Processing ", i, " th trajectory"

           ! Initialize union-find for all non-Li atoms
           do j = num_Li+1, nmedia
              parent(j) = j
           enddo

           ! Build connected components via union-find
           do j = num_Li+1, nmedia-1
              do k = j+1, nmedia

                 if (j < S_idx .and. k < S_idx) criterion = PP_dist
                 if (j < S_idx .and. k >= S_idx) criterion = PS_dist
                 if (j >= S_idx .and. k >= S_idx) criterion = SS_dist

                 call pbc_dist(i, j, k, dx, dy, dz, dist)

                 if (dist < criterion) then
                    root_j = j
                    do while (parent(root_j) /= root_j)
                       root_j = parent(root_j)
                    enddo
                    root_k = k
                    do while (parent(root_k) /= root_k)
                       root_k = parent(root_k)
                    enddo
                    if (root_j /= root_k) then
                       if (root_j < root_k) then
                          parent(root_k) = root_j
                       else
                          parent(root_j) = root_k
                       endif
                    endif
                 endif

              enddo
           enddo

           ! Path compression pass
           do j = num_Li+1, nmedia
              root_j = j
              do while (parent(root_j) /= root_j)
                 root_j = parent(root_j)
              enddo
              parent(j) = root_j
           enddo

           ! Count cluster sizes
           clst_siz = 0
           do j = num_Li+1, nmedia
              clst_siz(parent(j)) = clst_siz(parent(j)) + 1
           enddo

           ! Accumulate histogram over root atoms
           clst_idx = 0
           do j = num_Li+1, nmedia
              if (parent(j) == j) then
                 siz = clst_siz(j)
                 if (siz >= 1 .and. siz <= max_siz) then
                    counter(siz) = counter(siz) + 1
                    total_clst = total_clst + 1
                 endif
                 clst_idx = clst_idx + 1
              endif
           enddo

        enddo

        ! Output time-averaged cluster size distribution
        if (total_clst > 0) then
           do i = 1, max_siz
              if (counter(i) > 0) then
                 write(100,999,advance='no') i, counter(i), dble(counter(i))/dble(total_clst)
              endif
           enddo
        else
           write(100,'(A)') "# No clusters found"
        endif

        print*
        print*, "   Saved as 'CLF.out' file   "
        print*, "-----------------------------"
        print*, "          SUCCESS !"
        print*, "-----------------------------"
        print*

        close(50)
        close(100)

contains

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

end program cluster_fraction

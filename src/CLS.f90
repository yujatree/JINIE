program polyanion_cluster

! Variables Declaration-----------------------------------------------------------------

        use READER

        implicit none

        integer            :: S_idx, time_window, nroot, rnode, aidx
        double precision   :: PS_dist, PP_dist, SS_dist, criterion, &
                              dx, dy, dz, dist

        integer, allocatable :: parent(:), rlabel(:), clst_siz(:), atm_clst(:)

! Setting-------------------------------------------------------------------------------

        call dcd_reader

998     format(A,I5,A)

        open(50, file='/mnt/CALC_SSD/jina/CODES/JINIE/inp/cluster.inp', status='old')
        open(200, file='CLS.lammpstraj', status='unknown')

        print*
        print*, "-----------------------------"
        print*, "         CLS START !"
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
        print*, " Trajectories to process  :"
        write(*,'(20X,I7)') time_window

        S_idx = num_Li + num_P + 1

        allocate(parent(nmedia), rlabel(nmedia), clst_siz(nmedia), atm_clst(nmedia))
        atm_clst = 0
        clst_siz = 0
        nroot = 0

!---------------------------------------------------------------------------------------

        do i = total_traj-time_window+1, total_traj

           if (mod(i,100)==0) write(*,998) " Processing ", i, " th trajectory"

           ! Initialize union-find for non-Li atoms
           parent = 0
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

                 if (dist < criterion) call union_sets(parent, j, k)

              enddo
           enddo

           ! Relabel roots to sequential cluster IDs and count sizes
           rlabel = 0
           nroot = 0
           clst_siz = 0
           atm_clst = 0
           do j = num_Li+1, nmedia
              rnode = find_root(parent, j)
              if (rlabel(rnode) == 0) then
                 nroot = nroot + 1
                 rlabel(rnode) = nroot
              endif
              atm_clst(j) = rlabel(rnode)
              clst_siz(atm_clst(j)) = clst_siz(atm_clst(j)) + 1
           enddo

           ! Write LAMMPS trajectory frame
           write(200,'(A)') "ITEM: TIMESTEP"
           write(200,'(I8)') init_check + (i-1)*dt_check
           write(200,'(A)') "ITEM: NUMBER OF ATOMS"
           write(200,'(I8)') nmedia
           write(200,'(A)') "ITEM: BOX BOUNDS pp pp pp"
           write(200,'(E23.16," ",E23.16)') 0.d0, box(i,1)
           write(200,'(E23.16," ",E23.16)') 0.d0, box(i,3)
           write(200,'(E23.16," ",E23.16)') 0.d0, box(i,6)
           write(200,'(A)') "ITEM: ATOMS id type x y z cluster polytype"
           do j = 1, nmedia
              if (j <= num_Li) then
                 write(200,'(I6," ",I1," ",3(F12.6,1X),I6,I5)') &
                      j, 1, x(i,j), y(i,j), z(i,j), 0, 0
              else
                 if (j <= num_Li+num_P) then
                    aidx = 2
                 else
                    aidx = 3
                 endif
                 write(200,'(I6," ",I1," ",3(F12.6,1X),I6,I5)') &
                      j, aidx, x(i,j), y(i,j), z(i,j), atm_clst(j), clst_siz(atm_clst(j))
              endif
           enddo

        enddo

        print*
        print*, "   Saved as 'CLS.lammpstraj'"
        print*, "-----------------------------"
        print*, "          SUCCESS !"
        print*, "-----------------------------"
        print*

        close(50)
        close(200)

contains

!--------------------------------------------------------------------------------------------

   function find_root(par, idx) result(root)

        integer, intent(inout) :: par(:)
        integer, intent(in)    :: idx
        integer :: root, next, tmp

        root = idx
        do while (par(root) /= root)
           root = par(root)
        enddo

        next = idx
        do while (par(next) /= root)
           tmp = par(next)
           par(next) = root
           next = tmp
        enddo

   end function find_root

!--------------------------------------------------------------------------------------------

   subroutine union_sets(par, a, b)

        integer, intent(inout) :: par(:)
        integer, intent(in)    :: a, b
        integer :: ra, rb

        ra = find_root(par, a)
        rb = find_root(par, b)
        if (ra == rb) return

        if (ra < rb) then
           par(rb) = ra
        else
           par(ra) = rb
        endif

   end subroutine union_sets

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

end program polyanion_cluster

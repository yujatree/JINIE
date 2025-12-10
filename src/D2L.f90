program dcd2lammpstraj


   use READER
   
   call dcd_reader
 
   ! writing Trajectory Data---------------------------------------------------------------

   open(100, file='./out.lammpstraj', status='unknown')
   
   print*
   print*, "-----------------------------"
   print*, "         D2L START !"
   print*, "-----------------------------"
   print*

   do i=1, total_traj
      
      if (mod(i,10)==0) print*, i, 'th traj writing'
      write(100,'(A)') "ITEM: TIMESTEP"
      write(100,'(I8)') dt_check*i
      write(100,'(A)') "ITEM: NUMBER OF ATOMS"
      write(100,'(I8)') nmedia
      write(100,'(A)') "ITEM: BOX BOUNDS pp pp pp"
      write(100,'(E25.16," ",G0)') 0.d0, box(i,1)
      write(100,'(E25.16," ",G0)') 0.d0, box(i,3)
      write(100,'(E25.16," ",G0)') 0.d0, box(i,6)
      !write(100,'(E25.16," ",G0)') -0.5d0*box(i,1), 0.5d0*box(i,1)
      !write(100,'(E25.16," ",G0)') -0.5d0*box(i,3), 0.5d0*box(i,3)
      !write(100,'(E25.16," ",G0)') -0.5d0*box(i,6), 0.5d0*box(i,6)
      write(100,'(A)') "ITEM: ATOMS id type x y z"
   
      do j=1, nmedia
   
         if (j<=num_Li) then
            id=1
         else if (j<=num_Li+num_P) then
            id=2
         else
            id=3
         endif
   
         write(100,'(I4," ",I1," ",3(G0,1X))') j, id, x(i,j), y(i,j), z(i,j)
      enddo
   
   enddo
   
   close(100)
 
   print*
   print*, "-----------------------------"
   print*, "          SUCCESS !"
   print*, "-----------------------------"
   print*
   
   !---------------------------------------------------------------------------------------

end program dcd2lammpstraj

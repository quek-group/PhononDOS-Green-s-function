!------------------------
! Read and write .dyn file
!------------------------


	program main
	implicit none
	INTEGER :: i, j,k, num7,num8
	INTEGER :: n_elmt,n_atm, ibrav, n_atmnew
	REAL*8 :: para(6), num1, num2, num3, num4, num5, num6
	REAL*8 :: latt_vect(3:3)
	REAL*8, allocatable :: mass(:), coord(:,:)
	INTEGER, allocatable :: elmt_index(:), atm_index(:)
	COMPLEX*8, allocatable ::  fc_const(:,:,:,:)
	


! READ INPUT
	open(21, file = 'input.dat', status = 'unknown')
	read (21,*) n_elmt
	read (21,*) n_atm, n_atmnew
	allocate (atm_index( n_atmnew))
	read (21,*) (atm_index(i), i = 1,  n_atmnew)

	close (21)
!READ .dyn
	open(22, file="dyn.dat", status = 'unknown')
	allocate (fc_const(n_atm, n_atm, 3, 3))
	allocate (elmt_index(n_atm))
	allocate (coord(n_atm,3))

	do i = 1,n_atm
	  read (22,*) num7, num8, num1, num2, num3
	  elmt_index(i)= num8
	  coord(i,1) = num1
	  coord(i,2) = num2
	  coord(i,3) = num3
	enddo

	read(22,*)
	read(22,*)
	read(22,*)
	read(22,*)
	read(22,*)
	do i = 1, n_atm
	  do j = 1,n_atm
	    read(22,*)
	    read(22,*) num1, num2, num3, num4, num5, num6
	    fc_const(i, j, 1, 1) = CMPLX( num1, num2)
	    fc_const(i, j, 1, 2) = CMPLX( num3, num4)
	    fc_const(i, j, 1, 3) = CMPLX( num5, num6)
	    read(22,*) num1, num2, num3, num4, num5, num6
	    fc_const(i, j, 2, 1) = CMPLX( num1, num2)
	    fc_const(i, j, 2, 2) = CMPLX( num3, num4)
	    fc_const(i, j, 2, 3) = CMPLX( num5, num6)
	    read(22,*) num1, num2, num3, num4, num5, num6
	    fc_const(i, j, 3, 1) = CMPLX( num1, num2)
	    fc_const(i, j, 3, 2) = CMPLX( num3, num4)
	    fc_const(i, j, 3, 3) = CMPLX( num5, num6)
	  enddo
	enddo

	close(22)
!write FC

	open(23, file = "fc.dat", action = 'WRITE')
	1000 format(i3, 2X,i3)
	2000 format( 3X, 3(f10.8,2X,f10.8,2X))
	3000 format (2X, i3, 2X,i3, 3X,3(f12.10,2X))

	do i = 1,n_atmnew
	  write(23,3000) i, elmt_index(atm_index(i)), (coord(atm_index(i), k), k= 1,3)
	enddo


	do i = 1, n_atmnew
	  do j = 1, n_atmnew
	    write(23,1000) i,j
	    write(23,2000) (fc_const(atm_index(i), atm_index(j), 1, k), k=1,3)
	    write(23,2000) (fc_const(atm_index(i), atm_index(j), 2, k), k=1,3)
	    write(23,2000) (fc_const(atm_index(i), atm_index(j), 3, k), k=1,3)
	  enddo
	enddo

	close(23)
	end program main


	    


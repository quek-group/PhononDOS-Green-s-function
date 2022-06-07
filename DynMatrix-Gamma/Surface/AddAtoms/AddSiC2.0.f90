!------------------------
! Read and print .dyn file in order (atoms from bottom to top)
! add SiC layers at bottom
!------------------------


	program main
	implicit none
	INTEGER :: i, j,k, num7,num8,p,q,atm1,nint,n_cnct,n_diff,n_mt
	INTEGER :: n_elmt,n_atm, prt_n_atm,index_row,index_col,n
	INTEGER :: add_atm, n_bulk,n_left
	REAL*8 :: para(6), num1, num2, num3, num4, num5, num6,sum
	REAL*8 :: latt_vect(3:3), c
	REAL*8, allocatable :: mass(:), coord(:,:),coord_ori(:,:)
	INTEGER, allocatable :: elmt_index(:), atm_index(:),atm_sorted(:),elmt_sorted(:)
	COMPLEX*8, allocatable ::  fc_const(:,:,:,:),fc_const_ord(:,:,:,:),dyn_matrix(:,:),connection_matrix(:,:)
	COMPLEX*8, allocatable ::  Int_fc(:,:,:,:),Int_dyn(:,:),fc_const_prt(:,:,:,:),dyn_matrix_prt(:,:)
	COMPLEX*8, allocatable ::  fc_cnct(:,:,:,:),fc_cnct_d(:,:,:,:),fc_bulk(:,:,:,:),fc_wh(:,:,:,:)
        REAL*8, PARAMETER :: amu_to_ry = 1.660538782E-27/9.10938356E-31/2 
	LOGICAL :: asr


! READ INPUT
	open(21, file = 'input.dat', status = 'unknown')
	read (21,*) n_elmt
	read (21,*) n_atm
	read (21,*) n_mt
	read (21,*) add_atm
	read (21,*) n_cnct  != n_bulk

	close (21)
	prt_n_atm = n_atm + add_atm
	n_bulk = add_atm/12
	n_left = mod(add_atm,12)


!READ dyn.dat
	open(22, file="dyn.dat", status = 'unknown')
	allocate (fc_const(n_atm, n_atm, 3, 3))
	allocate (atm_index(n_atm))
	allocate (elmt_index(n_atm))
	allocate (coord(n_atm,3))
	allocate (coord_ori(n_atm,3))
	allocate (mass(n_elmt))
	
	do i = 1, n_elmt
	  read (22,*) num7, num1
	  mass(i) = num1
	enddo

	do i = 1,n_atm
	  read (22,*) num7, num8, num1, num2, num3
	  atm_index(i) = num7
	  elmt_index(i) = num8
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
!sort the atoms based on it's z coordinates from top to bottom
	allocate (atm_sorted(n_atm))
	allocate (elmt_sorted(n_atm))

	coord_ori = coord
	do i = 1,n_atm
	  atm_sorted(i) = maxloc(coord(:,3),1)
	  elmt_sorted(i) = elmt_index(atm_sorted(i))
	  write(6,*) atm_sorted(i),elmt_index(atm_sorted(i)),coord(atm_sorted(i),3)
	  coord(atm_sorted(i),3)=0
	enddo

!re-order the dynamic matrix


	allocate (fc_const_ord(n_atm, n_atm, 3, 3))
	do i = 1, n_atm
	  do j = 1,n_atm
	    do p = 1,3
	      do q = 1,3
	        fc_const_ord(i,j,p,q)=fc_const(atm_sorted(i),atm_sorted(j),p,q)
	      enddo
	    enddo
	  enddo
	enddo


! read force constant of connnection matrix dagger
	!n_diff = n_atm - n_cnct

	open(22, file="FC_3Nx3N-cnct_d-nint.dat", status = 'unknown')
	allocate (fc_cnct_d(n_cnct, n_cnct, 3, 3))
	fc_cnct_d = complex(0.d0, 0.d0)
	
	do i = 1, n_cnct
	  do j = 1,n_cnct
	    read(22,*)
	    read(22,*) num1, num2, num3, num4, num5, num6
	    fc_cnct_d(i, j, 1, 1) = CMPLX( num1, num2)
	    fc_cnct_d(i, j, 1, 2) = CMPLX( num3, num4)
	    fc_cnct_d(i, j, 1, 3) = CMPLX( num5, num6)
	    read(22,*) num1, num2, num3, num4, num5, num6
	    fc_cnct_d(i, j, 2, 1) = CMPLX( num1, num2)
	    fc_cnct_d(i, j, 2, 2) = CMPLX( num3, num4)
	    fc_cnct_d(i, j, 2, 3) = CMPLX( num5, num6)
	    read(22,*) num1, num2, num3, num4, num5, num6
	    fc_cnct_d(i, j, 3, 1) = CMPLX( num1, num2)
	    fc_cnct_d(i, j, 3, 2) = CMPLX( num3, num4)
	    fc_cnct_d(i, j, 3, 3) = CMPLX( num5, num6)
	  enddo
	enddo

! Define force constant of cnct
	allocate (fc_cnct(n_cnct,n_cnct,  3, 3))
	fc_cnct = complex(0.d0, 0.d0)
	do i = 1,n_cnct
	  do j = 1,n_cnct
	    do p = 1,3
	      do q = 1,3
	        fc_cnct(i, j, p, q) = CONJG(fc_cnct_d(j,i,q,p))
	      enddo
	    enddo
	  enddo
	enddo

! Read force constant of bulk matrix
	open(22, file="FC_3Nx3N-bulk-nint.dat", status = 'unknown')
	allocate (fc_bulk(n_cnct, n_cnct, 3, 3))
	fc_bulk = complex(0.d0, 0.d0)
	
	do i = 1, n_cnct
	  do j = 1,n_cnct
	    read(22,*)
	    read(22,*) num1, num2, num3, num4, num5, num6
	    fc_bulk(i, j, 1, 1) = CMPLX( num1, num2)
	    fc_bulk(i, j, 1, 2) = CMPLX( num3, num4)
	    fc_bulk(i, j, 1, 3) = CMPLX( num5, num6)
	    read(22,*) num1, num2, num3, num4, num5, num6
	    fc_bulk(i, j, 2, 1) = CMPLX( num1, num2)
	    fc_bulk(i, j, 2, 2) = CMPLX( num3, num4)
	    fc_bulk(i, j, 2, 3) = CMPLX( num5, num6)
	    read(22,*) num1, num2, num3, num4, num5, num6
	    fc_bulk(i, j, 3, 1) = CMPLX( num1, num2)
	    fc_bulk(i, j, 3, 2) = CMPLX( num3, num4)
	    fc_bulk(i, j, 3, 3) = CMPLX( num5, num6)
	  enddo
	enddo

! Build the whole force constant matrix
	n = n_atm + n_cnct*(n_bulk+1)
	allocate (fc_wh(n, n, 3, 3))
	fc_wh = complex(0.d0, 0.d0)

!	fc_wh = fc_wh + fc_const_ord	
	do i = 1,n_atm
	  do j = 1,n_atm
	    fc_wh(i,j,:,:) = fc_const_ord(i,j,:,:)
	  enddo
	enddo

	do k = 1, n_bulk+1
	  do i = n_atm+(k-1)*n_cnct+1,n_atm+k*n_cnct
	    do j = n_atm+(k-1)*n_cnct+1,n_atm+k*n_cnct
	      fc_wh(i,j,:,:) = fc_bulk(i-(n_atm+(k-1)*n_cnct), j-(n_atm+(k-1)*n_cnct),:,:)
	    enddo
	  enddo

	  do i = n_atm+(k-1)*n_cnct+1,n_atm+k*n_cnct
	    do j = n_atm+(k-2)*n_cnct+1,n_atm+(k-1)*n_cnct

	      if (j .gt. 0) then

	      fc_wh(i,j,:,:) = fc_cnct(i-(n_atm+(k-1)*n_cnct),j-(n_atm+(k-2)*n_cnct),:,:)
	      endif

	    enddo
	  enddo
	  do i = n_atm+(k-2)*n_cnct+1,n_atm+(k-1)*n_cnct
	    if (i .gt. 0) then
	    do j = n_atm+(k-1)*n_cnct+1,n_atm+k*n_cnct

	      fc_wh(i,j,:,:) = fc_cnct_d(i-(n_atm+(k-2)*n_cnct),j-(n_atm+(k-1)*n_cnct),:,:)
	    enddo
	    endif
	  enddo

	enddo  ! k

! remove interactions between metal and added atoms
	do i = 1, n_mt
	  do j = n_atm+1,n
	    fc_wh(i,j,:,:) =complex(0.d0, 0.d0)
	  enddo
	enddo
	do i = n_atm+1,n
	  do j = 1, n_mt
	    fc_wh(i,j,:,:) =complex(0.d0, 0.d0)
	  enddo
	enddo


! Write prt_n_atm(n_atm+add_atm) dyn.dat
	open(23, file = "fc.dat", action = 'WRITE')
	500 format (4X, i3, 4X, f18.10)
	1000 format(i3, 2X,i3)
	2000 format( 3X, 3(f10.8,2X,f10.8,2X))
	3000 format (2X, i3, 2X,i3, 3X,3(f14.10,2X))
	4000 format (2X,i3, 2X,i3,2X a,f14.10,2X)

	do i = 1, n_elmt
	  write(23,500) i, mass(i)
	enddo

	!write(23,*) "           1     10947.0833589781"
	!write(23,*) "           2     63548.6268941884"
	!write(23,*) "           3     25598.3672624059" 

	do i = 1,n_atm
	  write(23,3000) i, elmt_sorted(i), (coord_ori(atm_sorted(i),k), k= 1,3)
	enddo
	
	k = 1
	do i = n_atm+1, n_atm+add_atm
	  if (mod(k,2) .eq. 0) then
	    j = 1
	  else
	    j = 3
	  endif
	  !write(*,*) j
	  c = 1-REAL(i)/1000
	  write(23,4000) i, j, "0.0000000000  0.0000000000  ",c
	  k = k+1
	enddo


	write(23,*) "  "
	write(23,*) "     Dynamical  Matrix in cartesian axes"
	write(23,*) "  "
	write(23,*) "     q = (    0.000000000   0.000000000   0.000000000 ) "
	write(23,*) "  "

	do i = 1, prt_n_atm
	  do j = 1, prt_n_atm
	    write(23,1000) i,j
	    write(23,2000) (fc_wh(i, j, 1, k), k=1,3)
	    write(23,2000) (fc_wh(i, j, 2, k), k=1,3)
	    write(23,2000) (fc_wh(i, j, 3, k), k=1,3)
	  enddo
	enddo

	close(23)
	end program main




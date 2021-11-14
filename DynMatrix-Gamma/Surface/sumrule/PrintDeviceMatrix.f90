!------------------------
! Read and print .dyn file in order (atoms from bottom to top)
!------------------------


	program main
	implicit none
	INTEGER :: i, j,k, num7,num8,p,q,atm1,nint,n_cnct,n_diff
	INTEGER :: n_elmt,n_atm, prt_n_atm,index_row,index_col
	REAL*8 :: para(6), num1, num2, num3, num4, num5, num6,sum
	REAL*8 :: latt_vect(3:3)
	REAL*8, allocatable :: mass(:), coord(:,:)
	INTEGER, allocatable :: elmt_index(:), atm_index(:),atm_sorted(:),elmt_sorted(:)
	COMPLEX*8, allocatable ::  fc_const(:,:,:,:),fc_const_ord(:,:,:,:),dyn_matrix(:,:),connection_matrix(:,:)
	COMPLEX*8, allocatable ::  Int_fc(:,:,:,:),Int_dyn(:,:),fc_const_prt(:,:,:,:),dyn_matrix_prt(:,:)
	COMPLEX*8, allocatable ::  fc_cnct(:,:,:,:)
        REAL*8, PARAMETER :: amu_to_ry = 1.660538782E-27/9.10938356E-31/2 
	LOGICAL :: asr

! READ INPUT
	open(21, file = 'input.dat', status = 'unknown')
	read (21,*) n_elmt
	read (21,*) n_atm
	read (21,*) prt_n_atm
	read (21,*) asr
	read (21,*) n_cnct

	close (21)


!READ .dyn
	open(22, file="dyn.dat", status = 'unknown')
	allocate (fc_const(n_atm, n_atm, 3, 3))
	allocate (atm_index(n_atm))
	allocate (elmt_index(n_atm))
	allocate (coord(n_atm,3))
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

! read force constant of connnection matrix
	!n_diff = prt_n_atm - n_cnct

	open(22, file="FC_3Nx3N-cnct_d-nint.dat", status = 'unknown')
	allocate (fc_cnct(n_cnct, n_cnct, 3, 3))
	fc_cnct = complex(0.d0, 0.d0)
	

	do i = 1, n_cnct
	  do j = 1,n_cnct
	    read(22,*)
	    read(22,*) num1, num2, num3, num4, num5, num6
	    fc_cnct(i, j, 1, 1) = CMPLX( num1, num2)
	    fc_cnct(i, j, 1, 2) = CMPLX( num3, num4)
	    fc_cnct(i, j, 1, 3) = CMPLX( num5, num6)
	    read(22,*) num1, num2, num3, num4, num5, num6
	    fc_cnct(i, j, 2, 1) = CMPLX( num1, num2)
	    fc_cnct(i, j, 2, 2) = CMPLX( num3, num4)
	    fc_cnct(i, j, 2, 3) = CMPLX( num5, num6)
	    read(22,*) num1, num2, num3, num4, num5, num6
	    fc_cnct(i, j, 3, 1) = CMPLX( num1, num2)
	    fc_cnct(i, j, 3, 2) = CMPLX( num3, num4)
	    fc_cnct(i, j, 3, 3) = CMPLX( num5, num6)
	  enddo
	enddo
close(22)
!sort the atoms based on it's z coordinates from top to bottom
	allocate (atm_sorted(n_atm))
	allocate (elmt_sorted(n_atm))


	do i = 1,n_atm
	  atm_sorted(i) = maxloc(coord(:,3),1)
	  elmt_sorted(i) = elmt_index(atm_sorted(i))
	  write(6,*) atm_sorted(i),elmt_index(atm_sorted(i)),coord(atm_sorted(i),3)
	  coord(atm_sorted(i),3)=0
	enddo

!re-order the dynamic matrix
!Write the whole dynmic matrix

	allocate (fc_const_ord(n_atm, n_atm, 3, 3))
	allocate (dyn_matrix(3*n_atm, 3*n_atm))

	do i = 1, n_atm
	  do j = 1,n_atm
	    do p = 1,3
	      do q = 1,3
	        fc_const_ord(i,j,p,q)=fc_const(atm_sorted(i),atm_sorted(j),p,q)
	        index_row = p+3*(i-1)
	        index_col = q+3*(j-1)
	        dyn_matrix(index_row,index_col) = fc_const_ord(i,j,p,q)/ &
	          	sqrt(mass(elmt_sorted(i))*mass(elmt_sorted(j)))
!write (*,*) mass(elmt_sorted(i))
	      enddo
	    enddo
	  enddo
	enddo
!write(6,*) fc_const_ord(24,24,1,1),fc_const_ord(24,24,1,2)
!write(6,*) force_matrix(4,4),force_matrix(7,1),force_matrix(1,7),force_matrix(71,71)


!Write fc_con vs atoms

	!open(23, file = "dyn_vs_atoms1n.dat", action = 'WRITE')
	!atm1 = 1
	!do i = 1,n_atm
	!  write(23,*) atm1, i, &
	!              ABS(REAL(fc_const_ord(atm1,i,1,1))),ABS(REAL(fc_const_ord(atm1,i,2,2))),ABS(REAL(fc_const_ord(atm1,i,3,3))), &
	!              ABS(REAL(fc_const_ord(atm1,i,1,2))),ABS(REAL(fc_const_ord(atm1,i,2,3))),ABS(REAL(fc_const_ord(atm1,i,3,1)))
	!enddo
	!close(23)



!write the whole dynamic matrix
	!open(23, file = "dyn_3Nx3N.dat", action = 'WRITE')

	!do i = 1,3*n_atm
	!    write(23,*) (REAL(dyn_matrix(i,j)),AIMAG(dyn_matrix(i,j)),j = 1,3*n_atm)  !
	!enddo
	!close(23)



!write the prt_n_atm dynamic matrix
	open(23, file = "dyn_prt.dat", action = 'WRITE')

	do i = 1,3*prt_n_atm
	    write(23,*) (REAL(dyn_matrix(i,j)),AIMAG(dyn_matrix(i,j)),j = 1,3*prt_n_atm)  !
	enddo
	close(23)

	allocate (fc_const_prt(prt_n_atm, prt_n_atm, 3, 3))
	allocate (dyn_matrix_prt(3*prt_n_atm, 3*prt_n_atm))

!Hemitian

	  do p = 1,3
	    do q = 1,3
	      do i = 1,prt_n_atm
	        do j = i+1, 1,prt_n_atm
	          fc_const_ord(i,j,p,q) = (fc_const_ord(i,j,p,q) + CONJG(fc_const_ord(j,i,q,p)))/2
	          fc_const_ord(j,i,q,p) = CONJG(fc_const_ord(i,j,p,q))
	          enddo
	      enddo
	    enddo
	  enddo

! sumrule for prt part

	if (asr) then
	  do p = 1,3
	    do q = 1,3
	      do i = 1,n_atm
	        if (n_atm-i .ge. 4) then
	          sum = 0.d0
	          do j = 1,n_atm
	            sum = sum + fc_const_ord(i,j,p,q)
	            fc_const_prt(i,j,p,q) = fc_const_ord(i,j,p,q)
	          enddo
	          fc_const_prt(i,i,p,q) = fc_const_prt(i,i,p,q) - sum
	        else
	          sum = 0.d0
	          do j = 1,n_atm
	            sum = sum + fc_const_ord(i,j,p,q)
	            fc_const_prt(i,j,p,q) = fc_const_ord(i,j,p,q)
	          enddo
	          do k = 1,4
	             sum = sum + fc_cnct(n_cnct-n_atm+i,k,p,q)
	          enddo
	          fc_const_prt(i,i,p,q) = fc_const_prt(i,i,p,q) - sum
	          endif
	      enddo
	    enddo
	  enddo
	endif

	do i = 1, prt_n_atm
	  do j = 1,prt_n_atm
	    do p = 1,3
	      do q = 1,3
	        index_row = p+3*(i-1)
	        index_col = q+3*(j-1)
	        dyn_matrix_prt(index_row,index_col) = fc_const_prt(i,j,p,q)/ &
	          	sqrt(mass(elmt_sorted(i))*mass(elmt_sorted(j)))
	      enddo
	    enddo
	  enddo
	enddo

	open(23, file = "dyn_prt_sumrule.dat", action = 'WRITE')

	do i = 1,3*prt_n_atm
	    write(23,*) (REAL(dyn_matrix_prt(i,j)),AIMAG(dyn_matrix_prt(i,j)),j = 1,3*prt_n_atm)  !
	enddo
	close(23)





!!! write dynamic matrix of top SiC unitcell. Interactions between nearest nint atoms are considered.
	allocate (Int_fc(n_atm, n_atm, 3, 3))
	allocate (Int_dyn(3*n_atm, 3*n_atm))
	Int_fc(:,:,:,:) = complex(0.d0,0.d0)
	Int_dyn(:,:) = complex(0.d0,0.d0)
	nint = 3 

	do i = 1, n_atm
	  do j = i,i+nint
	  if (j .le. n_atm) then
	    do p = 1,3
	      do q = 1,3
	        Int_fc(i,j,p,q) = fc_const_ord(i,j,p,q)
	        index_row = p+3*(i-1)
	        index_col = q+3*(j-1)
	        Int_dyn(index_row,index_col) = Int_fc(i,j,p,q)/ &
	          	sqrt(mass(elmt_sorted(i))*mass(elmt_sorted(j)))
	      enddo
	    enddo
	    endif
	  enddo
	  do j = i-nint,i-1
	    if (j .gt. 0) then
	    do p = 1,3
	      do q = 1,3
	        Int_fc(i,j,p,q) = fc_const_ord(i,j,p,q)
	        index_row = p+3*(i-1)
	        index_col = q+3*(j-1)
!write(6,*) index_row,index_col
	        Int_dyn(index_row,index_col) = Int_fc(i,j,p,q)/ &
	          	sqrt(mass(elmt_sorted(i))*mass(elmt_sorted(j)))
!write(6,*) Int_dyn(22,1),Int_dyn(index_row,index_col)

	      enddo
	    enddo
	  endif
	  enddo

	enddo
!write(6,*) Int_fc(1,2,1,1),Int_fc(1,3,1,1),Int_fc(1,4,1,1),Int_fc(1,5,1,1)

	open(23, file = "dyn_3Nx3N-device-nint.dat", action = 'WRITE')
	do i = 1,3*prt_n_atm
	    write(23,*) (REAL(Int_dyn(i,j)),AIMAG(Int_dyn(i,j)),j = 1,3*prt_n_atm)  !/2
	enddo
	close(23)





	end program main



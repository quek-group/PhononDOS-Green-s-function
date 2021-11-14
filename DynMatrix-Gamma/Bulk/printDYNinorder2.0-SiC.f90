!------------------------
! Read and print .dyn file in order (atoms from bottom to top)
!------------------------


	program main
	implicit none
	INTEGER :: i, j,k, num7,num8,p,q,atm1,nint,k1,k2,l_sic,start
	INTEGER :: n_elmt,n_atm, ibrav,index_row,index_col
	REAL*8 :: para(6), num1, num2, num3, num4, num5, num6,sum
	REAL*8 :: latt_vect(3:3)
	REAL*8, allocatable :: mass(:), coord(:,:)
	INTEGER, allocatable :: elmt_index(:), atm_index(:),atm_sorted(:),elmt_sorted(:),layer(:,:)
	COMPLEX*8, allocatable ::  fc_const(:,:,:,:),fc_const_ord(:,:,:,:),dyn_matrix(:,:),connection_matrix(:,:)
	COMPLEX*8, allocatable ::  Int_fc(:,:,:,:),Int_dyn(:,:), Int_dyn_wh(:,:),cnct(:,:),cnct_d(:,:)
	COMPLEX*8, allocatable ::  Int_fc_cnct(:,:,:,:), Int_dyn_cnct(:,:),Int_fc_wh(:,:,:,:),Int_fc_cnct_d(:,:,:,:)
	COMPLEX*8, allocatable ::Int_fc_wh_order(:,:,:,:)
        REAL*8, PARAMETER :: amu_to_ry = 1.660538782E-27/9.10938356E-31/2 
	LOGICAL :: asr


! READ INPUT
	open(21, file = 'input.dat', status = 'unknown')
	read (21,*) n_elmt
	read (21,*) n_atm
	read (21,*) asr
	read (21,*) start

	close (21)
	l_sic = n_atm/2

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

!sort the atoms based on it's z coordinates from top to bottom
	allocate (atm_sorted(n_atm))
	allocate (elmt_sorted(n_atm))

	do i = 1,n_atm
	  atm_sorted(i) = maxloc(coord(:,3),1)
	  elmt_sorted(i) = elmt_index(atm_sorted(i))
	  !write(6,*) atm_sorted(i),elmt_index(atm_sorted(i)),coord(atm_sorted(i),3)
	  coord(atm_sorted(i),3)=0
	enddo

! Divide the atms into layers
	allocate (layer(l_sic,2))
	do i = 1,l_sic
	  layer (i,1) = (i-1)*2+1
	  layer (i,2) = (i-1)*2+2
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
	      enddo
	    enddo
	  enddo
	enddo
!write(6,*) fc_const_ord(24,24,1,1),fc_const_ord(24,24,1,2)
!write(6,*) force_matrix(4,4),force_matrix(7,1),force_matrix(1,7),force_matrix(71,71)

!Write fc_con vs atoms

	!open(23, file = "fcDm_vs_atoms1n.dat", action = 'WRITE')
	!atm1 = 1
	!do i = 1,n_atm
	!  write(23,*) atm1, i, &
	!              ABS(REAL(fc_const_ord(atm1,i,1,1))),ABS(REAL(fc_const_ord(atm1,i,2,2))),ABS(REAL(fc_const_ord(atm1,i,3,3))), &
	!              ABS(REAL(fc_const_ord(atm1,i,1,2))),ABS(REAL(fc_const_ord(atm1,i,2,3))),ABS(REAL(fc_const_ord(atm1,i,3,1)))
	!enddo
	!  write(23,*) "2  1",  &
	!  ABS(REAL(fc_const_ord(2,1,1,1))),ABS(REAL(fc_const_ord(2,1,2,2))),ABS(REAL(fc_const_ord(2,1,3,3))), &
	!  ABS(REAL(fc_const_ord(2,1,1,2))),ABS(REAL(fc_const_ord(2,1,2,3))),ABS(REAL(fc_const_ord(2,1,3,1)))

	!close(23)



	
!write the whole dynamic matrix
	!open(23, file = "dyn_3Nx3N.dat", action = 'WRITE')

	!do i = 1,3*n_atm
	!    write(23,*) (REAL(dyn_matrix(i,j)),AIMAG(dyn_matrix(i,j)),j = 1,3*n_atm)  !
	!enddo
	!close(23)

!write top SiC unitcell (bulk)


	!open(23, file = "dyn_3Nx3N-bulk.dat", action = 'WRITE')

	!do i = 1,3*n_atm
	!    write(23,*) (REAL(dyn_matrix(i,j)),AIMAG(dyn_matrix(i,j)),j = 1,3*n_atm)  !
	!enddo
	!close(23)


!write connection matrix

	!open(23, file = "dyn_3Nx3N-connection-12L.dat", action = 'WRITE')
	!allocate (connection_matrix(3*n_atm,3*n_atm))
	!connection_matrix(:,:)=complex(0.d0,0.d0)

	!do i = 1,3*n_atm
	!  do j = i,3*n_atm
	!    connection_matrix(i,j)=dyn_matrix(3*n_atm/2+i,j)
	!  enddo
	!enddo

	!do i = 1,3*n_atm/2
	!    write(23,*) (REAL(connection_matrix(i,j)),AIMAG(connection_matrix(i,j)),j = 1,3*n_atm/2)  !
	!enddo
	!close(23)




!!! write dynamic matrix of top SiC unitcell. Interactions between nearest nint atoms are considered.
	allocate (Int_fc(n_atm, n_atm, 3, 3))
	allocate (Int_dyn(3*n_atm, 3*n_atm))
	Int_fc(:,:,:,:) = complex(0.d0,0.d0)
	Int_dyn(:,:) = complex(0.d0,0.d0)
	nint = 2 

	do i = 1, l_sic
	  do k1 = 1,2
	    do j = i, i+nint
	      if (j .le. l_sic) then
	      do k2 = 1,2
	        do p = 1,3
	        do q = 1,3
	          Int_fc(layer(i,k1),layer(j,k2),p,q) = fc_const_ord(layer(i,k1),layer(j,k2),p,q)
	        enddo
	        enddo
	      enddo
	    endif
	    enddo
	    do j = i- nint, i-1
	if (j .gt. 0) then
	      do k2 = 1,2
	        do p = 1,3
	        do q = 1,3
	          Int_fc(layer(i,k1),layer(j,k2),p,q) = fc_const_ord(layer(i,k1),layer(j,k2),p,q)
	        enddo
	        enddo
	      enddo
	  endif
	    enddo
	  enddo
	enddo
	


!write connection matrix
	allocate (Int_fc_cnct(n_atm, n_atm, 3, 3))
	allocate (Int_fc_cnct_d(n_atm, n_atm, 3, 3))
	Int_fc_cnct = complex(0.d0,0.d0)
	Int_fc_cnct_d = complex(0.d0,0.d0)

	do i = 1,nint
	  do k1 = 1,2
  	    do j = i -nint+l_sic, l_sic
	      do k2 = 1,2
	        do p = 1,3
	        do q = 1,3
	          Int_fc_cnct(layer(i,k1),layer(j,k2),p,q) = fc_const_ord(layer(i,k1),layer(j,k2),p,q)
	          Int_fc_cnct_d(layer(j,k2),layer(i,k1),p,q) = fc_const_ord(layer(j,k2),layer(i,k1),p,q)
	        enddo
	        enddo
	      enddo
	     enddo
	   enddo
	enddo
        

!!!Hermitian and sumrule-simple to get whole dynamic matrix with nint SiC layers interaction
	allocate (Int_fc_wh(n_atm, n_atm, 3, 3))
	Int_fc_wh(:,:,:,:) = Int_fc(:,:,:,:) + Int_fc_cnct(:,:,:,:) + Int_fc_cnct_d(:,:,:,:)

	!Index symmetry Cij = Conjugate(Cji)
	  do p = 1,3
	    do q = 1,3
	      do i = 1,n_atm
	        do j = i+1, n_atm
	          Int_fc_wh(i,j,p,q) = (Int_fc_wh(i,j,p,q) + CONJG(Int_fc_wh(j,i,q,p)))/2
	          Int_fc_wh(j,i,q,p) = CONJG(Int_fc_wh(i,j,p,q))
	          enddo
	      enddo
	    enddo
	  enddo

	if (asr) then
	  do p = 1,3
	    do q = 1,3
	      do i = 1,n_atm
	        sum = 0.d0
	        do j = 1,n_atm
	          sum = sum + Int_fc_wh(i,j,p,q)
	        enddo
	        Int_fc_wh(i,i,p,q) = Int_fc_wh(i,i,p,q) - sum
	      enddo
	    enddo
	  enddo
	endif

!!! Reorder the atoms 

	allocate (Int_fc_wh_order(n_atm, n_atm, 3, 3))
	Int_fc_wh_order = complex(0.d0,0.d0)
	do i = 1,n_atm
	  do j = 1,n_atm
	    index_row = i - (start-1)
	    if (index_row .le. 0) then
	      index_row = index_row+n_atm
	    endif
	    index_col = j - (start-1)
	    if (index_col .le. 0) then
	      index_col = index_col+n_atm
	    endif
	    Int_fc_wh_order(index_row,index_col,:,:) = Int_fc_wh(i,j,:,:)
	  enddo
	 enddo

Int_fc_wh = Int_fc_wh_order

!!! Rewrite nint bulk and connection matrix
	Int_fc(:,:,:,:) = complex(0.d0,0.d0)
	Int_dyn(:,:) = complex(0.d0,0.d0)
	nint = 2 

	do i = 1, l_sic
	  do k1 = 1,2
	    do j = i, i+nint
	      if (j .le. l_sic) then
	      do k2 = 1,2
	        do p = 1,3
	        do q = 1,3
	          Int_fc(layer(i,k1),layer(j,k2),p,q) = Int_fc_wh(layer(i,k1),layer(j,k2),p,q)
	        enddo
	        enddo
	      enddo
	    endif
	    enddo
	    do j = i- nint, i-1
	if (j .gt. 0) then
	      do k2 = 1,2
	        do p = 1,3
	        do q = 1,3
	          Int_fc(layer(i,k1),layer(j,k2),p,q) = Int_fc_wh(layer(i,k1),layer(j,k2),p,q)
	        enddo
	        enddo
	      enddo
	  endif
	    enddo
	  enddo
	enddo

!write bulk force constant 
	open(23, file = "FC_3Nx3N-bulk-nint.dat", action = 'WRITE')
	do i = 1, n_atm
	  do j = 1,n_atm
	    write(23,*) i,j
	    do p = 1,3
	      write(23,*) (REAL(Int_fc(i,j,p,q)),AIMAG(Int_fc(i,j,p,q)),q = 1,3)  
	    enddo
	  enddo
	enddo
	close(23)

	
	do i = 1, n_atm
	  do j = 1,n_atm
	    do p = 1,3
	      do q = 1,3
	        index_row = p+3*(i-1)
	        index_col = q+3*(j-1)
	        Int_dyn(index_row,index_col) = Int_fc(i,j,p,q)/ &
	          	sqrt(mass(elmt_sorted(i))*mass(elmt_sorted(j)))
	      enddo
	    enddo
	  enddo
	enddo

	open(23, file = "dyn_3Nx3N-bulk-nint.dat", action = 'WRITE')
	do i = 1,3*n_atm
	    write(23,*) (REAL(Int_dyn(i,j)),AIMAG(Int_dyn(i,j)),j = 1,3*n_atm) 
	enddo
	close(23)

!bulk finish
!cncnt begin
	Int_fc_cnct = complex(0.d0,0.d0)
	Int_fc_cnct_d = complex(0.d0,0.d0)

	do i = 1,nint
	  do k1 = 1,2
  	    do j = i -nint+l_sic, l_sic
	      do k2 = 1,2
	        do p = 1,3
	        do q = 1,3
	          Int_fc_cnct(layer(i,k1),layer(j,k2),p,q) = Int_fc_wh(layer(i,k1),layer(j,k2),p,q)
	          Int_fc_cnct_d(layer(j,k2),layer(i,k1),p,q) = Int_fc_wh(layer(j,k2),layer(i,k1),p,q)
	        enddo
	        enddo
	      enddo
	     enddo
	   enddo
	enddo

!write cnct force constant 
	open(23, file = "FC_3Nx3N-cnct_d-nint.dat", action = 'WRITE')
	do i = 1, n_atm
	  do j = 1,n_atm
	    write(23,*) i,j
	    do p = 1,3
	      write(23,*) (REAL(Int_fc_cnct_d(i,j,p,q)),AIMAG(Int_fc_cnct_d(i,j,p,q)),q = 1,3)  
	    enddo
	  enddo
	enddo
	close(23)

	allocate (Int_dyn_cnct(3*n_atm, 3*n_atm))
	do i = 1, n_atm
	  do j = 1,n_atm
	    do p = 1,3
	      do q = 1,3
	        index_row = p+3*(i-1)
	        index_col = q+3*(j-1)
	        Int_dyn_cnct(index_row,index_col) = Int_fc_cnct(i,j,p,q)/ &
	          	sqrt(mass(elmt_sorted(i))*mass(elmt_sorted(j)))
	      enddo
	    enddo
	  enddo
	enddo
	open(23, file = "dyn_3Nx3N-connection-nint.dat", action = 'WRITE')

	do i = 1, 3*n_atm
	    write(23,*) (REAL(Int_dyn_cnct(i,j)),AIMAG(Int_dyn_cnct(i,j)),j = 1,3*n_atm)  !
	enddo
	close(23)


!write connection and bulk together
	allocate (Int_dyn_wh(3*n_atm, 3*n_atm))
	allocate (cnct(3*n_atm, 3*n_atm))
	allocate (cnct_d(3*n_atm, 3*n_atm))
	

	Int_dyn_wh = complex(0.d0,0.d0)

	do i = 1,3*n_atm
	  do j = 1, 3*n_atm
	    Int_dyn_wh(i,j) = Int_dyn(i,j)
	    cnct(i,j) = Int_dyn_cnct(i,j)
	  enddo
	enddo
	cnct_d = cnct
	call GetT(cnct_d,3*n_atm)
	do i = 1,3*n_atm
	  do j = 1, 3*n_atm
	    if (cnct(i,j) .ne. complex(0.d0,0.d0)) then
	      Int_dyn_wh(i,j) = cnct(i,j)
	    endif
	    if (cnct_d(i,j) .ne. complex(0.d0,0.d0)) then
	      Int_dyn_wh(i,j) = cnct_d(i,j)
	    endif
	  enddo
	enddo

	
	!open(23, file = "dyn_3Nx3N-whole-nint.dat", action = 'WRITE')

	!do i = 1,3*n_atm
	!    write(23,*) (REAL(Int_dyn_wh(i,j)),AIMAG(Int_dyn_wh(i,j)),j = 1,3*n_atm)  !
	!enddo
	!close(23)
	
	end program main
!-------------------------------------------------------
	subroutine GetT(mtx,n)

	COMPLEX*8 :: mtx(n,n), mtx_T(n,n)
	INTEGER :: n, i, j
	
	do i = 1,n
	  do j = 1,n
	    mtx_T(i,j) = mtx(j,i)
	  enddo
	enddo
	mtx = mtx_T
	end subroutine

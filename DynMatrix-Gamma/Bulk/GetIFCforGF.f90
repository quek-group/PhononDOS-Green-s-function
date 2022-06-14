!-----------
! Read the IFCs after crstal sumrule by QE and analyze.
!  Write scattering, bulk and connection
!-----------


	program main
	implicit none
	INTEGER :: i, j, k, p, q, a,b, s, k1, k2, at
	INTEGER :: nat, nat3, nsic, dimsic, nint, nelmt, l_sic, n, start,atm1
	INTEGER, allocatable :: elmt_index(:), layer(:,:), atm_sorted(:),elmt_sorted(:)
	REAL :: ima
	REAL, allocatable ::  ifcm_sorted(:,:,:,:)
	REAL*8, allocatable :: ifc(:,:), ifcm(:,:,:,:), dynm(:,:,:,:), dyn(:,:), mass(:), dyn_wh(:,:),dyn_wh2(:,:)
	REAL*8, allocatable :: Int_ifc(:,:,:,:), Int_dyn_cnct(:,:,:,:), Int_dyn_bulk(:,:,:,:)
	REAL*8, allocatable :: dyn_cnct(:,:),dyn_bulk(:,:),coord(:,:)
	CHARACTER(len = 200) :: filename 




! read input

	open (21, file = 'input.dat', status = 'unknown')
	read (21, *) nat
	read (21,*) nelmt
	read (21,*) nint
	read (21,*) start

	close(21)
	
	dimsic = 36
	nat3 = nat*3
	l_sic = nat/2
	ima = 0.d0
	n = 36
	at = 12

! Read mass and coord
	open (21, file = 'mass.dat', status = 'unknown')
	allocate (mass(nelmt))
	allocate (elmt_index(nat))
	allocate (coord (nat,3))
	do i = 1, nelmt
	  read (21,*) p, mass(i)
	enddo
	do i = 1,nat
	  read (21,*) p, elmt_index(i), coord(i,1), coord(i,2), coord(i,3)
	enddo


!sort the atoms based on it's z coordinates from top to bottom
	allocate (atm_sorted(nat))
	allocate (elmt_sorted(nat))

	do i = 1,nat
	  atm_sorted(i) = maxloc(coord(:,3),1)
	  elmt_sorted(i) = elmt_index(atm_sorted(i))
	  !write(6,*) atm_sorted(i),elmt_index(atm_sorted(i)),coord(atm_sorted(i),3)
	  coord(atm_sorted(i),3)=0
	enddo
	elmt_index = elmt_sorted

! read dyn
	open (22, file = 'IFC_sumrule.dat', status = 'unknown')
	allocate (ifc(nat3, nat3))
	do i = 1, nat3
	  do j = 1, nat3
	    read(22, *) ifc(i,j)
	  enddo
	enddo



! write ifc into a*b*i*j matrix and sort it
	allocate (ifcm(nat, nat, 3,3))
	do a = 1,nat
	  do b = 1,nat
	    do i = 1,3
	      do j = 1,3
	        ifcm(a,b,i,j)=ifc((a-1)*3+i, (b-1)*3+j)
	      enddo
	     enddo
	    enddo
	enddo
	allocate (ifcm_sorted(nat, nat, 3,3))
	do a = 1,nat
	  do b = 1,nat
	     ifcm_sorted(a,b,:,:)=ifcm(atm_sorted(a), atm_sorted(b), :,:)
	  enddo
	enddo

	ifcm = ifcm_sorted

!Write fc_con vs atoms

	open(28, file = "fcDm_vs_atoms1n.dat", action = 'WRITE')
	atm1 = 1
	do i = 1,nat
	  write(28,*) atm1, i, &
	              ABS(REAL(ifcm(atm1,i,1,1))),ABS(REAL(ifcm(atm1,i,2,2))),ABS(REAL(ifcm(atm1,i,3,3))), &
	              ABS(REAL(ifcm(atm1,i,1,2))),ABS(REAL(ifcm(atm1,i,2,3))),ABS(REAL(ifcm(atm1,i,3,1)))
	enddo

	close(28)


! Divide the atms into layers
	allocate (layer(l_sic,2))
	do i = 1,l_sic
	  layer (i,1) = (i-1)*2+1
	  layer (i,2) = (i-1)*2+2
	enddo


! write out the whole dynm
	allocate (dyn_wh(nat3, nat3))
	dyn_wh = 0.d0
	do a = 1, nat
	  do b = 1, nat
	    do i = 1,3
	      do j = 1,3
	      dyn_wh((a-1)*3+i, (b-1)*3+j) = ifcm(a,b,i,j)/&
	      sqrt (mass(elmt_index(a))*mass(elmt_index(b)))
	      enddo
	     enddo
	  enddo
	enddo

	open (32, file = 'whdyn_bulk.dat', action = 'WRITE')
	do i = 1, nat3
	  do j = 1, nat3
	    write (32, '(4x, e15.8,f15.8,  $)') dyn_wh(i,j), ima
	  enddo
	  write (32, *) '  '
	enddo
	close(32)

! Reorder IFC

	ifcm_sorted = 0.d0
	do a = 1,nat
	  do b = 1,nat
	    p = a - (start -1)
	    if (p .le. 0) then
	      p = p+nat
	    endif
	    q = b-(start -1)
	    if (q .le. 0) then
	      q = q+nat
	    endif
	    ifcm_sorted(p, q, :,:) = ifcm (a, b, :,:)
!write (*,*) p,q,a,b
	  enddo
	enddo
	ifcm = ifcm_sorted



! write ifcm with interactions between nearest nint atoms are considered - for bulk
	allocate (Int_ifc(nat, nat, 3,3))
	Int_ifc = 0.d0
	do i = 1, l_sic
	  do k1 = 1,2
	    do j = i, i+nint
	      if (j .le. l_sic) then
	      do k2 = 1,2
	        do p = 1,3
	        do q = 1,3
	          Int_ifc(layer(i,k1),layer(j,k2),p,q) = ifcm(layer(i,k1),layer(j,k2),p,q)
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
	          Int_ifc(layer(i,k1),layer(j,k2),p,q) = ifcm(layer(i,k1),layer(j,k2),p,q)
	        enddo
	        enddo
	      enddo
	    endif
	    enddo
	  enddo
	enddo


! write IFC
	open(23, file = "FC_3Nx3N-bulk-nint.dat", action = 'WRITE')
	do i = 1, at
	  do j = 1,at
	    write(23,*) i,j
	    do p = 1,3
	      do q = 1,3
	        write(23,'(4x, e15.8,f15.8,  $)') Int_ifc(i,j,p,q), ima
	      enddo
	      write (23, *) '  '
	    enddo
	  enddo
	enddo
	close(23)


	open(23, file = "FC_3Nx3N-cnct_d-nint.dat", action = 'WRITE')
	do i = 1, at
	  do j = at+1, 2*at
	    write(23,*) i,j-at
	    do p = 1,3
	      do q = 1,3
	        write(23,'(4x, e15.8,f15.8,  $)') Int_ifc(i,j,p,q), ima
	      enddo
	      write (23, *) '  '
	    enddo
	  enddo
	enddo
	close(23)


! convert ifc to dyn
	allocate (dynm(nat, nat, 3,3))
	dynm = 0.d0
	do a = 1, nat
	  do b = 1, nat
	    dynm(a,b,:,:) = Int_ifc(a,b,:,:)/&
	    sqrt (mass(elmt_index(a))*mass(elmt_index(b)))
	  enddo
	enddo

! write bulk matrix
	allocate (dyn_bulk(nat3, nat3))
	dyn_bulk = 0.d0
	do a = 1, nat
	  do b = 1, nat
	    do i = 1,3
	      do j = 1,3
	        dyn_bulk((a-1)*3+i, (b-1)*3+j) = dynm(a, b, i, j)
	      enddo
	    enddo
	   enddo
	enddo



! Write bulk region IFC: top_SiC

	open (32, file = 'dyn_3Nx3N-bulk.dat', action = 'WRITE')
	do i = 1, n
	  do j = 1, n
	    write (32, '(4x, e15.8,f15.8,  $)') dyn_bulk(i,j), ima
	  enddo
	  write (32, *) '  '
	enddo
	close(32)



! Write connection region IFC : between top and bottom SiC

	open (33, file = 'dyn_3Nx3N-connection.dat', action = 'WRITE')
	do i = n+1, 2*n
	  do j = 1, n
	    write (33, '(4x, e15.8, f15.8,$)') dyn_bulk(i,j), ima
	  enddo
	  write (33, *) '  '
	enddo
	close(33)

	end program main



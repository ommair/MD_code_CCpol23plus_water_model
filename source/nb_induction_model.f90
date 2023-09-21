module  nb_induction_model

  use, intrinsic :: iso_fortran_env
  use common_arrays
  use utility
  
  implicit none
  private

  public :: nb_induction

  public :: rspace_efield_stat, kspace_efield_stat, self_efield_stat

  public :: rspace_efield_idm, kspace_efield_idm, self_efield_idm

  public :: real_forces, self_forces, kspace_forces

  public :: iteridm

  public :: shortrange_efield_stat_damped,shortrange_efield_idm_damped, shortrange_forces_damped


contains

  subroutine nb_induction(xx,yy,zz,xxs,yys,zzs,boxi,qindpe)

    implicit none

    integer :: is,js,i,j,k,isteps
    real(8) :: xx(nom),yy(nom),zz(nom)
    real(8) :: xxs(nsite,nom),yys(nsite,nom),zzs(nsite,nom)
    real(8) :: boxi,idmlength
    real(8) :: qindpe    
	! DBG
	real(4) iniTime, deltaTime;
	real(4) iniTimeTot, deltaTimeTot;
	! DBG

	! DBG
	iniTimeTot = SECNDS(0.0);
	! DBG

     alphad = alpha

	! DBG
	iniTime = SECNDS(0.0);
	! DBG
    call rspace_efield_stat(xx,yy,zz,xxs,yys,zzs,alphad,boxi)
	! DBG
	deltaTime = SECNDS(iniTime);
!	write(dbgUnit, *) 'P#', MPIrank, '. rspace_efield_stat time	', deltaTime;
	iniTime = SECNDS(0.0);
	! DBG
    call kspace_efield_stat(xx,yy,zz,xxs,yys,zzs,kmax1,kmax2,kmax3,alphad,boxi)
    
	! DBG
	deltaTime = SECNDS(iniTime);
!	write(dbgUnit, *) 'P#', MPIrank, '. kspace_efield_stat time	', deltaTime;
	iniTime = SECNDS(0.0);
	! DBG
    call self_efield_stat(alphad)
	! DBG
	deltaTime = SECNDS(iniTime);
!	write(dbgUnit, *) 'P#', MPIrank, '. self_efield_stat time	', deltaTime;
	! DBG

    if ( induction_type .eq. 'damped' ) then
		! DBG
		iniTime = SECNDS(0.0);
		! DBG
       call shortrange_efield_stat_damped(xx,yy,zz,xxs,yys,zzs,boxi) ! output shtE0x, shtE0y, shtE0z
		! DBG
		deltaTime = SECNDS(iniTime);
!		write(dbgUnit, *) 'P#', MPIrank, '. shortrange_efield_stat_damped time	', deltaTime;
		! DBG
    end if

	! DBG
	iniTime = SECNDS(0.0);
	! DBG
    call iteridm(xx,yy,zz,xxs,yys,zzs,boxi,qindpe)  ! it includes all contribution from real, recip, self 
	! DBG
	deltaTime = SECNDS(iniTime);
!	write(dbgUnit, *) 'P#', MPIrank, '. iteridm time	', deltaTime;
	! DBG

	! DBG
	iniTime = SECNDS(0.0);
	! DBG
    call real_forces(xx,yy,zz,xxs,yys,zzs,boxi,alphad)
	! DBG
	deltaTime = SECNDS(iniTime);
!	write(dbgUnit, *) 'P#', MPIrank, '. real_forces time	', deltaTime;
	iniTime = SECNDS(0.0);
	! DBG
    call kspace_forces(xx,yy,zz,xxs,yys,zzs,kmax1,kmax2,kmax3,alphad,boxi)
	! DBG
	deltaTime = SECNDS(iniTime);
!	write(dbgUnit, *) 'P#', MPIrank, '. kspace_forces time	', deltaTime;
	iniTime = SECNDS(0.0);
	! DBG
    call self_forces(alphad)
	! DBG
	deltaTime = SECNDS(iniTime);
!	write(dbgUnit, *) 'P#', MPIrank, '. self_forces time	', deltaTime;
	iniTime = SECNDS(0.0);
	! DBG

    if ( induction_type .eq. 'damped' ) then
		! DBG
		iniTime = SECNDS(0.0);
		! DBG
       call shortrange_forces_damped(xx,yy,zz,xxs,yys,zzs,boxi) 
		! DBG
		deltaTime = SECNDS(iniTime);
!		write(dbgUnit, *) 'P#', MPIrank, '. shortrange_forces_damped time	', deltaTime;
		! DBG
    end if
	! DBG
	deltaTimeTot = SECNDS(iniTimeTot);
!	write(dbgUnit, *) 'P#', MPIrank, '. nb_induction time	', deltaTimeTot;
	call MPI_BARRIER(MPI_COMM_WORLD, MPIierr);
!	STOP;
	! DBG

  end subroutine nb_induction

  subroutine rspace_efield_stat(xx,yy,zz,xxs,yys,zzs,alp,boxi)
    implicit none
    real(8) :: xx(nom),yy(nom),zz(nom)
    real(8) :: xxs(nsite,nom),yys(nsite,nom),zzs(nsite,nom)
    real(8) :: alp,boxi

    integer :: is,js,i,j,k,isteps
    real(8) :: dist2,rccsq,dist,doti,dotj
    real(8) :: qis,qjs  ! charges on sites
    real(8) :: xcc,ycc,zcc,dx,dy,dz,ddx,ddy,ddz,d1i,d2i,d3i,d5i,d7i
    real(8) :: expar2
	! MPI
	INTEGER(4) MPIchunk_size_loc, iii;
	real(8), allocatable :: MPIbuff1(:), MPIbuff2(:);
	! DBG
	real(4) iniTime, deltaTime;
	! DBG
	! MPI

    do i=1,nm
       do is=1,ns
          rE0x(is,i) = 0.d0 ! initiating permanent charges fields
          rE0y(is,i) = 0.d0
          rE0z(is,i) = 0.d0
       end do
    end do

	! DBG
	iniTime = SECNDS(0.0);
	! DBG

!$omp parallel DEFAULT(SHARED) private(i,id)
!    do i=1,nm       ! loop over molecule i   starts
!	 do i = MPIrank + 1, nm, MPIcomm_size
	if (MPIcomm_size .GT. 1) then
		i = 1;
	else
		i = 0;
	endif
	do while (i .LE. nm)
	  if (MPIrank .EQ. 0) then
		if (MPIcomm_size .GT. 1) then
			call MPI_RECV(MPIaddr, 1, MPI_INTEGER, MPI_ANY_SOURCE, MPI_TAG_CommReady, MPI_COMM_WORLD, MPIstat, MPIierr);
			call MPI_SEND(i, 1, MPI_INTEGER, MPIaddr, MPI_TAG_SendData, MPI_COMM_WORLD, MPIierr);
		endif
		i = i + 1;
	  else
		if (MPIcomm_size .GT. 1) then
			call MPI_SEND(MPIrank, 1, MPI_INTEGER, 0, MPI_TAG_CommReady, MPI_COMM_WORLD, MPIierr);
			call MPI_RECV(i, 1, MPI_INTEGER, 0, MPI_TAG_SendData, MPI_COMM_WORLD, MPIstat, MPIierr);
		endif
	  endif
 	  if ((i .LE. (nm)) .AND. ((MPIrank .GE. 1) .OR. (MPIcomm_size .EQ. 1))) then
		!$omp do private(is,js,iii,j,k,isteps,dist2,rccsq,dist,doti,dotj,qis,qjs,xcc,ycc,zcc,dx)&
		!$omp& private(dy,dz,ddx,ddy,ddz,d1i,d2i,d3i,d5i,d7i,expar2) schedule(dynamic)
        do j=1,nm    ! loop over molecule j   starts

!          if (j.eq.i) goto 11
          if (j .NE. i) then

          dx = xx(i)-xx(j)
          dy = yy(i)-yy(j)
          dz = zz(i)-zz(j)

          xcc = dx - boxi*nint(dx/boxi)
          ycc = dy - boxi*nint(dy/boxi)
          zcc = dz - boxi*nint(dz/boxi)

          rccsq = xcc**2 + ycc**2 + zcc**2

          if (rccsq.lt.rcutsd) then  ! com-com cutoff check starts

             do is=ns-npol+1,ns

                do js=1,ns

                   qis = chgs(is,i)
                   qjs = chgs(js,j)

!                   if (qjs .ne. 0.d0) then

                      ddx = xcc + (xxs(is,i) - xxs(js,j))  ! separation between polarization
                      ddy = ycc + (yys(is,i) - yys(js,j))  ! center of A and charged sites of B
                      ddz = zcc + (zzs(is,i) - zzs(js,j))

                      dist2 = ddx**2 + ddy**2 + ddz**2

!                      if (dist2.lt.rcutsd) then ! cut off on  sites separation
!                      starts

                      dist = sqrt(dist2)

                      expar2 = exp(-(alp*dist)**2.d0)

                      d1i = erfc(alp*dist)/dist
                      d2i = 1.d0/dist2
                      d3i = d2i*(d1i + twosqpi*alp*expar2)
                     
!                      expar2 = exp(-(alp*dist)**2.d0)

                      !rE0x(is,i) = rE0x(is,i) + qjs*ddx*(erfc(alp*dist)*d3i+twosqpi*alp*expar2*d2i)
                      !rE0y(is,i) = rE0y(is,i) + qjs*ddy*(erfc(alp*dist)*d3i+twosqpi*alp*expar2*d2i)
                      !rE0z(is,i) = rE0z(is,i) + qjs*ddz*(erfc(alp*dist)*d3i+twosqpi*alp*expar2*d2i)

                      rE0x(is,i) = rE0x(is,i) + qjs*ddx*d3i
                      rE0y(is,i) = rE0y(is,i) + qjs*ddy*d3i
                      rE0z(is,i) = rE0z(is,i) + qjs*ddz*d3i

!                      end if  ! cut off on  sites separation ends

!                   end if

                end do

             end do

          end if     ! com-com cutoff check ends

          endif

!     11 continue 
        end do        ! loop over molecule j ends
		!$omp end do
	   endif
    end do           ! loop over molecule i ends
	if (MPIrank .EQ. 0) then
		do i = 1, MPIcomm_size - 1
			call MPI_RECV(MPIaddr, 1, MPI_INTEGER, MPI_ANY_SOURCE, MPI_TAG_CommReady, MPI_COMM_WORLD, MPIstat, MPIierr);
			call MPI_SEND(nm + 1, 1, MPI_INTEGER, MPIaddr, MPI_TAG_SendData, MPI_COMM_WORLD, MPIierr);
		enddo
	endif
!$omp end parallel
	! DBG
	deltaTime = SECNDS(iniTime);
!	write(dbgUnit, *) 'P#', MPIrank, '. rspace_efield_stat calc time	', deltaTime;
	! DBG
! MPI sync
	j = nm * ns * 3; !rE0x, rE0y, rE0z
	allocate(MPIbuff1(0:j - 1));
	allocate(MPIbuff2(0:j - 1));
	do k = 0, j - 1
		MPIbuff2(k) = 0.0D0;
	enddo
	k = 0;
	do i = 1, nm
		do is = 1, ns
			MPIbuff1(k) = rE0x(is, i); k = k + 1;
			MPIbuff1(k) = rE0y(is, i); k = k + 1;
			MPIbuff1(k) = rE0z(is, i); k = k + 1;
		enddo
	enddo
	call MPI_ALLREDUCE(MPIbuff1, MPIbuff2, j, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, MPIierr);
	k = 0;
	do i = 1, nm
		do is = 1, ns
			rE0x(is, i) = MPIbuff2(k); k = k + 1;
			rE0y(is, i) = MPIbuff2(k); k = k + 1;
			rE0z(is, i) = MPIbuff2(k); k = k + 1;
		enddo
	enddo
	deallocate(MPIbuff1);
	deallocate(MPIbuff2);
! MPI sync

	return;
  end subroutine rspace_efield_stat

  subroutine kspace_efield_stat2(xx,yy,zz,xxs,yys,zzs,lmax,mmax,nmax,alp,boxi)
    implicit none
    real(8) :: xx(nom),yy(nom),zz(nom)
    real(8) :: xxs(nsite,nom),yys(nsite,nom),zzs(nsite,nom)
    integer :: lmax,mmax,nmax 
    real(8) :: alp,boxi

    real(8) :: qkcpe,qvir_kc
    real(8) :: volmi 
    integer :: i,j,is,js,ii,k
    integer, parameter :: nqmax = nsite*nom  
    integer :: klmmax,klmmaxsq 
!    complex(8) :: el(nqmax,0:20)
!    complex(8) :: em(nqmax,-20:20),en(nqmax,-20:20)
 !   real(8) :: qi(nqmax),xi(nqmax),yi(nqmax),zi(nqmax)
    complex(8), allocatable :: expikr(:)
    complex(8), allocatable :: el(:, :)
    complex(8), allocatable :: em(:, :), en(:, :)
    real(8), allocatable :: qi(:), xi(:), yi(:), zi(:)
    real(8) :: qpe,rl,rm,rn,rksq,qforce,qefld
    complex(8) :: sumqex
    real(8), allocatable :: lkE0x(:, :), lkE0y(:, :), lkE0z(:, :)

    real(8), parameter :: twopi = 2.d0*pi
    integer :: noq
!    real(8) :: ak(ksqmax)
    real(8), allocatable :: ak(:)
    integer :: llim,mlim,nlim,mmin,nmin

    integer :: ll,mm,nn,kk,lmnp
	! MPI
	INTEGER(4) MPIchunk_size_loc, iii;
	real(8), allocatable :: MPIbuff1(:), MPIbuff2(:);
	! DBG
	real(4) iniTime, deltaTime;
	! DBG
	! MPI

    do i=1,nm
       do is=1,ns
          kE0x(is,i) = 0.d0 ! initiating permanent charges fields
          kE0y(is,i) = 0.d0
          kE0z(is,i) = 0.d0
       end do
    end do 

    allocate(el(nqmax,0:20))
    allocate(em(nqmax,-20:20))
    allocate(en(nqmax,-20:20))
    allocate(qi(nqmax))
    allocate(xi(nqmax))
    allocate(yi(nqmax))
    allocate(zi(nqmax))

!$omp parallel DEFAULT(SHARED) private(id,i,iii,j,is,js,ii,ll,mm,nn,kk,lmnp)&
!$omp& private(klmmax,klmmaxsq,llim,mlim,nlim,volmi,mmin,nmin,expikr,ak)&
!$omp& private(lkE0x,lkE0y,lkE0z)

    allocate(expikr(nqmax))
    allocate(ak(ksqmax))
    allocate(lkE0x(ns,nm))
    allocate(lkE0y(ns,nm))
    allocate(lkE0z(ns,nm))

    klmmax = min(dble(lmax),dble(mmax),dble(nmax))
    klmmaxsq = klmmax**2

    llim = lmax 
    mlim = mmax
    nlim = nmax

    volmi = boxi**3
    
    do i=1,nm
       do is=1,ns
          lkE0x(is,i) = 0.d0 ! initiating permanent charges fields
          lkE0y(is,i) = 0.d0
          lkE0z(is,i) = 0.d0
       end do
    end do 
!    collect quantities for charged sites

!$OMP MASTER
      i = 0
      do j = 1,nm
         do js = 1,ns
!            if (chgs(js,j).ne.0.0) then
               i = i+1
               qi(i) = chgs(js,j)
               xi(i) = xx(j) + xxs(js,j)
               yi(i) = yy(j) + yys(js,j)
               zi(i) = zz(j) + zzs(js,j)
!            endif
         enddo
      enddo
      noq = i
!$omp end MASTER
!$omp barrier

!$omp do
    do i = 1,noq
       el(i,0) = (1.0,0.0)
       em(i,0) = (1.0,0.0)
       en(i,0) = (1.0,0.0)
       el(i,1) = exp((0.,1.)*twopi*xi(i)/boxi)
       em(i,1) = exp((0.,1.)*twopi*yi(i)/boxi)
       en(i,1) = exp((0.,1.)*twopi*zi(i)/boxi)
    enddo
!$omp end do
!$omp barrier
 
!$omp MASTER
    do ll = 2,llim
         do i = 1,noq
            el(i,ll) = el(i,ll-1)*el(i,1)
         enddo
    enddo
    do mm = 2,mlim
        do i = 1,noq
            em(i,mm) = em(i,mm-1)*em(i,1)
        enddo
    enddo
    do nn = 2,nlim
        do i = 1,noq
            en(i,nn) = en(i,nn-1)*en(i,1)
        enddo
    enddo
   do mm = -mlim,-1
        do i = 1,noq
            em(i,mm) = conjg(em(i,-mm))
        enddo
    enddo
    do nn = -nlim,-1
        do i = 1,noq
            en(i,nn) = conjg(en(i,-nn))
        enddo
    enddo   
!$omp end MASTER
!$omp barrier

    mmin = 0
    nmin = 1

	! DBG
	iniTime = SECNDS(0.0);
	! DBG

!    do ll=0,llim
!	 do ll = MPIrank, llim, MPIcomm_size
	if (MPIcomm_size .GT. 1) then
		ll = 0;
	else
		ll = -1;
	endif
	do while (ll .LE. llim)
	  if (MPIrank .EQ. 0) then
		if (MPIcomm_size .GT. 1) then
			call MPI_RECV(MPIaddr, 1, MPI_INTEGER, MPI_ANY_SOURCE, MPI_TAG_CommReady, MPI_COMM_WORLD, MPIstat, MPIierr);
			call MPI_SEND(ll, 1, MPI_INTEGER, MPIaddr, MPI_TAG_SendData, MPI_COMM_WORLD, MPIierr);
		endif
		ll = ll + 1;
	  else
		if (MPIcomm_size .GT. 1) then
			call MPI_SEND(MPIrank, 1, MPI_INTEGER, 0, MPI_TAG_CommReady, MPI_COMM_WORLD, MPIierr);
			call MPI_RECV(ll, 1, MPI_INTEGER, 0, MPI_TAG_SendData, MPI_COMM_WORLD, MPIstat, MPIierr);
		endif
	  endif
 	  if ((ll .LE. (llim)) .AND. ((MPIrank .GE. 1) .OR. (MPIcomm_size .EQ. 1))) then
       
       rl = twopi*ll/boxi
   
       if (ll .GT. 0) mmin = -mlim
 
	   !$omp do private(rl,rm,rn,rksq,qefld,sumqex) schedule(dynamic)
       do mm=mmin,mlim

          rm = twopi*mm/boxi

          do nn=nmin,nlim
 
             rn = twopi*nn/boxi

             kk=ll**2+mm**2+nn**2     
             lmnp = ll+mm+nn 

!             if ((kk.le.klmmaxsq) .and. .not.(mod(lmnp,2).ne.0) ) then
             if ((kk.le.klmmaxsq)) then

                rksq = rl*rl + rm*rm + rn*rn

                ak(kk) = twopi/volmi*exp(-0.25d0*rksq/alp**2)/rksq

                ! form expikr for each charge
                do i = 1,noq
                   expikr(i) = el(i,ll)*em(i,mm)*en(i,nn)
                enddo      

             !form sum of qi*expikr

                sumqex = (0.0,0.0)
                do i = 1,noq
                   sumqex = sumqex+qi(i)*expikr(i)

!                    sumqex = sumqex+cmplx(qi(i),0.d0)*expikr(i)
                enddo

               ii=0
               do i=1,nm
                  do is=1,ns

                       ii=ii+1

                       qefld = -4.0*ak(kk)*aimag(sumqex*conjg(expikr(ii)))

            !            qefld = -4.0*ak(kk)*aimag(conjg(sumqex)*expikr(ii))                        
                      
                        lkE0x(is,i) = lkE0x(is,i) + rl*qefld
                        lkE0y(is,i) = lkE0y(is,i) + rm*qefld
                        lkE0z(is,i) = lkE0z(is,i) + rn*qefld
                
                  end do 
               end do

            end if

          end do
          nmin = -nlim

       end do
!       mmin = -mlim
	   !$omp end do
	  endif
    end do
	if (MPIrank .EQ. 0) then
		do i = 1, MPIcomm_size - 1
			call MPI_RECV(MPIaddr, 1, MPI_INTEGER, MPI_ANY_SOURCE, MPI_TAG_CommReady, MPI_COMM_WORLD, MPIstat, MPIierr);
			call MPI_SEND(llim + 1, 1, MPI_INTEGER, MPIaddr, MPI_TAG_SendData, MPI_COMM_WORLD, MPIierr);
		enddo
	endif
	! DBG
	deltaTime = SECNDS(iniTime);
!	write(dbgUnit, *) 'P#', MPIrank, '. kspace_efield_stat calc time	', deltaTime;
	! DBG

!$omp critical
    do i=1,nm
       do is=1,ns
          kE0x(is,i) = kE0x(is,i) + lkE0x(is,i)
          kE0y(is,i) = kE0y(is,i) + lkE0y(is,i)
          kE0z(is,i) = kE0z(is,i) + lkE0z(is,i)
       end do
    end do 

!$omp end critical

    deallocate(expikr)
    deallocate(ak)
    deallocate(lkE0x)
    deallocate(lkE0y)
    deallocate(lkE0z)

!$omp end parallel

    deallocate(el)
    deallocate(em)
    deallocate(en)
    deallocate(qi)
    deallocate(xi)
    deallocate(yi)
    deallocate(zi)
! MPI sync
	j = nm * ns * 3; !kE0x, kE0y, kE0z
	allocate(MPIbuff1(0:j - 1));
	allocate(MPIbuff2(0:j - 1));
	do k = 0, j - 1
		MPIbuff2(k) = 0.0D0;
	enddo
	k = 0;
	do i = 1, nm
		do is = 1, ns
			MPIbuff1(k) = kE0x(is, i); k = k + 1;
			MPIbuff1(k) = kE0y(is, i); k = k + 1;
			MPIbuff1(k) = kE0z(is, i); k = k + 1;
		enddo
	enddo
	call MPI_ALLREDUCE(MPIbuff1, MPIbuff2, j, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, MPIierr);
	k = 0;
	do i = 1, nm
		do is = 1, ns
			kE0x(is, i) = MPIbuff2(k); k = k + 1;
			kE0y(is, i) = MPIbuff2(k); k = k + 1;
			kE0z(is, i) = MPIbuff2(k); k = k + 1;
		enddo
	enddo
	deallocate(MPIbuff1);
	deallocate(MPIbuff2);
! MPI sync

	return;
  end subroutine kspace_efield_stat2


  subroutine kspace_efield_stat(xx,yy,zz,xxs,yys,zzs,lmax,mmax,nmax,alp,boxi)
    implicit none
    real(8) :: xx(nom),yy(nom),zz(nom)
    real(8) :: xxs(nsite,nom),yys(nsite,nom),zzs(nsite,nom)
    integer :: lmax,mmax,nmax 
    real(8) :: alp,boxi

    real(8) :: qkcpe,qvir_kc
    real(8) :: volmi 
    integer :: i,j,is,js,ii,k,li,mi,ni,mx
    integer, parameter :: nqmax = nsite*nom  
    integer :: klmmax,klmmaxsq 
!    complex(8) :: el(nqmax,0:20)
!    complex(8) :: em(nqmax,-20:20),en(nqmax,-20:20)
 !   real(8) :: qi(nqmax),xi(nqmax),yi(nqmax),zi(nqmax)
    complex(8), allocatable :: expikr(:)
    complex(8), allocatable :: el(:, :)
    complex(8), allocatable :: em(:, :), en(:, :)
    complex(8), allocatable :: elc(:, :),emc(:, :),enc(:, :)
    complex(8), allocatable :: els(:, :),ems(:, :),ens(:, :)
    real(8), allocatable :: qi(:), xi(:), yi(:), zi(:)
    real(8), allocatable :: ckr(:),skr(:),clm(:),slm(:),ckc(:),cks(:)
    real(8) :: ckcs,ckss
    real(8) :: qpe,rl,rm,rn,rksq,qforce,qefld
    complex(8) :: sumqex
    real(8), allocatable :: lkE0x(:, :), lkE0y(:, :), lkE0z(:, :)

    real(8), parameter :: twopi = 2.d0*pi
    integer :: noq
!    real(8) :: ak(ksqmax)
    real(8), allocatable :: ak(:)
    integer :: llim,mlim,nlim,mmin,nmin

    integer :: ll,mm,nn,kk,lmnp
        ! MPI
        INTEGER(4) MPIchunk_size_loc, iii;
        real(8), allocatable :: MPIbuff1(:), MPIbuff2(:);
        ! DBG
        real(4) iniTime, deltaTime;
        ! DBG
        ! MPI

    do i=1,nm
       do is=1,ns
          kE0x(is,i) = 0.d0 ! initiating permanent charges fields
          kE0y(is,i) = 0.d0
          kE0z(is,i) = 0.d0
       end do
    end do 

    allocate(el(nqmax,0:20))
    allocate(em(nqmax,-20:20))
    allocate(en(nqmax,-20:20))
    allocate(qi(nqmax))
    allocate(xi(nqmax))
    allocate(yi(nqmax))
    allocate(zi(nqmax))
    allocate(elc(nqmax,0:20))
    allocate(emc(nqmax,-20:20))
    allocate(enc(nqmax,-20:20))
    allocate(els(nqmax,0:20))
    allocate(ems(nqmax,-20:20))
    allocate(ens(nqmax,-20:20))
    allocate(ckr(nqmax))
    allocate(skr(nqmax))
    allocate(clm(nqmax))
    allocate(slm(nqmax))
    allocate(ckc(nqmax))
    allocate(cks(nqmax))

!$omp parallel DEFAULT(SHARED) private(id,i,iii,j,is,js,ii,ll,mm,nn,kk,lmnp)&
!$omp& private(klmmax,klmmaxsq,llim,mlim,nlim,volmi,mmin,nmin,expikr,ak)&
!$omp& private(lkE0x,lkE0y,lkE0z)

    allocate(expikr(nqmax))
    allocate(ak(ksqmax))
    allocate(lkE0x(ns,nm))
    allocate(lkE0y(ns,nm))
    allocate(lkE0z(ns,nm))

    klmmax = min(dble(lmax),dble(mmax),dble(nmax))
    klmmaxsq = klmmax**2

    llim = lmax 
    mlim = mmax
    nlim = nmax

    volmi = boxi**3
    
    do i=1,nm
       do is=1,ns
          lkE0x(is,i) = 0.d0 ! initiating permanent charges fields
          lkE0y(is,i) = 0.d0
          lkE0z(is,i) = 0.d0
       end do
    end do 
!    collect quantities for charged sites

!$OMP MASTER
      i = 0
      do j = 1,nm
         do js = 1,ns
!            if (chgs(js,j).ne.0.0) then
               i = i+1
               qi(i) = chgs(js,j)
               xi(i) = xx(j) + xxs(js,j)
               yi(i) = yy(j) + yys(js,j)
               zi(i) = zz(j) + zzs(js,j)
!            endif
         enddo
      enddo
      noq = i
!$omp end MASTER
!$omp barrier

!$omp do
    do ii = 1,noq
       elc(ii,0) = 1.d0
       emc(ii,0) = 1.d0
       enc(ii,0) = 1.d0

       els(ii,0) = 0.d0
       ems(ii,0) = 0.d0
       ens(ii,0) = 0.d0

       elc(ii,1) = cos(twopi*xi(ii)/boxi)
       emc(ii,1) = cos(twopi*yi(ii)/boxi)
       enc(ii,1) = cos(twopi*zi(ii)/boxi)
 
       els(ii,1) = sin(twopi*xi(ii)/boxi)
       ems(ii,1) = sin(twopi*yi(ii)/boxi)
       ens(ii,1) = sin(twopi*zi(ii)/boxi)
    enddo
!$omp end do
!$omp barrier
 
!$omp MASTER

    do ll=2,llim
       do ii=1,noq

          elc(ii,ll) = elc(ii,ll-1)*elc(ii,1) - els(ii,ll-1)*els(ii,1)
          emc(ii,ll) = emc(ii,ll-1)*emc(ii,1) - ems(ii,ll-1)*ems(ii,1)
          enc(ii,ll) = enc(ii,ll-1)*enc(ii,1) - ens(ii,ll-1)*ens(ii,1)

          els(ii,ll) = els(ii,ll-1)*elc(ii,1) + elc(ii,ll-1)*els(ii,1)
          ems(ii,ll) = ems(ii,ll-1)*emc(ii,1) + emc(ii,ll-1)*ems(ii,1)
          ens(ii,ll) = ens(ii,ll-1)*enc(ii,1) + enc(ii,ll-1)*ens(ii,1)

       end do
    end do   
  
!$omp end MASTER
!$omp barrier

    mmin = 0
    nmin = 1

        ! DBG
        iniTime = SECNDS(0.0);
        ! DBG

!    do ll=0,llim
!        do ll = MPIrank, llim, MPIcomm_size
        if (MPIcomm_size .GT. 1) then
                ll = 0;
        else
                ll = -1;
        endif
        do while (ll .LE. llim)
          if (MPIrank .EQ. 0) then
                if (MPIcomm_size .GT. 1) then
                        call MPI_RECV(MPIaddr, 1, MPI_INTEGER, MPI_ANY_SOURCE,MPI_TAG_CommReady, MPI_COMM_WORLD, MPIstat, MPIierr);
                        call MPI_SEND(ll, 1, MPI_INTEGER, MPIaddr,MPI_TAG_SendData, MPI_COMM_WORLD, MPIierr);
                endif
                ll = ll + 1;
          else
                if (MPIcomm_size .GT. 1) then
                        call MPI_SEND(MPIrank, 1, MPI_INTEGER, 0,MPI_TAG_CommReady, MPI_COMM_WORLD, MPIierr);
                        call MPI_RECV(ll, 1, MPI_INTEGER, 0, MPI_TAG_SendData,MPI_COMM_WORLD, MPIstat, MPIierr);
                endif
          endif
          if ((ll .LE. (llim)) .AND. ((MPIrank .GE. 1) .OR. (MPIcomm_size .EQ.1))) then
       
       rl = twopi*ll/boxi
  
       li = ll
 
       if (ll .GT. 0) mmin = -mlim
 
           !$omp do private(rl,rm,rn,rksq,qefld,sumqex) schedule(dynamic)
       do mm=mmin,mlim

          rm = twopi*mm/boxi

          mi = abs(mm)

          if (mm.ge.0) then
            do ii=1,noq

               clm(ii)=elc(ii,li)*emc(ii,mi)-els(ii,li)*ems(ii,mi)
               slm(ii)=els(ii,li)*emc(ii,mi)+ems(ii,mi)*elc(ii,li)

            end do
          else
            do ii=1,noq

               clm(ii)=elc(ii,li)*emc(ii,mi)+els(ii,li)*ems(ii,mi)
               slm(ii)=els(ii,li)*emc(ii,mi)-ems(ii,mi)*elc(ii,li)

            end do
          end if


          do nn=nmin,nlim

             ni = abs(nn)
 
             rn = twopi*nn/boxi

             kk=ll**2+mm**2+nn**2     
             lmnp = ll+mm+nn 

!             if ((kk.le.klmmaxsq) .and. .not.(mod(lmnp,2).ne.0) ) then
             if ((kk.le.klmmaxsq)) then

                rksq = rl*rl + rm*rm + rn*rn

                ak(kk) = twopi/volmi*exp(-0.25d0*rksq/alp**2)/rksq

               if (nn.ge.0) then

                   do ii=1,noq

                      !ckc(ii)=qi(ii)*(clm(ii)*enc(ii,ni)-slm(ii)*ens(ii,ni))
                      !cks(ii)=qi(ii)*(slm(ii)*enc(ii,ni)+clm(ii)*ens(ii,ni))

                     ckr(ii)=clm(ii)*enc(ii,ni)-slm(ii)*ens(ii,ni)
                     skr(ii)=slm(ii)*enc(ii,ni)+clm(ii)*ens(ii,ni)

                     ckc(ii)=ckr(ii)
                     cks(ii)=skr(ii)

                   end do

               else

                   do ii=1,noq

                    !  ckc(ii)=qi(ii)*(clm(ii)*enc(ii,ni)+slm(ii)*ens(ii,ni))
                    !  cks(ii)=qi(ii)*(slm(ii)*enc(ii,ni)-clm(ii)*ens(ii,ni))

                      ckr(ii)=clm(ii)*enc(ii,ni)+slm(ii)*ens(ii,ni)
                      skr(ii)=slm(ii)*enc(ii,ni)-clm(ii)*ens(ii,ni)

                      ckc(ii)=ckr(ii)
                      cks(ii)=skr(ii)

                   end do

               end if

              ckcs = 0.d0
              ckss = 0.d0
 
               do mx=1,noq

                  ckcs=ckcs+qi(mx)*ckc(mx)
                  ckss=ckss+qi(mx)*cks(mx)

               end do

               ii=0
               do i=1,nm
                  do is=1,ns

                       ii=ii+1

                       qefld = 4.0*ak(kk)*(cks(ii)*ckcs-ckc(ii)*ckss)
 
!                        qefld = -4.0*ak(kk)*aimag(conjg(sumqex)*expikr(ii))
                       
                        lkE0x(is,i) = lkE0x(is,i) + rl*qefld
                        lkE0y(is,i) = lkE0y(is,i) + rm*qefld
                        lkE0z(is,i) = lkE0z(is,i) + rn*qefld
                
                  end do 
               end do

            end if

          end do
          nmin = -nlim

       end do
!       mmin = -mlim
           !$omp end do
          endif
    end do
        if (MPIrank .EQ. 0) then
                do i = 1, MPIcomm_size - 1
                        call MPI_RECV(MPIaddr, 1, MPI_INTEGER, MPI_ANY_SOURCE,MPI_TAG_CommReady, MPI_COMM_WORLD, MPIstat, MPIierr);
                        call MPI_SEND(llim + 1, 1, MPI_INTEGER, MPIaddr,MPI_TAG_SendData, MPI_COMM_WORLD, MPIierr);
                enddo
        endif
        ! DBG
        deltaTime = SECNDS(iniTime);
!       write(dbgUnit, *) 'P#', MPIrank, '. kspace_efield_stat calc time
!       ', deltaTime;
        ! DBG

!$omp critical
    do i=1,nm
       do is=1,ns
          kE0x(is,i) = kE0x(is,i) + lkE0x(is,i)
          kE0y(is,i) = kE0y(is,i) + lkE0y(is,i)
          kE0z(is,i) = kE0z(is,i) + lkE0z(is,i)
       end do
    end do 

!$omp end critical

    deallocate(expikr)
    deallocate(ak)
    deallocate(lkE0x)
    deallocate(lkE0y)
    deallocate(lkE0z)

!$omp end parallel

    deallocate(el)
    deallocate(em)
    deallocate(en)
    deallocate(qi)
    deallocate(xi)
    deallocate(yi)
    deallocate(zi)
    deallocate(elc)
    deallocate(emc)
    deallocate(enc)
    deallocate(els)
    deallocate(ems)
    deallocate(ens)
    deallocate(ckr)
    deallocate(skr)
    deallocate(clm)
    deallocate(slm)
    deallocate(ckc)
    deallocate(cks)
! MPI sync
        j = nm * ns * 3; !kE0x, kE0y, kE0z
        allocate(MPIbuff1(0:j - 1));
        allocate(MPIbuff2(0:j - 1));
        do k = 0, j - 1
                MPIbuff2(k) = 0.0D0;
        enddo
        k = 0;
        do i = 1, nm
                do is = 1, ns
                        MPIbuff1(k) = kE0x(is, i); k = k + 1;
                        MPIbuff1(k) = kE0y(is, i); k = k + 1;
                        MPIbuff1(k) = kE0z(is, i); k = k + 1;
                enddo
        enddo
        call MPI_ALLREDUCE(MPIbuff1, MPIbuff2, j, MPI_DOUBLE_PRECISION, MPI_SUM,MPI_COMM_WORLD, MPIierr);
        k = 0;
        do i = 1, nm
                do is = 1, ns
                        kE0x(is, i) = MPIbuff2(k); k = k + 1;
                        kE0y(is, i) = MPIbuff2(k); k = k + 1;
                        kE0z(is, i) = MPIbuff2(k); k = k + 1;
                enddo
        enddo
        deallocate(MPIbuff1);
        deallocate(MPIbuff2);
! MPI sync

        return;
  end subroutine kspace_efield_stat

  subroutine self_efield_stat(alp)
    implicit none
    real(8) :: alp

    integer :: is,js,i,j,k
    real(8) :: expar2 
    real(8) :: drr,dist2
    real(8) :: qis,qjs
    real(8) :: ddx,ddy,ddz,d1i,d2i,d3i
	! MPI
	INTEGER(4) MPIchunk_size_loc, iii;
	real(8), allocatable :: MPIbuff1(:), MPIbuff2(:);
	! DBG
	real(4) iniTime, deltaTime;
	! DBG
	! MPI

    do i=1,nm
       do is=1,ns

       slfE0x(is,i) = 0.d0 ! initiating permanent charges fields
       slfE0y(is,i) = 0.d0
       slfE0z(is,i) = 0.d0

       end do
    end do

	! DBG
	iniTime = SECNDS(0.0);
	! DBG

!$omp parallel DEFAULT(SHARED) private(id,is,js,i,iii,j,k)
!    do i=1,nm
!	 do i = MPIrank + 1, nm, MPIcomm_size
	if (MPIcomm_size .GT. 1) then
		i = 1;
	else
		i = 0;
	endif
	do while (i .LE. nm)
	  if (MPIrank .EQ. 0) then
		if (MPIcomm_size .GT. 1) then
			call MPI_RECV(MPIaddr, 1, MPI_INTEGER, MPI_ANY_SOURCE, MPI_TAG_CommReady, MPI_COMM_WORLD, MPIstat, MPIierr);
			call MPI_SEND(i, 1, MPI_INTEGER, MPIaddr, MPI_TAG_SendData, MPI_COMM_WORLD, MPIierr);
		endif
		i = i + 1;
	  else
		if (MPIcomm_size .GT. 1) then
			call MPI_SEND(MPIrank, 1, MPI_INTEGER, 0, MPI_TAG_CommReady, MPI_COMM_WORLD, MPIierr);
			call MPI_RECV(i, 1, MPI_INTEGER, 0, MPI_TAG_SendData, MPI_COMM_WORLD, MPIstat, MPIierr);
		endif
	  endif
 	  if ((i .LE. (nm)) .AND. ((MPIrank .GE. 1) .OR. (MPIcomm_size .EQ. 1))) then
	   !$omp do private(expar2,drr,dist2,qis,qjs,ddx,ddy,ddz,d1i,d2i,d3i) schedule(dynamic)
       do is=ns-npol+1,ns 
          do js=1,ns-npol

             qjs = charge(js)

!             if (js.eq.is) goto 12

!              if (qjs .ne. 0.d0) then

              ddx = (x(i)+xs(is,i)) - (x(i)+xs(js,i)) 
              ddy = (y(i)+ys(is,i)) - (y(i)+ys(js,i)) 
              ddz = (z(i)+zs(is,i)) - (z(i)+zs(js,i)) 

              dist2 = ddx**2+ddy**2+ddz**2

              if ( dist2 .eq. 0.d0 ) goto 10

                 drr = sqrt(dist2)

                 d1i = 1.d0/drr
                 d2i = d1i**2
                 d3i = d2i*d1i

                 expar2 = exp(-(alp*drr)**2.d0)

                 slfE0x(is,i) = slfE0x(is,i) - qjs*ddx*(erf(alp*drr)*d3i-twosqpi*alp*expar2*d2i) 
                 slfE0y(is,i) = slfE0y(is,i) - qjs*ddy*(erf(alp*drr)*d3i-twosqpi*alp*expar2*d2i)
                 slfE0z(is,i) = slfE0z(is,i) - qjs*ddz*(erf(alp*drr)*d3i-twosqpi*alp*expar2*d2i)     

!                 write(*,*) i,is,slfE0x(is,i),slfE0y(is,i),slfE0z(is,i)

!              end if
          
              10 continue     
   
!          12 continue    
          end do

!           write(*,*) i,is,slfE0x(is,i),slfE0y(is,i),slfE0z(is,i)

       end do    
	   !$omp end do
	  endif
    end do
	if (MPIrank .EQ. 0) then
		do i = 1, MPIcomm_size - 1
			call MPI_RECV(MPIaddr, 1, MPI_INTEGER, MPI_ANY_SOURCE, MPI_TAG_CommReady, MPI_COMM_WORLD, MPIstat, MPIierr);
			call MPI_SEND(nm + 1, 1, MPI_INTEGER, MPIaddr, MPI_TAG_SendData, MPI_COMM_WORLD, MPIierr);
		enddo
	endif
!$omp end parallel
	! DBG
	deltaTime = SECNDS(iniTime);
!	write(dbgUnit, *) 'P#', MPIrank, '. self_efield_stat calc time	', deltaTime;
	! DBG
! MPI sync
	j = nm * ns * 3; !slfE0x, slfE0y, slfE0z
	allocate(MPIbuff1(0:j - 1));
	allocate(MPIbuff2(0:j - 1));
	do k = 0, j - 1
		MPIbuff2(k) = 0.0D0;
	enddo
	k = 0;
	do i = 1, nm
		do is = 1, ns
			MPIbuff1(k) = slfE0x(is, i); k = k + 1;
			MPIbuff1(k) = slfE0y(is, i); k = k + 1;
			MPIbuff1(k) = slfE0z(is, i); k = k + 1;
		enddo
	enddo
	call MPI_ALLREDUCE(MPIbuff1, MPIbuff2, j, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, MPIierr);
	k = 0;
	do i = 1, nm
		do is = 1, ns
			slfE0x(is, i) = MPIbuff2(k); k = k + 1;
			slfE0y(is, i) = MPIbuff2(k); k = k + 1;
			slfE0z(is, i) = MPIbuff2(k); k = k + 1;
		enddo
	enddo
	deallocate(MPIbuff1);
	deallocate(MPIbuff2);
! MPI sync

	return;
  end subroutine self_efield_stat


  subroutine surf_efield_stat(xx,yy,zz,xxs,yys,zzs,boxi)
    implicit none
    real(8) :: xx(nom),yy(nom),zz(nom)
    real(8) :: xxs(nsite,nom),yys(nsite,nom),zzs(nsite,nom)
    real(8) :: boxi

    integer :: is,js,i,j,k
    real(8) :: qis
    real(8) :: lsumqrx, lsumqry, lsumqrz
    real(8) :: sumqrx, sumqry, sumqrz
	real(8) :: volmi,fact

    volmi = boxi**3.d0

    fact = 4.0d0*pi/(3.0d0*volmi)

    do i=1,nm
       do is=1,ns
          srfE0x(is,i) = 0.d0 ! initiating permanent charges fields
          srfE0y(is,i) = 0.d0
          srfE0z(is,i) = 0.d0
       end do
    end do

    sumqrx = 0.0D0 ! GRU: these local variables have not been initialized!
    sumqry = 0.0D0
    sumqrz = 0.0D0

!$omp parallel DEFAULT(SHARED) private(id,is,js,i,j,k)&
!$omp& private(qis,lsumqrx, lsumqry, lsumqrz)

    lsumqrx = 0.0D0
    lsumqry = 0.0D0
    lsumqrz = 0.0D0

!$omp do
    do i=1,nm
       do is=1,ns

          qis = charge(is)

          if (qis .ne. 0.d0) then

             lsumqrx = lsumqrx + qis*(xx(i)+xxs(is,i))
             lsumqry = lsumqry + qis*(yy(i)+yys(is,i))
             lsumqrz = lsumqrz + qis*(zz(i)+zzs(is,i))

          end if

       end do
    end do
!$omp end do

!$omp critical
   sumqrx = sumqrx + lsumqrx
   sumqry = sumqry + lsumqry
   sumqrz = sumqrz + lsumqrz
!$omp end critical
!$omp barrier

!$omp do
    do i=1,nm
       do is=ns-npol+1,ns

          srfE0x(is,i) = srfE0x(is,i) - fact*sumqrx
          srfE0y(is,i) = srfE0y(is,i) - fact*sumqry
          srfE0z(is,i) = srfE0z(is,i) - fact*sumqrz

       end do
    end do
!$omp end do
!$omp end parallel

  end subroutine surf_efield_stat

  subroutine iteridm(xx,yy,zz,xxs,yys,zzs,boxi,qindpe)
    implicit none
    real(8) :: xx(nom),yy(nom),zz(nom)
    real(8) :: xxs(nsite,nom),yys(nsite,nom),zzs(nsite,nom)
    real(8) :: boxi,qindpe

    integer :: is,js,i,j,k,isteps
!    real(8) :: Etx(nsite,nom),Ety(nsite,nom),Etz(nsite,nom)  ! total efield 
!    real(8) :: shtEtx(nsite,nom),shtEty(nsite,nom),shtEtz(nsite,nom)  ! short range total efield
    real(8), allocatable :: Etx(:, :), Ety(:, :), Etz(:, :)  ! total efield 
    real(8), allocatable :: shtEtx(:, :), shtEty(:, :), shtEtz(:, :)  ! short range total efield
    real(8) :: efldxi,efldyi,efldzi,efldxj,efldyj,efldzj
    real(8) :: xpolc,ypolc,zpolc !separation between dipole centers
    real(8) :: dist2,rccsq,dist,doti,dotj
    real(8) :: qis,qjs  ! charges on sites
    real(8) :: xcc,ycc,zcc,dx,dy,dz,ddx,ddy,ddz,d1i,d2i,d3i,d5i,d7i
    real(8) :: polr,d3ii,d5ii,expar2,selffac
    real(8) :: thr_iter,change
    real(8) :: energy,engtst
    real(8), parameter :: maxit = 500
    real(8) :: polE1x,polE1y,polE1z
    real(8) lchange,lenergy

    real(8) :: muE0x,muE0y,muE0z,muEix,muEiy,muEiz
    real(8) :: mufE0x,mufE0y,mufE0z,mufEix,mufEiy,mufEiz
	! DBG
	real(4) iniTime, deltaTime;
	real(4) iniTime1, deltaTime1;
	real(4) iniTime2, deltaTime2;
	real(4) iniTime3, deltaTime3;
	real(4) iniTime4, deltaTime4;
	real(4) iniTimeTot, deltaTimeTot;
	real(4) syncTime1, syncTime2, syncTime3, syncTime4, tmpTime;
	! DBG

	! DBG
	deltaTimeTot = 0.0;
	deltaTime = 0.0;
	deltaTime1 = 0.0;
	deltaTime2 = 0.0;
	deltaTime3 = 0.0;
	deltaTime4 = 0.0;
	syncTime1 = 0.0;
	! DBG
    allocate(Etx(ns,nm))
    allocate(Ety(ns,nm))
    allocate(Etz(ns,nm))
    allocate(shtEtx(ns,nm))
    allocate(shtEty(ns,nm))
    allocate(shtEtz(ns,nm))

    do i=1,nm
       do is=1,ns

       idmx(is,i) = 0.d0 ! initiating induced dipoles
       idmy(is,i) = 0.d0
       idmz(is,i) = 0.d0

       E0x(is,i) = 0.d0 ! initiating permanent charges fields
       E0y(is,i) = 0.d0
       E0z(is,i) = 0.d0

       Eidmx(is,i) = 0.d0 ! initiating permanent charges fields
       Eidmy(is,i) = 0.d0
       Eidmz(is,i) = 0.d0  

       Etx(is,i) = 0.d0 ! initiating induced dipoles
       Ety(is,i) = 0.d0
       Etz(is,i) = 0.d0

       end do
    end do

!**************************************************************************
! Now iterate to calculate converged polarization energy        starts here

    thr_iter = 1.d-20
    change = 10.d0
    isteps = 0

    do while (change.gt.thr_iter.and.isteps.lt.maxit)  ! while loop starts

    energy = 0.0d0
    change = 0.d0 

	! DBG
	iniTime = SECNDS(0.0);
	iniTime1 = SECNDS(0.0);
	! DBG
    call rspace_efield_idm(xx,yy,zz,xxs,yys,zzs,alphad,boxi,tmpTime)                    ! output rEidmx,   rEidmy,   rEidmz
	! DBG
	deltaTime1 = deltaTime1 + SECNDS(iniTime1);
	syncTime1 = syncTime1 + tmpTime;
	iniTime2 = SECNDS(0.0);
	! DBG
    call kspace_efield_idm(xx,yy,zz,xxs,yys,zzs,kmax1,kmax2,kmax3,alphad,boxi)  ! output kEidmx,   kEidmy,   kEidmz
	! DBG
	deltaTime2 = deltaTime2 + SECNDS(iniTime2);
	iniTime3 = SECNDS(0.0);
	! DBG
    call self_efield_idm(alphad)                                                ! output slfEidmx, slfEidmy, slfEidmz 
	! DBG
	deltaTime3 = deltaTime3 + SECNDS(iniTime3);
	! DBG

    if ( induction_type .eq. 'damped' ) then
		! DBG
		iniTime4 = SECNDS(0.0);
		! DBG
       call shortrange_efield_idm_damped(xx,yy,zz,xxs,yys,zzs,boxi) ! output shtE0x, shtE0y, shtE0z
		! DBG
		deltaTime4 = deltaTime4 + SECNDS(iniTime4);
		! DBG
    end if 
	! DBG
	deltaTime = deltaTime + SECNDS(iniTime);
	iniTimeTot = SECNDS(0.0);
	! DBG

!$omp parallel DEFAULT(SHARED)&
!$omp& private(is,i,polr,polE1x,polE1y,polE1z,lchange,lenergy)
    lchange = 0.0D0
    lenergy = 0.0D0
!$omp do schedule(dynamic)
      do i=1,nm    ! loop over ith molecule starts

         do is=ns-npol+1,ns 

            E0x(is,i) = rE0x(is,i)+kE0x(is,i)+slfE0x(is,i) + shtE0x(is,i)
            E0y(is,i) = rE0y(is,i)+kE0y(is,i)+slfE0y(is,i) + shtE0y(is,i)
            E0z(is,i) = rE0z(is,i)+kE0z(is,i)+slfE0z(is,i) + shtE0z(is,i)  

            Eidmx(is,i) = rEidmx(is,i)+kEidmx(is,i)+slfEidmx(is,i) + shtEidmx(is,i)
            Eidmy(is,i) = rEidmy(is,i)+kEidmy(is,i)+slfEidmy(is,i) + shtEidmy(is,i)
            Eidmz(is,i) = rEidmz(is,i)+kEidmz(is,i)+slfEidmz(is,i) + shtEidmz(is,i)

            Etx(is,i) =  Eidmx(is,i) + E0x(is,i) ! total efield from Ewald sum 
            Ety(is,i) =  Eidmy(is,i) + E0y(is,i) 
            Etz(is,i) =  Eidmz(is,i) + E0z(is,i) 

!            write(*,*)i,is,E0x(is,i),E0y(is,i),E0z(is,i),Eidmx(is,i),Eidmy(is,i),Eidmz(is,i) 

            polr = apol(is,i)

            polE1x = polr* Etx(is,i) 
            polE1y = polr* Ety(is,i) 
            polE1z = polr* Etz(is,i)   
 
            lchange = (idmx(is,i)-polE1x)**2 + &
                      (idmy(is,i)-polE1y)**2 + &
                      (idmz(is,i)-polE1z)**2 + lchange

            idmx(is,i) = polE1x
            idmy(is,i) = polE1y
            idmz(is,i) = polE1z

            lenergy = -0.5d0*polr*(Etx(is,i)*E0x(is,i) + &
                                  Ety(is,i)*E0y(is,i) + &
                                  Etz(is,i)*E0z(is,i)) +  lenergy

!            energy = -0.5d0*(idmx(is,i)*E0x(is,i) + &
!                             idmy(is,i)*E0y(is,i) + &
!                             idmz(is,i)*E0z(is,i)) +  energy

         end do   ! loop over is 

      end do  ! loop over ith molecule ends  
!$omp end do
!$omp critical
      change = change + lchange
      energy = energy + lenergy
!$omp end critical
!$omp end parallel

      isteps = isteps + 1

		! DBG
		deltaTimeTot = deltaTimeTot + SECNDS(iniTimeTot);
		! DBG
    end do       ! while loop ends

    if (isteps.ge.maxit) then
        write (*,*) 'No convergence in indN_iter'
        write (*,'(a,g12.3)') 'energy change=',change
        write (*,'(a,g12.3)') 'thr_iter=',thr_iter
        stop
    end if

    qindpe = r4pie0*energy/418.4d0

!    write(*,*) qindpe
!    vir = vir - qindpe*418.4d0

!    write(*,*) 'vir_ind', -qindpe*418.4d0

!**************************************************************************
! calculation of converged polarization energy        ends here

    deallocate(Etx)
    deallocate(Ety)
    deallocate(Etz)
    deallocate(shtEtx)
    deallocate(shtEty)
    deallocate(shtEtz)
	! DBG
!	write(dbgUnit, *) 'P#', MPIrank, '. rspace_efield_idm calc time	', deltaTime1;
!	write(dbgUnit, *) 'P#', MPIrank, '. rspace_efield_idm sync time	', syncTime1;
!	write(dbgUnit, *) 'P#', MPIrank, '. kspace_efield_idm calc time	', deltaTime2;
!	write(dbgUnit, *) 'P#', MPIrank, '. self_efield_idm calc time	', deltaTime3;
!	write(dbgUnit, *) 'P#', MPIrank, '. shortrange_efield_idm_damped calc time	', deltaTime4;
!	write(dbgUnit, *) 'P#', MPIrank, '. iteridm calc1 time	', deltaTime;
!	write(dbgUnit, *) 'P#', MPIrank, '. iteridm calc2 time	', deltaTimeTot;
	! DBG

  end subroutine iteridm

  subroutine rspace_efield_idm(xx,yy,zz,xxs,yys,zzs,alp,boxi,rdeltaTime)
    implicit none
    real(8) :: xx(nom),yy(nom),zz(nom)
    real(8) :: xxs(nsite,nom),yys(nsite,nom),zzs(nsite,nom)
	real(8) :: alp,boxi
	! DBG
	real(4) rdeltaTime;
	! DBG

    integer :: is,js,i,j,k,isteps
    real(8) :: efldxi,efldyi,efldzi,efldxj,efldyj,efldzj
    real(8) :: dist2,rccsq,dist,doti,dotj,d3ii,d5ii,drr,expar2
    real(8) :: qis,qjs  ! charges on sites
    real(8) :: xcc,ycc,zcc,dx,dy,dz,ddx,ddy,ddz,d1i,d2i,d3i,d5i,d7i
	! MPI
	INTEGER(4) MPIchunk_size_loc, iii;
	real(8), allocatable :: MPIbuff1(:), MPIbuff2(:);
	! DBG
	real(4) iniTime, deltaTime;
	! DBG
	! MPI

    do i=1,nm
       do is=1,ns

       rEidmx(is,i) = 0.d0 ! initiating permanent charges fields
       rEidmy(is,i) = 0.d0
       rEidmz(is,i) = 0.d0

       end do
    end do
	! DBG
	iniTime = SECNDS(0.0);
	! DBG

!$omp parallel DEFAULT(SHARED) private(id,is,js,i,iii,j,k,isteps)
!    do i=1,nm   ! loop over ith molecule starts
    do i = MPIrank + 1, nm, MPIcomm_size
		!$omp do private(efldxi,efldyi,efldzi,efldxj,efldyj,efldzj,dist2,rccsq,dist,doti,dotj,d3ii,d5ii) &
		!$omp& private(drr,expar2,qis,qjs,xcc,ycc,zcc,dx,dy,dz,ddx,ddy,ddz,d1i,d2i,d3i,d5i,d7i) schedule(dynamic)
          do j=1,nm   ! loop over jth molecule starts

!             do is=ns-npol+1,ns

!             if (j.eq.i) goto 13 !cycle
             if (j .NE. i) then

                dx = xx(i)-xx(j)
                dy = yy(i)-yy(j)
                dz = zz(i)-zz(j)

                xcc = dx - boxi*nint(dx/boxi)
                ycc = dy - boxi*nint(dy/boxi)
                zcc = dz - boxi*nint(dz/boxi)

                rccsq = xcc**2 + ycc**2 + zcc**2

                if (rccsq.lt.rcutsd) then  ! com-com cutoff check starts

                   do is=ns-npol+1,ns

                   do js=ns-npol+1,ns    ! loop over sites of molecule j  Starts

                      ddx = xcc + xxs(is,i) - xxs(js,j)
                      ddy = ycc + yys(is,i) - yys(js,j)
                      ddz = zcc + zzs(is,i) - zzs(js,j)

                      dist2 = ddx**2 + ddy**2 + ddz**2

!                      if (dist2.lt.rcutsd) then ! cut off on  sites separation
!                      starts

                      drr = sqrt(dist2)
                      d1i = 1.d0/drr
                      d2i = d1i**2
                      d3i = d2i*d1i
                      d5i = d3i*d2i

                      ! dot product of induced dipole and separation vector
                      doti = idmx(is,i)*ddx+idmy(is,i)*ddy+idmz(is,i)*ddz
                      dotj = idmx(js,j)*ddx+idmy(js,j)*ddy+idmz(js,j)*ddz 

                      expar2 = exp(-(alp*drr)**2.d0)

                      d3ii = erfc(alp*drr)*d3i + twosqpi*alp*expar2*d2i
      
                      d5ii = erfc(alp*drr)*d5i + twosqpi*alp*expar2*d2i*d2i + &
                             2.d0*twosqpi*alp**(3.d0)*expar2*d2i/3.d0 

                      ! calculate the  efield of inducd dipole of molecule j at
                      ! induced dipole of molecule i

                      efldxi = 3.d0*d5ii*dotj*ddx - idmx(js,j)*d3ii
                      efldyi = 3.d0*d5ii*dotj*ddy - idmy(js,j)*d3ii
                      efldzi = 3.d0*d5ii*dotj*ddz - idmz(js,j)*d3ii

                      ! calculate total efield at induced dipole at molecule i
                      ! because of parmenant charges
                      !  and induced dipoles at molecule j

                      rEidmx(is,i) = rEidmx(is,i) + efldxi
                      rEidmy(is,i) = rEidmy(is,i) + efldyi
                      rEidmz(is,i) = rEidmz(is,i) + efldzi

!                       end if

                   end do     ! loop over sites of molecule j  ends

                   end do

                end if   ! com-com cutoff check ends

!            end do

             endif
!         13 continue
         end do      ! loop over jth molecule ends
		!$omp end do
    end do  ! loop over ith molecule ends
!$omp end parallel
	! DBG
	deltaTime = SECNDS(iniTime);
!	write(dbgUnit, *) 'P#', MPIrank, '. rspace_efield_idm calc time	', deltaTime;
	! DBG
! MPI sync
	! DBG
	iniTime = SECNDS(0.0);
	! DBG
	j = nm * ns * 3; !rEidmx, rEidmy, rEidmz
	allocate(MPIbuff1(0:j - 1));
	allocate(MPIbuff2(0:j - 1));
	do k = 0, j - 1
		MPIbuff2(k) = 0.0D0;
	enddo
	k = 0;
	do i = 1, nm
		do is = 1, ns
			MPIbuff1(k) = rEidmx(is, i); k = k + 1;
			MPIbuff1(k) = rEidmy(is, i); k = k + 1;
			MPIbuff1(k) = rEidmz(is, i); k = k + 1;
		enddo
	enddo
	call MPI_ALLREDUCE(MPIbuff1, MPIbuff2, j, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, MPIierr);
	k = 0;
	do i = 1, nm
		do is = 1, ns
			rEidmx(is, i) = MPIbuff2(k); k = k + 1;
			rEidmy(is, i) = MPIbuff2(k); k = k + 1;
			rEidmz(is, i) = MPIbuff2(k); k = k + 1;
		enddo
	enddo
	deallocate(MPIbuff1);
	deallocate(MPIbuff2);
	! DBG
	rdeltaTime = SECNDS(iniTime);
!	write(dbgUnit, *) 'P#', MPIrank, '. rspace_efield_idm sync time	', deltaTime;
	! DBG
! MPI sync

	return;
  end subroutine rspace_efield_idm

  subroutine rspace_efield_idm_v0(xx,yy,zz,xxs,yys,zzs,alp,boxi,rdeltaTime)
    implicit none
    real(8) :: xx(nom),yy(nom),zz(nom)
    real(8) :: xxs(nsite,nom),yys(nsite,nom),zzs(nsite,nom)
	real(8) :: alp,boxi
	! DBG
	real(4) rdeltaTime;
	! DBG

    integer :: is,js,i,j,k,isteps
    real(8) :: efldxi,efldyi,efldzi,efldxj,efldyj,efldzj
    real(8) :: dist2,rccsq,dist,doti,dotj,d3ii,d5ii,drr,expar2
    real(8) :: qis,qjs  ! charges on sites
    real(8) :: xcc,ycc,zcc,dx,dy,dz,ddx,ddy,ddz,d1i,d2i,d3i,d5i,d7i
	! MPI
	INTEGER(4) MPIchunk_size_loc, iii;
	real(8), allocatable :: MPIbuff1(:), MPIbuff2(:);
	! DBG
	real(4) iniTime, deltaTime;
	! DBG
	! MPI

    do i=1,nm
       do is=1,ns

       rEidmx(is,i) = 0.d0 ! initiating permanent charges fields
       rEidmy(is,i) = 0.d0
       rEidmz(is,i) = 0.d0

       end do
    end do
	! DBG
	iniTime = SECNDS(0.0);
	! DBG

!$omp parallel DEFAULT(SHARED) private(id,is,js,i,iii,j,k,isteps)
!    do i=1,nm   ! loop over ith molecule starts
    do i = MPIrank + 1, nm, MPIcomm_size
		!$omp do private(efldxi,efldyi,efldzi,efldxj,efldyj,efldzj,dist2,rccsq,dist,doti,dotj,d3ii,d5ii) &
		!$omp& private(drr,expar2,qis,qjs,xcc,ycc,zcc,dx,dy,dz,ddx,ddy,ddz,d1i,d2i,d3i,d5i,d7i) schedule(dynamic)
          do j=1,nm   ! loop over jth molecule starts

!             do is=ns-npol+1,ns

!             if (j.eq.i) goto 13 !cycle
             if (j .NE. i) then

                dx = xx(i)-xx(j)
                dy = yy(i)-yy(j)
                dz = zz(i)-zz(j)

                xcc = dx - boxi*nint(dx/boxi)
                ycc = dy - boxi*nint(dy/boxi)
                zcc = dz - boxi*nint(dz/boxi)

                rccsq = xcc**2 + ycc**2 + zcc**2

                if (rccsq.lt.rcutsd) then  ! com-com cutoff check starts

                   do is=ns-npol+1,ns

                   do js=ns-npol+1,ns    ! loop over sites of molecule j  Starts

                      ddx = xcc + xxs(is,i) - xxs(js,j)
                      ddy = ycc + yys(is,i) - yys(js,j)
                      ddz = zcc + zzs(is,i) - zzs(js,j)

                      dist2 = ddx**2 + ddy**2 + ddz**2

!                      if (dist2.lt.rcutsd) then ! cut off on  sites separation
!                      starts

                      drr = sqrt(dist2)
                      d1i = 1.d0/drr
                      d2i = d1i**2
                      d3i = d2i*d1i
                      d5i = d3i*d2i

                      ! dot product of induced dipole and separation vector
                      doti = idmx(is,i)*ddx+idmy(is,i)*ddy+idmz(is,i)*ddz
                      dotj = idmx(js,j)*ddx+idmy(js,j)*ddy+idmz(js,j)*ddz 

                      expar2 = exp(-(alp*drr)**2.d0)

                      d3ii = erfc(alp*drr)*d3i + twosqpi*alp*expar2*d2i
      
                      d5ii = erfc(alp*drr)*d5i + twosqpi*alp*expar2*d2i*d2i + &
                             2.d0*twosqpi*alp**(3.d0)*expar2*d2i/3.d0 

                      ! calculate the  efield of inducd dipole of molecule j at
                      ! induced dipole of molecule i

                      efldxi = 3.d0*d5ii*dotj*ddx - idmx(js,j)*d3ii
                      efldyi = 3.d0*d5ii*dotj*ddy - idmy(js,j)*d3ii
                      efldzi = 3.d0*d5ii*dotj*ddz - idmz(js,j)*d3ii

                      ! calculate total efield at induced dipole at molecule i
                      ! because of parmenant charges
                      !  and induced dipoles at molecule j

                      rEidmx(is,i) = rEidmx(is,i) + efldxi
                      rEidmy(is,i) = rEidmy(is,i) + efldyi
                      rEidmz(is,i) = rEidmz(is,i) + efldzi

!                       end if

                   end do     ! loop over sites of molecule j  ends

                   end do

                end if   ! com-com cutoff check ends

!            end do

             endif
!         13 continue
         end do      ! loop over jth molecule ends
		!$omp end do
    end do  ! loop over ith molecule ends
!$omp end parallel
	! DBG
	deltaTime = SECNDS(iniTime);
!	write(dbgUnit, *) 'P#', MPIrank, '. rspace_efield_idm calc time	', deltaTime;
	! DBG
! MPI sync
	! DBG
	iniTime = SECNDS(0.0);
	! DBG
	j = nm * ns * 3; !rEidmx, rEidmy, rEidmz
	allocate(MPIbuff1(0:j - 1));
	allocate(MPIbuff2(0:j - 1));
	do k = 0, j - 1
		MPIbuff2(k) = 0.0D0;
	enddo
	k = 0;
	do i = 1, nm
		do is = 1, ns
			MPIbuff1(k) = rEidmx(is, i); k = k + 1;
			MPIbuff1(k) = rEidmy(is, i); k = k + 1;
			MPIbuff1(k) = rEidmz(is, i); k = k + 1;
		enddo
	enddo
	call MPI_ALLREDUCE(MPIbuff1, MPIbuff2, j, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, MPIierr);
	k = 0;
	do i = 1, nm
		do is = 1, ns
			rEidmx(is, i) = MPIbuff2(k); k = k + 1;
			rEidmy(is, i) = MPIbuff2(k); k = k + 1;
			rEidmz(is, i) = MPIbuff2(k); k = k + 1;
		enddo
	enddo
	deallocate(MPIbuff1);
	deallocate(MPIbuff2);
	! DBG
	rdeltaTime = SECNDS(iniTime);
!	write(dbgUnit, *) 'P#', MPIrank, '. rspace_efield_idm sync time	', deltaTime;
	! DBG
! MPI sync

	return;
  end subroutine rspace_efield_idm_v0

  subroutine rspace_efield_idm_v1(xx,yy,zz,xxs,yys,zzs,alp,boxi)
    implicit none
    real(8) :: xx(nom),yy(nom),zz(nom)
    real(8) :: xxs(nsite,nom),yys(nsite,nom),zzs(nsite,nom)
	real(8) :: alp,boxi

    integer :: is,js,i,j,k,isteps
    real(8) :: efldxi,efldyi,efldzi,efldxj,efldyj,efldzj
    real(8) :: dist2,rccsq,dist,doti,dotj,d3ii,d5ii,drr,expar2
    real(8) :: qis,qjs  ! charges on sites
    real(8) :: xcc,ycc,zcc,dx,dy,dz,ddx,ddy,ddz,d1i,d2i,d3i,d5i,d7i
	! MPI
	INTEGER(4) MPIchunk_size_loc, iii;
	real(8), allocatable :: MPIbuff1(:), MPIbuff2(:);
	! DBG
	real(4) iniTime, deltaTime;
	! DBG
	! MPI

    do i=1,nm
       do is=1,ns

       rEidmx(is,i) = 0.d0 ! initiating permanent charges fields
       rEidmy(is,i) = 0.d0
       rEidmz(is,i) = 0.d0

       end do
    end do
	! DBG
	iniTime = SECNDS(0.0);
	! DBG

!$omp parallel DEFAULT(SHARED) private(id,is,js,i,iii,j,k,isteps)
!    do i=1,nm   ! loop over ith molecule starts
!	 do i = MPIrank + 1, nm, MPIcomm_size
	if (MPIcomm_size .GT. 1) then
		i = 1;
	else
		i = 0;
	endif
	do while (i .LE. nm)
	  if (MPIrank .EQ. 0) then
		if (MPIcomm_size .GT. 1) then
			call MPI_RECV(MPIaddr, 1, MPI_INTEGER, MPI_ANY_SOURCE, MPI_TAG_CommReady, MPI_COMM_WORLD, MPIstat, MPIierr);
			call MPI_SEND(i, 1, MPI_INTEGER, MPIaddr, MPI_TAG_SendData, MPI_COMM_WORLD, MPIierr);
		endif
		i = i + 1;
	  else
		if (MPIcomm_size .GT. 1) then
			call MPI_SEND(MPIrank, 1, MPI_INTEGER, 0, MPI_TAG_CommReady, MPI_COMM_WORLD, MPIierr);
			call MPI_RECV(i, 1, MPI_INTEGER, 0, MPI_TAG_SendData, MPI_COMM_WORLD, MPIstat, MPIierr);
		endif
	  endif
 	  if ((i .LE. (nm)) .AND. ((MPIrank .GE. 1) .OR. (MPIcomm_size .EQ. 1))) then
		!$omp do private(efldxi,efldyi,efldzi,efldxj,efldyj,efldzj,dist2,rccsq,dist,doti,dotj,d3ii,d5ii) &
		!$omp& private(drr,expar2,qis,qjs,xcc,ycc,zcc,dx,dy,dz,ddx,ddy,ddz,d1i,d2i,d3i,d5i,d7i) schedule(dynamic)
          do j=1,nm   ! loop over jth molecule starts

!             do is=ns-npol+1,ns

!             if (j.eq.i) goto 13 !cycle
             if (j .NE. i) then

                dx = xx(i)-xx(j)
                dy = yy(i)-yy(j)
                dz = zz(i)-zz(j)

                xcc = dx - boxi*nint(dx/boxi)
                ycc = dy - boxi*nint(dy/boxi)
                zcc = dz - boxi*nint(dz/boxi)

                rccsq = xcc**2 + ycc**2 + zcc**2

                if (rccsq.lt.rcutsd) then  ! com-com cutoff check starts

                   do is=ns-npol+1,ns

                   do js=ns-npol+1,ns    ! loop over sites of molecule j  Starts

                      ddx = xcc + xxs(is,i) - xxs(js,j)
                      ddy = ycc + yys(is,i) - yys(js,j)
                      ddz = zcc + zzs(is,i) - zzs(js,j)

                      dist2 = ddx**2 + ddy**2 + ddz**2

!                      if (dist2.lt.rcutsd) then ! cut off on  sites separation
!                      starts

                      drr = sqrt(dist2)
                      d1i = 1.d0/drr
                      d2i = d1i**2
                      d3i = d2i*d1i
                      d5i = d3i*d2i

                      ! dot product of induced dipole and separation vector
                      doti = idmx(is,i)*ddx+idmy(is,i)*ddy+idmz(is,i)*ddz
                      dotj = idmx(js,j)*ddx+idmy(js,j)*ddy+idmz(js,j)*ddz 

                      expar2 = exp(-(alp*drr)**2.d0)

                      d3ii = erfc(alp*drr)*d3i + twosqpi*alp*expar2*d2i
      
                      d5ii = erfc(alp*drr)*d5i + twosqpi*alp*expar2*d2i*d2i + &
                             2.d0*twosqpi*alp**(3.d0)*expar2*d2i/3.d0 

                      ! calculate the  efield of inducd dipole of molecule j at
                      ! induced dipole of molecule i

                      efldxi = 3.d0*d5ii*dotj*ddx - idmx(js,j)*d3ii
                      efldyi = 3.d0*d5ii*dotj*ddy - idmy(js,j)*d3ii
                      efldzi = 3.d0*d5ii*dotj*ddz - idmz(js,j)*d3ii

                      ! calculate total efield at induced dipole at molecule i
                      ! because of parmenant charges
                      !  and induced dipoles at molecule j

                      rEidmx(is,i) = rEidmx(is,i) + efldxi
                      rEidmy(is,i) = rEidmy(is,i) + efldyi
                      rEidmz(is,i) = rEidmz(is,i) + efldzi

!                       end if

                   end do     ! loop over sites of molecule j  ends

                   end do

                end if   ! com-com cutoff check ends

!            end do

             endif
!         13 continue
         end do      ! loop over jth molecule ends
		!$omp end do
	   endif
    end do  ! loop over ith molecule ends
	if (MPIrank .EQ. 0) then
		do i = 1, MPIcomm_size - 1
			call MPI_RECV(MPIaddr, 1, MPI_INTEGER, MPI_ANY_SOURCE, MPI_TAG_CommReady, MPI_COMM_WORLD, MPIstat, MPIierr);
			call MPI_SEND(nm + 1, 1, MPI_INTEGER, MPIaddr, MPI_TAG_SendData, MPI_COMM_WORLD, MPIierr);
		enddo
	endif
!$omp end parallel
	! DBG
	deltaTime = SECNDS(iniTime);
!	write(dbgUnit, *) 'P#', MPIrank, '. rspace_efield_idm calc time	', deltaTime;
	! DBG
! MPI sync
	j = nm * ns * 3; !rEidmx, rEidmy, rEidmz
	allocate(MPIbuff1(0:j - 1));
	allocate(MPIbuff2(0:j - 1));
	do k = 0, j - 1
		MPIbuff2(k) = 0.0D0;
	enddo
	k = 0;
	do i = 1, nm
		do is = 1, ns
			MPIbuff1(k) = rEidmx(is, i); k = k + 1;
			MPIbuff1(k) = rEidmy(is, i); k = k + 1;
			MPIbuff1(k) = rEidmz(is, i); k = k + 1;
		enddo
	enddo
	call MPI_ALLREDUCE(MPIbuff1, MPIbuff2, j, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, MPIierr);
	k = 0;
	do i = 1, nm
		do is = 1, ns
			rEidmx(is, i) = MPIbuff2(k); k = k + 1;
			rEidmy(is, i) = MPIbuff2(k); k = k + 1;
			rEidmz(is, i) = MPIbuff2(k); k = k + 1;
		enddo
	enddo
	deallocate(MPIbuff1);
	deallocate(MPIbuff2);
! MPI sync

	return;
  end subroutine rspace_efield_idm_v1

     
  subroutine kspace_efield_idm2(xx,yy,zz,xxs,yys,zzs,lmax,mmax,nmax,alp,boxi)
    implicit none
    real(8) :: xx(nom),yy(nom),zz(nom)
    real(8) :: xxs(nsite,nom),yys(nsite,nom),zzs(nsite,nom)
    integer :: lmax,mmax,nmax 
	real(8) :: alp,boxi

    real(8) :: qkcpe,qvir_kc
    real(8) :: volmi 
    integer :: i,j,is,js
!    integer, parameter :: nqmax = nsite*nom  
    integer :: nqmax ! = ns*nm ! This should be enough
    integer :: klmmax,klmmaxsq 
!    complex(8) :: expikr(nqmax),el(nqmax,0:20)
!    complex(8) :: em(nqmax,-20:20),en(nqmax,-20:20)
!    real(8) :: qi(nqmax),xi(nqmax),yi(nqmax),zi(nqmax)
!    real(8) :: dk,idmxi(nqmax),idmyi(nqmax),idmzi(nqmax)
    complex(8), allocatable :: expikr(:), el(:, :)
    complex(8), allocatable :: em(:, :), en(:, :)
    real(8), allocatable :: qi(:), xi(:), yi(:), zi(:)
    real(8), allocatable :: idmxi(:), idmyi(:), idmzi(:)
    real(8) :: dk,qpe,rl,rm,rn,rksq,qforce,muefld
    complex(8) :: sumqex
    real(8), allocatable :: lkEidmx(:, :), lkEidmy(:, :), lkEidmz(:, :)

    real(8), parameter :: twopi = 2.d0*pi
    integer :: noq
!    real(8) :: a(ksqmax)
    real(8), allocatable :: a(:)
    integer :: llim,mlim,nlim,mmin,nmin

    integer :: ll,mm,nn,kk,lmnp
	! MPI
	INTEGER(4) MPIchunk_size_loc, iii, k;
	INTEGER(4) :: MPIbuffI(0:2); ! ll, mm, nn
	real(8), allocatable :: MPIbuff1(:), MPIbuff2(:);
	! DBG
	real(4) iniTime, deltaTime;
	logical flag;
	! DBG
	! MPI

    nqmax = ns*nm

    allocate(el(nqmax,0:20))
    allocate(em(nqmax,-20:20))
    allocate(en(nqmax,-20:20))
    allocate(qi(nqmax))
    allocate(xi(nqmax))
    allocate(yi(nqmax))
    allocate(zi(nqmax))
    allocate(idmxi(nqmax))
    allocate(idmyi(nqmax))
    allocate(idmzi(nqmax))

    do i=1,nm
       do is=1,ns

       kEidmx(is,i) = 0.d0 ! initiating permanent charges fields
       kEidmy(is,i) = 0.d0
       kEidmz(is,i) = 0.d0

       end do
    end do
!$omp parallel DEFAULT(SHARED) private(id,i,iii,j,is,js,expikr,a,llim,mlim,nlim,mmin,nmin)&
!$omp private(ll,mm,nn,kk,lmnp,volmi,klmmaxsq,klmmax,lkEidmx,lkEidmy,lkEidmz)

   allocate(expikr(nqmax))
   allocate(a(ksqmax))
   allocate(lkEidmx(ns,nm))
   allocate(lkEidmy(ns,nm))
   allocate(lkEidmz(ns,nm))

    klmmax = min(dble(lmax),dble(mmax),dble(nmax))
    klmmaxsq = klmmax**2

    llim = lmax 
    mlim = mmax
    nlim = nmax

!    qpe = 0.d0 ! GRU: Not used local variable

    volmi = boxi**3
    
    do i=1,nm
       do is=1,ns
          lkEidmx(is,i) = 0.d0
          lkEidmy(is,i) = 0.d0
          lkEidmz(is,i) = 0.d0
       end do
    end do

!    collect quantities for charged sites

!$omp MASTER
      i = 0

      do j = 1,nm
         do js = 1,ns
               i = i+1
               idmxi(i) = idmx(js,j)
               idmyi(i) = idmy(js,j)
               idmzi(i) = idmz(js,j)
               xi(i) = xx(j) + xxs(js,j)
               yi(i) = yy(j) + yys(js,j)
               zi(i) = zz(j) + zzs(js,j)
         enddo
      enddo 

      noq = i
!$omp end MASTER
!$omp barrier

!$omp do
    do i = 1,noq
       el(i,0) = (1.0,0.0)
       em(i,0) = (1.0,0.0)
       en(i,0) = (1.0,0.0)
       el(i,1) = exp((0.,1.)*twopi*xi(i)/boxi)
       em(i,1) = exp((0.,1.)*twopi*yi(i)/boxi)
       en(i,1) = exp((0.,1.)*twopi*zi(i)/boxi)
    enddo
!$omp end do
!$omp barrier

!$omp MASTER
    do ll = 2,llim
         do i = 1,noq
            el(i,ll) = el(i,ll-1)*el(i,1)
         enddo
      enddo
      do mm = 2,mlim
         do i = 1,noq
            em(i,mm) = em(i,mm-1)*em(i,1)
         enddo
      enddo
      do nn = 2,nlim
         do i = 1,noq
            en(i,nn) = en(i,nn-1)*en(i,1)
         enddo
      enddo
      do mm = -mlim,-1
         do i = 1,noq
            em(i,mm) = conjg(em(i,-mm))
         enddo
      enddo
      do nn = -nlim,-1
         do i = 1,noq
            en(i,nn) = conjg(en(i,-nn))
         enddo
      enddo   
!$omp end MASTER
!$omp barrier

	! DBG
	iniTime = SECNDS(0.0);
	! DBG

    mmin = 0
    nmin = 1

!    do ll=0,llim
!     do ll = MPIrank, llim, MPIcomm_size
	ll = 0;
	mm = mmin;
	nn = nmin - 1;
!	if (MPIcomm_size .GT. 1) ll = -1;
	do while (ll .LE. llim)
	  if (MPIrank .EQ. 0) then
!		do ll = 0, llim
!			if (ll .GT. 0) mmin = -mlim
!			do mm=mmin,mlim
!				do nn=nmin,nlim
!				enddo
!				nmin = -nlim;
!			enddo
!		enddo
		! Linearization
		nn = nn + 1;
		if (nn .GT. nlim) then
			nmin = -nlim;
			nn = nmin;
			mm = mm + 1;
		endif
		if (mm .GT. mlim) then
			mmin = -mlim;
			mm = mmin;
			ll = ll + 1;
		endif
!		if (ll .GT. 0) mmin = -mlim;
		MPIbuffI(0) = ll;
		MPIbuffI(1) = mm;
		MPIbuffI(2) = nn;
		if (MPIcomm_size .GT. 1) then
			call MPI_RECV(MPIaddr, 1, MPI_INTEGER, MPI_ANY_SOURCE, MPI_TAG_CommReady, MPI_COMM_WORLD, MPIstat, MPIierr);
			call MPI_SEND(MPIbuffI, SIZE(MPIbuffI), MPI_INTEGER, MPIaddr, MPI_TAG_SendData, MPI_COMM_WORLD, MPIierr);
		endif
	  else
		call MPI_SEND(MPIrank, 1, MPI_INTEGER, 0, MPI_TAG_CommReady, MPI_COMM_WORLD, MPIierr);
		call MPI_RECV(MPIbuffI, SIZE(MPIbuffI), MPI_INTEGER, 0, MPI_TAG_SendData, MPI_COMM_WORLD, MPIstat, MPIierr);
		ll = MPIbuffI(0);
		mm = MPIbuffI(1);
		nn = MPIbuffI(2);
	  endif
 	  if ((ll .LE. (llim)) .AND. ((MPIrank .GE. 1) .OR. (MPIcomm_size .EQ. 1))) then
		! DBG
!		flag = .TRUE.;
		! DBG
       
       rl = twopi*ll/boxi
       rm = twopi*mm/boxi
       rn = twopi*nn/boxi
       kk=ll**2+mm**2+nn**2     
       lmnp = ll+mm+nn 

!       if ((kk.le.klmmaxsq) .and. .not.(mod(lmnp,2).ne.0) ) then 
       if ((kk.le.klmmaxsq)) then
			! DBG
			! if (flag) then
				! flag = .FALSE.;
				! write(500, '(5I6)') ll, mmin, mlim, nmin, nlim;
			! endif
			! write(501, '(5I6)') ll, mm, nn, kk;
			! DBG

			rksq = rl*rl + rm*rm + rn*rn

			a(kk) = twopi/volmi*exp(-0.25d0*rksq/alp**2)/rksq             

			! form expikr for each charge
			do i = 1,noq
			   expikr(i) = el(i,ll)*em(i,mm)*en(i,nn)
			enddo      

			 !form sum of (idm dot k)*expikr

			sumqex = (0.0,0.0)
			do i = 1,noq
			   dk = rl*idmxi(i) + rm*idmyi(i) + rn*idmzi(i)
			   sumqex = sumqex+dk*expikr(i)
           !                  sumqex = sumqex+conjg(cmplx(0.d0,dk)*expikr(i))
			enddo

            ! accumulate potential energy

			   i=0
			   do j=1,nm
				  do js=1,ns
			  !       if (chgs(js,j).ne.0.0) then
						i=i+1

                                  muefld = -4.0*a(kk)*sumqex*conjg(expikr(i))

           !!                        muefld = -4.0*a(kk)*aimag(sumqex*expikr(i))

						lkEidmx(js,j) = lkEidmx(js,j) + rl*muefld
						lkEidmy(js,j) = lkEidmy(js,j) + rm*muefld
						lkEidmz(js,j) = lkEidmz(js,j) + rn*muefld
					  
			  !       end if
				  end do 
			   end do
            end if
	  endif
    end do
	if (MPIrank .EQ. 0) then
		MPIbuffI(0) = llim + 1;
		MPIbuffI(1) = 0;
		MPIbuffI(2) = 0;
		do i = 1, MPIcomm_size - 2
			call MPI_RECV(MPIaddr, 1, MPI_INTEGER, MPI_ANY_SOURCE, MPI_TAG_CommReady, MPI_COMM_WORLD, MPIstat, MPIierr);
			call MPI_SEND(MPIbuffI, SIZE(MPIbuffI), MPI_INTEGER, MPIaddr, MPI_TAG_SendData, MPI_COMM_WORLD, MPIierr);
		enddo
	endif
	! DBG
	deltaTime = SECNDS(iniTime);
!	write(dbgUnit, *) 'P#', MPIrank, '. kspace_efield_idm calc time	', deltaTime;
!	call MPI_BARRIER(MPI_COMM_WORLD, MPIierr);
!	STOP;
	! DBG

!$omp critical
   do i=1,nm
      do is=1,ns
         kEidmx(is,i) = kEidmx(is,i) + lkEidmx(is,i)
         kEidmy(is,i) = kEidmy(is,i) + lkEidmy(is,i)
         kEidmz(is,i) = kEidmz(is,i) + lkEidmz(is,i)
      end do
   end do
!$omp end critical

   deallocate(expikr)
   deallocate(a)
   deallocate(lkEidmx)
   deallocate(lkEidmy)
   deallocate(lkEidmz)

!$omp end parallel

   deallocate(el)
   deallocate(em)
   deallocate(en)
   deallocate(qi)
   deallocate(xi)
   deallocate(yi)
   deallocate(zi)
   deallocate(idmxi)
   deallocate(idmyi)
   deallocate(idmzi)
! MPI sync
	j = nm * ns * 3; !kEidmx, kEidmy, kEidmz
	allocate(MPIbuff1(0:j - 1));
	allocate(MPIbuff2(0:j - 1));
	do k = 0, j - 1
		MPIbuff2(k) = 0.0D0;
	enddo
	k = 0;
	do i = 1, nm
		do is = 1, ns
			MPIbuff1(k) = kEidmx(is, i); k = k + 1;
			MPIbuff1(k) = kEidmy(is, i); k = k + 1;
			MPIbuff1(k) = kEidmz(is, i); k = k + 1;
		enddo
	enddo
	call MPI_ALLREDUCE(MPIbuff1, MPIbuff2, j, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, MPIierr);
	k = 0;
	do i = 1, nm
		do is = 1, ns
			kEidmx(is, i) = MPIbuff2(k); k = k + 1;
			kEidmy(is, i) = MPIbuff2(k); k = k + 1;
			kEidmz(is, i) = MPIbuff2(k); k = k + 1;
		enddo
	enddo
	deallocate(MPIbuff1);
	deallocate(MPIbuff2);
! MPI sync
	! DBG
	! if (MPIrank .EQ. 0) then
		! write(340, '(a)') '#';
		! write(340, '(a)') 'kEidmx';
		! write(340, '(1p3E16.6E3)') kEidmx;
		! write(340, '(a)') 'kEidmy';
		! write(340, '(1p3E16.6E3)') kEidmy;
		! write(340, '(a)') 'kEidmz';
		! write(340, '(1p3E16.6E3)') kEidmz;
	! endif
	! DBG

	return;
  end subroutine kspace_efield_idm2

  subroutine kspace_efield_idm(xx,yy,zz,xxs,yys,zzs,lmax,mmax,nmax,alp,boxi)
    implicit none
    real(8) :: xx(nom),yy(nom),zz(nom)
    real(8) :: xxs(nsite,nom),yys(nsite,nom),zzs(nsite,nom)
    integer :: lmax,mmax,nmax 
    real(8) :: alp,boxi

    real(8) :: qkcpe,qvir_kc
    real(8) :: volmi 
    integer :: i,j,is,js,ii,k,li,mi,ni,mx
    integer, parameter :: nqmax = nsite*nom  
    integer :: klmmax,klmmaxsq 
!    complex(8) :: el(nqmax,0:20)
!    complex(8) :: em(nqmax,-20:20),en(nqmax,-20:20)
 !   real(8) :: qi(nqmax),xi(nqmax),yi(nqmax),zi(nqmax)
    complex(8), allocatable :: expikr(:)
    complex(8), allocatable :: el(:, :)
    complex(8), allocatable :: em(:, :), en(:, :)
    complex(8), allocatable :: elc(:, :),emc(:, :),enc(:, :)
    complex(8), allocatable :: els(:, :),ems(:, :),ens(:, :)
    real(8), allocatable :: qi(:), xi(:), yi(:), zi(:)
    real(8), allocatable :: idmxi(:), idmyi(:), idmzi(:)
    real(8), allocatable :: ckr(:),skr(:),clm(:),slm(:),ckc(:),cks(:)
    real(8) :: dk,ckcs,ckss
    real(8) :: qpe,rl,rm,rn,rksq,qforce,qefld,muefld
    complex(8) :: sumqex
    real(8), allocatable :: lkEidmx(:, :), lkEidmy(:, :), lkEidmz(:, :)
    !real(8), allocatable :: lkE0x(:, :), lkE0y(:, :), lkE0z(:, :)

    real(8), parameter :: twopi = 2.d0*pi
    integer :: noq
!    real(8) :: ak(ksqmax)
    real(8), allocatable :: ak(:)
    integer :: llim,mlim,nlim,mmin,nmin

    integer :: ll,mm,nn,kk,lmnp
        ! MPI
        INTEGER(4) MPIchunk_size_loc, iii;
        real(8), allocatable :: MPIbuff1(:), MPIbuff2(:);
        ! DBG
        real(4) iniTime, deltaTime;
        ! DBG
        ! MPI

    do i=1,nm
       do is=1,ns

       kEidmx(is,i) = 0.d0 ! initiating permanent charges fields
       kEidmy(is,i) = 0.d0
       kEidmz(is,i) = 0.d0

       end do
    end do

    allocate(el(nqmax,0:20))
    allocate(em(nqmax,-20:20))
    allocate(en(nqmax,-20:20))
    allocate(qi(nqmax))
    allocate(xi(nqmax))
    allocate(yi(nqmax))
    allocate(zi(nqmax))
    allocate(idmxi(nqmax))
    allocate(idmyi(nqmax))
    allocate(idmzi(nqmax))
    allocate(elc(nqmax,0:20))
    allocate(emc(nqmax,-20:20))
    allocate(enc(nqmax,-20:20))
    allocate(els(nqmax,0:20))
    allocate(ems(nqmax,-20:20))
    allocate(ens(nqmax,-20:20))
    allocate(ckr(nqmax))
    allocate(skr(nqmax))
    allocate(clm(nqmax))
    allocate(slm(nqmax))
    allocate(ckc(nqmax))
    allocate(cks(nqmax))

!$omp parallel DEFAULT(SHARED) private(id,i,iii,j,is,js,ii,ll,mm,nn,kk,lmnp)&
!$omp& private(klmmax,klmmaxsq,llim,mlim,nlim,volmi,mmin,nmin,expikr,ak)&
!$omp& private(lkE0x,lkE0y,lkE0z)

    allocate(expikr(nqmax))
    allocate(ak(ksqmax))
    allocate(lkEidmx(ns,nm))
    allocate(lkEidmy(ns,nm))
    allocate(lkEidmz(ns,nm))
!    allocate(lkE0x(ns,nm))
!    allocate(lkE0y(ns,nm))
!    allocate(lkE0z(ns,nm))

    klmmax = min(dble(lmax),dble(mmax),dble(nmax))
    klmmaxsq = klmmax**2

    llim = lmax 
    mlim = mmax
    nlim = nmax

    volmi = boxi**3

    do i=1,nm
       do is=1,ns
          lkEidmx(is,i) = 0.d0
          lkEidmy(is,i) = 0.d0
          lkEidmz(is,i) = 0.d0
       end do
    end do
     
!    collect quantities for charged sites

!$OMP MASTER
      i = 0
      do j = 1,nm
         do js = 1,ns
!            if (chgs(js,j).ne.0.0) then
               i = i+1
               idmxi(i) = idmx(js,j)
               idmyi(i) = idmy(js,j)
               idmzi(i) = idmz(js,j)
               xi(i) = xx(j) + xxs(js,j)
               yi(i) = yy(j) + yys(js,j)
               zi(i) = zz(j) + zzs(js,j)
!            endif
         enddo
      enddo
      noq = i
!$omp end MASTER
!$omp barrier

!$omp do
    do ii = 1,noq
       elc(ii,0) = 1.d0
       emc(ii,0) = 1.d0
       enc(ii,0) = 1.d0

       els(ii,0) = 0.d0
       ems(ii,0) = 0.d0
       ens(ii,0) = 0.d0

       elc(ii,1) = cos(twopi*xi(ii)/boxi)
       emc(ii,1) = cos(twopi*yi(ii)/boxi)
       enc(ii,1) = cos(twopi*zi(ii)/boxi)
 
       els(ii,1) = sin(twopi*xi(ii)/boxi)
       ems(ii,1) = sin(twopi*yi(ii)/boxi)
       ens(ii,1) = sin(twopi*zi(ii)/boxi)
    enddo
!$omp end do
!$omp barrier
 
!$omp MASTER

    do ll=2,llim
       do ii=1,noq

          elc(ii,ll) = elc(ii,ll-1)*elc(ii,1) - els(ii,ll-1)*els(ii,1)
          emc(ii,ll) = emc(ii,ll-1)*emc(ii,1) - ems(ii,ll-1)*ems(ii,1)
          enc(ii,ll) = enc(ii,ll-1)*enc(ii,1) - ens(ii,ll-1)*ens(ii,1)

          els(ii,ll) = els(ii,ll-1)*elc(ii,1) + elc(ii,ll-1)*els(ii,1)
          ems(ii,ll) = ems(ii,ll-1)*emc(ii,1) + emc(ii,ll-1)*ems(ii,1)
          ens(ii,ll) = ens(ii,ll-1)*enc(ii,1) + enc(ii,ll-1)*ens(ii,1)

       end do
    end do   
  
!$omp end MASTER
!$omp barrier

    mmin = 0
    nmin = 1

        ! DBG
        iniTime = SECNDS(0.0);
        ! DBG

!    do ll=0,llim
!        do ll = MPIrank, llim, MPIcomm_size
        if (MPIcomm_size .GT. 1) then
                ll = 0;
        else
                ll = -1;
        endif
        do while (ll .LE. llim)
          if (MPIrank .EQ. 0) then
                if (MPIcomm_size .GT. 1) then
                        call MPI_RECV(MPIaddr, 1, MPI_INTEGER,MPI_ANY_SOURCE,MPI_TAG_CommReady, MPI_COMM_WORLD, MPIstat, MPIierr);
                        call MPI_SEND(ll, 1, MPI_INTEGER,MPIaddr,MPI_TAG_SendData, MPI_COMM_WORLD, MPIierr);
                endif
                ll = ll + 1;
          else
                if (MPIcomm_size .GT. 1) then
                        call MPI_SEND(MPIrank, 1, MPI_INTEGER,0,MPI_TAG_CommReady, MPI_COMM_WORLD, MPIierr);
                        call MPI_RECV(ll, 1, MPI_INTEGER, 0,MPI_TAG_SendData,MPI_COMM_WORLD, MPIstat, MPIierr);
                endif
          endif
          if ((ll .LE. (llim)) .AND. ((MPIrank .GE. 1) .OR. (MPIcomm_size .EQ.1))) then
       
       rl = twopi*ll/boxi
  
       li = ll
 
       if (ll .GT. 0) mmin = -mlim
 
           !$omp do private(rl,rm,rn,rksq,qefld,sumqex) schedule(dynamic)
       do mm=mmin,mlim

          rm = twopi*mm/boxi

          mi = abs(mm)

          if (mm.ge.0) then
            do ii=1,noq

               clm(ii)=elc(ii,li)*emc(ii,mi)-els(ii,li)*ems(ii,mi)
               slm(ii)=els(ii,li)*emc(ii,mi)+ems(ii,mi)*elc(ii,li)

            end do
          else
            do ii=1,noq

               clm(ii)=elc(ii,li)*emc(ii,mi)+els(ii,li)*ems(ii,mi)
               slm(ii)=els(ii,li)*emc(ii,mi)-ems(ii,mi)*elc(ii,li)

            end do
          end if


          do nn=nmin,nlim

             ni = abs(nn)
 
             rn = twopi*nn/boxi

             kk=ll**2+mm**2+nn**2     
             lmnp = ll+mm+nn 

!             if ((kk.le.klmmaxsq) .and. .not.(mod(lmnp,2).ne.0) ) then
             if ((kk.le.klmmaxsq)) then

                rksq = rl*rl + rm*rm + rn*rn

                ak(kk) = twopi/volmi*exp(-0.25d0*rksq/alp**2)/rksq

               if (nn.ge.0) then

                   do ii=1,noq

                      !ckc(ii)=qi(ii)*(clm(ii)*enc(ii,ni)-slm(ii)*ens(ii,ni))
                      !cks(ii)=qi(ii)*(slm(ii)*enc(ii,ni)+clm(ii)*ens(ii,ni))

                     ckr(ii)=clm(ii)*enc(ii,ni)-slm(ii)*ens(ii,ni)
                     skr(ii)=slm(ii)*enc(ii,ni)+clm(ii)*ens(ii,ni)

                     ckc(ii)=ckr(ii)
                     cks(ii)=skr(ii)

                   end do

               else

                   do ii=1,noq

                    !  ckc(ii)=qi(ii)*(clm(ii)*enc(ii,ni)+slm(ii)*ens(ii,ni))
                    !  cks(ii)=qi(ii)*(slm(ii)*enc(ii,ni)-clm(ii)*ens(ii,ni))

                      ckr(ii)=clm(ii)*enc(ii,ni)+slm(ii)*ens(ii,ni)
                      skr(ii)=slm(ii)*enc(ii,ni)-clm(ii)*ens(ii,ni)

                      ckc(ii)=ckr(ii)
                      cks(ii)=skr(ii)

                   end do

               end if

              ckcs = 0.d0
              ckss = 0.d0
 
               do mx=1,noq

                  dk = rl*idmxi(mx) + rm*idmyi(mx) + rn*idmzi(mx)

                  ckcs=ckcs+dk*ckc(mx)
                  ckss=ckss+dk*cks(mx)

               end do

               ii=0
               do i=1,nm
                  do is=1,ns

                       ii=ii+1

                       muefld = -4.0*ak(kk)*(ckc(ii)*ckcs+cks(ii)*ckss)
 
!                        qefld = -4.0*ak(kk)*aimag(conjg(sumqex)*expikr(ii))
                       
                       lkEidmx(is,i) = lkEidmx(is,i) + rl*muefld
                       lkEidmy(is,i) = lkEidmy(is,i) + rm*muefld
                       lkEidmz(is,i) = lkEidmz(is,i) + rn*muefld
                
                  end do 
               end do

            end if

          end do
          nmin = -nlim

       end do
!       mmin = -mlim
           !$omp end do
          endif
    end do
        if (MPIrank .EQ. 0) then
                do i = 1, MPIcomm_size - 1
                        call MPI_RECV(MPIaddr, 1, MPI_INTEGER,MPI_ANY_SOURCE,MPI_TAG_CommReady, MPI_COMM_WORLD, MPIstat, MPIierr);
                        call MPI_SEND(llim + 1, 1, MPI_INTEGER,MPIaddr,MPI_TAG_SendData, MPI_COMM_WORLD, MPIierr);
                enddo
        endif
        ! DBG
        deltaTime = SECNDS(iniTime);
!       write(dbgUnit, *) 'P#', MPIrank, '. kspace_efield_stat calc time
!       ', deltaTime;
        ! DBG

!$omp critical

   do i=1,nm
      do is=1,ns
         kEidmx(is,i) = kEidmx(is,i) + lkEidmx(is,i)
         kEidmy(is,i) = kEidmy(is,i) + lkEidmy(is,i)
         kEidmz(is,i) = kEidmz(is,i) + lkEidmz(is,i)
      end do
   end do

!$omp end critical

    deallocate(expikr)
    deallocate(ak)
    deallocate(lkEidmx)
    deallocate(lkEidmy)
    deallocate(lkEidmz)

!$omp end parallel

    deallocate(el)
    deallocate(em)
    deallocate(en)
    deallocate(qi)
    deallocate(xi)
    deallocate(yi)
    deallocate(zi)
    deallocate(idmxi)
    deallocate(idmyi)
    deallocate(idmzi)
    deallocate(elc)
    deallocate(emc)
    deallocate(enc)
    deallocate(els)
    deallocate(ems)
    deallocate(ens)
    deallocate(ckr)
    deallocate(skr)
    deallocate(clm)
    deallocate(slm)
    deallocate(ckc)
    deallocate(cks)
! MPI sync
        j = nm * ns * 3; 
        allocate(MPIbuff1(0:j - 1));
        allocate(MPIbuff2(0:j - 1));
        do k = 0, j - 1
                MPIbuff2(k) = 0.0D0;
        enddo
        k = 0;
        do i = 1, nm
                do is = 1, ns
                        MPIbuff1(k) = kEidmx(is, i); k = k + 1;
                        MPIbuff1(k) = kEidmy(is, i); k = k + 1;
                        MPIbuff1(k) = kEidmz(is, i); k = k + 1;    
                enddo
        enddo
        call MPI_ALLREDUCE(MPIbuff1, MPIbuff2, j, MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD, MPIierr);
        k = 0;
        do i = 1, nm
                do is = 1, ns
                        kEidmx(is, i) = MPIbuff2(k); k = k + 1;
                        kEidmy(is, i) = MPIbuff2(k); k = k + 1;
                        kEidmz(is, i) = MPIbuff2(k); k = k + 1;
                enddo
        enddo
        deallocate(MPIbuff1);
        deallocate(MPIbuff2);
! MPI sync

        return;
  end subroutine kspace_efield_idm


  subroutine self_efield_idm(alp)
    implicit none
	real(8) :: alp

    integer :: is,js,i,j,k
    real(8) :: expar2
    real(8) :: drr,dist2
    real(8) :: qis,qjs
    real(8) :: ddx,ddy,ddz,d1i,d2i,d3i,d5i
    real(8) :: doti,dotj,d3ii,d5ii
    real(8) :: selffac   ! self interaction factor
	! MPI
	INTEGER(4) MPIchunk_size_loc, iii;
	real(8), allocatable :: MPIbuff1(:), MPIbuff2(:);
	real(8), allocatable :: llslfEidmx(:, :), llslfEidmy(:, :), llslfEidmz(:, :);
	! DBG
	real(4) iniTime, deltaTime;
	! DBG
	! MPI

    do i=1,nm
       do is=1,ns
		   slfEidmx(is,i) = 0.d0 ! initiating permanent charges fields
		   slfEidmy(is,i) = 0.d0
		   slfEidmz(is,i) = 0.d0
       end do
    end do

!   point dipole self interaction

    selffac = 4.d0*alp**(3.d0)/(3.d0*sqrt(pi))
	
	allocate(llslfEidmx(1:ns, 1:nm));
	allocate(llslfEidmy(1:ns, 1:nm));
	allocate(llslfEidmz(1:ns, 1:nm));
    do i=1,nm
       do is=1,ns
		   llslfEidmx(is,i) = 0.0D0;
		   llslfEidmy(is,i) = 0.0D0;
		   llslfEidmz(is,i) = 0.0D0;
       end do
    end do

!$omp parallel DEFAULT(SHARED) private(id,is,js,i,j,k)
!$omp do
    do i=1,nm
       do is=ns-npol+1,ns
          slfEidmx(is,i) = slfEidmx(is,i) + selffac*idmx(is,i)
          slfEidmy(is,i) = slfEidmy(is,i) + selffac*idmy(is,i)
          slfEidmz(is,i) = slfEidmz(is,i) + selffac*idmz(is,i)

!          slfEidmx(is,i) = slfEidmx(is,i) - selffac*idmx(is,i)
!          slfEidmy(is,i) = slfEidmy(is,i) - selffac*idmy(is,i)
!          slfEidmz(is,i) = slfEidmz(is,i) - selffac*idmz(is,i)

          !write(*,*),slfEidmx(is,i),slfEidmy(is,i),slfEidmz(is,i) 
 
       end do
    end do  
!$omp end do
!$omp barrier

! molecular contribution 
	! DBG
	iniTime = SECNDS(0.0);
	! DBG

!    do i=1,nm
    do i = MPIrank + 1, nm, MPIcomm_size
	  !$omp do private(expar2,drr,dist2,qjs,ddx,ddy,ddz,d1i,d2i,d3i,d5i,dotj,d3ii,d5ii) schedule(dynamic)
       do is=ns-npol+1,ns
          do js=ns-npol+1,ns !is+1,ns

!             if (js.eq.is) goto 15
             if (js .NE. is) then

!              if (qjs .ne. 0.d0) then

              ddx = (x(i)+xs(is,i)) - (x(i)+xs(js,i))
              ddy = (y(i)+ys(is,i)) - (y(i)+ys(js,i))
              ddz = (z(i)+zs(is,i)) - (z(i)+zs(js,i))   

!              ddx = xs(is,i) - xs(js,i) 
!              ddy = ys(is,i) - ys(js,i) 
!              ddz = zs(is,i) - zs(js,i)  

              dist2 = ddx**2+ddy**2+ddz**2

              drr = sqrt(dist2)

              d1i = 1.d0/drr
              d2i = d1i**2
              d3i = d2i*d1i
              d5i = d3i*d2i  
            
              dotj = idmx(js,i)*ddx+idmy(js,i)*ddy+idmz(js,i)*ddz

              expar2 = exp(-(alp*drr)**2.d0)

              d3ii = erf(alp*drr)*d3i - twosqpi*alp*expar2*d2i

              d5ii = erf(alp*drr)*d5i - twosqpi*alp*expar2*d2i*d2i - &
                             2.d0*twosqpi*alp**(3.d0)*expar2*d2i/3.d0

              llslfEidmx(is,i) = llslfEidmx(is,i) - (3.d0*d5ii*dotj*ddx - idmx(js,i)*d3ii)
              llslfEidmy(is,i) = llslfEidmy(is,i) - (3.d0*d5ii*dotj*ddy - idmy(js,i)*d3ii)
              llslfEidmz(is,i) = llslfEidmz(is,i) - (3.d0*d5ii*dotj*ddz - idmz(js,i)*d3ii)

!              write(*,*) llslfEidmx(is,i),llslfEidmy(is,i),llslfEidmz(is,i)

!              end if

             endif
!          15 continue
          end do

!           write(*,*) i,is,slfEidmx(is,i),slfEidmy(is,i),slfEidmz(is,i)

       end do
	   !$omp end do
    end do
!$omp end parallel
	! DBG
	deltaTime = SECNDS(iniTime);
!	write(dbgUnit, *) 'P#', MPIrank, '. self_efield_idm calc time	', deltaTime;
	! DBG
! MPI sync
	j = nm * ns * 3; !slfEidmx, slfEidmy, slfEidmz
	allocate(MPIbuff1(0:j - 1));
	allocate(MPIbuff2(0:j - 1));
	do k = 0, j - 1
		MPIbuff2(k) = 0.0D0;
	enddo
	k = 0;
	do i = 1, nm
		do is = 1, ns
			MPIbuff1(k) = llslfEidmx(is, i); k = k + 1;
			MPIbuff1(k) = llslfEidmy(is, i); k = k + 1;
			MPIbuff1(k) = llslfEidmz(is, i); k = k + 1;
		enddo
	enddo
	call MPI_ALLREDUCE(MPIbuff1, MPIbuff2, j, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, MPIierr);
	k = 0;
	do i = 1, nm
		do is = 1, ns
			slfEidmx(is, i) = slfEidmx(is, i) + MPIbuff2(k); k = k + 1;
			slfEidmy(is, i) = slfEidmy(is, i) + MPIbuff2(k); k = k + 1;
			slfEidmz(is, i) = slfEidmz(is, i) + MPIbuff2(k); k = k + 1;
		enddo
	enddo
	deallocate(MPIbuff1);
	deallocate(MPIbuff2);
! MPI sync
	deallocate(llslfEidmx);
	deallocate(llslfEidmy);
	deallocate(llslfEidmz);

	return;

!        do i = 1, nm
!           do is = 1, ns
!              write(*,*) i,is,slfEidmx(is, i),slfEidmy(is, i),slfEidmz(is, i) 
!           end do
!        end do
           

  end subroutine self_efield_idm


  subroutine self_efield_idm_v1(alp)
    implicit none
	real(8) :: alp

    integer :: is,js,i,j,k
    real(8) :: expar2
    real(8) :: drr,dist2
    real(8) :: qis,qjs
    real(8) :: ddx,ddy,ddz,d1i,d2i,d3i,d5i
    real(8) :: doti,dotj,d3ii,d5ii
    real(8) :: selffac   ! self interaction factor
	! MPI
	INTEGER(4) MPIchunk_size_loc, iii;
	real(8), allocatable :: MPIbuff1(:), MPIbuff2(:);
	real(8), allocatable :: llslfEidmx(:, :), llslfEidmy(:, :), llslfEidmz(:, :);
	! DBG
	real(4) iniTime, deltaTime;
	! DBG
	! MPI

    do i=1,nm
       do is=1,ns
		   slfEidmx(is,i) = 0.d0 ! initiating permanent charges fields
		   slfEidmy(is,i) = 0.d0
		   slfEidmz(is,i) = 0.d0
       end do
    end do

!   point dipole self interaction

    selffac = 4.d0*alp**(3.d0)/(3.d0*sqrt(pi))
	
	allocate(llslfEidmx(1:ns, 1:nm));
	allocate(llslfEidmy(1:ns, 1:nm));
	allocate(llslfEidmz(1:ns, 1:nm));
    do i=1,nm
       do is=1,ns
		   llslfEidmx(is,i) = 0.0D0;
		   llslfEidmy(is,i) = 0.0D0;
		   llslfEidmz(is,i) = 0.0D0;
       end do
    end do

!$omp parallel DEFAULT(SHARED) private(id,is,js,i,j,k)
!$omp do
    do i=1,nm
       do is=ns-npol+1,ns
          slfEidmx(is,i) = slfEidmx(is,i) + selffac*idmx(is,i)
          slfEidmy(is,i) = slfEidmy(is,i) + selffac*idmy(is,i)
          slfEidmz(is,i) = slfEidmz(is,i) + selffac*idmz(is,i) 
       end do
    end do  
!$omp end do
!$omp barrier

! molecular contribution 
	! DBG
	iniTime = SECNDS(0.0);
	! DBG

!    do i=1,nm
	if (MPIcomm_size .GT. 1) then
		i = 1;
	else
		i = 0;
	endif
	do while (i .LE. nm)
	  if (MPIrank .EQ. 0) then
		if (MPIcomm_size .GT. 1) then
			call MPI_RECV(MPIaddr, 1, MPI_INTEGER, MPI_ANY_SOURCE, MPI_TAG_CommReady, MPI_COMM_WORLD, MPIstat, MPIierr);
			call MPI_SEND(i, 1, MPI_INTEGER, MPIaddr, MPI_TAG_SendData, MPI_COMM_WORLD, MPIierr);
		endif
		i = i + 1;
	  else
		if (MPIcomm_size .GT. 1) then
			call MPI_SEND(MPIrank, 1, MPI_INTEGER, 0, MPI_TAG_CommReady, MPI_COMM_WORLD, MPIierr);
			call MPI_RECV(i, 1, MPI_INTEGER, 0, MPI_TAG_SendData, MPI_COMM_WORLD, MPIstat, MPIierr);
		endif
	  endif
 	  if ((i .LE. (nm)) .AND. ((MPIrank .GE. 1) .OR. (MPIcomm_size .EQ. 1))) then
	  !$omp do private(expar2,drr,dist2,qjs,ddx,ddy,ddz,d1i,d2i,d3i,d5i,dotj,d3ii,d5ii) schedule(dynamic)
       do is=ns-npol+1,ns
          do js=ns-npol+1,ns !is+1,ns

!             if (js.eq.is) goto 15
             if (js .NE. is) then

!              if (qjs .ne. 0.d0) then

              ddx = xs(is,i) - xs(js,i) 
              ddy = ys(is,i) - ys(js,i) 
              ddz = zs(is,i) - zs(js,i)  

              dist2 = ddx**2+ddy**2+ddz**2

              drr = sqrt(dist2)

              d1i = 1.d0/drr
              d2i = d1i**2
              d3i = d2i*d1i
              d5i = d3i*d2i  
            
              dotj = idmx(js,i)*ddx+idmy(js,i)*ddy+idmz(js,i)*ddz

              expar2 = exp(-(alp*drr)**2.d0)

              d3ii = erf(alp*drr)*d3i - twosqpi*alp*expar2*d2i

              d5ii = erf(alp*drr)*d5i - twosqpi*alp*expar2*d2i*d2i - &
                             2.d0*twosqpi*alp**(3.d0)*expar2*d2i/3.d0

              llslfEidmx(is,i) = llslfEidmx(is,i) - (3.d0*d5ii*dotj*ddx - idmx(js,i)*d3ii)
              llslfEidmy(is,i) = llslfEidmy(is,i) - (3.d0*d5ii*dotj*ddy - idmy(js,i)*d3ii)
              llslfEidmz(is,i) = llslfEidmz(is,i) - (3.d0*d5ii*dotj*ddz - idmz(js,i)*d3ii)

!              end if

             endif
!          15 continue
          end do
       end do
	   !$omp end do
	  endif
    end do
	if (MPIrank .EQ. 0) then
		do i = 1, MPIcomm_size - 1
			call MPI_RECV(MPIaddr, 1, MPI_INTEGER, MPI_ANY_SOURCE, MPI_TAG_CommReady, MPI_COMM_WORLD, MPIstat, MPIierr);
			call MPI_SEND(nm + 1, 1, MPI_INTEGER, MPIaddr, MPI_TAG_SendData, MPI_COMM_WORLD, MPIierr);
		enddo
	endif
!$omp end parallel
	! DBG
	deltaTime = SECNDS(iniTime);
!	write(dbgUnit, *) 'P#', MPIrank, '. self_efield_idm calc time	', deltaTime;
	! DBG
! MPI sync
	j = nm * ns * 3; !slfEidmx, slfEidmy, slfEidmz
	allocate(MPIbuff1(0:j - 1));
	allocate(MPIbuff2(0:j - 1));
	do k = 0, j - 1
		MPIbuff2(k) = 0.0D0;
	enddo
	k = 0;
	do i = 1, nm
		do is = 1, ns
			MPIbuff1(k) = llslfEidmx(is, i); k = k + 1;
			MPIbuff1(k) = llslfEidmy(is, i); k = k + 1;
			MPIbuff1(k) = llslfEidmz(is, i); k = k + 1;
		enddo
	enddo
	call MPI_ALLREDUCE(MPIbuff1, MPIbuff2, j, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, MPIierr);
	k = 0;
	do i = 1, nm
		do is = 1, ns
			slfEidmx(is, i) = slfEidmx(is, i) + MPIbuff2(k); k = k + 1;
			slfEidmy(is, i) = slfEidmy(is, i) + MPIbuff2(k); k = k + 1;
			slfEidmz(is, i) = slfEidmz(is, i) + MPIbuff2(k); k = k + 1;
		enddo
	enddo
	deallocate(MPIbuff1);
	deallocate(MPIbuff2);
! MPI sync
	deallocate(llslfEidmx);
	deallocate(llslfEidmy);
	deallocate(llslfEidmz);

	return;
  end subroutine self_efield_idm_v1

  subroutine surf_efield_idm(boxi)
    implicit none
    real(8) :: boxi

    integer :: is,js,i,j,k
    real(8) :: qis
    real(8) :: lsumdx, lsumdy, lsumdz
	real(8) :: volmi,fact
    real(8) :: sumdx, sumdy, sumdz

    volmi = boxi**3.d0

    fact = 4.0d0*pi/(3.0d0*volmi)

    do i=1,nm
       do is=1,ns
          srfEidmx(is,i) = 0.d0 ! initiating permanent charges fields
          srfEidmy(is,i) = 0.d0
          srfEidmz(is,i) = 0.d0
       end do
    end do

    sumdx=0.d0
    sumdy=0.d0
    sumdz=0.d0

!$omp parallel DEFAULT(SHARED) private(is,js,i,j,k)&
!$omp& private(qis,lsumdx, lsumdy, lsumdz)

    lsumdx=0.d0
    lsumdy=0.d0
    lsumdz=0.d0

!$omp do
    do i=1,nm
       do is=ns-npol+1,ns
          lsumdx = lsumdx + idmx(is,i)
          lsumdy = lsumdy + idmy(is,i)
          lsumdz = lsumdz + idmz(is,i)
       end do
    end do
!$omp end do

!$omp critical

    sumdx = sumdx + lsumdx
    sumdy = sumdy + lsumdy
    sumdz = sumdz + lsumdz

!$omp end critical

!$omp do
    do i=1,nm
       do is=ns-npol+1,ns
          srfEidmx(is,i) = srfEidmx(is,i) - fact*sumdx
          srfEidmy(is,i) = srfEidmy(is,i) - fact*sumdy
          srfEidmz(is,i) = srfEidmz(is,i) - fact*sumdz
       end do
    end do
!$omp end do
!$omp end parallel

  end subroutine surf_efield_idm


!**************************************************************************
! ************* Now calculating forces ****************    starts here
 
  subroutine real_forces(xx,yy,zz,xxs,yys,zzs,boxi,alp)
    implicit none
    real(8) :: xx(nom),yy(nom),zz(nom)
    real(8) :: xxs(nsite,nom),yys(nsite,nom),zzs(nsite,nom)
    real(8) :: boxi,alp

    real(8) :: qvirial,lqvirial

    integer :: is,js,i,j,k,isteps
    real(8) :: dist,d1i,d2i,d3i,d4i,d5i,d6i,d7i,doti,dotj

    real(8) :: ef0gxx,ef0gyy,ef0gzz
    real(8) :: ef0gxy,ef0gyx,ef0gxz
    real(8) :: ef0gzx,ef0gyz,ef0gzy

    real(8) :: efigxx,efigyy,efigzz
    real(8) :: efigxy,efigyx,efigxz
    real(8) :: efigzx,efigyz,efigzy

    real(8) :: efi0gxx,efi0gyy,efi0gzz
    real(8) :: efi0gxy,efi0gyx,efi0gxz
    real(8) :: efi0gzx,efi0gyz,efi0gzy 

!    real(8) :: refgxx(nsite,nom),refgyy(nsite,nom),refgzz(nsite,nom)
!    real(8) :: refgxy(nsite,nom),refgyx(nsite,nom),refgxz(nsite,nom)
!    real(8) :: refgzx(nsite,nom),refgyz(nsite,nom),refgzy(nsite,nom)
!    real(8) :: refxt(nsite,nom),refyt(nsite,nom),refzt(nsite,nom)
    real(8), allocatable :: refgxx(:, :), refgyy(:, :), refgzz(:, :)
    real(8), allocatable :: refgxy(:, :), refgyx(:, :), refgxz(:, :)
    real(8), allocatable :: refgzx(:, :), refgyz(:, :), refgzy(:, :)
    real(8), allocatable :: refxt(:, :), refyt(:, :), refzt(:, :)
	! Shared versions
    real(8), allocatable :: srefgxx(:, :), srefgyy(:, :), srefgzz(:, :)
    real(8), allocatable :: srefgxy(:, :), srefgyx(:, :), srefgxz(:, :)
    real(8), allocatable :: srefgzx(:, :), srefgyz(:, :), srefgzy(:, :)
    real(8), allocatable :: srefxt(:, :), srefyt(:, :), srefzt(:, :)

    real(8), allocatable :: lrefxt(:, :), lrefyt(:, :), lrefzt(:, :)
    real(8), allocatable :: llrefxt(:, :), llrefyt(:, :), llrefzt(:, :)

    real(8) :: efldx,efldy,efldz

    real(8) :: fxss,fyss,fzss
    real(8) :: qis,qjs  ! charges on sites
    real(8) :: xcc,ycc,zcc,dx,dy,dz,ddx,ddy,ddz,dist2,rccsq
    real(8) :: expar2
    real(8) :: d3ii,d5ii,d7ii,term1,term2
	! MPI
	INTEGER(4) MPIchunk_size_loc, iii;
	real(8), allocatable :: MPIbuff1(:), MPIbuff2(:);
    real(8), allocatable :: llfxs(:, :), llfys(:, :), llfzs(:, :);
	! DBG
	real(4) iniTime, deltaTime;
	! DBG
	! MPI

    allocate(srefgxx(ns,nm))
    allocate(srefgyy(ns,nm))
    allocate(srefgzz(ns,nm))
    allocate(srefgxy(ns,nm))
    allocate(srefgyx(ns,nm))
    allocate(srefgxz(ns,nm))
    allocate(srefgzx(ns,nm))
    allocate(srefgyz(ns,nm))
    allocate(srefgzy(ns,nm))
    allocate(srefxt(ns,nm))
    allocate(srefyt(ns,nm))
    allocate(srefzt(ns,nm))

    allocate(lrefxt(ns,nm))
    allocate(lrefyt(ns,nm))
    allocate(lrefzt(ns,nm))

    allocate(llrefxt(ns,nm))
    allocate(llrefyt(ns,nm))
    allocate(llrefzt(ns,nm))

    qvirial = 0.d0 

    do i=1,nm
       do is=1,ns
 
          srefgxx(is,i) = 0.d0   ! total gradient of efield
          srefgyy(is,i) = 0.d0
          srefgzz(is,i) = 0.d0
          srefgxy(is,i) = 0.d0
          srefgyx(is,i) = 0.d0
          srefgxz(is,i) = 0.d0
          srefgzx(is,i) = 0.d0
          srefgyz(is,i) = 0.d0
          srefgzy(is,i) = 0.d0

          srefxt(is,i) = 0.d0
          srefyt(is,i) = 0.d0
          srefzt(is,i) = 0.d0

          lrefxt(is,i) = 0.d0
          lrefyt(is,i) = 0.d0
          lrefzt(is,i) = 0.d0

          llrefxt(is,i) = 0.d0
          llrefyt(is,i) = 0.d0
          llrefzt(is,i) = 0.d0

       end do
    end do

    allocate(llfxs(ns,nm));
    allocate(llfys(ns,nm));
    allocate(llfzs(ns,nm));
    do i = 1, nm
       do j = 1, ns
          llfxs(j, i) = 0.0D0;
          llfys(j, i) = 0.0D0;
          llfzs(j, i) = 0.0D0;

          Etotx(j, i) = 0.0D0;
          Etoty(j, i) = 0.0D0;
          Etotz(j, i) = 0.0D0;

       enddo
    enddo
	! DBG
	! if (MPIrank .EQ. 0) then
		! write(331, '(a)') '#';
		! write(331, '(a)') 'fxs';
		! write(331, '(1p3E16.6E3)') fxs;
		! write(331, '(a)') 'fys';
		! write(331, '(1p3E16.6E3)') fys;
		! write(331, '(a)') 'fzs';
		! write(331, '(1p3E16.6E3)') fzs;
	! endif
	! DBG

!$omp parallel DEFAULT(SHARED)&
!$omp& private(is,js,i,iii,j,k,isteps,id,dist,d1i,d2i,d3i,d4i,d5i,d6i,d7i,doti,dotj)&
!$omp& private(ef0gxx,ef0gyy,ef0gzz,ef0gxy,ef0gyx,ef0gxz,ef0gzx,ef0gyz,ef0gzy)&
!$omp& private(efigxx,efigyy,efigzz,efigxy,efigyx,efigxz,efigzx,efigyz,efigzy)&
!$omp& private(efi0gxx,efi0gyy,efi0gzz,efi0gxy,efi0gyx,efi0gxz,efi0gzx,efi0gyz,efi0gzy)&
!$omp& private(efldx,efldy,efldz,fxss,fyss,fzss,qis,qjs,xcc,ycc,zcc,dx,dy,dz)&
!$omp& private(ddx,ddy,ddz,dist2,rccsq,expar2,d3ii,d5ii,d7ii,term1,term2)&
!$omp& private(refgxx,refgyy,refgzz,refgxy,refgyx,refgxz,refgzx,refgyz,refgzy)&
!$omp& private(refxt,refyt,refzt,lqvirial)

   lqvirial = 0.0D0

   allocate(refgxx(ns,nm))
   allocate(refgyy(ns,nm))
   allocate(refgzz(ns,nm))
   allocate(refgxy(ns,nm))
   allocate(refgyx(ns,nm))
   allocate(refgxz(ns,nm))
   allocate(refgzx(ns,nm))
   allocate(refgyz(ns,nm))
   allocate(refgzy(ns,nm))
   allocate(refxt(ns,nm))
   allocate(refyt(ns,nm))
   allocate(refzt(ns,nm))

    do i=1,nm
       do is=1,ns
          refgxx(is,i) = 0.d0   ! total gradient of efield
          refgyy(is,i) = 0.d0
          refgzz(is,i) = 0.d0
          refgxy(is,i) = 0.d0
          refgyx(is,i) = 0.d0
          refgxz(is,i) = 0.d0
          refgzx(is,i) = 0.d0
          refgyz(is,i) = 0.d0
          refgzy(is,i) = 0.d0

          refxt(is,i) = 0.d0
          refyt(is,i) = 0.d0
          refzt(is,i) = 0.d0
      end do
    end do
	! DBG
	iniTime = SECNDS(0.0);
	! DBG
!    do i=1,nm-1   ! loop over ith molecule starts
!	 do i = MPIrank + 1, nm - 1, MPIcomm_size
	if (MPIcomm_size .GT. 1) then
		i = 1;
	else
		i = 0;
	endif
	do while (i .LE. nm - 1)
	  if (MPIrank .EQ. 0) then
		if (MPIcomm_size .GT. 1) then
			call MPI_RECV(MPIaddr, 1, MPI_INTEGER, MPI_ANY_SOURCE, MPI_TAG_CommReady, MPI_COMM_WORLD, MPIstat, MPIierr);
			call MPI_SEND(i, 1, MPI_INTEGER, MPIaddr, MPI_TAG_SendData, MPI_COMM_WORLD, MPIierr);
		endif
		i = i + 1;
	  else
		if (MPIcomm_size .GT. 1) then
			call MPI_SEND(MPIrank, 1, MPI_INTEGER, 0, MPI_TAG_CommReady, MPI_COMM_WORLD, MPIierr);
			call MPI_RECV(i, 1, MPI_INTEGER, 0, MPI_TAG_SendData, MPI_COMM_WORLD, MPIstat, MPIierr);
		endif
	  endif
 	  if ((i .LE. (nm - 1)) .AND. ((MPIrank .GE. 1) .OR. (MPIcomm_size .EQ. 1))) then

!       do is=ns-npol+1,ns

	!$omp do schedule(dynamic)
      do j=i+1,nm   ! loop over jth molecule starts 

!          if (j.eq.i) goto 13 !cycle

          dx = xx(i)-xx(j)
          dy = yy(i)-yy(j)
          dz = zz(i)-zz(j)

          xcc = dx - boxi*nint(dx/boxi)
          ycc = dy - boxi*nint(dy/boxi)
          zcc = dz - boxi*nint(dz/boxi)

          rccsq = xcc**2 + ycc**2 + zcc**2

          if (rccsq.lt.rcutsd) then  ! com-com cutoff check starts

             do is=1,ns   ! loop over sites of molecule i  Starts 

                do js=1,ns    ! loop over sites of molecule j  Starts

                   qis = chgs(is,i)
                   qjs = chgs(js,j)

                   ddx = xcc + (xxs(is,i) - xxs(js,j))  
                   ddy = ycc + (yys(is,i) - yys(js,j))  
                   ddz = zcc + (zzs(is,i) - zzs(js,j))

                   dist2 = ddx**2 + ddy**2 + ddz**2

!                   if (dist2.lt.rcutsd) then  ! site cutoff

                   dist = sqrt(dist2)
                   d1i = 1.d0/dist
                   d2i = d1i**2.d0
                   d3i = d2i*d1i
                   d4i = d2i*d2i
                   d5i = d3i*d2i
                   d6i = d3i*d3i
                   d7i = d5i*d2i 

                   expar2 = exp(-(alp*dist)**2.d0)

                   ! dot product of induced dipole and separation vector
                   doti = idmx(is,i)*ddx+idmy(is,i)*ddy+idmz(is,i)*ddz
                   dotj = idmx(js,j)*ddx+idmy(js,j)*ddy+idmz(js,j)*ddz         

                   ! gradient of efield of permanent charges w.r.t position of charge sites

                   d3ii = erfc(alp*dist)*d3i + twosqpi*alp*expar2*d2i 

                   d5ii = erfc(alp*dist)*d5i + twosqpi*alp*expar2*d4i &
                         + 2.d0*twosqpi*alp**(3.d0)*expar2*d2i/3.d0

                   d7ii = erfc(alp*dist)*d7i + twosqpi*alp*expar2*d6i &
                         + 2.d0*twosqpi*alp**(3.d0)*expar2*d4i/3.d0 &
                         + 4.d0*twosqpi*alp**(5.d0)*expar2*d2i/15.d0 

                   efldx = 3.d0*d5ii*dotj*ddx - idmx(js,j)*d3ii
                   efldy = 3.d0*d5ii*dotj*ddy - idmy(js,j)*d3ii
                   efldz = 3.d0*d5ii*dotj*ddz - idmz(js,j)*d3ii                    

                   refxt(is,i) = refxt(is,i) + efldx
                   refyt(is,i) = refyt(is,i) + efldy
                   refzt(is,i) = refzt(is,i) + efldz

                   ! test OI
                   lrefxt(is,i) = lrefxt(is,i) + efldx + qjs*d3ii*ddx
                   lrefyt(is,i) = lrefyt(is,i) + efldy + qjs*d3ii*ddy
                   lrefzt(is,i) = lrefzt(is,i) + efldz + qjs*d3ii*ddz             
                   ! test OI
 
                   ! virial contribution ! OI  
                   lqvirial = lqvirial - qis*(efldx*ddx + efldy*ddy + efldz*ddz)        
 
                  ! gradient of efield of permanent charges w.r.t position of pol. centers

                   ef0gxx = qjs*(-3.d0*d5ii*ddx*ddx + d3ii)    
                   ef0gyy = qjs*(-3.d0*d5ii*ddy*ddy + d3ii)   
                   ef0gzz = qjs*(-3.d0*d5ii*ddz*ddz + d3ii)  
                   ef0gxy = qjs*(-3.d0*d5ii*ddx*ddy) 
                   ef0gyx = qjs*(-3.d0*d5ii*ddy*ddx) 
                   ef0gxz = qjs*(-3.d0*d5ii*ddx*ddz) 
                   ef0gzx = qjs*(-3.d0*d5ii*ddz*ddx) 
                   ef0gyz = qjs*(-3.d0*d5ii*ddy*ddz)       
                   ef0gzy = qjs*(-3.d0*d5ii*ddz*ddy) 

                  ! gradient of efield of induced dipoles w.r.t position of pol. centers 

                   efigxx = -15.d0*d7ii*ddx*ddx*dotj + 3.d0*d5ii*( dotj + 2.d0*ddx*idmx(js,j) )
                   efigyy = -15.d0*d7ii*ddy*ddy*dotj + 3.d0*d5ii*( dotj + 2.d0*ddy*idmy(js,j) )
                   efigzz = -15.d0*d7ii*ddz*ddz*dotj + 3.d0*d5ii*( dotj + 2.d0*ddz*idmz(js,j) )
                   efigxy = -15.d0*d7ii*ddx*ddy*dotj + 3.d0*d5ii*( idmx(js,j)*ddy + idmy(js,j)*ddx )
                   efigyx = -15.d0*d7ii*ddy*ddx*dotj + 3.d0*d5ii*( idmy(js,j)*ddx + idmx(js,j)*ddy )
                   efigxz = -15.d0*d7ii*ddx*ddz*dotj + 3.d0*d5ii*( idmx(js,j)*ddz + idmz(js,j)*ddx )
                   efigzx = -15.d0*d7ii*ddz*ddx*dotj + 3.d0*d5ii*( idmz(js,j)*ddx + idmx(js,j)*ddz )
                   efigyz = -15.d0*d7ii*ddy*ddz*dotj + 3.d0*d5ii*( idmy(js,j)*ddz + idmz(js,j)*ddy )
                   efigzy = -15.d0*d7ii*ddz*ddy*dotj + 3.d0*d5ii*( idmz(js,j)*ddy + idmy(js,j)*ddz )

                   refgxx(is,i) = refgxx(is,i) + ef0gxx + efigxx
                   refgyy(is,i) = refgyy(is,i) + ef0gyy + efigyy 
                   refgzz(is,i) = refgzz(is,i) + ef0gzz + efigzz 
                   refgxy(is,i) = refgxy(is,i) + ef0gxy + efigxy  
                   refgyx(is,i) = refgyx(is,i) + ef0gyx + efigyx 
                   refgxz(is,i) = refgxz(is,i) + ef0gxz + efigxz 
                   refgzx(is,i) = refgzx(is,i) + ef0gzx + efigzx  
                   refgyz(is,i) = refgyz(is,i) + ef0gyz + efigyz 
                   refgzy(is,i) = refgzy(is,i) + ef0gzy + efigzy

                   ! virial contribution  ! OI

                   efi0gxx = ef0gxx + efigxx
                   efi0gyy = ef0gyy + efigyy
                   efi0gzz = ef0gzz + efigzz
                   efi0gxy = ef0gxy + efigxy
                   efi0gyx = ef0gyx + efigyx
                   efi0gxz = ef0gxz + efigxz
                   efi0gzx = ef0gzx + efigzx
                   efi0gyz = ef0gyz + efigyz
                   efi0gzy = ef0gzy + efigzy

                   lqvirial=lqvirial - (  (idmx(is,i)*efi0gxx + idmy(is,i)*efi0gyx + idmz(is,i)*efi0gzx)*ddx &
                                        + (idmx(is,i)*efi0gxy + idmy(is,i)*efi0gyy + idmz(is,i)*efi0gzy)*ddy &
                                        + (idmx(is,i)*efi0gxz + idmy(is,i)*efi0gyz + idmz(is,i)*efi0gzz)*ddz )

!                   lqvirial=lqvirial-( idmx(is,i)*efi0gxx*ddx+idmy(is,i)*efi0gyx*ddx+idmz(is,i)*efi0gzx*ddx &
!                                      +idmx(is,i)*efi0gxy*ddy+idmy(is,i)*efi0gyy*ddy+idmz(is,i)*efi0gzy*ddy &
!                                      +idmx(is,i)*efi0gxz*ddz+idmy(is,i)*efi0gyz*ddz+idmz(is,i)*efi0gzz*ddz)

!                   lqvirial=lqvirial-( idmx(is,i)*efi0gxx*ddx+idmx(is,i)*efi0gyx*ddy+idmx(is,i)*efi0gzx*ddz &
!                                      +idmy(is,i)*efi0gxy*ddx+idmy(is,i)*efi0gyy*ddy+idmy(is,i)*efi0gzy*ddz &
!                                      +idmz(is,i)*efi0gxz*ddx+idmz(is,i)*efi0gyz*ddy+idmz(is,i)*efi0gzz*ddz)

            !       write(*,*) term1xx,term2xx,term3xx,term4xx,term5xx

                    ! gradient of efield of permanent charges w.r.t position of
                    ! charge sites

                   efldx = 3.d0*d5ii*doti*ddx - idmx(is,i)*d3ii
                   efldy = 3.d0*d5ii*doti*ddy - idmy(is,i)*d3ii
                   efldz = 3.d0*d5ii*doti*ddz - idmz(is,i)*d3ii

                   refxt(js,j) = refxt(js,j) + efldx
                   refyt(js,j) = refyt(js,j) + efldy
                   refzt(js,j) = refzt(js,j) + efldz

                   ! test OI
                   lrefxt(js,j) = lrefxt(js,j) + efldx - qis*d3ii*ddx
                   lrefyt(js,j) = lrefyt(js,j) + efldy - qis*d3ii*ddy
                   lrefzt(js,j) = lrefzt(js,j) + efldz - qis*d3ii*ddz
                   ! test OI

                   ! gradient of efield of permanent charges in ith molecule at is site 

                   ef0gxx = qis*(-3.d0*d5ii*ddx*ddx + d3ii)
                   ef0gyy = qis*(-3.d0*d5ii*ddy*ddy + d3ii)
                   ef0gzz = qis*(-3.d0*d5ii*ddz*ddz + d3ii)
                   ef0gxy = qis*(-3.d0*d5ii*ddx*ddy)
                   ef0gyx = qis*(-3.d0*d5ii*ddy*ddx)
                   ef0gxz = qis*(-3.d0*d5ii*ddx*ddz)
                   ef0gzx = qis*(-3.d0*d5ii*ddz*ddx)
                   ef0gyz = qis*(-3.d0*d5ii*ddy*ddz)
                   ef0gzy = qis*(-3.d0*d5ii*ddz*ddy)

                   ! gradient of efield of induced dipoles in ith molecule at is site

                   efigxx =  15.d0*d7ii*ddx*ddx*doti + 3.d0*d5ii*( -doti - 2.d0*ddx*idmx(is,i) ) 
                   efigyy =  15.d0*d7ii*ddy*ddy*doti + 3.d0*d5ii*( -doti - 2.d0*ddy*idmy(is,i) ) 
                   efigzz =  15.d0*d7ii*ddz*ddz*doti + 3.d0*d5ii*( -doti - 2.d0*ddz*idmz(is,i) ) 
                   efigxy =  15.d0*d7ii*ddx*ddy*doti + 3.d0*d5ii*( -idmx(is,i)*ddy - idmy(is,i)*ddx ) 
                   efigyx =  15.d0*d7ii*ddy*ddx*doti + 3.d0*d5ii*( -idmy(is,i)*ddx - idmx(is,i)*ddy ) 
                   efigxz =  15.d0*d7ii*ddx*ddz*doti + 3.d0*d5ii*( -idmx(is,i)*ddz - idmz(is,i)*ddx ) 
                   efigzx =  15.d0*d7ii*ddz*ddx*doti + 3.d0*d5ii*( -idmz(is,i)*ddx - idmx(is,i)*ddz ) 
                   efigyz =  15.d0*d7ii*ddy*ddz*doti + 3.d0*d5ii*( -idmy(is,i)*ddz - idmz(is,i)*ddy ) 
                   efigzy =  15.d0*d7ii*ddz*ddy*doti + 3.d0*d5ii*( -idmz(is,i)*ddy - idmy(is,i)*ddz ) 

                   refgxx(js,j) = refgxx(js,j) + ef0gxx + efigxx
                   refgyy(js,j) = refgyy(js,j) + ef0gyy + efigyy
                   refgzz(js,j) = refgzz(js,j) + ef0gzz + efigzz
                   refgxy(js,j) = refgxy(js,j) + ef0gxy + efigxy
                   refgyx(js,j) = refgyx(js,j) + ef0gyx + efigyx
                   refgxz(js,j) = refgxz(js,j) + ef0gxz + efigxz
                   refgzx(js,j) = refgzx(js,j) + ef0gzx + efigzx
                   refgyz(js,j) = refgyz(js,j) + ef0gyz + efigyz
                   refgzy(js,j) = refgzy(js,j) + ef0gzy + efigzy

!                   end if   ! site off ends

                end do         ! loop over sites of molecule j  ends

             end do        ! loop over sites of molecule i  ends

          end if    ! com-com cutoff check ends

!       13 continue
         end do
		!$omp end do

!      end do
		endif

    end do   ! loop over ith molecule ends   
	if (MPIrank .EQ. 0) then
		do i = 1, MPIcomm_size - 1
			call MPI_RECV(MPIaddr, 1, MPI_INTEGER, MPI_ANY_SOURCE, MPI_TAG_CommReady, MPI_COMM_WORLD, MPIstat, MPIierr);
			call MPI_SEND(nm + 1, 1, MPI_INTEGER, MPIaddr, MPI_TAG_SendData, MPI_COMM_WORLD, MPIierr);
		enddo
	endif
	! DBG
	deltaTime = SECNDS(iniTime);
!	write(dbgUnit, *) 'P#', MPIrank, '. real_forces calc time	', deltaTime;
	! DBG

!$omp critical

    qvirial = qvirial + lqvirial

    do i=1,nm
       do is=1,ns
          srefgxx(is,i) = srefgxx(is,i) + refgxx(is,i)
          srefgyy(is,i) = srefgyy(is,i) + refgyy(is,i)
          srefgzz(is,i) = srefgzz(is,i) + refgzz(is,i)
          srefgxy(is,i) = srefgxy(is,i) + refgxy(is,i)
          srefgyx(is,i) = srefgyx(is,i) + refgyx(is,i)
          srefgxz(is,i) = srefgxz(is,i) + refgxz(is,i)
          srefgzx(is,i) = srefgzx(is,i) + refgzx(is,i)
          srefgyz(is,i) = srefgyz(is,i) + refgyz(is,i)
          srefgzy(is,i) = srefgzy(is,i) + refgzy(is,i)
          srefxt(is,i) = srefxt(is,i) + refxt(is,i)
          srefyt(is,i) = srefyt(is,i) + refyt(is,i)
          srefzt(is,i) = srefzt(is,i) + refzt(is,i)

          llrefxt(is,i) = llrefxt(is,i) + lrefxt(is,i) !*r4pie0
          llrefyt(is,i) = llrefyt(is,i) + lrefyt(is,i) !*r4pie0
          llrefzt(is,i) = llrefzt(is,i) + lrefzt(is,i) !*r4pie0

       enddo
    enddo
!$omp end critical
!$omp barrier


!                               tot    tot
! ... calculate forces: f  = q e    + u    efg             !!! check this
!                        a      a      b      ab

!$omp do
    do i=1,nm
       do is=1,ns
        
         qis = chgs(is,i)

         fxss = qis*srefxt(is,i) + idmx(is,i)*srefgxx(is,i) + idmy(is,i)*srefgyx(is,i) + idmz(is,i)*srefgzx(is,i)
         fyss = qis*srefyt(is,i) + idmx(is,i)*srefgxy(is,i) + idmy(is,i)*srefgyy(is,i) + idmz(is,i)*srefgzy(is,i)
         fzss = qis*srefzt(is,i) + idmx(is,i)*srefgxz(is,i) + idmy(is,i)*srefgyz(is,i) + idmz(is,i)*srefgzz(is,i)
     
         llfxs(is,i)= llfxs(is,i) + fxss*r4pie0
         llfys(is,i)= llfys(is,i) + fyss*r4pie0
         llfzs(is,i)= llfzs(is,i) + fzss*r4pie0

!         write(*,*) 'real force ',fxss,fyss,fzss 

       end do
    end do
!$omp end do
    deallocate(refgxx)
    deallocate(refgyy)
    deallocate(refgzz)
    deallocate(refgxy)
    deallocate(refgyx)
    deallocate(refgxz)
    deallocate(refgzx)
    deallocate(refgyz)
    deallocate(refgzy)
    deallocate(refxt)
    deallocate(refyt)
    deallocate(refzt)
!$omp end parallel

!******************************************************
    deallocate(srefgxx)
    deallocate(srefgyy)
    deallocate(srefgzz)
    deallocate(srefgxy)
    deallocate(srefgyx)
    deallocate(srefgxz)
    deallocate(srefgzx)
    deallocate(srefgyz)
    deallocate(srefgzy)
    deallocate(srefxt)
    deallocate(srefyt)
    deallocate(srefzt)
! MPI sync
	j = nm * ns * 3 + nm * ns * 3 + 1; !fxs, fys, fzs, Etotx,Etoty,Etotz, qvirial
	allocate(MPIbuff1(0:j - 1));
	allocate(MPIbuff2(0:j - 1));
	do k = 0, j - 1
		MPIbuff2(k) = 0.0D0;
	enddo
	k = 0;
	do i = 1, nm
		do is = 1, ns
			MPIbuff1(k) = llfxs(is, i); k = k + 1;
			MPIbuff1(k) = llfys(is, i); k = k + 1;
			MPIbuff1(k) = llfzs(is, i); k = k + 1;

                        MPIbuff1(k) = llrefxt(is, i); k = k + 1;
                        MPIbuff1(k) = llrefyt(is, i); k = k + 1;
                        MPIbuff1(k) = llrefzt(is, i); k = k + 1;
		enddo
	enddo
        MPIbuff1(k) = qvirial; k = k + 1;   ! OI
	call MPI_ALLREDUCE(MPIbuff1, MPIbuff2, j, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, MPIierr);
	k = 0;
	do i = 1, nm
		do is = 1, ns
			fxs(is, i) = fxs(is, i) + MPIbuff2(k); k = k + 1;
			fys(is, i) = fys(is, i) + MPIbuff2(k); k = k + 1;
			fzs(is, i) = fzs(is, i) + MPIbuff2(k); k = k + 1;

                        Etotx(is, i) = Etotx(is, i) + MPIbuff2(k); k = k + 1;
                        Etoty(is, i) = Etoty(is, i) + MPIbuff2(k); k = k + 1;
                        Etotz(is, i) = Etotz(is, i) + MPIbuff2(k); k = k + 1;

		enddo
	enddo
        qvirial = MPIbuff2(k); k = k + 1;   ! OI
	deallocate(MPIbuff1);
	deallocate(MPIbuff2);
! MPI sync
    deallocate(llfxs);
    deallocate(llfys);
    deallocate(llfzs);

    deallocate(lrefxt)
    deallocate(lrefyt)
    deallocate(lrefzt)
    deallocate(llrefxt)
    deallocate(llrefyt)
    deallocate(llrefzt)

	! DBG
	! if (MPIrank .EQ. 0) then
		! write(341, '(a)') '#';
		! write(341, '(a)') 'fxs';
		! write(341, '(1p3E16.6E3)') fxs;
		! write(341, '(a)') 'fys';
		! write(341, '(1p3E16.6E3)') fys;
		! write(341, '(a)') 'fzs';
		! write(341, '(1p3E16.6E3)') fzs;
	! endif
	! DBG

        vir_rind = r4pie0*qvirial     ! OI
        vir = vir + vir_rind   ! OI

!       write(*,*) 'vir_ind_real', r4pie0*qvirial  

	return;
  end subroutine real_forces

  subroutine kspace_forces2(xx,yy,zz,xxs,yys,zzs,lmax,mmax,nmax,alp,boxi)
    implicit none
    real(8) :: xx(nom),yy(nom),zz(nom)
    real(8) :: xxs(nsite,nom),yys(nsite,nom),zzs(nsite,nom)
    integer :: lmax,mmax,nmax
    real(8) :: alp,boxi

!    integer, parameter :: nqmax = nsite*nom
    integer :: nqmax
    real(8), parameter :: twopi = 2.d0*pi
    integer :: i,j,is,js
    integer :: klmmax,klmmaxsq
    integer :: llim,mlim,nlim,mmin,nmin
    integer :: ll,mm,nn,kk,lmnp
    integer :: noq
!    complex(8) :: el(nqmax,0:20), em(nqmax,-20:20),en(nqmax,-20:20)
!    real(8) :: qi(nqmax),xi(nqmax),yi(nqmax),zi(nqmax)
!    real(8) :: dk,idmxi(nqmax),idmyi(nqmax),idmzi(nqmax)
    complex(8), allocatable :: el(:, :), em(:, :), en(:, :)
    real(8), allocatable :: qi(:), xi(:), yi(:), zi(:)
    real(8), allocatable :: idmxi(:), idmyi(:), idmzi(:)

!    complex(8) :: expikr(nqmax)
!    real(8) :: a(ksqmax)
    real(8), allocatable :: a(:)
    complex(8), allocatable :: expikr(:)

    complex(8) :: sumqex,sumidmex
    real(8) :: qkcpe,qvir_kc,qis,qjs
    real(8) :: volmi
    real(8) :: dk,qpe,rl,rm,rn,rksq,qfac,mufac,mufaci

    real(8) :: ef0gxx,ef0gyy,ef0gzz
    real(8) :: ef0gxy,ef0gyx,ef0gxz
    real(8) :: ef0gzx,ef0gyz,ef0gzy

    real(8) :: efigxx,efigyy,efigzz
    real(8) :: efigxy,efigyx,efigxz
    real(8) :: efigzx,efigyz,efigzy

    real(8) :: efi0gxx,efi0gyy,efi0gzz
    real(8) :: efi0gxy,efi0gyx,efi0gxz
    real(8) :: efi0gzx,efi0gyz,efi0gzy

    real(8) :: qvirial,lqvirial

    real(8) :: efldx,efldy,efldz
    real(8) :: fxss,fyss,fzss

!    real(8) :: kefgxx(nsite,nom),kefgyy(nsite,nom),kefgzz(nsite,nom)
!    real(8) :: kefgxy(nsite,nom),kefgyx(nsite,nom),kefgxz(nsite,nom)
!    real(8) :: kefgzx(nsite,nom),kefgyz(nsite,nom),kefgzy(nsite,nom)
!    real(8) :: kefxt(nsite,nom),kefyt(nsite,nom),kefzt(nsite,nom)
    real(8), allocatable :: kefgxx(:, :), kefgyy(:, :), kefgzz(:, :)
    real(8), allocatable :: kefgxy(:, :),kefgyx(:, :),kefgxz(:, :)
    real(8), allocatable :: kefgzx(:, :),kefgyz(:, :),kefgzy(:, :)
    real(8), allocatable :: kefxt(:, :),kefyt(:, :),kefzt(:, :)
    ! Shared
    real(8), allocatable :: skefgxx(:, :), skefgyy(:, :), skefgzz(:, :)
    real(8), allocatable :: skefgxy(:, :),skefgyx(:, :),skefgxz(:, :)
    real(8), allocatable :: skefgzx(:, :),skefgyz(:, :),skefgzy(:, :)
    real(8), allocatable :: skefxt(:, :),skefyt(:, :),skefzt(:, :)
	! MPI
	INTEGER(4) MPIchunk_size_loc, iii, k;
	real(8), allocatable :: MPIbuff1(:), MPIbuff2(:);
    real(8), allocatable :: llfxs(:, :), llfys(:, :), llfzs(:, :);
	! DBG
	real(4) iniTime, deltaTime;
	! DBG
	! MPI

    nqmax = ns*nm

    allocate(el(nqmax,0:20))
    allocate(em(nqmax,-20:20))
    allocate(en(nqmax,-20:20))
    allocate(qi(nqmax))
    allocate(xi(nqmax))
    allocate(yi(nqmax))
    allocate(zi(nqmax))
    allocate(idmxi(nqmax))
    allocate(idmyi(nqmax))
    allocate(idmzi(nqmax))
    allocate(skefgxx(ns,nm))
    allocate(skefgyy(ns,nm))
    allocate(skefgzz(ns,nm))
    allocate(skefgxy(ns,nm))
    allocate(skefgyx(ns,nm))
    allocate(skefgxz(ns,nm))
    allocate(skefgzx(ns,nm))
    allocate(skefgyz(ns,nm))
    allocate(skefgzy(ns,nm))
    allocate(skefxt(ns,nm))
    allocate(skefyt(ns,nm))
    allocate(skefzt(ns,nm))

    qvirial = 0.d0

    do i=1,nm
       do is=1,ns

          skefgxx(is,i) = 0.d0   ! total gradient of efield
          skefgyy(is,i) = 0.d0
          skefgzz(is,i) = 0.d0
          skefgxy(is,i) = 0.d0
          skefgyx(is,i) = 0.d0
          skefgxz(is,i) = 0.d0
          skefgzx(is,i) = 0.d0
          skefgyz(is,i) = 0.d0
          skefgzy(is,i) = 0.d0

          skefxt(is,i) = 0.d0
          skefyt(is,i) = 0.d0
          skefzt(is,i) = 0.d0

       end do
    end do

    allocate(llfxs(ns,nm))
    allocate(llfys(ns,nm))
    allocate(llfzs(ns,nm))
    do i=1,nm
       do is=1,ns
          llfxs(is,i) = 0.0D0
          llfys(is,i) = 0.0D0
          llfzs(is,i) = 0.0D0
       enddo
    enddo

!$omp parallel DEFAULT(SHARED)&
!$omp& private(a, expikr)&
!$omp& private(i,iii,j,is,js,klmmax,klmmaxsq,llim,mlim,nlim,mmin,nmin)&
!$omp& private(ll,mm,nn,kk,lmnp,sumqex,sumidmex,qkcpe,qvir_kc,qis,qjs)&
!$omp& private(volmi,dk,qpe,rl,rm,rn,rksq,qfac,mufac,mufaci,ef0gxx,ef0gyy,ef0gzz)&
!$omp& private(ef0gxy,ef0gyx,ef0gxz,ef0gzx,ef0gyz,ef0gzy,efigxx,efigyy,efigzz)&
!$omp& private(efigxy,efigyx,efigxz,efigzx,efigyz,efigzy,efldx,efldy,efldz)&
!$omp& private(efi0gxx,efi0gyy,efi0gzz,efi0gxy,efi0gyx,efi0gxz,efi0gzx,efi0gyz,efi0gzy,lqvirial)&
!$omp& private(fxss,fyss,fzss,kefgxx,kefgyy,kefgzz,kefgxy,kefgyx,kefgxz)&
!$omp& private(kefgzx,kefgyz,kefgzy,kefxt,kefyt,kefzt)

    lqvirial = 0.d0

    allocate(expikr(nqmax))
    allocate(a(ksqmax))
    allocate(kefgxx(ns,nm))
    allocate(kefgyy(ns,nm))
    allocate(kefgzz(ns,nm))
    allocate(kefgxy(ns,nm))
    allocate(kefgyx(ns,nm))
    allocate(kefgxz(ns,nm))
    allocate(kefgzx(ns,nm))
    allocate(kefgyz(ns,nm))
    allocate(kefgzy(ns,nm))
    allocate(kefxt(ns,nm))
    allocate(kefyt(ns,nm))
    allocate(kefzt(ns,nm))

    klmmax = min(dble(lmax),dble(mmax),dble(nmax))
    klmmaxsq = klmmax**2

    llim = lmax
    mlim = mmax
    nlim = nmax
!    qpe = 0.d0 ! GRU: Not used local variable

    volmi = boxi**3

    do i=1,nm
       do is=1,ns

          kefgxx(is,i) = 0.d0   ! total gradient of efield
          kefgyy(is,i) = 0.d0
          kefgzz(is,i) = 0.d0
          kefgxy(is,i) = 0.d0
          kefgyx(is,i) = 0.d0
          kefgxz(is,i) = 0.d0
          kefgzx(is,i) = 0.d0
          kefgyz(is,i) = 0.d0
          kefgzy(is,i) = 0.d0

          kefxt(is,i) = 0.d0
          kefyt(is,i) = 0.d0
          kefzt(is,i) = 0.d0

       end do
    end do
!    collect quantities for charged sites

!$omp MASTER
      i = 0
      do j = 1,nm
         do js = 1,ns
               i = i+1
               qi(i) = chgs(js,j)
               idmxi(i) = idmx(js,j)
               idmyi(i) = idmy(js,j)
               idmzi(i) = idmz(js,j)
               xi(i) = xx(j) + xxs(js,j)
               yi(i) = yy(j) + yys(js,j)
               zi(i) = zz(j) + zzs(js,j)
         enddo
      enddo

      noq = i
!$omp end MASTER
!$omp barrier

!$omp do
    do i = 1,noq
       el(i,0) = (1.0,0.0)
       em(i,0) = (1.0,0.0)
       en(i,0) = (1.0,0.0)
       el(i,1) = exp((0.,1.)*twopi*xi(i)/boxi)
       em(i,1) = exp((0.,1.)*twopi*yi(i)/boxi)
       en(i,1) = exp((0.,1.)*twopi*zi(i)/boxi)
    enddo
!$omp end do
!$omp barrier

!$omp MASTER
    do ll = 2,llim
         do i = 1,noq
            el(i,ll) = el(i,ll-1)*el(i,1)
         enddo
      enddo
      do mm = 2,mlim
         do i = 1,noq
            em(i,mm) = em(i,mm-1)*em(i,1)
         enddo
      enddo
      do nn = 2,nlim
         do i = 1,noq
            en(i,nn) = en(i,nn-1)*en(i,1)
         enddo
      enddo
      do mm = -mlim,-1
         do i = 1,noq
            em(i,mm) = conjg(em(i,-mm))
         enddo
      enddo
      do nn = -nlim,-1
         do i = 1,noq
            en(i,nn) = conjg(en(i,-nn))
         enddo
      enddo
!$omp end MASTER
!$omp barrier

	! DBG
	iniTime = SECNDS(0.0);
	! DBG

    mmin = 0
    nmin = 1

!    do ll=0,llim
!	 do ll = MPIrank, llim, MPIcomm_size
	if (MPIcomm_size .GT. 1) then
		ll = 0;
	else
		ll = -1;
	endif
	do while (ll .LE. llim)
	  if (MPIrank .EQ. 0) then
		if (MPIcomm_size .GT. 1) then
			call MPI_RECV(MPIaddr, 1, MPI_INTEGER, MPI_ANY_SOURCE, MPI_TAG_CommReady, MPI_COMM_WORLD, MPIstat, MPIierr);
			call MPI_SEND(ll, 1, MPI_INTEGER, MPIaddr, MPI_TAG_SendData, MPI_COMM_WORLD, MPIierr);
		endif
		ll = ll + 1;
	  else
		if (MPIcomm_size .GT. 1) then
			call MPI_SEND(MPIrank, 1, MPI_INTEGER, 0, MPI_TAG_CommReady, MPI_COMM_WORLD, MPIierr);
			call MPI_RECV(ll, 1, MPI_INTEGER, 0, MPI_TAG_SendData, MPI_COMM_WORLD, MPIstat, MPIierr);
		endif
	  endif
 	  if ((ll .LE. (llim)) .AND. ((MPIrank .GE. 1) .OR. (MPIcomm_size .EQ. 1))) then

       rl = twopi*ll/boxi
   
       if (ll .GT. 0) mmin = -mlim

	   !$omp do schedule(dynamic)
       do mm=mmin,mlim

          rm = twopi*mm/boxi

          do nn=nmin,nlim

             rn = twopi*nn/boxi

             kk=ll**2+mm**2+nn**2
             lmnp = ll+mm+nn

!             if ((kk.le.klmmaxsq) .and. .not.(mod(lmnp,2).ne.0) ) then
             if ((kk.le.klmmaxsq)) then

                rksq = rl*rl + rm*rm + rn*rn

                a(kk) = twopi/volmi*exp(-0.25d0*rksq/alp**2)/rksq

                ! form expikr for each charge
                do i = 1,noq
                   expikr(i) = el(i,ll)*em(i,mm)*en(i,nn)
                enddo

                sumqex = (0.0,0.0)
                sumidmex = (0.0,0.0)

                do i = 1,noq

                   dk = rl*idmxi(i) + rm*idmyi(i) + rn*idmzi(i)

                   sumqex = sumqex + qi(i)*expikr(i)
                   sumidmex = sumidmex + dk*expikr(i)

                enddo

               ! accumulate virials

!               lqvirial = lqvirial - 4.d0*a(kk)*(aimag(sumidmex*conjg(sumqex))+sumidmex*conjg(sumidmex)  &
!                                  +0.5d0*(aimag(sumidmex*conjg(sumqex)) - aimag(conjg(sumidmex)*sumqex)  &
!                                  +sumidmex*conjg(sumidmex))*(1.d0-0.5d0*rksq/alp**2))

               lqvirial = lqvirial - 2.d0*a(kk)*(aimag(sumidmex*conjg(sumqex))+sumidmex*conjg(sumidmex)  &
                                  +0.5d0*(sumidmex*conjg(sumidmex))*(1.d0-0.5d0*rksq/alp**2))

!               lqvirial = lqvirial - 4.d0*a(kk)*(real(sumidmex)*real(sumqex)+aimag(sumidmex*aimag(sumidmex))  &
!                                  +0.5d0*(aimag(sumidmex)**2.d0)*(1.d0-0.5d0*rksq/alp**2))

               i=0
               do j=1,nm
                  do js=1,ns
                
                      i=i+1

                      ! gradient of efield of permanent charges w.r.t position
                      ! of charge sites

                      mufac = -4.0*a(kk)*sumidmex*conjg(expikr(i))

                      efldx = rl*mufac
                      efldy = rm*mufac
                      efldz = rn*mufac

                      kefxt(js,j) = kefxt(js,j) + efldx
                      kefyt(js,j) = kefyt(js,j) + efldy
                      kefzt(js,j) = kefzt(js,j) + efldz
             
                      ! gradient of efield of permanent charges w.r.t position
                      ! of pol. centers 
                      
                      qfac = 4.0*a(kk)*sumqex*conjg(expikr(i))         
                 
                      ef0gxx = rl*rl*qfac
                      ef0gyy = rm*rm*qfac
                      ef0gzz = rn*rn*qfac
                      ef0gxy = rl*rm*qfac
                      ef0gyx = rm*rl*qfac
                      ef0gxz = rl*rn*qfac
                      ef0gzx = rn*rl*qfac
                      ef0gyz = rm*rn*qfac
                      ef0gzy = rn*rm*qfac

                      ! gradient of efield of induced dipoles w.r.t position of
                      ! pol. centers

                      mufaci = -4.0*a(kk)*aimag(sumidmex*conjg(expikr(i))) 

                      efigxx = rl*rl*mufaci
                      efigyy = rm*rm*mufaci
                      efigzz = rn*rn*mufaci
                      efigxy = rl*rm*mufaci
                      efigyx = rm*rl*mufaci
                      efigxz = rl*rn*mufaci
                      efigzx = rn*rl*mufaci
                      efigyz = rm*rn*mufaci
                      efigzy = rn*rm*mufaci

                      kefgxx(js,j) = kefgxx(js,j) + ef0gxx + efigxx
                      kefgyy(js,j) = kefgyy(js,j) + ef0gyy + efigyy
                      kefgzz(js,j) = kefgzz(js,j) + ef0gzz + efigzz
                      kefgxy(js,j) = kefgxy(js,j) + ef0gxy + efigxy
                      kefgyx(js,j) = kefgyx(js,j) + ef0gyx + efigyx
                      kefgxz(js,j) = kefgxz(js,j) + ef0gxz + efigxz
                      kefgzx(js,j) = kefgzx(js,j) + ef0gzx + efigzx
                      kefgyz(js,j) = kefgyz(js,j) + ef0gyz + efigyz
                      kefgzy(js,j) = kefgzy(js,j) + ef0gzy + efigzy

                  end do
               end do

            end if

          end do
          nmin = -nlim

       end do
!       mmin = -mlim
	   !$omp end do
	  endif
    end do
	if (MPIrank .EQ. 0) then
		do i = 1, MPIcomm_size - 1
			call MPI_RECV(MPIaddr, 1, MPI_INTEGER, MPI_ANY_SOURCE, MPI_TAG_CommReady, MPI_COMM_WORLD, MPIstat, MPIierr);
			call MPI_SEND(llim + 1, 1, MPI_INTEGER, MPIaddr, MPI_TAG_SendData, MPI_COMM_WORLD, MPIierr);
		enddo
	endif
!$omp barrier
	! DBG
	deltaTime = SECNDS(iniTime);
!	write(dbgUnit, *) 'P#', MPIrank, '. kspace_forces calc time	', deltaTime;
	! DBG

!$omp critical

   qvirial = qvirial + lqvirial

   do i=1,nm
      do is=1,ns
         skefgxx(is,i) = skefgxx(is,i) + kefgxx(is,i)
         skefgyy(is,i) = skefgyy(is,i) + kefgyy(is,i)
         skefgzz(is,i) = skefgzz(is,i) + kefgzz(is,i)
         skefgxy(is,i) = skefgxy(is,i) + kefgxy(is,i)
         skefgyx(is,i) = skefgyx(is,i) + kefgyx(is,i)
         skefgxz(is,i) = skefgxz(is,i) + kefgxz(is,i)
         skefgzx(is,i) = skefgzx(is,i) + kefgzx(is,i)
         skefgyz(is,i) = skefgyz(is,i) + kefgyz(is,i)
         skefgzy(is,i) = skefgzy(is,i) + kefgzy(is,i)
         skefxt(is,i)  = skefxt(is,i)  + kefxt(is,i)
         skefyt(is,i)  = skefyt(is,i)  + kefyt(is,i)
         skefzt(is,i)  = skefzt(is,i)  + kefzt(is,i)
      end do
   end do
!$omp end critical
!$omp barrier

!                               tot    tot
! ... calculate forces: f  = q e    + u    efg             !!! check this
!                        a      a      b      ab

!$omp do
    do i=1,nm
       do is=1,ns

         qis = chgs(is,i)

         fxss = qis*skefxt(is,i) + idmx(is,i)*skefgxx(is,i) + idmy(is,i)*skefgyx(is,i) + idmz(is,i)*skefgzx(is,i)
         fyss = qis*skefyt(is,i) + idmx(is,i)*skefgxy(is,i) + idmy(is,i)*skefgyy(is,i) + idmz(is,i)*skefgzy(is,i)
         fzss = qis*skefzt(is,i) + idmx(is,i)*skefgxz(is,i) + idmy(is,i)*skefgyz(is,i) + idmz(is,i)*skefgzz(is,i)

         llfxs(is,i)= llfxs(is,i) + fxss*r4pie0
         llfys(is,i)= llfys(is,i) + fyss*r4pie0
         llfzs(is,i)= llfzs(is,i) + fzss*r4pie0

       end do
    end do
!$omp end do
    deallocate(expikr)
    deallocate(a)
    deallocate(kefgxx)
    deallocate(kefgyy)
    deallocate(kefgzz)
    deallocate(kefgxy)
    deallocate(kefgyx)
    deallocate(kefgxz)
    deallocate(kefgzx)
    deallocate(kefgyz)
    deallocate(kefgzy)
    deallocate(kefxt)
    deallocate(kefyt)
    deallocate(kefzt)
!$omp end parallel

    deallocate(el)
    deallocate(em)
    deallocate(en)
    deallocate(qi)
    deallocate(xi)
    deallocate(yi)
    deallocate(zi)
    deallocate(idmxi)
    deallocate(idmyi)
    deallocate(idmzi)
    deallocate(skefgxx)
    deallocate(skefgyy)
    deallocate(skefgzz)
    deallocate(skefgxy)
    deallocate(skefgyx)
    deallocate(skefgxz)
    deallocate(skefgzx)
    deallocate(skefgyz)
    deallocate(skefgzy)
    deallocate(skefxt)
    deallocate(skefyt)
    deallocate(skefzt)
! MPI sync
	j = nm * ns * 3 + 1 ; !fxs, fys, fzs, qvirial
	allocate(MPIbuff1(0:j - 1));
	allocate(MPIbuff2(0:j - 1));
	do k = 0, j - 1
		MPIbuff2(k) = 0.0D0;
	enddo
	k = 0;
	do i = 1, nm
		do is = 1, ns
			MPIbuff1(k) = llfxs(is, i); k = k + 1;
			MPIbuff1(k) = llfys(is, i); k = k + 1;
			MPIbuff1(k) = llfzs(is, i); k = k + 1;
		enddo
	enddo
        MPIbuff1(k) = qvirial; k = k + 1;   ! OI 
	call MPI_ALLREDUCE(MPIbuff1, MPIbuff2, j, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, MPIierr);
	k = 0;
	do i = 1, nm
		do is = 1, ns
			fxs(is, i) = fxs(is, i) + MPIbuff2(k); k = k + 1;
			fys(is, i) = fys(is, i) + MPIbuff2(k); k = k + 1;
			fzs(is, i) = fzs(is, i) + MPIbuff2(k); k = k + 1;
		enddo
	enddo
        qvirial = MPIbuff2(k); k = k + 1;   ! OI
	deallocate(MPIbuff1);
	deallocate(MPIbuff2);
! MPI sync
    deallocate(llfxs);
    deallocate(llfys);
    deallocate(llfzs);


    vir_kind = 2.d0*r4pie0*qvirial
    vir = vir + vir_kind

!    write(*,*) 'vir_ind_kspace', vir_kind

	return;
  end subroutine kspace_forces2


  subroutine kspace_forces(xx,yy,zz,xxs,yys,zzs,lmax,mmax,nmax,alp,boxi)
    implicit none
    real(8) :: xx(nom),yy(nom),zz(nom)
    real(8) :: xxs(nsite,nom),yys(nsite,nom),zzs(nsite,nom)
    integer :: lmax,mmax,nmax
    real(8) :: alp,boxi

!    integer, parameter :: nqmax = nsite*nom
    integer :: nqmax
    real(8), parameter :: twopi = 2.d0*pi
    integer :: i,j,is,js,ii,li,mi,ni,mx
    integer :: klmmax,klmmaxsq
    integer :: llim,mlim,nlim,mmin,nmin
    integer :: ll,mm,nn,kk,lmnp
    integer :: noq
!    complex(8) :: el(nqmax,0:20), em(nqmax,-20:20),en(nqmax,-20:20)
!    real(8) :: qi(nqmax),xi(nqmax),yi(nqmax),zi(nqmax)
!    real(8) :: dk,idmxi(nqmax),idmyi(nqmax),idmzi(nqmax)
    complex(8), allocatable :: el(:, :), em(:, :), en(:, :)
    real(8), allocatable :: qi(:), xi(:), yi(:), zi(:)
    real(8), allocatable :: idmxi(:), idmyi(:), idmzi(:)
    complex(8), allocatable :: elc(:, :),emc(:, :),enc(:, :)
    complex(8), allocatable :: els(:, :),ems(:, :),ens(:, :)
    real(8), allocatable :: ckr(:),skr(:),clm(:),slm(:),ckc(:),cks(:)
!    complex(8) :: expikr(nqmax)
!    real(8) :: a(ksqmax)
    real(8), allocatable :: a(:)
    complex(8), allocatable :: expikr(:)

    complex(8) :: sumqex,sumidmex
    real(8) :: qkcpe,qvir_kc,qis,qjs
    real(8) :: volmi,ckcqs,cksqs,ckcds,cksds
    real(8) :: dk,qpe,rl,rm,rn,rksq,qfac,mufac,mufaci

    real(8) :: ef0gxx,ef0gyy,ef0gzz
    real(8) :: ef0gxy,ef0gyx,ef0gxz
    real(8) :: ef0gzx,ef0gyz,ef0gzy

    real(8) :: efigxx,efigyy,efigzz
    real(8) :: efigxy,efigyx,efigxz
    real(8) :: efigzx,efigyz,efigzy

    real(8) :: efi0gxx,efi0gyy,efi0gzz
    real(8) :: efi0gxy,efi0gyx,efi0gxz
    real(8) :: efi0gzx,efi0gyz,efi0gzy

    real(8) :: qvirial,lqvirial

    real(8) :: efldx,efldy,efldz
    real(8) :: fxss,fyss,fzss

!    real(8) :: kefgxx(nsite,nom),kefgyy(nsite,nom),kefgzz(nsite,nom)
!    real(8) :: kefgxy(nsite,nom),kefgyx(nsite,nom),kefgxz(nsite,nom)
!    real(8) :: kefgzx(nsite,nom),kefgyz(nsite,nom),kefgzy(nsite,nom)
!    real(8) :: kefxt(nsite,nom),kefyt(nsite,nom),kefzt(nsite,nom)
    real(8), allocatable :: kefgxx(:, :), kefgyy(:, :), kefgzz(:, :)
    real(8), allocatable :: kefgxy(:, :),kefgyx(:, :),kefgxz(:, :)
    real(8), allocatable :: kefgzx(:, :),kefgyz(:, :),kefgzy(:, :)
    real(8), allocatable :: kefxt(:, :),kefyt(:, :),kefzt(:, :)
    ! Shared
    real(8), allocatable :: skefgxx(:, :), skefgyy(:, :), skefgzz(:, :)
    real(8), allocatable :: skefgxy(:, :),skefgyx(:, :),skefgxz(:, :)
    real(8), allocatable :: skefgzx(:, :),skefgyz(:, :),skefgzy(:, :)
    real(8), allocatable :: skefxt(:, :),skefyt(:, :),skefzt(:, :)

    real(8), allocatable :: lkefxt(:, :),lkefyt(:, :),lkefzt(:, :)
    real(8), allocatable :: llkefxt(:, :),llkefyt(:, :),llkefzt(:, :)
        ! MPI
        INTEGER(4) MPIchunk_size_loc, iii, k;
        real(8), allocatable :: MPIbuff1(:), MPIbuff2(:);
    real(8), allocatable :: llfxs(:, :), llfys(:, :), llfzs(:, :);
        ! DBG
        real(4) iniTime, deltaTime;
        ! DBG
        ! MPI

    nqmax = ns*nm

    allocate(el(nqmax,0:20))
    allocate(em(nqmax,-20:20))
    allocate(en(nqmax,-20:20))
    allocate(qi(nqmax))
    allocate(xi(nqmax))
    allocate(yi(nqmax))
    allocate(zi(nqmax))
    allocate(idmxi(nqmax))
    allocate(idmyi(nqmax))
    allocate(idmzi(nqmax))

    allocate(elc(nqmax,0:20))
    allocate(emc(nqmax,-20:20))
    allocate(enc(nqmax,-20:20))
    allocate(els(nqmax,0:20))
    allocate(ems(nqmax,-20:20))
    allocate(ens(nqmax,-20:20))
    allocate(ckr(nqmax))
    allocate(skr(nqmax))
    allocate(clm(nqmax))
    allocate(slm(nqmax))
    allocate(ckc(nqmax))
    allocate(cks(nqmax))

    allocate(skefgxx(ns,nm))
    allocate(skefgyy(ns,nm))
    allocate(skefgzz(ns,nm))
    allocate(skefgxy(ns,nm))
    allocate(skefgyx(ns,nm))
    allocate(skefgxz(ns,nm))
    allocate(skefgzx(ns,nm))
    allocate(skefgyz(ns,nm))
    allocate(skefgzy(ns,nm))
    allocate(skefxt(ns,nm))
    allocate(skefyt(ns,nm))
    allocate(skefzt(ns,nm))

    allocate(lkefxt(ns,nm))
    allocate(lkefyt(ns,nm))
    allocate(lkefzt(ns,nm))
    allocate(llkefxt(ns,nm))
    allocate(llkefyt(ns,nm))
    allocate(llkefzt(ns,nm))

    qvirial = 0.d0

    do i=1,nm
       do is=1,ns

          skefgxx(is,i) = 0.d0   ! total gradient of efield
          skefgyy(is,i) = 0.d0
          skefgzz(is,i) = 0.d0
          skefgxy(is,i) = 0.d0
          skefgyx(is,i) = 0.d0
          skefgxz(is,i) = 0.d0
          skefgzx(is,i) = 0.d0
          skefgyz(is,i) = 0.d0
          skefgzy(is,i) = 0.d0

          skefxt(is,i) = 0.d0
          skefyt(is,i) = 0.d0
          skefzt(is,i) = 0.d0

          lkefxt(is,i) = 0.d0
          lkefyt(is,i) = 0.d0
          lkefzt(is,i) = 0.d0

          llkefxt(is,i) = 0.d0
          llkefyt(is,i) = 0.d0
          llkefzt(is,i) = 0.d0

       end do
    end do

    allocate(llfxs(ns,nm))
    allocate(llfys(ns,nm))
    allocate(llfzs(ns,nm))
    do i=1,nm
       do is=1,ns
          llfxs(is,i) = 0.0D0
          llfys(is,i) = 0.0D0
          llfzs(is,i) = 0.0D0
       enddo
    enddo

!$omp parallel DEFAULT(SHARED)&
!$omp& private(a, expikr)&
!$omp& private(i,iii,j,is,js,klmmax,klmmaxsq,llim,mlim,nlim,mmin,nmin)&
!$omp& private(ll,mm,nn,kk,lmnp,sumqex,sumidmex,qkcpe,qvir_kc,qis,qjs)&
!$omp&
!private(volmi,dk,qpe,rl,rm,rn,rksq,qfac,mufac,mufaci,ef0gxx,ef0gyy,ef0gzz)&
!$omp& private(ef0gxy,ef0gyx,ef0gxz,ef0gzx,ef0gyz,ef0gzy,efigxx,efigyy,efigzz)&
!$omp& private(efigxy,efigyx,efigxz,efigzx,efigyz,efigzy,efldx,efldy,efldz)&
!$omp&
!private(efi0gxx,efi0gyy,efi0gzz,efi0gxy,efi0gyx,efi0gxz,efi0gzx,efi0gyz,efi0gzy,lqvirial)&
!$omp& private(fxss,fyss,fzss,kefgxx,kefgyy,kefgzz,kefgxy,kefgyx,kefgxz)&
!$omp& private(kefgzx,kefgyz,kefgzy,kefxt,kefyt,kefzt)

    lqvirial = 0.d0

    allocate(expikr(nqmax))
    allocate(a(ksqmax))
    allocate(kefgxx(ns,nm))
    allocate(kefgyy(ns,nm))
    allocate(kefgzz(ns,nm))
    allocate(kefgxy(ns,nm))
    allocate(kefgyx(ns,nm))
    allocate(kefgxz(ns,nm))
    allocate(kefgzx(ns,nm))
    allocate(kefgyz(ns,nm))
    allocate(kefgzy(ns,nm))
    allocate(kefxt(ns,nm))
    allocate(kefyt(ns,nm))
    allocate(kefzt(ns,nm))

    klmmax = min(dble(lmax),dble(mmax),dble(nmax))
    klmmaxsq = klmmax**2

    llim = lmax
    mlim = mmax
    nlim = nmax
!    qpe = 0.d0 ! GRU: Not used local variable

    volmi = boxi**3

    do i=1,nm
       do is=1,ns

          kefgxx(is,i) = 0.d0   ! total gradient of efield
          kefgyy(is,i) = 0.d0
          kefgzz(is,i) = 0.d0
          kefgxy(is,i) = 0.d0
          kefgyx(is,i) = 0.d0
          kefgxz(is,i) = 0.d0
          kefgzx(is,i) = 0.d0
          kefgyz(is,i) = 0.d0
          kefgzy(is,i) = 0.d0

          kefxt(is,i) = 0.d0
          kefyt(is,i) = 0.d0
          kefzt(is,i) = 0.d0

       end do
    end do
!    collect quantities for charged sites

!$omp MASTER
      i = 0
      do j = 1,nm
         do js = 1,ns
               i = i+1
               qi(i) = chgs(js,j)
               idmxi(i) = idmx(js,j)
               idmyi(i) = idmy(js,j)
               idmzi(i) = idmz(js,j)
               xi(i) = xx(j) + xxs(js,j)
               yi(i) = yy(j) + yys(js,j)
               zi(i) = zz(j) + zzs(js,j)
         enddo
      enddo

      noq = i
!$omp end MASTER
!$omp barrier

!$omp do
    do ii = 1,noq
       elc(ii,0) = 1.d0
       emc(ii,0) = 1.d0
       enc(ii,0) = 1.d0

       els(ii,0) = 0.d0
       ems(ii,0) = 0.d0
       ens(ii,0) = 0.d0

       elc(ii,1) = cos(twopi*xi(ii)/boxi)
       emc(ii,1) = cos(twopi*yi(ii)/boxi)
       enc(ii,1) = cos(twopi*zi(ii)/boxi)
 
       els(ii,1) = sin(twopi*xi(ii)/boxi)
       ems(ii,1) = sin(twopi*yi(ii)/boxi)
       ens(ii,1) = sin(twopi*zi(ii)/boxi)
    enddo

!$omp end do
!$omp barrier

!$omp MASTER

    do ll=2,llim
       do ii=1,noq

          elc(ii,ll) = elc(ii,ll-1)*elc(ii,1) - els(ii,ll-1)*els(ii,1)
          emc(ii,ll) = emc(ii,ll-1)*emc(ii,1) - ems(ii,ll-1)*ems(ii,1)
          enc(ii,ll) = enc(ii,ll-1)*enc(ii,1) - ens(ii,ll-1)*ens(ii,1)

          els(ii,ll) = els(ii,ll-1)*elc(ii,1) + elc(ii,ll-1)*els(ii,1)
          ems(ii,ll) = ems(ii,ll-1)*emc(ii,1) + emc(ii,ll-1)*ems(ii,1)
          ens(ii,ll) = ens(ii,ll-1)*enc(ii,1) + enc(ii,ll-1)*ens(ii,1)

       end do
    end do

!$omp end MASTER
!$omp barrier

        ! DBG
        iniTime = SECNDS(0.0);
        ! DBG

    mmin = 0
    nmin = 1

!    do ll=0,llim
!        do ll = MPIrank, llim, MPIcomm_size
        if (MPIcomm_size .GT. 1) then
                ll = 0;
        else
                ll = -1;
        endif
        do while (ll .LE. llim)
          if (MPIrank .EQ. 0) then
                if (MPIcomm_size .GT. 1) then
                        call MPI_RECV(MPIaddr, 1, MPI_INTEGER, MPI_ANY_SOURCE,MPI_TAG_CommReady, MPI_COMM_WORLD, MPIstat, MPIierr);
                        call MPI_SEND(ll, 1, MPI_INTEGER, MPIaddr,MPI_TAG_SendData, MPI_COMM_WORLD, MPIierr);
                endif
                ll = ll + 1;
          else
                if (MPIcomm_size .GT. 1) then
                        call MPI_SEND(MPIrank, 1, MPI_INTEGER, 0, MPI_TAG_CommReady, MPI_COMM_WORLD, MPIierr);
                        call MPI_RECV(ll, 1, MPI_INTEGER, 0, MPI_TAG_SendData, MPI_COMM_WORLD, MPIstat, MPIierr);
                endif
          endif
          if ((ll .LE. (llim)) .AND. ((MPIrank .GE. 1) .OR. (MPIcomm_size .EQ. 1))) then

       rl = twopi*ll/boxi
  
       li = ll
 
       if (ll .GT. 0) mmin = -mlim

           !$omp do schedule(dynamic)
       do mm=mmin,mlim

          rm = twopi*mm/boxi

          mi = abs(mm)

          if (mm.ge.0) then
            do ii=1,noq

               clm(ii)=elc(ii,li)*emc(ii,mi)-els(ii,li)*ems(ii,mi)
               slm(ii)=els(ii,li)*emc(ii,mi)+ems(ii,mi)*elc(ii,li)

            end do
          else
            do ii=1,noq

               clm(ii)=elc(ii,li)*emc(ii,mi)+els(ii,li)*ems(ii,mi)
               slm(ii)=els(ii,li)*emc(ii,mi)-ems(ii,mi)*elc(ii,li)

            end do
          end if 

          do nn=nmin,nlim

             ni = abs(nn)
             rn = twopi*nn/boxi

             kk=ll**2+mm**2+nn**2
             lmnp = ll+mm+nn

!             if ((kk.le.klmmaxsq) .and. .not.(mod(lmnp,2).ne.0) ) then
             if ((kk.le.klmmaxsq)) then

                rksq = rl*rl + rm*rm + rn*rn

                a(kk) = twopi/volmi*exp(-0.25d0*rksq/alp**2)/rksq

                if (nn.ge.0) then

                   do ii=1,noq

                      !ckc(ii)=qi(ii)*(clm(ii)*enc(ii,ni)-slm(ii)*ens(ii,ni))
                      !cks(ii)=qi(ii)*(slm(ii)*enc(ii,ni)+clm(ii)*ens(ii,ni))

                     ckr(ii)=clm(ii)*enc(ii,ni)-slm(ii)*ens(ii,ni)
                     skr(ii)=slm(ii)*enc(ii,ni)+clm(ii)*ens(ii,ni)

                     ckc(ii)=ckr(ii)
                     cks(ii)=skr(ii)

                   end do

                else

                   do ii=1,noq

                    !  ckc(ii)=qi(ii)*(clm(ii)*enc(ii,ni)+slm(ii)*ens(ii,ni))
                    !  cks(ii)=qi(ii)*(slm(ii)*enc(ii,ni)-clm(ii)*ens(ii,ni))

                      ckr(ii)=clm(ii)*enc(ii,ni)+slm(ii)*ens(ii,ni)
                      skr(ii)=slm(ii)*enc(ii,ni)-clm(ii)*ens(ii,ni)

                      ckc(ii)=ckr(ii)
                      cks(ii)=skr(ii)

                   end do

                end if

                ckcqs = 0.d0
                cksqs = 0.d0
                ckcds = 0.d0
                cksds = 0.d0

                do mx=1,noq

                  dk = rl*idmxi(mx) + rm*idmyi(mx) + rn*idmzi(mx)

                  ckcqs=ckcqs+qi(mx)*ckc(mx)
                  cksqs=cksqs+qi(mx)*cks(mx) 
                  ckcds=ckcds+dk*ckc(mx)
                  cksds=cksds+dk*cks(mx)

               end do

               ! accumulate virials

!               lqvirial = lqvirial -
!               4.d0*a(kk)*(aimag(sumidmex*conjg(sumqex))+sumidmex*conjg(sumidmex)
!               &
!                                  +0.5d0*(aimag(sumidmex*conjg(sumqex)) -
!                                  aimag(conjg(sumidmex)*sumqex)  &
!                                  +sumidmex*conjg(sumidmex))*(1.d0-0.5d0*rksq/alp**2))

               lqvirial = lqvirial -2.d0*a(kk)*((ckcds*cksqs-cksds*ckcqs)+(ckcds*ckcds+cksds*cksds)  &
                                  +0.5d0*(2.d0*(ckcds*cksqs-cksds*ckcqs)+(ckcds*ckcds+cksds*cksds))  &
                                          *(1.d0-0.5d0*rksq/alp**2))

!               lqvirial = lqvirial -
!               4.d0*a(kk)*(real(sumidmex)*real(sumqex)+aimag(sumidmex*aimag(sumidmex))
!               &
!                                  +0.5d0*(aimag(sumidmex)**2.d0)*(1.d0-0.5d0*rksq/alp**2))

               i=0
               do j=1,nm
                  do js=1,ns
                
                      i=i+1

                      ! gradient of efield of permanent charges w.r.t position
                      ! of charge sites

                      mufac = -4.0*a(kk)*(ckc(i)*ckcds+cks(i)*cksds)

                      efldx = rl*mufac
                      efldy = rm*mufac
                      efldz = rn*mufac

                      kefxt(js,j) = kefxt(js,j) + efldx
                      kefyt(js,j) = kefyt(js,j) + efldy
                      kefzt(js,j) = kefzt(js,j) + efldz

                      lkefxt(js,j) = lkefxt(js,j) + efldx-4.0*a(kk)*rl*(ckc(i)*cksqs-cks(i)*ckcqs) 
                      lkefyt(js,j) = lkefyt(js,j) + efldy-4.0*a(kk)*rm*(ckc(i)*cksqs-cks(i)*ckcqs)
                      lkefzt(js,j) = lkefzt(js,j) + efldz-4.0*a(kk)*rn*(ckc(i)*cksqs-cks(i)*ckcqs)
             
                      ! gradient of efield of permanent charges w.r.t position
                      ! of pol. centers 
                      
                      qfac = -4.0*a(kk)*(-ckc(i)*ckcqs-cks(i)*cksqs) 
                 
                      ef0gxx = rl*rl*qfac
                      ef0gyy = rm*rm*qfac
                      ef0gzz = rn*rn*qfac
                      ef0gxy = rl*rm*qfac
                      ef0gyx = rm*rl*qfac
                      ef0gxz = rl*rn*qfac
                      ef0gzx = rn*rl*qfac
                      ef0gyz = rm*rn*qfac
                      ef0gzy = rn*rm*qfac

                      ! gradient of efield of induced dipoles w.r.t position of
                      ! pol. centers

                      mufaci = -4.0*a(kk)*(ckc(i)*cksds-cks(i)*ckcds)

                      efigxx = rl*rl*mufaci
                      efigyy = rm*rm*mufaci
                      efigzz = rn*rn*mufaci
                      efigxy = rl*rm*mufaci
                      efigyx = rm*rl*mufaci
                      efigxz = rl*rn*mufaci
                      efigzx = rn*rl*mufaci
                      efigyz = rm*rn*mufaci
                      efigzy = rn*rm*mufaci

                      kefgxx(js,j) = kefgxx(js,j) + ef0gxx + efigxx
                      kefgyy(js,j) = kefgyy(js,j) + ef0gyy + efigyy
                      kefgzz(js,j) = kefgzz(js,j) + ef0gzz + efigzz
                      kefgxy(js,j) = kefgxy(js,j) + ef0gxy + efigxy
                      kefgyx(js,j) = kefgyx(js,j) + ef0gyx + efigyx
                      kefgxz(js,j) = kefgxz(js,j) + ef0gxz + efigxz
                      kefgzx(js,j) = kefgzx(js,j) + ef0gzx + efigzx
                      kefgyz(js,j) = kefgyz(js,j) + ef0gyz + efigyz
                      kefgzy(js,j) = kefgzy(js,j) + ef0gzy + efigzy

                  end do
               end do

            end if

          end do
          nmin = -nlim

       end do
!       mmin = -mlim
           !$omp end do
          endif
    end do
        if (MPIrank .EQ. 0) then
                do i = 1, MPIcomm_size - 1
                        call MPI_RECV(MPIaddr, 1, MPI_INTEGER, MPI_ANY_SOURCE, MPI_TAG_CommReady, MPI_COMM_WORLD, MPIstat, MPIierr);
                        call MPI_SEND(llim + 1, 1, MPI_INTEGER, MPIaddr, MPI_TAG_SendData, MPI_COMM_WORLD, MPIierr);
                enddo
        endif
!$omp barrier
        ! DBG
        deltaTime = SECNDS(iniTime);
!       write(dbgUnit, *) 'P#', MPIrank, '. kspace_forces calc time     ',
!       deltaTime;
        ! DBG

!$omp critical

   qvirial = qvirial + lqvirial

   do i=1,nm
      do is=1,ns
         skefgxx(is,i) = skefgxx(is,i) + kefgxx(is,i)
         skefgyy(is,i) = skefgyy(is,i) + kefgyy(is,i)
         skefgzz(is,i) = skefgzz(is,i) + kefgzz(is,i)
         skefgxy(is,i) = skefgxy(is,i) + kefgxy(is,i)
         skefgyx(is,i) = skefgyx(is,i) + kefgyx(is,i)
         skefgxz(is,i) = skefgxz(is,i) + kefgxz(is,i)
         skefgzx(is,i) = skefgzx(is,i) + kefgzx(is,i)
         skefgyz(is,i) = skefgyz(is,i) + kefgyz(is,i)
         skefgzy(is,i) = skefgzy(is,i) + kefgzy(is,i)
         skefxt(is,i)  = skefxt(is,i)  + kefxt(is,i)
         skefyt(is,i)  = skefyt(is,i)  + kefyt(is,i)
         skefzt(is,i)  = skefzt(is,i)  + kefzt(is,i)

         llkefxt(is,i)  = llkefxt(is,i)  + lkefxt(is,i)!*r4pie0
         llkefyt(is,i)  = llkefyt(is,i)  + lkefyt(is,i)!*r4pie0
         llkefzt(is,i)  = llkefzt(is,i)  + lkefzt(is,i)!*r4pie0

      end do
   end do
!$omp end critical
!$omp barrier

!                               tot    tot
! ... calculate forces: f  = q e    + u    efg             !!! check this
!                        a      a      b      ab

!$omp do
    do i=1,nm
       do is=1,ns

         qis = chgs(is,i)

         fxss = qis*skefxt(is,i) + idmx(is,i)*skefgxx(is,i) + idmy(is,i)*skefgyx(is,i) + idmz(is,i)*skefgzx(is,i)
         fyss = qis*skefyt(is,i) + idmx(is,i)*skefgxy(is,i) + idmy(is,i)*skefgyy(is,i) + idmz(is,i)*skefgzy(is,i)
         fzss = qis*skefzt(is,i) + idmx(is,i)*skefgxz(is,i) + idmy(is,i)*skefgyz(is,i) + idmz(is,i)*skefgzz(is,i)

         llfxs(is,i)= llfxs(is,i) + fxss*r4pie0
         llfys(is,i)= llfys(is,i) + fyss*r4pie0
         llfzs(is,i)= llfzs(is,i) + fzss*r4pie0

       end do
    end do
!$omp end do
    deallocate(expikr)
    deallocate(a)
    deallocate(kefgxx)
    deallocate(kefgyy)
    deallocate(kefgzz)
    deallocate(kefgxy)
    deallocate(kefgyx)
    deallocate(kefgxz)
    deallocate(kefgzx)
    deallocate(kefgyz)
    deallocate(kefgzy)
    deallocate(kefxt)
    deallocate(kefyt)
    deallocate(kefzt)
!$omp end parallel

    deallocate(el)
    deallocate(em)
    deallocate(en)
    deallocate(qi)
    deallocate(xi)
    deallocate(yi)
    deallocate(zi)
    deallocate(idmxi)
    deallocate(idmyi)
    deallocate(idmzi)
    deallocate(skefgxx)
    deallocate(skefgyy)
    deallocate(skefgzz)
    deallocate(skefgxy)
    deallocate(skefgyx)
    deallocate(skefgxz)
    deallocate(skefgzx)
    deallocate(skefgyz)
    deallocate(skefgzy)
    deallocate(skefxt)
    deallocate(skefyt)
    deallocate(skefzt)
    deallocate(elc)
    deallocate(emc)
    deallocate(enc)
    deallocate(els)
    deallocate(ems)
    deallocate(ens)
    deallocate(ckr)
    deallocate(skr)
    deallocate(clm)
    deallocate(slm)
    deallocate(ckc)
    deallocate(cks)

! MPI sync
        j = nm * ns * 3 + nm * ns * 3 + 1 ; !fxs, fys, fzs, Etotx,Etoty,Etotz, qvirial
        allocate(MPIbuff1(0:j - 1));
        allocate(MPIbuff2(0:j - 1));
        do k = 0, j - 1
                MPIbuff2(k) = 0.0D0;
        enddo
        k = 0;
        do i = 1, nm
                do is = 1, ns
                        MPIbuff1(k) = llfxs(is, i); k = k + 1;
                        MPIbuff1(k) = llfys(is, i); k = k + 1;
                        MPIbuff1(k) = llfzs(is, i); k = k + 1;
                   
                        MPIbuff1(k) = llkefxt(is, i); k = k + 1;
                        MPIbuff1(k) = llkefyt(is, i); k = k + 1;
                        MPIbuff1(k) = llkefzt(is, i); k = k + 1;

                enddo
        enddo
        MPIbuff1(k) = qvirial; k = k + 1;   ! OI 
        call MPI_ALLREDUCE(MPIbuff1, MPIbuff2, j, MPI_DOUBLE_PRECISION, MPI_SUM,MPI_COMM_WORLD, MPIierr);
        k = 0;
        do i = 1, nm
                do is = 1, ns
                        fxs(is, i) = fxs(is, i) + MPIbuff2(k); k = k + 1;
                        fys(is, i) = fys(is, i) + MPIbuff2(k); k = k + 1;
                        fzs(is, i) = fzs(is, i) + MPIbuff2(k); k = k + 1;

                        Etotx(is, i) = Etotx(is, i) + MPIbuff2(k); k = k + 1;
                        Etoty(is, i) = Etoty(is, i) + MPIbuff2(k); k = k + 1;
                        Etotz(is, i) = Etotz(is, i) + MPIbuff2(k); k = k + 1;

                enddo
        enddo
        qvirial = MPIbuff2(k); k = k + 1;   ! OI
        deallocate(MPIbuff1);
        deallocate(MPIbuff2);
! MPI sync
    deallocate(llfxs);
    deallocate(llfys);
    deallocate(llfzs);

    deallocate(lkefxt)
    deallocate(lkefyt)
    deallocate(lkefzt)
    deallocate(llkefxt)
    deallocate(llkefyt)
    deallocate(llkefzt)

    vir_kind = 2.d0*r4pie0*qvirial
    vir = vir + vir_kind

!    write(*,*) 'vir_ind_kspace', vir_kind

        return;
  end subroutine kspace_forces

  subroutine self_forces(alp)
    implicit none
    real(8) :: alp

    integer :: is,js,i,j,k
    real(8) :: expar2
    real(8) :: drr,dist2
    real(8) :: qis,qjs
    real(8) :: ddx,ddy,ddz,d1i,d2i,d3i,d4i,d5i,d6i,d7i,d5ii,d3ii,d7ii

    real(8) :: ef0gxx,ef0gyy,ef0gzz
    real(8) :: ef0gxy,ef0gyx,ef0gxz
    real(8) :: ef0gzx,ef0gyz,ef0gzy

    real(8) :: efigxx,efigyy,efigzz
    real(8) :: efigxy,efigyx,efigxz
    real(8) :: efigzx,efigyz,efigzy

    real(8) :: efi0gxx,efi0gyy,efi0gzz
    real(8) :: efi0gxy,efi0gyx,efi0gxz
    real(8) :: efi0gzx,efi0gyz,efi0gzy

    real(8) :: qvirial,lqvirial

    real(8) :: efldx,efldy,efldz
    real(8) :: fxss,fyss,fzss
    real(8) :: term1,term2,dotis,dotjs

    real(8) :: selffac

!    real(8) :: sefgxx(nsite,nom),sefgyy(nsite,nom),sefgzz(nsite,nom)
!    real(8) :: sefgxy(nsite,nom),sefgyx(nsite,nom),sefgxz(nsite,nom)
!    real(8) :: sefgzx(nsite,nom),sefgyz(nsite,nom),sefgzy(nsite,nom)
!    real(8) :: sefxt(nsite,nom),sefyt(nsite,nom),sefzt(nsite,nom)
    real(8), allocatable :: sefgxx(:, :),sefgyy(:, :),sefgzz(:, :)
    real(8), allocatable :: sefgxy(:, :),sefgyx(:, :),sefgxz(:, :)
    real(8), allocatable :: sefgzx(:, :),sefgyz(:, :),sefgzy(:, :)
    real(8), allocatable :: sefxt(:, :),sefyt(:, :),sefzt(:, :)

    real(8), allocatable :: lsefxt(:, :),lsefyt(:, :),lsefzt(:, :)
    real(8), allocatable :: llsefxt(:, :),llsefyt(:, :),llsefzt(:, :)
	! MPI
	INTEGER(4) MPIchunk_size_loc, iii;
	real(8), allocatable :: MPIbuff1(:), MPIbuff2(:);
    real(8), allocatable :: llfxs(:, :), llfys(:, :), llfzs(:, :);
    real(8), allocatable :: qlfxs(:, :), qlfys(:, :), qlfzs(:, :);

	! DBG
	real(4) iniTime, deltaTime;
	! DBG
	! MPI

    allocate(sefgxx(nsite,nom))
    allocate(sefgyy(nsite,nom))
    allocate(sefgzz(nsite,nom))
    allocate(sefgxy(nsite,nom))
    allocate(sefgyx(nsite,nom))
    allocate(sefgxz(nsite,nom))
    allocate(sefgzx(nsite,nom))
    allocate(sefgyz(nsite,nom))
    allocate(sefgzy(nsite,nom))
    allocate(sefxt(nsite,nom))
    allocate(sefyt(nsite,nom))
    allocate(sefzt(nsite,nom))

    allocate(lsefxt(nsite,nom))
    allocate(lsefyt(nsite,nom))
    allocate(lsefzt(nsite,nom))
    allocate(llsefxt(nsite,nom))
    allocate(llsefyt(nsite,nom))
    allocate(llsefzt(nsite,nom))

    qvirial = 0.d0

    do i=1,nm
       do is=1,ns

          sefgxx(is,i) = 0.d0   ! total gradient of efield
          sefgyy(is,i) = 0.d0
          sefgzz(is,i) = 0.d0
          sefgxy(is,i) = 0.d0
          sefgyx(is,i) = 0.d0
          sefgxz(is,i) = 0.d0
          sefgzx(is,i) = 0.d0
          sefgyz(is,i) = 0.d0
          sefgzy(is,i) = 0.d0

          sefxt(is,i) = 0.d0
          sefyt(is,i) = 0.d0
          sefzt(is,i) = 0.d0

          lsefxt(is,i) = 0.d0
          lsefyt(is,i) = 0.d0
          lsefzt(is,i) = 0.d0

          llsefxt(is,i) = 0.d0
          llsefyt(is,i) = 0.d0
          llsefzt(is,i) = 0.d0

       end do
    end do


    !   point self force

    selffac = 4.d0*alp**(3.d0)/(3.d0*sqrt(pi))
!    write(*,*) selffac

    allocate(llfxs(ns,nm));
    allocate(llfys(ns,nm));
    allocate(llfzs(ns,nm));

    allocate(qlfxs(ns,nm));
    allocate(qlfys(ns,nm));
    allocate(qlfzs(ns,nm));


    do i = 1, nm
       do j = 1, ns
          llfxs(j, i) = 0.0D0;
          llfys(j, i) = 0.0D0;
          llfzs(j, i) = 0.0D0;

          qlfxs(j, i) = 0.0D0;
          qlfys(j, i) = 0.0D0;
          qlfzs(j, i) = 0.0D0;
       enddo
    enddo

!$omp parallel DEFAULT(SHARED)&
!$omp& private(is,js,i,iii,j,k,id)&
!$omp& private(expar2,drr,dist2,qis,qjs,ddx,ddy,ddz,d1i,d2i,d3i,d4i,d5i,d6i,d7i,d5ii,d3ii,d7ii)&
!$omp& private(ef0gxx,ef0gyy,ef0gzz,ef0gxy,ef0gyx,ef0gxz,ef0gzx,ef0gyz,ef0gzy,efigxx,efigyy,efigzz)&
!$omp& private(efi0gxx,efi0gyy,efi0gzz,efi0gxy,efi0gyx,efi0gxz,efi0gzx,efi0gyz,efi0gzy,lqvirial)&
!$omp& private(efigxy,efigyx,efigxz,efigzx,efigyz,efigzy,efldx,efldy,efldz,fxss,fyss,fzss)&
!$omp& private(term1,term2,dotis,dotjs)
!$omp do

    lqvirial = 0.d0

    do i=1,nm
       do is=1,ns

          qis = chgs(is,i)

          sefxt(is,i) = sefxt(is,i) + selffac*idmx(is,i)
          sefyt(is,i) = sefyt(is,i) + selffac*idmy(is,i)
          sefzt(is,i) = sefzt(is,i) + selffac*idmz(is,i)

          sefgxx(is,i) = sefgxx(is,i) - selffac*qis
          sefgyy(is,i) = sefgyy(is,i) - selffac*qis
          sefgzz(is,i) = sefgzz(is,i) - selffac*qis  

          lsefxt(is,i) = lsefxt(is,i) + selffac*idmx(is,i)
          lsefyt(is,i) = lsefyt(is,i) + selffac*idmy(is,i)
          lsefzt(is,i) = lsefzt(is,i) + selffac*idmz(is,i)

!          sefxt(is,i) = sefxt(is,i) - selffac*idmx(is,i)
!          sefyt(is,i) = sefyt(is,i) - selffac*idmy(is,i)
!          sefzt(is,i) = sefzt(is,i) - selffac*idmz(is,i)

!          sefgxx(is,i) = sefgxx(is,i) + selffac*qis
!          sefgyy(is,i) = sefgyy(is,i) + selffac*qis
!          sefgzz(is,i) = sefgzz(is,i) + selffac*qis

!          lsefxt(is,i) = lsefxt(is,i)  - selffac*idmx(is,i)
!          lsefyt(is,i) = lsefyt(is,i)  - selffac*idmy(is,i)
!          lsefzt(is,i) = lsefzt(is,i)  - selffac*idmz(is,i)

!          write(*,*) i,is,sefxt(is,i),sefyt(is,i),sefzt(is,i)
!          write(*,*) i,is,sefgxx(is,i),sefgyy(is,i),sefgzz(is,i)

       end do
    end do
!$omp end do
!$omp barrier

   ! molecular self force
	! DBG
	iniTime = SECNDS(0.0);
	! DBG

!    do i=1,nm
!	 do i = MPIrank + 1, nm, MPIcomm_size
	if (MPIcomm_size .GT. 1) then
		i = 1;
	else
		i = 0;
	endif
	do while (i .LE. nm)
	  if (MPIrank .EQ. 0) then
		if (MPIcomm_size .GT. 1) then
			call MPI_RECV(MPIaddr, 1, MPI_INTEGER, MPI_ANY_SOURCE, MPI_TAG_CommReady, MPI_COMM_WORLD, MPIstat, MPIierr);
			call MPI_SEND(i, 1, MPI_INTEGER, MPIaddr, MPI_TAG_SendData, MPI_COMM_WORLD, MPIierr);
		endif
		i = i + 1;
	  else
		if (MPIcomm_size .GT. 1) then
			call MPI_SEND(MPIrank, 1, MPI_INTEGER, 0, MPI_TAG_CommReady, MPI_COMM_WORLD, MPIierr);
			call MPI_RECV(i, 1, MPI_INTEGER, 0, MPI_TAG_SendData, MPI_COMM_WORLD, MPIstat, MPIierr);
		endif
	  endif
 	  if ((i .LE. (nm)) .AND. ((MPIrank .GE. 1) .OR. (MPIcomm_size .EQ. 1))) then
	   !$omp do schedule(dynamic)
       do is=1,ns-1 
          do js=is+1,ns

             qis = chgs(is,i)
             qjs = chgs(js,i)

!              if (qis*qjs .eq. 0.d0) then

!              ddx = xs(is,i) - xs(js,i) !pcoord(1,is)-pcoord(1,js)
!              ddy = ys(is,i) - ys(js,i) !pcoord(2,is)-pcoord(2,js)
!              ddz = zs(is,i) - zs(js,i) !pcoord(3,is)-pcoord(3,js)

              ddx = (x(i)+xs(is,i)) - (x(i)+xs(js,i))
              ddy = (y(i)+ys(is,i)) - (y(i)+ys(js,i))
              ddz = (z(i)+zs(is,i)) - (z(i)+zs(js,i))

              dist2 = ddx**2+ddy**2+ddz**2

!              if (dist2.eq.0.d0) goto 12
              if (dist2 .NE. 0.0D0) then

              drr = sqrt(dist2)

              d1i = 1.d0/drr
              d2i = d1i**2
              d3i = d2i*d1i
              d4i = d2i*d2i
              d5i = d3i*d2i
              d6i = d3i*d3i
              d7i = d6i*d1i 

              expar2 = exp(-(alp*drr)**2.d0)
              
              ! dot product of induced dipole and separation vector
              dotis = idmx(is,i)*ddx+idmy(is,i)*ddy+idmz(is,i)*ddz
              dotjs = idmx(js,i)*ddx+idmy(js,i)*ddy+idmz(js,i)*ddz

              d3ii = erf(alp*drr)*d3i - twosqpi*alp*expar2*d2i
 
              d5ii = erf(alp*drr)*d5i - twosqpi*alp*expar2*d4i &
                     - (2.d0/3.d0)*twosqpi*alp**(3.d0)*expar2*d2i

              d7ii = erf(alp*drr)*d7i - twosqpi*alp*expar2*d6i &
                     - (2.d0/3.d0)*twosqpi*alp**(3.d0)*expar2*d4i &
                     - (4.d0/15.d0)*twosqpi*alp**(5.d0)*expar2*d2i

              efldx = - (3.d0*d5ii*dotjs*ddx - idmx(js,i)*d3ii)
              efldy = - (3.d0*d5ii*dotjs*ddy - idmy(js,i)*d3ii)
              efldz = - (3.d0*d5ii*dotjs*ddz - idmz(js,i)*d3ii)

              sefxt(is,i) = sefxt(is,i) + efldx
              sefyt(is,i) = sefyt(is,i) + efldy
              sefzt(is,i) = sefzt(is,i) + efldz

              lsefxt(is,i) = lsefxt(is,i) + efldx - qjs*d3ii*ddx
              lsefyt(is,i) = lsefyt(is,i) + efldy - qjs*d3ii*ddy
              lsefzt(is,i) = lsefzt(is,i) + efldz - qjs*d3ii*ddz

              !write(*,*) i,is,js,sefxt(is,i),sefyt(is,i),sefzt(is,i)

              ! virial contribution ! OI
              lqvirial = lqvirial - qis*(efldx*ddx + efldy*ddy + efldz*ddz)         

              ! gradient of efield of permanent charges w.r.t position of pol. centers

              ef0gxx = -qjs*(-3.d0*d5ii*ddx*ddx + d3ii)
              ef0gyy = -qjs*(-3.d0*d5ii*ddy*ddy + d3ii)
              ef0gzz = -qjs*(-3.d0*d5ii*ddz*ddz + d3ii)
              ef0gxy = -qjs*(-3.d0*d5ii*ddx*ddy)
              ef0gyx = -qjs*(-3.d0*d5ii*ddy*ddx)
              ef0gxz = -qjs*(-3.d0*d5ii*ddx*ddz)
              ef0gzx = -qjs*(-3.d0*d5ii*ddz*ddx)
              ef0gyz = -qjs*(-3.d0*d5ii*ddy*ddz)
              ef0gzy = -qjs*(-3.d0*d5ii*ddz*ddy)

              efigxx = -(-15.d0*d7ii*ddx*ddx*dotjs + 3.d0*d5ii*(dotjs + 2.d0*ddx*idmx(js,i)) )
              efigyy = -(-15.d0*d7ii*ddy*ddy*dotjs + 3.d0*d5ii*(dotjs + 2.d0*ddy*idmy(js,i)) )
              efigzz = -(-15.d0*d7ii*ddz*ddz*dotjs + 3.d0*d5ii*(dotjs + 2.d0*ddz*idmz(js,i)) )
              efigxy = -(-15.d0*d7ii*ddx*ddy*dotjs + 3.d0*d5ii*(idmx(js,i)*ddy + idmy(js,i)*ddx) )
              efigyx = -(-15.d0*d7ii*ddy*ddx*dotjs + 3.d0*d5ii*(idmy(js,i)*ddx + idmx(js,i)*ddy) )
              efigxz = -(-15.d0*d7ii*ddx*ddz*dotjs + 3.d0*d5ii*(idmx(js,i)*ddz + idmz(js,i)*ddx) )
              efigzx = -(-15.d0*d7ii*ddz*ddx*dotjs + 3.d0*d5ii*(idmz(js,i)*ddx + idmx(js,i)*ddz) )
              efigyz = -(-15.d0*d7ii*ddy*ddz*dotjs + 3.d0*d5ii*(idmy(js,i)*ddz + idmz(js,i)*ddy) )
              efigzy = -(-15.d0*d7ii*ddz*ddy*dotjs + 3.d0*d5ii*(idmz(js,i)*ddy + idmy(js,i)*ddz) )

              sefgxx(is,i) = sefgxx(is,i) + ef0gxx + efigxx
              sefgyy(is,i) = sefgyy(is,i) + ef0gyy + efigyy
              sefgzz(is,i) = sefgzz(is,i) + ef0gzz + efigzz
              sefgxy(is,i) = sefgxy(is,i) + ef0gxy + efigxy
              sefgyx(is,i) = sefgyx(is,i) + ef0gyx + efigyx
              sefgxz(is,i) = sefgxz(is,i) + ef0gxz + efigxz
              sefgzx(is,i) = sefgzx(is,i) + ef0gzx + efigzx
              sefgyz(is,i) = sefgyz(is,i) + ef0gyz + efigyz
              sefgzy(is,i) = sefgzy(is,i) + ef0gzy + efigzy

!              write(*,*) i,is,js,sefgxx(is,i),sefgyy(is,i),sefgzz(is,i)

              ! virial contribution  ! OI

              efi0gxx = ef0gxx + efigxx
              efi0gyy = ef0gyy + efigyy
              efi0gzz = ef0gzz + efigzz
              efi0gxy = ef0gxy + efigxy
              efi0gyx = ef0gyx + efigyx
              efi0gxz = ef0gxz + efigxz
              efi0gzx = ef0gzx + efigzx
              efi0gyz = ef0gyz + efigyz
              efi0gzy = ef0gzy + efigzy

              lqvirial=lqvirial - (  (idmx(is,i)*efi0gxx + idmy(is,i)*efi0gyx + idmz(is,i)*efi0gzx)*ddx &
                                   + (idmx(is,i)*efi0gxy + idmy(is,i)*efi0gyy + idmz(is,i)*efi0gzy)*ddy &
                                   + (idmx(is,i)*efi0gxz + idmy(is,i)*efi0gyz + idmz(is,i)*efi0gzz)*ddz ) 

!              lqvirial=lqvirial-( idmx(is,i)*efi0gxx*ddx+idmy(is,i)*efi0gyx*ddx+idmz(is,i)*efi0gzx*ddx &
!                                 +idmx(is,i)*efi0gxy*ddy+idmy(is,i)*efi0gyy*ddy+idmz(is,i)*efi0gzy*ddy &
!                                 +idmx(is,i)*efi0gxz*ddz+idmy(is,i)*efi0gyz*ddz+idmz(is,i)*efi0gzz*ddz) 

!              lqvirial=lqvirial-( idmx(is,i)*efi0gxx*ddx+idmx(is,i)*efi0gyx*ddy+idmx(is,i)*efi0gzx*ddz &
!                                 +idmy(is,i)*efi0gxy*ddx+idmy(is,i)*efi0gyy*ddy+idmy(is,i)*efi0gzy*ddz &
!                                 +idmz(is,i)*efi0gxz*ddx+idmz(is,i)*efi0gyz*ddy+idmz(is,i)*efi0gzz*ddz)

             ! gradient of efield of permanent charges w.r.t position of
             ! charge sites

              efldx = -(3.d0*d5ii*dotis*ddx - idmx(is,i)*d3ii)
              efldy = -(3.d0*d5ii*dotis*ddy - idmy(is,i)*d3ii)
              efldz = -(3.d0*d5ii*dotis*ddz - idmz(is,i)*d3ii)

              sefxt(js,i) = sefxt(js,i) + efldx
              sefyt(js,i) = sefyt(js,i) + efldy
              sefzt(js,i) = sefzt(js,i) + efldz

              lsefxt(js,i) = lsefxt(js,i) + efldx + qis*d3ii*ddx 
              lsefyt(js,i) = lsefyt(js,i) + efldy + qis*d3ii*ddy
              lsefzt(js,i) = lsefzt(js,i) + efldz + qis*d3ii*ddz

              ! gradient of efield of permanent charges in ith molecule at
              ! is site

              ef0gxx = -qis*(-3.d0*d5ii*ddx*ddx + d3ii)
              ef0gyy = -qis*(-3.d0*d5ii*ddy*ddy + d3ii)
              ef0gzz = -qis*(-3.d0*d5ii*ddz*ddz + d3ii)
              ef0gxy = -qis*(-3.d0*d5ii*ddx*ddy)
              ef0gyx = -qis*(-3.d0*d5ii*ddy*ddx)
              ef0gxz = -qis*(-3.d0*d5ii*ddx*ddz)
              ef0gzx = -qis*(-3.d0*d5ii*ddz*ddx)
              ef0gyz = -qis*(-3.d0*d5ii*ddy*ddz)
              ef0gzy = -qis*(-3.d0*d5ii*ddz*ddy)

              ! gradient of efield of induced dipoles in ith molecule at is
              ! site

              efigxx = -(15.d0*d7ii*ddx*ddx*dotis + 3.d0*d5ii*(-dotis - 2.d0*ddx*idmx(is,i)) )
              efigyy = -(15.d0*d7ii*ddy*ddy*dotis + 3.d0*d5ii*(-dotis - 2.d0*ddy*idmy(is,i)) )
              efigzz = -(15.d0*d7ii*ddz*ddz*dotis + 3.d0*d5ii*(-dotis - 2.d0*ddz*idmz(is,i)) )
              efigxy = -(15.d0*d7ii*ddx*ddy*dotis + 3.d0*d5ii*(-idmx(is,i)*ddy - idmy(is,i)*ddx) )
              efigyx = -(15.d0*d7ii*ddy*ddx*dotis + 3.d0*d5ii*(-idmy(is,i)*ddx - idmx(is,i)*ddy) )
              efigxz = -(15.d0*d7ii*ddx*ddz*dotis + 3.d0*d5ii*(-idmx(is,i)*ddz - idmz(is,i)*ddx) )
              efigzx = -(15.d0*d7ii*ddz*ddx*dotis + 3.d0*d5ii*(-idmz(is,i)*ddx - idmx(is,i)*ddz) )
              efigyz = -(15.d0*d7ii*ddy*ddz*dotis + 3.d0*d5ii*(-idmy(is,i)*ddz - idmz(is,i)*ddy) )
              efigzy = -(15.d0*d7ii*ddz*ddy*dotis + 3.d0*d5ii*(-idmz(is,i)*ddy - idmy(is,i)*ddz) )

              sefgxx(js,i) = sefgxx(js,i) + ef0gxx + efigxx
              sefgyy(js,i) = sefgyy(js,i) + ef0gyy + efigyy
              sefgzz(js,i) = sefgzz(js,i) + ef0gzz + efigzz
              sefgxy(js,i) = sefgxy(js,i) + ef0gxy + efigxy
              sefgyx(js,i) = sefgyx(js,i) + ef0gyx + efigyx
              sefgxz(js,i) = sefgxz(js,i) + ef0gxz + efigxz
              sefgzx(js,i) = sefgzx(js,i) + ef0gzx + efigzx
              sefgyz(js,i) = sefgyz(js,i) + ef0gyz + efigyz
              sefgzy(js,i) = sefgzy(js,i) + ef0gzy + efigzy

!              write(*,*) ef0gxx, efigxx
 
              endif
!              12 continue 
!!              end if

          end do
       end do
	   !$omp end do
	  endif
    end do
	if (MPIrank .EQ. 0) then
		do i = 1, MPIcomm_size - 1
			call MPI_RECV(MPIaddr, 1, MPI_INTEGER, MPI_ANY_SOURCE, MPI_TAG_CommReady, MPI_COMM_WORLD, MPIstat, MPIierr);
			call MPI_SEND(nm + 1, 1, MPI_INTEGER, MPIaddr, MPI_TAG_SendData, MPI_COMM_WORLD, MPIierr);
		enddo
	endif
!$omp barrier
	! DBG
	deltaTime = SECNDS(iniTime);
!	write(dbgUnit, *) 'P#', MPIrank, '. self_forces calc time	', deltaTime;
	! DBG

!$omp critical

    qvirial = qvirial + lqvirial

!$omp end critical
!$omp barrier

!$omp do

    do i=1,nm
       do is=1,ns

         llsefxt(is,i) = llsefxt(is,i) + lsefxt(is,i)!*r4pie0
         llsefyt(is,i) = llsefyt(is,i) + lsefyt(is,i)!*r4pie0
         llsefzt(is,i) = llsefzt(is,i) + lsefzt(is,i)!*r4pie0

!         write(*,*) llfxs(is,i),llfys(is,i),llfzs(is,i)

       end do
    end do


    do i=1,nm
       do is=1,ns

         qis = chgs(is,i)

         fxss = qis*sefxt(is,i) + idmx(is,i)*sefgxx(is,i) + idmy(is,i)*sefgyx(is,i) + idmz(is,i)*sefgzx(is,i)
         fyss = qis*sefyt(is,i) + idmx(is,i)*sefgxy(is,i) + idmy(is,i)*sefgyy(is,i) + idmz(is,i)*sefgzy(is,i)
         fzss = qis*sefzt(is,i) + idmx(is,i)*sefgxz(is,i) + idmy(is,i)*sefgyz(is,i) + idmz(is,i)*sefgzz(is,i)

         llfxs(is,i)= llfxs(is,i) + fxss*r4pie0
         llfys(is,i)= llfys(is,i) + fyss*r4pie0
         llfzs(is,i)= llfzs(is,i) + fzss*r4pie0

!         write(*,*) llfxs(is,i),llfys(is,i),llfzs(is,i) 

       end do
    end do    
!$omp end do
!$omp end parallel

    deallocate(sefgxx)
    deallocate(sefgyy)
    deallocate(sefgzz)
    deallocate(sefgxy)
    deallocate(sefgyx)
    deallocate(sefgxz)
    deallocate(sefgzx)
    deallocate(sefgyz)
    deallocate(sefgzy)
    deallocate(sefxt)
    deallocate(sefyt)
    deallocate(sefzt)
! MPI sync
	j =  nm * ns * 3 + nm * ns * 3 + 1; !fxs, fys, fzs,Etotx,Etoty,Etotz, qvirial
	allocate(MPIbuff1(0:j - 1));
	allocate(MPIbuff2(0:j - 1));
	do k = 0, j - 1
		MPIbuff2(k) = 0.0D0;
	enddo
	k = 0;
	do i = 1, nm
		do is = 1, ns
			MPIbuff1(k) = llfxs(is, i); k = k + 1;
			MPIbuff1(k) = llfys(is, i); k = k + 1;
			MPIbuff1(k) = llfzs(is, i); k = k + 1;

                        MPIbuff1(k) = llsefxt(is, i); k = k + 1;
                        MPIbuff1(k) = llsefyt(is, i); k = k + 1;
                        MPIbuff1(k) = llsefzt(is, i); k = k + 1;

		enddo
	enddo
        MPIbuff1(k) = qvirial; k = k + 1;   ! OI
	call MPI_ALLREDUCE(MPIbuff1, MPIbuff2, j, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, MPIierr);
	k = 0;
	do i = 1, nm
		do is = 1, ns
			fxs(is, i) = fxs(is, i) + MPIbuff2(k); k = k + 1;
			fys(is, i) = fys(is, i) + MPIbuff2(k); k = k + 1;
			fzs(is, i) = fzs(is, i) + MPIbuff2(k); k = k + 1;

                        Etotx(is, i) = Etotx(is, i) + MPIbuff2(k); k = k + 1;
                        Etoty(is, i) = Etoty(is, i) + MPIbuff2(k); k = k + 1;
                        Etotz(is, i) = Etotz(is, i) + MPIbuff2(k); k = k + 1;

		enddo
	enddo
        qvirial = MPIbuff2(k); k = k + 1; 
	deallocate(MPIbuff1);
	deallocate(MPIbuff2);
! MPI sync
    deallocate(llfxs);
    deallocate(llfys);
    deallocate(llfzs);

    deallocate(lsefxt);
    deallocate(lsefyt);
    deallocate(lsefzt);
    deallocate(llsefxt);
    deallocate(llsefyt);
    deallocate(llsefzt);

!    deallocate(qlfxs);
!    deallocate(qlfys);
!    deallocate(qlfzs);

    vir_selfind = r4pie0*qvirial
    vir = vir + vir_selfind

!    write(*,*) 'vir_ind_self', vir_selfind !! qlfxs(ns,nm) 

!    do i = 1, nm
!       do is = 1, ns 
!          write(*,*) i,is,fxs(is, i),fys(is, i),fzs(is, i)
!       end do
!    end do
    
	return;
  end subroutine self_forces

  subroutine surf_forces(xx,yy,zz,xxs,yys,zzs,boxi)
    implicit none
    real(8) :: xx(nom),yy(nom),zz(nom)
    real(8) :: xxs(nsite,nom),yys(nsite,nom),zzs(nsite,nom)
	real(8) :: boxi

    integer :: is,js,i,j,k
    real(8) :: qis,volmi,fact
    real(8) :: sumqrx,sumqry,sumqrz,sumdx,sumdy,sumdz
    real(8) :: fxss,fyss,fzss
!    real(8) :: srefxt(nsite,nom),srefyt(nsite,nom),srefzt(nsite,nom)
    real(8), allocatable :: srefxt(:, :), srefyt(:, :), srefzt(:, :)

    allocate(srefxt(nsite,nom))
    allocate(srefyt(nsite,nom))
    allocate(srefzt(nsite,nom))

    do i=1,nm
       do is=1,ns
          srefxt(is,i) = 0.d0
          srefyt(is,i) = 0.d0
          srefzt(is,i) = 0.d0
       end do
    end do

! GRU: Never used local variables
!	sumqrx = 0.d0
!	sumqry = 0.d0
!	sumqrz = 0.d0
! GRU: Never used local variables

   sumdx = 0.d0
   sumdy = 0.d0
   sumdz = 0.d0

    volmi = boxi**3.d0
    fact = 4.0d0*pi/(3.0d0*volmi)

!$omp parallel DEFAULT(SHARED)&
!$omp& private(is,js,i,j,k,qis,fxss,fyss,fzss)
!$omp master
    do i=1,nm
       do is=1,ns
          qis = charge(is)
!          sumqrxA(i) = sumqrxA(i) + qis*(xx(i)+xxs(is,i))
!          sumqryA(i) = sumqryA(i) + qis*(yy(i)+yys(is,i))
!          sumqrzA(i) = sumqrzA(i) + qis*(zz(i)+zzs(is,i))

          sumdx = sumdx + idmx(is,i)
          sumdy = sumdy + idmy(is,i)
          sumdz = sumdz + idmz(is,i)
       end do
    end do
!$omp end master
!$omp barrier

!$omp do
    do i=1,nm
       do is=1,ns
          srefxt(is,i) = srefxt(is,i) - fact*sumdx
          srefyt(is,i) = srefyt(is,i) - fact*sumdy
          srefzt(is,i) = srefzt(is,i) - fact*sumdz
       end do
    end do
!$omp end do
!$omp barrier

!$omp do
    do i=1,nm
       do is=1,ns
         qis = chgs(is,i)

         fxss = qis*srefxt(is,i)
         fyss = qis*srefyt(is,i)
         fzss = qis*srefzt(is,i)

         fxs(is,i)= fxs(is,i) + fxss*r4pie0
         fys(is,i)= fys(is,i) + fyss*r4pie0
         fzs(is,i)= fzs(is,i) + fzss*r4pie0
       end do
    end do
!$omp end do
!$omp end parallel

    deallocate(srefxt)
    deallocate(srefyt)
    deallocate(srefzt)

  end subroutine surf_forces
    
! *************** Short Range Damped Induction terms ***************************

  subroutine shortrange_efield_stat_damped(xx,yy,zz,xxs,yys,zzs,boxi)
    implicit none
    real(8) :: xx(nom),yy(nom),zz(nom)
    real(8) :: xxs(nsite,nom),yys(nsite,nom),zzs(nsite,nom)
    real(8) :: boxi

    integer :: is,js,i,j,k,isteps
    real(8) :: efldx,efldy,efldz
    real(8) :: dist2,rccsq,dist,doti,dotj
    real(8) :: qis,qjs  ! charges on sites
    real(8) :: xcc,ycc,zcc,dx,dy,dz,ddx,ddy,ddz,d1i,d2i,d3i,d5i,d7i
    real(8) :: polr,f2,f3  
    real(8) :: qindpe
    real(8) :: expdr,detr

!    real(8), parameter :: delta2 = 1.65d0/bohr2a  ! 1/Bohr to 1/Ang
!    real(8), parameter :: delta3 = 1.55d0/bohr2a
	! MPI
     INTEGER(4) MPIchunk_size_loc, iii;
     real(8), allocatable :: MPIbuff1(:), MPIbuff2(:);
	! DBG
     real(4) iniTime, deltaTime;
	! DBG
	! MPI

    do i=1,nm
       do is=1,ns
          shtE0x(is,i)= 0.d0 
          shtE0y(is,i)= 0.d0
          shtE0z(is,i)= 0.d0
       end do
    end do

! ************* Calculating Electric from permanant charges ****************
! starts here
	! DBG
	iniTime = SECNDS(0.0);
	! DBG

!$omp parallel DEFAULT(SHARED)&
!$omp& private(is,js,i,iii,j,k,isteps)&
!$omp& private(efldx,efldy,efldz,dist2,rccsq,dist,doti,dotj,qis,qjs)&
!$omp& private(xcc,ycc,zcc,dx,dy,dz,ddx,ddy,ddz,d1i,d2i,d3i,d5i,d7i)&
!$omp& private(polr,f2,f3,qindpe,expdr,detr)
!    do i=1,nm       ! loop over molecule i   starts
!	 do i = MPIrank + 1, nm, MPIcomm_size
	if (MPIcomm_size .GT. 1) then
		i = 1;
	else
		i = 0;
	endif
	do while (i .LE. nm)
	  if (MPIrank .EQ. 0) then
		if (MPIcomm_size .GT. 1) then
			call MPI_RECV(MPIaddr, 1, MPI_INTEGER, MPI_ANY_SOURCE, MPI_TAG_CommReady, MPI_COMM_WORLD, MPIstat, MPIierr);
			call MPI_SEND(i, 1, MPI_INTEGER, MPIaddr, MPI_TAG_SendData, MPI_COMM_WORLD, MPIierr);
		endif
		i = i + 1;
	  else
		if (MPIcomm_size .GT. 1) then
			call MPI_SEND(MPIrank, 1, MPI_INTEGER, 0, MPI_TAG_CommReady, MPI_COMM_WORLD, MPIierr);
			call MPI_RECV(i, 1, MPI_INTEGER, 0, MPI_TAG_SendData, MPI_COMM_WORLD, MPIstat, MPIierr);
		endif
	  endif
 	  if ((i .LE. (nm)) .AND. ((MPIrank .GE. 1) .OR. (MPIcomm_size .EQ. 1))) then
		!$omp do schedule(dynamic)
        do j=1,nm    ! loop over molecule j   starts

!          if (j.eq.i) goto 17
          if (j .NE. i) then

          dx = xx(i)-xx(j)
          dy = yy(i)-yy(j)
          dz = zz(i)-zz(j)

          xcc = dx - boxi*nint(dx/boxi)
          ycc = dy - boxi*nint(dy/boxi)
          zcc = dz - boxi*nint(dz/boxi)

          rccsq = xcc**2 + ycc**2 + zcc**2

          if (rccsq.lt.rcutsd) then 
!          if (rccsq.lt.3.0d0**2) then
   
             do is=ns-npol+1,ns

                do js=1,ns

                   qis = chgs(is,i)
                   qjs = chgs(js,j)

!                   if (qjs .ne. 0.d0) then

                      ddx = xcc + (xxs(is,i) - xxs(js,j)) 
                      ddy = ycc + (yys(is,i) - yys(js,j)) 
                      ddz = zcc + (zzs(is,i) - zzs(js,j))

                      dist2 = ddx**2 + ddy**2 + ddz**2

!                      if (dist2.lt.rcutsd) then ! cut off on  sites separation
!                      starts

                      dist = sqrt(dist2)
                      d1i = 1.d0/dist
                      d2i = d1i**2
                      d3i = d2i*d1i

                      f2 = tt(2,delta2,dist)
                     ! f2=erf(alpha*dist)

                      ! components of efield in internal units (multiply
                      ! r4pie0)! Ei= qi* vecr_rij/rij^3

                      shtE0x(is,i)= shtE0x(is,i) + (f2 - 1.d0)*qjs*ddx*d3i
                      shtE0y(is,i)= shtE0y(is,i) + (f2 - 1.d0)*qjs*ddy*d3i
                      shtE0z(is,i)= shtE0z(is,i) + (f2 - 1.d0)*qjs*ddz*d3i     

!                       end if  ! cut off on  sites separation ends

!                   end if

                end do

             end do

          end if     ! com-com cutoff check ends

          endif
!     17 continue
        end do        ! loop over molecule j ends
		!$omp end do
	   endif
    end do           ! loop over molecule i ends
	if (MPIrank .EQ. 0) then
		do i = 1, MPIcomm_size - 1
			call MPI_RECV(MPIaddr, 1, MPI_INTEGER, MPI_ANY_SOURCE, MPI_TAG_CommReady, MPI_COMM_WORLD, MPIstat, MPIierr);
			call MPI_SEND(nm + 1, 1, MPI_INTEGER, MPIaddr, MPI_TAG_SendData, MPI_COMM_WORLD, MPIierr);
		enddo
	endif
!$omp end parallel
	! DBG
	deltaTime = SECNDS(iniTime);
!	write(dbgUnit, *) 'P#', MPIrank, '. shortrange_efield_stat_damped calc time	', deltaTime;
	! DBG

! ************* Calculating Electric from permanant charges ends ***********
! MPI sync
	j = nm * ns * 3; !shtE0x, shtE0y, shtE0z
	allocate(MPIbuff1(0:j - 1));
	allocate(MPIbuff2(0:j - 1));
	do k = 0, j - 1
		MPIbuff2(k) = 0.0D0;
	enddo
	k = 0;
	do i = 1, nm
		do is = 1, ns
			MPIbuff1(k) = shtE0x(is, i); k = k + 1;
			MPIbuff1(k) = shtE0y(is, i); k = k + 1;
			MPIbuff1(k) = shtE0z(is, i); k = k + 1;
		enddo
	enddo
	call MPI_ALLREDUCE(MPIbuff1, MPIbuff2, j, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, MPIierr);
	k = 0;
	do i = 1, nm
		do is = 1, ns
			shtE0x(is, i) = MPIbuff2(k); k = k + 1;
			shtE0y(is, i) = MPIbuff2(k); k = k + 1;
			shtE0z(is, i) = MPIbuff2(k); k = k + 1;
		enddo
	enddo
	deallocate(MPIbuff1);
	deallocate(MPIbuff2);
! MPI sync

	return;
  end subroutine shortrange_efield_stat_damped


! *************** Short Range Damped Induction terms ***************************

  subroutine shortrange_efield_idm_damped(xx,yy,zz,xxs,yys,zzs,boxi)
    implicit none
    real(8) :: xx(nom),yy(nom),zz(nom)
    real(8) :: xxs(nsite,nom),yys(nsite,nom),zzs(nsite,nom)
    real(8) :: boxi

    integer :: is,js,i,j,k,isteps
    real(8) :: efldx,efldy,efldz
    real(8) :: dist2,rccsq,dist,doti,dotj
    real(8) :: qis,qjs  ! charges on sites
    real(8) :: xcc,ycc,zcc,dx,dy,dz,ddx,ddy,ddz,d1i,d2i,d3i,d5i,d7i
    real(8) :: polr,f2,f3  
    real(8) :: expdr,detr
    real(8) :: qindpe
!    real(8), parameter :: delta2 = 1.65d0/bohr2a ! 1/Bohr to 1/Ang
!    real(8), parameter :: delta3 = 1.55d0/bohr2a
	! MPI
	INTEGER(4) MPIchunk_size_loc, iii;
	real(8), allocatable :: MPIbuff1(:), MPIbuff2(:);
	! DBG
	real(4) iniTime, deltaTime;
	! DBG
	! MPI

    do i=1,nm
       do is=1,ns
          shtEidmx(is,i)= 0.d0
          shtEidmy(is,i)= 0.d0
          shtEidmz(is,i)= 0.d0
       end do
    end do

!**************************************************************************
! Now Calculating Electric from point induced dipoles
	! DBG
	iniTime = SECNDS(0.0);
	! DBG
!$omp parallel DEFAULT(SHARED)&
!$omp& private(is,js,i,iii,j,k,isteps)&
!$omp& private(efldx,efldy,efldz,dist2,rccsq,dist,doti,dotj,qis,qjs)&
!$omp& private(xcc,ycc,zcc,dx,dy,dz,ddx,ddy,ddz,d1i,d2i,d3i,d5i,d7i)&
!$omp& private(polr,f2,f3,expdr,detr,qindpe)
!    do i=1,nm   ! loop over ith molecule starts
    do i = MPIrank + 1, nm, MPIcomm_size
	   !$omp do schedule(dynamic)
       do j=1,nm   ! loop over jth molecule starts

!          if (j.eq.i) goto 18 !cycle
          if (j .NE. i) then

             dx = xx(i)-xx(j)
             dy = yy(i)-yy(j)
             dz = zz(i)-zz(j)

             xcc = dx - boxi*nint(dx/boxi)
             ycc = dy - boxi*nint(dy/boxi)
             zcc = dz - boxi*nint(dz/boxi)

             rccsq = xcc**2 + ycc**2 + zcc**2

             if (rccsq.lt.rcutsd) then  ! com-com cutoff check starts
!             if (rccsq.lt.3.0d0**2) then
 
                do is=ns-npol+1,ns

                   do js=ns-npol+1,ns    ! loop over sites of molecule j  Starts

                      ddx = xcc + xxs(is,i) - xxs(js,j)
                      ddy = ycc + yys(is,i) - yys(js,j)
                      ddz = zcc + zzs(is,i) - zzs(js,j)

                      dist2 = ddx**2 + ddy**2 + ddz**2

!                      if (dist2.lt.rcutsd) then ! cut off on  sites separation
!                      starts

                      dist = sqrt(dist2)
                      d1i = 1.d0/dist
                      d2i = d1i**2
                      d3i = d2i*d1i
                      d5i = d3i*d2i

                      f3 =tt(3,delta3,dist)
                    !  f3=erf(alpha*dist)
                                            
                      ! dot product of induced dipole and separation vector
                      doti = idmx(is,i)*ddx+idmy(is,i)*ddy+idmz(is,i)*ddz
                      dotj = idmx(js,j)*ddx+idmy(js,j)*ddy+idmz(js,j)*ddz  
     
                      ! calculate the  efield of inducd dipole of molecule j at
                      ! induced dipole of molecule i

                      efldx = (f3 - 1.d0)*(3.d0*d5i*dotj*ddx - idmx(js,j)*d3i)
                      efldy = (f3 - 1.d0)*(3.d0*d5i*dotj*ddy - idmy(js,j)*d3i)
                      efldz = (f3 - 1.d0)*(3.d0*d5i*dotj*ddz - idmz(js,j)*d3i)

                      ! calculate total efield at induced dipole at molecule i
                      ! because of parmenant charges and induced dipoles at
                      ! molecule j

                      shtEidmx(is,i)= shtEidmx(is,i) + efldx
                      shtEidmy(is,i)= shtEidmy(is,i) + efldy
                      shtEidmz(is,i)= shtEidmz(is,i) + efldz 

!                      qindpe = qindpe - 0.5d0*(f3-1.d0)*(efldx*idmx(is,i)+efldy*idmy(is,i)+efldz*idmz(is,i))

!                      end if

                   end do     ! loop over sites js  ends

                end do     ! loop over sites is  ends

             end if   ! com-com cutoff check ends

             endif
!       18 continue
       end do      ! loop over jth molecule ends
	   !$omp end do
    end do  ! loop over ith molecule ends
!$omp end parallel
	! DBG
	deltaTime = SECNDS(iniTime);
!	write(dbgUnit, *) 'P#', MPIrank, '. shortrange_efield_idm_damped calc time	', deltaTime;
	! DBG
!**************************************************************************
! calculation of converged polarization energy        ends here
! MPI sync
	j = nm * ns * 3; !shtEidmx, shtEidmx, shtEidmx
	allocate(MPIbuff1(0:j - 1));
	allocate(MPIbuff2(0:j - 1));
	do k = 0, j - 1
		MPIbuff2(k) = 0.0D0;
	enddo
	k = 0;
	do i = 1, nm
		do is = 1, ns
			MPIbuff1(k) = shtEidmx(is, i); k = k + 1;
			MPIbuff1(k) = shtEidmy(is, i); k = k + 1;
			MPIbuff1(k) = shtEidmz(is, i); k = k + 1;
		enddo
	enddo
	call MPI_ALLREDUCE(MPIbuff1, MPIbuff2, j, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, MPIierr);
	k = 0;
	do i = 1, nm
		do is = 1, ns
			shtEidmx(is, i) = MPIbuff2(k); k = k + 1;
			shtEidmy(is, i) = MPIbuff2(k); k = k + 1;
			shtEidmz(is, i) = MPIbuff2(k); k = k + 1;
		enddo
	enddo
	deallocate(MPIbuff1);
	deallocate(MPIbuff2);
! MPI sync

	return;
  end subroutine shortrange_efield_idm_damped

  subroutine shortrange_efield_idm_damped_v1(xx,yy,zz,xxs,yys,zzs,boxi)
    implicit none
    real(8) :: xx(nom),yy(nom),zz(nom)
    real(8) :: xxs(nsite,nom),yys(nsite,nom),zzs(nsite,nom)
    real(8) :: boxi

    integer :: is,js,i,j,k,isteps
    real(8) :: efldx,efldy,efldz
    real(8) :: dist2,rccsq,dist,doti,dotj
    real(8) :: qis,qjs  ! charges on sites
    real(8) :: xcc,ycc,zcc,dx,dy,dz,ddx,ddy,ddz,d1i,d2i,d3i,d5i,d7i
    real(8) :: polr,f2,f3  
    real(8) :: expdr,detr
    real(8) :: qindpe
!    real(8), parameter :: delta2 = 1.65d0/bohr2a ! 1/Bohr to 1/Ang
!    real(8), parameter :: delta3 = 1.55d0/bohr2a
	! MPI
	INTEGER(4) MPIchunk_size_loc, iii;
	real(8), allocatable :: MPIbuff1(:), MPIbuff2(:);
	! DBG
	real(4) iniTime, deltaTime;
	! DBG
	! MPI

    do i=1,nm
       do is=1,ns
          shtEidmx(is,i)= 0.d0
          shtEidmy(is,i)= 0.d0
          shtEidmz(is,i)= 0.d0
       end do
    end do

!**************************************************************************
! Now Calculating Electric from point induced dipoles
	! DBG
	iniTime = SECNDS(0.0);
	! DBG
!$omp parallel DEFAULT(SHARED)&
!$omp& private(is,js,i,iii,j,k,isteps)&
!$omp& private(efldx,efldy,efldz,dist2,rccsq,dist,doti,dotj,qis,qjs)&
!$omp& private(xcc,ycc,zcc,dx,dy,dz,ddx,ddy,ddz,d1i,d2i,d3i,d5i,d7i)&
!$omp& private(polr,f2,f3,expdr,detr,qindpe)
!    do i=1,nm   ! loop over ith molecule starts
!	 do i = MPIrank + 1, nm, MPIcomm_size
	if (MPIcomm_size .GT. 1) then
		i = 1;
	else
		i = 0;
	endif
	do while (i .LE. nm)
	  if (MPIrank .EQ. 0) then
		if (MPIcomm_size .GT. 1) then
			call MPI_RECV(MPIaddr, 1, MPI_INTEGER, MPI_ANY_SOURCE, MPI_TAG_CommReady, MPI_COMM_WORLD, MPIstat, MPIierr);
			call MPI_SEND(i, 1, MPI_INTEGER, MPIaddr, MPI_TAG_SendData, MPI_COMM_WORLD, MPIierr);
		endif
		i = i + 1;
	  else
		if (MPIcomm_size .GT. 1) then
			call MPI_SEND(MPIrank, 1, MPI_INTEGER, 0, MPI_TAG_CommReady, MPI_COMM_WORLD, MPIierr);
			call MPI_RECV(i, 1, MPI_INTEGER, 0, MPI_TAG_SendData, MPI_COMM_WORLD, MPIstat, MPIierr);
		endif
	  endif
 	  if ((i .LE. (nm)) .AND. ((MPIrank .GE. 1) .OR. (MPIcomm_size .EQ. 1))) then
	   !$omp do schedule(dynamic)
       do j=1,nm   ! loop over jth molecule starts

!          if (j.eq.i) goto 18 !cycle
          if (j .NE. i) then

             dx = xx(i)-xx(j)
             dy = yy(i)-yy(j)
             dz = zz(i)-zz(j)

             xcc = dx - boxi*nint(dx/boxi)
             ycc = dy - boxi*nint(dy/boxi)
             zcc = dz - boxi*nint(dz/boxi)

             rccsq = xcc**2 + ycc**2 + zcc**2

             if (rccsq.lt.rcutsd) then  ! com-com cutoff check starts
!             if (rccsq.lt.3.0d0**2) then
 
                do is=ns-npol+1,ns

                   do js=ns-npol+1,ns    ! loop over sites of molecule j  Starts

                      ddx = xcc + xxs(is,i) - xxs(js,j)
                      ddy = ycc + yys(is,i) - yys(js,j)
                      ddz = zcc + zzs(is,i) - zzs(js,j)

                      dist2 = ddx**2 + ddy**2 + ddz**2

!                      if (dist2.lt.rcutsd) then ! cut off on  sites separation
!                      starts

                      dist = sqrt(dist2)
                      d1i = 1.d0/dist
                      d2i = d1i**2
                      d3i = d2i*d1i
                      d5i = d3i*d2i

                      f3 =tt(3,delta3,dist)
                    !  f3=erf(alpha*dist)
                                            
                      ! dot product of induced dipole and separation vector
                      doti = idmx(is,i)*ddx+idmy(is,i)*ddy+idmz(is,i)*ddz
                      dotj = idmx(js,j)*ddx+idmy(js,j)*ddy+idmz(js,j)*ddz  
     
                      ! calculate the  efield of inducd dipole of molecule j at
                      ! induced dipole of molecule i

                      efldx = (f3 - 1.d0)*(3.d0*d5i*dotj*ddx - idmx(js,j)*d3i)
                      efldy = (f3 - 1.d0)*(3.d0*d5i*dotj*ddy - idmy(js,j)*d3i)
                      efldz = (f3 - 1.d0)*(3.d0*d5i*dotj*ddz - idmz(js,j)*d3i)

                      ! calculate total efield at induced dipole at molecule i
                      ! because of parmenant charges and induced dipoles at
                      ! molecule j

                      shtEidmx(is,i)= shtEidmx(is,i) + efldx
                      shtEidmy(is,i)= shtEidmy(is,i) + efldy
                      shtEidmz(is,i)= shtEidmz(is,i) + efldz 

!                      qindpe = qindpe - 0.5d0*(f3-1.d0)*(efldx*idmx(is,i)+efldy*idmy(is,i)+efldz*idmz(is,i))

!                      end if

                   end do     ! loop over sites js  ends

                end do     ! loop over sites is  ends

             end if   ! com-com cutoff check ends

             endif
!       18 continue
       end do      ! loop over jth molecule ends
	   !$omp end do
	  endif
    end do  ! loop over ith molecule ends
	if (MPIrank .EQ. 0) then
		do i = 1, MPIcomm_size - 1
			call MPI_RECV(MPIaddr, 1, MPI_INTEGER, MPI_ANY_SOURCE, MPI_TAG_CommReady, MPI_COMM_WORLD, MPIstat, MPIierr);
			call MPI_SEND(nm + 1, 1, MPI_INTEGER, MPIaddr, MPI_TAG_SendData, MPI_COMM_WORLD, MPIierr);
		enddo
	endif
!$omp end parallel
	! DBG
	deltaTime = SECNDS(iniTime);
!	write(dbgUnit, *) 'P#', MPIrank, '. shortrange_efield_idm_damped calc time	', deltaTime;
	! DBG
!**************************************************************************
! calculation of converged polarization energy        ends here
! MPI sync
	j = nm * ns * 3; !shtEidmx, shtEidmx, shtEidmx
	allocate(MPIbuff1(0:j - 1));
	allocate(MPIbuff2(0:j - 1));
	do k = 0, j - 1
		MPIbuff2(k) = 0.0D0;
	enddo
	k = 0;
	do i = 1, nm
		do is = 1, ns
			MPIbuff1(k) = shtEidmx(is, i); k = k + 1;
			MPIbuff1(k) = shtEidmy(is, i); k = k + 1;
			MPIbuff1(k) = shtEidmz(is, i); k = k + 1;
		enddo
	enddo
	call MPI_ALLREDUCE(MPIbuff1, MPIbuff2, j, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, MPIierr);
	k = 0;
	do i = 1, nm
		do is = 1, ns
			shtEidmx(is, i) = MPIbuff2(k); k = k + 1;
			shtEidmy(is, i) = MPIbuff2(k); k = k + 1;
			shtEidmz(is, i) = MPIbuff2(k); k = k + 1;
		enddo
	enddo
	deallocate(MPIbuff1);
	deallocate(MPIbuff2);
! MPI sync

	return;
  end subroutine shortrange_efield_idm_damped_v1


  subroutine shortrange_forces_damped(xx,yy,zz,xxs,yys,zzs,boxi)
    implicit none
    real(8) :: xx(nom),yy(nom),zz(nom)
    real(8) :: xxs(nsite,nom),yys(nsite,nom),zzs(nsite,nom)
	real(8) :: boxi

    integer :: is,js,i,j,k,isteps
    real(8) :: dist,d1i,d2i,d3i,d5i,d7i,doti,dotj,f2,f3,df2,df3

    real(8) :: ef0gxx,ef0gyy,ef0gzz
    real(8) :: ef0gxy,ef0gyx,ef0gxz
    real(8) :: ef0gzx,ef0gyz,ef0gzy

    real(8) :: efigxx,efigyy,efigzz
    real(8) :: efigxy,efigyx,efigxz
    real(8) :: efigzx,efigyz,efigzy
    real(8) :: efldx,efldy,efldz

    real(8) :: efi0gxx,efi0gyy,efi0gzz
    real(8) :: efi0gxy,efi0gyx,efi0gxz
    real(8) :: efi0gzx,efi0gyz,efi0gzy

    real(8) :: qvirial,lqvirial

    real(8) :: fxss,fyss,fzss
    real(8) :: qis,qjs  ! charges on sites
    real(8) :: xcc,ycc,zcc,dx,dy,dz,ddx,ddy,ddz,dist2,rccsq
    real(8) :: expdr,detr,expar2,derf

	real(8) :: cst

!    real(8) :: efgxx(nsite,nom),efgyy(nsite,nom),efgzz(nsite,nom)
!    real(8) :: efgxy(nsite,nom),efgyx(nsite,nom),efgxz(nsite,nom)
!    real(8) :: efgzx(nsite,nom),efgyz(nsite,nom),efgzy(nsite,nom)
!    real(8) :: efxt(nsite,nom),efyt(nsite,nom),efzt(nsite,nom)
    real(8), allocatable :: efgxx(:, :),efgyy(:, :),efgzz(:, :)
    real(8), allocatable :: efgxy(:, :),efgyx(:, :),efgxz(:, :)
    real(8), allocatable :: efgzx(:, :),efgyz(:, :),efgzy(:, :)
    real(8), allocatable :: efxt(:, :),efyt(:, :),efzt(:, :)
! Shared
    real(8), allocatable :: sefgxx(:, :),sefgyy(:, :),sefgzz(:, :)
    real(8), allocatable :: sefgxy(:, :),sefgyx(:, :),sefgxz(:, :)
    real(8), allocatable :: sefgzx(:, :),sefgyz(:, :),sefgzy(:, :)
    real(8), allocatable :: sefxt(:, :),sefyt(:, :),sefzt(:, :)

    real(8), allocatable :: lshefxt(:, :), lshefyt(:, :), lshefzt(:, :)
    real(8), allocatable :: llshefxt(:, :), llshefyt(:, :), llshefzt(:, :)

!    real(8), parameter :: delta2 = 1.65d0/bohr2a  ! 1/Bohr to 1/Ang
!    real(8), parameter :: delta3 = 1.55d0/bohr2a
	! MPI
	INTEGER(4) MPIchunk_size_loc, iii;
	real(8), allocatable :: MPIbuff1(:), MPIbuff2(:);
    real(8), allocatable :: llfxs(:, :), llfys(:, :), llfzs(:, :);
	! DBG
	real(4) iniTime, deltaTime;
	! DBG
	! MPI


    allocate(sefgxx(ns,nm))
    allocate(sefgyy(ns,nm))
    allocate(sefgzz(ns,nm))
    allocate(sefgxy(ns,nm))
    allocate(sefgyx(ns,nm))
    allocate(sefgxz(ns,nm))
    allocate(sefgzx(ns,nm))
    allocate(sefgyz(ns,nm))
    allocate(sefgzy(ns,nm))
    allocate(sefxt(ns,nm))
    allocate(sefyt(ns,nm))
    allocate(sefzt(ns,nm))

    allocate(lshefxt(ns,nm))
    allocate(lshefyt(ns,nm))
    allocate(lshefzt(ns,nm))

    allocate(llshefxt(ns,nm))
    allocate(llshefyt(ns,nm))
    allocate(llshefzt(ns,nm)) 

    cst = 1.0d0  ! it should be 1.d0 ideally
    qvirial = 0.d0

    do i=1,nm
       do is=1,ns
          sefgxx(is,i) = 0.d0   ! total gradient of sefield
          sefgyy(is,i) = 0.d0
          sefgzz(is,i) = 0.d0
          sefgxy(is,i) = 0.d0
          sefgyx(is,i) = 0.d0
          sefgxz(is,i) = 0.d0
          sefgzx(is,i) = 0.d0
          sefgyz(is,i) = 0.d0
          sefgzy(is,i) = 0.d0
          sefxt(is,i) = 0.d0
          sefyt(is,i) = 0.d0
          sefzt(is,i) = 0.d0

          lshefxt(is,i) = 0.d0
          lshefyt(is,i) = 0.d0
          lshefzt(is,i) = 0.d0

          llshefxt(is,i) = 0.d0
          llshefyt(is,i) = 0.d0
          llshefzt(is,i) = 0.d0  

       end do
    end do

    allocate(llfxs(ns,nm));
    allocate(llfys(ns,nm));
    allocate(llfzs(ns,nm));
    do i = 1, nm
       do j = 1, ns
          llfxs(j, i) = 0.0D0;
          llfys(j, i) = 0.0D0;
          llfzs(j, i) = 0.0D0;
       enddo
    enddo
	! DBG
	! if (MPIrank .EQ. 0) then
		! write(332, '(a)') '#';
		! write(332, '(a)') 'fxs';
		! write(332, '(1p3E16.6E3)') fxs;
		! write(332, '(a)') 'fys';
		! write(332, '(1p3E16.6E3)') fys;
		! write(332, '(a)') 'fzs';
		! write(332, '(1p3E16.6E3)') fzs;
	! endif
	! DBG

!$omp parallel DEFAULT(SHARED)&
!$omp& private(is,js,i,iii,j,k,isteps)&
!$omp& private(dist,d1i,d2i,d3i,d5i,d7i,doti,dotj,f2,f3,df2,df3,ef0gxx,ef0gyy,ef0gzz)&
!$omp& private(ef0gxy,ef0gyx,ef0gxz,ef0gzx,ef0gyz,ef0gzy,efigxx,efigyy,efigzz)&
!$omp& private(efigxy,efigyx,efigxz,efigzx,efigyz,efigzy,efldx,efldy,efldz)&
!$omp& private(efi0gxx,efi0gyy,efi0gzz,efi0gxy,efi0gyx,efi0gxz,efi0gzx,efi0gyz,efi0gzy,lqvirial)&
!$omp& private(fxss,fyss,fzss,qis,qjs,xcc,ycc,zcc,dx,dy,dz,ddx,ddy,ddz,dist2,rccsq)&
!$omp& private(expdr,detr,expar2,derf,efgxx,efgyy,efgzz,efgxy,efgyx,efgxz,efgzx)&
!$omp& private(efgyz,efgzy,efxt,efyt,efzt)

   lqvirial = 0.d0

   allocate(efgxx(ns,nm))
   allocate(efgyy(ns,nm))
   allocate(efgzz(ns,nm))
   allocate(efgxy(ns,nm))
   allocate(efgyx(ns,nm))
   allocate(efgxz(ns,nm))
   allocate(efgzx(ns,nm))
   allocate(efgyz(ns,nm))
   allocate(efgzy(ns,nm))
   allocate(efxt(ns,nm))
   allocate(efyt(ns,nm))
   allocate(efzt(ns,nm))

    do i=1,nm
       do is=1,ns
          efgxx(is,i) = 0.d0   ! total gradient of efield
          efgyy(is,i) = 0.d0
          efgzz(is,i) = 0.d0
          efgxy(is,i) = 0.d0
          efgyx(is,i) = 0.d0
          efgxz(is,i) = 0.d0
          efgzx(is,i) = 0.d0
          efgyz(is,i) = 0.d0
          efgzy(is,i) = 0.d0
          efxt(is,i) = 0.d0
          efyt(is,i) = 0.d0
          efzt(is,i) = 0.d0
       end do
    end do

	! DBG
	iniTime = SECNDS(0.0);
	! DBG
!    do i=1,nm-1   ! loop over ith molecule starts
!	 do i = MPIrank + 1, nm - 1, MPIcomm_size
	if (MPIcomm_size .GT. 1) then
		i = 1;
	else
		i = 0;
	endif
	do while (i .LE. nm - 1)
	  if (MPIrank .EQ. 0) then
		if (MPIcomm_size .GT. 1) then
			call MPI_RECV(MPIaddr, 1, MPI_INTEGER, MPI_ANY_SOURCE, MPI_TAG_CommReady, MPI_COMM_WORLD, MPIstat, MPIierr);
			call MPI_SEND(i, 1, MPI_INTEGER, MPIaddr, MPI_TAG_SendData, MPI_COMM_WORLD, MPIierr);
		endif
		i = i + 1;
	  else
		if (MPIcomm_size .GT. 1) then
			call MPI_SEND(MPIrank, 1, MPI_INTEGER, 0, MPI_TAG_CommReady, MPI_COMM_WORLD, MPIierr);
			call MPI_RECV(i, 1, MPI_INTEGER, 0, MPI_TAG_SendData, MPI_COMM_WORLD, MPIstat, MPIierr);
		endif
	  endif
 	  if ((i .LE. (nm - 1)) .AND. ((MPIrank .GE. 1) .OR. (MPIcomm_size .EQ. 1))) then
	   !$omp do schedule(dynamic)
       do j=i+1,nm   ! loop over jth molecule starts

          dx = xx(i)-xx(j)
          dy = yy(i)-yy(j)
          dz = zz(i)-zz(j)

          xcc = dx - boxi*nint(dx/boxi)
          ycc = dy - boxi*nint(dy/boxi)
          zcc = dz - boxi*nint(dz/boxi)

          rccsq = xcc**2 + ycc**2 + zcc**2


          if (rccsq.lt.rcutsd) then  ! com-com cutoff check starts

!           if (rccsq.lt.3.0d0**2) then 

             do is=1,ns   ! loop over sites of molecule i  Starts

                do js=1,ns    ! loop over sites of molecule j  Starts

                   qis = chgs(is,i)
                   qjs = chgs(js,j)

                   ddx = xcc + (xxs(is,i) - xxs(js,j))
                   ddy = ycc + (yys(is,i) - yys(js,j))
                   ddz = zcc + (zzs(is,i) - zzs(js,j))

                   dist2 = ddx**2 + ddy**2 + ddz**2

!                   if (dist2.lt.rcutsd) then  ! site cutoff

                   dist = sqrt(dist2)
                   d1i = 1.d0/dist
                   d2i = d1i**2.d0
                   d3i = d2i*d1i
                   d5i = d3i*d2i
                   d7i = d5i*d2i

                   f2 = tt(2,delta2,dist)
                   f3 = tt(3,delta3,dist)
                   df2 = dtt(2,delta2,dist)
                   df3 = dtt(3,delta3,dist)

                   ! test starts
                   !  f2=erf(alpha*dist)
                   !  f3=erf(alpha*dist) 
                   !  expar2 = exp(-(alpha*dist)**2.d0)
                   !  derf = twosqpi*alpha*expar2
                   !  df2 = derf
                   !  df3= derf 
                   ! test ends

                   ! dot product of induced dipole and separation vector
                   doti = idmx(is,i)*ddx+idmy(is,i)*ddy+idmz(is,i)*ddz
                   dotj = idmx(js,j)*ddx+idmy(js,j)*ddy+idmz(js,j)*ddz

                   ! gradient of efield of permanent charges w.r.t position of
                   ! charge sites

                   efldx = (f2-1.d0)*(3.d0*d5i*dotj*ddx-idmx(js,j)*d3i) - dotj*d3i*df2*d1i*ddx
                   efldy = (f2-1.d0)*(3.d0*d5i*dotj*ddy-idmy(js,j)*d3i) - dotj*d3i*df2*d1i*ddy
                   efldz = (f2-1.d0)*(3.d0*d5i*dotj*ddz-idmz(js,j)*d3i) - dotj*d3i*df2*d1i*ddz

                   efxt(is,i) = efxt(is,i) + efldx
                   efyt(is,i) = efyt(is,i) + efldy
                   efzt(is,i) = efzt(is,i) + efldz

                   ! test OI
                   lshefxt(is,i) = lshefxt(is,i) + efldx + (f2-1.d0)*qjs*d3i*ddx
                   lshefyt(is,i) = lshefyt(is,i) + efldy + (f2-1.d0)*qjs*d3i*ddy
                   lshefzt(is,i) = lshefzt(is,i) + efldz + (f2-1.d0)*qjs*d3i*ddz
                   ! test OI

                   ! virial contribution ! OI
                   lqvirial = lqvirial - qis*(efldx*ddx + efldy*ddy + efldz*ddz)

                  ! gradient of efield of permanent charges w.r.t position of
                  ! pol. centers

                   ef0gxx = (f2-1.d0)*qjs*(-3.d0*d5i*ddx*ddx+d3i) + qjs*ddx*d3i*df2*ddx*d1i 
                   ef0gyy = (f2-1.d0)*qjs*(-3.d0*d5i*ddy*ddy+d3i) + qjs*ddy*d3i*df2*ddy*d1i
                   ef0gzz = (f2-1.d0)*qjs*(-3.d0*d5i*ddz*ddz+d3i) + qjs*ddz*d3i*df2*ddz*d1i
                   ef0gxy = (f2-1.d0)*qjs*(-3.d0*d5i*ddx*ddy)     + qjs*ddy*d3i*df2*ddx*d1i
                   ef0gyx = (f2-1.d0)*qjs*(-3.d0*d5i*ddy*ddx)     + qjs*ddx*d3i*df2*ddy*d1i
                   ef0gxz = (f2-1.d0)*qjs*(-3.d0*d5i*ddx*ddz)     + qjs*ddz*d3i*df2*ddx*d1i
                   ef0gzx = (f2-1.d0)*qjs*(-3.d0*d5i*ddz*ddx)     + qjs*ddx*d3i*df2*ddz*d1i
                   ef0gyz = (f2-1.d0)*qjs*(-3.d0*d5i*ddy*ddz)     + qjs*ddz*d3i*df2*ddy*d1i
                   ef0gzy = (f2-1.d0)*qjs*(-3.d0*d5i*ddz*ddy)     + qjs*ddy*d3i*df2*ddz*d1i

                  ! gradient of efield of induced dipoles w.r.t position of pol.
                  ! centers

                   efigxx =(f3-1.d0)*(-15.d0*dotj*d7i*ddx*ddx+3.d0*dotj*d5i+2.d0*3.d0*d5i*idmx(js,j)*ddx) &
                            + cst*(3.d0*ddx*dotj*d5i-idmx(js,j)*d3i)*df3*ddx*d1i

                   efigyy =(f3-1.d0)*(-15.d0*dotj*d7i*ddy*ddy+3.d0*dotj*d5i+2.d0*3.d0*d5i*idmy(js,j)*ddy) &
                            + cst*(3.d0*ddy*dotj*d5i-idmy(js,j)*d3i)*df3*ddy*d1i

                   efigzz =(f3-1.d0)*(-15.d0*dotj*d7i*ddz*ddz+3.d0*dotj*d5i+2.d0*3.d0*d5i*idmz(js,j)*ddz) &
                            + cst*(3.d0*ddz*dotj*d5i-idmz(js,j)*d3i)*df3*ddz*d1i

                   efigxy =(f3-1.d0)*(-15.d0*dotj*d7i*ddx*ddy+3.d0*d5i*idmx(js,j)*ddy+3.d0*d5i*idmy(js,j)*ddx) &
                            + cst*(3.d0*ddy*dotj*d5i-idmy(js,j)*d3i)*df3*ddx*d1i 

                   efigyx =(f3-1.d0)*(-15.d0*dotj*d7i*ddy*ddx+3.d0*d5i*idmy(js,j)*ddx+3.d0*d5i*idmx(js,j)*ddy) &
                            + cst*(3.d0*ddx*dotj*d5i-idmx(js,j)*d3i)*df3*ddy*d1i

                   efigxz =(f3-1.d0)*(-15.d0*dotj*d7i*ddx*ddz+3.d0*d5i*idmx(js,j)*ddz+3.d0*d5i*idmz(js,j)*ddx) &
                            + cst*(3.d0*ddz*dotj*d5i-idmz(js,j)*d3i)*df3*ddx*d1i

                   efigzx =(f3-1.d0)*(-15.d0*dotj*d7i*ddz*ddx+3.d0*d5i*idmz(js,j)*ddx+3.d0*d5i*idmx(js,j)*ddz) &
                            + cst*(3.d0*ddx*dotj*d5i-idmx(js,j)*d3i)*df3*ddz*d1i

                   efigyz =(f3-1.d0)*(-15.d0*dotj*d7i*ddy*ddz+3.d0*d5i*idmy(js,j)*ddz+3.d0*d5i*idmz(js,j)*ddy) &
                            + cst*(3.d0*ddz*dotj*d5i-idmz(js,j)*d3i)*df3*ddy*d1i 

                   efigzy =(f3-1.d0)*(-15.d0*dotj*d7i*ddz*ddy+3.d0*d5i*idmz(js,j)*ddy+3.d0*d5i*idmy(js,j)*ddz) &
                            + cst*(3.d0*ddy*dotj*d5i-idmy(js,j)*d3i)*df3*ddz*d1i 

                   efgxx(is,i) = efgxx(is,i) + ef0gxx +efigxx
                   efgyy(is,i) = efgyy(is,i) + ef0gyy +efigyy
                   efgzz(is,i) = efgzz(is,i) + ef0gzz +efigzz
                   efgxy(is,i) = efgxy(is,i) + ef0gxy +efigxy
                   efgyx(is,i) = efgyx(is,i) + ef0gyx +efigyx
                   efgxz(is,i) = efgxz(is,i) + ef0gxz +efigxz
                   efgzx(is,i) = efgzx(is,i) + ef0gzx +efigzx
                   efgyz(is,i) = efgyz(is,i) + ef0gyz +efigyz
                   efgzy(is,i) = efgzy(is,i) + ef0gzy +efigzy

                   ! virial contribution  ! OI

                   efi0gxx = ef0gxx + efigxx
                   efi0gyy = ef0gyy + efigyy
                   efi0gzz = ef0gzz + efigzz
                   efi0gxy = ef0gxy + efigxy
                   efi0gyx = ef0gyx + efigyx
                   efi0gxz = ef0gxz + efigxz
                   efi0gzx = ef0gzx + efigzx
                   efi0gyz = ef0gyz + efigyz
                   efi0gzy = ef0gzy + efigzy

                   lqvirial=lqvirial-( idmx(is,i)*efi0gxx*ddx+idmy(is,i)*efi0gxy*ddx+idmz(is,i)*efi0gxz*ddx &
                                      +idmx(is,i)*efi0gyx*ddy+idmy(is,i)*efi0gyy*ddy+idmz(is,i)*efi0gyz*ddy &
                                      +idmx(is,i)*efi0gzx*ddz+idmy(is,i)*efi0gzy*ddz+idmz(is,i)*efi0gzz*ddz)

                    ! gradient of efield of permanent charges w.r.t position of
                    ! charge sites

                   efldx = (f2-1.d0)*(3.d0*d5i*doti*ddx-idmx(is,i)*d3i) - doti*d3i*df2*d1i*ddx 
                   efldy = (f2-1.d0)*(3.d0*d5i*doti*ddy-idmy(is,i)*d3i) - doti*d3i*df2*d1i*ddy
                   efldz = (f2-1.d0)*(3.d0*d5i*doti*ddz-idmz(is,i)*d3i) - doti*d3i*df2*d1i*ddz

                   efxt(js,j) = efxt(js,j) + efldx
                   efyt(js,j) = efyt(js,j) + efldy
                   efzt(js,j) = efzt(js,j) + efldz

                   ! test OI
                   lshefxt(js,j) = lshefxt(js,j) + efldx - (f2-1.d0)*qis*d3i*ddx
                   lshefyt(js,j) = lshefyt(js,j) + efldy - (f2-1.d0)*qis*d3i*ddy
                   lshefzt(js,j) = lshefzt(js,j) + efldz - (f2-1.d0)*qis*d3i*ddz
                   ! test OI

                   ! gradient of efield of permanent charges in ith molecule at
                   ! is site

                   ef0gxx = (f2-1.d0)*qis*(-3.d0*d5i*ddx*ddx + d3i) + qis*ddx*d3i*df2*ddx*d1i
                   ef0gyy = (f2-1.d0)*qis*(-3.d0*d5i*ddy*ddy + d3i) + qis*ddy*d3i*df2*ddy*d1i
                   ef0gzz = (f2-1.d0)*qis*(-3.d0*d5i*ddz*ddz + d3i) + qis*ddz*d3i*df2*ddz*d1i
                   ef0gxy = (f2-1.d0)*qis*(-3.d0*d5i*ddx*ddy)       + qis*ddy*d3i*df2*ddx*d1i
                   ef0gyx = (f2-1.d0)*qis*(-3.d0*d5i*ddy*ddx)       + qis*ddx*d3i*df2*ddy*d1i 
                   ef0gxz = (f2-1.d0)*qis*(-3.d0*d5i*ddx*ddz)       + qis*ddz*d3i*df2*ddx*d1i
                   ef0gzx = (f2-1.d0)*qis*(-3.d0*d5i*ddz*ddx)       + qis*ddx*d3i*df2*ddz*d1i 
                   ef0gyz = (f2-1.d0)*qis*(-3.d0*d5i*ddy*ddz)       + qis*ddz*d3i*df2*ddy*d1i
                   ef0gzy = (f2-1.d0)*qis*(-3.d0*d5i*ddz*ddy)       + qis*ddy*d3i*df2*ddz*d1i

                   ! gradient of efield of induced dipoles in ith molecule at is
                   ! site   !!! check sign

                   efigxx = -(f3-1.d0)*(-15.d0*doti*d7i*ddx*ddx+3.d0*doti*d5i+2.d0*3.d0*d5i*idmx(is,i)*ddx) &
                              - cst*(3.d0*ddx*doti*d5i-idmx(is,i)*d3i)*df3*ddx*d1i

                   efigyy = -(f3-1.d0)*(-15.d0*doti*d7i*ddy*ddy+3.d0*doti*d5i+2.d0*3.d0*d5i*idmy(is,i)*ddy) &
                              - cst*(3.d0*ddy*doti*d5i-idmy(is,i)*d3i)*df3*ddy*d1i
  
                   efigzz = -(f3-1.d0)*(-15.d0*doti*d7i*ddz*ddz+3.d0*doti*d5i+2.d0*3.d0*d5i*idmz(is,i)*ddz) &
                              - cst*(3.d0*ddz*doti*d5i-idmz(is,i)*d3i)*df3*ddz*d1i

                   efigxy = -(f3-1.d0)*(-15.d0*doti*d7i*ddx*ddy+3.d0*d5i*idmx(is,i)*ddy+3.d0*d5i*idmy(is,i)*ddx) &
                              - cst*(3.d0*ddy*doti*d5i-idmy(is,i)*d3i)*df3*ddx*d1i

                   efigyx = -(f3-1.d0)*(-15.d0*doti*d7i*ddy*ddx+3.d0*d5i*idmy(is,i)*ddx+3.d0*d5i*idmx(is,i)*ddy) &
                              - cst*(3.d0*ddx*doti*d5i-idmx(is,i)*d3i)*df3*ddy*d1i
                 
                   efigxz = -(f3-1.d0)*(-15.d0*doti*d7i*ddx*ddz+3.d0*d5i*idmx(is,i)*ddz+3.d0*d5i*idmz(is,i)*ddx) &
                              - cst*(3.d0*ddz*doti*d5i-idmz(is,i)*d3i)*df3*ddx*d1i

                   efigzx = -(f3-1.d0)*(-15.d0*doti*d7i*ddz*ddx+3.d0*d5i*idmz(is,i)*ddx+3.d0*d5i*idmx(is,i)*ddz) &
                              - cst*(3.d0*ddx*doti*d5i-idmx(is,i)*d3i)*df3*ddz*d1i 

                   efigyz = -(f3-1.d0)*(-15.d0*doti*d7i*ddy*ddz+3.d0*d5i*idmy(is,i)*ddz+3.d0*d5i*idmz(is,i)*ddy) &
                              - cst*(3.d0*ddz*doti*d5i-idmz(is,i)*d3i)*df3*ddy*d1i

                   efigzy = -(f3-1.d0)*(-15.d0*doti*d7i*ddz*ddy+3.d0*d5i*idmz(is,i)*ddy+3.d0*d5i*idmy(is,i)*ddz) &
                              - cst*(3.d0*ddy*doti*d5i-idmy(is,i)*d3i)*df3*ddz*d1i


                   !write(*,*) (f3-1.d0)*(-15.d0*doti*d7i*ddz*ddy+3.d0*d5i*idmz(is,i)*ddy+3.d0*d5i*idmy(is,i)*ddz), & 
                   !           (3.d0*ddy*doti*d5i-idmy(is,i)*d3i)*df3*ddz*d1i 

                   efgxx(js,j) = efgxx(js,j) + ef0gxx +efigxx
                   efgyy(js,j) = efgyy(js,j) + ef0gyy +efigyy
                   efgzz(js,j) = efgzz(js,j) + ef0gzz +efigzz
                   efgxy(js,j) = efgxy(js,j) + ef0gxy +efigxy
                   efgyx(js,j) = efgyx(js,j) + ef0gyx +efigyx
                   efgxz(js,j) = efgxz(js,j) + ef0gxz +efigxz
                   efgzx(js,j) = efgzx(js,j) + ef0gzx +efigzx
                   efgyz(js,j) = efgyz(js,j) + ef0gyz +efigyz
                   efgzy(js,j) = efgzy(js,j) + ef0gzy +efigzy

!                   end if   ! site off ends

                end do         ! loop over sites of molecule j  ends

             end do        ! loop over sites of molecule i  ends

          end if    ! com-com cutoff check ends

        end do ! loop over jth molecule ends
		!$omp end do
	   endif
    end do   ! loop over ith molecule ends
	if (MPIrank .EQ. 0) then
		do i = 1, MPIcomm_size - 1
			call MPI_RECV(MPIaddr, 1, MPI_INTEGER, MPI_ANY_SOURCE, MPI_TAG_CommReady, MPI_COMM_WORLD, MPIstat, MPIierr);
			call MPI_SEND(nm + 1, 1, MPI_INTEGER, MPIaddr, MPI_TAG_SendData, MPI_COMM_WORLD, MPIierr);
		enddo
	endif
!$omp barrier
	! DBG
	deltaTime = SECNDS(iniTime);
!	write(dbgUnit, *) 'P#', MPIrank, '. shortrange_forces_damped calc time	', deltaTime;
	! DBG

!$omp critical
    qvirial = qvirial + lqvirial

    do i=1,nm
       do is=1,ns
          sefgxx(is,i) = sefgxx(is,i) + efgxx(is,i)
          sefgyy(is,i) = sefgyy(is,i) + efgyy(is,i)
          sefgzz(is,i) = sefgzz(is,i) + efgzz(is,i)
          sefgxy(is,i) = sefgxy(is,i) + efgxy(is,i)
          sefgyx(is,i) = sefgyx(is,i) + efgyx(is,i)
          sefgxz(is,i) = sefgxz(is,i) + efgxz(is,i)
          sefgzx(is,i) = sefgzx(is,i) + efgzx(is,i)
          sefgyz(is,i) = sefgyz(is,i) + efgyz(is,i)
          sefgzy(is,i) = sefgzy(is,i) + efgzy(is,i)
          sefxt(is,i)  = sefxt(is,i)  + efxt(is,i)
          sefyt(is,i)  = sefyt(is,i)  + efyt(is,i)
          sefzt(is,i)  = sefzt(is,i)  + efzt(is,i)

          llshefxt(is,i) = llshefxt(is,i) + lshefxt(is,i) !*r4pie0
          llshefyt(is,i) = llshefyt(is,i) + lshefyt(is,i) !*r4pie0
          llshefzt(is,i) = llshefzt(is,i) + lshefzt(is,i) !*r4pie0

       end do
    end do
!$omp end critical
!$omp barrier
!                               tot    tot
! ... calculate forces: f  = q e    + u    efg             !!! check this
!                        a      a      b      ab

!$omp do
    do i=1,nm
       do is=1,ns

         qis = chgs(is,i)

!         fxss = qis*sefxt(is,i) + idmx(is,i)*sefgxx(is,i) + idmy(is,i)*sefgyx(is,i) + idmz(is,i)*sefgzx(is,i)
!         fyss = qis*sefyt(is,i) + idmx(is,i)*sefgxy(is,i) + idmy(is,i)*sefgyy(is,i) + idmz(is,i)*sefgzy(is,i)
!         fzss = qis*sefzt(is,i) + idmx(is,i)*sefgxz(is,i) + idmy(is,i)*sefgyz(is,i) + idmz(is,i)*sefgzz(is,i)

         fxss = qis*sefxt(is,i) + idmx(is,i)*sefgxx(is,i) + idmy(is,i)*sefgxy(is,i) + idmz(is,i)*sefgxz(is,i)
         fyss = qis*sefyt(is,i) + idmx(is,i)*sefgyx(is,i) + idmy(is,i)*sefgyy(is,i) + idmz(is,i)*sefgyz(is,i)
         fzss = qis*sefzt(is,i) + idmx(is,i)*sefgzx(is,i) + idmy(is,i)*sefgzy(is,i) + idmz(is,i)*sefgzz(is,i)

         llfxs(is,i)= llfxs(is,i) + fxss*r4pie0
         llfys(is,i)= llfys(is,i) + fyss*r4pie0
         llfzs(is,i)= llfzs(is,i) + fzss*r4pie0

       end do
    end do
!$omp end do

    deallocate(efgxx)
    deallocate(efgyy)
    deallocate(efgzz)
    deallocate(efgxy)
    deallocate(efgyx)
    deallocate(efgxz)
    deallocate(efgzx)
    deallocate(efgyz)
    deallocate(efgzy)
    deallocate(efxt)
    deallocate(efyt)
    deallocate(efzt)
!$omp end parallel

!******************************************************
! ************* Now calculating forces ****************    ends here
    deallocate(sefgxx)
    deallocate(sefgyy)
    deallocate(sefgzz)
    deallocate(sefgxy)
    deallocate(sefgyx)
    deallocate(sefgxz)
    deallocate(sefgzx)
    deallocate(sefgyz)
    deallocate(sefgzy)
    deallocate(sefxt)
    deallocate(sefyt)
    deallocate(sefzt)

! MPI sync
	j = nm * ns * 3 + nm * ns * 3 + 1; !fxs, fys, fzs, qvirial
	allocate(MPIbuff1(0:j - 1));
	allocate(MPIbuff2(0:j - 1));
	do k = 0, j - 1
		MPIbuff2(k) = 0.0D0;
	enddo
	k = 0;
	do i = 1, nm
		do is = 1, ns
			MPIbuff1(k) = llfxs(is, i); k = k + 1;
			MPIbuff1(k) = llfys(is, i); k = k + 1;
			MPIbuff1(k) = llfzs(is, i); k = k + 1;

                        MPIbuff1(k) = llshefxt(is, i); k = k + 1;
                        MPIbuff1(k) = llshefyt(is, i); k = k + 1;
                        MPIbuff1(k) = llshefzt(is, i); k = k + 1;

		enddo
	enddo
        MPIbuff1(k) = qvirial; k = k + 1;   ! OI
	call MPI_ALLREDUCE(MPIbuff1, MPIbuff2, j, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, MPIierr);
	k = 0;
	do i = 1, nm
		do is = 1, ns
			fxs(is, i) = fxs(is, i) + MPIbuff2(k); k = k + 1;
			fys(is, i) = fys(is, i) + MPIbuff2(k); k = k + 1;
			fzs(is, i) = fzs(is, i) + MPIbuff2(k); k = k + 1;

                        Etotx(is, i) = Etotx(is, i) + MPIbuff2(k); k = k + 1;
                        Etoty(is, i) = Etoty(is, i) + MPIbuff2(k); k = k + 1;
                        Etotz(is, i) = Etotz(is, i) + MPIbuff2(k); k = k + 1;

		enddo
	enddo
        qvirial = MPIbuff2(k); k = k + 1;   ! OI
	deallocate(MPIbuff1);
	deallocate(MPIbuff2);
! MPI sync
    deallocate(llfxs);
    deallocate(llfys);
    deallocate(llfzs);
	! DBG
	! if (MPIrank .EQ. 0) then
		! write(342, '(a)') '#';
		! write(342, '(a)') 'fxs';
		! write(342, '(1p3E16.6E3)') fxs;
		! write(342, '(a)') 'fys';
		! write(342, '(1p3E16.6E3)') fys;
		! write(342, '(a)') 'fzs';
		! write(342, '(1p3E16.6E3)') fzs;
	! endif
	! DBG

        vir_shtind = r4pie0*qvirial
        vir = vir + vir_shtind

	return;
  end subroutine shortrange_forces_damped

end module nb_induction_model

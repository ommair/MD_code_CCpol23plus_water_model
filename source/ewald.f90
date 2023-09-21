module ewald

  use, intrinsic :: iso_fortran_env
  use common_arrays

  implicit none
  private

  public :: ewald_setup,ewald_self_correction,rspace_ewald
  public :: kspace_ewald,kspace_ewald2,coul_pot


contains

! compute alpha and the kmax

  subroutine ewald_setup(ewps,alp,k1,k2,k3,boxi) 
! (input) ewps is precision of ewald sum as input
! (output) alp is alpha value calculted using ewps and rcut
! (output) max number reciprocal space vectors(k1,k2,k3) required

    implicit none

    real(8) :: boxi
    real(8) :: tol,tol1
    real(8) :: ewps, eprec
    real(8) :: alp
    integer :: k1,k2,k3  
    real(8) :: fac   
   
    eprec = min(abs(ewps),0.5d0)
    tol=sqrt(abs(log(eprec*rcut)))
    alp=sqrt(abs(log(eprec*rcut*tol)))/rcut
    tol1=sqrt(-log(eprec*rcut*(2.d0*tol*alp)**2))

    fac = 1.d0 ! for cubic boundary

    k1=nint(0.25d0+fac*boxi*alp*tol1/pi)
    k2=nint(0.25d0+fac*boxi*alp*tol1/pi)
    k3=nint(0.25d0+fac*boxi*alp*tol1/pi)    
    

  end subroutine ewald_setup

  subroutine ewald_self_correction(alp,boxi,qscpe)
    implicit none

    !real(8), parameter :: twosqpi= 2.0d0/sqrt(pi)
    integer :: i,j,is,js
    real(8) :: qscpe,qvir_selfp,qvir_selfi,qvir_self,qvirial
    real(8) :: e_selfp ! Point self-energy (of 1 molecule)
    real(8) :: e_selfi ! Intramolecular self-energy (of 1 molecule)
    real(8) :: sumq,sumq2,sumff,lsumff,lqvirial
    real(8) :: alp,boxi,volmi
    real(8) ::  chgprd  ! charge product
    real(8) :: drx,dry,drz,ddrx,ddry,ddrz
    real(8) :: rrx,rry,rrz,drr
    real(8) :: sfij,fs     ! force form self energy term (point plus intramolecular)
!    real(8) :: sfx(nsite),sfy(nsite),sfz(nsite)
    real(8) :: fxsij,fysij,fzsij

    qvirial = 0.d0
! calculating point self energy for single molecule

    e_selfp =0.d0
    sumq=0.d0
    sumq2=0.d0
    volmi = boxi**3

    do i=1,nm
    do is=1, ns
!       if(charge(is).ne.0.d0) then

          sumq = sumq + charge(is)
          sumq2 = sumq2 + charge(is)*charge(is)  
       
!       end if                      
    end do     
    end do

!   background term for non-neutral systems
    
    e_selfp  =  -pi * (sumq/alp)**2 / (2.0*volmi)

    ! self energy term
    ! q**2*alpha/((4 pi eps0)*(sqrt(pi)) see dlpoly classic manual(page 47, eq 2.196)

    e_selfp = e_selfp - alp*sumq2/sqrt(pi)

!    qvir_selfp =  - e_selfp

! calculating intra-molecular corrections for energy for single molecule

     sumff = 0.d0

!$omp parallel DEFAULT(SHARED)&
!$omp& private(i,is,js,chgprd,drx,dry,drz,drr,fs)&
!$omp& private(fxsij,fysij,fzsij,lqvirial,lsumff)

     lqvirial = 0.0D0
     lsumff = 0.0D0

!$omp do schedule(dynamic)
    do i=1,nm
       do is=1,ns-1
          do js=is+1,ns
             chgprd = charge(is)*charge(js)
             if(chgprd.ne.0.d0) then

              drx = (x(i)+xs(is,i)) - (x(i)+xs(js,i))
              dry = (y(i)+ys(is,i)) - (y(i)+ys(js,i))
              drz = (z(i)+zs(is,i)) - (z(i)+zs(js,i))

!              drx = xs(is,i) - xs(js,i) 
!              dry = ys(is,i) - ys(js,i) 
!              drz = zs(is,i) - zs(js,i)  

               drr = drx**2+dry**2+drz**2
               drr = sqrt(drr)
      
               lsumff = lsumff + erf(alp*drr)*chgprd/drr
 
               fs = -r4pie0*chgprd*( erf(alp*drr)/drr**3-twosqpi*alp*exp(-(alp*drr)**2)/drr**2 )

               fxsij = drx*fs
               fysij = dry*fs
               fzsij = drz*fs
             
               fxs(is,i)=fxs(is,i) + fxsij 
               fys(is,i)=fys(is,i) + fysij 
               fzs(is,i)=fzs(is,i) + fzsij

               fxs(js,i)=fxs(js,i) - fxsij
               fys(js,i)=fys(js,i) - fysij
               fzs(js,i)=fzs(js,i) - fzsij 

               ! virial = - rab * fab 

               lqvirial = lqvirial -(drx*fxsij + dry*fysij + drz*fzsij)
!               qvirial = qvirial - fs*drr**2

             end if
          end do
       end do
    end do
!$omp end do

!$omp atomic
   sumff = sumff + lsumff
!$omp atomic
   qvirial = qvirial + lqvirial

!$omp end parallel
 
    e_selfi = -sumff

!    qvir_selfi = -e_selfi 

    qscpe = r4pie0*(e_selfp + e_selfi)

    qscpe = qscpe/418.4d0  ! convert to kcal from internal units

!    qvir_self = -qscpe*418.4d0 !qvirial
    qvir_self = qvirial

    ! accumalation of virial

    vir = vir + qvir_self  ! in internal units

!!    vir = vir + qvir_selfi+qvir_selfp

!    write(*,*)'qvir_self ',qvir_self

  end subroutine ewald_self_correction


  subroutine rspace_ewald(xx,yy,zz,xxs,yys,zzs,boxi,alp,qrcpe)
    implicit none
    real(8) :: xx(nom),yy(nom),zz(nom)
    real(8) :: xxs(nsite,nom),yys(nsite,nom),zzs(nsite,nom)
    real(8) :: boxi,alp,qrcpe

    integer :: i,j,is,js,k
    real(8) :: qvir_rc
    real(8) :: wij
    real(8) :: xcc,ycc,zcc
    real(8) :: dx,dy,dz
    real(8) :: rsx,rsy,rsz
    real(8) :: rccsq,rssq,drr
    real(8) :: qis,qjs,chgprd
    real(8) :: fc,cfij
    real(8) :: fxsij,fysij,fzsij

    real(8) :: qvirial,lqrcpe,lqvirial
    ! Shared
    real(8), allocatable :: lfxs(:, :), lfys(:, :), lfzs(:, :)
	! MPI
	INTEGER(4) MPIchunk_size_loc, iii;
	real(8), allocatable :: MPIbuff1(:), MPIbuff2(:);
    real(8), allocatable :: llfxs(:, :), llfys(:, :), llfzs(:, :);
	! DBG
	real(4) iniTime, deltaTime;
	! DBG
	! MPI

    qvirial = 0.d0
    qrcpe =0.d0

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

!$omp parallel DEFAULT(SHARED)&
!$omp& private(i,iii,j,is,js,k)&
!$omp& private(qvir_rc,wij,xcc,ycc,zcc,dx,dy,dz,rsx,rsy,rsz,rccsq,rssq)&
!$omp& private(drr,qis,qjs,chgprd,fc,cfij,fxsij,fysij,fzsij,lfxs,lfys,lfzs,lqrcpe,lqvirial)

    lqrcpe = 0.0D0
    lqvirial = 0.0D0

    allocate(lfxs(ns,nm))
    allocate(lfys(ns,nm))
    allocate(lfzs(ns,nm))

    do j=1,nm
       do js=1,ns
          lfxs(js,j) = 0.0D0
          lfys(js,j) = 0.0D0
          lfzs(js,j) = 0.0D0
       end do
    end do

	! DBG
	iniTime = SECNDS(0.0);
	! DBG
!     do i=1,nm-1
!	  do i = MPIrank + 1, nm - 1, MPIcomm_size
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
        do j=i+1, nm

           dx = xx(i)-xx(j)
           dy = yy(i)-yy(j)
           dz = zz(i)-zz(j)

           xcc = dx - boxi*nint(dx/boxi)
           ycc = dy - boxi*nint(dy/boxi)
           zcc = dz - boxi*nint(dz/boxi)

           rccsq=xcc**2 + ycc**2 + zcc**2

           if (rccsq.lt.rcutsd) then   ! com cut off starts

              do is=1,ns

                  do js=1,ns
                 
                  qis = chgs(is,i)
                  qjs = chgs(js,j) 
                  
                  rsx = xcc + xxs(is,i)-xxs(js,j)
                  rsy = ycc + yys(is,i)-yys(js,j)
                  rsz = zcc + zzs(is,i)-zzs(js,j) 
                  
                  rssq = rsx**2 + rsy**2 + rsz**2

!                  if (rssq.le.rcutsd) then       ! site cut off starts 
                          
                  drr = sqrt(rssq)

                  lqrcpe = lqrcpe + r4pie0*qis*qjs*erfc(alp*drr)/drr

                  fc = r4pie0*qis*qjs*(erfc(alp*drr)/drr+twosqpi*alp*exp(-(alp*drr)**2))

                  cfij = fc/rssq

                  fxsij = rsx*cfij
                  fysij = rsy*cfij
                  fzsij = rsz*cfij 
      
                  lfxs(is,i)=lfxs(is,i)+fxsij
                  lfys(is,i)=lfys(is,i)+fysij
                  lfzs(is,i)=lfzs(is,i)+fzsij

                  lfxs(js,j)=lfxs(js,j)-fxsij
                  lfys(js,j)=lfys(js,j)-fysij
                  lfzs(js,j)=lfzs(js,j)-fzsij

                  ! virial = - rab * fab

                  lqvirial = lqvirial -(rsx*fxsij + rsy*fysij + rsz*fzsij) 
         
!                  end if    ! site cut off ends
 
                 end do 

              end do
              
           end if         ! com cut off ends
           
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
	! DBG
	deltaTime = SECNDS(iniTime);
!	write(dbgUnit, *) 'P#', MPIrank, '. rspace_ewald calc time	', deltaTime;
	! DBG
!$omp critical
     qvirial = qvirial + lqvirial
     qrcpe = qrcpe + lqrcpe

     do j=1,nm
        do js=1,ns
           llfxs(js,j) = llfxs(js,j) + lfxs(js,j)
           llfys(js,j) = llfys(js,j) + lfys(js,j)
           llfzs(js,j) = llfzs(js,j) + lfzs(js,j)
        end do
     end do
!$omp end critical

    deallocate(lfxs)
    deallocate(lfys)
    deallocate(lfzs)
!$omp end parallel
! MPI sync
	j = nm * ns * 3 + 2; !fxs, fys, fzs, qvirial, qrcpe
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
	MPIbuff1(k) = qvirial; k = k + 1;
	MPIbuff1(k) = qrcpe; k = k + 1;
	call MPI_ALLREDUCE(MPIbuff1, MPIbuff2, j, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, MPIierr);
	k = 0;
	do i = 1, nm
		do is = 1, ns
			fxs(is, i) = fxs(is, i) + MPIbuff2(k); k = k + 1;
			fys(is, i) = fys(is, i) + MPIbuff2(k); k = k + 1;
			fzs(is, i) = fzs(is, i) + MPIbuff2(k); k = k + 1;
		enddo
	enddo
	qvirial = MPIbuff2(k); k = k + 1;
	qrcpe = MPIbuff2(k); k = k + 1;
	deallocate(MPIbuff1);
	deallocate(MPIbuff2);
! MPI sync
    deallocate(llfxs);
    deallocate(llfys);
    deallocate(llfzs);

!     qvir_rc = -qrcpe

     qrcpe = qrcpe/418.4d0       ! real space energy for coulomb in internal units

!     qvir_rc = -qrcpe*418.4d0 !qvirial
     qvir_rc = qvirial

     vir = vir + qvir_rc

!!     vir = vir - qrcpe*418.4d0 

!     write(*,*)'qvir_rc ',qvir_rc

  end subroutine rspace_ewald


  subroutine kspace_ewald2(xx,yy,zz,xxs,yys,zzs,lmax,mmax,nmax,alp,boxi,qkcpe)
    implicit none
    real(8) :: xx(nom),yy(nom),zz(nom)
    real(8) :: xxs(nsite,nom),yys(nsite,nom),zzs(nsite,nom)
    integer :: lmax,mmax,nmax 
    real(8) :: alp,boxi,qkcpe

!    integer, parameter :: nqmax = nsite*nom  
    integer :: nqmax !ns*nm
    real(8), parameter :: twopi = 2.d0*pi

    integer :: i,j,is,js
    integer :: ll,mm,nn,kk,lmnp
    integer :: mmin,nmin

    integer :: llim,mlim,nlim
    integer :: klmmax,klmmaxsq 
    integer :: noq

    real(8) :: qvirial,lqvirial 
    real(8) :: qpe,lqpe,rl,rm,rn,rksq,qforce

    real(8) :: volmi 
    real(8) :: qvir_kc
!    real(8) :: qi(nqmax),xi(nqmax),yi(nqmax),zi(nqmax)
    real(8), allocatable :: qi(:),xi(:),yi(:),zi(:)
!    real(8) :: a(ksqmax)
    real(8), allocatable :: a(:)
    complex(8) :: sumqex
!    complex(8) :: expikr(nqmax)
    complex(8), allocatable :: expikr(:)
!    complex(8) :: el(nqmax,0:20),em(nqmax,-20:20),en(nqmax,-20:20)
    complex(8), allocatable :: el(:, :),em(:, :),en(:, :)

   ! Shared
    real(8), allocatable :: lfxs(:, :), lfys(:, :), lfzs(:, :)
	! MPI
	INTEGER(4) MPIchunk_size_loc, iii, k;
	real(8), allocatable :: MPIbuff1(:), MPIbuff2(:);
    real(8), allocatable :: llfxs(:, :), llfys(:, :), llfzs(:, :);
	! DBG
	real(4) iniTime, deltaTime;
	! DBG
	! MPI

    nqmax = ns*nm

    allocate(qi(nqmax))
    allocate(xi(nqmax))
    allocate(yi(nqmax))
    allocate(zi(nqmax))
    allocate(el(nqmax,0:20))
    allocate(em(nqmax,-20:20))
    allocate(en(nqmax,-20:20))

    klmmax = min(dble(lmax),dble(mmax),dble(nmax))
    klmmaxsq = klmmax**2

    llim = lmax 
    mlim = mmax
    nlim = nmax

    qpe = 0.d0
    qvirial = 0.d0

    volmi = boxi**3
    
!    collect quantities for charged sites
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

!$omp parallel DEFAULT(SHARED)&
!$omp& private(i,iii,j,is,js,ll,mm,nn,kk,lmnp,mmin,nmin)&
!$omp& private(rl,rm,rn,rksq,qforce,sumqex,lqpe,lqvirial)&
!$omp& private(a,expikr,lfxs,lfys,lfzs)

    lqpe = 0.0D0
    lqvirial = 0.0D0

    allocate(a(ksqmax))
    allocate(expikr(nqmax))
    allocate(lfxs(ns,nm))
    allocate(lfys(ns,nm))
    allocate(lfzs(ns,nm))

    do j=1,nm
       do js=1,ns
          lfxs(js,j) = 0.0D0
          lfys(js,j) = 0.0D0
          lfzs(js,j) = 0.0D0
       end do
    end do

!$OMP MASTER
      i = 0
      do j = 1,nm
         do js = 1,ns
            if (chgs(js,j).ne.0.0) then
               i = i+1
               qi(i) = chgs(js,j)
               xi(i) = xx(j) + xxs(js,j)
               yi(i) = yy(j) + yys(js,j)
               zi(i) = zz(j) + zzs(js,j)
            endif
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
 
	   !$omp do schedule(dynamic)
       do mm=mmin,mlim

          rm = twopi*mm/boxi

          do nn=nmin,nlim
 
             rn = twopi*nn/boxi

             kk=ll**2+mm**2+nn**2     
             lmnp = ll+mm+nn 

             if ((kk.le.klmmaxsq) .and. .not.(mod(lmnp,2).ne.0) ) then
!             if ((kk.le.klmmaxsq)) then

                rksq = rl*rl + rm*rm + rn*rn

                a(kk) = twopi/volmi*exp(-0.25d0*rksq/alp**2)/rksq             

            
                ! form expikr for each charge
                do i = 1,noq
                   expikr(i) = el(i,ll)*em(i,mm)*en(i,nn)
                enddo      

             !form sum of qi*expikr

                sumqex = (0.0,0.0)
                do i = 1,noq
                   sumqex = sumqex+qi(i)*expikr(i)
                enddo

               ! accumulate potential energy

               lqpe = lqpe + a(kk)*conjg(sumqex)*sumqex
               lqvirial = lqvirial - (1.d0-0.5d0*rksq/alp**2)*a(kk)*conjg(sumqex)*sumqex

               i=0
               do j=1,nm
                  do js=1,ns
                     if (chgs(js,j).ne.0.0) then
                        i=i+1
                        qforce = -4.0*r4pie0*a(kk)*qi(i)*aimag(sumqex*conjg(expikr(i)))
                      
                        lfxs(js,j) = lfxs(js,j) + rl*qforce
                        lfys(js,j) = lfys(js,j) + rm*qforce
                        lfzs(js,j) = lfzs(js,j) + rn*qforce

!                        lqvirial = lqvirial - (rl*rl*qforce+rm*rm*qforce+rn*rn*qforce)

!                        write(*,*)qvirial !qforce

                     end if
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
!	write(dbgUnit, *) 'P#', MPIrank, '. kspace_ewald calc time	', deltaTime;
	! DBG
!$omp critical
   qpe = qpe + lqpe
   qvirial = qvirial + lqvirial

   do j=1,nm
      do js=1,ns
         llfxs(js,j) = llfxs(js,j) + lfxs(js,j)
         llfys(js,j) = llfys(js,j) + lfys(js,j)
         llfzs(js,j) = llfzs(js,j) + lfzs(js,j)
      end do
   end do
!$omp end critical

   deallocate(a)
   deallocate(expikr)
   deallocate(lfxs)
   deallocate(lfys)
   deallocate(lfzs)
!$omp end parallel
! MPI sync
	j = nm * ns * 3 + 2; !fxs, fys, fzs, qpe, qvirial
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

        MPIbuff1(k) = qvirial; k = k + 1;
	MPIbuff1(k) = qpe; k = k + 1;
	call MPI_ALLREDUCE(MPIbuff1, MPIbuff2, j, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, MPIierr);
	k = 0;
	do i = 1, nm
		do is = 1, ns
			fxs(is, i) = fxs(is, i) + MPIbuff2(k); k = k + 1;
			fys(is, i) = fys(is, i) + MPIbuff2(k); k = k + 1;
			fzs(is, i) = fzs(is, i) + MPIbuff2(k); k = k + 1;
		enddo
	enddo
        qvirial = MPIbuff2(k); k = k + 1;
	qpe = MPIbuff2(k); k = k + 1;
	deallocate(MPIbuff1);
	deallocate(MPIbuff2);
! MPI sync
    deallocate(llfxs);
    deallocate(llfys);
    deallocate(llfzs);

!    write(*,*) qforce*rl,qforce*rm,qforce*rn

    qkcpe = 2.d0*r4pie0*qpe/418.4d0

!    qvir_kc = -qkcpe*418.4d0 !2.d0*r4pie0*qvirial 
    qvir_kc =2.d0*r4pie0*qvirial

   vir = vir + qvir_kc

!   vir = vir - qkcpe*418.4d0

!   write(*,*) 'qvir_kc ',qvir_kc

   deallocate(qi)
   deallocate(xi)
   deallocate(yi)
   deallocate(zi)
   deallocate(el)
   deallocate(em)
   deallocate(en)

  end subroutine kspace_ewald2



  subroutine kspace_ewald(xx,yy,zz,xxs,yys,zzs,lmax,mmax,nmax,alp,boxi,qkcpe)
    implicit none
    real(8) :: xx(nom),yy(nom),zz(nom)
    real(8) :: xxs(nsite,nom),yys(nsite,nom),zzs(nsite,nom)
    integer :: lmax,mmax,nmax 
    real(8) :: alp,boxi,qkcpe

!    integer, parameter :: nqmax = nsite*nom  
    integer :: nqmax !ns*nm
    real(8), parameter :: twopi = 2.d0*pi

    integer :: i,j,is,js,ii,li,mi,ni,mx
    integer :: ll,mm,nn,kk,lmnp
    integer :: mmin,nmin

    integer :: llim,mlim,nlim
    integer :: klmmax,klmmaxsq 
    integer :: noq

    real(8) :: qvirial,lqvirial 
    real(8) :: qpe,lqpe,rl,rm,rn,rksq,qforce

    real(8) :: volmi 
    real(8) :: qvir_kc
!    real(8) :: qi(nqmax),xi(nqmax),yi(nqmax),zi(nqmax)
    real(8), allocatable :: qi(:),xi(:),yi(:),zi(:)
    real(8), allocatable :: ckr(:),skr(:),clm(:),slm(:),ckc(:),cks(:)
    real(8) :: ckcs,ckss
!    real(8) :: a(ksqmax)
    real(8), allocatable :: a(:)
    complex(8) :: sumqex
!    complex(8) :: expikr(nqmax)
    complex(8), allocatable :: expikr(:)
!    complex(8) :: el(nqmax,0:20),em(nqmax,-20:20),en(nqmax,-20:20)
!    complex(8), allocatable :: el(:, :),em(:, :),en(:, :)
    complex(8), allocatable :: elc(:, :),emc(:, :),enc(:, :)
    complex(8), allocatable :: els(:, :),ems(:, :),ens(:, :)
   ! Shared
    real(8), allocatable :: lfxs(:, :), lfys(:, :), lfzs(:, :)
        ! MPI
        INTEGER(4) MPIchunk_size_loc, iii, k;
        real(8), allocatable :: MPIbuff1(:), MPIbuff2(:);
    real(8), allocatable :: llfxs(:, :), llfys(:, :), llfzs(:, :);
        ! DBG
        real(4) iniTime, deltaTime;
        ! DBG
        ! MPI

    nqmax = ns*nm

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

    klmmax = min(dble(lmax),dble(mmax),dble(nmax))
    klmmaxsq = klmmax**2

    llim = lmax 
    mlim = mmax
    nlim = nmax

    qpe = 0.d0
    qvirial = 0.d0

    volmi = boxi**3
    
!    collect quantities for charged sites
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

!$omp parallel DEFAULT(SHARED)&
!$omp& private(i,iii,j,is,js,ll,mm,nn,kk,lmnp,mmin,nmin)&
!$omp& private(rl,rm,rn,rksq,qforce,sumqex,lqpe,lqvirial)&
!$omp& private(a,expikr,lfxs,lfys,lfzs)

    lqpe = 0.0D0
    lqvirial = 0.0D0

    allocate(a(ksqmax))
!    allocate(expikr(nqmax))
    allocate(lfxs(ns,nm))
    allocate(lfys(ns,nm))
    allocate(lfzs(ns,nm))

    do j=1,nm
       do js=1,ns
          lfxs(js,j) = 0.0D0
          lfys(js,j) = 0.0D0
          lfzs(js,j) = 0.0D0
       end do
    end do

!$OMP MASTER
      i = 0
      do j = 1,nm
         do js = 1,ns
            if (chgs(js,j).ne.0.0) then
               i = i+1
               qi(i) = chgs(js,j)
               xi(i) = xx(j) + xxs(js,j)
               yi(i) = yy(j) + yys(js,j)
               zi(i) = zz(j) + zzs(js,j)
            endif
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
 
           !$omp do schedule(dynamic)
       do mm=mmin,mlim

          mi = abs(mm)

          rm = twopi*mm/boxi

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

                     ckc(ii)=qi(ii)*ckr(ii)
                     cks(ii)=qi(ii)*skr(ii)

                   end do

               else

                   do ii=1,noq

                    !  ckc(ii)=qi(ii)*(clm(ii)*enc(ii,ni)+slm(ii)*ens(ii,ni))
                    !  cks(ii)=qi(ii)*(slm(ii)*enc(ii,ni)-clm(ii)*ens(ii,ni))

                      ckr(ii)=clm(ii)*enc(ii,ni)+slm(ii)*ens(ii,ni)
                      skr(ii)=slm(ii)*enc(ii,ni)-clm(ii)*ens(ii,ni)

                      ckc(ii)=qi(ii)*ckr(ii)
                      cks(ii)=qi(ii)*skr(ii)

                   end do

               end if                                

               ckcs = 0.d0
               ckss = 0.d0
 
               do mx=1,noq

                  ckcs=ckcs+ckc(mx)
                  ckss=ckss+cks(mx)

               end do

  

                ! accumulate potetial energy

               lqpe = lqpe + a(kk)*(ckcs*ckcs+ckss*ckss)

               lqvirial = lqvirial -(1.d0-0.5d0*rksq/alp**2)*a(kk)*(ckcs*ckcs+ckss*ckss)

               i=0
               do j=1,nm
                  do js=1,ns
                     if (chgs(js,j).ne.0.0) then
                        i=i+1
                        qforce = 4.0*r4pie0*a(kk)*(cks(i)*ckcs-ckc(i)*ckss)
                      
                        lfxs(js,j) = lfxs(js,j) + rl*qforce
                        lfys(js,j) = lfys(js,j) + rm*qforce
                        lfzs(js,j) = lfzs(js,j) + rn*qforce

!                        lqvirial = lqvirial -
!                        (rl*rl*qforce+rm*rm*qforce+rn*rn*qforce)

!                        write(*,*)qvirial !qforce

                     end if
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
!       write(dbgUnit, *) 'P#', MPIrank, '. kspace_ewald calc time      ',
!       deltaTime;
        ! DBG
!$omp critical
   qpe = qpe + lqpe
   qvirial = qvirial + lqvirial

   do j=1,nm
      do js=1,ns
         llfxs(js,j) = llfxs(js,j) + lfxs(js,j)
         llfys(js,j) = llfys(js,j) + lfys(js,j)
         llfzs(js,j) = llfzs(js,j) + lfzs(js,j)
      end do
   end do
!$omp end critical

   deallocate(a)
!   deallocate(expikr)
   deallocate(lfxs)
   deallocate(lfys)
   deallocate(lfzs)
!$omp end parallel
! MPI sync
        j = nm * ns * 3 + 2; !fxs, fys, fzs, qpe, qvirial
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

        MPIbuff1(k) = qvirial; k = k + 1;
        MPIbuff1(k) = qpe; k = k + 1;
        call MPI_ALLREDUCE(MPIbuff1, MPIbuff2, j, MPI_DOUBLE_PRECISION, MPI_SUM,MPI_COMM_WORLD, MPIierr);
        k = 0;
        do i = 1, nm
                do is = 1, ns
                        fxs(is, i) = fxs(is, i) + MPIbuff2(k); k = k + 1;
                        fys(is, i) = fys(is, i) + MPIbuff2(k); k = k + 1;
                        fzs(is, i) = fzs(is, i) + MPIbuff2(k); k = k + 1;
                enddo
        enddo
        qvirial = MPIbuff2(k); k = k + 1;
        qpe = MPIbuff2(k); k = k + 1;
        deallocate(MPIbuff1);
        deallocate(MPIbuff2);
! MPI sync
    deallocate(llfxs);
    deallocate(llfys);
    deallocate(llfzs);

!    write(*,*) qforce*rl,qforce*rm,qforce*rn

    qkcpe = 2.d0*r4pie0*qpe/418.4d0

!    qvir_kc = -qkcpe*418.4d0 !2.d0*r4pie0*qvirial 
    qvir_kc =2.d0*r4pie0*qvirial

   vir = vir + qvir_kc

!   vir = vir - qkcpe*418.4d0

!   write(*,*) 'qvir_kc ',qvir_kc

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

  end subroutine kspace_ewald



! NO EWALD Coulumb Potential and  forces (shifted potential)

  subroutine coul_pot(xx,yy,zz,xxs,yys,zzs,boxi,alp,qrcpe)
    implicit none
    real(8) :: xx(nom),yy(nom),zz(nom)
    real(8) :: xxs(nsite,nom),yys(nsite,nom),zzs(nsite,nom)
    real(8) :: boxi,alp,qrcpe

    integer :: i,j,is,js,k
    real(8) :: qvir_rc
    real(8) :: wij
    real(8) :: xcc,ycc,zcc
    real(8) :: dx,dy,dz
    real(8) :: rsx,rsy,rsz
    real(8) :: rccsq,rssq,drr
    real(8) :: qis,qjs,chgprd
    real(8) :: fc,cfij
    real(8) :: fxsij,fysij,fzsij

    real(8) :: qvirial,lqrcpe,lqvirial
    ! Shared
    real(8), allocatable :: lfxs(:, :), lfys(:, :), lfzs(:, :)
        ! MPI
        INTEGER(4) MPIchunk_size_loc, iii;
        real(8), allocatable :: MPIbuff1(:), MPIbuff2(:);
    real(8), allocatable :: llfxs(:, :), llfys(:, :), llfzs(:, :);
        ! DBG
        real(4) iniTime, deltaTime;
        ! DBG
        ! MPI

    qvirial = 0.d0
    qrcpe =0.d0

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

!$omp parallel DEFAULT(SHARED)&
!$omp& private(i,iii,j,is,js,k)&
!$omp& private(qvir_rc,wij,xcc,ycc,zcc,dx,dy,dz,rsx,rsy,rsz,rccsq,rssq)&
!$omp&
!private(drr,qis,qjs,chgprd,fc,cfij,fxsij,fysij,fzsij,lfxs,lfys,lfzs,lqrcpe,lqvirial)

    lqrcpe = 0.0D0
    lqvirial = 0.0D0

    allocate(lfxs(ns,nm))
    allocate(lfys(ns,nm))
    allocate(lfzs(ns,nm))

    do j=1,nm
       do js=1,ns
          lfxs(js,j) = 0.0D0
          lfys(js,j) = 0.0D0
          lfzs(js,j) = 0.0D0
       end do
    end do

        ! DBG
        iniTime = SECNDS(0.0);
        ! DBG
!     do i=1,nm-1
!         do i = MPIrank + 1, nm - 1, MPIcomm_size
        if (MPIcomm_size .GT. 1) then
                i = 1;
        else
                i = 0;
        endif
        do while (i .LE. nm - 1)
          if (MPIrank .EQ. 0) then
                if (MPIcomm_size .GT. 1) then
                        call MPI_RECV(MPIaddr, 1, MPI_INTEGER, MPI_ANY_SOURCE,MPI_TAG_CommReady, MPI_COMM_WORLD, MPIstat, MPIierr);
                        call MPI_SEND(i, 1, MPI_INTEGER, MPIaddr,MPI_TAG_SendData, MPI_COMM_WORLD, MPIierr);
                endif
                i = i + 1;
          else
                if (MPIcomm_size .GT. 1) then
                        call MPI_SEND(MPIrank, 1, MPI_INTEGER, 0,MPI_TAG_CommReady, MPI_COMM_WORLD, MPIierr);
                        call MPI_RECV(i, 1, MPI_INTEGER, 0, MPI_TAG_SendData,MPI_COMM_WORLD, MPIstat, MPIierr);
                endif
          endif
          if ((i .LE. (nm - 1)) .AND. ((MPIrank .GE. 1) .OR. (MPIcomm_size .EQ.1))) then
                !$omp do schedule(dynamic)
        do j=i+1, nm

           dx = xx(i)-xx(j)
           dy = yy(i)-yy(j)
           dz = zz(i)-zz(j)

           xcc = dx - boxi*nint(dx/boxi)
           ycc = dy - boxi*nint(dy/boxi)
           zcc = dz - boxi*nint(dz/boxi)

           rccsq=xcc**2 + ycc**2 + zcc**2

           if (rccsq.lt.rcutsd) then   ! com cut off starts

              do is=1,ns

                  do js=1,ns

                  qis = chgs(is,i)
                  qjs = chgs(js,j)

                  rsx = xcc + xxs(is,i)-xxs(js,j)
                  rsy = ycc + yys(is,i)-yys(js,j)
                  rsz = zcc + zzs(is,i)-zzs(js,j)

                  rssq = rsx**2 + rsy**2 + rsz**2

!                  if (rssq.le.rcutsd) then       ! site cut off starts

                  drr = sqrt(rssq)

                  lqrcpe = lqrcpe + r4pie0*qis*qjs/drr

                  fc = r4pie0*qis*qjs/drr

                  cfij = fc/rssq

                  fxsij = rsx*cfij
                  fysij = rsy*cfij
                  fzsij = rsz*cfij

                  lfxs(is,i)=lfxs(is,i)+fxsij
                  lfys(is,i)=lfys(is,i)+fysij
                  lfzs(is,i)=lfzs(is,i)+fzsij

                  lfxs(js,j)=lfxs(js,j)-fxsij
                  lfys(js,j)=lfys(js,j)-fysij
                  lfzs(js,j)=lfzs(js,j)-fzsij

                  ! virial = - rab * fab

                  lqvirial = lqvirial -(rsx*fxsij + rsy*fysij + rsz*fzsij)


!                  end if    ! site cut off ends

                 end do

              end do

           end if         ! com cut off ends

        end do
                !$omp end do
           endif
     end do
        if (MPIrank .EQ. 0) then
                do i = 1, MPIcomm_size - 1
                        call MPI_RECV(MPIaddr, 1, MPI_INTEGER, MPI_ANY_SOURCE,MPI_TAG_CommReady, MPI_COMM_WORLD, MPIstat, MPIierr);
                        call MPI_SEND(nm + 1, 1, MPI_INTEGER, MPIaddr,MPI_TAG_SendData, MPI_COMM_WORLD, MPIierr);
                enddo
        endif
        ! DBG
        deltaTime = SECNDS(iniTime);
!       write(dbgUnit, *) 'P#', MPIrank, '. rspace_ewald calc time      ',
!       deltaTime;
        ! DBG
!$omp critical
     qvirial = qvirial + lqvirial
     qrcpe = qrcpe + lqrcpe

     do j=1,nm
        do js=1,ns
           llfxs(js,j) = llfxs(js,j) + lfxs(js,j)
           llfys(js,j) = llfys(js,j) + lfys(js,j)
           llfzs(js,j) = llfzs(js,j) + lfzs(js,j)
        end do
     end do
!$omp end critical

    deallocate(lfxs)
    deallocate(lfys)
    deallocate(lfzs)
!$omp end parallel
! MPI sync
        j = nm * ns * 3 + 2; !fxs, fys, fzs, qvirial, qrcpe
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
        MPIbuff1(k) = qvirial; k = k + 1;
        MPIbuff1(k) = qrcpe; k = k + 1;
        call MPI_ALLREDUCE(MPIbuff1, MPIbuff2, j, MPI_DOUBLE_PRECISION, MPI_SUM,MPI_COMM_WORLD, MPIierr);
        k = 0;
        do i = 1, nm
                do is = 1, ns
                        fxs(is, i) = fxs(is, i) + MPIbuff2(k); k = k + 1;
                        fys(is, i) = fys(is, i) + MPIbuff2(k); k = k + 1;
                        fzs(is, i) = fzs(is, i) + MPIbuff2(k); k = k + 1;
                enddo
        enddo
        qvirial = MPIbuff2(k); k = k + 1;
        qrcpe = MPIbuff2(k); k = k + 1;
        deallocate(MPIbuff1);
        deallocate(MPIbuff2);
! MPI sync
    deallocate(llfxs);
    deallocate(llfys);
    deallocate(llfzs);

!     qvir_rc = -qrcpe

     qrcpe = qrcpe/418.4d0       ! real space energy for coulomb in internal units

!     qvir_rc = -qrcpe*418.4d0 !qvirial
     qvir_rc = qvirial

     vir = vir + qvir_rc

!     vir_test = vir_test - qrcpe*418.4d0

!     write(*,*) 'ewald_rspace: ',qrcpe,qvir_rc/418.4d0
!     !vir,vir_test,qvir_rc,-qrcpe*418.4d0

  end subroutine coul_pot

end module ewald

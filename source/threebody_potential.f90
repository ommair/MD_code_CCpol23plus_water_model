! Looks good, but MPI processes are unbalanced.
! Use of Master thread, discard existing OpenMP
module  threebody_potential

  use, intrinsic :: iso_fortran_env
  use common_arrays
  use utility

  implicit none
  private

  public :: threebody_term
!  public :: FS2_term, FS3_term

contains  

  subroutine threebody_term(xxx,yyy,zzz,xxxs,yyys,zzzs,boxii,qfs2,qfs3)

    implicit none
    real(8) :: xxx(nom),yyy(nom),zzz(nom)
    real(8) :: xxxs(nsite,nom),yyys(nsite,nom),zzzs(nsite,nom)
    real(8) :: boxii
    real(8) :: qfs2,qfs3

    integer :: is,js,i,j,k,isteps
!    real(8) :: xx(nom),yy(nom),zz(nom)
!    real(8) :: xxs(nsite,nom),yys(nsite,nom),zzs(nsite,nom)
    real(8), allocatable :: xx(:), yy(:), zz(:)
    real(8), allocatable :: xxs(:, :), yys(:, :), zzs(:, :)
    real(8) :: boxi

    real(8) :: tol3,eprec3,alph3
	! DBG
	real(4) iniTime, deltaTime;
	real(4) iniTimeTot, deltaTimeTot;
	! DBG

	! DBG
	iniTimeTot = SECNDS(0.0);
	! DBG
	allocate(xx(nom));
	allocate(yy(nom));
	allocate(zz(nom));
	allocate(xxs(nsite,nom));
	allocate(yys(nsite,nom));
	allocate(zzs(nsite,nom));

!    eprec3 = min(abs(ewp),0.5d0)
!    tol3=sqrt(abs(log(eprec3*rcut3)))
!    alpha3=sqrt(abs(log(eprec3*rcut3*tol3)))/rcut3  ! calculating alpha for 3B interaction for a givein rcut3

    do i=1,nm

       xx(i) = xxx(i)/bohr2a
       yy(i) = yyy(i)/bohr2a
       zz(i) = zzz(i)/bohr2a

       do is=1,ns
          xxs(is,i) = (xxxs(is,i))/bohr2a
          yys(is,i) = (yyys(is,i))/bohr2a
          zzs(is,i) = (zzzs(is,i))/bohr2a
       end do

    end do

    boxi = boxii/bohr2a
    alpha3 = alpha
    alph3 = alpha3*bohr2a

	! DBG
	iniTime = SECNDS(0.0);
	! DBG
    call FS2_term(xx,yy,zz,xxs,yys,zzs,boxi,alph3,qfs2)
	! DBG
	deltaTime = SECNDS(iniTime);
!	write(dbgUnit, *) 'P#', MPIrank, '. FS2 time	', deltaTime;
	iniTime = SECNDS(0.0);
	! DBG
    call FS3_term(xx,yy,zz,xxs,yys,zzs,boxi,alph3,qfs3)
	! DBG
	deltaTime = SECNDS(iniTime);
!	write(dbgUnit, *) 'P#', MPIrank, '. FS3 time	', deltaTime;
	! DBG

	deallocate(xx);
	deallocate(yy);
	deallocate(zz);
	deallocate(xxs);
	deallocate(yys);
	deallocate(zzs);
	! DBG
	deltaTimeTot = SECNDS(iniTimeTot);
!	write(dbgUnit, *) 'P#', MPIrank, '. threebody_term time	', deltaTimeTot;
	call MPI_BARRIER(MPI_COMM_WORLD, MPIierr);
!	STOP;
	! DBG
	return;
  end subroutine threebody_term

  subroutine FS2_term(xx,yy,zz,xxs,yys,zzs,boxi,alp3,qfs2)
    implicit none
    real(8) :: xx(nom),yy(nom),zz(nom)
    real(8) :: xxs(nsite,nom),yys(nsite,nom),zzs(nsite,nom)
    real(8) :: boxi
    real(8) :: alp3,qfs2,lqvirial,qvirial,fs2_vir

!    real(8) :: xxx(nom),yyy(nom),zzz(nom)
!    real(8) :: xxxs(nsite,nom),yyys(nsite,nom),zzzs(nsite,nom)
!    real(8) :: charge3b(nsite)
    real(8), allocatable :: xxx(:), yyy(:), zzz(:);
    real(8), allocatable :: xxxs(:, :), yyys(:, :), zzzs(:, :);
    real(8), allocatable :: charge3b(:);
    real(8) :: boxii  ! charges on sites
    real(8) :: alph3

    integer :: k
    real(8) qqfs2,qqfs20,efs2
    integer :: is,js,ks,i,j,isteps,l
    real(8) :: dsij2,dski2,dskj2,dsij,dski,dskj
    real(8) :: dsij1i,rc2_3b

    real(8) :: dij2,dki2,dkj2
    real(8) :: dij,dij1i,dij3i
    real(8) :: dkj,dkj1i,dkj3i
    real(8) :: dki,dki1i,dki3i

!    real(8) :: qis,qjs,qks,boxi  ! charges on sites
    real(8) :: qis,qjs,qks  ! charges on sites
    real(8) :: eqij,eqijsum 

    real(8) :: dxij,dyij,dzij,dxki,dyki,dzki,dxkj,dykj,dzkj
    real(8) :: xij,yij,zij,xki,yki,zki,xkj,ykj,zkj

    real(8) :: xsij,ysij,zsij
    real(8) :: xqi,yqi,zqi,xqj,yqj,zqj          ! position of additional charges w.r.t com
    real(8) :: dxkqi,dykqi,dzkqi,dxkqj,dykqj,dzkqj     ! sep of site k and xqi ...
    real(8) :: dxkci,dykci,dzkci,dxkcj,dykcj,dzkcj     ! sep of site k and com i ...

    real(8) :: dkqi2,dkqj2
    real(8) :: dkqi,dkqj,dkqi1i,dkqj1i,dkqi2i,dkqi3i,dkqj2i,dkqj3i

    real(8) :: dkci2,dkcj2
    real(8) :: dkci,dkcj,dkci1i,dkcj1i,dkci2i,dkci3i,dkcj2i,dkcj3i
    real(8) :: dotkqi_ij,dotkqj_ij             
    real(8) :: term0,termxy1,termxy2    
    real(8) :: gradxqij,gradyqij,gradzqij
    real(8) :: torqxija,torqyija,torqzija
    real(8) :: torqxijb,torqyijb,torqzijb  

    real(8) :: fackci,fackcj,fackqi,fackqj

!    real(8) :: qqfs2,qfs2,alp3

    real(8) :: fxi1,fyi1,fzi1,fxi2,fyi2,fzi2,fxi3,fyi3,fzi3,fxi4,fyi4,fzi4
    real(8) :: fxj1,fyj1,fzj1,fxj2,fyj2,fzj2,fxj3,fyj3,fzj3,fxj4,fyj4,fzj4
    real(8) :: fxk,fyk,fzk

    real(8) :: vpxi,vpyi,vpzi 
    real(8) :: vpxj,vpyj,vpzj

    real(8) :: qfxt, qfyt, qfzt

    real(8) :: erf_ki_kj, derf_ki,derf_kj

    !   long range damping as used by Robert starts
    real(8), parameter :: epsprim = 1.d-10
    real(8), parameter :: gamdmp = 1.62823d0*bohr2a
    real(8), parameter :: delta3 = 10.d0 !log(1.d0/epsprim)/log(2.d0)
    real(8) :: r0dmp,flrki,flrkj,dflrki,dflrkj,dflrkiv,dflrkjv
    real(8) :: pom,pomkj,pomki
    real(8) :: termx1,termy1,termz1,termx2,termy2,termz2
    !   long range damping as used by Robert ends

! Shared
    real(8), allocatable :: lfxs(:,:),lfys(:,:),lfzs(:,:)
    real(8), allocatable :: lfxfs2(:),lfyfs2(:),lfzfs2(:)

    real(8), parameter :: eta3 =0.48056977605478400d0  !0.693231401521011**2.d0  ! eta OO  HH OH
    real(8), parameter :: beta3 = 0.85836606579445000d0 !/bohr2a !2.d0*0.429183032897225  ! beta OO HH OH
    real(8), parameter :: g3 = 0.400000000000000 !*bohr2a 
	! MPI
	INTEGER(4) MPIchunk_size_loc, iii;
	real(8), allocatable :: MPIbuff1(:), MPIbuff2(:);
    real(8), allocatable :: llfxs(:, :), llfys(:, :), llfzs(:, :);
	! DBG
!	real(4) iniTime, deltaTime;
	! DBG
	! MPI

	allocate(charge3b(nsite));

    rc2_3b = rcutsd3*bohr2a**(-2.d0)

    do i=1,ns
       charge3b(i) = 0.d0
    end do

    charge3b(1) = 0.258504607653328389d0
    charge3b(2) = 0.564050362569283648d0
    charge3b(3) = 0.564050362569283648d0
    charge3b(4) = -0.693302666395947953d0
    charge3b(5) = -0.693302666395947953d0

    do i=1,nm
       fxfs2(i) = 0.d0
       fyfs2(i) = 0.d0
       fzfs2(i) = 0.d0 

       txfs2(i) = 0.d0
       tyfs2(i) = 0.d0
       tzfs2(i) = 0.d0 
    end do
 
    qqfs2 = 0.d0
    qqfs20 = 0.d0

    qvirial = 0.d0
    lqvirial = 0.d0 

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

!   Molecule k is the spectator interacting electrostatically
!   with the exchange quadrupole of molecules i and j.
  	! DBG
	! if (MPIrank .EQ. 0) then
		! write(202, '(a)') '#';
		! write(202, '(a)') 'fxs';
		! write(202, '(1p3E16.6E3)') fxs;
		! write(202, '(a)') 'fys';
		! write(202, '(1p3E16.6E3)') fys;
		! write(202, '(a)') 'fzs';
		! write(202, '(1p3E16.6E3)') fzs;
		! write(202, '(a)') 'fxfs2';
		! write(202, '(1p3E16.6E3)') fxfs2;
		! write(202, '(a)') 'fyfs2';
		! write(202, '(1p3E16.6E3)') fyfs2;
		! write(202, '(a)') 'fzfs2';
		! write(202, '(1p3E16.6E3)') fzfs2;
		! write(202, '(a)') 'qfs2';
		! write(202, '(1p3E16.6E3)') qfs2;
	! endif
	! DBG

!$omp parallel DEFAULT(SHARED)&
!$omp& private(k,id,is,js,ks,i,iii,j,isteps,l,dsij2,dski2,dskj2,dsij,dski,dskj,dsij1i,dij2,dki2) &
!$omp& private(dkj2,dij,dij1i,dij3i,dkj,dkj1i,dkj3i,dki,dki1i,dki3i)&
!$omp& private(qis,qjs,qks,eqij,eqijsum,dxij,dyij,dzij,dxki,dyki,dzki)&
!$omp& private(dxkj,dykj,dzkj,xij,yij,zij,xki,yki,zki,xkj,ykj,zkj,xsij)&
!$omp& private(ysij,zsij,xqi,yqi,zqi,xqj,yqj,zqj,dxkqi,dykqi,dzkqi,dxkqj)&
!$omp& private(dykqj,dzkqj,dxkci,dykci,dzkci,dxkcj,dykcj,dzkcj,dkqi2,dkqj2)&
!$omp& private(dkqi,dkqj,dkqi1i,dkqj1i,dkqi2i,dkqi3i,dkqj2i,dkqj3i,dkci2,dkcj2)&
!$omp& private(dkci,dkcj,dkci1i,dkcj1i,dkci2i,dkci3i,dkcj2i,dkcj3i,dotkqi_ij)&
!$omp& private(dotkqj_ij,term0,gradxqij,gradyqij,gradzqij,torqxija,torqyija)&
!$omp& private(torqzija,torqxijb,torqyijb,torqzijb,fackci,fackcj,fackqi)&
!$omp& private(fackqj,fxi1,fyi1,fzi1,fxi2,fyi2,fzi2,fxi3,fyi3,fzi3,fxi4)&
!$omp& private(fyi4,fzi4,fxj1,fyj1,fzj1,fxj2,fyj2,fzj2,fxj3,fyj3,fzj3,fxj4)&
!$omp& private(fyj4,fzj4,fxk,fyk,fzk,vpxi,vpyi,vpzi,vpxj,vpyj,vpzj,qfxt,qfyt)&
!$omp& private(qfzt,erf_ki_kj,derf_ki,derf_kj,lfxs,lfys,lfzs,lfxfs2,lfyfs2,lfzfs2)&
!$omp& private(r0dmp,flrki,flrkj,dflrki,dflrkj,pom,pomkj,pomki)&
!$omp& private(termx1,termy1,termz1,termx2,termy2,termz2)&
!$omp& private(xxx,yyy,zzz,xxxs,yyys,zzzs)

    allocate(lfxs(ns,nm))
    allocate(lfys(ns,nm))
    allocate(lfzs(ns,nm))
    allocate(lfxfs2(nm)) 
    allocate(lfyfs2(nm)) 
    allocate(lfzfs2(nm)) 
	allocate(xxx(nom));
	allocate(yyy(nom));
	allocate(zzz(nom));
	allocate(xxxs(nsite, nom));
	allocate(yyys(nsite, nom));
	allocate(zzzs(nsite, nom));

    do i=1,nm
       lfxfs2(i) = 0.d0
       lfyfs2(i) = 0.d0
       lfzfs2(i) = 0.d0 
       do j=1,ns
          lfxs(j,i) = 0.0D0
          lfys(j,i) = 0.0D0
          lfzs(j,i) = 0.0D0
       enddo
    enddo

!$omp do schedule(dynamic) reduction(+:qqfs2,qqfs20)
!    do k=1,nm
	 do k = MPIrank + 1, nm, MPIcomm_size
       do i=1,nm-1
           if (k .NE. i) then
          do j=i+1,nm
             if (j .NE. k) then

!              write(*,*) 'k, i, j: ',k,i,j 

              dxij = xx(i)-xx(j)   ! com sep of  i - j
              dyij = yy(i)-yy(j)
              dzij = zz(i)-zz(j)

              dxki = xx(k)-xx(i)   ! com sep of  k - i
              dyki = yy(k)-yy(i)
              dzki = zz(k)-zz(i)

              dxkj = xx(k)-xx(j)   ! com sep of  k - j
              dykj = yy(k)-yy(j)
              dzkj = zz(k)-zz(j)
 
              xij = dxij - boxi*nint(dxij/boxi)  ! minimum image convenion
              yij = dyij - boxi*nint(dyij/boxi)
              zij = dzij - boxi*nint(dzij/boxi)

              xki = dxki - boxi*nint(dxki/boxi)
              yki = dyki - boxi*nint(dyki/boxi)
              zki = dzki - boxi*nint(dzki/boxi)

              xkj = dxkj - boxi*nint(dxkj/boxi)
              ykj = dykj - boxi*nint(dykj/boxi)
              zkj = dzkj - boxi*nint(dzkj/boxi)

              dij2 = xij**2 + yij**2 + zij**2   ! rij2
              dki2 = xki**2 + yki**2 + zki**2   ! rki2
              dkj2 = xkj**2 + ykj**2 + zkj**2   ! rkj2

              if ( (dij2 .lt. rc2_3b) .and. &
                   (dki2 .lt. rc2_3b) .and. &
                   (dkj2 .lt. rc2_3b) ) then  ! com-com cutoff  on i-j, i-k, j-k starts

!                    r0dmp = 8.81785d0/bohr2a
                    r0dmp = 8.36785d0/bohr2a

                    ! long range damping ****  dkj2 dki2
                    pom = exp(gamdmp*(sqrt(dki2)-r0dmp))  ! ks and com i+
                    pomki = 1.d0/(1.d0+pom)
                    flrki = pomki**delta3
                    dflrki = -gamdmp*delta3*pom*pomki*flrki
                    dflrkiv = -gamdmp*delta3*pom*pomki*sqrt(dki2)

                    pom = exp(gamdmp*(sqrt(dkj2)-r0dmp))  ! ks and com j+
                    pomkj = 1.d0/(1.d0+pom)
                    flrkj = pomkj**delta3
                    dflrkj = -gamdmp*delta3*pom*pomkj*flrkj 
                    dflrkjv = -gamdmp*delta3*pom*pomkj*sqrt(dkj2)

                    !write(*,*) flrki,flrkj,dflrki,dflrkj,gamdmp,r0dmp,pom,delta3,pomki
 
                    do ks=1, 5  ! loop ks

                       qks= charge3b(ks) !chgs(ks,k)

!                       if (qks .ne. 0.d0 ) then

                      ! Calculate the "exchange charge", i.e., eta*sum(exp(-beta*r_ab))
                    
                          do is=1,3    ! loop is         only for OO, OH, HH

                             do js=1,3 ! loop js
                           
                                xsij = xij + (xxs(is,i) - xxs(js,j))
                                ysij = yij + (yys(is,i) - yys(js,j))
                                zsij = zij + (zzs(is,i) - zzs(js,j))     

                                dsij2 = xsij**2 + ysij**2 + zsij**2 

                                dsij = sqrt(dsij2)
                                dsij1i = 1.d0/dsij 
                               
                                eqij = eta3*exp(-beta3*dsij)

                                ! Calculate the positions of exchange charges 

                                dij  = sqrt(dij2) ! separation between i-j
                                dij1i = 1.d0/dij  ! inverse separation between i-j
                                dij3i = dij1i**3

                                dki  = sqrt(dki2) ! separation between k-i
                                dki1i = 1.d0/dki  ! inverse separation between k-i
                                dki3i = dki1i**3 

                                dkj  = sqrt(dkj2) ! separation between k-j
                                dkj1i = 1.d0/dkj  ! inverse separation between k-j
                                dkj3i = dkj1i**3

                                xqi = xx(i) + g3*xij*dij1i  ! comp of pos of i- charges
                                yqi = yy(i) + g3*yij*dij1i
                                zqi = zz(i) + g3*zij*dij1i

                                xqj = xx(j) - g3*xij*dij1i  ! comp of pos of i- charges
                                yqj = yy(j) - g3*yij*dij1i
                                zqj = zz(j) - g3*zij*dij1i
                                
                               ! separation of site ks and additional charges i- and j-
  
                                dxkqi = (xx(k)+xxs(ks,k)) - xqi       ! sep b/w charge at ks and i-
                                dykqi = (yy(k)+yys(ks,k)) - yqi
                                dzkqi = (zz(k)+zzs(ks,k)) - zqi

                                dkqi2 = dxkqi**2 + dykqi**2 + dzkqi**2
                                dkqi = sqrt(dkqi2)
                                dkqi1i = 1.d0/dkqi
                                dkqi2i = dkqi1i**2
                                dkqi3i = dkqi1i**3

                                dxkqj = (xx(k)+xxs(ks,k)) - xqj       ! sep b/w charge at ks and j-
                                dykqj = (yy(k)+yys(ks,k)) - yqj
                                dzkqj = (zz(k)+zzs(ks,k)) - zqj

                                dkqj2 = dxkqj**2 + dykqj**2 + dzkqj**2
                                dkqj = sqrt(dkqj2)
                                dkqj1i = 1.d0/dkqj
                                dkqj2i = dkqj1i**2
                                dkqj3i = dkqj1i**3

                                ! separation of site ks and additional charges i+ and j+

                                dxkci = (xx(k)+xxs(ks,k)) - xx(i)     ! sep ks and com i+
                                dykci = (yy(k)+yys(ks,k)) - yy(i)
                                dzkci = (zz(k)+zzs(ks,k)) - zz(i)

                                dkci2 = dxkci**2 + dykci**2 + dzkci**2
                                dkci = sqrt(dkci2)
                                dkci1i = 1.d0/dkci
                                dkci2i = dkci1i**2
                                dkci3i = dkci1i**3

                                dxkcj = (xx(k)+xxs(ks,k)) - xx(j)     ! sep ks and com j+
                                dykcj = (yy(k)+yys(ks,k)) - yy(j)
                                dzkcj = (zz(k)+zzs(ks,k)) - zz(j)

                                dkcj2 = dxkcj**2 + dykcj**2 + dzkcj**2
                                dkcj = sqrt(dkcj2)
                                dkcj1i = 1.d0/dkcj
                                dkcj2i = dkcj1i**2
                                dkcj3i = dkcj1i**3

                               ! *** FS2 energy contribution *** 
                         
                               term0 = (dkci1i + dkcj1i - dkqi1i - dkqj1i) 
                               qqfs20 = qqfs20 + qks*eqij*term0
                               qqfs2 = qqfs2 + qks*eqij*term0*flrki*flrkj

                               !!! virial part start
                               termxy1 =(dxkqi*xij+dykqi*yij+dzkqi*zij)*g3*dij1i*dkqi3i 
                               termxy2 =(dxkqj*xij+dykqj*yij+dzkqj*zij)*g3*dij1i*dkqj3i

                               lqvirial = lqvirial + ((-beta3*dsij - 1.d0) + dflrkiv + dflrkjv)*   &
                                                      qks*eqij*term0*flrki*flrkj                   &
                                                   + qks*eqij*(termxy1-termxy2)*flrki*flrkj   

!                               lqvirial = lqvirial + ((-beta3*dsij - 1.d0) )*   &
!                                                      qks*eqij*term0*flrki*flrkj                   &
!                                                   + qks*eqij*(termxy1-termxy2)*flrki*flrkj
 
                               !!! virial part ends
 
                               ! *** FS2 Forces contribution ***

                               ! scalar product vec_ij * vec_kiq and vec_ij * vec_kjq 
                               dotkqi_ij = xij*dxkqi + yij*dykqi + zij*dzkqi
                               dotkqj_ij = xij*dxkqj + yij*dykqj + zij*dzkqj                 

                              ! Contribution to derivatives over coordinates of spectator molecule (k)

                               fxk = qks*eqij*(-dxkci*dkci3i-dxkcj*dkcj3i+dxkqi*dkqi3i+dxkqj*dkqj3i)
                               fyk = qks*eqij*(-dykci*dkci3i-dykcj*dkcj3i+dykqi*dkqi3i+dykqj*dkqj3i)
                               fzk = qks*eqij*(-dzkci*dkci3i-dzkcj*dkcj3i+dzkqi*dkqi3i+dzkqj*dkqj3i)

                               lfxfs2(k) = lfxfs2(k) - qks*eqij*term0*dflrki*xki*flrkj/sqrt(dki2) &
                                                     - qks*eqij*term0*dflrkj*xkj*flrki/sqrt(dkj2)

                               lfyfs2(k) = lfyfs2(k) - qks*eqij*term0*dflrki*yki*flrkj/sqrt(dki2) &
                                                     - qks*eqij*term0*dflrkj*ykj*flrki/sqrt(dkj2)

                               lfzfs2(k) = lfzfs2(k) - qks*eqij*term0*dflrki*zki*flrkj/sqrt(dki2) &
                                                     - qks*eqij*term0*dflrkj*zkj*flrki/sqrt(dkj2) 

                              ! Contribution to derivatives over coordinates of molecule (i)

                               fxi1 = qks*eqij*(-beta3)*term0*xsij*dsij1i   
                               fyi1 = qks*eqij*(-beta3)*term0*ysij*dsij1i 
                               fzi1 = qks*eqij*(-beta3)*term0*zsij*dsij1i

                               fxi2 = qks*eqij*(dxkci*dkci3i)
                               fyi2 = qks*eqij*(dykci*dkci3i)
                               fzi2 = qks*eqij*(dzkci*dkci3i)

                               fxi3 = qks*eqij*dkqi3i*(-(1.d0+g3*dij1i)*dxkqi + g3*dotkqi_ij*xij*dij3i)
                               fyi3 = qks*eqij*dkqi3i*(-(1.d0+g3*dij1i)*dykqi + g3*dotkqi_ij*yij*dij3i)
                               fzi3 = qks*eqij*dkqi3i*(-(1.d0+g3*dij1i)*dzkqi + g3*dotkqi_ij*zij*dij3i)

                               fxi4 = qks*eqij*dkqj3i*(g3*dxkqj*dij1i - g3*dotkqj_ij*xij*dij3i)
                               fyi4 = qks*eqij*dkqj3i*(g3*dykqj*dij1i - g3*dotkqj_ij*yij*dij3i)
                               fzi4 = qks*eqij*dkqj3i*(g3*dzkqj*dij1i - g3*dotkqj_ij*zij*dij3i)

                              ! Contribution to derivatives over coordinates of molecule (j)

                               fxj1 = qks*eqij*beta3*term0*xsij*dsij1i
                               fyj1 = qks*eqij*beta3*term0*ysij*dsij1i
                               fzj1 = qks*eqij*beta3*term0*zsij*dsij1i

                               fxj2 = qks*eqij*(dxkcj*dkcj3i)
                               fyj2 = qks*eqij*(dykcj*dkcj3i)
                               fzj2 = qks*eqij*(dzkcj*dkcj3i)

                               fxj3 = qks*eqij*dkqi3i*(g3*dxkqi*dij1i - g3*dotkqi_ij*xij*dij3i)
                               fyj3 = qks*eqij*dkqi3i*(g3*dykqi*dij1i - g3*dotkqi_ij*yij*dij3i)
                               fzj3 = qks*eqij*dkqi3i*(g3*dzkqi*dij1i - g3*dotkqi_ij*zij*dij3i)

                               fxj4 = qks*eqij*dkqj3i*(-(1.d0+g3*dij1i)*dxkqj + g3*dotkqj_ij*xij*dij3i)
                               fyj4 = qks*eqij*dkqj3i*(-(1.d0+g3*dij1i)*dykqj + g3*dotkqj_ij*yij*dij3i)
                               fzj4 = qks*eqij*dkqj3i*(-(1.d0+g3*dij1i)*dzkqj + g3*dotkqj_ij*zij*dij3i)

                               lfxs(ks,k)= lfxs(ks,k) - fxk*flrki*flrkj  
                               lfys(ks,k)= lfys(ks,k) - fyk*flrki*flrkj 
                               lfzs(ks,k)= lfzs(ks,k) - fzk*flrki*flrkj                     
                               lfxs(is,i)= lfxs(is,i) - fxi1*flrki*flrkj 
                               lfys(is,i)= lfys(is,i) - fyi1*flrki*flrkj 
                               lfzs(is,i)= lfzs(is,i) - fzi1*flrki*flrkj 
                               lfxs(js,j)= lfxs(js,j) - fxj1*flrki*flrkj 
                               lfys(js,j)= lfys(js,j) - fyj1*flrki*flrkj 
                               lfzs(js,j)= lfzs(js,j) - fzj1*flrki*flrkj 

                               lfxfs2(i) = lfxfs2(i) - (fxi2+fxi3+fxi4)*flrki*flrkj &
                                                     + qks*eqij*term0*dflrki*xki*flrkj/sqrt(dki2)  

                               lfyfs2(i) = lfyfs2(i) - (fyi2+fyi3+fyi4)*flrki*flrkj &
                                                     + qks*eqij*term0*dflrki*yki*flrkj/sqrt(dki2)

                               lfzfs2(i) = lfzfs2(i) - (fzi2+fzi3+fzi4)*flrki*flrkj &
                                                     + qks*eqij*term0*dflrki*zki*flrkj/sqrt(dki2)

                               lfxfs2(j) = lfxfs2(j) - (fxj2+fxj3+fxj4)*flrki*flrkj &
                                                     + qks*eqij*term0*flrki*dflrkj*xkj/sqrt(dkj2)  

                               lfyfs2(j) = lfyfs2(j) - (fyj2+fyj3+fyj4)*flrki*flrkj &
                                                     + qks*eqij*term0*flrki*dflrkj*ykj/sqrt(dkj2)

                               lfzfs2(j) = lfzfs2(j) - (fzj2+fzj3+fzj4)*flrki*flrkj &
                                                     + qks*eqij*term0*flrki*dflrkj*zkj/sqrt(dkj2)

                           end do ! loop js
                          end do    ! loop is

!                       end if 

                    end do  ! loop ks

          !          lfxfs2(k) = lfxfs2(k) - qqfs20*dflrki*xki*flrkj/sqrt(dki2) &
          !                                - qqfs20*dflrkj*xkj*flrki/sqrt(dkj2)

          !          lfyfs2(k) = lfyfs2(k) - qqfs20*dflrki*yki*flrkj/sqrt(dki2) &
          !                                - qqfs20*dflrkj*ykj*flrki/sqrt(dkj2)
  
          !          lfzfs2(k) = lfzfs2(k) - qqfs20*dflrki*zki*flrkj/sqrt(dki2) &
          !                                - qqfs20*dflrkj*zkj*flrki/sqrt(dkj2)

          !          lfxfs2(i) = lfxfs2(i) + qqfs20*dflrki*xki*flrkj/sqrt(dki2)
          !          lfyfs2(i) = lfyfs2(i) + qqfs20*dflrki*yki*flrkj/sqrt(dki2)
          !          lfzfs2(i) = lfzfs2(i) + qqfs20*dflrki*zki*flrkj/sqrt(dki2)
          !          lfxfs2(j) = lfxfs2(j) + qqfs20*flrki*dflrkj*xkj/sqrt(dkj2)
          !          lfyfs2(j) = lfyfs2(j) + qqfs20*flrki*dflrkj*ykj/sqrt(dkj2)
          !          lfzfs2(j) = lfzfs2(j) + qqfs20*flrki*dflrkj*zkj/sqrt(dkj2)                                   

              end if   ! com-com cutoff on i-j, i-k, j-k ends

!          12 continue       
           endif !if (j .NE. k) then
          end do  ! k loop

!       11 continue
        endif !if (i .NE. k) then
       end do ! j loop
    end do    ! i loop
!$omp end do

!$omp critical
    do i=1,nm
       fxfs2(i) = fxfs2(i) + lfxfs2(i)*h2kcal*418.4d0/bohr2a
       fyfs2(i) = fyfs2(i) + lfyfs2(i)*h2kcal*418.4d0/bohr2a
       fzfs2(i) = fzfs2(i) + lfzfs2(i)*h2kcal*418.4d0/bohr2a
       do is=1,ns
          llfxs(is,i) = llfxs(is,i) + lfxs(is,i)*h2kcal*418.4d0/bohr2a
          llfys(is,i) = llfys(is,i) + lfys(is,i)*h2kcal*418.4d0/bohr2a
          llfzs(is,i) = llfzs(is,i) + lfzs(is,i)*h2kcal*418.4d0/bohr2a
       enddo
    enddo
!$omp end critical

    deallocate(lfxs)
    deallocate(lfys)
    deallocate(lfzs)
    deallocate(lfxfs2)
    deallocate(lfyfs2)
    deallocate(lfzfs2)
	deallocate(xxx);
	deallocate(yyy);
	deallocate(zzz);
	deallocate(xxxs);
	deallocate(yyys);
	deallocate(zzzs);
!$omp end parallel

!    qfs2 = qqfs2*r4pie0/418.4d0

     qvirial = qvirial + lqvirial
     qfs2 = qqfs2*h2kcal
!     write(*,*) qfs2
	deallocate(charge3b);

! MPI sync
	j = nm * ns * 3 + nm * 3 + 2; !fxs, fys, fzs, fxfs3, fyfs3, fzfs3,qfs2,qvirial
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
		MPIbuff1(k) = fxfs2(i); k = k + 1;
		MPIbuff1(k) = fyfs2(i); k = k + 1;
		MPIbuff1(k) = fzfs2(i); k = k + 1;
	enddo
	MPIbuff1(k) = qfs2; k = k + 1;
        MPIbuff1(k) = qvirial; k = k + 1;

	call MPI_ALLREDUCE(MPIbuff1, MPIbuff2, j, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, MPIierr);
	k = 0;
	do i = 1, nm
		do is = 1, ns
			fxs(is, i) = fxs(is, i) + MPIbuff2(k); k = k + 1;
			fys(is, i) = fys(is, i) + MPIbuff2(k); k = k + 1;
			fzs(is, i) = fzs(is, i) + MPIbuff2(k); k = k + 1;
		enddo
		fxfs2(i) = MPIbuff2(k); k = k + 1;
		fyfs2(i) = MPIbuff2(k); k = k + 1;
		fzfs2(i) = MPIbuff2(k); k = k + 1;
	enddo
	qfs2 = MPIbuff2(k); k = k + 1;
        qvirial = MPIbuff2(k); k = k + 1;
	deallocate(MPIbuff1);
	deallocate(MPIbuff2);
! MPI sync
    deallocate(llfxs);
    deallocate(llfys);
    deallocate(llfzs);
	! DBG
	! if (MPIrank .EQ. 0) then
		! write(302, '(a)') '#';
		! write(302, '(a)') 'fxs';
		! write(302, '(1p3E16.6E3)') fxs;
		! write(302, '(a)') 'fys';
		! write(302, '(1p3E16.6E3)') fys;
		! write(302, '(a)') 'fzs';
		! write(302, '(1p3E16.6E3)') fzs;
		! write(302, '(a)') 'fxfs2';
		! write(302, '(1p3E16.6E3)') fxfs2;
		! write(302, '(a)') 'fyfs2';
		! write(302, '(1p3E16.6E3)') fyfs2;
		! write(302, '(a)') 'fzfs2';
		! write(302, '(1p3E16.6E3)') fzfs2;
		! write(302, '(a)') 'qfs2';
		! write(302, '(1p3E16.6E3)') qfs2;
	! endif
	! DBG

    fs2_vir = qvirial*h2kcal*418.4d0

    vir = vir + fs2_vir

    return;
end subroutine FS2_term

subroutine FS3_term(xx,yy,zz,xxs,yys,zzs,boxi,alp3,qfs3)
    implicit none
    real(8) :: xx(nom),yy(nom),zz(nom)
    real(8) :: xxs(nsite,nom),yys(nsite,nom),zzs(nsite,nom)
    real(8) :: boxi
    real(8) :: qfs3,alp3,lqvirial,qvirial

!    real(8) :: xxx(nom),yyy(nom),zzz(nom)
!    real(8) :: xxxs(nsite,nom),yyys(nsite,nom),zzzs(nsite,nom)
    real(8), allocatable :: xxx(:), yyy(:), zzz(:)
    real(8), allocatable :: xxxs(:, :), yyys(:, :), zzzs(:, :)
    real(8) :: alph3
    real(8) :: boxii
    real(8) :: pomij,ppomij,fdmpij,dfdmpij,r_av,rc3
    real(8) :: pomik,ppomik,fdmpik,dfdmpik
    real(8) :: pomjk,ppomjk,fdmpjk,dfdmpjk
    real(8) :: dfdmpijv,dfdmpikv,dfdmpjkv

    integer :: is,js,ks,i,j,k,ii,idx
    real(8) :: qis,qjs,rc2_3b
    
    real(8) :: dxij,dyij,dzij,dxik,dyik,dzik,dxjk,dyjk,dzjk
    real(8) :: xij,yij,zij,xik,yik,zik,xjk,yjk,zjk

    real(8) :: dij2,dik2,djk2

    real(8) :: xsij,ysij,zsij,dsij2,dsij,dsij1i
    real(8) :: xsik,ysik,zsik,dsik2,dsik,dsik1i
    real(8) :: xsjk,ysjk,zsjk,dsjk2,dsjk,dsjk1i

    real(8) :: dampsij,gsij,r0sij,bsij
    real(8) :: dampsik,gsik,r0sik,bsik
    real(8) :: dampsjk,gsjk,r0sjk,bsjk

    real(8) :: eksij,eksik,eksjk,eks
    
    real(8) :: qqfs3, lqqfs3,fs3_vir
! Shared
    real(8), allocatable :: lfxs(:,:),lfys(:,:),lfzs(:,:)
    real(8), allocatable :: lfxfs3(:),lfyfs3(:),lfzfs3(:)

    integer, parameter :: ns3=17, lmax=1
    integer, parameter :: nzero = ns3*ns3*ns3*(lmax+1)*(lmax+1)*(lmax+1)
!    integer :: ind_lin (0:lmax,0:lmax,0:lmax,ns3,ns3,ns3)
!    integer :: ind_lin_2 (nzero)
    integer, allocatable :: ind_lin(:,:,:,:,:,:)
!    integer, pointer :: ind_lin_2(:) => NULL()
    integer, allocatable :: ind_lin_2(:)
!    equivalence (ind_lin,ind_lin_2)
    logical :: init
    integer :: iat,jat,kat,ip,jp,kp,new,nlin,itab,iind,nl
!    real(8) :: aj(5000),term1,term2,term3
	real(8) :: term1,term2,term3
    real(8), allocatable :: aj(:) ! (5000)

    real(8) :: der_fij,der_fik,der_fjk,der_fijv,der_fikv,der_fjkv
    real(8) :: fsij,fsik,fsjk

    real(8) :: erf_ij_ik_jk, der_erfij,der_erfjk,der_erfik
    real(8), parameter :: alp =1.62823d0*bohr2a !10.0d0*bohr2a !9.21d0*bohr2a
    real(8), parameter :: epsprim = 1.d-10
!    real(8), parameter :: gamdmp = 10.d0*bohr2a
    real(8), parameter :: delta3 = 10.d0 !log(1.d0/epsprim)/log(2.d0)
	! MPI
	INTEGER(4) MPIchunk_size_loc, iii;
	real(8), allocatable :: MPIbuff1(:), MPIbuff2(:);
    real(8), allocatable :: llfxs(:,:),llfys(:,:),llfzs(:,:)
	! DBG
	real(4) iniTime, deltaTime;
	! DBG
	! MPI

    do i=1,nm
       fxfs3(i) = 0.d0
       fyfs3(i) = 0.d0
       fzfs3(i) = 0.d0

       txfs3(i) = 0.d0
       tyfs3(i) = 0.d0
       tzfs3(i) = 0.d0
    end do

    lqvirial = 0.d0
    qvirial = 0.d0
    qqfs3 = 0.0D0
    ! *** converting positions in AU ***
!    boxi = box/bohr2a
    rc2_3b = rcutsd3*bohr2a**(-2.d0)
!    rc3 = (rcut3+diameter)*bohr2a**(-1.d0)
!    rc3 =9.0d0/bohr2a

!    alp3 = alpha*bohr2a

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
!$omp& private(id,is,js,ks,i,j,k,ii,iii,idx,qis,qjs,dxij,dyij,dzij,dxik,dyik,dzik,dxjk,dyjk,dzjk,ind_lin,ind_lin_2,aj)&
!$omp& private(xij,yij,zij,xik,yik,zik,xjk,yjk,zjk,dij2,dik2,djk2,xsij,ysij,zsij,dsij2,dsij,dsij1i)&
!$omp& private(xsik,ysik,zsik,dsik2,dsik,dsik1i,xsjk,ysjk,zsjk,dsjk2,dsjk,dsjk1i,dampsij,gsij,r0sij,bsij)&
!$omp& private(dampsik,gsik,r0sik,bsik,dampsjk,gsjk,r0sjk,bsjk,eksij,eksik,eksjk,eks)&
!$omp& private(init,iat,jat,kat,ip,jp,kp,new,nlin,itab,iind,nl,term1,term2,term3,der_fij,der_fik,der_fjk)&
!$omp& private(fsij,fsik,fsjk,erf_ij_ik_jk, der_erfij,der_erfjk,der_erfik,lfxs,lfys,lfzs,lfxfs3,lfyfs3,lfzfs3)&
!$omp& private(r_av,pomij,ppomij,fdmpij,dfdmpij,pomik,ppomik,fdmpik,dfdmpik,pomjk,ppomjk,fdmpjk,dfdmpjk)&
!$omp& private(xxx,yyy,zzz,xxxs,yyys,zzzs,lqqfs3)


    ! *** converting positions in AU ***
    ! *** ends here ***
   
!    allocate(ind_lin(0:lmax,0:lmax,0:lmax,ns3,ns3,ns3))
    allocate(ind_lin_2(nzero))
    allocate(aj(5000))

    allocate(lfxs(ns,nm))
    allocate(lfys(ns,nm))
    allocate(lfzs(ns,nm))
    allocate(lfxfs3(nm))
    allocate(lfyfs3(nm))
    allocate(lfzfs3(nm))
	allocate(xxx(nom));
	allocate(yyy(nom));
	allocate(zzz(nom));
	allocate(xxxs(nsite,nom));
	allocate(yyys(nsite,nom));
	allocate(zzzs(nsite,nom));
	
	lqqfs3 = 0.0D0;

    do i=1,nm
       lfxfs3(i) = 0.d0
       lfyfs3(i) = 0.d0
       lfzfs3(i) = 0.d0
       do is=1,ns
          lfxs(is,i) = 0.0D0
          lfys(is,i) = 0.0D0
          lfzs(is,i) = 0.0D0
       enddo
    enddo
	! DBG
	iniTime = SECNDS(0.0);
	! DBG

!    do i=1,nm-2
!	 do i = MPIrank + 1, nm - 2, MPIcomm_size
	if (MPIcomm_size .GT. 1) then
		i = 1;
	else
		i = 0;
	endif
	do while (i .LE. nm - 2)
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
 	  if ((i .LE. (nm - 2)) .AND. ((MPIrank .GE. 1) .OR. (MPIcomm_size .EQ. 1))) then
       !$omp do schedule(dynamic)
       do j=i+1,nm-1
          do k=j+1,nm

             dxij = xx(i)-xx(j)   ! com sep of  i-j
             dyij = yy(i)-yy(j)
             dzij = zz(i)-zz(j)

             dxik = xx(i)-xx(k)   ! com sep of  i-k
             dyik = yy(i)-yy(k)
             dzik = zz(i)-zz(k)

             dxjk = xx(j)-xx(k)   ! com sep of  j-k
             dyjk = yy(j)-yy(k)
             dzjk = zz(j)-zz(k)

             xij = dxij - boxi*nint(dxij/boxi)  ! minimum image convenion
             yij = dyij - boxi*nint(dyij/boxi)
             zij = dzij - boxi*nint(dzij/boxi)

             xik = dxik - boxi*nint(dxik/boxi)
             yik = dyik - boxi*nint(dyik/boxi)
             zik = dzik - boxi*nint(dzik/boxi)

             xjk = dxjk - boxi*nint(dxjk/boxi)
             yjk = dyjk - boxi*nint(dyjk/boxi)
             zjk = dzjk - boxi*nint(dzjk/boxi)

             dij2 = xij**2 + yij**2 + zij**2   ! rij2
             dik2 = xik**2 + yik**2 + zik**2   ! rik2
             djk2 = xjk**2 + yjk**2 + zjk**2   ! rjk2

             if ( (dij2 .lt. rc2_3b) .and. &
                  (dik2 .lt. rc2_3b) .and. &
                  (djk2 .lt. rc2_3b) ) then  ! com-com cutoff  on i-j, i-k, j-k starts

             !   rc3 = 8.81785d0/bohr2a 
               rc3 = 8.36785d0/bohr2a

                pomij = exp(alp*(sqrt(dij2)-rc3))
                ppomij = 1.d0/(1.d0+pomij)
                fdmpij = ppomij**delta3
                dfdmpij = -alp*delta3*pomij*ppomij*fdmpij !*(1/3.d0)
                dfdmpijv = -alp*delta3*pomij*ppomij*sqrt(dij2) 
 
                pomik = exp(alp*(sqrt(dik2)-rc3))
                ppomik = 1.d0/(1.d0+pomik)
                fdmpik = ppomik**delta3
                dfdmpik = -alp*delta3*pomik*ppomik*fdmpik
                dfdmpikv = -alp*delta3*pomik*ppomik*sqrt(dik2)

                pomjk = exp(alp*(sqrt(djk2)-rc3))
                ppomjk = 1.d0/(1.d0+pomjk)
                fdmpjk = ppomjk**delta3
                dfdmpjk = -alp*delta3*pomjk*ppomjk*fdmpjk 
                dfdmpjkv = -alp*delta3*pomjk*ppomjk*sqrt(djk2)

                !write(*,*)fdmpij,dfdmpij,fdmpik,dfdmpik,fdmpjk,dfdmpjk

                ! **************
                init = .true. ! GRU:This variable is useless.
                
                if(init) then
                  init = .false.
!                  ind_lin = 0
                ind_lin_2 = 0
                  new = 0  
 
                  do iat=1,ns3
                  do jat=1,ns3
                  do kat=1,ns3
                  do ip=0,lmax
                  do jp=0,lmax
                  do kp=0,lmax

                     idx = (kp + 1)
                     idx = idx + jp * (lmax + 1)
                     idx = idx + ip * (lmax + 1)**2
                     idx = idx + (kat - 1) * (lmax + 1)**3
                     idx = idx + (jat - 1) * (lmax + 1)**3 * ns3
                     idx = idx + (iat - 1) * (lmax + 1)**3 * ns3**2

!                     if (ind_lin(kp,jp,ip,kat,jat,iat).eq.0) then
                     if (ind_lin_2(idx) .EQ. 0) then
                        new=new+1
!                        call fill(ind_lin, iat,jat,kat, ip,jp,kp, new)
                        call fill(ind_lin_2, iat,jat,kat, ip,jp,kp, new)
!			ind_lin_2 = PACK(ind_lin, .TRUE.)
                     end if

                  end do
                  end do
                  end do
                  end do
                  end do
                  end do              

                end if 

                nlin = new
                aj = 0.d0 
                itab = 0
                !***************

                do is=1,ns3
                   do js=1,ns3
                      do ks=1,ns3

                         bsij  = bS3(is,js)          ! sites is-js 
                         gsij  = gS3(is,js)
                         r0sij = r0S3(is,js)

                         xsij = xij+(xxs(is,i) - xxs(js,j))
                         ysij = yij+(yys(is,i) - yys(js,j))
                         zsij = zij+(zzs(is,i) - zzs(js,j))

                         dsij2 = xsij**2 + ysij**2 + zsij**2

                         dsij = sqrt(dsij2)
                         dsij1i = 1.d0/dsij ! GRU: this variable isn't used

                         dampsij = 1.d0 /(1.d0 +exp(-gsij*(dsij-r0sij)))

                         bsik  = bS3(is,ks)          ! sites is-ks 
                         gsik  = gS3(is,ks)
                         r0sik = r0S3(is,ks)

                         xsik = xik+(xxs(is,i) - xxs(ks,k))
                         ysik = yik+(yys(is,i) - yys(ks,k))
                         zsik = zik+(zzs(is,i) - zzs(ks,k))

                         dsik2 = xsik**2 + ysik**2 + zsik**2

                         dsik = sqrt(dsik2)
                         dsik1i = 1.d0/dsik

                         dampsik = 1.d0 /(1.d0 +exp(-gsik*(dsik-r0sik)))

                         bsjk  = bS3(js,ks)          ! sites js-ks
                         gsjk  = gS3(js,ks)
                         r0sjk = r0S3(js,ks)

                         xsjk = xjk+(xxs(js,j) - xxs(ks,k))
                         ysjk = yjk+(yys(js,j) - yys(ks,k))
                         zsjk = zjk+(zzs(js,j) - zzs(ks,k))

                         dsjk2 = xsjk**2 + ysjk**2 + zsjk**2

                         dsjk = sqrt(dsjk2)
                         dsjk1i = 1.d0/dsjk

                         dampsjk = 1.d0 /(1.d0 +exp(-gsjk*(dsjk-r0sjk)))
 
                         eksij = exp(-bsij*dsij)*dampsij 
                         eksik = exp(-bsik*dsik)*dampsik
                         eksjk = exp(-bsjk*dsjk)*dampsjk  

                         eks = eksij*eksik*eksjk

                         der_fij = -bsij + gsij*exp(-gsij*(dsij-r0sij))*dampsij
                         der_fik = -bsik + gsik*exp(-gsik*(dsik-r0sik))*dampsik
                         der_fjk = -bsjk + gsjk*exp(-gsjk*(dsjk-r0sjk))*dampsjk       

                         ! virial terms
                         der_fijv = dampsij*exp(-gsij*(dsij-r0sij))*gsij*dsij-bsij*dsij
                         der_fikv = dampsik*exp(-gsik*(dsik-r0sik))*gsik*dsik-bsik*dsik
                         der_fjkv = dampsjk*exp(-gsjk*(dsjk-r0sjk))*gsjk*dsjk-bsjk*dsjk

                         do ip=0,lmax

                            term1 = eks*dsjk**ip

                            fsjk = der_fjk/dsjk + ip/dsjk**(2.d0) 

                         do jp=0,lmax

                            term2 = term1*dsik**jp

                            fsik = der_fik/dsik + jp/dsik**(2.d0)
                            
                         do kp=0,lmax

                            itab = itab + 1 

                            iind = ind_lin_2(itab) 

                            term3 = term2*dsij**kp
                            
                            fsij = der_fij/dsij + kp/dsij**(2.d0)

!                            aj(iind) = aj(iind) + term3

                            lqqfs3 = lqqfs3 + cs3(iind)*term3*fdmpij*fdmpik*fdmpjk 
   
                            lqvirial = lqvirial + ((der_fijv+der_fikv+der_fjkv+ip+jp+kp) &
                                                   +dfdmpijv + dfdmpikv + dfdmpjkv)      & 
                                                  *cs3(iind)*term3*fdmpij*fdmpik*fdmpjk

!                            lqvirial = lqvirial + (der_fijv+der_fikv+der_fjkv+ip+jp+kp) &
!                                                   *cs3(iind)*term3*fdmpij*fdmpik*fdmpjk

                            ! molecule i

                            lfxs(is,i)= lfxs(is,i) - (cs3(iind)*term3*(fsij*xsij + fsik*xsik)*fdmpij*fdmpik*fdmpjk )
                            lfys(is,i)= lfys(is,i) - (cs3(iind)*term3*(fsij*ysij + fsik*ysik)*fdmpij*fdmpik*fdmpjk ) 
                            lfzs(is,i)= lfzs(is,i) - (cs3(iind)*term3*(fsij*zsij + fsik*zsik)*fdmpij*fdmpik*fdmpjk )

                            lfxfs3(i) = lfxfs3(i) - (cs3(iind)*term3*(dfdmpij*fdmpik*fdmpjk*xij/sqrt(dij2) &
                                                                     +fdmpij*dfdmpik*fdmpjk*xik/sqrt(dik2)))

                            lfyfs3(i) = lfyfs3(i) - (cs3(iind)*term3*(dfdmpij*fdmpik*fdmpjk*yij/sqrt(dij2) &
                                                                     +fdmpij*dfdmpik*fdmpjk*yik/sqrt(dik2)))

                            lfzfs3(i) = lfzfs3(i) - (cs3(iind)*term3*(dfdmpij*fdmpik*fdmpjk*zij/sqrt(dij2) &
                                                                     +fdmpij*dfdmpik*fdmpjk*zik/sqrt(dik2)))

                            ! molecule j

                            lfxs(js,j)= lfxs(js,j) - (cs3(iind)*term3*(-fsij*xsij + fsjk*xsjk)*fdmpij*fdmpik*fdmpjk )
                            lfys(js,j)= lfys(js,j) - (cs3(iind)*term3*(-fsij*ysij + fsjk*ysjk)*fdmpij*fdmpik*fdmpjk )
                            lfzs(js,j)= lfzs(js,j) - (cs3(iind)*term3*(-fsij*zsij + fsjk*zsjk)*fdmpij*fdmpik*fdmpjk ) 

                            lfxfs3(j) = lfxfs3(j) - (cs3(iind)*term3*(-dfdmpij*fdmpik*fdmpjk*xij/sqrt(dij2) &
                                                                      +fdmpij*fdmpik*dfdmpjk*xjk/sqrt(djk2)))

                            lfyfs3(j) = lfyfs3(j) - (cs3(iind)*term3*(-dfdmpij*fdmpik*fdmpjk*yij/sqrt(dij2) &
                                                                      +fdmpij*fdmpik*dfdmpjk*yjk/sqrt(djk2)))

                            lfzfs3(j) = lfzfs3(j) - (cs3(iind)*term3*(-dfdmpij*fdmpik*fdmpjk*zij/sqrt(dij2) &
                                                                      +fdmpij*fdmpik*dfdmpjk*zjk/sqrt(djk2)))

                            ! molecule k

                            lfxs(ks,k)= lfxs(ks,k) - (cs3(iind)*term3*(-fsik*xsik - fsjk*xsjk)*fdmpij*fdmpik*fdmpjk ) 
                            lfys(ks,k)= lfys(ks,k) - (cs3(iind)*term3*(-fsik*ysik - fsjk*ysjk)*fdmpij*fdmpik*fdmpjk )
                            lfzs(ks,k)= lfzs(ks,k) - (cs3(iind)*term3*(-fsik*zsik - fsjk*zsjk)*fdmpij*fdmpik*fdmpjk )

                            lfxfs3(k) = lfxfs3(k) - (cs3(iind)*term3*(-fdmpij*dfdmpik*fdmpjk*xik/sqrt(dik2) &
                                                                      -fdmpij*fdmpik*dfdmpjk*xjk/sqrt(djk2)))

                            lfyfs3(k) = lfyfs3(k) - (cs3(iind)*term3*(-fdmpij*dfdmpik*fdmpjk*yik/sqrt(dik2) &
                                                                      -fdmpij*fdmpik*dfdmpjk*yjk/sqrt(djk2)))

                            lfzfs3(k) = lfzfs3(k) - (cs3(iind)*term3*(-fdmpij*dfdmpik*fdmpjk*zik/sqrt(dik2) &
                                                                      -fdmpij*fdmpik*dfdmpjk*zjk/sqrt(djk2)))  

                         end do 
                         end do
                         end do 

                         !end if ! check whether r_avg is with in rcut3 !! ends

                      end do  ! ks
                   end do     ! js
                end do        ! is 

!                do nl=1,nlin
!                qqfs3 = qqfs3 + cs3(nl)*aj(nl)
!                end do

            end if  ! com-com cutoff ends

          end do             ! k
       end do                ! j
       !$omp end do
	  endif
    enddo                   ! i
	if (MPIrank .EQ. 0) then
		do i = 1, MPIcomm_size - 1
			call MPI_RECV(MPIaddr, 1, MPI_INTEGER, MPI_ANY_SOURCE, MPI_TAG_CommReady, MPI_COMM_WORLD, MPIstat, MPIierr);
			call MPI_SEND(nm + 1, 1, MPI_INTEGER, MPIaddr, MPI_TAG_SendData, MPI_COMM_WORLD, MPIierr);
		enddo
	endif

	! DBG
	deltaTime = SECNDS(iniTime);
!	write(dbgUnit, *) 'P#', MPIrank, '. FS3 calc time	', deltaTime;
	! DBG

!$omp critical
    do i=1,nm
       fxfs3(i) = fxfs3(i) + lfxfs3(i)*h2kcal*418.4d0/bohr2a
       fyfs3(i) = fyfs3(i) + lfyfs3(i)*h2kcal*418.4d0/bohr2a
       fzfs3(i) = fzfs3(i) + lfzfs3(i)*h2kcal*418.4d0/bohr2a
       do is=1,ns
          llfxs(is,i) = llfxs(is,i) + lfxs(is,i)*h2kcal*418.4d0/bohr2a
          llfys(is,i) = llfys(is,i) + lfys(is,i)*h2kcal*418.4d0/bohr2a
          llfzs(is,i) = llfzs(is,i) + lfzs(is,i)*h2kcal*418.4d0/bohr2a
       enddo
    enddo
        qvirial = qvirial + lqvirial
	qqfs3 = qqfs3 + lqqfs3;
!$omp end critical

!    deallocate(ind_lin)
    deallocate(ind_lin_2)
    deallocate(aj)
    deallocate(lfxs)
    deallocate(lfys)
    deallocate(lfzs)
    deallocate(lfxfs3)
    deallocate(lfyfs3)
    deallocate(lfzfs3)
	deallocate(xxx);
	deallocate(yyy);
	deallocate(zzz);
	deallocate(xxxs);
	deallocate(yyys);
	deallocate(zzzs);
!$omp end parallel

    qfs3 = h2kcal*qqfs3
!    write(*,*) qfs3
! MPI sync
	j = nm * ns * 3 + nm * 3 + 2; !fxs, fys, fzs, fxfs3, fyfs3, fzfs3,qfs3,qvirial
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
		MPIbuff1(k) = fxfs3(i); k = k + 1;
		MPIbuff1(k) = fyfs3(i); k = k + 1;
		MPIbuff1(k) = fzfs3(i); k = k + 1;
	enddo
	MPIbuff1(k) = qfs3; k = k + 1;
        MPIbuff1(k) = qvirial; k = k + 1;
	call MPI_ALLREDUCE(MPIbuff1, MPIbuff2, j, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, MPIierr);
	k = 0;
	do i = 1, nm
		do is = 1, ns
			fxs(is, i) = fxs(is, i) + MPIbuff2(k); k = k + 1;
			fys(is, i) = fys(is, i) + MPIbuff2(k); k = k + 1;
			fzs(is, i) = fzs(is, i) + MPIbuff2(k); k = k + 1;
		enddo
		fxfs3(i) = MPIbuff2(k); k = k + 1;
		fyfs3(i) = MPIbuff2(k); k = k + 1;
		fzfs3(i) = MPIbuff2(k); k = k + 1;
	enddo
	qfs3 = MPIbuff2(k); k = k + 1;
        qvirial = MPIbuff2(k); k = k + 1;
	deallocate(MPIbuff1);
	deallocate(MPIbuff2);
! MPI sync
    deallocate(llfxs);
    deallocate(llfys);
    deallocate(llfzs);
	! DBG
	! if (MPIrank .EQ. 0) then
		! write(303, '(a)') '#';
		! write(303, '(a)') 'fxs';
		! write(303, '(1p3E16.6E3)') fxs;
		! write(303, '(a)') 'fys';
		! write(303, '(1p3E16.6E3)') fys;
		! write(303, '(a)') 'fzs';
		! write(303, '(1p3E16.6E3)') fzs;
		! write(303, '(a)') 'fxfs3';
		! write(303, '(1p3E16.6E3)') fxfs3;
		! write(303, '(a)') 'fyfs3';
		! write(303, '(1p3E16.6E3)') fyfs3;
		! write(303, '(a)') 'fzfs3';
		! write(303, '(1p3E16.6E3)') fzfs3;
		! write(303, '(a)') 'qfs3';
		! write(303, '(1p3E16.6E3)') qfs3;
	! endif
	! DBG

    fs3_vir = qvirial*h2kcal*418.4d0 
    vir = vir + fs3_vir 

    return;
  end subroutine FS3_term

  subroutine fill (ntab,iat,jat,kat,ip,jp,kp,nr)

    implicit none

    integer :: iat,jat,kat,ip,jp,kp,nr,is
    integer :: ifrst,ilast,jfrst,jlast,kfrst,klast 
    integer :: ii,jj,kk
    logical :: ilook,jlook,klook 
    integer, parameter :: nsite3=17, lmax=1
    integer :: ntab(0:lmax,0:lmax,0:lmax,nsite3,nsite3,nsite3)   
    integer :: ntype(nsite3)

    data ntype /1,2,2,3,3,4,4,4,4,5,5,5,5,6,6,6,6/  ! must be blocks with consecutive numbers

    ilook=.true.
    jlook=.true.
    klook=.true.

    do is=1,nsite3
         if ((ntype(iat).eq.ntype(is)).and.ilook) then
            ifrst=is
            ilook=.false.
         end if
         if ((ntype(jat).eq.ntype(is)).and.jlook) then
            jfrst=is
            jlook=.false.
         end if
         if ((ntype(kat).eq.ntype(is)).and.klook) then
            kfrst=is
            klook=.false.
         end if
         if (ntype(iat).eq.ntype(is)) then
            if (is.eq.nsite3) then
               ilast=is
            else
               if (ntype(iat).ne.ntype(is+1)) ilast=is
            end if
         end if
         if (ntype(jat).eq.ntype(is)) then
            if (is.eq.nsite3) then
               jlast=is
            else
               if (ntype(jat).ne.ntype(is+1)) jlast=is
            end if
         end if
         if (ntype(kat).eq.ntype(is)) then
            if (is.eq.nsite3) then
               klast=is
            else
               if (ntype(kat).ne.ntype(is+1)) klast=is
            end if
         end if
    end do

    do ii=ifrst,ilast
    do jj=jfrst,jlast
    do kk=kfrst,klast
         ntab(kp,jp,ip, kk,jj,ii)=nr
         ntab(jp,kp,ip, jj,kk,ii)=nr
         ntab(ip,kp,jp, ii,kk,jj)=nr
         ntab(kp,ip,jp, kk,ii,jj)=nr
         ntab(ip,jp,kp, ii,jj,kk)=nr
         ntab(jp,ip,kp, jj,ii,kk)=nr
    end do
    end do
    end do  

  end subroutine fill

end module threebody_potential

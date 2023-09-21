module vdw_ccpol23plus

  use, intrinsic :: iso_fortran_env
  use common_arrays
  use utility 
  use ewald
  
  implicit none
  private

  public :: ccpol23plus, ccpol23plus_lrc

contains

  subroutine ccpol23plus(xx,yy,zz,xxs,yys,zzs,boxi,qqvdwe)
    implicit none
    real(8) :: xx(nom),yy(nom),zz(nom)
    real(8) :: xxs(nsite,nom),yys(nsite,nom),zzs(nsite,nom)
    real(8) :: boxi,qqvdwe

    integer :: i,j,is,js,k
    real(8) :: wij,qis,qjs
    real(8) :: xcc,ycc,zcc
    real(8) :: dx,dy,dz
    real(8) :: rccsq
    real(8) :: rssq,rsx,rsy,rsz
    real(8) :: vfij
    real(8) :: sqr
    real(8) :: eks,u_exp,du_exp
    real(8) :: u_ele, f1,df1, du_ele
    real(8) :: u_ind,f6,f8,f10,df6,df8,df10
    real(8) :: r6,r8,r10,du_ind
    real(8) :: fdamp
    real(8) :: fxsij,fysij,fzsij
    ! Shared
    real(8), allocatable :: lfxs(:, :), lfys(:, :), lfzs(:, :)

    real(8) :: qvdwe, qvirial, qvir_vdw
    real(8) :: lqvdwe, lqvirial,lvdwlrc,lvirlrc
    real(8), parameter :: gamma = 20.d0/bohr2a
    real(8), parameter :: sqr0 = 1.4d0*bohr2a  
	! MPI
	INTEGER(4) MPIchunk_size_loc, iii;
	real(8), allocatable :: MPIbuff1(:), MPIbuff2(:);
	! DBG
	real(4) iniTime, deltaTime;
	! DBG
	! MPI

    qvdwe=0.d0
    qvirial=0.d0 ! GRU: this variable seems to be not used

    vir = 0.d0  ! global variable for virial
!    vir_test = 0.d0


    do  j=1,nm
      do js=1,ns
         fxs(js,j)=0.0
         fys(js,j)=0.0
         fzs(js,j)=0.0
      end do
    end do

!$omp parallel DEFAULT(SHARED)&
!$omp& private(i,iii,j,is,js,k)&
!$omp& private(wij,qis,qjs,xcc,ycc,zcc,dx,dy,dz,rccsq)&
!$omp& private(rssq,rsx,rsy,rsz,vfij,sqr,eks,u_exp,du_exp)&
!$omp& private(u_ele, f1,df1, du_ele,u_ind,f6,f8,f10,df6,df8,df10)&
!$omp& private(r6,r8,r10,du_ind,fdamp,lfxs,lfys,lfzs,lqvdwe, lqvirial)

    lqvdwe=0.d0
    lqvirial=0.d0

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
!reduction(+:qvdwe,qvirial)
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

          if (rccsq.lt.rcutsd) then   ! rc com

             do is=1,ns

                do js=1,ns
                 
                  rsx = xcc + xxs(is,i)-xxs(js,j)
                  rsy = ycc + yys(is,i)-yys(js,j)
                  rsz = zcc + zzs(is,i)-zzs(js,j) 
                  
                  qis = chgs(is,i)
                  qjs = chgs(js,j)

                  rssq = rsx**2 + rsy**2 + rsz**2
            
!                  sqr = sqrt(rssq)

!                  if (rssq.le.rcutsd) then    ! rc sites
 
                     sqr = sqrt(rssq) ! square root of r^2

                     do k=1,npi
                        if( (an(is) .eq. an1(k) .and. an(js) .eq. an2(k)) .or. &
                          (an(js) .eq. an1(k) .and. an(is) .eq. an2(k)) )then         
 
!             ccpol8s Energy, Virial and Forces                    
                          ! exponential (excahnge - repulsion) 
                          eks = exp(-beta(k)*sqr)  !  exp(-beta r)
                          u_exp = eks*(c0(k)+c1(k)*sqr+c2(k)*sqr**2+c3(k)*sqr**3)

                          !if (u_exp .lt. 0.d0) then

                          fdamp = 1.d0/( 1.d0+exp(-gamma*(sqr-sqr0))) ! short range damping function
                          u_exp = u_exp * fdamp   ! u_exp in atomic units       

                          !end if    

                          ! short range coulumb contribution
                          f1 = tt(1, d1(k), sqr)
                          u_ele = r4pie0*(f1-1.d0)*qis*qjs/sqr

                          ! induction-dispersion
                          f6 = tt(6, d6(k), sqr) 
                          f8 = tt(8, d8(k), sqr)
                          f10 = tt(10, d10(k), sqr) 
                          r6 = rssq*rssq*rssq
                          r8 = rssq*r6
                          r10 = rssq*r8     
                          u_ind = - f6*c6(k)/r6 - f8*c8(k)/r8 - f10*c10(k)/r10 ! u_ind in atomic units
                          
                          lqvdwe = lqvdwe + u_exp + u_ele + u_ind 

                          ! drivative of echange repulsion 
                          du_exp = - beta(k)*u_exp*fdamp &
                                  + eks*fdamp*(c1(k)+2.d0*c2(k)*sqr+3.d0*c3(k)*sqr**2) &
                                  + gamma*exp(-gamma*(sqr-sqr0))*fdamp**2*u_exp 
                                   
                          ! derivative of electrostatic part (du_ele/drij)
                          df1=dtt(1, d1(k),sqr)
                          du_ele=-r4pie0*qis*qjs*((f1-1.d0)/sqr**2 - df1/sqr)

                          
                          ! derivative of asymptotic part of the dispersion interaction (du_ind/drij)
                          df6=dtt(6, d6(k), sqr)
                          df8=dtt(8, d8(k), sqr)
                          df10=dtt(10, d10(k), sqr) 
                          du_ind = - df6*c6(k)/r6 + 6.d0*f6*c6(k)/sqr**7 &
                                  - df8*c8(k)/r8 + 8.d0*f8*c8(k)/sqr**9 & 
                                  - df10*c10(k)/r10 + 10.d0*f10*c10(k)/sqr**11 

                          wij = du_exp+du_ele+du_ind

!                          lqvirial = lqvirial + wij

                          vfij = -wij/sqr

!              Accumalating focces

                          fxsij = rsx*vfij
                          fysij = rsy*vfij
                          fzsij = rsz*vfij
      
                          lfxs(is,i)=lfxs(is,i)+ fxsij
                          lfys(is,i)=lfys(is,i)+ fysij
                          lfzs(is,i)=lfzs(is,i)+ fzsij

                          lfxs(js,j)=lfxs(js,j)- fxsij
                          lfys(js,j)=lfys(js,j)- fysij
                          lfzs(js,j)=lfzs(js,j)- fzsij

                  ! virial = - rab * fab

                           lqvirial = lqvirial - (rsx*fxsij + rsy*fysij + rsz*fzsij)

                        end if
                     end do

!                  end if    ! rc sites
 
                end do 

             end do
              
          end if            ! rc com      
           
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
!	write(dbgUnit, *) 'P#', MPIrank, '. ccpol23plus calc time	', deltaTime;
	! DBG
!$omp critical
     qvdwe = qvdwe + lqvdwe
     qvirial = qvirial + lqvirial

     do j=1,nm
        do js=1,ns
           fxs(js,j) = fxs(js,j) + lfxs(js,j)
           fys(js,j) = fys(js,j) + lfys(js,j)
           fzs(js,j) = fzs(js,j) + lfzs(js,j)
        end do
     end do
!$omp end critical

     deallocate(lfxs)
     deallocate(lfys)
     deallocate(lfzs)
!$omp end parallel

! MPI sync
	j = nm * ns * 3 + 2; !fxs, fys, fzs, qvdwe, qvirial
	allocate(MPIbuff1(0:j - 1));
	allocate(MPIbuff2(0:j - 1));
	do k = 0, j - 1
		MPIbuff2(k) = 0.0D0;
	enddo
	k = 0;
	do i = 1, nm
		do is = 1, ns
			MPIbuff1(k) = fxs(is, i); k = k + 1;
			MPIbuff1(k) = fys(is, i); k = k + 1;
			MPIbuff1(k) = fzs(is, i); k = k + 1;
		enddo
	enddo
	MPIbuff1(k) = qvdwe; k = k + 1;
	MPIbuff1(k) = qvirial; k = k + 1;
	call MPI_ALLREDUCE(MPIbuff1, MPIbuff2, j, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, MPIierr);
	k = 0;
	do i = 1, nm
		do is = 1, ns
			fxs(is, i) = MPIbuff2(k); k = k + 1;
			fys(is, i) = MPIbuff2(k); k = k + 1;
			fzs(is, i) = MPIbuff2(k); k = k + 1;
		enddo
	enddo
	qvdwe = MPIbuff2(k); k = k + 1;
	qvirial = MPIbuff2(k); k = k + 1;
	deallocate(MPIbuff1);
	deallocate(MPIbuff2);
! MPI sync

     lvdwlrc = vdwlrc
     lvirlrc = virlrc
     call ccpol23plus_lrc(nm,box,rcut,lvdwlrc,lvirlrc)
     vdwlrc = lvdwlrc/418.4d0
     virlrc = lvirlrc/418.4d0

!     call ccpol23plus_lrc(nm,box,rcut,vdwlrc,virlrc)
!     vdwlrc = vdwlrc/418.4d0
!     virlrc = virlrc/418.4d0

     qvdwe = qvdwe/418.4d0
     qqvdwe = qvdwe !+ vdwlrc  

!     write(*,*) qvirial

     qvir_vdw = qvirial !+ lvirlrc

     vir = vir + qvir_vdw

!     vir_test = vir_test+qvir_vdw

!     write(*,*) 'vdw: ', vir,vir_test,qvir_vdw


  end subroutine ccpol23plus

  subroutine ccpol23plus_lrc(n,bx,rc,ulrc,wlrc)
    implicit none
    integer :: n
    real(8) :: bx,rc
    real(8) :: ulrc ! pot eng correction
    real(8) :: wlrc ! virial correct

    integer :: m
    real(8) :: rct
    real(8) :: uedis,uele,uin6,uin8,uin10
    real(8) :: wedis,wele,win6,win8,win10
    real(8) :: lulrc, lwlrc

    real(8) :: vol,dens
    real(8) :: ulrcs,wlrcs
    character(len=5) :: lrcs ! GRU: Local variable never used

    vol    = bx**3.d0
    dens   = real(n)/vol      

    ! ** calculates long range corrections 
! dlpoly_classical manual 1.10 page 30 equation 2.105

    ulrcs=0.d0 ! GRU: Local variable never used.
    wlrcs=0.d0 ! GRU: Local variable never used.

    ulrc=0.d0
    wlrc=0.d0

! GRU: parallel writing into a file will require synchronization, comment unless really necessary
!    write(5000,*)'a ',' b ','     ulrc(kcal/mol)        ','     plrc(katm)', &
!              '                ulrc_nc            ', '    plrc_nc'
!$omp parallel DEFAULT(SHARED)&
!$omp& private(m,rct,uedis,uele,uin6,uin8,uin10,wedis,wele,win6,win8,win10)&
!$omp& private(lulrc,lwlrc)
    lulrc = 0.0d0
    lwlrc = 0.0d0

!$omp do schedule(dynamic)
    do m=1,npi

      if(abs(beta(m)) > 0.d0) then
!      if(beta(m) .ne. 0.d0) then
      uedis =(120*c3(m)+24*beta(m)*(c2(m)+5*c3(m)*rc)+6*beta(m)**2*(c1(m)+2*rc*(2*c2(m) &
               +5*c3(m)*rc))+beta(m)**4*rc*(2*c0(m)+rc*(3*c1(m)+4*c2(m)*rc+5*c3(m)*rc**2)) &
               +2*beta(m)**3*(c0(m)+rc*(3*c1(m)+6*c2(m)*rc+10*c3(m)*rc**2))+beta(m)**5*rc**2 &
               *(c0(m)+rc*(c1(m)+rc*(c2(m)+c3(m)*rc))))/(beta(m)**6*exp(beta(m)*rc))
      else

      uedis = 0.d0

      end if

      if(abs(d1(m)) > 0.d0) then        ! r4pie0*
!      if(d1(m) .ne. 0.d0) then
      uele = - r4pie0*((qa(m)*qb(m)*(3 + d1(m)*rc*(3+d1(m)*rc)))/(d1(m)**2*exp(d1(m)*rc)))

      else
      uele = 0.d0
      end if

      uin6 = (c6(m)*(240 - 240*exp(d6(m)*rc) + d6(m)*rc*(240 + d6(m)*rc*(120+d6(m)*rc*(38+ &
              d6(m)*rc*(8 + d6(m)*rc))))))/(720.*exp(d6(m)*rc)*rc**3)
!      write(*,*) uin6

      uin8 =(c8(m)*(8064-8064*exp(d8(m)*rc)+d8(m)*rc*(8064+d8(m)*rc*(4032+d8(m)*rc*(1344+  &
             d8(m)*rc*(336+d8(m)*rc*(66+d8(m)*rc*(10+d8(m)*rc))))))))/(40320.* &
             exp(d8(m)*rc)*rc**5)
!      write(*,*) uin8

      uin10 = -(c10(m)*(-518400+518400*exp(d10(m)*rc)-d10(m)*rc*(518400+d10(m)*rc*   &
              (259200+d10(m)*rc*(86400 + d10(m)*rc*(21600 + d10(m)*rc*(4320 + &
              d10(m)*rc*(720+d10(m)*rc*(102+d10(m)*rc*(12+d10(m)*rc)))))))))) &
              /(3.6288e6*exp(d10(m)*rc)*rc**7)
!      write(*,*) uin10

!       write(*,*) an1(m),an2(m),(uedis+uele+uin6+uin8+uin10)/418.4d0
       lulrc = lulrc + uedis + uele + uin6 + uin8 + uin10
!       write(*,*) an1(m),an2(m),ulrc/418.4d0
   ! conversion factor for pressure from internal units to katm is 0.163882576
! 0.163882576 * WLRC

! following is LRC for pressure, definition is taken fron dlpoly_classical
! manual 1.10 page 30 equation 2.106 divide by 3*volume to get pressure
! corrections

       if(abs(beta(m)) > 0.d0) then
       wedis= -((360*c3(m)+72*beta(m)*(c2(m)+5*c3(m)*rc)+18*beta(m)**2*(c1(m) &
              + 2*rc*(2*c2(m) +5*c3(m)*rc))+3*beta(m)**4*rc*(2*c0(m)+rc*(3*c1(m) &
              + 4*c2(m)*rc+5*c3(m)*rc**2))+6*beta(m)**3*(c0(m)+rc*(3*c1(m) + &
               6*c2(m)*rc+10*c3(m)*rc**2))+3*beta(m)**5*rc**2*(c0(m)+rc*(c1(m) &
              + rc*(c2(m)+c3(m)*rc)))+beta(m)**6*rc**3*(c0(m)+rc*(c1(m) + &
              rc*(c2(m)+c3(m)*rc))))/(beta(m)**6*exp(beta(m)*rc)))
       else

       wedis = 0.d0

       end if

      if(abs(d1(m)) > 0.d0) then
      wele=r4pie0*(qa(m)*qb(m)*(9+d1(m)*rc*(9+d1(m)*rc*(4+d1(m)*rc))))/(d1(m)**2*exp(d1(m)*rc))

      else

      wele = 0.d0

      end if

      win6 = -(c6(m)*(1440-1440*exp(d6(m)*rc)+d6(m)*rc*(1440+d6(m)*rc*(720 &
             + d6(m)*rc*(234+d6(m)*rc*(54+d6(m)*rc*(9+d6(m)*rc))))))) &
             /(720.*exp(d6(m)*rc)*rc**3)

      win8 = -(c8(m)*(64512-64512*exp(d8(m)*rc)+d8(m)*rc*(64512+d8(m)*rc*(32256 &
            + d8(m)*rc*(10752+d8(m)*rc*(2688+d8(m)*rc*(534+d8(m)*rc* &
             (86 + d8(m)*rc*(11 + d8(m)*rc)))))))))/(40320.*exp(d8(m)*rc)*rc**5)      

      win10 = (c10(m)*(-5184000+5184000*exp(d10(m)*rc)-d10(m)*rc*(5184000 + &
             d10(m)*rc*(2592000+d10(m)*rc*(864000+d10(m)*rc*(216000 &
             + d10(m)*rc*(43200+d10(m)*rc*(7200+d10(m)*rc*(1026 &
             + d10(m)*rc*(126+d10(m)*rc*(13+d10(m)*rc)))))))))))/ &
                (3.6288e6*exp(d10(m)*rc)*rc**7)

       lwlrc =  lwlrc + wedis+wele+win6+win8+win10

!	   write(5000,*) an1(m),an2(m),&
!                     (2.d0*pi*n**2/vol)*(uedis+uele+uin6+uin8+uin10)/418.4d0, &
!                     (0.163882576*(-(2.d0*pi*n**2/vol)* &
!                     (wedis+wele+win6+win8+win10))/(3.d0*vol))/418.4d0, &
!                     (uedis+uele+uin6+uin8+uin10)/418.4d0,&
!                     (wedis+wele+win6+win8+win10)/418.4d0

    end do
!$omp end do

!$omp critical
    wlrc = wlrc + lwlrc
    ulrc = ulrc + lulrc
 !$omp end critical
!$omp end parallel

    ulrc = (2.d0*pi*n**2/vol) * ulrc
    wlrc = 0.163882576*(-(2.d0*pi*n**2/vol) * wlrc)/(3.d0*vol)

  end subroutine ccpol23plus_lrc 

end module vdw_ccpol23plus 

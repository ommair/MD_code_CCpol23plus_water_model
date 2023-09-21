module vdw_ccpol8s_omo

  use, intrinsic :: iso_fortran_env
  use common_arrays
  use utility 
  use ewald
  
  implicit none
  private

  public :: ccpol8s_omo, ccpol8s_omo_lrc

contains

  subroutine ccpol8s_omo(xx,yy,zz,xxs,yys,zzs,boxi,qqvdwe)

    implicit none
 
    real(8) :: xx(nom),yy(nom),zz(nom)
    real(8) :: xxs(nsite,nom),yys(nsite,nom),zzs(nsite,nom)
    real(8) :: boxi   

    integer :: i,j,is,js,k
    real(8) :: qqvdwe,qvdwe, qvirial,qvir_vdw
    real(8) :: wij,qis,qjs
    real(8) :: xcc,ycc,zcc
    real(8) :: dx,dy,dz
    real(8) :: rccsq
    real(8) :: rssq, rsx,rsy,rsz
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
    real(8) :: lqvdwe, lqvirial,lvdwlrc,lvirlrc
	! MPI
	INTEGER(4) MPIchunk_size_loc, iii;
	real(8), allocatable :: MPIbuff1(:), MPIbuff2(:);
	! DBG
	real(4) iniTime, deltaTime;
	! DBG
	! MPI

    qvdwe=0.d0
    qvirial=0.d0

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
!$omp& private(i,j,is,js,k)&
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

          if (rccsq.lt.rcutsd) then  ! rc com

             do is=1,ns

                do js=1,ns
                 
                  rsx = xcc + xxs(is,i)-xxs(js,j)
                  rsy = ycc + yys(is,i)-yys(js,j)
                  rsz = zcc + zzs(is,i)-zzs(js,j) 
                  
                  qis = chgs(is,i)
                  qjs = chgs(js,j)

                  rssq = rsx**2 + rsy**2 + rsz**2
!                  sqr = sqrt(rssq)

!                  if (rssq.lt.rcutsd) then  ! rc sites
 
                     sqr = sqrt(rssq) ! square root of r^2

                     do k=1,npi
                        if( (an(is) .eq. an1(k) .and. an(js) .eq. an2(k)) .or. &
                          (an(js) .eq. an1(k) .and. an(is) .eq. an2(k)) )then           
 
!             ccpol8s_omo Energy, Virial and Forces                   
         
                          ! exponential (excahnge - repulsion) 
                          eks = exp(-beta(k)*sqr)  !  exp(-beta r)
                          u_exp = eks*(c0(k) + c1(k)*sqr + c2(k)*sqr**2)    

                          ! induction-dispersion
                          f6 = tt(6, domo(k), sqr) 
                          f8 = tt(8, domo(k), sqr) 
                          r6 = rssq*rssq*rssq
                          r8 = rssq*r6
                          u_ind =  f6*c6(k)/r6 + f8*c8(k)/r8   !! check sign if required
                          
                          lqvdwe = lqvdwe + u_exp + u_ind

        !                  write(*,*) u_exp , u_ind 

                          ! drivative of echange repulsion (du_exp/drij) 
                          du_exp=-beta(k)*u_exp + eks*(c1(k)+2.d0*c2(k)*sqr)
                          
                          ! derivative of asymptotic part of the dispersion interaction (du_ind/drij)
                          df6=dtt(6, domo(k), sqr)
                          df8=dtt(8, domo(k), sqr)
                          du_ind =  df6*c6(k)/r6 - 6.d0*f6*c6(k)/sqr**7 &
                                  + df8*c8(k)/r8 - 8.d0*f8*c8(k)/sqr**9 

!                          wij = du_exp+du_ele+du_ind ! GRU: du_ele is not assigned.
                          wij = du_exp+du_ind ! GRU: du_ele is not assigned.

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

                           qvirial = qvirial + wij*sqr

!                           qvirial = qvirial - (  rsx*fxsij &
!                                                + rsy*fysij &
!                                                + rsz*fzsij)

                        end if
                     end do

!                  end if   ! rc sites
 
                end do 

             end do
              
          end if         ! rc com
           
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
!	write(dbgUnit, *) 'P#', MPIrank, '. ccpol8s_omo calc time	', deltaTime;
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
     call ccpol8s_omo_lrc(nm,box,rcut,lvdwlrc,lvirlrc)
     vdwlrc = lvdwlrc/418.4d0
     virlrc = lvirlrc*prsunt/(-3.d0*volm) !/418.4d0
!     qvdwe = qvdwe/418.4d0
!     vdwlrc = vdwlrc/418.4d0
!     virlrc = virlrc/418.4d0

     qvdwe = qvdwe/418.4d0
     qqvdwe = qvdwe !+ vdwlrc  

!     write(*,*) qvirial

     qvir_vdw = qvirial !+ lvirlrc

     vir = vir + qvir_vdw

!     vir_test = vir_test+qvir_vdw

!     write(*,*) 'vdw: ', vir,qvir_vdw     
 
	return;
  end subroutine ccpol8s_omo

  subroutine ccpol8s_omo_lrc(n,bx,rc,ulrc,wlrc)

    implicit none
    integer :: n,m
    real(8) :: ulrc ! pot eng correction
    real(8) :: wlrc ! virial correct
    real(8) :: vol, dens, bx
    real(8) :: rc,rct
    real(8) :: uedis,uin6,uin8,uin10
    real(8) :: wedis,win6,win8,win10
    real(8) :: lulrc, lwlrc

    character(len=5) :: lrcs
    real(8) :: ulrcs,wlrcs
 

    vol    = bx**3.d0
    dens   = real(n)/vol      

    ! ** calculates long range corrections ** for LJ POTENTIAL
! dlpoly_classical manual 1.10 page 30 equation 2.105

!analytical expressions of lrc integrals are calculted using mathemtica
    ulrcs=0.d0
    wlrcs=0.d0

    ulrc=0.d0
    wlrc=0.d0

!    write(5000,*)'a ',' b ','     ulrc(kcal/mol)        ','     plrc(katm)', &
!              '                ulrc_nc            ', '    plrc_nc'
!$omp parallel DEFAULT(SHARED)&
!$omp& private(m,uedis,uin6,uin8,uin10,wedis,win6,win8,win10)&
!$omp& private(lulrc,lwlrc)
    lulrc = 0.0d0
    lwlrc = 0.0d0

!$omp do schedule(dynamic)
    do m=1,npi

       uedis = (24*c2(m) + 6*beta(m)*(c1(m) + 4*c2(m)*rc) + beta(m)**4*rc**2*(c0(m)    &      
               + rc*(c1(m)+c2(m)*rc))+ 2*beta(m)**2*(c0(m) + 3*rc*(c1(m) + 2*c2(m)*rc))     &
               + beta(m)**3*rc*(2*c0(m) +rc*(3*c1(m)+4*c2(m)*rc)))/(beta(m)**5*exp(beta(m)*rc)) 

       uin6 = -(c6(m)*(240 - 240*exp(domo(m)*rc)+ domo(m)*rc*(240 + domo(m)*rc*(120     &
               + domo(m)*rc*(38 + domo(m)*rc*(8 +domo(m)*rc))))))/(720.*exp(domo(m)*rc)*rc**3)  

       uin8 = -(c8(m)*(8064 - 8064*exp(domo(m)*rc) + domo(m)*rc*(8064 + domo(m)*rc*(4032 &
               + domo(m)*rc*(1344+domo(m)*rc*(336 + domo(m)*rc*(66 + domo(m)*rc*(10 &
               + domo(m)*rc))))))))/(40320.*exp(domo(m)*rc)*rc**5)

       lulrc = lulrc + uedis + uin6 + uin8

   ! conversion factor for pressure from internal units to katm is 0.163882576
! 0.163882576 * WLRC

! following is LRC for pressure, definition is taken fron dlpoly_classical
! manual 1.10 page 30 equation 2.106 divide by 3*volume to get pressure
! corrections 
       wedis = -((72*c2(m)+18*beta(m)*(c1(m)+4*c2(m)*rc)+3*beta(m)**4*rc**2*(c0(m) &   
                    +rc*(c1(m)+c2(m)*rc)) +beta(m)**5*rc**3*(c0(m)+rc*(c1(m)+c2(m)*rc)) &
                    +6*beta(m)**2*(c0(m)+3*rc*(c1(m)+2*c2(m)*rc))+3*beta(m)**3*rc*(2*c0(m)  &
                    +rc*(3*c1(m)+4*c2(m)*rc)))/(beta(m)**5*exp(beta(m)*rc)))

       win6 = (c6(m)*(1440-1440*exp(domo(m)*rc)+domo(m)*rc*(1440+domo(m)*rc*(720+domo(m)*rc*(234  &
              +domo(m)*rc*(54+domo(m)*rc*(9 + domo(m)*rc)))))))/(720.*exp(domo(m)*rc)*rc**3)          
                 
       win8 = +(c8(m)*(64512-64512*exp(domo(m)*rc)+domo(m)*rc*(64512+domo(m)*rc*(32256       &
              +domo(m)*rc*(10752+domo(m)*rc*(2688+domo(m)*rc*(534+domo(m)*rc*(86+domo(m)*rc  & 
              *(11+domo(m)*rc)))))))))/(40320.*exp(domo(m)*rc)*rc**5)
       
       lwlrc = lwlrc + wedis + win6 + win8  

!       write(5000,*) an1(m),an2(m),&
!                     (2.d0*pi*n**2/vol)*(uedis+uin6+uin8)/418.4d0, &
!                     (0.163882576*(-(2.d0*pi*n**2/vol)* &
!                     (wedis+win6+win8))/(3.d0*vol))/418.4d0, &
!                     (uedis+uin6+uin8)/418.4d0,&
!                     (wedis+win6+win8)/418.4d0
 
    end do 
!$omp end do

!$omp critical
 wlrc = wlrc + lwlrc
 ulrc = ulrc + lulrc
 !$omp end critical
!$omp end parallel      

   ulrc = (2.d0*pi*n**2/vol) * ulrc
   wlrc = -(2.d0*pi*n**2/vol) * wlrc

!   wlrc = 0.163882576*(-(2.d0*pi*n**2/vol) * wlrc)/(3.d0*vol)

  end subroutine ccpol8s_omo_lrc 


end module vdw_ccpol8s_omo 

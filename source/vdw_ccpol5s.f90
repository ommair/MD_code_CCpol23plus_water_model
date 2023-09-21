module vdw_ccpol5s

  use, intrinsic :: iso_fortran_env
  use common_arrays
  use utility 
  use ewald
  
  implicit none
  private

  public :: ccpol5s, ccpol5s_lrc

contains

  subroutine ccpol5s(xx,yy,zz,xxs,yys,zzs,boxi,qvdwe)

    implicit none

    real(8) :: xx(nom),yy(nom),zz(nom)
    real(8) :: xxs(nsite,nom),yys(nsite,nom),zzs(nsite,nom)
    real(8) :: boxi, qvdwe

    integer :: i,j,is,js,k
    real(8) :: qvirial,qvir_vdw
    real(8) :: wij
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
    real(8) :: fxsij,fysij,fzsij
!    real(8) :: tt
    real(8) :: qis,qjs
	! MPI
    real(8) :: lvdwlrc,lvirlrc
	INTEGER(4) MPIchunk_size_loc, iii;
	real(8), allocatable :: MPIbuff1(:), MPIbuff2(:);
	! DBG
	real(4) iniTime, deltaTime;
	! DBG
	! MPI

    qvdwe=0.d0
    qvirial=0.d0

    vir = 0.d0  ! global variable for virial
    vir_test = 0.d0

    do  j=1,nm
      do js=1,ns
         fxs(js,j)=0.0
         fys(js,j)=0.0
         fzs(js,j)=0.0
      end do
    end do

	! DBG
	iniTime = SECNDS(0.0);
	! DBG

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

!                  if (rssq.le.rcutsd) then  ! rc sites
 
                     sqr = sqrt(rssq)

                     do k=1,npi
                        if( (an(is) .eq. an1(k) .and. an(js) .eq. an2(k)) .or. &
                          (an(js) .eq. an1(k) .and. an(is) .eq. an2(k)) )then  
 
                          !qis = chgs(is,i)
                          !qjs = chgs(js,j)              
 
!             ccpol5s Energy, Virial and Forces                   
                          !write(*,*) an1(k),an2(k), i,is,j,js,qis, qjs 
                          ! exponential (excahnge - repulsion) 
                          
                          eks = exp(-beta(k)*sqr)  !  exp(-beta r)
                          u_exp = eks*(c0(k)+c1(k)*sqr+c2(k)*sqr**2+c3(k)*sqr**3)
                          !u_exp = eks*expalp(k)*(1.d0 + a1(k)*sqr + a2(k)*sqr**2 + a3(k)*sqr**3)    

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
                          u_ind = - f6*c6(k)/r6 - f8*c8(k)/r8 - f10*c10(k)/r10

                         !write(*,'(1x,4i3,A,A,1p,14e15.8)')i,is,j,js,an1(k),an2(k),beta(k),aalp(k),a1(k),a2(k),a3(k), &
                         !          c6(k),c8(k),c10(k),d1(k),d6(k),d8(k),d10(k),qis,qjs
                          !write(*,'(1x,4i3,A,A,1p,14e15.8)')i,is,j,js,an(is),an(js),qis,qjs
                          
                          qvdwe = qvdwe + u_exp+u_ele+u_ind 

                          ! drivative of exchange repulsion 
                          du_exp =-beta(k)*u_exp+eks*(c1(k)+2.d0*c2(k)*sqr+3.d0*c3(k)*sqr**2)
                          !du_exp=-beta(k)*u_exp+eks*expalp(k)*(a1(k)+2.d0*a2(k)*sqr+3.d0*a3(k)*sqr**2)

                          ! derivative of electrostatic part
                          df1=dtt(1, d1(k),sqr)
                          du_ele=-r4pie0*qis*qjs*((f1-1.d0)/sqr**2 - df1/sqr)
                          
                          ! derivative of asymptotic part of the dispersion interaction
                          df6=dtt(6, d6(k), sqr)
                          df8=dtt(8, d8(k), sqr)
                          df10=dtt(10, d10(k), sqr) 
                          du_ind= - df6*c6(k)/r6 + 6.d0*f6*c6(k)/sqr**7 &
                                  - df8*c8(k)/r8 + 8.d0*f8*c8(k)/sqr**9 & 
                                  - df10*c10(k)/r10 + 10.d0*f10*c10(k)/sqr**11 

                          wij = du_exp+du_ele+du_ind

!                          qvirial = qvirial - wij

                          vfij = -wij/sqr

!              Accumalating focces

                          fxsij = rsx*vfij
                          fysij = rsy*vfij
                          fzsij = rsz*vfij
      
                          fxs(is,i)=fxs(is,i)+ fxsij
                          fys(is,i)=fys(is,i)+ fysij
                          fzs(is,i)=fzs(is,i)+ fzsij

                          fxs(js,j)=fxs(js,j)- fxsij
                          fys(js,j)=fys(js,j)- fysij
                          fzs(js,j)=fzs(js,j)- fzsij

                  ! virial = - rab * fab

                           qvirial = qvirial - (  rsx*fxsij &
                                                + rsy*fysij &
                                                + rsz*fzsij)

!                            qvirial = qvirial + wij*sqr

                        end if
                     end do

!                  end if    ! rc sites
 
                end do 

             end do
              
          end if      ! rc com      
           
        end do
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
!	write(dbgUnit, *) 'P#', MPIrank, '. ccpol5s calc time	', deltaTime;
	! DBG
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
     call ccpol5s_lrc(nm,box,rcut,lvdwlrc,lvirlrc)
     vdwlrc = lvdwlrc/418.4d0
     virlrc = lvirlrc/418.4d0

     qvdwe = qvdwe/418.4d0

     qvdwe = qvdwe + vdwlrc

     qvir_vdw = qvirial !+ lvirlrc

     vir = vir + qvir_vdw

     vir_test = vir_test+qvir_vdw

!     write(*,*) 'vdw: ', vir,vir_test,qvir_vdw
 
	return;
  end subroutine ccpol5s

  subroutine ccpol5s_lrc(n,bx,rc,ulrc,wlrc)

    implicit none
    integer :: n,m
    real(8) :: ulrc ! pot eng correction
    real(8) :: wlrc ! virial correct
    real(8) :: vol, dens, bx
    real(8) :: rc
    real(8) :: uedis,uele,uin6,uin8,uin10
    real(8) :: wedis,wele,win6,win8,win10

    character(len=5) :: lrcs
    real(8) :: ulrcs,wlrcs

    vol    = bx**3.d0
    dens   = real(n)/vol      

    ! ** calculates long range corrections
! dlpoly_classical manual 1.10 page 30 equation 2.105
    ulrcs=0.d0
    wlrcs=0.d0

    ulrc=0.d0
    wlrc=0.d0

    write(5000,*)'a ',' b ','     ulrc(kcal/mol)        ','     plrc(katm)', &
              '                ulrc_nc            ', '    plrc_nc' 
    do m=1,npi

      if(abs(beta(m)) > 0.d0) then 
!      if(beta(m) .ne. 0.d0) then
      uedis = (120*c3(m)+24*beta(m)*(c2(m)+5*c3(m)*rc)+6*beta(m)**2*(c1(m)+2*rc*(2*c2(m) &
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
       ulrc = ulrc + uedis + uele + uin6 + uin8 + uin10 
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
      wele= r4pie0*(qa(m)*qb(m)*(9+d1(m)*rc*(9+d1(m)*rc*(4+d1(m)*rc))))/(d1(m)**2*exp(d1(m)*rc))
   
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

       wlrc =  wlrc + wedis+wele+win6+win8+win10


       write(5000,*) an1(m),an2(m),&
                     (2.d0*pi*n**2/vol)*(uedis+uele+uin6+uin8+uin10)/418.4d0, &
                     (0.163882576*(-(2.d0*pi*n**2/vol)* &
                     (wedis+wele+win6+win8+win10))/(3.d0*vol))/418.4d0, &
                     (uedis+uele+uin6+uin8+uin10)/418.4d0,&
                     (wedis+wele+win6+win8+win10)/418.4d0  

    end do       

    ulrc = (2.d0*pi*n**2/vol) * ulrc
    wlrc = 0.163882576*(-(2.d0*pi*n**2/vol) * wlrc)/(3.d0*vol)

  end subroutine ccpol5s_lrc 

end module vdw_ccpol5s 

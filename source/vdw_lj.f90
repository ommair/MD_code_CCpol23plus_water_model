module vdw_lj

  use, intrinsic :: iso_fortran_env
  use common_arrays
  use ewald
  
  implicit none
  private

  public :: lj, lj_lrc

contains

  subroutine lj(xx,yy,zz,xxs,yys,zzs,boxi,qvdwe)

    implicit none

    real(8) :: xx(nom),yy(nom),zz(nom)
    real(8) :: xxs(nsite,nom),yys(nsite,nom),zzs(nsite,nom)
    real(8) :: boxi

    integer :: i,j,is,js,k
    real(8) :: qvdwe, qvirial, qvir_vdw
    real(8) :: wij,qis,qjs
    real(8) :: xcc,ycc,zcc
    real(8) :: dx,dy,dz
    real(8) :: rccsq
    real(8) :: rssq, rsx,rsy,rsz
    real(8) :: sr2,sr6,sr12,vfij
    real(8) :: fxsij,fysij,fzsij
    integer :: ibox
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
!    vir_test = 0.d0 

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

          if (rccsq.lt.rcutsd) then   ! com cut off starts

             ! centre-centre rdf histogram

            ! ibox=nint(conrdf*sqrt(rccsq))+0.5
            ! write(*,*) 'ibox ',ibox 
            ! ibox=min0(ibox,maxint)
            ! krdf(ibox,1)=krdf(ibox,1)+1

             do is=1,ns

                do js=1,ns
                 
                  rsx = xcc + xxs(is,i)-xxs(js,j)
                  rsy = ycc + yys(is,i)-yys(js,j)
                  rsz = zcc + zzs(is,i)-zzs(js,j) 
                  
                  qis = chgs(is,i)
                  qjs = chgs(js,j)

                  rssq = rsx**2 + rsy**2 + rsz**2

!                  if (rssq.le.rcutsd) then       ! site cut off starts 

                  !!*     site-site rdf histogram   
                                        
                  !ibox=nint(conrdf*sqrt(rssq)+0.5)
                  !ibox=min0(ibox,maxint)
                  !krdf(ibox,3)=krdf(ibox,3)+1                                           

                     do k=1,npi
                        if( (an(is) .eq. an1(k) .and. an(js) .eq. an2(k)) .or. &
                          (an(js) .eq. an1(k) .and. an(is) .eq. an2(k)) )then

               
!                 Lennard-Jones Energy, Virial and Forces      
             
                           sr2 = sig(k)**2 / rssq
                           sr6=sr2**3

                           qvdwe = qvdwe + 4.d0*eps(k)*sr6*(sr6 - 1.d0)
        
                           wij = 24.d0*eps(k)*sr6*(2.d0*sr6 - 1.d0)

!                           qvirial = qvirial - wij  ! virial vdw
                           
                           vfij = wij/rssq

!                 Accumalating focces

                           fxsij = rsx*vfij
                           fysij = rsy*vfij
                           fzsij = rsz*vfij
      
                           fxs(is,i)=fxs(is,i)+fxsij
                           fys(is,i)=fys(is,i)+fysij
                           fzs(is,i)=fzs(is,i)+fzsij

                           fxs(js,j)=fxs(js,j)-fxsij
                           fys(js,j)=fys(js,j)-fysij
                           fzs(js,j)=fzs(js,j)-fzsij

                  ! virial = - rab * fab

                           qvirial = qvirial - (  rsx*fxsij &
                                                + rsy*fysij &
                                                + rsz*fzsij) 

                        end if       
                     end do


!                  end if    ! site cut off ends
 
                end do 

             end do
              
          end if         ! com cut off ends
           
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
!	write(dbgUnit, *) 'P#', MPIrank, '. lj calc time	', deltaTime;
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
     call lj_lrc(nm,box,rcut,lvdwlrc,lvirlrc)
     vdwlrc = lvdwlrc/418.4d0
     virlrc = lvirlrc*prsunt/(-3.d0*volm) 

     qvdwe = qvdwe/418.4d0 
     qvdwe = qvdwe + vdwlrc

!     write(*,*) qvirial

     qvir_vdw = qvirial + lvirlrc 

     vir = vir + qvir_vdw

!     vir_test = vir_test+qvir_vdw

!     write(*,*) 'vdw: ',qvir_vdw/418.4d0 !vir/418.4d0,vir_test/418.4d0,qvir_vdw/418.4d0

	return;
  end subroutine lj

  subroutine lj_lrc(n,bx,rc,ulrc,wlrc)

    implicit none
    integer :: n,m
    real(8) :: ulrc ! pot eng correction
    real(8) :: wlrc ! virial correct
    real(8) :: vol, dens, bx
    real(8) :: rc

    vol    = bx**3.d0
    dens   = real(n)/vol      

    ! ** calculates long range corrections ** for LJ POTENTIAL
! dlpoly_classical manual 1.10 page 30 equation 2.105

    ulrc=0.d0
    wlrc=0.d0
    do m=1,npi
       ulrc = ulrc+ (8.d0/9.d0)*pi*dens*real(n)*eps(m)*(sig(m)**12/rc**9) - &
           (8.d0/3.d0)*pi*dens*real(n)*eps(m)*(sig(m)**6/rc**3)

   ! conversion factor for pressure from internal units to katm is 0.163882576
! 0.163882576 * WLRC

! following is LRC for pressure, definition is taken fron dlpoly_classical
! manual 1.10 page 30 equation 2.106 divide by 3*volume to get pressure
! corrections for LJ potential

!       wlrc =  wlrc + (32.d0/9.d0)*pi*dens**2*eps(m)*(sig(m)**12/rc**9) - &
!            (16.d0/3.d0)*pi*dens**2*eps(m)*(sig(m)**6/rc**3)

       wlrc =  wlrc - ((32.d0/3.d0)*pi*dens*real(n)*eps(m)*(sig(m)**12/rc**9) - &
            (16.d0/1.d0)*pi*dens*real(n)*eps(m)*(sig(m)**6/rc**3))

    end do       

!    write(*,*) 'virlrc, pressure ', wlrc,prsunt*wlrc/(-3.d0*vol)

  end subroutine lj_lrc 

end module vdw_lj 

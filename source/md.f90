program md

! use, intrinsic :: omp_lib

 use common_arrays
 use utility
 use initialization
 use molecular_sites
 use vdw_lj
 use vdw_sapt5s 
 use vdw_ccpol5s 
 use vdw_ccpol8s
 use vdw_ccpol8s_omo
 use vdw_ccpol23plus 
 use evolve
 use statistics 
 

  implicit none

  integer :: i,j,is,js,k
  real(4) iniTime, deltaTime, tStart, tFinish
  real(4) iniTimeTotal, deltaTimeTotal, tStartTotal, tFinishTotal
  real(8) :: sumvir_test

! MPI
  CALL MPI_Init(MPIierr);
  CALL MPI_Comm_size(MPI_COMM_WORLD, MPIcomm_size, MPIierr);
  CALL MPI_Comm_rank(MPI_COMM_WORLD, MPIrank, MPIierr);
  CALL MPI_BARRIER(MPI_COMM_WORLD, MPIierr);
! MPI

! Scalar
	OMP_NumOfThreads = 1;
	id = 0;
! Scalar
! OpenMP
!	call OMP_SET_NUM_THREADS(OMP_NumOfThreads);
!	!$omp parallel DEFAULT(SHARED) private(id)
!	id = OMP_GET_THREAD_NUM();
! OpenMP
  write(*, '(a,I6,a,I6,a,I6,a,I6)') 'Process #', MPIrank, ' of', MPIcomm_size, '; thread #', id, ' of', OMP_NumOfThreads;
  !$omp barrier
  !$omp end parallel
  ! Useful environment variable: OMP_STACKSIZE
  if (MPIrank .EQ. 0) then
	  iniTimeTotal = SECNDS(0.0)
	  call CPU_TIME(tStartTotal)
  endif

  tr = 0.d0  ! restarting time in ps
  stepr = 0  ! restarting time step
  laststp = 0 ! last time step of previous simulation
  ! GRU
  vir_vdw = 0.0D0
  vir_rc = 0.0D0
  vir_kc = 0.0D0
  vir_self = 0.0D0
  ! GRU

  if (MPIrank .EQ. 0) then
	  open(100, file='OUTPUT')

	  write(100,*)
	  write(100,*) 'runninng water simulation'
	  write(100,*)


	  iniTime = SECNDS(0.0)
	  call CPU_TIME(tStart)
  endif

  call start

  if (MPIrank .EQ. 0) then
	  deltaTime = SECNDS(iniTime)
	  call CPU_TIME(tFinish)
	  write(*, *) 'Init:', deltaTime
	  write(*, *) 'Init CPU:', tFinish - tStart
	  write(401, *) 'Init:', deltaTime
	  write(401, *) 'Init CPU:', tFinish - tStart

	  write(100,*)
	  write(100,'(1x,A,e16.6,A)') 'vdw LRC: ',vdwlrc,' kcal/mol'
	  write(100,'(1x,A,e16.6,A)') 'Pres LRC: ',virlrc,' katm' !virlrc*prsunt
  endif

  t = tr          ! update restarting time in ps
  step = stepr    ! update restarting time step
  laststp = stepr ! save last time step of previous simulation

  if (MPIrank .EQ. 0) then
	  write(100,*)  
	  write(100,*) '------------------------------------------------------------------------------------------&
					------------------------------------------------------------------------------------------------'
	  write(100,*)"      STEP", "     TIME(ps)", "       SYS_TEMP(K)", "     KE(kcal/mol)", &
                      "   VDW(kcal/mol)","  COUL(kcal/mol)","  IND(kcal/mol)", "    3B(kcal/mol)", & 
                      "  TOTPE(kcal/mol)", "  TOTE(kcal/mol)", "   VIR(kcal/mol)", "   Pres(katm)", "   Volm(A^3)"
	  write(100,*)'------------------------------------------------------------------------------------------&
					------------------------------------------------------------------------------------------------'
   endif

  stpvdwpe = vdwpe
  stpcpe = scpe+rcpe+kcpe
  stpindpe = indpe
  stp3bpe = fs2pe+fs3pe
  stppe = stpvdwpe+stpcpe+stpindpe+stp3bpe
  stptotke = stptke+stprke
  stpte = stptotke + stppe

!  sumvir_test = 0.d0
!  do i=1,nm
!  do js=1,ns
!  vir = vir + (fxs(js,i)*xs(js,i)+fys(js,i)*ys(js,i)+fzs(js,i)*zs(js,i))
!  end do 
!  end do 

!  stpvir = vir/418.4d0

!   pres = (2.d0*tke - vir)/(3.d0*volm) !+ virlrc
!  pres = (2.d0*stptke - stpvir)/(3.d0*volm) !+ virlrc
!  stpprs = prsunt*pres

!  stpprs = prsunt*(2.d0*tke - vir)/(3.d0*volm)

!  write(1000,*) vir_vdw/418.4d0-3.d0*volm*virlrc/418.4d0,vir_rc+vir_kc+vir_self
    if (MPIrank .EQ. 0) write(100,'(1x,i8,1p,13e16.8)') step,t,stptmp,stptotke,stpvdwpe,stpcpe,stpindpe,&
                                                        stp3bpe,stppe,stpte,stpvir,stpprs,stpvolm

!  write(*,*) step,stptke,stprke,stptotke,stpttp,stprtp,stptmp

   stepr=stepr+1

! MPI
  MPIchunk_size = INT(nm / MPIcomm_size);
! MPI

!**************************************
! Header of History file 
  if (MPIrank .EQ. 0) then
	  open(500, file='HISTORY.xyz')
	  write(500,*) 'Created by Ommair'
	  write(500,'(3i10)') 0, 1, nm*3 
  endif
!**************************************
  if (MPIrank .EQ. 0) then
	  iniTime = SECNDS(0.0)
	  call CPU_TIME(tStart)
  endif

  do step=stepr,nsteps  ! start md time loop

    t = t+dt

    call motion

    if (step <= nequil) then
       if (ensemble.eq.'nve') then
!        call equilibration(nm,temp0,stpttp,vx,vy,vz,stprtp,jx,jy,jz) ! GRU DBG: vy->vz
       call equilibrationGRU 
       end if
    end if

    stpvdwpe = vdwpe
    stpcpe = scpe+rcpe+kcpe
    stpindpe = indpe
    stp3bpe = fs2pe+fs3pe
    stppe = stpvdwpe+stpcpe+stpindpe+stp3bpe
    stptotke = stptke+stprke
    stpte = stptotke+stppe

!    stpvir = vir_vdw+vir_rc+vir_kc+vir_self
!    stpprs = prsunt*((2.d0/3.d0)*stptke - (1.d0/3.d0)*stpvir)/volm

!     stpvir = vir/418.4d0
!     stpprs = prsunt*pres/418.4d0

!    write(*,*) step,stptke,stprke,stptotke,stpttp,stprtp,stptmp

    if (MPIrank .EQ. 0) write(100,'(1x,i8,1p,13e16.8)')step,t,stptmp,stptotke,stpvdwpe,stpcpe,stpindpe,&
                                                       stp3bpe,stppe,stpte,stpvir,stpprs,stpvolm
       

    if (step.gt.nequil) then
       call tot_sum   

!*****************************************************
! Creating HISTORY FILE
       if ((mod(step,hist_freq).eq.0) .AND. (MPIrank .EQ. 0)) then
          write(500,'(a8,4i10,f12.6)') 'timestep',step,nm*3,0,1,dt
          write(500,'(3g12.4)') box,0.d0,0.d0 
          write(500,'(3g12.4)') 0.d0,box,0.d0
          write(500,'(3g12.4)') 0.d0,0.d0,box
          k=0
          do j=1,nm
             do js=1,3 !ns
                k=k+1
                write(500,'(a8,i10,2f12.6)')an(js),k,massites(js),charge(js) 
                write(500,'(3e12.4)') x(j)+xs(js,j),y(j)+ys(js,j),z(j)+zs(js,j)
             end do
          end do
!********************************************************
       end if

    end if

! save info for restart simulation
	if (MPIrank .EQ. 0) then
		open(600, file='REVCON')
		write(600,*) 'Created by Ommair'
		write(600,*) step,t,box
		do j=1,nm
		   write(600,*)x(j),y(j),z(j),vx(j),vy(j),vz(j),jx(j),jy(j),jz(j),q0(j),q1(j),q2(j),q3(j)
	!       write(600,'(1p6E19.9E3)')x(j),y(j),z(j),vx(j),vy(j),vz(j),jx(j),jy(j),jz(j),q0(j),q1(j),q2(j),q3(j)
		end do
		close(600) 
	endif

  end do    ! ends time loop

  call averages_fluctuations

  if (MPIrank .EQ. 0) then
	deltaTime = SECNDS(iniTime)
	call CPU_TIME(tFinish)
	write(401, *) 'Main calc', deltaTime
	write(401, *) 'Main calc CPU:', tFinish - tStart
	write(*, *) 'Main calc:', deltaTime
	write(*, *) 'Main calc CPU:', tFinish - tStart


	write(100,*)'------------------------------------------------------------------------------------------&
				------------------------------------------------------------------------------------------------'
	write(100,'(24x,13A)')"       SYS_TEMP(K) ", "    KE(kcal/mol)","   VDW(kcal/mol) " ," COUL(kcal/mol)",&
        "  IND(kcal/mol)","    3B(kcal/mol)","  TOTPE(kcal/mol) ", " TOTE(kcal/mol)", "   VIR(kcal/mol)", &
                             "   Pres(katm)",  "   Volm(A^3)"
	write(100,*)'------------------------------------------------------------------------------------------&
				------------------------------------------------------------------------------------------------'

	write(100,'(A,1x,i6,1x,A,1p,13e16.8)')'Average over',navge,'steps',avtmp,avtotke,avvdwpe,avcpe,avindpe,av3bpe,&
                                             avpe,avte,avvir,avprs,avvolm  
	write(100,'(A,9x,1p,13e16.8)') 'rms fluctuations',flctmp,flctotke,flcvdwpe,flccpe,flcindpe,flc3bpe,flcpe,flcte,& 
                                                     flcvir,flcprs,flcvolm
	close(500)
	close(100)

	deltaTimeTotal = SECNDS(iniTimeTotal)
	call CPU_TIME(tFinishTotal)
	write(*, *) 'Total time:', deltaTimeTotal
	write(*, *) 'Total time CPU:', tFinishTotal - tStartTotal
	write(401, *) 'Total time:', deltaTimeTotal
	write(401, *) 'Total time CPU:', tFinishTotal - tStartTotal
  endif

! MPI
  CALL MPI_BARRIER(MPI_COMM_WORLD, MPIierr);
  CALL MPI_Finalize(MPIierr);
! MPI

end program md

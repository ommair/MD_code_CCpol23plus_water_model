module initialization

  use, intrinsic :: iso_fortran_env
  use common_arrays
  use utility
  use molecular_sites
  use statistics
  use vdw_lj
  use vdw_sapt5s 
  use vdw_ccpol5s
  use vdw_ccpol8s
  use vdw_ccpol8s_omo
  use vdw_ccpol23plus
  use ewald
  use nb_induction_model
  use threebody_potential 

  implicit none  

  private
  public :: start, control_parameters, molecule, shift_to_com
  public :: principle_moment_inertia,fcc_positions, initial_vel
  public :: ini_quaternions,random_positions
contains

  subroutine start
       
    implicit none
  
    integer :: i,j,is,js,k,m,ics,is3,nlin0,ii3
	! DBG
	real(4) iniTime, deltaTime;
	real(4) iniTimeTot, deltaTimeTot;
	real(4) iniTimePart, deltaTimePart;
	! DBG

	! DBG
	iniTimeTot = SECNDS(0.0);
	! DBG

	! DBG
	iniTimePart = SECNDS(0.0);
	! DBG
    call control_parameters     ! simulation parameters
    call molecule               ! reads molecular geometery(site's masses,charges,coords)
    call shift_to_com           ! shifiting molecule to its com
    call principle_moment_inertia  ! calculaing principle moment of inertia       
	! DBG
	deltaTimePart = SECNDS(iniTimePart);
!	write(dbgUnit, *) 'P#', MPIrank, '. start part1 time	', deltaTimePart;
	! DBG

	! DBG
	iniTimePart = SECNDS(0.0);
	! DBG
    if (restart) then

       ! read com v, J, rcom and rsites
       open(600, file='REVCON_eq')
       read(600,*) 
       read(600,*) stepr,tr,box
       do j=1,nm

          read(600,*)x(j),y(j),z(j),vx(j),vy(j),vz(j),jx(j),jy(j),jz(j),q0(j),q1(j),q2(j),q3(j)

       end do
       close(600)

      ! creating images
       do i=1,nm
          x(i) = x(i) - box*nint(x(i)/box)
          y(i) = y(i) - box*nint(y(i)/box)
          z(i) = z(i) - box*nint(z(i)/box)
       end do 

!      call sites

      call tkinetic_eng (nm, totm, vx, vy, vz, tke)
      call trans_temp (nm, tke, ttemp)
      stptke = tke/418.4d0
      stpttp = ttemp 

      call rkinetic_eng (nm,jx,jy,jz,rke)
      call rot_temp (nm, rke, rtemp)
      stprke = rke/418.4d0
      stprtp = rtemp

      call system_temp(nm,tke,rke,stptmp)

      call sites

!      write(*,*) stptke,stprke,stptke+stprke,stpttp,stprtp,stptmp

    else

      if (initial_config .eq. 'fcc_config') then
         call fcc_positions
      else if (initial_config .eq. 'ran_config') then
         call random_positions
      end if

      call ini_quaternions ! assigns initial quaternions to com

      call initial_vel     ! assigns initial linear and angular velocities to com

      call sites           ! creates sites for each com using quaternions  

      call system_mass
      call system_com
      call system_vcom

!      write(*,*) nm*totm,sysmass,xcom,ycom,zcom,vxcom,vycom,vzcom 

    end if
	! DBG
	deltaTimePart = SECNDS(iniTimePart);
!	write(dbgUnit, *) 'P#', MPIrank, '. start part2 time	', deltaTimePart;
	! DBG

    call zero            ! zero all arrays

	! DBG
	iniTimePart = SECNDS(0.0);
	! DBG
  if (MPIrank .EQ. 0) then
    write(100,*) 
    write(100,*) 'contorl parameters '
    write(100,*)


    write(100,*) 'potenial name ', pot_name
    write(100,*) 'Thermodynamical ensemble ', ensemble
    write(100,*) 'relaxation time taut and taup ', taut, taup
    write(100,*) 'random velocities and ang. momentum ', rand_vcom
    write(100,*) 'random quaternions ', rand_quat
    write(100,*)'desired temperature: ',temp0,' K'
    write(100,*)'desired pressure: ',pres0*1000.d0,' atm'
    write(100,*)'number of steps: ',nsteps
    write(100,*)'step size: ',dt
    write(100,*)'number of equilibration steps: ',nequil
  
    if (induction) then
       write(100,*)'induction model: ', induction_type, ind_model, method
    end if

    if (nonadditive_3B) then
       write(100,*)'Three Body Potential included'
    end if 
  endif
    
    volm = box**3
    rcutsd = rcut**2
    rcutsd3 = rcut3**2

    pres0 = pres0/prsunt

!    write(*,*) 'pres0',pres0

    volmt = volm
    boxt = box 
  
  if (MPIrank .EQ. 0) then
    write(100,*) 'length of box in each directon: ', box, ' angstrom' 
    write(100,*) 'volume of the simulation box: ', volm, ' angstrom**3'
    write(100,*) 'cutoff radius: ', rcut, ' angstrom'
    write(100,*) 'cutoff radius for 3B interactions: ', rcut3, ' angstrom'
    write(100,*) 'molecular mass: ', totm, ' amu'
 
    write(100,*)
    write(100,*)'****molecular diameter****'
  endif
    diameter = 2.d0*sqrt(maxval(sum(pcoord**2,dim=1)))

!    do is=1,ns
!    write(*,*) sqrt(pcoord(1,is)**2+pcoord(2,is)**2+pcoord(3,is)**2)  
!    end do

  if (MPIrank .EQ. 0) then
    write(100,*) diameter   

    write(100,*)
    write(100,*) 'field file '
    write(100,*)
    write(100,*) system_name
    write(100,*) 'No of Molecules in the Box', nm
    write(100,*) 'No of Sites per Molecule', ns
    write(100,*) 'No of pol. centers per Molecule', npol
    do is=1, ns
       write(100,*)an(is),(coords(j,is),j=1,3),massites(is),charge(is),polarizability(is)
    end do

    if (induction) then
       if (ind_model .eq. 'physical_dipole') then
          write(100,*) '    spring const       ', '        core-site # ',  'shell-site # '
          do ics=1,npol
             write(100,*) spring_k(ics)/418.4d0,icore(ics),ishell(ics)
          end do
       end if
    end if
       
    write(100,*)
    write(100,*)'no of vdw paramters',npi
    if (pot_name.eq.'tip4p' ) then
        do k=1,npi
           write(100,*) an1(k),an2(k),eps(k)/418.4d0,sig(k)
        end do
    else if (pot_name.eq.'tip5p' ) then
        do k=1,npi
           write(100,*) an1(k),an2(k),eps(k)/418.4d0,sig(k)
        end do
    else if (pot_name.eq.'sapt5s' ) then
        do k=1,npi
           write(100,'(1x,A,A,1p,14e15.8)') an1(k),an2(k),beta(k),c0(k)/418.4d0, &
                     c1(k)/418.4d0,c2(k)/418.4d0,c3(k)/418.4d0, &
                     c6(k)/418.4d0,c8(k)/418.4d0,c10(k)/418.4d0, &
                     d1(k),d6(k),d8(k),d10(k),qa(k),qb(k)
       end do    
    else if (pot_name.eq.'ccpol5s' ) then   
       do k=1,npi
          write(100,'(1x,A,A,1p,14e15.8)')an1(k),an2(k),beta(k),expalp(k), &
                     a1(k),a2(k),a3(k), &
                     c6(k)/418.4d0,c8(k)/418.4d0,c10(k)/418.4d0, &
                     d1(k),d6(k),d8(k),d10(k),qa(k),qb(k)
       end do
    else if (pot_name.eq.'ccpol8s' ) then
       do k=1,npi
          write(100,'(1x,A,A,1p,14e15.8)') an1(k),an2(k),beta(k),c0(k)/418.4d0, &
                     c1(k)/418.4d0,c2(k)/418.4d0,c3(k)/418.4d0, &
                     c6(k)/418.4d0,c8(k)/418.4d0,c10(k)/418.4d0, &
                     d1(k),d6(k),d8(k),d10(k),qa(k),qb(k)
       end do

    else if (pot_name.eq.'ccpol8s_omo' ) then
       do k=1,npi
          write(100,'(1x,A,A,1p,12e15.8)') an1(k),an2(k),beta(k),c0(k)/418.4d0, &
                   c1(k)/418.4d0,c2(k)/418.4d0,c6(k)/418.4d0,c8(k)/418.4d0,domo(k)
       end do
   
     else if (pot_name.eq.'ccpol23plus' ) then
       do k=1,npi
          write(100,'(1x,A,A,1p,14e15.8)') an1(k),an2(k),beta(k),c0(k)/418.4d0,&
                     c1(k)/418.4d0,c2(k)/418.4d0,c3(k)/418.4d0, &
                     c6(k)/418.4d0,c8(k)/418.4d0,c10(k)/418.4d0, &
                     d1(k),d6(k),d8(k),d10(k),qa(k),qb(k)
       end do


    end if 
  endif
	! DBG
	deltaTimePart = SECNDS(iniTimePart);
!	write(dbgUnit, *) 'P#', MPIrank, '. start part3 time	', deltaTimePart;
	! DBG

!!    write(100,*)
!!    write(100,*) 'shifting sites w.r.t com of molecule '
!!    write(100,*)
!!    do is=1, ns
!!       write(100,*) (coords(j,is),j=1,3)
!!    end do

!!    write(100,*)
!!    write(100,*) 'pcoords of molecule '
!!    write(100,*)
!!    do is=1, ns
!!       write(100,*) (pcoord(j,is),j=1,3)
!!    end do


!!    write(100,*)
!!    write(100,*) 'com coords of molecules on fcc '
!!    write(100,*)
!!    do i=1, nm
!!       write(100,*) x(i), y(i), z(i)
!!    end do

	! DBG
	iniTimePart = SECNDS(0.0);
	! DBG
  if (MPIrank .EQ. 0) write(100,*)
!    write(100,*) 'rescaled initial kinetic energy: ', stptke, ' kcal/mol'
!    write(100,*) 'temperature from initial rescaled velocities: ',stpttp
!    write(100,*)
!    write(100,*) 'rescaled initial rot. kinetic energy: ', stprke, ' kcal/mol'
!    write(100,*) 'rot. temperature from initial rescaled ang. momentum: ',stprtp

!    call zero

!   seting up Ewald summation
    call ewald_setup(ewp,alpha,kmax1,kmax2,kmax3,box)

!    call sites

    ! testing forces numerical vs analytical

    if(check_forces) then

!      call forces_check(hsize)
      call pressure_check(hsize)
      stop

    end if
	! DBG
	deltaTimePart = SECNDS(iniTimePart);
!	write(dbgUnit, *) 'P#', MPIrank, '. start part4 time	', deltaTimePart;
	! DBG
    
	! DBG
	iniTime = SECNDS(0.0);
	! DBG
    if (pot_name.eq.'tip4p' ) then
        call lj(x,y,z,xs,ys,zs,box,vdwpe)
    else if (pot_name.eq.'tip5p') then 
        call lj(x,y,z,xs,ys,zs,box,vdwpe)
    else if (pot_name.eq.'sapt5s') then
        call sapt5s(x,y,z,xs,ys,zs,box,vdwpe)
    else if (pot_name.eq.'ccpol5s') then
        call ccpol5s(x,y,z,xs,ys,zs,box,vdwpe)
    else if (pot_name.eq.'ccpol8s') then
        call ccpol8s(x,y,z,xs,ys,zs,box,vdwpe)
    else if (pot_name.eq.'ccpol8s_omo') then
        call ccpol8s_omo(x,y,z,xs,ys,zs,box,vdwpe)
    else if (pot_name.eq.'ccpol23plus' ) then
       call ccpol23plus(x,y,z,xs,ys,zs,box,vdwpe)
		! DBG
		deltaTime = SECNDS(iniTime);
!		write(dbgUnit, *) 'P#', MPIrank, '. ccpol23plus time	', deltaTime;
		! DBG
    end if


	! DBG
	iniTimePart = SECNDS(0.0);
	! DBG
!  Ewald summation for charges

!   call coul_pot(x,y,z,xs,ys,zs,box,alpha,rcpe)

   call ewald_self_correction(alpha,box,scpe)
   call rspace_ewald(x,y,z,xs,ys,zs,box,alpha,rcpe)
   call kspace_ewald(x,y,z,xs,ys,zs,kmax1,kmax2,kmax3,alpha,box,kcpe)
	! DBG
	deltaTimePart = SECNDS(iniTimePart);
!	write(dbgUnit, *) 'P#', MPIrank, '. start part5 time	', deltaTimePart;
	! DBG

	! DBG
	iniTimePart = SECNDS(0.0);
	! DBG
!  calling N-Body Iterative Induction model
   if (induction) then
      if(ind_model .eq. 'point_dipole') then

        call nb_induction(x,y,z,xs,ys,zs,box,indpe) 

      end if
   end if   

! calling 3B-potenital (FS2 and FS3 terms)

  if (nonadditive_3B) then 

     call threebody_term(x,y,z,xs,ys,zs,box,fs2pe,fs3pe) 

  end if

!!! pressure calculation start
  vircom=0.d0
  do i=1,nm
  do js=1,ns
  vircom = vircom + (fxs(js,i)*xs(js,i)+fys(js,i)*ys(js,i)+fzs(js,i)*zs(js,i))
!  vir_test = vir_test + (fxs(js,i)*xs(js,i)+fys(js,i)*ys(js,i)+fzs(js,i)*zs(js,i))
  end do
  end do

  stpvir = vir/418.4d0 + vircom/418.4d0

   pres = (2.d0*tke - vir - vircom)/(3.d0*volm) !+ virlrc
!   pres = (3.d0*boltz*(2.d0*nm-1.d0)*stptmp)/volm - (vir + vircom)/(3.d0*volm)

  stpprs = prsunt*pres !/418.4d0
 
  stpvolm = volm

!  write(*,*)'pres,volm,vir,vircom,vir+vircom,tke,rke',pres,volm,vir,vircom,vir+vircom,tke,rke  
!!! pressure calculation ends

	! DBG
	deltaTimePart = SECNDS(iniTimePart);
!	write(dbgUnit, *) 'P#', MPIrank, '. start part6 time	', deltaTimePart;
	! DBG

  ! testing forces numerical vs analytical  

!  if(.true.) then
!  call forces_check(0.000001d0) 
!  end if

  
    ! writing in output file
  if (MPIrank .EQ. 0) then
    write(100,*) 'Ewald Setup'
    write(100,*) 'precision required for Ewald sum ', ewp
    write(100,*) 'alpha for Ewald sum ', alpha
    write(100,*) 'kmax1, kmax2, kmax3 ', kmax1,kmax2,kmax3
    write(100,*) 'self interaction coul energy correction ', scpe, ' kcal/mol'
    write(100,*) 'rspace coul energy ', rcpe, ' kcal/mol'
    write(100,*) 'kspace coul energy ',kcpe,' kcal/mol'
    write(100,*) 'spring energy ',spge,' kcal/mol' ! GRU: spge is not calculated at all;
    write(100,*) 'alpha for Induction Model ', alphad
    write(100,*) 'alpha for 3B interactions', alpha3, 'if rcut3b = ',rcut3
    write(100,*)
!*************************************************************
!  Creating Config file for DLPOLY
!*************************************************************
    open(400, file='config_dlpoly.xyz')
    write(400,*)'Config file of ',system_name, 'for ',nm,'molecules'
    write(400,'(3i10)')0,1,nm*ns
    write(400,'(3f20.12)') box,0.d0,0.d0
    write(400,'(3f20.12)') 0.d0,box,0.d0
    write(400,'(3f20.12)') 0.d0,0.d0,box
    m=0
    do j=1,nm
       do js=1,ns
         m=m+1
         write(400,'(A,i10)')an(js), m
         write(400,'(3f20.12)') x(j)+xs(js,j), y(j)+ys(js,j),z(j)+zs(js,j)
!         write(400,*) j,js,x(j)+xs(js,j), y(j)+ys(js,j),z(j)+zs(js,j)  

       end do
    end do
    close(400)
  endif
!***************************************************************    
	! DBG
	deltaTimeTot = SECNDS(iniTimeTot);
!	write(dbgUnit, *) 'P#', MPIrank, '. start time	', deltaTimeTot;
	call MPI_BARRIER(MPI_COMM_WORLD, MPIierr);
!	STOP;
	! DBG

  end subroutine start

  subroutine control_parameters

    implicit none
  
    open(unit=10, file='CONTROL', form='formatted')

 
    read(10,*) pot_name   ! name of potential required 
    read(10,*) ensemble, taut, taup   ! thermodynamical ensemble (nve,nvt_hov,npt_hov)
    read(10,*) initial_config ! selects initial config (fcc or random)
    read(10,*) temp0      ! desired temperature in Kelvin K
    read(10,*) pres0      ! desired pressure in katm
    read(10,*) nsteps     ! number of steps
    read(10,*) dt         ! in pico-seconds (ps)
    read(10,*) nequil     ! number of equilibration steps
    read(10,*) rcut       ! cutoff radius
    read(10,*) hist_freq  ! write in history file with frquency 'hist_freq'
    read(10,*) rand_vcom  ! random velocities
    read(10,*) rand_quat  ! random quaternions
    read(10,*) induction, ind_model, induction_type, method  ! T or F induction model, point_dipole or physical_dipole , damped or undamped, (cutoff,GSF,RB,Ewald)
    read(10,*) nonadditive_3B  ! T or F
    read(10,*) rcut3
    read(10,*) restart  ! T or F (restart simulation using previous info saved)
    read(10,*) check_forces, hsize   ! check force (T of F), size of step 

  close(10)

  end subroutine control_parameters

  subroutine molecule

    implicit none

    integer :: i,j,is,js,k,ics 
    integer :: is3,nlin0,ii3

    open(unit=20, file='FIELD',form='formatted')
    read(20,*) system_name
    read(20,*) nm  ! number of molecules
    read(20,*) ns ! number of sites in a molecule
    read(20,*) npol ! number of pol centers in a molecule
    read(20,*)

    do is=1,ns
       read(20,*)an(is),(coords(j,is),j=1,3),massites(is),charge(is),polarizability(is)
    end do

!   ************ if physical dipole is in use ****************
!    read(20,*)
!    do ics=1,npol
!       read(20,*) spring_k(ics),icore(ics),ishell(ics)
!    end do

    ! converting spring const into internal units

!    do ics=1,npol
!       spring_k(ics) = spring_k(ics)*418.4d0
!    end do
!   *********************************************************


    read(20,*)
    read(20,*) npi
    read(20,*)

    if(pot_name.eq.'tip4p') then
       do k=1,npi
          read(20,*) an1(k),an2(k),eps(k),sig(k)
       end do 

       do k=1,npi
          eps(k)=eps(k)*418.4d0 ! converting into internal units of energy
       end do
    else if(pot_name.eq.'tip5p') then
       do k=1,npi
          read(20,*) an1(k),an2(k),eps(k),sig(k)
       end do

       do k=1,npi
          eps(k)=eps(k)*418.4d0 ! converting into internal units of energy
       end do
    else if (pot_name.eq.'sapt5s') then
       do k=1,npi 
          read(20,*) an1(k),an2(k),beta(k),c0(k),c1(k),c2(k),c3(k), &
               c6(k),c8(k),c10(k),d1(k),d6(k),d8(k),d10(k),qa(k),qb(k)

          !write(*,*) an1(k),an2(k),qa(k),qb(k)    
       end do

       ! sapt5s converting into internal units of energy
       do k=1,npi
          c0(k)=c0(k)*418.4d0
          c1(k)=c1(k)*418.4d0
          c2(k)=c2(k)*418.4d0
          c3(k)=c3(k)*418.4d0
          c6(k)=c6(k)*418.4d0
          c8(k)=c8(k)*418.4d0
          c10(k)=c10(k)*418.4d0
       end do
    else if (pot_name.eq.'ccpol5s') then
       do k=1,npi
          !read(20,*) an1(k),an2(k),expalp(k),beta(k),a1(k),a2(k),a3(k)
          read(20,*) an1(k),an2(k),beta(k),expalp(k),a1(k),a2(k),a3(k), &
                c6(k),c8(k),c10(k),d1(k),d6(k),d8(k),d10(k),qa(k),qb(k)
       end do
       ! ccpol5s converting into internal units of energy
       do k=1,npi    
       !   expalp(k)=expalp(k)*418.4d0
          !a1(k)=a1(k)*418.4d0
          !a2(k)=a2(k)*418.4d0
          !a3(k)=a3(k)*418.4d0
          c0(k)=expalp(k)*418.4d0
          c1(k)=a1(k)*expalp(k)*418.4d0
          c2(k)=a2(k)*expalp(k)*418.4d0
          c3(k)=a3(k)*expalp(k)*418.4d0
          c6(k)=c6(k)*418.4d0
          c8(k)=c8(k)*418.4d0
          c10(k)=c10(k)*418.4d0
       end do

    else if (pot_name.eq.'ccpol8s') then
       do k=1,npi
          !read(20,*) an1(k),an2(k),expalp(k),beta(k),a1(k),a2(k),a3(k)
          read(20,*) an1(k),an2(k),beta(k),c0(k),c1(k),c2(k),c3(k), &
                c6(k),c8(k),c10(k),d1(k),d6(k),d8(k),d10(k),qa(k),qb(k)
       end do

       ! ccpol8s converting into internal units of energy
       do k=1,npi
          c0(k)=c0(k)*418.4d0
          c1(k)=c1(k)*418.4d0
          c2(k)=c2(k)*418.4d0
          c3(k)=c3(k)*418.4d0
          c6(k)=c6(k)*418.4d0
          c8(k)=c8(k)*418.4d0
          c10(k)=c10(k)*418.4d0
       end do 

    else if (pot_name.eq.'ccpol8s_omo') then
       do k=1,npi
          read(20,*) an1(k),an2(k),beta(k),c0(k),c1(k),c2(k), &
                     c6(k),c8(k),domo(k)
       end do

       ! ccpol8s converting into internal units of energy
       do k=1,npi
          c0(k)=c0(k)*418.4d0
          c1(k)=c1(k)*418.4d0
          c2(k)=c2(k)*418.4d0
          c6(k)=c6(k)*418.4d0
          c8(k)=c8(k)*418.4d0
       end do 

    else if (pot_name.eq.'ccpol23plus') then
       do k=1,npi
          !read(20,*) an1(k),an2(k),expalp(k),beta(k),a1(k),a2(k),a3(k)
          read(20,*) an1(k),an2(k),beta(k),c0(k),c1(k),c2(k),c3(k), &
                c6(k),c8(k),c10(k),d1(k),d6(k),d8(k),d10(k),qa(k),qb(k)
       end do
       ! ccpol8s converting into internal units of energy
       do k=1,npi
          c0(k)=c0(k)*418.4d0
          c1(k)=c1(k)*418.4d0
          c2(k)=c2(k)*418.4d0
          c3(k)=c3(k)*418.4d0
          c6(k)=c6(k)*418.4d0
          c8(k)=c8(k)*418.4d0
          c10(k)=c10(k)*418.4d0
       end do   

    else
       stop
	   if (MPIrank .EQ. 0) write(100,*) 'no potential is selected'
    end if

    close(20)

! reading 3body parameters starts 

   if (nonadditive_3B) then

   open(unit=30, file='data3b_full', form='formatted')

    read(30,*)

    do is=1,17
    do js=1,17
       read(30,*) bS3(is,js),gS3(is,js),r0S3(is,js),buf1(is),buf2(js)
       !write(*,*) bS3(is,js),gS3(is,js),r0S3(is,js)
    end do
    end do


!   reading non-linear parameters fs3
    read(30,*)
    read(30,*) nlin0
    do is3=1,nlin0
       read(30,*) ii3, cs3(is3)
    end do

   close(30)

!    do is=1,17
!    do js=1,17
!      read(30,*) bS3(is,js),gS3(is,js),r0S3(is,js)
!       write(*,*) bS3(is,js),gS3(is,js),r0S3(is,js),buf1(is),buf2(js)
!    end do
!    end do

  end if

! reading 3body parameters ends  

    totm = 0.d0     ! initaiting total mass of the molecule
    natom = 0      ! initiating number of massive centers

    do i=1,ns
       if (massites(i).ne.0.d0) then
         natom = natom + 1            ! total number of massive sites
         totm  = totm + massites(i)    ! total mass of a molecule
       end if
    end do

!    write(*,*) totm

  end subroutine molecule

  subroutine shift_to_com

    implicit none

    integer :: i,j,is,js
    real(8) :: rx, ry, rz  ! com of coords

    rx = 0.d0
    ry = 0.d0
    rz = 0.d0

    do j=1, ns  ! loop on sites 
          rx = rx + massites(j)*coords(1,j)/totm
          ry = ry + massites(j)*coords(2,j)/totm
          rz = rz + massites(j)*coords(3,j)/totm
    end do

!    write(100,*) '*** COM ***'
!    write(100,*) rx, ry, rz    

    ! sites coordinates relative to com coords

    do j=1, ns
       coords(1,j) = coords(1,j) - rx
       coords(2,j) = coords(2,j) - ry
       coords(3,j) = coords(3,j) - rz
    end do

  end subroutine shift_to_com

  subroutine principle_moment_inertia

    implicit none

    integer :: i,j,is,js
    real(8) :: rmi(3,3), rot(3,3), size
    real(8) :: rxx,ryy,rzz
    real(8) :: rotxyz,aa1

    do i=1,3
       do j=1,3
          rmi(i,j) = 0.d0   ! zero all components
       end do
    end do


   do i=1,3
       do j=1,3
          do is=1, ns
             rmi(i,j) = rmi(i,j) - massites(is)*coords(i,is)*coords(j,is)
             if (i.eq.j) rmi(i,j) = rmi(i,j) + massites(is)* & 
                         (coords(1,is)**2+coords(2,is)**2+coords(3,is)**2 )
          end do
       end do
    end do
 

!   write(*,*) '*****rot. inertia tensor of molecule before jacobi*****'
!   do i=1,3
!      write(*,*) (rmi(i,j),j=1,3)
!   end do

   call jacobi(rmi,rot,3)  

   rmi(1,2)=rmi(2,1)
   rmi(1,3)=rmi(3,1)
   rmi(2,3)=rmi(3,2) 

!   write(*,*)

!   write(*,*) '*****rot. inertia tensor of molecule after jacobi*****'
!   do i=1,3
!      write(*,*) (rmi(i,j),j=1,3)
!   end do

!   write(*,*) '*****rot. mat*****'
!   do i=1,3
!      write(*,*) (rot(i,j),j=1,3)
!   end do

   size=sqrt(rmi(1,1)**2+rmi(2,2)**2+rmi(3,3)**2)/3.0
      if(rmi(1,1)/size.lt.1.0e-8.or.   &
        rmi(2,2)/size.lt.1.0e-8.or.    &
        rmi(3,3)/size.lt.1.0e-8) then 
        stop
		if (MPIrank .EQ. 0) write(100,*)'molecule off-digonal inertia error '
      end if
!  assigning principle moment of inertia

   pmi(1) = rmi(1,1)
   pmi(2) = rmi(2,2)
   pmi(3) = rmi(3,3) 
  


   do is=1,ns
      rxx = coords(1,is)
      ryy = coords(2,is)
      rzz = coords(3,is)

!  multiplication of rotation matrix with site corrdinates to transform into
!  principle corrdinates

     coords(1,is) = rxx*rot(1,1)+ryy*rot(2,1)+rzz*rot(3,1) 
     coords(2,is) = rxx*rot(1,2)+ryy*rot(2,2)+rzz*rot(3,2)
     coords(3,is) = rxx*rot(1,3)+ryy*rot(2,3)+rzz*rot(3,3)

   end do 

 

!   writee(*,*) pmi(1),pmi(2),pmi(3) 

   do is=1,ns

     pcoord(1,is) = coords(1,is)
     pcoord(2,is) = coords(2,is)
     pcoord(3,is) = coords(3,is) 

   end do  

  
  end subroutine principle_moment_inertia
  
  subroutine fcc_positions

    implicit none

    integer :: i,j,is,js
    integer :: ix,iy,iz,aa
    integer :: nc ! ! number of fcc unit cells in each coordinate direction; n=4*nc**3
    real(8) :: ndens ! number density =  na* mass_density / molarmass[pera^3]
    real(8) :: cell ! ! unit cell
    real(8) :: r(3,nm)   ! components and no of particles
    real(8), parameter :: mdens = 1.d0 ! mass density       ! g/cm^3
    real(8), parameter :: na = 6.02214076e+23 ! per mol

    ! ! sets up the fcc lattice: four molecules per unit cell
    real, dimension(3,4), parameter :: r_fcc = reshape ( [ &
         & 0.25, 0.25, 0.25, &
         & 0.25, 0.75, 0.75, &
         & 0.75, 0.75, 0.25, &
         & 0.75, 0.25, 0.75 ], [3,4] ) ! positions in unit cell
    ndens = (na * mdens)/totm
    ndens = ndens* 1.e-24
    box  = ( real(nm) / ndens ) ** ( 1.d0/3.d0 )
    nc = nint(real(nm/4)**(1.0/3.0) )
!    write(*,*) nm,4*nc**3    
    if (nm .ne. 4*nc**3) then
        write(*,*)' invalid number of molecules for fcc lattice '
        stop
    end if
    cell = box / real(nc) ! unit cell
    hbox = box / 2.0      ! half box length
    i = 0
    do iz = 0, nc-1  ! begin triple loop over unit cell indices
         do iy = 0, nc-1
            do ix = 0, nc-1
               do aa=1,4  ! begin loop over atoms in unit cell
                  i=i+1
                  r(1:3,i) = r_fcc(:,aa) +real([ix,iy,iz])  ! in range 0..real(nc)
                  r(1:3,i) = r(1:3,i) * cell            ! in range 0..real(nc)
               end do   ! end loop over atoms in unit cell
            end do
         end do
    end do    ! end triple loop over unit cell indices
   
    do i=1,nm
         x(i) = r(1,i)
         y(i) = r(2,i)
         z(i) = r(3,i)
    end do
! creating images
  do i=1,nm
         x(i) = x(i) - box*nint(x(i)/box)
         y(i) = y(i) - box*nint(y(i)/box)
         z(i) = z(i) - box*nint(z(i)/box)
   end do

  end subroutine fcc_positions

  subroutine random_positions

    implicit none

    integer :: i,j,is,js
    real(8) :: ndens ! number density =  na* mass_density / molarmass[pera^3]
    real(8) :: r(3,nm)   ! components and no of particles
    real(8), parameter :: mdens = 1.d0 ! mass density       ! g/cm^3
    real(8), parameter :: na = 6.02214076e+23 ! per mol

    ! ! sets up the random positions
    ndens = (na * mdens)/totm
    ndens = ndens* 1.e-24
    box  = ( real(nm) / ndens ) ** ( 1.d0/3.d0 )

    hbox = box / 2.0      ! half box length
   
    do i=1,nm


         call random_number(r(1:3,i))

         r(1:3,i) = r(1:3,i)*box
!         r(1:3,i) = (r(1:3,i)-0.5d0)*box 

    end do  

    do i=1,nm
         x(i) = r(1,i)
         y(i) = r(2,i)
         z(i) = r(3,i)
    end do
! creating images
   do i=1,nm
         x(i) = x(i) - box*nint(x(i)/box)
         y(i) = y(i) - box*nint(y(i)/box)
         z(i) = z(i) - box*nint(z(i)/box)
   end do

  end subroutine random_positions

!  initializing velocities and angular momentum
  subroutine initial_vel

    implicit none

    integer :: i,j,is,js
    real(8) :: totmx, totmy, totmz
    real(8) :: scfac                 ! velocity scale factor
    real(8) :: initke,inirke                ! initial ke
    real(8) :: initemp,inirtemp               ! initial temperature
    real(8) :: rtemp                 ! rotational temp
    real(8) :: qjj(3),rott(1:9)
    real(8) :: jtx(nom),jty(nom),jtz(nom)


    if (rand_vcom) then
	   if (MPIrank .EQ. 0) write(100,*) 'generating random velocities'
       scfac = sqrt(boltz*temp0/totm)

       call gauss(nm,vx,vy,vz) ! randomly generated velocities

       do i=1,nm
          vx(i) = scfac*vx(i)
          vy(i) = scfac*vy(i)
          vz(i) = scfac*vz(i)
       end do
     
 

      totmx = 0.d0
      totmy = 0.d0
      totmz = 0.d0

      do j=1,nm      
        totmx = totmx + vx(j) 
        totmy = totmy + vy(j)
        totmz = totmz + vz(j) 
      end do

      do j=1,nm
        vx(j) = vx(j) - totmx/real(nm)
        vy(j) = vy(j) - totmy/real(nm)
        vz(j) = vz(j) - totmz/real(nm)
      end do  

      call gauss(nm,jtx,jty,jtz) ! randomly generated angular momentum

     do i=1,nm

      ! current rotation matrix

!      call quat_to_rot(q0(i),q1(i),q2(i),q3(i),rott)

       ! calculate angular momentum in principal frame

      jx(i)=jtx(i)!*rott(1)+jty(i)*rott(4)+jtz(i)*rott(7)
      jy(i)=jtx(i)!*rott(2)+jty(i)*rott(5)+jtz(i)*rott(8)
      jz(i)=jtx(i)!*rott(3)+jty(i)*rott(6)+jtz(i)*rott(9)

!      write(*,*) jx(i),jtx(i)
 
    end do

   end if


   stptke = 0.d0
   stpttp = 0.d0

   call tkinetic_eng (nm, totm, vx, vy, vz, initke)
   call trans_temp (nm, initke, ttemp)

!   stptke = initke/418.4d0
!   write(*,*)
!   write(*,*) 'kinetic energy before rescaling: ', stptke, ' kcal/mol'
!   write(*,*) 'temperature from before rescaling: ', ttemp, ' k'  
       
   call rescale_velocities (nm, temp0, ttemp, vx, vy, vz)
   call tkinetic_eng (nm, totm, vx, vy, vz, tke)
   call trans_temp (nm, tke, ttemp)
   stptke = tke/418.4d0
   stpttp = ttemp

!   write(*,*)
!   write(*,*) 'rescaled kinetic energy: ', stptke, ' kcal/mol'
!   write(*,*) 'temperature from initial rescaled velocities: ',stpttp

! set angular momentum to zero

!   call gauss(nm,jx,jy,jz)


   call rkinetic_eng (nm,jx,jy,jz,inirke)
   call rot_temp (nm, inirke, rtemp)
!  write(*,*)'rotational ini. ke ', inirke/418.4d0, ' kcal/mol'
!  write(*,*)'rotational ini. temp ', rtemp

   call rescale_velocities (nm, temp0, rtemp, jx, jy, jz)
   call rkinetic_eng (nm,jx,jy,jz,rke)
   call rot_temp (nm, rke, rtemp)
! write(*,*) 'rescaled rot. kinetic energy: ', rke/418.4d0, ' kcal/mol'
! write(*,*) 'temperature from initial rescaled rot. velocities: ',rtemp

   stprke = rke/418.4d0 
   stprtp = rtemp
    
   call system_temp(nm,tke,rke,stptmp)

  end subroutine initial_vel
  
  subroutine ini_quaternions
    
    implicit none
   
    integer :: i,j,is,js
!    real(8) :: theta,phi,psi
    real(8) :: conorm 
    integer :: seed
    real(8) :: qrr(3) 
   


    if (rand_quat) then
       seed = 96873683

  ! set quaternions to give random orientiation
       call random_seed(seed)

       do j=1,nm   ! loop over molecules
!       write(*,*) quat(1,j),quat(2,j),quat(3,j),quat(4,j)

!       theta = pi*rand(seed)
!       phi = 2.0*pi*rand(seed)
!       psi = 2.0*pi*rand(seed)

!       q0(j) = cos(0.5*theta) * cos(0.5*(psi+phi))
!       q1(j) = sin(0.5*theta) * sin(0.5*(psi-phi))
!       q2(j) = sin(0.5*theta) * cos(0.5*(psi-phi))
!       q3(j) = cos(0.5*theta) * sin(0.5*(psi+phi))

          call random_number0(q0(j))
          call random_number0(q1(j))
          call random_number0(q2(j))
          call random_number0(q3(j))

          conorm=sqrt(q0(j)**2+q1(j)**2+q2(j)**2+q3(j)**2)
          q0(j)=q0(j)/conorm
          q1(j)=q1(j)/conorm
          q2(j)=q2(j)/conorm
          q3(j)=q3(j)/conorm
       end do  
 
    end if
   end subroutine ini_quaternions

   subroutine forces_check(h)

     implicit none

     integer :: i,j,is,js
     real(8) :: h
     real(8) :: qfxt,qfyt,qfzt
     real(8) :: fxa(nom),fya(nom),fza(nom) ! analytical force comp
     real(8) :: nfs2h0,nfs3h0,nvdwind0
     real(8) :: nfs2hp,nfs3hp,nvdwindp
     real(8) :: nfs2hm,nfs3hm,nvdwindm
     real(8) :: ntot0,ntotp,ntotm

    open(1001, file='numerical_vs_analytical_forces')

     fxs = 0.d0
     fys = 0.d0
     fzs = 0.d0

     call ccpol23plus(x,y,z,xs,ys,zs,box,vdwpe)

!     call ewald_self_correction(alpha,box,scpe)
!     call rspace_ewald(x,y,z,xs,ys,zs,box,alpha,rcpe)
!     call kspace_ewald(x,y,z,xs,ys,zs,kmax1,kmax2,kmax3,alpha,box,kcpe)

!     call nb_induction(x,y,z,xs,ys,zs,box,indpe)
 
!     call threebody_term(x,y,z,xs,ys,zs,box,fs2pe,fs3pe)

     nvdwind0 = vdwpe+scpe+rcpe+kcpe+indpe
     nfs2h0 = fs2pe  ! energy at h=0
     nfs3h0 = fs3pe  ! energy at h=0

     ntot0 = nfs2h0+nfs3h0+nvdwind0

     do j=1,nm
      qfxt=0.0
      qfyt=0.0
      qfzt=0.0
    do js=1,ns
      qfxt=qfxt+fxs(js,j)
      qfyt=qfyt+fys(js,j)
      qfzt=qfzt+fzs(js,j)
    end do
    fxa(j)=qfxt
    fya(j)=qfyt
    fza(j)=qfzt

    if (nonadditive_3B) then

         fxa(j) = fxa(j) + fxfs2(j) + fxfs3(j)
         fya(j) = fya(j) + fyfs2(j) + fyfs3(j)
         fza(j) = fza(j) + fzfs2(j) + fzfs3(j)

      end if

    end do

    do j=1,nm
    fxa(j) = fxa(j)/418.4d0
    fya(j) = fya(j)/418.4d0
    fza(j) = fza(j)/418.4d0
    end do

    write(1001,*) 'Analytical Forces on COM of Mol'
    write(1001,*) '      Mol. no      ', '      fx        ', '             fy        ', '               fz   '
    do j=1,3
       write(1001,*) j,fxa(j),fya(j),fza(j)
    end do

     fxs = 0.d0
     fys = 0.d0
     fzs = 0.d0

     call zero
     call fcc_positions
     call sites

!     h = 0.000001d0
     x(1) = x(1) + h

     call ccpol23plus(x,y,z,xs,ys,zs,box,vdwpe)

!     call ewald_self_correction(alpha,box,scpe)
!     call rspace_ewald(x,y,z,xs,ys,zs,box,alpha,rcpe)
!     call kspace_ewald(x,y,z,xs,ys,zs,kmax1,kmax2,kmax3,alpha,box,kcpe)

!     call nb_induction(x,y,z,xs,ys,zs,box,indpe)

!     call threebody_term(x,y,z,xs,ys,zs,box,fs2pe,fs3pe)

     nvdwindp = vdwpe+scpe+rcpe+kcpe+indpe
     nfs2hp = fs2pe  ! energy at h= +h
     nfs3hp = fs3pe  ! energy at h= +h

     ntotp = nfs2hp+nfs3hp+nvdwindp

     fxs = 0.d0
     fys = 0.d0
     fzs = 0.d0

     call zero
     call fcc_positions   
     call sites

     x(1) = x(1) - h

     call ccpol23plus(x,y,z,xs,ys,zs,box,vdwpe)

!     call ewald_self_correction(alpha,box,scpe)
!     call rspace_ewald(x,y,z,xs,ys,zs,box,alpha,rcpe)
!     call kspace_ewald(x,y,z,xs,ys,zs,kmax1,kmax2,kmax3,alpha,box,kcpe)

!     call nb_induction(x,y,z,xs,ys,zs,box,indpe)  

!     call threebody_term(x,y,z,xs,ys,zs,box,fs2pe,fs3pe)   

     nvdwindm = vdwpe+scpe+rcpe+kcpe+indpe
     nfs2hm = fs2pe  ! energy at h= -h
     nfs3hm = fs3pe  ! energy at h= -h

     ntotm = nfs2hm+nfs3hm+nvdwindm

     write(1001,*) 'Numerical Forces'
     write(1001,*) '      U(r)           ','           U(r+h)   ', '                  U(r-h)   ', & 
                   '            (U(r+h) - U(r-h))/2h           ', '% diff'
     write(1001,*) ntot0,ntotp, ntotm, -((ntotp)-(ntotm))/(2.d0*h), &
                   (-((ntotp)-(ntotm))/(2.d0*h) - fxa(1))*100/(-((ntotp)-(ntotm))/(2.d0*h)) 
        
    close(1001)

   end subroutine forces_check

   subroutine pressure_check(h)

     implicit none

     integer :: i,j,is,js
     real(8) :: h,qpres,volmi,boxi,ltke
     real(8) :: qfxt,qfyt,qfzt
     real(8) :: qvir,qvircom,lvir ! analytical virial
     real(8) :: nfs2h0,nfs3h0,nvdwind0
     real(8) :: nfs2hp,nfs3hp,nvdwindp
     real(8) :: nfs2hm,nfs3hm,nvdwindm
     real(8) :: ntot0,ntotp,ntotm

     open(1002, file='numerical_vs_analytical_pressure')
     
     qvircom = 0.d0  
     vir =0.d0
     lvir = 0.d0

     boxi=box
     volmi = boxi**3.d0

     call ccpol23plus(x,y,z,xs,ys,zs,boxi,vdwpe)
!     call ewald_self_correction(alpha,boxi,scpe)
!     call rspace_ewald(x,y,z,xs,ys,zs,boxi,alpha,rcpe)
!     call kspace_ewald(x,y,z,xs,ys,zs,kmax1,kmax2,kmax3,alpha,boxi,kcpe)
!     call nb_induction(x,y,z,xs,ys,zs,boxi,indpe)
!     call threebody_term(x,y,z,xs,ys,zs,boxi,fs2pe,fs3pe)

     nvdwind0 = vdwpe+scpe+rcpe+kcpe+indpe
     nfs2h0 = fs2pe  ! energy at h=0
     nfs3h0 = fs3pe  ! energy at h=0

     ntot0 = nfs2h0+nfs3h0+nvdwind0

     do i=1,nm
     do js=1,ns
        lvir = lvir + fxs(js,i)*xs(js,i)+fys(js,i)*ys(js,i)+fzs(js,i)*zs(js,i)
     end do
     end do

     qvircom = qvircom + lvir
     qvir = vir + qvircom
 
     ltke = tke
     call tkinetic_eng (nm, totm, vx, vy, vz, ltke)
     tke = ltke

     qpres = (2.d0*tke - qvir)/(3.d0*volmi)

     qpres = prsunt*qpres

    write(1002,*) 'Analytical pressure'
    write(1002,*)'energy','box   ', 'vol   ', 'pres ' 
    write(1002,*) ntot0,boxi,volmi,qpres

!    qvircom = 0.d0
!    vir =0.d0
!    lvir = 0.d0
    
    boxi=box+h
    volmi = boxi**3.d0

    call ccpol23plus(x,y,z,xs,ys,zs,boxi,vdwpe)
!    call ewald_self_correction(alpha,boxi,scpe)
!    call rspace_ewald(x,y,z,xs,ys,zs,boxi,alpha,rcpe)
!    call kspace_ewald(x,y,z,xs,ys,zs,kmax1,kmax2,kmax3,alpha,boxi,kcpe)
!    call nb_induction(x,y,z,xs,ys,zs,boxi,indpe)
!    call threebody_term(x,y,z,xs,ys,zs,boxi,fs2pe,fs3pe)

    nvdwindp = vdwpe+scpe+rcpe+kcpe+indpe
    nfs2hp = fs2pe  ! energy at h= +h
    nfs3hp = fs3pe  ! energy at h= +h

    ntotp = nfs2hp+nfs3hp+nvdwindp

!    qvircom = 0.d0
!    vir =0.d0
!    lvir = 0.d0

    boxi=box-h
    volmi = boxi**3.d0

    call ccpol23plus(x,y,z,xs,ys,zs,boxi,vdwpe)
!    call ewald_self_correction(alpha,boxi,scpe)
!    call rspace_ewald(x,y,z,xs,ys,zs,boxi,alpha,rcpe)
!    call kspace_ewald(x,y,z,xs,ys,zs,kmax1,kmax2,kmax3,alpha,boxi,kcpe)
!    call nb_induction(x,y,z,xs,ys,zs,boxi,indpe)
!    call threebody_term(x,y,z,xs,ys,zs,boxi,fs2pe,fs3pe)

    nvdwindm = vdwpe+scpe+rcpe+kcpe+indpe
    nfs2hm = fs2pe  ! energy at h= -h
    nfs3hm = fs3pe  ! energy at h= -h

    ntotm = nfs2hm+nfs3hm+nvdwindm

    write(1002,*) 'Numerical Pressure'
    write(1002,*) '      U(V)           ','           U(V+delV)   ', ' U(V-delV)   ', &
                   '            (U(V+delV) - U(V-delV))/2delV           ', '% diff'
    write(1002,*) ntot0,ntotp,ntotm,prsunt*(2.d0*tke-((ntotp)-(ntotm))/(2.d0*h**3.d0))/(3.d0*box**3.d0)

!    write(*,*) ((ntotp)-(ntotm))/(2.d0*h**3.d0),vir,lvir,vir+lvir 

    close(1002)
 
   end subroutine pressure_check


end module initialization

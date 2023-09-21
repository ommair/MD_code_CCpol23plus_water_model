module evolve
! OpenMP version
  
  use, intrinsic :: iso_fortran_env
  use common_arrays
  use molecular_sites
  use utility
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

  public :: motion  
  public :: nosquish
  public :: nve,nvt_hov,nvt_b,npt_hov

contains

  subroutine motion

    implicit none

    real(8) :: jjj(3)

    if (ensemble.eq.'nve') then

!        call nve(dt)
        call nve_2(dt)


    else if (ensemble.eq.'nvt_h') then
 
!        call nvt_hov(dt)
        call nvt_hov_2(dt) 

    else if (ensemble.eq.'nvt_b') then 

!        call nvt_b(dt)
        call nvt_b_2(dt)

    else if (ensemble.eq.'npt_h') then

!        call npt_hov(dt)
!        call npt_hov_2(dt)
        call npt_hov_3(dt) 

    else if (ensemble.eq.'npt_b') then 
  
        call npt_b(dt)

    else 
       stop
       write(*,*) 'No ensemble  selected'

     end if  

  end subroutine motion

  subroutine nve(ttt)
   implicit none
   real(8) :: ttt

   integer :: i,j,is,js
   real(8) :: qtke,qrke,lvir,qvir
   real(8) :: qfx,qfy,qfz
   real(8) :: qtx,qty,qtz, tqx,tqy,tqz
!   real(8) :: qvx,qvy,qvz
   real(8) :: qx,qy,qz
   real(8) :: qjj(3),omg(3)
   real(8) :: qq(0:3), rott(1:9)
   real(8) :: qt0,qt1,qt2,qt3
   real(8) :: opx,opy,opz
   real(8) :: lvdwpe,lscpe,lrcpe,lkcpe,lindpe,lfs2pe,lfs3pe
   real(8) :: ltke,lttemp,lrke,lrtemp,lstptmp
   real(8) :: qq0,qq1,qq2,qq3,pp0,pp1,pp2,pp3
!   real(8) :: q0t(nom),q1t(nom),q2t(nom),q3t(nom)
!   real(8) :: p0(nom),p1(nom),p2(nom),p3(nom)
   real(8), allocatable :: p0(:),p1(:),p2(:),p3(:)
   real(8) :: cconorm
	! DBG
	real(4) iniTime, deltaTime;
	real(4) iniTimeTot, deltaTimeTot;
	! DBG

	! DBG
	iniTimeTot = SECNDS(0.0);
	! DBG

   allocate(p0(nom))
   allocate(p1(nom))
   allocate(p2(nom))
   allocate(p3(nom))

   qtke = 0.d0
   qrke = 0.d0

!$omp parallel DEFAULT(SHARED) private(id,lvir,qq0,qq1,qq2,qq3,pp0,pp1,pp2,pp3)
   lvir = 0.0D0
!$omp do private(qfx,qfy,qfz,qtx,qty,qtz,opx,opy,opz) schedule(dynamic)
   do i=1,nm

      qfx=0.0
      qfy=0.0
      qfz=0.0
      qtx=0.0
      qty=0.0
      qtz=0.0

      do js=1,ns
         qfx=qfx+fxs(js,i)
         qfy=qfy+fys(js,i)
         qfz=qfz+fzs(js,i)
         qtx=qtx+ys(js,i)*fzs(js,i)-zs(js,i)*fys(js,i)
         qty=qty+zs(js,i)*fxs(js,i)-xs(js,i)*fzs(js,i)
         qtz=qtz+xs(js,i)*fys(js,i)-ys(js,i)*fxs(js,i)

         ! CONVERT FROM SITE-SITE TO MOLECULE-MOLECULE VIRIAL

         lvir = lvir + fxs(js,i)*xs(js,i)+fys(js,i)*ys(js,i)+fzs(js,i)*zs(js,i)   

      end do

      fx(i)=qfx
      fy(i)=qfy
      fz(i)=qfz
      tx(i)=qtx
      ty(i)=qty
      tz(i)=qtz

      if (nonadditive_3B) then

         fx(i) = fx(i) + fxfs2(i) + fxfs3(i)
         fy(i) = fy(i) + fyfs2(i) + fyfs3(i)
         fz(i) = fz(i) + fzfs2(i) + fzfs3(i)

         tx(i)= tx(i) + txfs2(i) + txfs3(i)
         ty(i)= ty(i) + tyfs2(i) + tyfs3(i)
         tz(i)= tz(i) + tzfs2(i) + tzfs3(i)

      end if 

 !  end do

   !   *************  Rigid body motion ****************************
!     operations common to first and second stages

!   do i=1,nm

      ! calculate quaternion momenta at start of time step
      opx = jx(i)
      opy = jy(i)
      opz = jz(i)

      p0(i)=2.0d0*(-q1(i)*opx-q2(i)*opy-q3(i)*opz)
      p1(i)=2.0d0*( q0(i)*opx-q3(i)*opy+q2(i)*opz)
      p2(i)=2.0d0*( q3(i)*opx+q0(i)*opy-q1(i)*opz)
      p3(i)=2.0d0*(-q2(i)*opx+q1(i)*opy+q0(i)*opz)
   end do
!$omp end do

!$omp atomic
   vircom = vircom + lvir

!$omp barrier

!$omp do private(rott,tqx,tqy,tqz,qt0,qt1,qt2,qt3) schedule(dynamic)
   do i=1,nm
      
      ! current rotation matrix

      call quat_to_rot(q0(i),q1(i),q2(i),q3(i),rott)

       ! calculate torque in principal frame

      tqx=tx(i)*rott(1)+ty(i)*rott(4)+tz(i)*rott(7)
      tqy=tx(i)*rott(2)+ty(i)*rott(5)+tz(i)*rott(8)
      tqz=tx(i)*rott(3)+ty(i)*rott(6)+tz(i)*rott(9)

       ! calculate quaternion torques

      qt0=2.0d0*(-q1(i)*tqx-q2(i)*tqy-q3(i)*tqz)
      qt1=2.0d0*( q0(i)*tqx-q3(i)*tqy+q2(i)*tqz)
      qt2=2.0d0*( q3(i)*tqx+q0(i)*tqy-q1(i)*tqz)
      qt3=2.0d0*(-q2(i)*tqx+q1(i)*tqy+q0(i)*tqz)


      ! calculate quaternion momenta at start of time step
!      opx = jx(i)
!      opy = jy(i)
!      opz = jz(i)

!      p0(i)=2.0d0*(-q1(i)*opx-q2(i)*opy-q3(i)*opz)
!      p1(i)=2.0d0*( q0(i)*opx-q3(i)*opy+q2(i)*opz)
!      p2(i)=2.0d0*( q3(i)*opx+q0(i)*opy-q1(i)*opz)
!      p3(i)=2.0d0*(-q2(i)*opx+q1(i)*opy+q0(i)*opz)

       ! update quaternion momenta by 1/2 time step

      p0(i)=p0(i)+qt0*0.5d0*ttt
      p1(i)=p1(i)+qt1*0.5d0*ttt
      p2(i)=p2(i)+qt2*0.5d0*ttt
      p3(i)=p3(i)+qt3*0.5d0*ttt

      ! update centre of mass velocity by 1/2 time step

      vx(i) = vx(i) + 0.5d0*fx(i)*ttt/totm
      vy(i) = vy(i) + 0.5d0*fy(i)*ttt/totm
      vz(i) = vz(i) + 0.5d0*fz(i)*ttt/totm

!      write(101,*) vx(i),vy(i),vz(i)

   end do
!$omp end do
    ! move centre of mass by full time step

!$omp do schedule(dynamic)
   do i=1,nm
      x(i)= x(i) + ttt * vx(i)
      y(i)= y(i) + ttt * vy(i)
      z(i)= z(i) + ttt * vz(i)

      ! periodic boundary conditions

!     x(i) = x(i) - box*nint(x(i)/box)
!     y(i) = y(i) - box*nint(y(i)/box)
!     z(i) = z(i) - box*nint(z(i)/box)
 !  end do

      ! rotate rigid groups: nosquish algorithm -update q to full tstep and
      ! amend p and get new rotation matrix

!   do i=1,nm   

   qq0 = q0(i)
   qq1 = q1(i)
   qq2 = q2(i)
   qq3 = q3(i)
   pp0 = p0(i)
   pp1 = p1(i)
   pp2 = p2(i)
   pp3 = p3(i)

   call nosquish(ttt,qq0,qq1,qq2,qq3,pp0,pp1,pp2,pp3) 

   q0(i) = qq0
   q1(i) = qq1
   q2(i) = qq2
   q3(i) = qq3
   p0(i) = pp0
   p1(i) = pp1
   p2(i) = pp2
   p3(i) = pp3
   end do
!$omp end do

!   call sites

!$omp do private(opx,opy,opz)
   do i=1,nm
!      qq(0) = q0(i)
!      qq(1) = q1(i)
!      qq(2) = q2(i)
!      qq(3) = q3(i)

!      call quat_to_rot(qq,rott)

       ! update COM velocities and angular momentum to 1/2 tstep

      opx=0.5d0*(-q1(i)*p0(i)+q0(i)*p1(i)+ &
                 q3(i)*p2(i)-q2(i)*p3(i))
      opy=0.5d0*(-q2(i)*p0(i)-q3(i)*p1(i)+ &
                 q0(i)*p2(i)+q1(i)*p3(i))
      opz=0.5d0*(-q3(i)*p0(i)+q2(i)*p1(i)- &
                 q1(i)*p2(i)+q0(i)*p3(i))

      jx(i) = opx
      jy(i) = opy
      jz(i) = opz
   end do
!$omp end do
!$omp end parallel
!      vx(i) = vx(i) + 0.5d0*fx(i)*tt/totm
!      vy(i) = vy(i) + 0.5d0*fy(i)*tt/totm
!      vz(i) = vz(i) + 0.5d0*fz(i)*tt/totm

     ! first stage of velocity verlet algorithm

     ! move centre of mass by full time step

!      x(i)= x(i) + tt * vx(i)
!      y(i)= y(i) + tt * vy(i)
!      z(i)= z(i) + tt * vz(i)

      ! periodic boundary conditions

!    do i=1,nm

!      x(i) = x(i) - box*nint(x(i)/box)
!      y(i) = y(i) - box*nint(y(i)/box)
!      z(i) = z(i) - box*nint(z(i)/box)

!   end do

! end of first stage of velocity verlet algorithm

   call sites

   lvdwpe = vdwpe

	! DBG
	iniTime = SECNDS(0.0);
	! DBG
   if (pot_name.eq.'tip4p' ) then 
       call lj(x,y,z,xs,ys,zs,box,lvdwpe)
   else if (pot_name.eq.'tip5p' ) then
       call lj(x,y,z,xs,ys,zs,box,lvdwpe)
   else if (pot_name.eq.'sapt5s' ) then
       call sapt5s(x,y,z,xs,ys,zs,box,lvdwpe)
   else if (pot_name.eq.'ccpol5s' ) then
       call ccpol5s(x,y,z,xs,ys,zs,box,lvdwpe)
   else if (pot_name.eq.'ccpol8s' ) then
       call ccpol8s(x,y,z,xs,ys,zs,box,lvdwpe)       
   else if (pot_name.eq.'ccpol8s_omo' ) then
       call ccpol8s_omo(x,y,z,xs,ys,zs,box,lvdwpe)       
   else if (pot_name.eq.'ccpol23plus' ) then
       call ccpol23plus(x,y,z,xs,ys,zs,box,lvdwpe)
		! DBG
		deltaTime = SECNDS(iniTime);
!		write(dbgUnit, *) 'P#', MPIrank, '. ccpol23plus time	', deltaTime;
		! DBG
   end if

   vdwpe = lvdwpe

   lscpe = scpe

   call ewald_self_correction(alpha,box,lscpe)

   scpe = lscpe
   lrcpe = rcpe

   call rspace_ewald(x,y,z,xs,ys,zs,box,alpha,lrcpe)
  
   rcpe = lrcpe
   lkcpe = kcpe

   call kspace_ewald(x,y,z,xs,ys,zs,kmax1,kmax2,kmax3,alpha,box,lkcpe)

   kcpe = lkcpe

!  calling N-Body Iterative Induction model

   if (induction) then
      if(ind_model .eq. 'point_dipole') then

        lindpe = indpe
        
        call nb_induction(x,y,z,xs,ys,zs,box,lindpe)

        indpe = lindpe

      end if
   end if 

! calling 3B-potenital (FS2 and FS3 terms)

  if (nonadditive_3B) then

     lfs2pe = fs2pe
     lfs3pe = fs3pe

     call threebody_term(x,y,z,xs,ys,zs,box,lfs2pe,lfs3pe)

     fs2pe = lfs2pe
     fs3pe = lfs3pe

  end if  

   ! second stage of velocity verlet algorithm

!$omp parallel DEFAULT(SHARED)&
!$omp& private(id,lvir,qq0,qq1,qq2,qq3,pp0,pp1,pp2,pp3)
   lvir = 0.0D0
!$omp do private(j,qfx,qfy,qfz,qtx,qty,qtz) schedule(dynamic)
   do j=1,nm

      qfx=0.0
      qfy=0.0
      qfz=0.0
      qtx=0.0
      qty=0.0
      qtz=0.0

      do js=1,ns
         qfx=qfx+fxs(js,j)
         qfy=qfy+fys(js,j)
         qfz=qfz+fzs(js,j)
         qtx=qtx+ys(js,j)*fzs(js,j)-zs(js,j)*fys(js,j)
         qty=qty+zs(js,j)*fxs(js,j)-xs(js,j)*fzs(js,j)
         qtz=qtz+xs(js,j)*fys(js,j)-ys(js,j)*fxs(js,j)

         ! CONVERT FROM SITE-SITE TO MOLECULE-MOLECULE VIRIAL

         lvir = lvir + fxs(js,j)*xs(js,j)+fys(js,j)*ys(js,j)+fzs(js,j)*zs(js,j)

      end do

      fx(j)=qfx
      fy(j)=qfy
      fz(j)=qfz
      tx(j)=qtx
      ty(j)=qty
      tz(j)=qtz

      if (nonadditive_3B) then

         fx(j) = fx(j) + fxfs2(j) + fxfs3(j)
         fy(j) = fy(j) + fyfs2(j) + fyfs3(j)
         fz(j) = fz(j) + fzfs2(j) + fzfs3(j)

         tx(j)= tx(j) + txfs2(j) + txfs3(j)
         ty(j)= ty(j) + tyfs2(j) + tyfs3(j)
         tz(j)= tz(j) + tzfs2(j) + tzfs3(j)

      end if 

   end do
!$omp end do
!$omp atomic

   vircom = vircom + lvir

!$omp barrier

!$omp do private(i,opx,opy,opz) schedule(dynamic)
   do i=1,nm

      ! calculate quaternion momenta at start of time step
      opx = jx(i)
      opy = jy(i)
      opz = jz(i)

      p0(i)=2.0d0*(-q1(i)*opx-q2(i)*opy-q3(i)*opz)
      p1(i)=2.0d0*( q0(i)*opx-q3(i)*opy+q2(i)*opz)
      p2(i)=2.0d0*( q3(i)*opx+q0(i)*opy-q1(i)*opz)
      p3(i)=2.0d0*(-q2(i)*opx+q1(i)*opy+q0(i)*opz)
   end do
!$omp end do

!$omp barrier

!$omp do private(i,tqx,tqy,tqz,rott,qt0,qt1,qt2,qt3) schedule(dynamic)
   do i=1,nm

       ! current rotation matrix

      call quat_to_rot(q0(i),q1(i),q2(i),q3(i),rott)

       ! calculate torque in principal frame

      tqx=tx(i)*rott(1)+ty(i)*rott(4)+tz(i)*rott(7)
      tqy=tx(i)*rott(2)+ty(i)*rott(5)+tz(i)*rott(8)
      tqz=tx(i)*rott(3)+ty(i)*rott(6)+tz(i)*rott(9)

       ! calculate quaternion torques

      qt0=2.0d0*(-q1(i)*tqx-q2(i)*tqy-q3(i)*tqz)
      qt1=2.0d0*( q0(i)*tqx-q3(i)*tqy+q2(i)*tqz)
      qt2=2.0d0*( q3(i)*tqx+q0(i)*tqy-q1(i)*tqz)
      qt3=2.0d0*(-q2(i)*tqx+q1(i)*tqy+q0(i)*tqz)

      ! calculate quaternion momenta at start of time step
!      opx = jx(i)
!      opy = jy(i)
!      opz = jz(i)

!      p0(i)=2.0d0*(-q1(i)*opx-q2(i)*opy-q3(i)*opz)
!      p1(i)=2.0d0*( q0(i)*opx-q3(i)*opy+q2(i)*opz)
!      p2(i)=2.0d0*( q3(i)*opx+q0(i)*opy-q1(i)*opz)
!      p3(i)=2.0d0*(-q2(i)*opx+q1(i)*opy+q0(i)*opz)

       ! update quaternion momenta by 1/2 time step

      p0(i)=p0(i)+qt0*0.5d0*ttt
      p1(i)=p1(i)+qt1*0.5d0*ttt
      p2(i)=p2(i)+qt2*0.5d0*ttt
      p3(i)=p3(i)+qt3*0.5d0*ttt

      !  update Rgid body angular momentum & COM velocities to full step

!      opx=0.5d0*(-q1(i)*p0(i)+q0(i)*p1(i)+ &
!                 q3(i)*p2(i)-q2(i)*p3(i))
!      opy=0.5d0*(-q2(i)*p0(i)-q3(i)*p1(i)+ &
!                 q0(i)*p2(i)+q1(i)*p3(i))
!      opz=0.5d0*(-q3(i)*p0(i)+q2(i)*p1(i)- &
!                 q1(i)*p2(i)+q0(i)*p3(i))

!      jx(i) = opx
!      jy(i) = opy
!      jz(i) = opz

      vx(i) = vx(i) + 0.5d0*fx(i)*ttt/totm
      vy(i) = vy(i) + 0.5d0*fy(i)*ttt/totm
      vz(i) = vz(i) + 0.5d0*fz(i)*ttt/totm

   end do
!$omp end do

!$omp barrier

!$omp do private(i,opx,opy,opz) schedule(dynamic)
   do i=1,nm
!      qq(0) = q0(i)
!      qq(1) = q1(i)
!      qq(2) = q2(i)
!      qq(3) = q3(i)

!      call quat_to_rot(qq,rott)

       ! update COM velocities and angular momentum to 1/2 tstep

      opx=0.5d0*(-q1(i)*p0(i)+q0(i)*p1(i)+ &
                 q3(i)*p2(i)-q2(i)*p3(i))
      opy=0.5d0*(-q2(i)*p0(i)-q3(i)*p1(i)+ &
                 q0(i)*p2(i)+q1(i)*p3(i))
      opz=0.5d0*(-q3(i)*p0(i)+q2(i)*p1(i)- &
                 q1(i)*p2(i)+q0(i)*p3(i))

      jx(i) = opx
      jy(i) = opy
      jz(i) = opz

!      call quat_to_rot(q0(i),q1(i),q2(i),q3(i),rott)

   end do
!$omp end do
!$omp end parallel
 
   ltke = tke

   call tkinetic_eng (nm, totm, vx, vy, vz, ltke)

   tke = ltke
   lttemp = ttemp

   call trans_temp (nm, tke, lttemp)

   ttemp = lttemp

   stptke = tke/418.4
   stpttp = ttemp

   lrke = rke

   call rkinetic_eng (nm, jx, jy, jz, lrke)

   rke = lrke
   stprke = rke/418.4

   lrtemp = rtemp
   call rot_temp (nm, rke, lrtemp)
   rtemp = lrtemp
   stprtp = rtemp

   lstptmp = stptmp
   call system_temp(nm,tke,rke,lstptmp)
   stptmp = lstptmp


   qvir = vir+vircom

!  calculate system pressure

   pres = (2.d0*ltke - qvir)/(3.d0*volm)

!   write(*,*) box,box**3.d0,volm,pres,vir,ltke

   stpvolm = volm

   stpvir = qvir/418.4d0
   stpprs = prsunt*pres

!   write(*,*) box,box**3.d0,stpvolm,stpprs,stpvir,stptmp

   deallocate(p0)
   deallocate(p1)
   deallocate(p2)
   deallocate(p3)
	! DBG
	deltaTimeTot = SECNDS(iniTimeTot);
!	write(dbgUnit, *) 'P#', MPIrank, '. nve time	', deltaTimeTot;
	call MPI_BARRIER(MPI_COMM_WORLD, MPIierr);
!	STOP;
	! DBG

  end subroutine nve  


  subroutine nve_2(ttt)
   implicit none
   real(8) :: ttt

   integer :: i,j,is,js
   real(8) :: qtke,qrke,lvir,qvir
   real(8) :: qfx,qfy,qfz
   real(8) :: qtx,qty,qtz, tqx,tqy,tqz
!   real(8) :: qvx,qvy,qvz
   real(8) :: qx,qy,qz
   real(8) :: qjj(3),omg(3)
   real(8) :: qq(0:3), rott(1:9)
   real(8) :: qt0,qt1,qt2,qt3
   real(8) :: opx,opy,opz
   real(8) :: lvdwpe,lscpe,lrcpe,lkcpe,lindpe,lfs2pe,lfs3pe
   real(8) :: ltke,lttemp,lrke,lrtemp,lstptmp
   real(8) :: qq0,qq1,qq2,qq3,pp0,pp1,pp2,pp3
!   real(8) :: q0t(nom),q1t(nom),q2t(nom),q3t(nom)
!   real(8) :: p0(nom),p1(nom),p2(nom),p3(nom)
   real(8), allocatable :: p0(:),p1(:),p2(:),p3(:)
   real(8) :: cconorm
        ! DBG
        real(4) iniTime, deltaTime;
        real(4) iniTimeTot, deltaTimeTot;
        ! DBG

        ! DBG
        iniTimeTot = SECNDS(0.0);
        ! DBG

   allocate(p0(nom))
   allocate(p1(nom))
   allocate(p2(nom))
   allocate(p3(nom))

   qtke = 0.d0
   qrke = 0.d0

   lvir = 0.0D0

   do i=1,nm

      qfx=0.0
      qfy=0.0
      qfz=0.0
      qtx=0.0
      qty=0.0
      qtz=0.0

      do js=1,ns
         qfx=qfx+fxs(js,i)
         qfy=qfy+fys(js,i)
         qfz=qfz+fzs(js,i)
         qtx=qtx+ys(js,i)*fzs(js,i)-zs(js,i)*fys(js,i)
         qty=qty+zs(js,i)*fxs(js,i)-xs(js,i)*fzs(js,i)
         qtz=qtz+xs(js,i)*fys(js,i)-ys(js,i)*fxs(js,i)

         ! CONVERT FROM SITE-SITE TO MOLECULE-MOLECULE VIRIAL

         lvir = lvir + fxs(js,i)*xs(js,i)+fys(js,i)*ys(js,i)+fzs(js,i)*zs(js,i)   

      end do

      fx(i)=qfx
      fy(i)=qfy
      fz(i)=qfz
      tx(i)=qtx
      ty(i)=qty
      tz(i)=qtz

      if (nonadditive_3B) then

         fx(i) = fx(i) + fxfs2(i) + fxfs3(i)
         fy(i) = fy(i) + fyfs2(i) + fyfs3(i)
         fz(i) = fz(i) + fzfs2(i) + fzfs3(i)

         tx(i)= tx(i) + txfs2(i) + txfs3(i)
         ty(i)= ty(i) + tyfs2(i) + tyfs3(i)
         tz(i)= tz(i) + tzfs2(i) + tzfs3(i)

      end if 

   end do

   !   *************  Rigid body motion ****************************
!     operations common to first and second stages

   do i=1,nm

      ! calculate quaternion momenta at start of time step
      opx = jx(i)
      opy = jy(i)
      opz = jz(i)

      p0(i)=2.0d0*(-q1(i)*opx-q2(i)*opy-q3(i)*opz)
      p1(i)=2.0d0*( q0(i)*opx-q3(i)*opy+q2(i)*opz)
      p2(i)=2.0d0*( q3(i)*opx+q0(i)*opy-q1(i)*opz)
      p3(i)=2.0d0*(-q2(i)*opx+q1(i)*opy+q0(i)*opz)
   end do

   vircom = vircom + lvir

   do i=1,nm
      
      ! current rotation matrix

      call quat_to_rot(q0(i),q1(i),q2(i),q3(i),rott)

       ! calculate torque in principal frame

      tqx=tx(i)*rott(1)+ty(i)*rott(4)+tz(i)*rott(7)
      tqy=tx(i)*rott(2)+ty(i)*rott(5)+tz(i)*rott(8)
      tqz=tx(i)*rott(3)+ty(i)*rott(6)+tz(i)*rott(9)

       ! calculate quaternion torques

      qt0=2.0d0*(-q1(i)*tqx-q2(i)*tqy-q3(i)*tqz)
      qt1=2.0d0*( q0(i)*tqx-q3(i)*tqy+q2(i)*tqz)
      qt2=2.0d0*( q3(i)*tqx+q0(i)*tqy-q1(i)*tqz)
      qt3=2.0d0*(-q2(i)*tqx+q1(i)*tqy+q0(i)*tqz)

      ! update quaternion momenta by 1/2 time step

      p0(i)=p0(i)+qt0*0.5d0*ttt
      p1(i)=p1(i)+qt1*0.5d0*ttt
      p2(i)=p2(i)+qt2*0.5d0*ttt
      p3(i)=p3(i)+qt3*0.5d0*ttt

      ! update centre of mass velocity by 1/2 time step

      vx(i) = vx(i) + 0.5d0*fx(i)*ttt/totm
      vy(i) = vy(i) + 0.5d0*fy(i)*ttt/totm
      vz(i) = vz(i) + 0.5d0*fz(i)*ttt/totm


   end do

    ! move centre of mass by full time step

   do i=1,nm

      x(i)= x(i) + ttt * vx(i)
      y(i)= y(i) + ttt * vy(i)
      z(i)= z(i) + ttt * vz(i)

   end do
 ! rotate rigid groups: nosquish algorithm -update q to full tstep and amend p and get new rotation matrix

   do i=1,nm   

      qq0 = q0(i)
      qq1 = q1(i)
      qq2 = q2(i)
      qq3 = q3(i)
      pp0 = p0(i)
      pp1 = p1(i)
      pp2 = p2(i)
      pp3 = p3(i)

      call nosquish(ttt,qq0,qq1,qq2,qq3,pp0,pp1,pp2,pp3) 

      q0(i) = qq0
      q1(i) = qq1
      q2(i) = qq2
      q3(i) = qq3
      p0(i) = pp0
      p1(i) = pp1
      p2(i) = pp2
      p3(i) = pp3

   end do

! end of first stage of velocity verlet algorithm

   call sites

   lvdwpe = vdwpe

        ! DBG
        iniTime = SECNDS(0.0);
        ! DBG
   if (pot_name.eq.'tip4p' ) then 
       call lj(x,y,z,xs,ys,zs,box,lvdwpe)
   else if (pot_name.eq.'tip5p' ) then
       call lj(x,y,z,xs,ys,zs,box,lvdwpe)
   else if (pot_name.eq.'sapt5s' ) then
       call sapt5s(x,y,z,xs,ys,zs,box,lvdwpe)
   else if (pot_name.eq.'ccpol5s' ) then
       call ccpol5s(x,y,z,xs,ys,zs,box,lvdwpe)
   else if (pot_name.eq.'ccpol8s' ) then
       call ccpol8s(x,y,z,xs,ys,zs,box,lvdwpe)       
   else if (pot_name.eq.'ccpol8s_omo' ) then
       call ccpol8s_omo(x,y,z,xs,ys,zs,box,lvdwpe)       
   else if (pot_name.eq.'ccpol23plus' ) then
       call ccpol23plus(x,y,z,xs,ys,zs,box,lvdwpe)
                ! DBG
                deltaTime = SECNDS(iniTime);
!               write(dbgUnit, *) 'P#', MPIrank, '. ccpol23plus time    ',
!               deltaTime;
                ! DBG
   end if

   vdwpe = lvdwpe

   lscpe = scpe

   call ewald_self_correction(alpha,box,lscpe)

   scpe = lscpe
   lrcpe = rcpe

   call rspace_ewald(x,y,z,xs,ys,zs,box,alpha,lrcpe)
  
   rcpe = lrcpe
   lkcpe = kcpe

   call kspace_ewald(x,y,z,xs,ys,zs,kmax1,kmax2,kmax3,alpha,box,lkcpe)

   kcpe = lkcpe

!  calling N-Body Iterative Induction model

   if (induction) then
      if(ind_model .eq. 'point_dipole') then

        lindpe = indpe
        
        call nb_induction(x,y,z,xs,ys,zs,box,lindpe)

        indpe = lindpe

      end if
   end if 

! calling 3B-potenital (FS2 and FS3 terms)

  if (nonadditive_3B) then

     lfs2pe = fs2pe
     lfs3pe = fs3pe

     call threebody_term(x,y,z,xs,ys,zs,box,lfs2pe,lfs3pe)

     fs2pe = lfs2pe
     fs3pe = lfs3pe

  end if  

   ! second stage of velocity verlet algorithm


   lvir = 0.0D0
   do j=1,nm

      qfx=0.0
      qfy=0.0
      qfz=0.0
      qtx=0.0
      qty=0.0
      qtz=0.0

      do js=1,ns
         qfx=qfx+fxs(js,j)
         qfy=qfy+fys(js,j)
         qfz=qfz+fzs(js,j)
         qtx=qtx+ys(js,j)*fzs(js,j)-zs(js,j)*fys(js,j)
         qty=qty+zs(js,j)*fxs(js,j)-xs(js,j)*fzs(js,j)
         qtz=qtz+xs(js,j)*fys(js,j)-ys(js,j)*fxs(js,j)

         ! CONVERT FROM SITE-SITE TO MOLECULE-MOLECULE VIRIAL

         lvir = lvir + fxs(js,j)*xs(js,j)+fys(js,j)*ys(js,j)+fzs(js,j)*zs(js,j)

      end do

      fx(j)=qfx
      fy(j)=qfy
      fz(j)=qfz
      tx(j)=qtx
      ty(j)=qty
      tz(j)=qtz

      if (nonadditive_3B) then

         fx(j) = fx(j) + fxfs2(j) + fxfs3(j)
         fy(j) = fy(j) + fyfs2(j) + fyfs3(j)
         fz(j) = fz(j) + fzfs2(j) + fzfs3(j)

         tx(j)= tx(j) + txfs2(j) + txfs3(j)
         ty(j)= ty(j) + tyfs2(j) + tyfs3(j)
         tz(j)= tz(j) + tzfs2(j) + tzfs3(j)

      end if 

   end do

   vircom = vircom + lvir

   do i=1,nm

       ! current rotation matrix

      call quat_to_rot(q0(i),q1(i),q2(i),q3(i),rott)

       ! calculate torque in principal frame

      tqx=tx(i)*rott(1)+ty(i)*rott(4)+tz(i)*rott(7)
      tqy=tx(i)*rott(2)+ty(i)*rott(5)+tz(i)*rott(8)
      tqz=tx(i)*rott(3)+ty(i)*rott(6)+tz(i)*rott(9)

       ! calculate quaternion torques

      qt0=2.0d0*(-q1(i)*tqx-q2(i)*tqy-q3(i)*tqz)
      qt1=2.0d0*( q0(i)*tqx-q3(i)*tqy+q2(i)*tqz)
      qt2=2.0d0*( q3(i)*tqx+q0(i)*tqy-q1(i)*tqz)
      qt3=2.0d0*(-q2(i)*tqx+q1(i)*tqy+q0(i)*tqz)


     ! update quaternion momenta by 1/2 time step

      p0(i)=p0(i)+qt0*0.5d0*ttt
      p1(i)=p1(i)+qt1*0.5d0*ttt
      p2(i)=p2(i)+qt2*0.5d0*ttt
      p3(i)=p3(i)+qt3*0.5d0*ttt

      !  update RB COM velocities to full step

      vx(i) = vx(i) + 0.5d0*fx(i)*ttt/totm
      vy(i) = vy(i) + 0.5d0*fy(i)*ttt/totm
      vz(i) = vz(i) + 0.5d0*fz(i)*ttt/totm

   end do

   do i=1,nm

       ! update COM velocities and angular momentum to 1/2 tstep

      opx=0.5d0*(-q1(i)*p0(i)+q0(i)*p1(i)+ &
                 q3(i)*p2(i)-q2(i)*p3(i))
      opy=0.5d0*(-q2(i)*p0(i)-q3(i)*p1(i)+ &
                 q0(i)*p2(i)+q1(i)*p3(i))
      opz=0.5d0*(-q3(i)*p0(i)+q2(i)*p1(i)- &
                 q1(i)*p2(i)+q0(i)*p3(i))

      jx(i) = opx
      jy(i) = opy
      jz(i) = opz
   end do

   do i=1,nm

      x(i) = x(i) - box*nint(x(i)/box)
      y(i) = y(i) - box*nint(y(i)/box)
      z(i) = z(i) - box*nint(z(i)/box)

   end do

 
   ltke = tke

   call tkinetic_eng (nm, totm, vx, vy, vz, ltke)

   tke = ltke
   lttemp = ttemp

   call trans_temp (nm, tke, lttemp)

   ttemp = lttemp

   stptke = tke/418.4
   stpttp = ttemp

   lrke = rke

   call rkinetic_eng (nm, jx, jy, jz, lrke)

   rke = lrke
   stprke = rke/418.4

   lrtemp = rtemp
   call rot_temp (nm, rke, lrtemp)
   rtemp = lrtemp
   stprtp = rtemp

   lstptmp = stptmp
   call system_temp(nm,tke,rke,lstptmp)
   stptmp = lstptmp

   qvir = vir+vircom

!  calculate system pressure

   pres = (2.d0*ltke - qvir)/(3.d0*volm)

!   write(*,*) box,box**3.d0,volm,pres,vir,ltke

   stpvolm = volm

   stpvir = qvir/418.4d0
   stpprs = prsunt*pres

!   write(*,*) box,box**3.d0,stpvolm,stpprs,stpvir,stptmp

   deallocate(p0)
   deallocate(p1)
   deallocate(p2)
   deallocate(p3)
        ! DBG
        deltaTimeTot = SECNDS(iniTimeTot);
!       write(dbgUnit, *) 'P#', MPIrank, '. nve time    ', deltaTimeTot;
        call MPI_BARRIER(MPI_COMM_WORLD, MPIierr);
!       STOP;
        ! DBG

  end subroutine nve_2

  subroutine nvt_b(ttt)
   implicit none
   real(8) :: ttt

   integer :: i,j,is,js
   real(8) :: qfx,qfy,qfz
   real(8) :: qtx,qty,qtz, tqx,tqy,tqz
   real(8) :: qvx,qvy,qvz
   real(8) :: qx,qy,qz
   real(8) :: qjj(3),omg(3)
   real(8) :: qq(0:3), rott(1:9)
   real(8) :: qt0,qt1,qt2,qt3
   real(8) :: opx,opy,opz
   real(8) :: cconorm,chit,lvir,qvir
   real(8) :: lvdwpe,lscpe,lrcpe,lkcpe,lindpe,lfs2pe,lfs3pe
   real(8) :: ltke,lttemp,lrke,lrtemp,lstptmp
   real(8) :: qq0,qq1,qq2,qq3,pp0,pp1,pp2,pp3
!   real(8) :: q0t(nom),q1t(nom),q2t(nom),q3t(nom)
!   real(8) :: p0(nom),p1(nom),p2(nom),p3(nom)
   real(8), allocatable :: p0(:),p1(:),p2(:),p3(:)

   real(8) :: qtke,qrke ! GRU: Unused local variables
   real(8),parameter :: taut = 0.5d0
	! DBG
	real(4) iniTime, deltaTime;
	real(4) iniTimeTot, deltaTimeTot;
	! DBG

	! DBG
	iniTimeTot = SECNDS(0.0);
	! DBG

   allocate(p0(nom))
   allocate(p1(nom))
   allocate(p2(nom))
   allocate(p3(nom))

   qtke = 0.d0
   qrke = 0.d0

!$omp parallel DEFAULT(SHARED) private(id,i,j,is,js,lvir)&
!$omp& private(qfx,qfy,qfz,qtx,qty,qtz, tqx,tqy,tqz)&
!$omp& private(qvx,qvy,qvz,qx,qy,qz,qjj,omg,qq,rott)&
!$omp& private(qt0,qt1,qt2,qt3,opx,opy,opz,cconorm,chit)&
!$omp& private(qq0,qq1,qq2,qq3,pp0,pp1,pp2,pp3)

   lvir = 0.0D0

!$omp do schedule(dynamic)
   do j=1,nm

      qfx=0.0
      qfy=0.0
      qfz=0.0
      qtx=0.0
      qty=0.0
      qtz=0.0

      do js=1,ns
         qfx=qfx+fxs(js,j)
         qfy=qfy+fys(js,j)
         qfz=qfz+fzs(js,j)
         qtx=qtx+ys(js,j)*fzs(js,j)-zs(js,j)*fys(js,j)
         qty=qty+zs(js,j)*fxs(js,j)-xs(js,j)*fzs(js,j)
         qtz=qtz+xs(js,j)*fys(js,j)-ys(js,j)*fxs(js,j)

         ! CONVERT FROM SITE-SITE TO MOLECULE-MOLECULE VIRIAL

         lvir = lvir + fxs(js,j)*xs(js,j)+fys(js,j)*ys(js,j)+fzs(js,j)*zs(js,j)  

      end do

      fx(j)=qfx
      fy(j)=qfy
      fz(j)=qfz
      tx(j)=qtx
      ty(j)=qty
      tz(j)=qtz

      if (nonadditive_3B) then

         fx(j) = fx(j) + fxfs2(j) + fxfs3(j)
         fy(j) = fy(j) + fyfs2(j) + fyfs3(j)
         fz(j) = fz(j) + fzfs2(j) + fzfs3(j)

         tx(j)= tx(j) + txfs2(j) + txfs3(j)
         ty(j)= ty(j) + tyfs2(j) + tyfs3(j)
         tz(j)= tz(j) + tzfs2(j) + tzfs3(j)

      end if 

!   end do

   !   *************  Rigid body motion ****************************
!     operations common to first and second stages

!   do i=1,nm

      ! calculate quaternion momenta at start of time step
      opx = jx(j)
      opy = jy(j)
      opz = jz(j)

      p0(j)=2.0d0*(-q1(j)*opx-q2(j)*opy-q3(j)*opz)
      p1(j)=2.0d0*( q0(j)*opx-q3(j)*opy+q2(j)*opz)
      p2(j)=2.0d0*( q3(j)*opx+q0(j)*opy-q1(j)*opz)
      p3(j)=2.0d0*(-q2(j)*opx+q1(j)*opy+q0(j)*opz)
   end do
!$omp end do
!$omp atomic

   vircom = vircom + lvir

!$omp barrier

!$omp do schedule(dynamic)
   do i=1,nm
      
      ! current rotation matrix

      call quat_to_rot(q0(i),q1(i),q2(i),q3(i),rott)

       ! calculate torque in principal frame

      tqx=tx(i)*rott(1)+ty(i)*rott(4)+tz(i)*rott(7)
      tqy=tx(i)*rott(2)+ty(i)*rott(5)+tz(i)*rott(8)
      tqz=tx(i)*rott(3)+ty(i)*rott(6)+tz(i)*rott(9)

       ! calculate quaternion torques

      qt0=2.0d0*(-q1(i)*tqx-q2(i)*tqy-q3(i)*tqz)
      qt1=2.0d0*( q0(i)*tqx-q3(i)*tqy+q2(i)*tqz)
      qt2=2.0d0*( q3(i)*tqx+q0(i)*tqy-q1(i)*tqz)
      qt3=2.0d0*(-q2(i)*tqx+q1(i)*tqy+q0(i)*tqz)


      ! calculate quaternion momenta at start of time step
!      opx = jx(i)
!      opy = jy(i)
!      opz = jz(i)

!      p0(i)=2.0d0*(-q1(i)*opx-q2(i)*opy-q3(i)*opz)
!      p1(i)=2.0d0*( q0(i)*opx-q3(i)*opy+q2(i)*opz)
!      p2(i)=2.0d0*( q3(i)*opx+q0(i)*opy-q1(i)*opz)
!      p3(i)=2.0d0*(-q2(i)*opx+q1(i)*opy+q0(i)*opz)

       ! update quaternion momenta by 1/2 time step

      p0(i)=p0(i)+qt0*0.5d0*ttt
      p1(i)=p1(i)+qt1*0.5d0*ttt
      p2(i)=p2(i)+qt2*0.5d0*ttt
      p3(i)=p3(i)+qt3*0.5d0*ttt

      ! update centre of mass velocity by 1/2 time step

      vx(i) = vx(i) + 0.5d0*fx(i)*ttt/totm
      vy(i) = vy(i) + 0.5d0*fy(i)*ttt/totm
      vz(i) = vz(i) + 0.5d0*fz(i)*ttt/totm

   end do
!$omp end do
!$omp barrier
    ! move centre of mass by full time step

!$omp do schedule(dynamic)
   do i=1,nm
      x(i)= x(i) + ttt * vx(i)
      y(i)= y(i) + ttt * vy(i)
      z(i)= z(i) + ttt * vz(i)

      ! periodic boundary conditions

!     x(i) = x(i) - box*nint(x(i)/box)
!     y(i) = y(i) - box*nint(y(i)/box)
!     z(i) = z(i) - box*nint(z(i)/box)
!   end do

      ! rotate rigid groups: nosquish algorithm -update q to full tstep and
      ! amend p and get new rotation matrix

 !  do i=1,nm  
 
      qq0 = q0(i)
      qq1 = q1(i)
      qq2 = q2(i)
      qq3 = q3(i)
      pp0 = p0(i)
      pp1 = p1(i)
      pp2 = p2(i)
      pp3 = p3(i)

      call nosquish(ttt,qq0,qq1,qq2,qq3,pp0,pp1,pp2,pp3) 

      q0(i) = qq0
      q1(i) = qq1
      q2(i) = qq2
      q3(i) = qq3
      p0(i) = pp0
      p1(i) = pp1
      p2(i) = pp2
      p3(i) = pp3

   end do

!$omp end do
!$omp barrier

!   call sites

!$omp do schedule(dynamic)
   do i=1,nm
!      qq(0) = q0(i)
!      qq(1) = q1(i)
!      qq(2) = q2(i)
!      qq(3) = q3(i)

!      call quat_to_rot(qq,rott)

       ! update COM velocities and angular momentum to 1/2 tstep

      opx=0.5d0*(-q1(i)*p0(i)+q0(i)*p1(i)+ &
                 q3(i)*p2(i)-q2(i)*p3(i))
      opy=0.5d0*(-q2(i)*p0(i)-q3(i)*p1(i)+ &
                 q0(i)*p2(i)+q1(i)*p3(i))
      opz=0.5d0*(-q3(i)*p0(i)+q2(i)*p1(i)- &
                 q1(i)*p2(i)+q0(i)*p3(i))

      jx(i) = opx
      jy(i) = opy
      jz(i) = opz
   end do
!$omp end do
!$omp end parallel

!      vx(i) = vx(i) + 0.5d0*fx(i)*tt/totm
!      vy(i) = vy(i) + 0.5d0*fy(i)*tt/totm
!      vz(i) = vz(i) + 0.5d0*fz(i)*tt/totm

     ! first stage of velocity verlet algorithm

     ! move centre of mass by full time step

!      x(i)= x(i) + tt * vx(i)
!      y(i)= y(i) + tt * vy(i)
!      z(i)= z(i) + tt * vz(i)

      ! periodic boundary conditions

!    do i=1,nm

!      x(i) = x(i) - box*nint(x(i)/box)
!      y(i) = y(i) - box*nint(y(i)/box)
!      z(i) = z(i) - box*nint(z(i)/box)

!   end do

! end of first stage of velocity verlet algorithm

   call sites

   lvdwpe = vdwpe

	! DBG
	iniTime = SECNDS(0.0);
	! DBG
   if (pot_name.eq.'tip4p' ) then 
       call lj(x,y,z,xs,ys,zs,box,lvdwpe)
   else if (pot_name.eq.'tip5p' ) then
       call lj(x,y,z,xs,ys,zs,box,lvdwpe)
   else if (pot_name.eq.'sapt5s' ) then
       call sapt5s(x,y,z,xs,ys,zs,box,lvdwpe)
   else if (pot_name.eq.'ccpol5s' ) then
       call ccpol5s(x,y,z,xs,ys,zs,box,lvdwpe)
   else if (pot_name.eq.'ccpol8s' ) then
       call ccpol8s(x,y,z,xs,ys,zs,box,lvdwpe)       
   else if (pot_name.eq.'ccpol8s_omo' ) then
       call ccpol8s_omo(x,y,z,xs,ys,zs,box,lvdwpe)       
   else if (pot_name.eq.'ccpol23plus' ) then
       call ccpol23plus(x,y,z,xs,ys,zs,box,lvdwpe)          
		! DBG
		deltaTime = SECNDS(iniTime);
!		write(dbgUnit, *) 'P#', MPIrank, '. ccpol23plus time	', deltaTime;
		! DBG
   end if

   vdwpe = lvdwpe

   lscpe = scpe
   call ewald_self_correction(alpha,box,lscpe)
   scpe = lscpe

   lrcpe = rcpe
   call rspace_ewald(x,y,z,xs,ys,zs,box,alpha,lrcpe)
   rcpe = lrcpe

   lkcpe = kcpe
   call kspace_ewald(x,y,z,xs,ys,zs,kmax1,kmax2,kmax3,alpha,box,lkcpe)
   kcpe = lkcpe

!  calling N-Body Iterative Induction model
   if (induction) then
      if(ind_model .eq. 'point_dipole') then

        lindpe = indpe
        call nb_induction(x,y,z,xs,ys,zs,box,lindpe)
        indpe = lindpe

      end if
   end if 

! calling 3B-potenital (FS2 and FS3 terms)

  if (nonadditive_3B) then

     lfs2pe = fs2pe
     lfs3pe = fs3pe
     call threebody_term(x,y,z,xs,ys,zs,box,lfs2pe,lfs3pe)
     fs2pe = lfs2pe
     fs3pe = lfs3pe

  end if  

   ! second stage of velocity verlet algorithm

!$omp parallel DEFAULT(SHARED) private(id,i,j,is,js,lvir)&
!$omp& private(qfx,qfy,qfz,qtx,qty,qtz, tqx,tqy,tqz)&
!$omp& private(qvx,qvy,qvz,qx,qy,qz,qjj,omg,qq,rott)&
!$omp& private(qt0,qt1,qt2,qt3,opx,opy,opz,cconorm,chit)&
!$omp& private(qq0,qq1,qq2,qq3,pp0,pp1,pp2,pp3)

   lvir = 0.0D0

!$omp do schedule(dynamic)
   do j=1,nm

      qfx=0.0
      qfy=0.0
      qfz=0.0
      qtx=0.0
      qty=0.0
      qtz=0.0

      do js=1,ns
         qfx=qfx+fxs(js,j)
         qfy=qfy+fys(js,j)
         qfz=qfz+fzs(js,j)
         qtx=qtx+ys(js,j)*fzs(js,j)-zs(js,j)*fys(js,j)
         qty=qty+zs(js,j)*fxs(js,j)-xs(js,j)*fzs(js,j)
         qtz=qtz+xs(js,j)*fys(js,j)-ys(js,j)*fxs(js,j)

         ! CONVERT FROM SITE-SITE TO MOLECULE-MOLECULE VIRIAL

         lvir = lvir + fxs(js,j)*xs(js,j)+fys(js,j)*ys(js,j)+fzs(js,j)*zs(js,j)

      end do

      fx(j)=qfx
      fy(j)=qfy
      fz(j)=qfz
      tx(j)=qtx
      ty(j)=qty
      tz(j)=qtz

      if (nonadditive_3B) then

         fx(j) = fx(j) + fxfs2(j) + fxfs3(j)
         fy(j) = fy(j) + fyfs2(j) + fyfs3(j)
         fz(j) = fz(j) + fzfs2(j) + fzfs3(j)

         tx(j)= tx(j) + txfs2(j) + txfs3(j)
         ty(j)= ty(j) + tyfs2(j) + tyfs3(j)
         tz(j)= tz(j) + tzfs2(j) + tzfs3(j)

      end if 

   end do
!$omp end do
!$omp atomic
   vircom = vircom + lvir

!$omp barrier


!$omp do schedule(dynamic)
   do i=1,nm

      ! calculate quaternion momenta at start of time step
      opx = jx(i)
      opy = jy(i)
      opz = jz(i)

      p0(i)=2.0d0*(-q1(i)*opx-q2(i)*opy-q3(i)*opz)
      p1(i)=2.0d0*( q0(i)*opx-q3(i)*opy+q2(i)*opz)
      p2(i)=2.0d0*( q3(i)*opx+q0(i)*opy-q1(i)*opz)
      p3(i)=2.0d0*(-q2(i)*opx+q1(i)*opy+q0(i)*opz)
   end do
!$omp end do
!$omp barrier

!$omp do schedule(dynamic)
   do i=1,nm

       ! current rotation matrix

      call quat_to_rot(q0(i),q1(i),q2(i),q3(i),rott)

       ! calculate torque in principal frame

      tqx=tx(i)*rott(1)+ty(i)*rott(4)+tz(i)*rott(7)
      tqy=tx(i)*rott(2)+ty(i)*rott(5)+tz(i)*rott(8)
      tqz=tx(i)*rott(3)+ty(i)*rott(6)+tz(i)*rott(9)

       ! calculate quaternion torques

      qt0=2.0d0*(-q1(i)*tqx-q2(i)*tqy-q3(i)*tqz)
      qt1=2.0d0*( q0(i)*tqx-q3(i)*tqy+q2(i)*tqz)
      qt2=2.0d0*( q3(i)*tqx+q0(i)*tqy-q1(i)*tqz)
      qt3=2.0d0*(-q2(i)*tqx+q1(i)*tqy+q0(i)*tqz)

      ! calculate quaternion momenta at start of time step
!      opx = jx(i)
!      opy = jy(i)
!      opz = jz(i)

!      p0(i)=2.0d0*(-q1(i)*opx-q2(i)*opy-q3(i)*opz)
!      p1(i)=2.0d0*( q0(i)*opx-q3(i)*opy+q2(i)*opz)
!      p2(i)=2.0d0*( q3(i)*opx+q0(i)*opy-q1(i)*opz)
!      p3(i)=2.0d0*(-q2(i)*opx+q1(i)*opy+q0(i)*opz)

       ! update quaternion momenta by 1/2 time step

      p0(i)=p0(i)+qt0*0.5d0*ttt
      p1(i)=p1(i)+qt1*0.5d0*ttt
      p2(i)=p2(i)+qt2*0.5d0*ttt
      p3(i)=p3(i)+qt3*0.5d0*ttt

      !  update Rgid body angular momentum & COM velocities to full step

!      opx=0.5d0*(-q1(i)*p0(i)+q0(i)*p1(i)+ &
!                 q3(i)*p2(i)-q2(i)*p3(i))
!      opy=0.5d0*(-q2(i)*p0(i)-q3(i)*p1(i)+ &
!                 q0(i)*p2(i)+q1(i)*p3(i))
!      opz=0.5d0*(-q3(i)*p0(i)+q2(i)*p1(i)- &
!                 q1(i)*p2(i)+q0(i)*p3(i))

!      jx(i) = opx
!      jy(i) = opy
!      jz(i) = opz

      vx(i) = vx(i) + 0.5d0*fx(i)*ttt/totm
      vy(i) = vy(i) + 0.5d0*fy(i)*ttt/totm
      vz(i) = vz(i) + 0.5d0*fz(i)*ttt/totm

   end do
!$omp end do
!$omp barrier

!$omp do schedule(dynamic)
   do i=1,nm
!      qq(0) = q0(i)
!      qq(1) = q1(i)
!      qq(2) = q2(i)
!      qq(3) = q3(i)

!      call quat_to_rot(qq,rott)

       ! update COM velocities and angular momentum to 1/2 tstep

      opx=0.5d0*(-q1(i)*p0(i)+q0(i)*p1(i)+ &
                 q3(i)*p2(i)-q2(i)*p3(i))
      opy=0.5d0*(-q2(i)*p0(i)-q3(i)*p1(i)+ &
                 q0(i)*p2(i)+q1(i)*p3(i))
      opz=0.5d0*(-q3(i)*p0(i)+q2(i)*p1(i)- &
                 q1(i)*p2(i)+q0(i)*p3(i))

      jx(i) = opx
      jy(i) = opy
      jz(i) = opz
   end do
!$omp end do
!$omp end parallel

! apply Berendsen thermostat - taut is the relaxation time

   ltke = tke
   call tkinetic_eng (nm, totm, vx, vy, vz, ltke)
   tke = ltke

   lttemp = ttemp
   call trans_temp (nm, tke, lttemp)
   ttemp = lttemp

   stptke = tke/418.4
   stpttp = ttemp

   lrke = rke
   call rkinetic_eng (nm, jx, jy, jz, lrke)
   rke = lrke
   stprke = rke/418.4

   lrtemp = rtemp
   call rot_temp (nm, rke, lrtemp)
   rtemp = lrtemp
   stprtp = rtemp

   lstptmp = stptmp
   call system_temp(nm,tke,rke,lstptmp)
   stptmp = lstptmp

   write(*,*) temp0,stptmp

   chit = sqrt(1.d0+ ttt/taut*(temp0/stptmp - 1.d0)) 

! thermostat velocities

   do i=1,nm

      vx(i) = chit*vx(i) 
      vy(i) = chit*vy(i) 
      vz(i) = chit*vz(i)

      jx(i) = chit*jx(i)
      jy(i) = chit*jy(i)
      jz(i) = chit*jz(i)


   end do
   
   ltke = tke
   call tkinetic_eng (nm, totm, vx, vy, vz, ltke)
   tke = ltke

   lttemp = ttemp
   call trans_temp (nm, tke, lttemp)
   ttemp = lttemp

   stptke = tke/418.4
   stpttp = ttemp

   lrke = rke
   call rkinetic_eng (nm, jx, jy, jz, lrke)
   rke = lrke
   stprke = rke/418.4

   lrtemp = rtemp
   call rot_temp (nm, rke, lrtemp)
   rtemp = lrtemp
   stprtp = rtemp

   lstptmp = stptmp
   call system_temp(nm,tke,rke,lstptmp)
   stptmp = lstptmp

!   write(*,*) temp0,stptmp

   qvir = vir+vircom

!  calculate system pressure

   pres = (2.d0*ltke - qvir)/(3.d0*volm)

!   write(*,*) box,box**3.d0,volm,pres,vir,ltke

   stpvolm = volm

   stpvir = qvir/418.4d0
   stpprs = prsunt*pres

   deallocate(p0)
   deallocate(p1)
   deallocate(p2)
   deallocate(p3)
	! DBG
	deltaTimeTot = SECNDS(iniTimeTot);
!	write(dbgUnit, *) 'P#', MPIrank, '. nvt_b time	', deltaTimeTot;
	call MPI_BARRIER(MPI_COMM_WORLD, MPIierr);
!	STOP;
	! DBG

  end subroutine nvt_b  


  subroutine nvt_b_2(ttt)
   implicit none
   real(8) :: ttt

   integer :: i,j,is,js
   real(8) :: qfx,qfy,qfz
   real(8) :: qtx,qty,qtz, tqx,tqy,tqz
   real(8) :: qvx,qvy,qvz
   real(8) :: qx,qy,qz
   real(8) :: qjj(3),omg(3)
   real(8) :: qq(0:3), rott(1:9)
   real(8) :: qt0,qt1,qt2,qt3
   real(8) :: opx,opy,opz
   real(8) :: cconorm,chit,lvir,qvir
   real(8) :: lvdwpe,lscpe,lrcpe,lkcpe,lindpe,lfs2pe,lfs3pe
   real(8) :: ltke,lttemp,lrke,lrtemp,lstptmp
   real(8) :: qq0,qq1,qq2,qq3,pp0,pp1,pp2,pp3
!   real(8) :: q0t(nom),q1t(nom),q2t(nom),q3t(nom)
!   real(8) :: p0(nom),p1(nom),p2(nom),p3(nom)
   real(8), allocatable :: p0(:),p1(:),p2(:),p3(:)

   real(8) :: qtke,qrke ! GRU: Unused local variables
   real(8),parameter :: taut = 0.5d0
        ! DBG
        real(4) iniTime, deltaTime;
        real(4) iniTimeTot, deltaTimeTot;
        ! DBG

        ! DBG
        iniTimeTot = SECNDS(0.0);
        ! DBG

   allocate(p0(nom))
   allocate(p1(nom))
   allocate(p2(nom))
   allocate(p3(nom))

   qtke = 0.d0
   qrke = 0.d0

!$omp parallel DEFAULT(SHARED) private(id,i,j,is,js,lvir)&
!$omp& private(qfx,qfy,qfz,qtx,qty,qtz, tqx,tqy,tqz)&
!$omp& private(qvx,qvy,qvz,qx,qy,qz,qjj,omg,qq,rott)&
!$omp& private(qt0,qt1,qt2,qt3,opx,opy,opz,cconorm,chit)&
!$omp& private(qq0,qq1,qq2,qq3,pp0,pp1,pp2,pp3)

   lvir = 0.0D0

!$omp do schedule(dynamic)
   do j=1,nm

      qfx=0.0
      qfy=0.0
      qfz=0.0
      qtx=0.0
      qty=0.0
      qtz=0.0

      do js=1,ns
         qfx=qfx+fxs(js,j)
         qfy=qfy+fys(js,j)
         qfz=qfz+fzs(js,j)
         qtx=qtx+ys(js,j)*fzs(js,j)-zs(js,j)*fys(js,j)
         qty=qty+zs(js,j)*fxs(js,j)-xs(js,j)*fzs(js,j)
         qtz=qtz+xs(js,j)*fys(js,j)-ys(js,j)*fxs(js,j)

         ! CONVERT FROM SITE-SITE TO MOLECULE-MOLECULE VIRIAL

         lvir = lvir + fxs(js,j)*xs(js,j)+fys(js,j)*ys(js,j)+fzs(js,j)*zs(js,j)  

      end do

      fx(j)=qfx
      fy(j)=qfy
      fz(j)=qfz
      tx(j)=qtx
      ty(j)=qty
      tz(j)=qtz

      if (nonadditive_3B) then

         fx(j) = fx(j) + fxfs2(j) + fxfs3(j)
         fy(j) = fy(j) + fyfs2(j) + fyfs3(j)
         fz(j) = fz(j) + fzfs2(j) + fzfs3(j)

         tx(j)= tx(j) + txfs2(j) + txfs3(j)
         ty(j)= ty(j) + tyfs2(j) + tyfs3(j)
         tz(j)= tz(j) + tzfs2(j) + tzfs3(j)

      end if 

!   end do

   !   *************  Rigid body motion ****************************
!     operations common to first and second stages

!   do i=1,nm

      ! calculate quaternion momenta at start of time step
      opx = jx(j)
      opy = jy(j)
      opz = jz(j)

      p0(j)=2.0d0*(-q1(j)*opx-q2(j)*opy-q3(j)*opz)
      p1(j)=2.0d0*( q0(j)*opx-q3(j)*opy+q2(j)*opz)
      p2(j)=2.0d0*( q3(j)*opx+q0(j)*opy-q1(j)*opz)
      p3(j)=2.0d0*(-q2(j)*opx+q1(j)*opy+q0(j)*opz)
   end do
!$omp end do
!$omp atomic

   vircom = vircom + lvir

!$omp barrier

!$omp do schedule(dynamic)
   do i=1,nm
      
      ! current rotation matrix

      call quat_to_rot(q0(i),q1(i),q2(i),q3(i),rott)

       ! calculate torque in principal frame

      tqx=tx(i)*rott(1)+ty(i)*rott(4)+tz(i)*rott(7)
      tqy=tx(i)*rott(2)+ty(i)*rott(5)+tz(i)*rott(8)
      tqz=tx(i)*rott(3)+ty(i)*rott(6)+tz(i)*rott(9)

       ! calculate quaternion torques

      qt0=2.0d0*(-q1(i)*tqx-q2(i)*tqy-q3(i)*tqz)
      qt1=2.0d0*( q0(i)*tqx-q3(i)*tqy+q2(i)*tqz)
      qt2=2.0d0*( q3(i)*tqx+q0(i)*tqy-q1(i)*tqz)
      qt3=2.0d0*(-q2(i)*tqx+q1(i)*tqy+q0(i)*tqz)


      ! calculate quaternion momenta at start of time step
!      opx = jx(i)
!      opy = jy(i)
!      opz = jz(i)

!      p0(i)=2.0d0*(-q1(i)*opx-q2(i)*opy-q3(i)*opz)
!      p1(i)=2.0d0*( q0(i)*opx-q3(i)*opy+q2(i)*opz)
!      p2(i)=2.0d0*( q3(i)*opx+q0(i)*opy-q1(i)*opz)
!      p3(i)=2.0d0*(-q2(i)*opx+q1(i)*opy+q0(i)*opz)

       ! update quaternion momenta by 1/2 time step

      p0(i)=p0(i)+qt0*0.5d0*ttt
      p1(i)=p1(i)+qt1*0.5d0*ttt
      p2(i)=p2(i)+qt2*0.5d0*ttt
      p3(i)=p3(i)+qt3*0.5d0*ttt

      ! update centre of mass velocity by 1/2 time step

      vx(i) = vx(i) + 0.5d0*fx(i)*ttt/totm
      vy(i) = vy(i) + 0.5d0*fy(i)*ttt/totm
      vz(i) = vz(i) + 0.5d0*fz(i)*ttt/totm

   end do
!$omp end do
!$omp barrier
    ! move centre of mass by full time step

!$omp do schedule(dynamic)
   do i=1,nm
      x(i)= x(i) + ttt * vx(i)
      y(i)= y(i) + ttt * vy(i)
      z(i)= z(i) + ttt * vz(i)

      ! periodic boundary conditions

!     x(i) = x(i) - box*nint(x(i)/box)
!     y(i) = y(i) - box*nint(y(i)/box)
!     z(i) = z(i) - box*nint(z(i)/box)
!   end do

      ! rotate rigid groups: nosquish algorithm -update q to full tstep and
      ! amend p and get new rotation matrix

 !  do i=1,nm  
 
      qq0 = q0(i)
      qq1 = q1(i)
      qq2 = q2(i)
      qq3 = q3(i)
      pp0 = p0(i)
      pp1 = p1(i)
      pp2 = p2(i)
      pp3 = p3(i)

      call nosquish(ttt,qq0,qq1,qq2,qq3,pp0,pp1,pp2,pp3) 

      q0(i) = qq0
      q1(i) = qq1
      q2(i) = qq2
      q3(i) = qq3
      p0(i) = pp0
      p1(i) = pp1
      p2(i) = pp2
      p3(i) = pp3

   end do

!$omp end do
!$omp barrier

!   call sites

!$omp do schedule(dynamic)
!!   do i=1,nm
!      qq(0) = q0(i)
!      qq(1) = q1(i)
!      qq(2) = q2(i)
!      qq(3) = q3(i)

!      call quat_to_rot(qq,rott)

       ! update COM velocities and angular momentum to 1/2 tstep

!!      opx=0.5d0*(-q1(i)*p0(i)+q0(i)*p1(i)+ &
!!                 q3(i)*p2(i)-q2(i)*p3(i))
!      opy=0.5d0*(-q2(i)*p0(i)-q3(i)*p1(i)+ &
!                 q0(i)*p2(i)+q1(i)*p3(i))
!      opz=0.5d0*(-q3(i)*p0(i)+q2(i)*p1(i)- &
!                 q1(i)*p2(i)+q0(i)*p3(i))

!      jx(i) = opx
!      jy(i) = opy
!      jz(i) = opz
!   end do

!$omp end do
!$omp end parallel

!      vx(i) = vx(i) + 0.5d0*fx(i)*tt/totm
!      vy(i) = vy(i) + 0.5d0*fy(i)*tt/totm
!      vz(i) = vz(i) + 0.5d0*fz(i)*tt/totm

     ! first stage of velocity verlet algorithm

     ! move centre of mass by full time step

!      x(i)= x(i) + tt * vx(i)
!      y(i)= y(i) + tt * vy(i)
!      z(i)= z(i) + tt * vz(i)

      ! periodic boundary conditions

!    do i=1,nm

!      x(i) = x(i) - box*nint(x(i)/box)
!      y(i) = y(i) - box*nint(y(i)/box)
!      z(i) = z(i) - box*nint(z(i)/box)

!   end do

! end of first stage of velocity verlet algorithm

   call sites

   lvdwpe = vdwpe

        ! DBG
        iniTime = SECNDS(0.0);
        ! DBG
   if (pot_name.eq.'tip4p' ) then 
       call lj(x,y,z,xs,ys,zs,box,lvdwpe)
   else if (pot_name.eq.'tip5p' ) then
       call lj(x,y,z,xs,ys,zs,box,lvdwpe)
   else if (pot_name.eq.'sapt5s' ) then
       call sapt5s(x,y,z,xs,ys,zs,box,lvdwpe)
   else if (pot_name.eq.'ccpol5s' ) then
       call ccpol5s(x,y,z,xs,ys,zs,box,lvdwpe)
   else if (pot_name.eq.'ccpol8s' ) then
       call ccpol8s(x,y,z,xs,ys,zs,box,lvdwpe)       
   else if (pot_name.eq.'ccpol8s_omo' ) then
       call ccpol8s_omo(x,y,z,xs,ys,zs,box,lvdwpe)       
   else if (pot_name.eq.'ccpol23plus' ) then
       call ccpol23plus(x,y,z,xs,ys,zs,box,lvdwpe)          
                ! DBG
                deltaTime = SECNDS(iniTime);
!               write(dbgUnit, *) 'P#', MPIrank, '. ccpol23plus time    ',
!               deltaTime;
                ! DBG
   end if

   vdwpe = lvdwpe

   lscpe = scpe
   call ewald_self_correction(alpha,box,lscpe)
   scpe = lscpe

   lrcpe = rcpe
   call rspace_ewald(x,y,z,xs,ys,zs,box,alpha,lrcpe)
   rcpe = lrcpe

   lkcpe = kcpe
   call kspace_ewald(x,y,z,xs,ys,zs,kmax1,kmax2,kmax3,alpha,box,lkcpe)
   kcpe = lkcpe

!  calling N-Body Iterative Induction model
   if (induction) then
      if(ind_model .eq. 'point_dipole') then

        lindpe = indpe
        call nb_induction(x,y,z,xs,ys,zs,box,lindpe)
        indpe = lindpe

      end if
   end if 

! calling 3B-potenital (FS2 and FS3 terms)

  if (nonadditive_3B) then

     lfs2pe = fs2pe
     lfs3pe = fs3pe
     call threebody_term(x,y,z,xs,ys,zs,box,lfs2pe,lfs3pe)
     fs2pe = lfs2pe
     fs3pe = lfs3pe

  end if  

   ! second stage of velocity verlet algorithm

!$omp parallel DEFAULT(SHARED) private(id,i,j,is,js,lvir)&
!$omp& private(qfx,qfy,qfz,qtx,qty,qtz, tqx,tqy,tqz)&
!$omp& private(qvx,qvy,qvz,qx,qy,qz,qjj,omg,qq,rott)&
!$omp& private(qt0,qt1,qt2,qt3,opx,opy,opz,cconorm,chit)&
!$omp& private(qq0,qq1,qq2,qq3,pp0,pp1,pp2,pp3)

   lvir = 0.0D0

!$omp do schedule(dynamic)
   do j=1,nm

      qfx=0.0
      qfy=0.0
      qfz=0.0
      qtx=0.0
      qty=0.0
      qtz=0.0

      do js=1,ns
         qfx=qfx+fxs(js,j)
         qfy=qfy+fys(js,j)
         qfz=qfz+fzs(js,j)
         qtx=qtx+ys(js,j)*fzs(js,j)-zs(js,j)*fys(js,j)
         qty=qty+zs(js,j)*fxs(js,j)-xs(js,j)*fzs(js,j)
         qtz=qtz+xs(js,j)*fys(js,j)-ys(js,j)*fxs(js,j)

         ! CONVERT FROM SITE-SITE TO MOLECULE-MOLECULE VIRIAL

         lvir = lvir + fxs(js,j)*xs(js,j)+fys(js,j)*ys(js,j)+fzs(js,j)*zs(js,j)

      end do

      fx(j)=qfx
      fy(j)=qfy
      fz(j)=qfz
      tx(j)=qtx
      ty(j)=qty
      tz(j)=qtz

      if (nonadditive_3B) then

         fx(j) = fx(j) + fxfs2(j) + fxfs3(j)
         fy(j) = fy(j) + fyfs2(j) + fyfs3(j)
         fz(j) = fz(j) + fzfs2(j) + fzfs3(j)

         tx(j)= tx(j) + txfs2(j) + txfs3(j)
         ty(j)= ty(j) + tyfs2(j) + tyfs3(j)
         tz(j)= tz(j) + tzfs2(j) + tzfs3(j)

      end if 

   end do
!$omp end do
!$omp atomic
   vircom = vircom + lvir

!$omp barrier


!$omp do schedule(dynamic)
!   do i=1,nm

      ! calculate quaternion momenta at start of time step
!      opx = jx(i)
!      opy = jy(i)
!      opz = jz(i)

!      p0(i)=2.0d0*(-q1(i)*opx-q2(i)*opy-q3(i)*opz)
!      p1(i)=2.0d0*( q0(i)*opx-q3(i)*opy+q2(i)*opz)
!      p2(i)=2.0d0*( q3(i)*opx+q0(i)*opy-q1(i)*opz)
!      p3(i)=2.0d0*(-q2(i)*opx+q1(i)*opy+q0(i)*opz)
!   end do
!$omp end do
!$omp barrier

!$omp do schedule(dynamic)
   do i=1,nm

       ! current rotation matrix

      call quat_to_rot(q0(i),q1(i),q2(i),q3(i),rott)

       ! calculate torque in principal frame

      tqx=tx(i)*rott(1)+ty(i)*rott(4)+tz(i)*rott(7)
      tqy=tx(i)*rott(2)+ty(i)*rott(5)+tz(i)*rott(8)
      tqz=tx(i)*rott(3)+ty(i)*rott(6)+tz(i)*rott(9)

       ! calculate quaternion torques

      qt0=2.0d0*(-q1(i)*tqx-q2(i)*tqy-q3(i)*tqz)
      qt1=2.0d0*( q0(i)*tqx-q3(i)*tqy+q2(i)*tqz)
      qt2=2.0d0*( q3(i)*tqx+q0(i)*tqy-q1(i)*tqz)
      qt3=2.0d0*(-q2(i)*tqx+q1(i)*tqy+q0(i)*tqz)

      ! calculate quaternion momenta at start of time step
!      opx = jx(i)
!      opy = jy(i)
!      opz = jz(i)

!      p0(i)=2.0d0*(-q1(i)*opx-q2(i)*opy-q3(i)*opz)
!      p1(i)=2.0d0*( q0(i)*opx-q3(i)*opy+q2(i)*opz)
!      p2(i)=2.0d0*( q3(i)*opx+q0(i)*opy-q1(i)*opz)
!      p3(i)=2.0d0*(-q2(i)*opx+q1(i)*opy+q0(i)*opz)

       ! update quaternion momenta by 1/2 time step

      p0(i)=p0(i)+qt0*0.5d0*ttt
      p1(i)=p1(i)+qt1*0.5d0*ttt
      p2(i)=p2(i)+qt2*0.5d0*ttt
      p3(i)=p3(i)+qt3*0.5d0*ttt

      !  update Rgid body angular momentum & COM velocities to full step

!      opx=0.5d0*(-q1(i)*p0(i)+q0(i)*p1(i)+ &
!                 q3(i)*p2(i)-q2(i)*p3(i))
!      opy=0.5d0*(-q2(i)*p0(i)-q3(i)*p1(i)+ &
!                 q0(i)*p2(i)+q1(i)*p3(i))
!      opz=0.5d0*(-q3(i)*p0(i)+q2(i)*p1(i)- &
!                 q1(i)*p2(i)+q0(i)*p3(i))

!      jx(i) = opx
!      jy(i) = opy
!      jz(i) = opz

      vx(i) = vx(i) + 0.5d0*fx(i)*ttt/totm
      vy(i) = vy(i) + 0.5d0*fy(i)*ttt/totm
      vz(i) = vz(i) + 0.5d0*fz(i)*ttt/totm

   end do
!$omp end do
!$omp barrier

!$omp do schedule(dynamic)
   do i=1,nm
!      qq(0) = q0(i)
!      qq(1) = q1(i)
!      qq(2) = q2(i)
!      qq(3) = q3(i)

!      call quat_to_rot(qq,rott)

       ! update COM velocities and angular momentum to 1/2 tstep

      opx=0.5d0*(-q1(i)*p0(i)+q0(i)*p1(i)+ &
                 q3(i)*p2(i)-q2(i)*p3(i))
      opy=0.5d0*(-q2(i)*p0(i)-q3(i)*p1(i)+ &
                 q0(i)*p2(i)+q1(i)*p3(i))
      opz=0.5d0*(-q3(i)*p0(i)+q2(i)*p1(i)- &
                 q1(i)*p2(i)+q0(i)*p3(i))

      jx(i) = opx
      jy(i) = opy
      jz(i) = opz
   end do
!$omp end do
!$omp end parallel

! apply Berendsen thermostat - taut is the relaxation time

   ltke = tke
   call tkinetic_eng (nm, totm, vx, vy, vz, ltke)
   tke = ltke

   lttemp = ttemp
   call trans_temp (nm, tke, lttemp)
   ttemp = lttemp

   stptke = tke/418.4
   stpttp = ttemp

   lrke = rke
   call rkinetic_eng (nm, jx, jy, jz, lrke)
   rke = lrke
   stprke = rke/418.4

   lrtemp = rtemp
   call rot_temp (nm, rke, lrtemp)
   rtemp = lrtemp
   stprtp = rtemp

   lstptmp = stptmp
   call system_temp(nm,tke,rke,lstptmp)
   stptmp = lstptmp

   write(*,*) temp0,stptmp

   chit = sqrt(1.d0+ ttt/taut*(temp0/stptmp - 1.d0)) 

! thermostat velocities

   do i=1,nm

      vx(i) = chit*vx(i) 
      vy(i) = chit*vy(i) 
      vz(i) = chit*vz(i)

      jx(i) = chit*jx(i)
      jy(i) = chit*jy(i)
      jz(i) = chit*jz(i)


   end do
   
   ltke = tke
   call tkinetic_eng (nm, totm, vx, vy, vz, ltke)
   tke = ltke

   lttemp = ttemp
   call trans_temp (nm, tke, lttemp)
   ttemp = lttemp

   stptke = tke/418.4
   stpttp = ttemp

   lrke = rke
   call rkinetic_eng (nm, jx, jy, jz, lrke)
   rke = lrke
   stprke = rke/418.4

   lrtemp = rtemp
   call rot_temp (nm, rke, lrtemp)
   rtemp = lrtemp
   stprtp = rtemp

   lstptmp = stptmp
   call system_temp(nm,tke,rke,lstptmp)
   stptmp = lstptmp

!   write(*,*) temp0,stptmp

   qvir = vir+vircom

!  calculate system pressure

   pres = (2.d0*ltke - qvir)/(3.d0*volm)

!   write(*,*) box,box**3.d0,volm,pres,vir,ltke

   stpvolm = volm

   stpvir = qvir/418.4d0
   stpprs = prsunt*pres

   deallocate(p0)
   deallocate(p1)
   deallocate(p2)
   deallocate(p3)
        ! DBG
        deltaTimeTot = SECNDS(iniTimeTot);
!       write(dbgUnit, *) 'P#', MPIrank, '. nvt_b time  ', deltaTimeTot;
        call MPI_BARRIER(MPI_COMM_WORLD, MPIierr);
!       STOP;
        ! DBG

  end subroutine nvt_b_2

  subroutine nvt_hov(ttt)
   implicit none
   real(8) :: ttt

   integer :: i,j,is,js
   real(8) :: qfx,qfy,qfz
   real(8) :: qtx,qty,qtz, tqx,tqy,tqz
   real(8) :: qvx,qvy,qvz
   real(8) :: qx,qy,qz
   real(8) :: qjj(3),omg(3)
   real(8) :: qq(0:3), rott(1:9)
   real(8) :: qt0,qt1,qt2,qt3
   real(8) :: opx,opy,opz
   real(8) :: cconorm,lvir,qvir
   real(8) :: lvdwpe,lscpe,lrcpe,lkcpe,lindpe,lfs2pe,lfs3pe
   real(8) :: ltke,lttemp,lrke,lrtemp,lstptmp
   real(8) :: qq0,qq1,qq2,qq3,pp0,pp1,pp2,pp3
!   real(8) :: q0t(nom),q1t(nom),q2t(nom),q3t(nom)
!   real(8) :: p0(nom),p1(nom),p2(nom),p3(nom)
    real(8), allocatable :: p0(:),p1(:),p2(:),p3(:)

! Nose Hoover parameters
   real(8),parameter :: taut = 0.5d0 ! ps (time constant)
   real(8) :: qmass,hstep,qtke,qrke,qsystmp,dof,chit,conint
   real(8) :: consv  
	! DBG
	real(4) iniTime, deltaTime;
	real(4) iniTimeTot, deltaTimeTot;
	! DBG

	! DBG
	iniTimeTot = SECNDS(0.0);
	! DBG

   allocate(p0(nm))
   allocate(p1(nm))
   allocate(p2(nm))
   allocate(p3(nm))

   ! GRU
   chit = 0.0D0
   conint = 0.0D0
   ! GRU

!$omp parallel DEFAULT(SHARED) private(id,i,j,is,js,lvir)&
!$omp& private(qfx,qfy,qfz,qtx,qty,qtz, tqx,tqy,tqz)&
!$omp& private(qvx,qvy,qvz,qx,qy,qz,qjj,omg,qq,rott)&
!$omp& private(qt0,qt1,qt2,qt3,opx,opy,opz,cconorm,qmass)&
!$omp& private(hstep,qtke,qrke,qsystmp,dof,chit,conint,consv)&
!$omp& private(qq0,qq1,qq2,qq3,pp0,pp1,pp2,pp3)

   lvir = 0.0D0

!$omp do schedule(dynamic)

   do j=1,nm

      qfx=0.0
      qfy=0.0
      qfz=0.0
      qtx=0.0
      qty=0.0
      qtz=0.0

      do js=1,ns
         qfx=qfx+fxs(js,j)
         qfy=qfy+fys(js,j)
         qfz=qfz+fzs(js,j)
         qtx=qtx+ys(js,j)*fzs(js,j)-zs(js,j)*fys(js,j)
         qty=qty+zs(js,j)*fxs(js,j)-xs(js,j)*fzs(js,j)
         qtz=qtz+xs(js,j)*fys(js,j)-ys(js,j)*fxs(js,j)

         ! CONVERT FROM SITE-SITE TO MOLECULE-MOLECULE VIRIAL

         lvir = lvir + fxs(js,j)*xs(js,j)+fys(js,j)*ys(js,j)+fzs(js,j)*zs(js,j)

      end do

      fx(j)=qfx
      fy(j)=qfy
      fz(j)=qfz
      tx(j)=qtx
      ty(j)=qty
      tz(j)=qtz

      if (nonadditive_3B) then

         fx(j) = fx(j) + fxfs2(j) + fxfs3(j)
         fy(j) = fy(j) + fyfs2(j) + fyfs3(j)
         fz(j) = fz(j) + fzfs2(j) + fzfs3(j)

         tx(j)= tx(j) + txfs2(j) + txfs3(j)
         ty(j)= ty(j) + tyfs2(j) + tyfs3(j)
         tz(j)= tz(j) + tzfs2(j) + tzfs3(j)


      end if

   end do
!$omp end do
!$omp atomic

   vircom = vircom + lvir

!$omp barrier
!$omp end parallel

!!$omp master
!  timestep parameters  
   hstep = 0.5d0*ttt

  dof = 3.d0*(2.d0*nm-1.d0) ! total dof transl + rot

!   dof = 3.d0*(nm-1.d0) ! total dof transl + rot 
!   nose-hoover inertia parameter
   qmass=dof*boltz*temp0*taut**2

   call nvtqscl(hstep,qmass,dof,taut,chit,conint)

!!$omp end master
!!$omp barrier

!$omp parallel DEFAULT(SHARED) private(id,i,j,is,js,lvir)&
!$omp& private(qfx,qfy,qfz,qtx,qty,qtz, tqx,tqy,tqz)&
!$omp& private(qvx,qvy,qvz,qx,qy,qz,qjj,omg,qq,rott)&
!$omp& private(qt0,qt1,qt2,qt3,opx,opy,opz,cconorm,qmass)&
!$omp& private(hstep,qtke,qrke,qsystmp,dof,chit,conint,consv)&
!$omp& private(qq0,qq1,qq2,qq3,pp0,pp1,pp2,pp3)

!   *************  Rigid body motion ****************************
!     operations common to first and second stages

!$omp do schedule(dynamic)
   do i=1,nm

      ! calculate quaternion momenta at start of time step
      opx = jx(i)
      opy = jy(i)
      opz = jz(i)

      p0(i)=2.0d0*(-q1(i)*opx-q2(i)*opy-q3(i)*opz)
      p1(i)=2.0d0*( q0(i)*opx-q3(i)*opy+q2(i)*opz)
      p2(i)=2.0d0*( q3(i)*opx+q0(i)*opy-q1(i)*opz)
      p3(i)=2.0d0*(-q2(i)*opx+q1(i)*opy+q0(i)*opz)
   end do
!$omp end do
!$omp barrier

!$omp do schedule(dynamic)
   do i=1,nm
      
      ! current rotation matrix

      call quat_to_rot(q0(i),q1(i),q2(i),q3(i),rott)

       ! calculate torque in principal frame

      tqx=tx(i)*rott(1)+ty(i)*rott(4)+tz(i)*rott(7)
      tqy=tx(i)*rott(2)+ty(i)*rott(5)+tz(i)*rott(8)
      tqz=tx(i)*rott(3)+ty(i)*rott(6)+tz(i)*rott(9)

       ! calculate quaternion torques

      qt0=2.0d0*(-q1(i)*tqx-q2(i)*tqy-q3(i)*tqz)
      qt1=2.0d0*( q0(i)*tqx-q3(i)*tqy+q2(i)*tqz)
      qt2=2.0d0*( q3(i)*tqx+q0(i)*tqy-q1(i)*tqz)
      qt3=2.0d0*(-q2(i)*tqx+q1(i)*tqy+q0(i)*tqz)


      ! calculate quaternion momenta at start of time step
!      opx = jx(i)
!      opy = jy(i)
!      opz = jz(i)

!      p0(i)=2.0d0*(-q1(i)*opx-q2(i)*opy-q3(i)*opz)
!      p1(i)=2.0d0*( q0(i)*opx-q3(i)*opy+q2(i)*opz)
!      p2(i)=2.0d0*( q3(i)*opx+q0(i)*opy-q1(i)*opz)
!      p3(i)=2.0d0*(-q2(i)*opx+q1(i)*opy+q0(i)*opz)

       ! update quaternion momenta by 1/2 time step

      p0(i)=p0(i)+qt0*0.5d0*ttt
      p1(i)=p1(i)+qt1*0.5d0*ttt
      p2(i)=p2(i)+qt2*0.5d0*ttt
      p3(i)=p3(i)+qt3*0.5d0*ttt

      ! update centre of mass velocity by 1/2 time step

      vx(i) = vx(i) + 0.5d0*fx(i)*ttt/totm
      vy(i) = vy(i) + 0.5d0*fy(i)*ttt/totm
      vz(i) = vz(i) + 0.5d0*fz(i)*ttt/totm

   end do  
!$omp end do
!$omp barrier
    ! move centre of mass by full time step
 
!$omp do schedule(dynamic)
   do i=1,nm
      x(i)= x(i) + ttt * vx(i)
      y(i)= y(i) + ttt * vy(i)
      z(i)= z(i) + ttt * vz(i)

      ! periodic boundary conditions

!     x(i) = x(i) - box*nint(x(i)/box)
!     y(i) = y(i) - box*nint(y(i)/box)
!     z(i) = z(i) - box*nint(z(i)/box)
!   end do

      ! rotate rigid groups: nosquish algorithm -update q to full tstep and
      ! amend p and get new rotation matrix

!   do i=1,nm

      qq0 = q0(i)
      qq1 = q1(i)
      qq2 = q2(i)
      qq3 = q3(i)
      pp0 = p0(i)
      pp1 = p1(i)
      pp2 = p2(i)
      pp3 = p3(i)

      call nosquish(ttt,qq0,qq1,qq2,qq3,pp0,pp1,pp2,pp3) 

      q0(i) = qq0
      q1(i) = qq1
      q2(i) = qq2
      q3(i) = qq3
      p0(i) = pp0
      p1(i) = pp1
      p2(i) = pp2
      p3(i) = pp3

   end do
!$omp end do
!$omp barrier

!   call sites

!$omp do schedule(dynamic)
   do i=1,nm
!      qq(0) = q0(i)
!      qq(1) = q1(i)
!      qq(2) = q2(i)
!      qq(3) = q3(i)

!      call quat_to_rot(qq,rott)

       ! update COM velocities and angular momentum to 1/2 tstep

      opx=0.5d0*(-q1(i)*p0(i)+q0(i)*p1(i)+ &
                 q3(i)*p2(i)-q2(i)*p3(i))
      opy=0.5d0*(-q2(i)*p0(i)-q3(i)*p1(i)+ &
                 q0(i)*p2(i)+q1(i)*p3(i))
      opz=0.5d0*(-q3(i)*p0(i)+q2(i)*p1(i)- &
                 q1(i)*p2(i)+q0(i)*p3(i))

      jx(i) = opx
      jy(i) = opy
      jz(i) = opz
   end do
!$omp end do
!$omp end parallel

!      vx(i) = vx(i) + 0.5d0*fx(i)*tt/totm
!      vy(i) = vy(i) + 0.5d0*fy(i)*tt/totm
!      vz(i) = vz(i) + 0.5d0*fz(i)*tt/totm

     ! first stage of velocity verlet algorithm

     ! move centre of mass by full time step

!      x(i)= x(i) + tt * vx(i)
!      y(i)= y(i) + tt * vy(i)
!      z(i)= z(i) + tt * vz(i)

      ! periodic boundary conditions

!    do i=1,nm

!      x(i) = x(i) - box*nint(x(i)/box)
!      y(i) = y(i) - box*nint(y(i)/box)
!      z(i) = z(i) - box*nint(z(i)/box)

!   end do

! end of first stage of velocity verlet algorithm

   call sites

   lvdwpe = vdwpe

	! DBG
	iniTime = SECNDS(0.0);
	! DBG
   if (pot_name.eq.'tip4p' ) then 
       call lj(x,y,z,xs,ys,zs,box,lvdwpe)
   else if (pot_name.eq.'tip5p' ) then
       call lj(x,y,z,xs,ys,zs,box,lvdwpe)
   else if (pot_name.eq.'sapt5s' ) then
       call sapt5s(x,y,z,xs,ys,zs,box,lvdwpe)
   else if (pot_name.eq.'ccpol5s' ) then
       call ccpol5s(x,y,z,xs,ys,zs,box,lvdwpe)
   else if (pot_name.eq.'ccpol8s' ) then
       call ccpol8s(x,y,z,xs,ys,zs,box,lvdwpe)       
   else if (pot_name.eq.'ccpol8s_omo' ) then
       call ccpol8s_omo(x,y,z,xs,ys,zs,box,lvdwpe)       
   else if (pot_name.eq.'ccpol23plus' ) then
       call ccpol23plus(x,y,z,xs,ys,zs,box,lvdwpe)
		! DBG
		deltaTime = SECNDS(iniTime);
!		write(dbgUnit, *) 'P#', MPIrank, '. ccpol23plus time	', deltaTime;
		! DBG
   end if

   vdwpe = lvdwpe

   lscpe = scpe
   call ewald_self_correction(alpha,box,lscpe)
   scpe = lscpe

   lrcpe = rcpe
   call rspace_ewald(x,y,z,xs,ys,zs,box,alpha,lrcpe)
   rcpe = lrcpe

   lkcpe = kcpe
   call kspace_ewald(x,y,z,xs,ys,zs,kmax1,kmax2,kmax3,alpha,box,lkcpe)
   kcpe = lkcpe

!  calling N-Body Iterative Induction model
   if (induction) then
      if(ind_model .eq. 'point_dipole') then

        lindpe = indpe
        call nb_induction(x,y,z,xs,ys,zs,box,lindpe)
        indpe = lindpe

      end if
   end if
 
! calling 3B-potenital (FS2 and FS3 terms)

  if (nonadditive_3B) then

     lfs2pe = fs2pe
     lfs3pe = fs3pe
     call threebody_term(x,y,z,xs,ys,zs,box,lfs2pe,lfs3pe)
     fs2pe = lfs2pe
     fs3pe = lfs3pe

  end if

! second stage of velocity verlet algorithm

!$omp parallel DEFAULT(SHARED) private(id,i,j,is,js,lvir)&
!$omp& private(qfx,qfy,qfz,qtx,qty,qtz, tqx,tqy,tqz)&
!$omp& private(qvx,qvy,qvz,qx,qy,qz,qjj,omg,qq,rott)&
!$omp& private(qt0,qt1,qt2,qt3,opx,opy,opz,cconorm,qmass)&
!$omp& private(hstep,qtke,qrke,qsystmp,dof,chit,conint,consv)&
!$omp& private(qq0,qq1,qq2,qq3,pp0,pp1,pp2,pp3)

   lvir = 0.0D0

!$omp do schedule(dynamic)
   do j=1,nm

      qfx=0.0
      qfy=0.0
      qfz=0.0
      qtx=0.0
      qty=0.0
      qtz=0.0

      do js=1,ns
         qfx=qfx+fxs(js,j)
         qfy=qfy+fys(js,j)
         qfz=qfz+fzs(js,j)
         qtx=qtx+ys(js,j)*fzs(js,j)-zs(js,j)*fys(js,j)
         qty=qty+zs(js,j)*fxs(js,j)-xs(js,j)*fzs(js,j)
         qtz=qtz+xs(js,j)*fys(js,j)-ys(js,j)*fxs(js,j)

         ! CONVERT FROM SITE-SITE TO MOLECULE-MOLECULE VIRIAL

         lvir = lvir + fxs(js,j)*xs(js,j)+fys(js,j)*ys(js,j)+fzs(js,j)*zs(js,j)

      end do

      fx(j)=qfx
      fy(j)=qfy
      fz(j)=qfz
      tx(j)=qtx
      ty(j)=qty
      tz(j)=qtz

      if (nonadditive_3B) then

         fx(j) = fx(j) + fxfs2(j) + fxfs3(j)
         fy(j) = fy(j) + fyfs2(j) + fyfs3(j)
         fz(j) = fz(j) + fzfs2(j) + fzfs3(j)

         tx(j)= tx(j) + txfs2(j) + txfs3(j)
         ty(j)= ty(j) + tyfs2(j) + tyfs3(j)
         tz(j)= tz(j) + tzfs2(j) + tzfs3(j)
      end if

   end do
!$omp end do
!$omp atomic

   vircom = vircom + lvir

!$omp barrier

!$omp do schedule(dynamic)
   do i=1,nm

      ! calculate quaternion momenta at start of time step
      opx = jx(i)
      opy = jy(i)
      opz = jz(i)

      p0(i)=2.0d0*(-q1(i)*opx-q2(i)*opy-q3(i)*opz)
      p1(i)=2.0d0*( q0(i)*opx-q3(i)*opy+q2(i)*opz)
      p2(i)=2.0d0*( q3(i)*opx+q0(i)*opy-q1(i)*opz)
      p3(i)=2.0d0*(-q2(i)*opx+q1(i)*opy+q0(i)*opz)
   end do
!$omp end do
!$omp barrier

!$omp do schedule(dynamic)
   do i=1,nm

       ! current rotation matrix

      call quat_to_rot(q0(i),q1(i),q2(i),q3(i),rott)

       ! calculate torque in principal frame

      tqx=tx(i)*rott(1)+ty(i)*rott(4)+tz(i)*rott(7)
      tqy=tx(i)*rott(2)+ty(i)*rott(5)+tz(i)*rott(8)
      tqz=tx(i)*rott(3)+ty(i)*rott(6)+tz(i)*rott(9)

       ! calculate quaternion torques

      qt0=2.0d0*(-q1(i)*tqx-q2(i)*tqy-q3(i)*tqz)
      qt1=2.0d0*( q0(i)*tqx-q3(i)*tqy+q2(i)*tqz)
      qt2=2.0d0*( q3(i)*tqx+q0(i)*tqy-q1(i)*tqz)
      qt3=2.0d0*(-q2(i)*tqx+q1(i)*tqy+q0(i)*tqz)

      ! calculate quaternion momenta at start of time step
!      opx = jx(i)
!      opy = jy(i)
!      opz = jz(i)

!      p0(i)=2.0d0*(-q1(i)*opx-q2(i)*opy-q3(i)*opz)
!      p1(i)=2.0d0*( q0(i)*opx-q3(i)*opy+q2(i)*opz)
!      p2(i)=2.0d0*( q3(i)*opx+q0(i)*opy-q1(i)*opz)
!      p3(i)=2.0d0*(-q2(i)*opx+q1(i)*opy+q0(i)*opz)

       ! update quaternion momenta by 1/2 time step

      p0(i)=p0(i)+qt0*0.5d0*ttt
      p1(i)=p1(i)+qt1*0.5d0*ttt
      p2(i)=p2(i)+qt2*0.5d0*ttt
      p3(i)=p3(i)+qt3*0.5d0*ttt

      !  update Rgid body angular momentum & COM velocities to full step

!      opx=0.5d0*(-q1(i)*p0(i)+q0(i)*p1(i)+ &
!                 q3(i)*p2(i)-q2(i)*p3(i))
!      opy=0.5d0*(-q2(i)*p0(i)-q3(i)*p1(i)+ &
!                 q0(i)*p2(i)+q1(i)*p3(i))
!      opz=0.5d0*(-q3(i)*p0(i)+q2(i)*p1(i)- &
!                 q1(i)*p2(i)+q0(i)*p3(i))

!      jx(i) = opx
!      jy(i) = opy
!      jz(i) = opz

      vx(i) = vx(i) + 0.5d0*fx(i)*ttt/totm
      vy(i) = vy(i) + 0.5d0*fy(i)*ttt/totm
      vz(i) = vz(i) + 0.5d0*fz(i)*ttt/totm

   end do
!$omp end do
!$omp barrier

!$omp do schedule(dynamic)
   do i=1,nm
!      qq(0) = q0(i)
!      qq(1) = q1(i)
!      qq(2) = q2(i)
!      qq(3) = q3(i)

!      call quat_to_rot(qq,rott)

       ! update COM velocities and angular momentum to 1/2 tstep

      opx=0.5d0*(-q1(i)*p0(i)+q0(i)*p1(i)+ &
                 q3(i)*p2(i)-q2(i)*p3(i))
      opy=0.5d0*(-q2(i)*p0(i)-q3(i)*p1(i)+ &
                 q0(i)*p2(i)+q1(i)*p3(i))
      opz=0.5d0*(-q3(i)*p0(i)+q2(i)*p1(i)- &
                 q1(i)*p2(i)+q0(i)*p3(i))

      jx(i) = opx
      jy(i) = opy
      jz(i) = opz
   end do
!$omp end do
!$omp end parallel

! apply thermostat for second stage and calculate kinetic energy

   call nvtqscl(hstep,qmass,dof,taut,chit,conint)

! conserved quantity less kinetic and potential energy terms

   consv=conint+0.5d0*qmass*chit**2

   ltke = tke
   call tkinetic_eng (nm, totm, vx, vy, vz, ltke)
   tke = ltke
   lttemp = ttemp
   call trans_temp (nm, tke, lttemp)
   ttemp = lttemp
   stptke = tke/418.4
   stpttp = ttemp

   lrke = rke
   call rkinetic_eng (nm, jx, jy, jz, lrke)
   rke = lrke
   stprke = rke/418.4

   lrtemp = rtemp
   call rot_temp (nm, rke, lrtemp)
   rtemp = lrtemp
   stprtp = rtemp

   lstptmp = stptmp
   call system_temp(nm,tke,rke,lstptmp)
   stptmp = lstptmp

   qvir = vir+vircom

!  calculate system pressure

   pres = (2.d0*ltke - qvir)/(3.d0*volm)

!   write(*,*) box,box**3.d0,volm,pres,vir,ltke

   stpvolm = volm

   stpvir = qvir/418.4d0
   stpprs = prsunt*pres

   deallocate(p0)
   deallocate(p1)
   deallocate(p2)
   deallocate(p3)
	! DBG
	deltaTimeTot = SECNDS(iniTimeTot);
!	write(dbgUnit, *) 'P#', MPIrank, '. nvt_hov time	', deltaTimeTot;
	call MPI_BARRIER(MPI_COMM_WORLD, MPIierr);
!	STOP;
	! DBG

  end subroutine nvt_hov  


  subroutine nvt_hov_2(ttt)
   implicit none
   real(8) :: ttt

   integer :: i,j,is,js
   real(8) :: qfx,qfy,qfz
   real(8) :: qtx,qty,qtz, tqx,tqy,tqz
   real(8) :: qvx,qvy,qvz
   real(8) :: qx,qy,qz
   real(8) :: qjj(3),omg(3)
   real(8) :: qq(0:3), rott(1:9)
   real(8) :: qt0,qt1,qt2,qt3
   real(8) :: opx,opy,opz
   real(8) :: cconorm,lvir,qvir
   real(8) :: lvdwpe,lscpe,lrcpe,lkcpe,lindpe,lfs2pe,lfs3pe
   real(8) :: ltke,lttemp,lrke,lrtemp,lstptmp
   real(8) :: qq0,qq1,qq2,qq3,pp0,pp1,pp2,pp3
!   real(8) :: q0t(nom),q1t(nom),q2t(nom),q3t(nom)
!   real(8) :: p0(nom),p1(nom),p2(nom),p3(nom)
    real(8), allocatable :: p0(:),p1(:),p2(:),p3(:)

! Nose Hoover parameters
   real(8),parameter :: taut = 0.5d0 ! ps (time constant)
   real(8) :: qmass,hstep,qtke,qrke,qsystmp,dof,chit,conint
   real(8) :: consv  
        ! DBG
        real(4) iniTime, deltaTime;
        real(4) iniTimeTot, deltaTimeTot;
        ! DBG

        ! DBG
        iniTimeTot = SECNDS(0.0);
        ! DBG

   allocate(p0(nm))
   allocate(p1(nm))
   allocate(p2(nm))
   allocate(p3(nm))

   ! GRU
   chit = 0.0D0
   conint = 0.0D0
   ! GRU

!$omp parallel DEFAULT(SHARED) private(id,i,j,is,js,lvir)&
!$omp& private(qfx,qfy,qfz,qtx,qty,qtz, tqx,tqy,tqz)&
!$omp& private(qvx,qvy,qvz,qx,qy,qz,qjj,omg,qq,rott)&
!$omp& private(qt0,qt1,qt2,qt3,opx,opy,opz,cconorm,qmass)&
!$omp& private(hstep,qtke,qrke,qsystmp,dof,chit,conint,consv)&
!$omp& private(qq0,qq1,qq2,qq3,pp0,pp1,pp2,pp3)

   lvir = 0.0D0

!$omp do schedule(dynamic)

   do j=1,nm

      qfx=0.0
      qfy=0.0
      qfz=0.0
      qtx=0.0
      qty=0.0
      qtz=0.0

      do js=1,ns
         qfx=qfx+fxs(js,j)
         qfy=qfy+fys(js,j)
         qfz=qfz+fzs(js,j)
         qtx=qtx+ys(js,j)*fzs(js,j)-zs(js,j)*fys(js,j)
         qty=qty+zs(js,j)*fxs(js,j)-xs(js,j)*fzs(js,j)
         qtz=qtz+xs(js,j)*fys(js,j)-ys(js,j)*fxs(js,j)

         ! CONVERT FROM SITE-SITE TO MOLECULE-MOLECULE VIRIAL

         lvir = lvir + fxs(js,j)*xs(js,j)+fys(js,j)*ys(js,j)+fzs(js,j)*zs(js,j)

      end do

      fx(j)=qfx
      fy(j)=qfy
      fz(j)=qfz
      tx(j)=qtx
      ty(j)=qty
      tz(j)=qtz

      if (nonadditive_3B) then

         fx(j) = fx(j) + fxfs2(j) + fxfs3(j)
         fy(j) = fy(j) + fyfs2(j) + fyfs3(j)
         fz(j) = fz(j) + fzfs2(j) + fzfs3(j)

         tx(j)= tx(j) + txfs2(j) + txfs3(j)
         ty(j)= ty(j) + tyfs2(j) + tyfs3(j)
         tz(j)= tz(j) + tzfs2(j) + tzfs3(j)


      end if

   end do
!$omp end do
!$omp atomic

   vircom = vircom + lvir

!$omp barrier
!$omp end parallel

!!$omp master
!  timestep parameters  
   hstep = 0.5d0*ttt

   dof = 3.d0*(2.d0*nm-1.d0) ! total dof transl + rot 
!   dof = 3.d0*(nm-1.d0) ! total dof transl + rot
!   nose-hoover inertia parameter
   qmass=dof*boltz*temp0*taut**2

   call nvtqscl(hstep,qmass,dof,taut,chit,conint)

!!$omp end master
!!$omp barrier

!$omp parallel DEFAULT(SHARED) private(id,i,j,is,js,lvir)&
!$omp& private(qfx,qfy,qfz,qtx,qty,qtz, tqx,tqy,tqz)&
!$omp& private(qvx,qvy,qvz,qx,qy,qz,qjj,omg,qq,rott)&
!$omp& private(qt0,qt1,qt2,qt3,opx,opy,opz,cconorm,qmass)&
!$omp& private(hstep,qtke,qrke,qsystmp,dof,chit,conint,consv)&
!$omp& private(qq0,qq1,qq2,qq3,pp0,pp1,pp2,pp3)

!   *************  Rigid body motion ****************************
!     operations common to first and second stages

!$omp do schedule(dynamic)
   do i=1,nm

      ! calculate quaternion momenta at start of time step
      opx = jx(i)
      opy = jy(i)
      opz = jz(i)

      p0(i)=2.0d0*(-q1(i)*opx-q2(i)*opy-q3(i)*opz)
      p1(i)=2.0d0*( q0(i)*opx-q3(i)*opy+q2(i)*opz)
      p2(i)=2.0d0*( q3(i)*opx+q0(i)*opy-q1(i)*opz)
      p3(i)=2.0d0*(-q2(i)*opx+q1(i)*opy+q0(i)*opz)
   end do
!$omp end do
!$omp barrier

!$omp do schedule(dynamic)
   do i=1,nm
      
      ! current rotation matrix

      call quat_to_rot(q0(i),q1(i),q2(i),q3(i),rott)

       ! calculate torque in principal frame

      tqx=tx(i)*rott(1)+ty(i)*rott(4)+tz(i)*rott(7)
      tqy=tx(i)*rott(2)+ty(i)*rott(5)+tz(i)*rott(8)
      tqz=tx(i)*rott(3)+ty(i)*rott(6)+tz(i)*rott(9)

       ! calculate quaternion torques

      qt0=2.0d0*(-q1(i)*tqx-q2(i)*tqy-q3(i)*tqz)
      qt1=2.0d0*( q0(i)*tqx-q3(i)*tqy+q2(i)*tqz)
      qt2=2.0d0*( q3(i)*tqx+q0(i)*tqy-q1(i)*tqz)
      qt3=2.0d0*(-q2(i)*tqx+q1(i)*tqy+q0(i)*tqz)


      ! calculate quaternion momenta at start of time step
!      opx = jx(i)
!      opy = jy(i)
!      opz = jz(i)

!      p0(i)=2.0d0*(-q1(i)*opx-q2(i)*opy-q3(i)*opz)
!      p1(i)=2.0d0*( q0(i)*opx-q3(i)*opy+q2(i)*opz)
!      p2(i)=2.0d0*( q3(i)*opx+q0(i)*opy-q1(i)*opz)
!      p3(i)=2.0d0*(-q2(i)*opx+q1(i)*opy+q0(i)*opz)

       ! update quaternion momenta by 1/2 time step

      p0(i)=p0(i)+qt0*0.5d0*ttt
      p1(i)=p1(i)+qt1*0.5d0*ttt
      p2(i)=p2(i)+qt2*0.5d0*ttt
      p3(i)=p3(i)+qt3*0.5d0*ttt

      ! update centre of mass velocity by 1/2 time step

      vx(i) = vx(i) + 0.5d0*fx(i)*ttt/totm
      vy(i) = vy(i) + 0.5d0*fy(i)*ttt/totm
      vz(i) = vz(i) + 0.5d0*fz(i)*ttt/totm

   end do  
!$omp end do
!$omp barrier
    ! move centre of mass by full time step
 
!$omp do schedule(dynamic)
   do i=1,nm
      x(i)= x(i) + ttt * vx(i)
      y(i)= y(i) + ttt * vy(i)
      z(i)= z(i) + ttt * vz(i)

      ! periodic boundary conditions

!     x(i) = x(i) - box*nint(x(i)/box)
!     y(i) = y(i) - box*nint(y(i)/box)
!     z(i) = z(i) - box*nint(z(i)/box)
!   end do

      ! rotate rigid groups: nosquish algorithm -update q to full tstep and
      ! amend p and get new rotation matrix

!   do i=1,nm

      qq0 = q0(i)
      qq1 = q1(i)
      qq2 = q2(i)
      qq3 = q3(i)
      pp0 = p0(i)
      pp1 = p1(i)
      pp2 = p2(i)
      pp3 = p3(i)

      call nosquish(ttt,qq0,qq1,qq2,qq3,pp0,pp1,pp2,pp3) 

      q0(i) = qq0
      q1(i) = qq1
      q2(i) = qq2
      q3(i) = qq3
      p0(i) = pp0
      p1(i) = pp1
      p2(i) = pp2
      p3(i) = pp3

   end do
!$omp end do
!$omp barrier

!   call sites

!$omp do schedule(dynamic)
!   do i=1,nm
!      qq(0) = q0(i)
!      qq(1) = q1(i)
!      qq(2) = q2(i)
!      qq(3) = q3(i)

!      call quat_to_rot(qq,rott)

       ! update COM velocities and angular momentum to 1/2 tstep

!      opx=0.5d0*(-q1(i)*p0(i)+q0(i)*p1(i)+ &
!                 q3(i)*p2(i)-q2(i)*p3(i))
!      opy=0.5d0*(-q2(i)*p0(i)-q3(i)*p1(i)+ &
!                 q0(i)*p2(i)+q1(i)*p3(i))
!      opz=0.5d0*(-q3(i)*p0(i)+q2(i)*p1(i)- &
!                 q1(i)*p2(i)+q0(i)*p3(i))

!      jx(i) = opx
!      jy(i) = opy
!      jz(i) = opz
!   end do
!$omp end do
!$omp end parallel

!      vx(i) = vx(i) + 0.5d0*fx(i)*tt/totm
!      vy(i) = vy(i) + 0.5d0*fy(i)*tt/totm
!      vz(i) = vz(i) + 0.5d0*fz(i)*tt/totm

     ! first stage of velocity verlet algorithm

     ! move centre of mass by full time step

!      x(i)= x(i) + tt * vx(i)
!      y(i)= y(i) + tt * vy(i)
!      z(i)= z(i) + tt * vz(i)

      ! periodic boundary conditions

!    do i=1,nm

!      x(i) = x(i) - box*nint(x(i)/box)
!      y(i) = y(i) - box*nint(y(i)/box)
!      z(i) = z(i) - box*nint(z(i)/box)

!   end do

! end of first stage of velocity verlet algorithm

   call sites

   lvdwpe = vdwpe

        ! DBG
        iniTime = SECNDS(0.0);
        ! DBG
   if (pot_name.eq.'tip4p' ) then 
       call lj(x,y,z,xs,ys,zs,box,lvdwpe)
   else if (pot_name.eq.'tip5p' ) then
       call lj(x,y,z,xs,ys,zs,box,lvdwpe)
   else if (pot_name.eq.'sapt5s' ) then
       call sapt5s(x,y,z,xs,ys,zs,box,lvdwpe)
   else if (pot_name.eq.'ccpol5s' ) then
       call ccpol5s(x,y,z,xs,ys,zs,box,lvdwpe)
   else if (pot_name.eq.'ccpol8s' ) then
       call ccpol8s(x,y,z,xs,ys,zs,box,lvdwpe)       
   else if (pot_name.eq.'ccpol8s_omo' ) then
       call ccpol8s_omo(x,y,z,xs,ys,zs,box,lvdwpe)       
   else if (pot_name.eq.'ccpol23plus' ) then
       call ccpol23plus(x,y,z,xs,ys,zs,box,lvdwpe)
                ! DBG
                deltaTime = SECNDS(iniTime);
!               write(dbgUnit, *) 'P#', MPIrank, '. ccpol23plus time    ',
!               deltaTime;
                ! DBG
   end if

   vdwpe = lvdwpe

   lscpe = scpe
   call ewald_self_correction(alpha,box,lscpe)
   scpe = lscpe

   lrcpe = rcpe
   call rspace_ewald(x,y,z,xs,ys,zs,box,alpha,lrcpe)
   rcpe = lrcpe

   lkcpe = kcpe
   call kspace_ewald(x,y,z,xs,ys,zs,kmax1,kmax2,kmax3,alpha,box,lkcpe)
   kcpe = lkcpe

!  calling N-Body Iterative Induction model
   if (induction) then
      if(ind_model .eq. 'point_dipole') then

        lindpe = indpe
        call nb_induction(x,y,z,xs,ys,zs,box,lindpe)
        indpe = lindpe

      end if
   end if
 
! calling 3B-potenital (FS2 and FS3 terms)

  if (nonadditive_3B) then

     lfs2pe = fs2pe
     lfs3pe = fs3pe
     call threebody_term(x,y,z,xs,ys,zs,box,lfs2pe,lfs3pe)
     fs2pe = lfs2pe
     fs3pe = lfs3pe

  end if

! second stage of velocity verlet algorithm

!$omp parallel DEFAULT(SHARED) private(id,i,j,is,js,lvir)&
!$omp& private(qfx,qfy,qfz,qtx,qty,qtz, tqx,tqy,tqz)&
!$omp& private(qvx,qvy,qvz,qx,qy,qz,qjj,omg,qq,rott)&
!$omp& private(qt0,qt1,qt2,qt3,opx,opy,opz,cconorm,qmass)&
!$omp& private(hstep,qtke,qrke,qsystmp,dof,chit,conint,consv)&
!$omp& private(qq0,qq1,qq2,qq3,pp0,pp1,pp2,pp3)

   lvir = 0.0D0

!$omp do schedule(dynamic)
   do j=1,nm

      qfx=0.0
      qfy=0.0
      qfz=0.0
      qtx=0.0
      qty=0.0
      qtz=0.0

      do js=1,ns
         qfx=qfx+fxs(js,j)
         qfy=qfy+fys(js,j)
         qfz=qfz+fzs(js,j)
         qtx=qtx+ys(js,j)*fzs(js,j)-zs(js,j)*fys(js,j)
         qty=qty+zs(js,j)*fxs(js,j)-xs(js,j)*fzs(js,j)
         qtz=qtz+xs(js,j)*fys(js,j)-ys(js,j)*fxs(js,j)

         ! CONVERT FROM SITE-SITE TO MOLECULE-MOLECULE VIRIAL

         lvir = lvir + fxs(js,j)*xs(js,j)+fys(js,j)*ys(js,j)+fzs(js,j)*zs(js,j)

      end do

      fx(j)=qfx
      fy(j)=qfy
      fz(j)=qfz
      tx(j)=qtx
      ty(j)=qty
      tz(j)=qtz

      if (nonadditive_3B) then

         fx(j) = fx(j) + fxfs2(j) + fxfs3(j)
         fy(j) = fy(j) + fyfs2(j) + fyfs3(j)
         fz(j) = fz(j) + fzfs2(j) + fzfs3(j)

         tx(j)= tx(j) + txfs2(j) + txfs3(j)
         ty(j)= ty(j) + tyfs2(j) + tyfs3(j)
         tz(j)= tz(j) + tzfs2(j) + tzfs3(j)
      end if

   end do
!$omp end do
!$omp atomic

   vircom = vircom + lvir

!$omp barrier

!$omp do schedule(dynamic)
!   do i=1,nm

      ! calculate quaternion momenta at start of time step
!      opx = jx(i)
!      opy = jy(i)
!      opz = jz(i)

!      p0(i)=2.0d0*(-q1(i)*opx-q2(i)*opy-q3(i)*opz)
!      p1(i)=2.0d0*( q0(i)*opx-q3(i)*opy+q2(i)*opz)
!      p2(i)=2.0d0*( q3(i)*opx+q0(i)*opy-q1(i)*opz)
!      p3(i)=2.0d0*(-q2(i)*opx+q1(i)*opy+q0(i)*opz)
!   end do
!$omp end do
!$omp barrier

!$omp do schedule(dynamic)
   do i=1,nm

       ! current rotation matrix

      call quat_to_rot(q0(i),q1(i),q2(i),q3(i),rott)

       ! calculate torque in principal frame

      tqx=tx(i)*rott(1)+ty(i)*rott(4)+tz(i)*rott(7)
      tqy=tx(i)*rott(2)+ty(i)*rott(5)+tz(i)*rott(8)
      tqz=tx(i)*rott(3)+ty(i)*rott(6)+tz(i)*rott(9)

       ! calculate quaternion torques

      qt0=2.0d0*(-q1(i)*tqx-q2(i)*tqy-q3(i)*tqz)
      qt1=2.0d0*( q0(i)*tqx-q3(i)*tqy+q2(i)*tqz)
      qt2=2.0d0*( q3(i)*tqx+q0(i)*tqy-q1(i)*tqz)
      qt3=2.0d0*(-q2(i)*tqx+q1(i)*tqy+q0(i)*tqz)

      ! calculate quaternion momenta at start of time step
!      opx = jx(i)
!      opy = jy(i)
!      opz = jz(i)

!      p0(i)=2.0d0*(-q1(i)*opx-q2(i)*opy-q3(i)*opz)
!      p1(i)=2.0d0*( q0(i)*opx-q3(i)*opy+q2(i)*opz)
!      p2(i)=2.0d0*( q3(i)*opx+q0(i)*opy-q1(i)*opz)
!      p3(i)=2.0d0*(-q2(i)*opx+q1(i)*opy+q0(i)*opz)

       ! update quaternion momenta by 1/2 time step

      p0(i)=p0(i)+qt0*0.5d0*ttt
      p1(i)=p1(i)+qt1*0.5d0*ttt
      p2(i)=p2(i)+qt2*0.5d0*ttt
      p3(i)=p3(i)+qt3*0.5d0*ttt

      !  update Rgid body angular momentum & COM velocities to full step

!      opx=0.5d0*(-q1(i)*p0(i)+q0(i)*p1(i)+ &
!                 q3(i)*p2(i)-q2(i)*p3(i))
!      opy=0.5d0*(-q2(i)*p0(i)-q3(i)*p1(i)+ &
!                 q0(i)*p2(i)+q1(i)*p3(i))
!      opz=0.5d0*(-q3(i)*p0(i)+q2(i)*p1(i)- &
!                 q1(i)*p2(i)+q0(i)*p3(i))

!      jx(i) = opx
!      jy(i) = opy
!      jz(i) = opz

      vx(i) = vx(i) + 0.5d0*fx(i)*ttt/totm
      vy(i) = vy(i) + 0.5d0*fy(i)*ttt/totm
      vz(i) = vz(i) + 0.5d0*fz(i)*ttt/totm

   end do
!$omp end do
!$omp barrier

!$omp do schedule(dynamic)
   do i=1,nm
!      qq(0) = q0(i)
!      qq(1) = q1(i)
!      qq(2) = q2(i)
!      qq(3) = q3(i)

!      call quat_to_rot(qq,rott)

       ! update COM velocities and angular momentum to 1/2 tstep

      opx=0.5d0*(-q1(i)*p0(i)+q0(i)*p1(i)+ &
                 q3(i)*p2(i)-q2(i)*p3(i))
      opy=0.5d0*(-q2(i)*p0(i)-q3(i)*p1(i)+ &
                 q0(i)*p2(i)+q1(i)*p3(i))
      opz=0.5d0*(-q3(i)*p0(i)+q2(i)*p1(i)- &
                 q1(i)*p2(i)+q0(i)*p3(i))

      jx(i) = opx
      jy(i) = opy
      jz(i) = opz
   end do
!$omp end do
!$omp end parallel

! apply thermostat for second stage and calculate kinetic energy

   call nvtqscl(hstep,qmass,dof,taut,chit,conint)

! conserved quantity less kinetic and potential energy terms

   consv=conint+0.5d0*qmass*chit**2

   ltke = tke
   call tkinetic_eng (nm, totm, vx, vy, vz, ltke)
   tke = ltke
   lttemp = ttemp
   call trans_temp (nm, tke, lttemp)
   ttemp = lttemp
   stptke = tke/418.4
   stpttp = ttemp

   lrke = rke
   call rkinetic_eng (nm, jx, jy, jz, lrke)
   rke = lrke
   stprke = rke/418.4

   lrtemp = rtemp
   call rot_temp (nm, rke, lrtemp)
   rtemp = lrtemp
   stprtp = rtemp

   lstptmp = stptmp
   call system_temp(nm,tke,rke,lstptmp)
   stptmp = lstptmp

!   write(*,*) temp0,stptmp

   qvir = vir+vircom

!  calculate system pressure

   pres = (2.d0*ltke - qvir)/(3.d0*volm)

!   write(*,*) box,box**3.d0,volm,pres,vir,ltke

   stpvolm = volm

   stpvir = qvir/418.4d0
   stpprs = prsunt*pres

   deallocate(p0)
   deallocate(p1)
   deallocate(p2)
   deallocate(p3)
        ! DBG
        deltaTimeTot = SECNDS(iniTimeTot);
!       write(dbgUnit, *) 'P#', MPIrank, '. nvt_hov time        ', deltaTimeTot;
        call MPI_BARRIER(MPI_COMM_WORLD, MPIierr);
!       STOP;
        ! DBG

  end subroutine nvt_hov_2


  subroutine npt_hov(ttt)
   implicit none
   real(8) :: ttt

   integer :: i,j,is,js,icyc
   real(8) :: qfx,qfy,qfz
   real(8) :: qtx,qty,qtz, tqx,tqy,tqz
   real(8) :: qvx,qvy,qvz
   real(8) :: qx,qy,qz
   real(8) :: qjj(3),omg(3)
   real(8) :: qq(0:3), rott(1:9)
   real(8) :: qt0,qt1,qt2,qt3
   real(8) :: opx,opy,opz
   real(8) :: cconorm,lvir,qvir,qpres
   real(8) :: lvdwpe,lscpe,lrcpe,lkcpe,lindpe,lfs2pe,lfs3pe
   real(8) :: ltke,lttemp,lrke,lrtemp,lstptmp,scale
   real(8) :: qq0,qq1,qq2,qq3,pp0,pp1,pp2,pp3
!   real(8) :: q0t(nom),q1t(nom),q2t(nom),q3t(nom)
!   real(8) :: p0(nom),p1(nom),p2(nom),p3(nom)
    real(8), allocatable :: p0(:),p1(:),p2(:),p3(:)

    real(8), allocatable :: p0t(:),p1t(:),p2t(:),p3t(:)
    real(8), allocatable :: q0tt(:),q1tt(:),q2tt(:),q3tt(:)
    real(8), allocatable :: x0t(:),y0t(:),z0t(:) 
    real(8), allocatable :: vx0t(:),vy0t(:),vz0t(:)
    real(8), allocatable :: jx0t(:),jy0t(:),jz0t(:)

! Nose Hoover parameters
   real(8),parameter :: taut = 0.05d0 ! ps (time constant)
   real(8),parameter :: taup = 0.05d0 ! ps (time constant check!!)
   integer, parameter :: ncyc = 5 
   real(8) :: qmass,hstep,qtke,qrke,qsystmp,dof,chit,conint
   real(8) :: pmass,consv,chip,chitt,chipt,conintt,qstep,fstep  
	! DBG
	real(4) iniTime, deltaTime;
	real(4) iniTimeTot, deltaTimeTot;
	! DBG

	! DBG
	iniTimeTot = SECNDS(0.0);
	! DBG

   allocate(p0(nm),p1(nm),p2(nm),p3(nm))
   allocate(q0tt(nm),q1tt(nm),q2tt(nm),q3tt(nm))
   allocate(x0t(nm),y0t(nm),z0t(nm))
   allocate(vx0t(nm),vy0t(nm),vz0t(nm))
   allocate(jx0t(nm),jy0t(nm),jz0t(nm))

   ! GRU
   chip = 0.0D0
   chit = 0.0D0
   conint = 0.0D0
   ! GRU

!$omp parallel DEFAULT(SHARED) private(id,i,j,is,js,lvir)&
!$omp& private(qfx,qfy,qfz,qtx,qty,qtz, tqx,tqy,tqz)&
!$omp& private(qvx,qvy,qvz,qx,qy,qz,qjj,omg,qq,rott)&
!$omp& private(qt0,qt1,qt2,qt3,opx,opy,opz,cconorm,qmass,pmass)&
!$omp& private(hstep,qtke,qrke,qsystmp,dof,chit,conint,consv,boxt,volmt)&
!$omp& private(qq0,qq1,qq2,qq3,pp0,pp1,pp2,pp3,chip,chitt,chipt,conintt)

   lvir = 0.0D0

!$omp do schedule(dynamic)

   do j=1,nm

      qfx=0.0
      qfy=0.0
      qfz=0.0
      qtx=0.0
      qty=0.0
      qtz=0.0

      do js=1,ns
         qfx=qfx+fxs(js,j)
         qfy=qfy+fys(js,j)
         qfz=qfz+fzs(js,j)
         qtx=qtx+ys(js,j)*fzs(js,j)-zs(js,j)*fys(js,j)
         qty=qty+zs(js,j)*fxs(js,j)-xs(js,j)*fzs(js,j)
         qtz=qtz+xs(js,j)*fys(js,j)-ys(js,j)*fxs(js,j)

         ! CONVERT FROM SITE-SITE TO MOLECULE-MOLECULE VIRIAL

         lvir = lvir + fxs(js,j)*xs(js,j)+fys(js,j)*ys(js,j)+fzs(js,j)*zs(js,j)

      end do

      fx(j)=qfx
      fy(j)=qfy
      fz(j)=qfz
      tx(j)=qtx
      ty(j)=qty
      tz(j)=qtz

      if (nonadditive_3B) then

         fx(j) = fx(j) + fxfs2(j) + fxfs3(j)
         fy(j) = fy(j) + fyfs2(j) + fyfs3(j)
         fz(j) = fz(j) + fzfs2(j) + fzfs3(j)

         tx(j)= tx(j) + txfs2(j) + txfs3(j)
         ty(j)= ty(j) + tyfs2(j) + tyfs3(j)
         tz(j)= tz(j) + tzfs2(j) + tzfs3(j)


      end if

   end do
!$omp end do
!$omp atomic

   vircom = vircom + lvir

!$omp barrier
!$omp end parallel

!!$omp master
!  timestep parameters  
   hstep = 0.5d0*ttt
   qstep = 0.25d0*ttt/dble(ncyc)
   fstep = 0.5d0*ttt/dble(ncyc)

! temporary saved variables

   volm=volm0
   volmt=volm
   box=box0
   boxt = box
   chipt = chip
   chitt = chit
   conintt = conint
   qvir = vir

   do j=1,nm

      q0tt(j) = q0(j)
      q1tt(j) = q1(j)
      q2tt(j) = q2(j)
      q3tt(j) = q3(j)
      x0t(j)  = x(j)
      y0t(j)  = y(j)
      z0t(j)  = z(j)
      vx0t(j) = vx(j)
      vy0t(j) = vy(j)
      vz0t(j) = vz(j)
      jx0t(j) = jx(j)
      jy0t(j) = jy(j)
      jz0t(j) = jz(j) 

   end do 

   dof = 3.d0*(2.d0*nm-1.d0) ! total dof transl + rot 
!   nose-hoover inertia parameter
   qmass=dof*boltz*temp0*taut**2
   pmass=dof*boltz*temp0*taup**2

   do i=icyc,ncyc

!  npt thermostat
   call nptqscl_t(qstep,qmass,pmass,dof,taut,chip,chit,conint)
!  npt barostat
   call nptqscl_p(fstep,box,pmass,dof,chit,chip,qvir,qpres,conint) 
!  npt thermostat
   call nptqscl_t(qstep,qmass,pmass,dof,taut,chip,chit,conint)

   end do

!   scale cell vectors - isotropic

!    volm=box**3.d0  

    scale=(volm/volm0)**(1.d0/3.d0)

    box=box0*scale
!    volm=box**3.d0

!!$omp end master
!!$omp barrier

!$omp parallel DEFAULT(SHARED) private(id,i,j,is,js,lvir)&
!$omp& private(qfx,qfy,qfz,qtx,qty,qtz, tqx,tqy,tqz)&
!$omp& private(qvx,qvy,qvz,qx,qy,qz,qjj,omg,qq,rott)&
!$omp& private(qt0,qt1,qt2,qt3,opx,opy,opz,cconorm,qmass)&
!$omp& private(hstep,qtke,qrke,qsystmp,dof,chit,conint,consv,boxt,volmt)&
!$omp& private(qq0,qq1,qq2,qq3,pp0,pp1,pp2,pp3,chip,chitt,chipt,conintt)

!   *************  Rigid body motion ****************************
!     operations common to first and second stages

!$omp do schedule(dynamic)
   do i=1,nm

      ! calculate quaternion momenta at start of time step
      opx = jx(i)
      opy = jy(i)
      opz = jz(i)

      p0(i)=2.0d0*(-q1(i)*opx-q2(i)*opy-q3(i)*opz)
      p1(i)=2.0d0*( q0(i)*opx-q3(i)*opy+q2(i)*opz)
      p2(i)=2.0d0*( q3(i)*opx+q0(i)*opy-q1(i)*opz)
      p3(i)=2.0d0*(-q2(i)*opx+q1(i)*opy+q0(i)*opz)
   end do
!$omp end do
!$omp barrier

!$omp do schedule(dynamic)
   do i=1,nm
      
      ! current rotation matrix

      call quat_to_rot(q0(i),q1(i),q2(i),q3(i),rott)

       ! calculate torque in principal frame

      tqx=tx(i)*rott(1)+ty(i)*rott(4)+tz(i)*rott(7)
      tqy=tx(i)*rott(2)+ty(i)*rott(5)+tz(i)*rott(8)
      tqz=tx(i)*rott(3)+ty(i)*rott(6)+tz(i)*rott(9)

       ! calculate quaternion torques

      qt0=2.0d0*(-q1(i)*tqx-q2(i)*tqy-q3(i)*tqz)
      qt1=2.0d0*( q0(i)*tqx-q3(i)*tqy+q2(i)*tqz)
      qt2=2.0d0*( q3(i)*tqx+q0(i)*tqy-q1(i)*tqz)
      qt3=2.0d0*(-q2(i)*tqx+q1(i)*tqy+q0(i)*tqz)


      ! calculate quaternion momenta at start of time step
!      opx = jx(i)
!      opy = jy(i)
!      opz = jz(i)

!      p0(i)=2.0d0*(-q1(i)*opx-q2(i)*opy-q3(i)*opz)
!      p1(i)=2.0d0*( q0(i)*opx-q3(i)*opy+q2(i)*opz)
!      p2(i)=2.0d0*( q3(i)*opx+q0(i)*opy-q1(i)*opz)
!      p3(i)=2.0d0*(-q2(i)*opx+q1(i)*opy+q0(i)*opz)

       ! update quaternion momenta by 1/2 time step

      p0(i)=p0(i)+qt0*0.5d0*ttt
      p1(i)=p1(i)+qt1*0.5d0*ttt
      p2(i)=p2(i)+qt2*0.5d0*ttt
      p3(i)=p3(i)+qt3*0.5d0*ttt

      ! update centre of mass velocity by 1/2 time step

      vx(i) = vx(i) + 0.5d0*fx(i)*ttt/totm
      vy(i) = vy(i) + 0.5d0*fy(i)*ttt/totm
      vz(i) = vz(i) + 0.5d0*fz(i)*ttt/totm

   end do  
!$omp end do
!$omp barrier
    ! move centre of mass by full time step
 
!$omp do schedule(dynamic)
   do i=1,nm
      do is=1,ns
      x(i)= x(i) + ttt * (vx(i)+chip*xs(i,is))
      y(i)= y(i) + ttt * (vy(i)+chip*ys(i,is))
      z(i)= z(i) + ttt * (vz(i)+chip*zs(i,is))
      end do
 
      ! periodic boundary conditions

!!     x(i) = x(i) - box*nint(x(i)/box)
!!     y(i) = y(i) - box*nint(y(i)/box)
!!     z(i) = z(i) - box*nint(z(i)/box)
!   end do

      ! rotate rigid groups: nosquish algorithm -update q to full tstep and
      ! amend p and get new rotation matrix

!   do i=1,nm

      qq0 = q0(i)
      qq1 = q1(i)
      qq2 = q2(i)
      qq3 = q3(i)
      pp0 = p0(i)
      pp1 = p1(i)
      pp2 = p2(i)
      pp3 = p3(i)

      call nosquish(ttt,qq0,qq1,qq2,qq3,pp0,pp1,pp2,pp3) 

      q0(i) = qq0
      q1(i) = qq1
      q2(i) = qq2
      q3(i) = qq3
      p0(i) = pp0
      p1(i) = pp1
      p2(i) = pp2
      p3(i) = pp3

   end do
!$omp end do
!$omp barrier

!      volm=volm0
!      box = box0
!      chip = chipt
!      chit = chitt
!      conint = conintt

!   do j=1,nm

!      q0(j) = q0tt(j)
!      q1(j) = q1tt(j)
!      q2(j) = q2tt(j)
!      q3(j) = q3tt(j)
!      x(j)  = x0t(j)
!      y(j)  = y0t(j)
!      z(j)  = z0t(j)
!      vx(j) = vx0t(j)
!      vy(j) = vy0t(j)
!      vz(j) = vz0t(j)
!      jx(j) = jx0t(j)
!      jy(j) = jy0t(j)
!      jz(j) = jz0t(j)

!   end do

!   call sites

!$omp do schedule(dynamic)
   do i=1,nm
!      qq(0) = q0(i)
!      qq(1) = q1(i)
!      qq(2) = q2(i)
!      qq(3) = q3(i)

!      call quat_to_rot(qq,rott)

       ! update COM velocities and angular momentum to 1/2 tstep

      opx=0.5d0*(-q1(i)*p0(i)+q0(i)*p1(i)+ &
                 q3(i)*p2(i)-q2(i)*p3(i))
      opy=0.5d0*(-q2(i)*p0(i)-q3(i)*p1(i)+ &
                 q0(i)*p2(i)+q1(i)*p3(i))
      opz=0.5d0*(-q3(i)*p0(i)+q2(i)*p1(i)- &
                 q1(i)*p2(i)+q0(i)*p3(i))

      jx(i) = opx
      jy(i) = opy
      jz(i) = opz
   end do

!$omp end do
!$omp end parallel

!      vx(i) = vx(i) + 0.5d0*fx(i)*tt/totm
!      vy(i) = vy(i) + 0.5d0*fy(i)*tt/totm
!      vz(i) = vz(i) + 0.5d0*fz(i)*tt/totm

     ! first stage of velocity verlet algorithm

     ! move centre of mass by full time step

!      x(i)= x(i) + tt * vx(i)
!      y(i)= y(i) + tt * vy(i)
!      z(i)= z(i) + tt * vz(i)

      ! periodic boundary conditions

!    do i=1,nm

!      x(i) = x(i) - box*nint(x(i)/box)
!      y(i) = y(i) - box*nint(y(i)/box)
!      z(i) = z(i) - box*nint(z(i)/box)

!   end do

! end of first stage of velocity verlet algorithm

   call sites

   lvdwpe = vdwpe

	! DBG
	iniTime = SECNDS(0.0);
	! DBG
   if (pot_name.eq.'tip4p' ) then 
       call lj(x,y,z,xs,ys,zs,box,lvdwpe)
   else if (pot_name.eq.'tip5p' ) then
       call lj(x,y,z,xs,ys,zs,box,lvdwpe)
   else if (pot_name.eq.'sapt5s' ) then
       call sapt5s(x,y,z,xs,ys,zs,box,lvdwpe)
   else if (pot_name.eq.'ccpol5s' ) then
       call ccpol5s(x,y,z,xs,ys,zs,box,lvdwpe)
   else if (pot_name.eq.'ccpol8s' ) then
       call ccpol8s(x,y,z,xs,ys,zs,box,lvdwpe)       
   else if (pot_name.eq.'ccpol8s_omo' ) then
       call ccpol8s_omo(x,y,z,xs,ys,zs,box,lvdwpe)       
   else if (pot_name.eq.'ccpol23plus' ) then
       call ccpol23plus(x,y,z,xs,ys,zs,box,lvdwpe)
		! DBG
		deltaTime = SECNDS(iniTime);
!		write(dbgUnit, *) 'P#', MPIrank, '. ccpol23plus time	', deltaTime;
		! DBG
   end if

   vdwpe = lvdwpe

   lscpe = scpe
   call ewald_self_correction(alpha,box,lscpe)
   scpe = lscpe

   lrcpe = rcpe
   call rspace_ewald(x,y,z,xs,ys,zs,box,alpha,lrcpe)
   rcpe = lrcpe

   lkcpe = kcpe
   call kspace_ewald(x,y,z,xs,ys,zs,kmax1,kmax2,kmax3,alpha,box,lkcpe)
   kcpe = lkcpe

!  calling N-Body Iterative Induction model
   if (induction) then
      if(ind_model .eq. 'point_dipole') then

        lindpe = indpe
        call nb_induction(x,y,z,xs,ys,zs,box,lindpe)
        indpe = lindpe

      end if
   end if
 
! calling 3B-potenital (FS2 and FS3 terms)

  if (nonadditive_3B) then

     lfs2pe = fs2pe
     lfs3pe = fs3pe
     call threebody_term(x,y,z,xs,ys,zs,box,lfs2pe,lfs3pe)
     fs2pe = lfs2pe
     fs3pe = lfs3pe

  end if

! second stage of velocity verlet algorithm

!$omp parallel DEFAULT(SHARED) private(id,i,j,is,js,lvir)&
!$omp& private(qfx,qfy,qfz,qtx,qty,qtz, tqx,tqy,tqz)&
!$omp& private(qvx,qvy,qvz,qx,qy,qz,qjj,omg,qq,rott)&
!$omp& private(qt0,qt1,qt2,qt3,opx,opy,opz,cconorm,qmass)&
!$omp& private(hstep,qtke,qrke,qsystmp,dof,chit,conint,consv)&
!$omp& private(qq0,qq1,qq2,qq3,pp0,pp1,pp2,pp3)

   lvir = 0.0D0

!$omp do schedule(dynamic)
   do j=1,nm

      qfx=0.0
      qfy=0.0
      qfz=0.0
      qtx=0.0
      qty=0.0
      qtz=0.0

      do js=1,ns
         qfx=qfx+fxs(js,j)
         qfy=qfy+fys(js,j)
         qfz=qfz+fzs(js,j)
         qtx=qtx+ys(js,j)*fzs(js,j)-zs(js,j)*fys(js,j)
         qty=qty+zs(js,j)*fxs(js,j)-xs(js,j)*fzs(js,j)
         qtz=qtz+xs(js,j)*fys(js,j)-ys(js,j)*fxs(js,j)

         ! CONVERT FROM SITE-SITE TO MOLECULE-MOLECULE VIRIAL

         lvir = lvir + fxs(js,j)*xs(js,j)+fys(js,j)*ys(js,j)+fzs(js,j)*zs(js,j)

      end do

      fx(j)=qfx
      fy(j)=qfy
      fz(j)=qfz
      tx(j)=qtx
      ty(j)=qty
      tz(j)=qtz

      if (nonadditive_3B) then

         fx(j) = fx(j) + fxfs2(j) + fxfs3(j)
         fy(j) = fy(j) + fyfs2(j) + fyfs3(j)
         fz(j) = fz(j) + fzfs2(j) + fzfs3(j)

         tx(j)= tx(j) + txfs2(j) + txfs3(j)
         ty(j)= ty(j) + tyfs2(j) + tyfs3(j)
         tz(j)= tz(j) + tzfs2(j) + tzfs3(j)
      end if

   end do
!$omp end do
!$omp atomic

   vircom = vircom + lvir

!$omp barrier

!$omp do schedule(dynamic)
   do i=1,nm

      ! calculate quaternion momenta at start of time step
      opx = jx(i)
      opy = jy(i)
      opz = jz(i)

      p0(i)=2.0d0*(-q1(i)*opx-q2(i)*opy-q3(i)*opz)
      p1(i)=2.0d0*( q0(i)*opx-q3(i)*opy+q2(i)*opz)
      p2(i)=2.0d0*( q3(i)*opx+q0(i)*opy-q1(i)*opz)
      p3(i)=2.0d0*(-q2(i)*opx+q1(i)*opy+q0(i)*opz)
   end do
!$omp end do
!$omp barrier

!$omp do schedule(dynamic)
   do i=1,nm

       ! current rotation matrix

      call quat_to_rot(q0(i),q1(i),q2(i),q3(i),rott)

       ! calculate torque in principal frame

      tqx=tx(i)*rott(1)+ty(i)*rott(4)+tz(i)*rott(7)
      tqy=tx(i)*rott(2)+ty(i)*rott(5)+tz(i)*rott(8)
      tqz=tx(i)*rott(3)+ty(i)*rott(6)+tz(i)*rott(9)

       ! calculate quaternion torques

      qt0=2.0d0*(-q1(i)*tqx-q2(i)*tqy-q3(i)*tqz)
      qt1=2.0d0*( q0(i)*tqx-q3(i)*tqy+q2(i)*tqz)
      qt2=2.0d0*( q3(i)*tqx+q0(i)*tqy-q1(i)*tqz)
      qt3=2.0d0*(-q2(i)*tqx+q1(i)*tqy+q0(i)*tqz)

      ! calculate quaternion momenta at start of time step
!      opx = jx(i)
!      opy = jy(i)
!      opz = jz(i)

!      p0(i)=2.0d0*(-q1(i)*opx-q2(i)*opy-q3(i)*opz)
!      p1(i)=2.0d0*( q0(i)*opx-q3(i)*opy+q2(i)*opz)
!      p2(i)=2.0d0*( q3(i)*opx+q0(i)*opy-q1(i)*opz)
!      p3(i)=2.0d0*(-q2(i)*opx+q1(i)*opy+q0(i)*opz)

       ! update quaternion momenta by 1/2 time step

      p0(i)=p0(i)+qt0*0.5d0*ttt
      p1(i)=p1(i)+qt1*0.5d0*ttt
      p2(i)=p2(i)+qt2*0.5d0*ttt
      p3(i)=p3(i)+qt3*0.5d0*ttt

      !  update Rgid body angular momentum & COM velocities to full step

!      opx=0.5d0*(-q1(i)*p0(i)+q0(i)*p1(i)+ &
!                 q3(i)*p2(i)-q2(i)*p3(i))
!      opy=0.5d0*(-q2(i)*p0(i)-q3(i)*p1(i)+ &
!                 q0(i)*p2(i)+q1(i)*p3(i))
!      opz=0.5d0*(-q3(i)*p0(i)+q2(i)*p1(i)- &
!                 q1(i)*p2(i)+q0(i)*p3(i))

!      jx(i) = opx
!      jy(i) = opy
!      jz(i) = opz

      vx(i) = vx(i) + 0.5d0*fx(i)*ttt/totm
      vy(i) = vy(i) + 0.5d0*fy(i)*ttt/totm
      vz(i) = vz(i) + 0.5d0*fz(i)*ttt/totm

   end do
!$omp end do
!$omp barrier

!$omp do schedule(dynamic)
   do i=1,nm
!      qq(0) = q0(i)
!      qq(1) = q1(i)
!      qq(2) = q2(i)
!      qq(3) = q3(i)

!      call quat_to_rot(qq,rott)

       ! update COM velocities and angular momentum to 1/2 tstep

      opx=0.5d0*(-q1(i)*p0(i)+q0(i)*p1(i)+ &
                 q3(i)*p2(i)-q2(i)*p3(i))
      opy=0.5d0*(-q2(i)*p0(i)-q3(i)*p1(i)+ &
                 q0(i)*p2(i)+q1(i)*p3(i))
      opz=0.5d0*(-q3(i)*p0(i)+q2(i)*p1(i)- &
                 q1(i)*p2(i)+q0(i)*p3(i))

      jx(i) = opx
      jy(i) = opy
      jz(i) = opz
   end do
!$omp end do
!$omp end parallel

   qvir = vir

   do i=icyc,ncyc

!  npt thermostat
   call nptqscl_t(qstep,qmass,pmass,dof,taut,chip,chit,conint)
!  npt barostat
   call nptqscl_p(fstep,box,pmass,dof,chit,chip,qvir,qpres,conint)
!  npt thermostat
   call nptqscl_t(qstep,qmass,pmass,dof,taut,chip,chit,conint)

   end do

!   scale cell vectors - isotropic

!    volm=box**3.d0
    scale=(volm/volm0)**(1.d0/3.d0)

    box=box0*scale

!    volm=box**3.d0


! conserved quantity less kinetic and potential energy terms

   consv=conint+0.5d0*qmass*chit**2+pres0*volm+0.5d0*pmass*chip**2.d0

   ltke = tke
   call tkinetic_eng (nm, totm, vx, vy, vz, ltke)
   tke = ltke
   lttemp = ttemp
   call trans_temp (nm, tke, lttemp)
   ttemp = lttemp
   stptke = tke/418.4
   stpttp = ttemp

   lrke = rke
   call rkinetic_eng (nm, jx, jy, jz, lrke)
   rke = lrke
   stprke = rke/418.4

   lrtemp = rtemp
   call rot_temp (nm, rke, lrtemp)
   rtemp = lrtemp
   stprtp = rtemp

   lstptmp = stptmp
   call system_temp(nm,tke,rke,lstptmp)
   stptmp = lstptmp

   stpvir = vir/418.4d0
   pres = (2.d0*tke - vir)/(3.d0*volm)
   stpprs = prsunt*pres

   stpvolm = volm

!   write(*,*) volm0,volm,box0,box

   deallocate(p0,p1,p2,p3)
   deallocate(q0tt,q1tt,q2tt,q3tt)
   deallocate(x0t,y0t,z0t)
   deallocate(vx0t,vy0t,vz0t)
   deallocate(jx0t,jy0t,jz0t)


	! DBG
	deltaTimeTot = SECNDS(iniTimeTot);
!	write(dbgUnit, *) 'P#', MPIrank, '. nvt_hov time	', deltaTimeTot;
	call MPI_BARRIER(MPI_COMM_WORLD, MPIierr);
!	STOP;
	! DBG

  end subroutine npt_hov


  subroutine npt_hov_2(ttt)
   implicit none
   real(8) :: ttt

   integer :: i,j,is,js,icyc
   real(8) :: qfx,qfy,qfz
   real(8) :: qtx,qty,qtz, tqx,tqy,tqz
   real(8) :: qvx,qvy,qvz
   real(8) :: qx,qy,qz
   real(8) :: qjj(3),omg(3)
   real(8) :: qq(0:3), rott(1:9)
   real(8) :: qt0,qt1,qt2,qt3
   real(8) :: opx,opy,opz
   real(8) :: cconorm,lvir,qvir,qpres
   real(8) :: lvdwpe,lscpe,lrcpe,lkcpe,lindpe,lfs2pe,lfs3pe
   real(8) :: ltke,lttemp,lrke,lrtemp,lstptmp,scale
   real(8) :: qq0,qq1,qq2,qq3,pp0,pp1,pp2,pp3
!   real(8) :: q0t(nom),q1t(nom),q2t(nom),q3t(nom)
!   real(8) :: p0(nom),p1(nom),p2(nom),p3(nom)
    real(8), allocatable :: p0(:),p1(:),p2(:),p3(:)

    real(8), allocatable :: p0t(:),p1t(:),p2t(:),p3t(:)
    real(8), allocatable :: q0tt(:),q1tt(:),q2tt(:),q3tt(:)
    real(8), allocatable :: x0t(:),y0t(:),z0t(:)
    real(8), allocatable :: vx0t(:),vy0t(:),vz0t(:)
    real(8), allocatable :: jx0t(:),jy0t(:),jz0t(:)

! Nose Hoover parameters
   real(8),parameter :: taut = 0.5d0 ! ps (time constant)
   real(8),parameter :: taup = 0.5d0 ! ps (time constant check!!)
   integer, parameter :: ncyc = 5
   real(8) :: qmass,hstep,qtke,qrke,qsystmp,dof,chit,conint,volmi
   real(8) :: pmass,consv,chip,chiti,chipi,coninti,qstep,fstep,boxi
        ! DBG
        real(4) iniTime, deltaTime;
        real(4) iniTimeTot, deltaTimeTot;
        ! DBG

        ! DBG
        iniTimeTot = SECNDS(0.0);
        ! DBG

   allocate(p0(nm),p1(nm),p2(nm),p3(nm))
   allocate(q0tt(nm),q1tt(nm),q2tt(nm),q3tt(nm))
   allocate(x0t(nm),y0t(nm),z0t(nm))
   allocate(vx0t(nm),vy0t(nm),vz0t(nm))
   allocate(jx0t(nm),jy0t(nm),jz0t(nm))

   call system_mass

!   call system_vcom

   ! GRU
   chip = 0.0D0
   chit = 0.0D0
   conint = 0.0D0
   ! GRU

!$omp parallel DEFAULT(SHARED) private(id,i,j,is,js,lvir)&
!$omp& private(qfx,qfy,qfz,qtx,qty,qtz, tqx,tqy,tqz)&
!$omp& private(qvx,qvy,qvz,qx,qy,qz,qjj,omg,qq,rott)&
!$omp& private(qt0,qt1,qt2,qt3,opx,opy,opz,cconorm,qmass,pmass)&
!$omp& private(hstep,qtke,qrke,qsystmp,dof,chit,conint,consv,boxt,volmt)&
!$omp& private(qq0,qq1,qq2,qq3,pp0,pp1,pp2,pp3,chip,chitt,chipt,conintt)

   lvir = 0.0D0
   vircom = 0.d0

!$omp do schedule(dynamic)

   do j=1,nm

      qfx=0.0
      qfy=0.0
      qfz=0.0
      qtx=0.0
      qty=0.0
      qtz=0.0

      do js=1,ns
         qfx=qfx+fxs(js,j)
         qfy=qfy+fys(js,j)
         qfz=qfz+fzs(js,j)
         qtx=qtx+ys(js,j)*fzs(js,j)-zs(js,j)*fys(js,j)
         qty=qty+zs(js,j)*fxs(js,j)-xs(js,j)*fzs(js,j)
         qtz=qtz+xs(js,j)*fys(js,j)-ys(js,j)*fxs(js,j)

         ! CONVERT FROM SITE-SITE TO MOLECULE-MOLECULE VIRIAL

         lvir = lvir + fxs(js,j)*xs(js,j)+fys(js,j)*ys(js,j)+fzs(js,j)*zs(js,j)

      end do

      fx(j)=qfx
      fy(j)=qfy
      fz(j)=qfz
      tx(j)=qtx
      ty(j)=qty
      tz(j)=qtz

      if (nonadditive_3B) then

         fx(j) = fx(j) + fxfs2(j) + fxfs3(j)
         fy(j) = fy(j) + fyfs2(j) + fyfs3(j)
         fz(j) = fz(j) + fzfs2(j) + fzfs3(j)

         tx(j)= tx(j) + txfs2(j) + txfs3(j)
         ty(j)= ty(j) + tyfs2(j) + tyfs3(j)
         tz(j)= tz(j) + tzfs2(j) + tzfs3(j)


      end if

   end do
!$omp end do
!$omp atomic

   vircom = vircom + lvir

!$omp barrier
!$omp end parallel

!!$omp master
!  timestep parameters
   hstep = 0.5d0*ttt
   qstep = 0.25d0*ttt/dble(ncyc)
   fstep = 0.5d0*ttt/dble(ncyc)

! temporary saved variables

   box0=box
   volm0=box0**3.d0

   boxi = box
   volmi = boxi**3.d0

   chipi = chip
   chiti = chit
   coninti = conint

!   volm=volm0
!   volmt=volm
!   box=box0
!   boxt = box
!   chipt = chip
!   chitt = chit
!   conintt = conint

   qvir = vir+vircom

!   do j=1,nm

!      q0tt(j) = q0(j)
!      q1tt(j) = q1(j)
!      q2tt(j) = q2(j)
!      q3tt(j) = q3(j)
!      x0t(j)  = x(j)
!      y0t(j)  = y(j)
!      z0t(j)  = z(j)
!      vx0t(j) = vx(j)
!      vy0t(j) = vy(j)
!      vz0t(j) = vz(j)
!      jx0t(j) = jx(j)
!      jy0t(j) = jy(j)
!      jz0t(j) = jz(j)

!   end do

   dof = 3.d0*(2.d0*nm-1.d0) ! total dof transl + rot
!   nose-hoover inertia parameter
   qmass=dof*boltz*temp0*taut**2
   pmass=dof*boltz*temp0*taup**2

   do i=icyc,ncyc

!  npt thermostat
   call nptqscl_t(qstep,qmass,pmass,dof,taut,chip,chit,conint)
!  npt barostat
   call nptqscl_p(fstep,boxi,pmass,dof,chit,chip,qvir,qpres,conint)
!  npt thermostat
   call nptqscl_t(qstep,qmass,pmass,dof,taut,chip,chit,conint)

   end do

!   scale cell vectors - isotropic

!    volm=box**3.d0

    volmi=boxi**3.d0
    scale=(volmi/volm0)**(1.d0/3.d0)


!    scale=(volm/volm0)**(1.d0/3.d0)

!    box=box0*scale
!    volm=box**3.d0

!!$omp end master
!!$omp barrier

!$omp parallel DEFAULT(SHARED) private(id,i,j,is,js,lvir)&
!$omp& private(qfx,qfy,qfz,qtx,qty,qtz, tqx,tqy,tqz)&
!$omp& private(qvx,qvy,qvz,qx,qy,qz,qjj,omg,qq,rott)&
!$omp& private(qt0,qt1,qt2,qt3,opx,opy,opz,cconorm,qmass)&
!$omp& private(hstep,qtke,qrke,qsystmp,dof,chit,conint,consv,boxt,volmt)&
!$omp& private(qq0,qq1,qq2,qq3,pp0,pp1,pp2,pp3,chip,chitt,chipt,conintt)

!   *************  Rigid body motion ****************************
!     operations common to first and second stages

!$omp do schedule(dynamic)
   do i=1,nm

      ! calculate quaternion momenta at start of time step
      opx = jx(i)
      opy = jy(i)
      opz = jz(i)

      p0(i)=2.0d0*(-q1(i)*opx-q2(i)*opy-q3(i)*opz)
      p1(i)=2.0d0*( q0(i)*opx-q3(i)*opy+q2(i)*opz)
      p2(i)=2.0d0*( q3(i)*opx+q0(i)*opy-q1(i)*opz)
      p3(i)=2.0d0*(-q2(i)*opx+q1(i)*opy+q0(i)*opz)
   end do
!$omp end do
!$omp barrier

   call system_vcom

!$omp do schedule(dynamic)
   do i=1,nm

      ! current rotation matrix

      call quat_to_rot(q0(i),q1(i),q2(i),q3(i),rott)

       ! calculate torque in principal frame

      tqx=tx(i)*rott(1)+ty(i)*rott(4)+tz(i)*rott(7)
      tqy=tx(i)*rott(2)+ty(i)*rott(5)+tz(i)*rott(8)
      tqz=tx(i)*rott(3)+ty(i)*rott(6)+tz(i)*rott(9)

       ! calculate quaternion torques

      qt0=2.0d0*(-q1(i)*tqx-q2(i)*tqy-q3(i)*tqz)
      qt1=2.0d0*( q0(i)*tqx-q3(i)*tqy+q2(i)*tqz)
      qt2=2.0d0*( q3(i)*tqx+q0(i)*tqy-q1(i)*tqz)
      qt3=2.0d0*(-q2(i)*tqx+q1(i)*tqy+q0(i)*tqz)


      ! calculate quaternion momenta at start of time step
!      opx = jx(i)
!      opy = jy(i)
!      opz = jz(i)

!      p0(i)=2.0d0*(-q1(i)*opx-q2(i)*opy-q3(i)*opz)
!      p1(i)=2.0d0*( q0(i)*opx-q3(i)*opy+q2(i)*opz)
!      p2(i)=2.0d0*( q3(i)*opx+q0(i)*opy-q1(i)*opz)
!      p3(i)=2.0d0*(-q2(i)*opx+q1(i)*opy+q0(i)*opz)

       ! update quaternion momenta by 1/2 time step

      p0(i)=p0(i)+qt0*0.5d0*ttt
      p1(i)=p1(i)+qt1*0.5d0*ttt
      p2(i)=p2(i)+qt2*0.5d0*ttt
      p3(i)=p3(i)+qt3*0.5d0*ttt

      ! update centre of mass velocity by 1/2 time step

      vx(i) = vx(i) - vxcom
      vy(i) = vy(i) - vycom
      vz(i) = vz(i) - vzcom

      vx(i) = vx(i) + 0.5d0*fx(i)*ttt/totm
      vy(i) = vy(i) + 0.5d0*fy(i)*ttt/totm
      vz(i) = vz(i) + 0.5d0*fz(i)*ttt/totm

   end do
!$omp end do
!$omp barrier
    ! move centre of mass by full time step

    call system_com

!$omp do schedule(dynamic)
   do i=1,nm

      x(i)= x(i) + ttt * (vx(i)+chip*(x(i)-xcom))
      y(i)= y(i) + ttt * (vy(i)+chip*(y(i)-ycom))
      z(i)= z(i) + ttt * (vz(i)+chip*(z(i)-zcom))

!      x(i)= x(i) + ttt * vx(i)
!      y(i)= y(i) + ttt * vy(i)
!      z(i)= z(i) + ttt * vz(i)

      ! periodic boundary conditions

!!     x(i) = x(i) - box*nint(x(i)/box)
!!     y(i) = y(i) - box*nint(y(i)/box)
!!     z(i) = z(i) - box*nint(z(i)/box)
  end do

      ! rotate rigid groups: nosquish algorithm -update q to full tstep and
      ! amend p and get new rotation matrix

   do i=1,nm

      qq0 = q0(i)
      qq1 = q1(i)
      qq2 = q2(i)
      qq3 = q3(i)
      pp0 = p0(i)
      pp1 = p1(i)
      pp2 = p2(i)
      pp3 = p3(i)

      call nosquish(ttt,qq0,qq1,qq2,qq3,pp0,pp1,pp2,pp3)

      q0(i) = qq0
      q1(i) = qq1
      q2(i) = qq2
      q3(i) = qq3
      p0(i) = pp0
      p1(i) = pp1
      p2(i) = pp2
      p3(i) = pp3

   end do
!$omp end do
!$omp barrier

!      volm=volm0
!      box = box0
!      chip = chipt
!      chit = chitt
!      conint = conintt

!   do j=1,nm

!      q0(j) = q0tt(j)
!      q1(j) = q1tt(j)
!      q2(j) = q2tt(j)
!      q3(j) = q3tt(j)
!      x(j)  = x0t(j)
!      y(j)  = y0t(j)
!      z(j)  = z0t(j)
!      vx(j) = vx0t(j)
!      vy(j) = vy0t(j)
!      vz(j) = vz0t(j)
!      jx(j) = jx0t(j)
!      jy(j) = jy0t(j)
!      jz(j) = jz0t(j)

!   end do

!   call sites

!$omp do schedule(dynamic)
!   do i=1,nm
!      qq(0) = q0(i)
!      qq(1) = q1(i)
!      qq(2) = q2(i)
!      qq(3) = q3(i)

!      call quat_to_rot(qq,rott)

       ! update COM velocities and angular momentum to 1/2 tstep

!      opx=0.5d0*(-q1(i)*p0(i)+q0(i)*p1(i)+ &
!                 q3(i)*p2(i)-q2(i)*p3(i))
!      opy=0.5d0*(-q2(i)*p0(i)-q3(i)*p1(i)+ &
!                 q0(i)*p2(i)+q1(i)*p3(i))
!      opz=0.5d0*(-q3(i)*p0(i)+q2(i)*p1(i)- &
!                 q1(i)*p2(i)+q0(i)*p3(i))

!      jx(i) = opx
!      jy(i) = opy
!      jz(i) = opz
!   end do

!$omp end do
!$omp end parallel

!      vx(i) = vx(i) + 0.5d0*fx(i)*tt/totm
!      vy(i) = vy(i) + 0.5d0*fy(i)*tt/totm
!      vz(i) = vz(i) + 0.5d0*fz(i)*tt/totm

     ! first stage of velocity verlet algorithm

     ! move centre of mass by full time step

!      x(i)= x(i) + tt * vx(i)
!      y(i)= y(i) + tt * vy(i)
!      z(i)= z(i) + tt * vz(i)

      ! periodic boundary conditions

!    do i=1,nm

!      x(i) = x(i) - box*nint(x(i)/box)
!      y(i) = y(i) - box*nint(y(i)/box)
!      z(i) = z(i) - box*nint(z(i)/box)

!   end do

! end of first stage of velocity verlet algorithm

   call sites

   lvdwpe = vdwpe

        ! DBG
        iniTime = SECNDS(0.0);
        ! DBG
   if (pot_name.eq.'tip4p' ) then
       call lj(x,y,z,xs,ys,zs,box,lvdwpe)
   else if (pot_name.eq.'tip5p' ) then
       call lj(x,y,z,xs,ys,zs,box,lvdwpe)
   else if (pot_name.eq.'sapt5s' ) then
       call sapt5s(x,y,z,xs,ys,zs,box,lvdwpe)
   else if (pot_name.eq.'ccpol5s' ) then
       call ccpol5s(x,y,z,xs,ys,zs,box,lvdwpe)
   else if (pot_name.eq.'ccpol8s' ) then
       call ccpol8s(x,y,z,xs,ys,zs,box,lvdwpe)
   else if (pot_name.eq.'ccpol8s_omo' ) then
       call ccpol8s_omo(x,y,z,xs,ys,zs,box,lvdwpe)
   else if (pot_name.eq.'ccpol23plus' ) then
       call ccpol23plus(x,y,z,xs,ys,zs,box,lvdwpe)
                ! DBG
                deltaTime = SECNDS(iniTime);
!               write(dbgUnit, *) 'P#', MPIrank, '. ccpol23plus time    ',
!               deltaTime;
                ! DBG
   end if

   vdwpe = lvdwpe

   lscpe = scpe
   call ewald_self_correction(alpha,box,lscpe)
   scpe = lscpe

   lrcpe = rcpe
   call rspace_ewald(x,y,z,xs,ys,zs,box,alpha,lrcpe)
   rcpe = lrcpe

   lkcpe = kcpe
   call kspace_ewald(x,y,z,xs,ys,zs,kmax1,kmax2,kmax3,alpha,box,lkcpe)
   kcpe = lkcpe

!  calling N-Body Iterative Induction model
   if (induction) then
      if(ind_model .eq. 'point_dipole') then

        lindpe = indpe
        call nb_induction(x,y,z,xs,ys,zs,box,lindpe)
        indpe = lindpe

      end if
   end if

! calling 3B-potenital (FS2 and FS3 terms)

  if (nonadditive_3B) then

     lfs2pe = fs2pe
     lfs3pe = fs3pe
     call threebody_term(x,y,z,xs,ys,zs,box,lfs2pe,lfs3pe)
     fs2pe = lfs2pe
     fs3pe = lfs3pe

  end if

! second stage of velocity verlet algorithm

!$omp parallel DEFAULT(SHARED) private(id,i,j,is,js,lvir)&
!$omp& private(qfx,qfy,qfz,qtx,qty,qtz, tqx,tqy,tqz)&
!$omp& private(qvx,qvy,qvz,qx,qy,qz,qjj,omg,qq,rott)&
!$omp& private(qt0,qt1,qt2,qt3,opx,opy,opz,cconorm,qmass)&
!$omp& private(hstep,qtke,qrke,qsystmp,dof,chit,conint,consv)&
!$omp& private(qq0,qq1,qq2,qq3,pp0,pp1,pp2,pp3)

   lvir = 0.0D0
   vircom = 0.d0

!$omp do schedule(dynamic)
   do j=1,nm

      qfx=0.0
      qfy=0.0
      qfz=0.0
      qtx=0.0
      qty=0.0
      qtz=0.0

      do js=1,ns
         qfx=qfx+fxs(js,j)
         qfy=qfy+fys(js,j)
         qfz=qfz+fzs(js,j)
         qtx=qtx+ys(js,j)*fzs(js,j)-zs(js,j)*fys(js,j)
         qty=qty+zs(js,j)*fxs(js,j)-xs(js,j)*fzs(js,j)
         qtz=qtz+xs(js,j)*fys(js,j)-ys(js,j)*fxs(js,j)

         ! CONVERT FROM SITE-SITE TO MOLECULE-MOLECULE VIRIAL

         lvir = lvir + fxs(js,j)*xs(js,j)+fys(js,j)*ys(js,j)+fzs(js,j)*zs(js,j)

      end do

      fx(j)=qfx
      fy(j)=qfy
      fz(j)=qfz
      tx(j)=qtx
      ty(j)=qty
      tz(j)=qtz

      if (nonadditive_3B) then

         fx(j) = fx(j) + fxfs2(j) + fxfs3(j)
         fy(j) = fy(j) + fyfs2(j) + fyfs3(j)
         fz(j) = fz(j) + fzfs2(j) + fzfs3(j)

         tx(j)= tx(j) + txfs2(j) + txfs3(j)
         ty(j)= ty(j) + tyfs2(j) + tyfs3(j)
         tz(j)= tz(j) + tzfs2(j) + tzfs3(j)
      end if

   end do
!$omp end do
!$omp atomic

   vircom = vircom + lvir

!$omp barrier

!$omp do schedule(dynamic)
!   do i=1,nm

      ! calculate quaternion momenta at start of time step
!      opx = jx(i)
!      opy = jy(i)
!      opz = jz(i)

!      p0(i)=2.0d0*(-q1(i)*opx-q2(i)*opy-q3(i)*opz)
!      p1(i)=2.0d0*( q0(i)*opx-q3(i)*opy+q2(i)*opz)
!      p2(i)=2.0d0*( q3(i)*opx+q0(i)*opy-q1(i)*opz)
!      p3(i)=2.0d0*(-q2(i)*opx+q1(i)*opy+q0(i)*opz)
!   end do
!$omp end do
!$omp barrier

!$omp do schedule(dynamic)
   do i=1,nm

       ! current rotation matrix

      call quat_to_rot(q0(i),q1(i),q2(i),q3(i),rott)

       ! calculate torque in principal frame

      tqx=tx(i)*rott(1)+ty(i)*rott(4)+tz(i)*rott(7)
      tqy=tx(i)*rott(2)+ty(i)*rott(5)+tz(i)*rott(8)
      tqz=tx(i)*rott(3)+ty(i)*rott(6)+tz(i)*rott(9)

       ! calculate quaternion torques

      qt0=2.0d0*(-q1(i)*tqx-q2(i)*tqy-q3(i)*tqz)
      qt1=2.0d0*( q0(i)*tqx-q3(i)*tqy+q2(i)*tqz)
      qt2=2.0d0*( q3(i)*tqx+q0(i)*tqy-q1(i)*tqz)
      qt3=2.0d0*(-q2(i)*tqx+q1(i)*tqy+q0(i)*tqz)

      ! calculate quaternion momenta at start of time step
!      opx = jx(i)
!      opy = jy(i)
!      opz = jz(i)

!      p0(i)=2.0d0*(-q1(i)*opx-q2(i)*opy-q3(i)*opz)
!      p1(i)=2.0d0*( q0(i)*opx-q3(i)*opy+q2(i)*opz)
!      p2(i)=2.0d0*( q3(i)*opx+q0(i)*opy-q1(i)*opz)
!      p3(i)=2.0d0*(-q2(i)*opx+q1(i)*opy+q0(i)*opz)

       ! update quaternion momenta by 1/2 time step

      p0(i)=p0(i)+qt0*0.5d0*ttt
      p1(i)=p1(i)+qt1*0.5d0*ttt
      p2(i)=p2(i)+qt2*0.5d0*ttt
      p3(i)=p3(i)+qt3*0.5d0*ttt

      !  update Rgid body angular momentum & COM velocities to full step

!      opx=0.5d0*(-q1(i)*p0(i)+q0(i)*p1(i)+ &
!                 q3(i)*p2(i)-q2(i)*p3(i))
!      opy=0.5d0*(-q2(i)*p0(i)-q3(i)*p1(i)+ &
!                 q0(i)*p2(i)+q1(i)*p3(i))
!      opz=0.5d0*(-q3(i)*p0(i)+q2(i)*p1(i)- &
!                 q1(i)*p2(i)+q0(i)*p3(i))

!      jx(i) = opx
!      jy(i) = opy
!      jz(i) = opz

      vx(i) = vx(i) + 0.5d0*fx(i)*ttt/totm
      vy(i) = vy(i) + 0.5d0*fy(i)*ttt/totm
      vz(i) = vz(i) + 0.5d0*fz(i)*ttt/totm

   end do
!$omp end do
!$omp barrier

!$omp do schedule(dynamic)
   do i=1,nm
!      qq(0) = q0(i)
!      qq(1) = q1(i)
!      qq(2) = q2(i)
!      qq(3) = q3(i)

!      call quat_to_rot(qq,rott)

       ! update COM velocities and angular momentum to 1/2 tstep

      opx=0.5d0*(-q1(i)*p0(i)+q0(i)*p1(i)+ &
                 q3(i)*p2(i)-q2(i)*p3(i))
      opy=0.5d0*(-q2(i)*p0(i)-q3(i)*p1(i)+ &
                 q0(i)*p2(i)+q1(i)*p3(i))
      opz=0.5d0*(-q3(i)*p0(i)+q2(i)*p1(i)- &
                 q1(i)*p2(i)+q0(i)*p3(i))

      jx(i) = opx
      jy(i) = opy
      jz(i) = opz
   end do
!$omp end do
!$omp end parallel

   qvir = vir+vircom

   do i=icyc,ncyc

!  npt thermostat
   call nptqscl_t(qstep,qmass,pmass,dof,taut,chip,chit,conint)
!  npt barostat
   call nptqscl_p(fstep,boxi,pmass,dof,chit,chip,qvir,qpres,conint)
!  npt thermostat
   call nptqscl_t(qstep,qmass,pmass,dof,taut,chip,chit,conint)

   end do

!   scale cell vectors - isotropic

    volmi=boxi**3.d0
    scale=(volmi/volm0)**(1.d0/3.d0)

!    box=box*scale

    box=boxi
    volm=box**3.d0

!    volm=box**3.d0
!    scale=(volm/volm0)**(1.d0/3.d0)

!    box=box0*scale

!    volm=box**3.d0


! conserved quantity less kinetic and potential energy terms

   consv=conint+0.5d0*qmass*chit**2+pres0*volm+0.5d0*pmass*chip**2.d0

   ltke = tke
   call tkinetic_eng (nm, totm, vx, vy, vz, ltke)
   tke = ltke
   lttemp = ttemp
   call trans_temp (nm, tke, lttemp)
   ttemp = lttemp
   stptke = tke/418.4
   stpttp = ttemp

   lrke = rke
   call rkinetic_eng (nm, jx, jy, jz, lrke)
   rke = lrke
   stprke = rke/418.4

   lrtemp = rtemp
   call rot_temp (nm, rke, lrtemp)
   rtemp = lrtemp
   stprtp = rtemp

   lstptmp = stptmp
   call system_temp(nm,tke,rke,lstptmp)
   stptmp = lstptmp

   qvir = vir+vircom

!  calculate system pressure

   pres = (2.d0*ltke - qvir)/(3.d0*volm)

!   write(*,*) box,box**3.d0,volm,pres,vir,ltke

   stpvolm = volm

   stpvir = qvir/418.4d0
   stpprs = prsunt*pres

!   write(*,*) temp0,stptmp!volm0,volm,box0,box

   deallocate(p0,p1,p2,p3)
   deallocate(q0tt,q1tt,q2tt,q3tt)
   deallocate(x0t,y0t,z0t)
   deallocate(vx0t,vy0t,vz0t)
   deallocate(jx0t,jy0t,jz0t)


        ! DBG
        deltaTimeTot = SECNDS(iniTimeTot);
!       write(dbgUnit, *) 'P#', MPIrank, '. nvt_hov time        ', deltaTimeTot;
        call MPI_BARRIER(MPI_COMM_WORLD, MPIierr);
!       STOP;
        ! DBG

  end subroutine npt_hov_2


  subroutine npt_hov_3(ttt)
   implicit none
   real(8) :: ttt

   integer :: i,j,is,js,icyc
   real(8) :: qfx,qfy,qfz
   real(8) :: qtx,qty,qtz, tqx,tqy,tqz
   real(8) :: qvx,qvy,qvz,ctx,cty,ctz
   real(8) :: qx,qy,qz
   real(8) :: qjj(3),omg(3)
   real(8) :: qq(0:3), rott(1:9)
   real(8) :: qt0,qt1,qt2,qt3
   real(8) :: opx,opy,opz
   real(8) :: cconorm,lvir,qvir,qpres
   real(8) :: lvdwpe,lscpe,lrcpe,lkcpe,lindpe,lfs2pe,lfs3pe
   real(8) :: ltke,lttemp,lrke,lrtemp,lstptmp,scale
   real(8) :: qq0,qq1,qq2,qq3,pp0,pp1,pp2,pp3,vtx,vty,vtz
!   real(8) :: q0t(nom),q1t(nom),q2t(nom),q3t(nom)
!   real(8) :: p0(nom),p1(nom),p2(nom),p3(nom)
    real(8), allocatable :: p0(:),p1(:),p2(:),p3(:)

    real(8), allocatable :: p0t(:),p1t(:),p2t(:),p3t(:)
    real(8), allocatable :: q0tt(:),q1tt(:),q2tt(:),q3tt(:)
    real(8), allocatable :: x0t(:),y0t(:),z0t(:)
    real(8), allocatable :: vx0t(:),vy0t(:),vz0t(:)
    real(8), allocatable :: jx0t(:),jy0t(:),jz0t(:)

! Nose Hoover parameters
   real(8),parameter :: taut = 0.5d0 ! ps (time constant)
   real(8),parameter :: taup = 0.5d0 ! ps (time constant check!!)
   integer, parameter :: ncyc = 5
   real(8) :: qmass,hstep,qtke,qrke,qsystmp,dof,chit,conint,volmi
   real(8) :: pmass,consv,chip,chiti,chipi,coninti,qstep,fstep,boxi
        ! DBG
        real(4) iniTime, deltaTime;
        real(4) iniTimeTot, deltaTimeTot;
        ! DBG

        ! DBG
        iniTimeTot = SECNDS(0.0);
        ! DBG

   allocate(p0(nm),p1(nm),p2(nm),p3(nm))
   allocate(q0tt(nm),q1tt(nm),q2tt(nm),q3tt(nm))
   allocate(x0t(nm),y0t(nm),z0t(nm))
   allocate(vx0t(nm),vy0t(nm),vz0t(nm))
   allocate(jx0t(nm),jy0t(nm),jz0t(nm))

   call system_mass

   ! GRU
   chip = 0.0D0
   chit = 0.0D0
   conint = 0.0D0
   ! GRU

   lvir = 0.0D0
   vircom = 0.d0

   do j=1,nm

      qfx=0.0
      qfy=0.0
      qfz=0.0
      qtx=0.0
      qty=0.0
      qtz=0.0

      do js=1,ns
         qfx=qfx+fxs(js,j)
         qfy=qfy+fys(js,j)
         qfz=qfz+fzs(js,j)
         qtx=qtx+ys(js,j)*fzs(js,j)-zs(js,j)*fys(js,j)
         qty=qty+zs(js,j)*fxs(js,j)-xs(js,j)*fzs(js,j)
         qtz=qtz+xs(js,j)*fys(js,j)-ys(js,j)*fxs(js,j)

         ! CONVERT FROM SITE-SITE TO MOLECULE-MOLECULE VIRIAL

         lvir = lvir + fxs(js,j)*xs(js,j)+fys(js,j)*ys(js,j)+fzs(js,j)*zs(js,j)

      end do

      fx(j)=qfx
      fy(j)=qfy
      fz(j)=qfz
      tx(j)=qtx
      ty(j)=qty
      tz(j)=qtz

      if (nonadditive_3B) then

         fx(j) = fx(j) + fxfs2(j) + fxfs3(j)
         fy(j) = fy(j) + fyfs2(j) + fyfs3(j)
         fz(j) = fz(j) + fzfs2(j) + fzfs3(j)

         tx(j)= tx(j) + txfs2(j) + txfs3(j)
         ty(j)= ty(j) + tyfs2(j) + tyfs3(j)
         tz(j)= tz(j) + tzfs2(j) + tzfs3(j)


      end if

   end do

   vircom = vircom + lvir

!  timestep parameters
   hstep = 0.5d0*ttt
   qstep = 0.25d0*ttt !/dble(ncyc)
   fstep = 0.5d0*ttt  !/dble(ncyc)

! temporary saved variables

   box0=box
   volm0=box0**3.d0

   write(*,*) 'Before Scaling '
   write(*,*) 'box0 ',box0
   write(*,*) 'volm0 ',volm0

   boxi = box
   volmi = boxi**3.d0

   chipi = chip
   chiti = chit
   coninti = conint

   qvir = vir + vircom

   dof = 3.d0*(2.d0*nm-1.d0) ! total dof transl + rot
!   nose-hoover inertia parameter
   qmass=dof*boltz*temp0*taut**2
   pmass=dof*boltz*temp0*taup**2

!   do i=icyc,ncyc

!*****************************************************************************
!******************  npt thermostat *****************************
!*****************************************************************************

!!   call nptqscl_t(qstep,qmass,pmass,dof,taut,chip,chit,conint)

   !   calculate kinetic energy and inst temp
   call tkinetic_eng (nm, totm, vx, vy, vz, qtke)
   call rkinetic_eng (nm, jx, jy, jz, qrke)
   call system_temp(nm,qtke,qrke,qsystmp)

!  update chit to 1/2 step
   chit=chit + 0.5d0*qstep*dof*boltz*(qsystmp-temp0)/qmass + &
               0.5d0*qstep*(pmass*chip**2.d0-boltz*temp0)/qmass

!  thermostat scale parameter
   scale=exp(-qstep*chit)   

!  thermostat rigid body velocities

   do i=1,nm

      vx(i) = scale*vx(i)
      vy(i) = scale*vy(i)
      vz(i) = scale*vz(i)

      jx(i) = scale*jx(i)
      jy(i) = scale*jy(i)
      jz(i) = scale*jz(i)

   end do

!  scale kinetic energy

   qtke = qtke*scale**2.d0
   qrke = qrke*scale**2.d0
   call system_temp(nm,qtke,qrke,qsystmp)

!  update chi to full step

   conint=conint+qstep*chit*(qmass/taut**2+boltz*temp0)

!  update chit to full step
   chit=chit + 0.5d0*qstep*dof*boltz*(qsystmp-temp0)/qmass + &
               0.5d0*qstep*(pmass*chip**2.d0-boltz*temp0)/qmass

   write(*,*) 'NPT thermostat'
   write(*,*) 'boxi ',boxi
   write(*,*) 'volmi ',volmi

!*****************************************************************************
!***********************  npt barostat ***********************************
!*****************************************************************************

!!   call nptqscl_p(fstep,boxi,pmass,dof,chit,chip,qvir,qpres,conint)

!   calculate kinetic energy and inst temp

   call tkinetic_eng (nm, totm, vx, vy, vz, qtke)

!   call rkinetic_eng (nm, jx, jy, jz, qrke)
!   call system_temp(nm,qtke,qrke,qsystmp)

   volmi = boxi**3

   qpres = (2.d0*qtke - qvir)/(3.d0*volmi)

!  update chit to 1/2 step
   chip=chip + 0.5d0*hstep*(3.d0*volmi*(qpres-pres0)/pmass - chip*chit)

!  barostat scale parameter
   scale=exp(-hstep*chip)

!  thermostat rigid body velocities

   do i=1,nm

      vx(i) = scale*vx(i)
      vy(i) = scale*vy(i)
      vz(i) = scale*vz(i)

!      jx(i) = scale*jx(i)
!      jy(i) = scale*jy(i)
!      jz(i) = scale*jz(i)

   end do

!  scale kinetic energy

   qtke = qtke*scale**2.d0

!   qrke = qrke*scale**2                   ! check

!   call system_temp(nm,qtke,qrke,qsystmp) ! check

!   volmi=boxi**3.d0

   volmi=volmi*exp(3.d0*hstep*chip)

!  update chip to full step
   boxi = volmi**(1.d0/3.d0)

!!   boxi=boxi*exp(3.d0*ttt*chip)
!!   volmi=boxi**3.d0

   qpres = (2.d0*qtke - qvir)/(3.d0*volmi)

!  update chit to 1/2 step
   chip=chip + 0.5d0*hstep*(3.d0*volmi*(qpres-pres0)/pmass - chip*chit)

   write(*,*) 'NPT barostat' 
   write(*,*) 'boxi ',boxi
   write(*,*) 'volmi ',volmi

!*****************************************************************************
!***********************  npt thermostat ********************************
!*****************************************************************************
!!   call nptqscl_t(qstep,qmass,pmass,dof,taut,chip,chit,conint)

   !   calculate kinetic energy and inst temp
   call tkinetic_eng (nm, totm, vx, vy, vz, qtke)
   call rkinetic_eng (nm, jx, jy, jz, qrke)
   call system_temp(nm,qtke,qrke,qsystmp)

!  update chit to 1/2 step
   chit=chit + 0.5d0*qstep*dof*boltz*(qsystmp-temp0)/qmass + &
               0.5d0*qstep*(pmass*chip**2.d0-boltz*temp0)/qmass

!  thermostat scale parameter
   scale=exp(-qstep*chit)

!  thermostat rigid body velocities

   do i=1,nm

      vx(i) = scale*vx(i)
      vy(i) = scale*vy(i)
      vz(i) = scale*vz(i)

      jx(i) = scale*jx(i)
      jy(i) = scale*jy(i)
      jz(i) = scale*jz(i)

   end do

!  scale kinetic energy

   qtke = qtke*scale**2.d0
   qrke = qrke*scale**2.d0
   call system_temp(nm,qtke,qrke,qsystmp)

!  update chi to full step

   conint=conint+qstep*chit*(qmass/taut**2+boltz*temp0)

!  update chit to full step
   chit=chit + 0.5d0*qstep*dof*boltz*(qsystmp-temp0)/qmass + &
               0.5d0*qstep*(pmass*chip**2.d0-boltz*temp0)/qmass

   write(*,*) 'NPT thermostat'
   write(*,*) 'boxi ',boxi
   write(*,*) 'volmi ',volmi

!*****************************************************************************
!*****************************************************************************
!   end do

!   scale cell vectors - isotropic

    volmi=boxi**3.d0
    scale=(volmi/volm0)**(1.d0/3.d0)

   write(*,*) 'After Scaling'
   write(*,*) 'boxi ',boxi
   write(*,*) 'volmi ',volmi

!   *************  Rigid body motion ****************************
!     operations common to first and second stages

   do i=1,nm

      ! calculate quaternion momenta at start of time step
      opx = jx(i)
      opy = jy(i)
      opz = jz(i)

      p0(i)=2.0d0*(-q1(i)*opx-q2(i)*opy-q3(i)*opz)
      p1(i)=2.0d0*( q0(i)*opx-q3(i)*opy+q2(i)*opz)
      p2(i)=2.0d0*( q3(i)*opx+q0(i)*opy-q1(i)*opz)
      p3(i)=2.0d0*(-q2(i)*opx+q1(i)*opy+q0(i)*opz)
   end do

   do i=1,nm

      ! current rotation matrix

      call quat_to_rot(q0(i),q1(i),q2(i),q3(i),rott)

       ! calculate torque in principal frame

      tqx=tx(i)*rott(1)+ty(i)*rott(4)+tz(i)*rott(7)
      tqy=tx(i)*rott(2)+ty(i)*rott(5)+tz(i)*rott(8)
      tqz=tx(i)*rott(3)+ty(i)*rott(6)+tz(i)*rott(9)

       ! calculate quaternion torques

      qt0=2.0d0*(-q1(i)*tqx-q2(i)*tqy-q3(i)*tqz)
      qt1=2.0d0*( q0(i)*tqx-q3(i)*tqy+q2(i)*tqz)
      qt2=2.0d0*( q3(i)*tqx+q0(i)*tqy-q1(i)*tqz)
      qt3=2.0d0*(-q2(i)*tqx+q1(i)*tqy+q0(i)*tqz)

       ! update quaternion momenta by 1/2 time step

      p0(i)=p0(i)+qt0*0.5d0*ttt
      p1(i)=p1(i)+qt1*0.5d0*ttt
      p2(i)=p2(i)+qt2*0.5d0*ttt
      p3(i)=p3(i)+qt3*0.5d0*ttt

   end do

      ! update centre of mass velocity by 1/2 time step

      ! remove system centre of mass velocity  
      call system_vcom 

   do i=1,nm 

      vx(i) = vx(i) - vxcom
      vy(i) = vy(i) - vycom
      vz(i) = vz(i) - vzcom

      vx(i) = vx(i) + 0.5d0*fx(i)*ttt/totm
      vy(i) = vy(i) + 0.5d0*fy(i)*ttt/totm
      vz(i) = vz(i) + 0.5d0*fz(i)*ttt/totm

   end do

    ! move centre of mass by full time step

   call system_com

   do i=1,nm

      ctx = x(i)-xcom
      cty = y(i)-ycom
      ctz = z(i)-zcom

      x(i)= x(i) + ttt * (vx(i)+chip*ctx)
      y(i)= y(i) + ttt * (vy(i)+chip*cty)
      z(i)= z(i) + ttt * (vz(i)+chip*ctz)

!      x(i)= x(i) + ttt * (vx(i))
!      y(i)= y(i) + ttt * (vy(i))
!      z(i)= z(i) + ttt * (vz(i))
 
  end do

      ! rotate rigid groups: nosquish algorithm -update q to full tstep and
      ! amend p and get new rotation matrix

   do i=1,nm

      qq0 = q0(i)
      qq1 = q1(i)
      qq2 = q2(i)
      qq3 = q3(i)
      pp0 = p0(i)
      pp1 = p1(i)
      pp2 = p2(i)
      pp3 = p3(i)

      call nosquish(ttt,qq0,qq1,qq2,qq3,pp0,pp1,pp2,pp3)

      q0(i) = qq0
      q1(i) = qq1
      q2(i) = qq2
      q3(i) = qq3
      p0(i) = pp0
      p1(i) = pp1
      p2(i) = pp2
      p3(i) = pp3

   end do

! end of first stage of velocity verlet algorithm

   call sites

   lvdwpe = vdwpe

        ! DBG
        iniTime = SECNDS(0.0);
        ! DBG
   if (pot_name.eq.'tip4p' ) then
       call lj(x,y,z,xs,ys,zs,box,lvdwpe)
   else if (pot_name.eq.'tip5p' ) then
       call lj(x,y,z,xs,ys,zs,box,lvdwpe)
   else if (pot_name.eq.'sapt5s' ) then
       call sapt5s(x,y,z,xs,ys,zs,box,lvdwpe)
   else if (pot_name.eq.'ccpol5s' ) then
       call ccpol5s(x,y,z,xs,ys,zs,box,lvdwpe)
   else if (pot_name.eq.'ccpol8s' ) then
       call ccpol8s(x,y,z,xs,ys,zs,box,lvdwpe)
   else if (pot_name.eq.'ccpol8s_omo' ) then
       call ccpol8s_omo(x,y,z,xs,ys,zs,box,lvdwpe)
   else if (pot_name.eq.'ccpol23plus' ) then
       call ccpol23plus(x,y,z,xs,ys,zs,box,lvdwpe)
                ! DBG
                deltaTime = SECNDS(iniTime);
!               write(dbgUnit, *) 'P#', MPIrank, '. ccpol23plus time    ',
!               deltaTime;
                ! DBG
   end if

   vdwpe = lvdwpe

   lscpe = scpe
   call ewald_self_correction(alpha,box,lscpe)
   scpe = lscpe

   lrcpe = rcpe
   call rspace_ewald(x,y,z,xs,ys,zs,box,alpha,lrcpe)
   rcpe = lrcpe

   lkcpe = kcpe
   call kspace_ewald(x,y,z,xs,ys,zs,kmax1,kmax2,kmax3,alpha,box,lkcpe)
   kcpe = lkcpe

!  calling N-Body Iterative Induction model
   if (induction) then
      if(ind_model .eq. 'point_dipole') then

        lindpe = indpe
        call nb_induction(x,y,z,xs,ys,zs,box,lindpe)
        indpe = lindpe

      end if
   end if

! calling 3B-potenital (FS2 and FS3 terms)

  if (nonadditive_3B) then

     lfs2pe = fs2pe
     lfs3pe = fs3pe
     call threebody_term(x,y,z,xs,ys,zs,box,lfs2pe,lfs3pe)
     fs2pe = lfs2pe
     fs3pe = lfs3pe

  end if

! second stage of velocity verlet algorithm

   lvir = 0.0D0
   vircom = 0.d0

!$omp do schedule(dynamic)
   do j=1,nm

      qfx=0.0
      qfy=0.0
      qfz=0.0
      qtx=0.0
      qty=0.0
      qtz=0.0

      do js=1,ns
         qfx=qfx+fxs(js,j)
         qfy=qfy+fys(js,j)
         qfz=qfz+fzs(js,j)
         qtx=qtx+ys(js,j)*fzs(js,j)-zs(js,j)*fys(js,j)
         qty=qty+zs(js,j)*fxs(js,j)-xs(js,j)*fzs(js,j)
         qtz=qtz+xs(js,j)*fys(js,j)-ys(js,j)*fxs(js,j)

         ! CONVERT FROM SITE-SITE TO MOLECULE-MOLECULE VIRIAL

         lvir = lvir + fxs(js,j)*xs(js,j)+fys(js,j)*ys(js,j)+fzs(js,j)*zs(js,j)

      end do

      fx(j)=qfx
      fy(j)=qfy
      fz(j)=qfz
      tx(j)=qtx
      ty(j)=qty
      tz(j)=qtz

      if (nonadditive_3B) then

         fx(j) = fx(j) + fxfs2(j) + fxfs3(j)
         fy(j) = fy(j) + fyfs2(j) + fyfs3(j)
         fz(j) = fz(j) + fzfs2(j) + fzfs3(j)

         tx(j)= tx(j) + txfs2(j) + txfs3(j)
         ty(j)= ty(j) + tyfs2(j) + tyfs3(j)
         tz(j)= tz(j) + tzfs2(j) + tzfs3(j)
      end if

   end do

   vircom = vircom + lvir

   do i=1,nm

       ! current rotation matrix

      call quat_to_rot(q0(i),q1(i),q2(i),q3(i),rott)

       ! calculate torque in principal frame

      tqx=tx(i)*rott(1)+ty(i)*rott(4)+tz(i)*rott(7)
      tqy=tx(i)*rott(2)+ty(i)*rott(5)+tz(i)*rott(8)
      tqz=tx(i)*rott(3)+ty(i)*rott(6)+tz(i)*rott(9)

       ! calculate quaternion torques

      qt0=2.0d0*(-q1(i)*tqx-q2(i)*tqy-q3(i)*tqz)
      qt1=2.0d0*( q0(i)*tqx-q3(i)*tqy+q2(i)*tqz)
      qt2=2.0d0*( q3(i)*tqx+q0(i)*tqy-q1(i)*tqz)
      qt3=2.0d0*(-q2(i)*tqx+q1(i)*tqy+q0(i)*tqz)

       ! update quaternion momenta by 1/2 time step

      p0(i)=p0(i)+qt0*0.5d0*ttt
      p1(i)=p1(i)+qt1*0.5d0*ttt
      p2(i)=p2(i)+qt2*0.5d0*ttt
      p3(i)=p3(i)+qt3*0.5d0*ttt

      !  update Rgid body angular momentum & COM velocities to full step

      vx(i) = vx(i) + 0.5d0*fx(i)*ttt/totm
      vy(i) = vy(i) + 0.5d0*fy(i)*ttt/totm
      vz(i) = vz(i) + 0.5d0*fz(i)*ttt/totm

   end do

   do i=1,nm

       ! update COM velocities and angular momentum to 1/2 tstep

      opx=0.5d0*(-q1(i)*p0(i)+q0(i)*p1(i)+ &
                 q3(i)*p2(i)-q2(i)*p3(i))
      opy=0.5d0*(-q2(i)*p0(i)-q3(i)*p1(i)+ &
                 q0(i)*p2(i)+q1(i)*p3(i))
      opz=0.5d0*(-q3(i)*p0(i)+q2(i)*p1(i)- &
                 q1(i)*p2(i)+q0(i)*p3(i))

      jx(i) = opx
      jy(i) = opy
      jz(i) = opz
   end do

   qvir = vir + vircom

!   do i=icyc,ncyc

!*****************************************************************************
!******************  npt thermostat *****************************
!*****************************************************************************

!!   call nptqscl_t(qstep,qmass,pmass,dof,taut,chip,chit,conint)

   !   calculate kinetic energy and inst temp
   call tkinetic_eng (nm, totm, vx, vy, vz, qtke)
   call rkinetic_eng (nm, jx, jy, jz, qrke)
   call system_temp(nm,qtke,qrke,qsystmp)

!  update chit to 1/2 step
   chit=chit + 0.5d0*qstep*dof*boltz*(qsystmp-temp0)/qmass + &
               0.5d0*qstep*(pmass*chip**2.d0-boltz*temp0)/qmass

!  thermostat scale parameter
   scale=exp(-qstep*chit)

!  thermostat rigid body velocities

   do i=1,nm

      vx(i) = scale*vx(i)
      vy(i) = scale*vy(i)
      vz(i) = scale*vz(i)

      jx(i) = scale*jx(i)
      jy(i) = scale*jy(i)
      jz(i) = scale*jz(i)

   end do

!  scale kinetic energy

   qtke = qtke*scale**2.d0
   qrke = qrke*scale**2.d0
   call system_temp(nm,qtke,qrke,qsystmp)

!  update chi to full step

   conint=conint+qstep*chit*(qmass/taut**2+boltz*temp0)

!  update chit to full step
   chit=chit + 0.5d0*qstep*dof*boltz*(qsystmp-temp0)/qmass + &
               0.5d0*qstep*(pmass*chip**2.d0-boltz*temp0)/qmass

!*****************************************************************************
!***********************  npt barostat ***********************************
!*****************************************************************************
!!   call nptqscl_p(fstep,boxi,pmass,dof,chit,chip,qvir,qpres,conint)

!   calculate kinetic energy and inst temp

   call tkinetic_eng (nm, totm, vx, vy, vz, qtke)

!   call rkinetic_eng (nm, jx, jy, jz, qrke)
!   call system_temp(nm,qtke,qrke,qsystmp)

   volmi = boxi**3

   qpres = (2.d0*qtke - qvir)/(3.d0*volmi)

!  update chit to 1/2 step
   chip=chip + 0.5d0*hstep*(3.d0*volmi*(qpres-pres0)/pmass - chip*chit)

!  barostat scale parameter
   scale=exp(-hstep*chip)

!  thermostat rigid body velocities

   do i=1,nm

      vx(i) = scale*vx(i)
      vy(i) = scale*vy(i)
      vz(i) = scale*vz(i)

!      jx(i) = scale*jx(i)
!      jy(i) = scale*jy(i)
!      jz(i) = scale*jz(i)

   end do

!  scale kinetic energy

   qtke = qtke*scale**2.d0

!   qrke = qrke*scale**2                   ! check

!   call system_temp(nm,qtke,qrke,qsystmp) ! check

!   volmi=boxi**3.d0

   volmi=volmi*exp(3.d0*hstep*chip)

!  update chip to full step
   boxi = volmi**(1.d0/3.d0)

!!   boxi=boxi*exp(3.d0*ttt*chip)
!!   volmi=boxi**3.d0

   qpres = (2.d0*qtke - qvir)/(3.d0*volmi)

!  update chit to 1/2 step
   chip=chip + 0.5d0*hstep*(3.d0*volmi*(qpres-pres0)/pmass - chip*chit)

!*****************************************************************************
!***********************  npt thermostat ********************************
!*****************************************************************************
!!   call nptqscl_t(qstep,qmass,pmass,dof,taut,chip,chit,conint)

   !   calculate kinetic energy and inst temp
   call tkinetic_eng (nm, totm, vx, vy, vz, qtke)
   call rkinetic_eng (nm, jx, jy, jz, qrke)
   call system_temp(nm,qtke,qrke,qsystmp)

!  update chit to 1/2 step
   chit=chit + 0.5d0*qstep*dof*boltz*(qsystmp-temp0)/qmass + &
               0.5d0*qstep*(pmass*chip**2.d0-boltz*temp0)/qmass

!  thermostat scale parameter
   scale=exp(-qstep*chit)

!  thermostat rigid body velocities

   do i=1,nm

      vx(i) = scale*vx(i)
      vy(i) = scale*vy(i)
      vz(i) = scale*vz(i)

      jx(i) = scale*jx(i)
      jy(i) = scale*jy(i)
      jz(i) = scale*jz(i)

   end do

!  scale kinetic energy

   qtke = qtke*scale**2.d0
   qrke = qrke*scale**2.d0
   call system_temp(nm,qtke,qrke,qsystmp)

!  update chi to full step

   conint=conint+qstep*chit*(qmass/taut**2+boltz*temp0)

!  update chit to full step
   chit=chit + 0.5d0*qstep*dof*boltz*(qsystmp-temp0)/qmass + &
               0.5d0*qstep*(pmass*chip**2.d0-boltz*temp0)/qmass

!*****************************************************************************
!*****************************************************************************
!   end do

!   scale cell vectors - isotropic

    volmi=boxi**3.d0
    scale=(volmi/volm0)**(1.d0/3.d0)

    write(*,*) 'boxi ',boxi
    write(*,*) 'volmi ',volmi

!    box=box*scale

    box=boxi
    volm=box**3.d0

    write(*,*) 'Final volm '
    write(*,*) 'volm ',volm

! conserved quantity less kinetic and potential energy terms

   consv=conint+0.5d0*qmass*chit**2+pres0*volm+0.5d0*pmass*chip**2.d0

   ltke = tke
   call tkinetic_eng (nm, totm, vx, vy, vz, ltke)
   tke = ltke
   lttemp = ttemp
   call trans_temp (nm, tke, lttemp)
   ttemp = lttemp
   stptke = tke/418.4
   stpttp = ttemp

   lrke = rke
   call rkinetic_eng (nm, jx, jy, jz, lrke)
   rke = lrke
   stprke = rke/418.4

   lrtemp = rtemp
   call rot_temp (nm, rke, lrtemp)
   rtemp = lrtemp
   stprtp = rtemp

   lstptmp = stptmp
   call system_temp(nm,tke,rke,lstptmp)
   stptmp = lstptmp

   qvir = vir+vircom

!  calculate system pressure

   pres = (2.d0*ltke - qvir)/(3.d0*volm)

!   write(*,*) box,box**3.d0,volm,pres,vir,ltke

   stpvolm = volm

   stpvir = qvir/418.4d0
   stpprs = prsunt*pres

!   write(*,*) temp0,stptmp!volm0,volm,box0,box

   deallocate(p0,p1,p2,p3)
   deallocate(q0tt,q1tt,q2tt,q3tt)
   deallocate(x0t,y0t,z0t)
   deallocate(vx0t,vy0t,vz0t)
   deallocate(jx0t,jy0t,jz0t)


        ! DBG
        deltaTimeTot = SECNDS(iniTimeTot);
!       write(dbgUnit, *) 'P#', MPIrank, '. nvt_hov time        ', deltaTimeTot;
        call MPI_BARRIER(MPI_COMM_WORLD, MPIierr);
!       STOP;
        ! DBG

  end subroutine npt_hov_3

  subroutine npt_b(ttt)
   implicit none
   real(8) :: ttt

   integer :: i,j,is,js
   real(8) :: qfx,qfy,qfz
   real(8) :: qtx,qty,qtz, tqx,tqy,tqz
   real(8) :: qvx,qvy,qvz
   real(8) :: qx,qy,qz
   real(8) :: qjj(3),omg(3)
   real(8) :: qq(0:3), rott(1:9)
   real(8) :: qt0,qt1,qt2,qt3
   real(8) :: opx,opy,opz
   real(8) :: cconorm,chit,lvir,qvir,chip,scale
   real(8) :: lvdwpe,lscpe,lrcpe,lkcpe,lindpe,lfs2pe,lfs3pe
   real(8) :: ltke,lttemp,lrke,lrtemp,lstptmp
   real(8) :: qq0,qq1,qq2,qq3,pp0,pp1,pp2,pp3
   real(8) :: qtke,qrke,qsystmp,volmi,boxi,qpres
!   real(8) :: q0t(nom),q1t(nom),q2t(nom),q3t(nom)
!   real(8) :: p0(nom),p1(nom),p2(nom),p3(nom)
   real(8), allocatable :: p0(:),p1(:),p2(:),p3(:)
   real(8), allocatable :: q0tt(:),q1tt(:),q2tt(:),q3tt(:)
   real(8), allocatable :: x0t(:),y0t(:),z0t(:)
   real(8), allocatable :: vx0t(:),vy0t(:),vz0t(:)
   real(8), allocatable :: jx0t(:),jy0t(:),jz0t(:)


!isothermal compressibility (betta) set to that of liquid water
!     = 0.007372 dlpoly units
   real(8),parameter :: betta=0.007372d0
   real(8),parameter :: taut = 0.5d0
   real(8),parameter :: taup = 0.5d0
        ! DBG
        real(4) iniTime, deltaTime;
        real(4) iniTimeTot, deltaTimeTot;
        ! DBG

        ! DBG
        iniTimeTot = SECNDS(0.0);
        ! DBG

   allocate(p0(nm),p1(nm),p2(nm),p3(nm))
   allocate(q0tt(nm),q1tt(nm),q2tt(nm),q3tt(nm))
   allocate(x0t(nm),y0t(nm),z0t(nm))
   allocate(vx0t(nm),vy0t(nm),vz0t(nm))
   allocate(jx0t(nm),jy0t(nm),jz0t(nm))

!   volmt = volm
!   boxt = box

   chip = 0.0D0
   chit = 0.0D0

   qtke = 0.d0
   qrke = 0.d0

!$omp parallel DEFAULT(SHARED) private(id,i,j,is,js,lvir)&
!$omp& private(qfx,qfy,qfz,qtx,qty,qtz, tqx,tqy,tqz)&
!$omp& private(qvx,qvy,qvz,qx,qy,qz,qjj,omg,qq,rott)&
!$omp& private(qt0,qt1,qt2,qt3,opx,opy,opz,cconorm,chit)&
!$omp& private(qq0,qq1,qq2,qq3,pp0,pp1,pp2,pp3)

   lvir = 0.0D0
   vircom = 0.d0

! store initial values
!   do j=1,nm

!      q0tt(j) = q0(j)
!      q1tt(j) = q1(j)
!      q2tt(j) = q2(j)
!      q3tt(j) = q3(j)
!      x0t(j)  = x(j)
!      y0t(j)  = y(j)
!      z0t(j)  = z(j)
!      vx0t(j) = vx(j)
!      vy0t(j) = vy(j)
!      vz0t(j) = vz(j)
!      jx0t(j) = jx(j)
!      jy0t(j) = jy(j)
!      jz0t(j) = jz(j)

!   end do



!$omp do schedule(dynamic)
   do j=1,nm

      qfx=0.0
      qfy=0.0
      qfz=0.0
      qtx=0.0
      qty=0.0
      qtz=0.0

      do js=1,ns
         qfx=qfx+fxs(js,j)
         qfy=qfy+fys(js,j)
         qfz=qfz+fzs(js,j)
         qtx=qtx+ys(js,j)*fzs(js,j)-zs(js,j)*fys(js,j)
         qty=qty+zs(js,j)*fxs(js,j)-xs(js,j)*fzs(js,j)
         qtz=qtz+xs(js,j)*fys(js,j)-ys(js,j)*fxs(js,j)

         ! CONVERT FROM SITE-SITE TO MOLECULE-MOLECULE VIRIAL

         lvir = lvir + fxs(js,j)*xs(js,j)+fys(js,j)*ys(js,j)+fzs(js,j)*zs(js,j)
!          vir_test=vir_test+fxs(js,j)*xs(js,j)+fys(js,j)*ys(js,j)+fzs(js,j)*zs(js,j) 

      end do

      fx(j)=qfx
      fy(j)=qfy
      fz(j)=qfz
      tx(j)=qtx
      ty(j)=qty
      tz(j)=qtz

      if (nonadditive_3B) then

         fx(j) = fx(j) + fxfs2(j) + fxfs3(j)
         fy(j) = fy(j) + fyfs2(j) + fyfs3(j)
         fz(j) = fz(j) + fzfs2(j) + fzfs3(j)

         tx(j)= tx(j) + txfs2(j) + txfs3(j)
         ty(j)= ty(j) + tyfs2(j) + tyfs3(j)
         tz(j)= tz(j) + tzfs2(j) + tzfs3(j)

      end if

   end do
!$omp end do
!$omp atomic

   vircom = vircom + lvir

!$omp barrier

!$omp do schedule(dynamic)
   do i=1,nm

      opx = jx(i)
      opy = jy(i)
      opz = jz(i)

      p0(i)=2.0d0*(-q1(i)*opx-q2(i)*opy-q3(i)*opz)
      p1(i)=2.0d0*( q0(i)*opx-q3(i)*opy+q2(i)*opz)
      p2(i)=2.0d0*( q3(i)*opx+q0(i)*opy-q1(i)*opz)
      p3(i)=2.0d0*(-q2(i)*opx+q1(i)*opy+q0(i)*opz)

! current rotation matrix

  end do

!!      call quat_to_rot(q0(i),q1(i),q2(i),q3(i),rott)

   do i=1,nm

      call quat_to_rot(q0(i),q1(i),q2(i),q3(i),rott)

! calculate torque in principal frame

      tqx=tx(i)*rott(1)+ty(i)*rott(4)+tz(i)*rott(7)
      tqy=tx(i)*rott(2)+ty(i)*rott(5)+tz(i)*rott(8)
      tqz=tx(i)*rott(3)+ty(i)*rott(6)+tz(i)*rott(9)

! calculate quaternion torques

      qt0=2.0d0*(-q1(i)*tqx-q2(i)*tqy-q3(i)*tqz)
      qt1=2.0d0*( q0(i)*tqx-q3(i)*tqy+q2(i)*tqz)
      qt2=2.0d0*( q3(i)*tqx+q0(i)*tqy-q1(i)*tqz)
      qt3=2.0d0*(-q2(i)*tqx+q1(i)*tqy+q0(i)*tqz)

! calculate quaternion momenta at start of time step

!      opx = jx(i)
!      opy = jy(i)
!      opz = jz(i)

!      p0(i)=2.0d0*(-q1(i)*opx-q2(i)*opy-q3(i)*opz)
!      p1(i)=2.0d0*( q0(i)*opx-q3(i)*opy+q2(i)*opz)
!      p2(i)=2.0d0*( q3(i)*opx+q0(i)*opy-q1(i)*opz)
!      p3(i)=2.0d0*(-q2(i)*opx+q1(i)*opy+q0(i)*opz)

! update quaternion momenta by 1/2 time step

      p0(i)=p0(i)+qt0*0.5d0*ttt
      p1(i)=p1(i)+qt1*0.5d0*ttt
      p2(i)=p2(i)+qt2*0.5d0*ttt
      p3(i)=p3(i)+qt3*0.5d0*ttt

! update centre of mass velocity by 1/2 time step

      vx(i) = vx(i) + 0.5d0*fx(i)*ttt/totm
      vy(i) = vy(i) + 0.5d0*fy(i)*ttt/totm
      vz(i) = vz(i) + 0.5d0*fz(i)*ttt/totm

   end do
!$omp end do
!$omp barrier

! store initial values
!!   do j=1,nm

!!      q0tt(j) = q0(j)
!!      q1tt(j) = q1(j)
!!      q2tt(j) = q2(j)
!!      q3tt(j) = q3(j)
!!      x0t(j)  = x(j)
!!      y0t(j)  = y(j)
!!      z0t(j)  = z(j)
!!      vx0t(j) = vx(j)
!!      vy0t(j) = vy(j)
!!      vz0t(j) = vz(j)
!!      jx0t(j) = jx(j)
!!      jy0t(j) = jy(j)
!!      jz0t(j) = jz(j)

!!   end do

    ! move centre of mass by full time step

   qvir = vir

   call tkinetic_eng (nm, totm, vx, vy, vz, qtke)
   call rkinetic_eng (nm, jx, jy, jz, qrke)
   call system_temp(nm,qtke,qrke,qsystmp)

!  calculate system pressure
!   boxi=box0
!   volmi = volm0 !boxi**3

   box0=box
   volm0=box0**3.d0

   boxi=box0
   volmi = boxi**3.d0

   qpres = (2.d0*qtke - qvir-vircom)/(3.d0*volmi)
!    qpres = (2.d0*qtke - qvir)/(3.d0*volmi)

!   write(*,*) 'qpres,volmi,qtke,qrke,qtke+qrke,qvir' 
!   write(*,*) qpres,volmi,qtke,qrke,qtke+qrke,qvir 

!  apply Berendsen barostat

   chip = 1.d0+betta*ttt*(qpres-pres0)/taup
   scale = chip**(1.d0/3.d0)

   volmi = volmi*chip
   volm = volmi

   boxi = scale*boxi
!   boxi = volmi**(1.d0/3.d0)
   box = boxi

!   write(*,*) 'chip,volm,betta,taup,qpres,pres0'
!   write(*,*) chip,volm,betta,taup,qpres,pres0

!   volmi = volmi*chip
!   boxi = boxi*scale

!   write(*,*) box,box**3.d0,volm,chip,scale,betta,ttt,pres0,qpres,taup

!$omp do schedule(dynamic)

   do i=1,nm

! move centre of mass by full time step

      x(i)= scale*x(i) + ttt * vx(i)
      y(i)= scale*y(i) + ttt * vy(i)
      z(i)= scale*z(i) + ttt * vz(i)

! rotate rigid groups: nosquish algorithm -update q to full tstep and get new
! rotation matrix


      qq0 = q0(i)
      qq1 = q1(i)
      qq2 = q2(i)
      qq3 = q3(i)
      pp0 = p0(i)
      pp1 = p1(i)
      pp2 = p2(i)
      pp3 = p3(i)

      call nosquish(ttt,qq0,qq1,qq2,qq3,pp0,pp1,pp2,pp3)

      q0(i) = qq0
      q1(i) = qq1
      q2(i) = qq2
      q3(i) = qq3
      p0(i) = pp0
      p1(i) = pp1
      p2(i) = pp2
      p3(i) = pp3

!      call quat_to_rot(q0(i),q1(i),q2(i),q3(i),rott)

   end do

!!      call quat_to_rot(q0(i),q1(i),q2(i),q3(i),rott)

!$omp end do
!$omp barrier

!   call sites

!$omp do schedule(dynamic)
!   do i=1,nm

       ! update COM velocities and angular momentum to 1/2 tstep

!      opx=0.5d0*(-q1(i)*p0(i)+q0(i)*p1(i)+ &
!                 q3(i)*p2(i)-q2(i)*p3(i))
!      opy=0.5d0*(-q2(i)*p0(i)-q3(i)*p1(i)+ &
!                 q0(i)*p2(i)+q1(i)*p3(i))
!      opz=0.5d0*(-q3(i)*p0(i)+q2(i)*p1(i)- &
!                 q1(i)*p2(i)+q0(i)*p3(i))

!      jx(i) = opx
!      jy(i) = opy
!      jz(i) = opz

! move centre of mass by full time step

!      x(i)= scale*x0t(i) + ttt * vx(i)
!      y(i)= scale*y0t(i) + ttt * vy(i)
!      z(i)= scale*z0t(i) + ttt * vz(i)

! restore initial conditions

!      q0(j) = q0tt(j)
!      q1(j) = q1tt(j)
!      q2(j) = q2tt(j)
!      q3(j) = q3tt(j)


!   end do
!$omp end do
!$omp end parallel

! update volume: chip=scale^3

!     volmi=chip*volmi

! end of first stage of velocity verlet algorithm

   call sites

   lvdwpe = vdwpe

        ! DBG
        iniTime = SECNDS(0.0);
        ! DBG
   if (pot_name.eq.'tip4p' ) then
       call lj(x,y,z,xs,ys,zs,box,lvdwpe)
   else if (pot_name.eq.'tip5p' ) then
       call lj(x,y,z,xs,ys,zs,box,lvdwpe)
   else if (pot_name.eq.'sapt5s' ) then
       call sapt5s(x,y,z,xs,ys,zs,box,lvdwpe)
   else if (pot_name.eq.'ccpol5s' ) then
       call ccpol5s(x,y,z,xs,ys,zs,box,lvdwpe)
   else if (pot_name.eq.'ccpol8s' ) then
       call ccpol8s(x,y,z,xs,ys,zs,box,lvdwpe)
   else if (pot_name.eq.'ccpol8s_omo' ) then
       call ccpol8s_omo(x,y,z,xs,ys,zs,box,lvdwpe)
   else if (pot_name.eq.'ccpol23plus' ) then
       call ccpol23plus(x,y,z,xs,ys,zs,box,lvdwpe)
                ! DBG
                deltaTime = SECNDS(iniTime);
!               write(dbgUnit, *) 'P#', MPIrank, '. ccpol23plus time    ',
!               deltaTime;
                ! DBG
   end if

   vdwpe = lvdwpe
   
!   lrcpe = rcpe
!   call coul_pot(x,y,z,xs,ys,zs,box,alpha,lrcpe)
!   rcpe = lrcpe
 
   lscpe = scpe
   call ewald_self_correction(alpha,box,lscpe)
   scpe = lscpe

   lrcpe = rcpe
   call rspace_ewald(x,y,z,xs,ys,zs,box,alpha,lrcpe)
   rcpe = lrcpe

   lkcpe = kcpe
   call kspace_ewald(x,y,z,xs,ys,zs,kmax1,kmax2,kmax3,alpha,box,lkcpe)
   kcpe = lkcpe

!  calling N-Body Iterative Induction model
   if (induction) then
      if(ind_model .eq. 'point_dipole') then

        lindpe = indpe
        call nb_induction(x,y,z,xs,ys,zs,box,lindpe)
        indpe = lindpe

      end if
   end if

! calling 3B-potenital (FS2 and FS3 terms)

  if (nonadditive_3B) then

     lfs2pe = fs2pe
     lfs3pe = fs3pe
     call threebody_term(x,y,z,xs,ys,zs,box,lfs2pe,lfs3pe)
     fs2pe = lfs2pe
     fs3pe = lfs3pe

  end if

   ! second stage of velocity verlet algorithm

!$omp parallel DEFAULT(SHARED) private(id,i,j,is,js,lvir)&
!$omp& private(qfx,qfy,qfz,qtx,qty,qtz, tqx,tqy,tqz)&
!$omp& private(qvx,qvy,qvz,qx,qy,qz,qjj,omg,qq,rott)&
!$omp& private(qt0,qt1,qt2,qt3,opx,opy,opz,cconorm,chit)&
!$omp& private(qq0,qq1,qq2,qq3,pp0,pp1,pp2,pp3)

   lvir = 0.0D0
   vircom = 0.d0

!$omp do schedule(dynamic)
   do j=1,nm

      qfx=0.0
      qfy=0.0
      qfz=0.0
      qtx=0.0
      qty=0.0
      qtz=0.0

      do js=1,ns
         qfx=qfx+fxs(js,j)
         qfy=qfy+fys(js,j)
         qfz=qfz+fzs(js,j)
         qtx=qtx+ys(js,j)*fzs(js,j)-zs(js,j)*fys(js,j)
         qty=qty+zs(js,j)*fxs(js,j)-xs(js,j)*fzs(js,j)
         qtz=qtz+xs(js,j)*fys(js,j)-ys(js,j)*fxs(js,j)

         ! CONVERT FROM SITE-SITE TO MOLECULE-MOLECULE VIRIAL

         lvir = lvir + fxs(js,j)*xs(js,j)+fys(js,j)*ys(js,j)+fzs(js,j)*zs(js,j)
!          vir_test=vir_test+fxs(js,j)*xs(js,j)+fys(js,j)*ys(js,j)+fzs(js,j)*zs(js,j)

      end do

      fx(j)=qfx
      fy(j)=qfy
      fz(j)=qfz
      tx(j)=qtx
      ty(j)=qty
      tz(j)=qtz

      if (nonadditive_3B) then

         fx(j) = fx(j) + fxfs2(j) + fxfs3(j)
         fy(j) = fy(j) + fyfs2(j) + fyfs3(j)
         fz(j) = fz(j) + fzfs2(j) + fzfs3(j)

         tx(j)= tx(j) + txfs2(j) + txfs3(j)
         ty(j)= ty(j) + tyfs2(j) + tyfs3(j)
         tz(j)= tz(j) + tzfs2(j) + tzfs3(j)

      end if

   end do
!$omp end do
!$omp atomic
   vircom = vircom + lvir

!$omp barrier


!$omp do schedule(dynamic)

!!   do i=1,nm

      ! calculate quaternion momenta at start of time step
!!      opx = jx(i)
!!      opy = jy(i)
!!      opz = jz(i)

!!      p0(i)=2.0d0*(-q1(i)*opx-q2(i)*opy-q3(i)*opz)
!!      p1(i)=2.0d0*( q0(i)*opx-q3(i)*opy+q2(i)*opz)
!!      p2(i)=2.0d0*( q3(i)*opx+q0(i)*opy-q1(i)*opz)
!!      p3(i)=2.0d0*(-q2(i)*opx+q1(i)*opy+q0(i)*opz)
!!   end do

!$omp end do
!$omp barrier

!$omp do schedule(dynamic)

   do i=1,nm

! current rotation matrix

      call quat_to_rot(q0(i),q1(i),q2(i),q3(i),rott)

! calculate torque in principal frame

      tqx=tx(i)*rott(1)+ty(i)*rott(4)+tz(i)*rott(7)
      tqy=tx(i)*rott(2)+ty(i)*rott(5)+tz(i)*rott(8)
      tqz=tx(i)*rott(3)+ty(i)*rott(6)+tz(i)*rott(9)

! calculate quaternion torques

      qt0=2.0d0*(-q1(i)*tqx-q2(i)*tqy-q3(i)*tqz)
      qt1=2.0d0*( q0(i)*tqx-q3(i)*tqy+q2(i)*tqz)
      qt2=2.0d0*( q3(i)*tqx+q0(i)*tqy-q1(i)*tqz)
      qt3=2.0d0*(-q2(i)*tqx+q1(i)*tqy+q0(i)*tqz)

! calculate quaternion momenta at half time step

!      opx = jx(i)
!      opy = jy(i)
!      opz = jz(i)

!      p0(i)=2.0d0*(-q1(i)*opx-q2(i)*opy-q3(i)*opz)
!      p1(i)=2.0d0*( q0(i)*opx-q3(i)*opy+q2(i)*opz)
!      p2(i)=2.0d0*( q3(i)*opx+q0(i)*opy-q1(i)*opz)
!      p3(i)=2.0d0*(-q2(i)*opx+q1(i)*opy+q0(i)*opz)

! update quaternion momenta full step

      p0(i)=p0(i)+qt0*0.5d0*ttt
      p1(i)=p1(i)+qt1*0.5d0*ttt
      p2(i)=p2(i)+qt2*0.5d0*ttt
      p3(i)=p3(i)+qt3*0.5d0*ttt

!  update Rgid body angular momentum & COM velocities to full step

!      opx=0.5d0*(-q1(i)*p0(i)+q0(i)*p1(i)+ &
!                 q3(i)*p2(i)-q2(i)*p3(i))
!      opy=0.5d0*(-q2(i)*p0(i)-q3(i)*p1(i)+ &
!                 q0(i)*p2(i)+q1(i)*p3(i))
!      opz=0.5d0*(-q3(i)*p0(i)+q2(i)*p1(i)- &
!                 q1(i)*p2(i)+q0(i)*p3(i))

!      jx(i) = opx
!      jy(i) = opy
!      jz(i) = opz

      vx(i) = vx(i) + 0.5d0*fx(i)*ttt/totm
      vy(i) = vy(i) + 0.5d0*fy(i)*ttt/totm
      vz(i) = vz(i) + 0.5d0*fz(i)*ttt/totm

   end do
!$omp end do
!$omp barrier

!$omp do schedule(dynamic)

   do i=1,nm
!      qq(0) = q0(i)
!      qq(1) = q1(i)
!      qq(2) = q2(i)
!      qq(3) = q3(i)

!      call quat_to_rot(qq,rott)

       ! update COM velocities and angular momentum to 1/2 tstep

      opx=0.5d0*(-q1(i)*p0(i)+q0(i)*p1(i)+ &
                 q3(i)*p2(i)-q2(i)*p3(i))
      opy=0.5d0*(-q2(i)*p0(i)-q3(i)*p1(i)+ &
                 q0(i)*p2(i)+q1(i)*p3(i))
      opz=0.5d0*(-q3(i)*p0(i)+q2(i)*p1(i)- &
                 q1(i)*p2(i)+q0(i)*p3(i))

      jx(i) = opx
      jy(i) = opy
      jz(i) = opz
   end do

!$omp end do
!$omp end parallel

! apply Berendsen thermostat - taut is the relaxation time

   ltke = tke
   call tkinetic_eng (nm, totm, vx, vy, vz, ltke)
   tke = ltke

   lttemp = ttemp
   call trans_temp (nm, tke, lttemp)
   ttemp = lttemp

   stptke = tke/418.4
   stpttp = ttemp

   lrke = rke
   call rkinetic_eng (nm, jx, jy, jz, lrke)
   rke = lrke
   stprke = rke/418.4

   lrtemp = rtemp
   call rot_temp (nm, rke, lrtemp)
   rtemp = lrtemp
   stprtp = rtemp

   lstptmp = stptmp
   call system_temp(nm,tke,rke,lstptmp)
   stptmp = lstptmp

!   write(*,*) temp0,stptmp

   chit = sqrt(1.d0+ ttt/taut*(temp0/stptmp - 1.d0))

!   write(1234,*) 'tke,rke,tke+rke,taut'
!   write(1234,*) tke,rke,tke+rke,taut,vir_test

! thermostat velocities

   do i=1,nm

      vx(i) = chit*vx(i)
      vy(i) = chit*vy(i)
      vz(i) = chit*vz(i)

      jx(i) = chit*jx(i)
      jy(i) = chit*jy(i)
      jz(i) = chit*jz(i)


   end do

   ltke = tke
   call tkinetic_eng (nm, totm, vx, vy, vz, ltke)
   tke = ltke

   lttemp = ttemp
   call trans_temp (nm, tke, lttemp)
   ttemp = lttemp

   stptke = tke/418.4
   stpttp = ttemp

   lrke = rke
   call rkinetic_eng (nm, jx, jy, jz, lrke)
   rke = lrke
   stprke = rke/418.4

   lrtemp = rtemp
   call rot_temp (nm, rke, lrtemp)
   rtemp = lrtemp
   stprtp = rtemp

   lstptmp = stptmp
   call system_temp(nm,tke,rke,lstptmp)
   stptmp = lstptmp

   qvir = vir

!  calculate system pressure
!   box=boxi
!   volmi = boxi**3
!   volm = volmi

   qpres = (2.d0*tke -qvir-vircom)/(3.d0*volm)

!   write(*,*) box,box**3.d0,volm,pres,vir,ltke

   pres = qpres
   stpvolm = volm

   stpvir = qvir/418.4d0+vircom/418.4d0
   stpprs = prsunt*pres

!   write(*,*) 'box,box**3.d0,volm,pres,vir,stptmp'
!   write(*,*) box,box**3.d0,volm,pres,vir,stptmp

   deallocate(p0,p1,p2,p3)
   deallocate(q0tt,q1tt,q2tt,q3tt)
   deallocate(x0t,y0t,z0t)
   deallocate(vx0t,vy0t,vz0t)
   deallocate(jx0t,jy0t,jz0t)

        ! DBG
        deltaTimeTot = SECNDS(iniTimeTot);
!       write(dbgUnit, *) 'P#', MPIrank, '. nvt_b time  ', deltaTimeTot;
        call MPI_BARRIER(MPI_COMM_WORLD, MPIierr);
!       STOP;
        ! DBG

  end subroutine npt_b


  subroutine nosquish(tstep,qq0,qq1,qq2,qq3,pp0,pp1,pp2,pp3)
    implicit none
	real(8) :: tstep
    real(8) :: qq0,qq1,qq2,qq3
    real(8) :: pp0,pp1,pp2,pp3
    integer :: m
    real(8) zetax,zetay,zetaz,cs,sn,trstep

    integer, parameter :: mrot=10
    real(8), parameter :: ov4=0.25d0
    real(8), parameter :: ov8=0.125d0


    real(8) qn1(0:3),pq2(0:3),qn2(0:3),pq3(0:3)
    real(8) qn3(0:3),pq4(0:3)

!    rotational time step

    trstep=tstep/dble(mrot)

    do m=1,mrot

       zetaz=ov8/pmi(3)*trstep* &
         (-pp0*qq3+pp1*qq2- &
         pp2*qq1+pp3*qq0)

       cs=cos(zetaz)
       sn=sin(zetaz)

       qn1(0)=cs*qq0-sn*qq3
       qn1(1)=cs*qq1+sn*qq2
       qn1(2)=cs*qq2-sn*qq1
       qn1(3)=cs*qq3+sn*qq0
       pq2(0)=cs*pp0-sn*pp3
       pq2(1)=cs*pp1+sn*pp2
       pq2(2)=cs*pp2-sn*pp1
       pq2(3)=cs*pp3+sn*pp0

       zetay=ov8/pmi(2)*trstep* &
         (-pq2(0)*qn1(2)-pq2(1)*qn1(3)+ &
         pq2(2)*qn1(0)+pq2(3)*qn1(1))
       cs=cos(zetay)
       sn=sin(zetay)
       qn2(0)=cs*qn1(0)-sn*qn1(2)
       qn2(1)=cs*qn1(1)-sn*qn1(3)
       qn2(2)=cs*qn1(2)+sn*qn1(0)
       qn2(3)=cs*qn1(3)+sn*qn1(1)
       pq3(0)=cs*pq2(0)-sn*pq2(2)
       pq3(1)=cs*pq2(1)-sn*pq2(3)
       pq3(2)=cs*pq2(2)+sn*pq2(0)
       pq3(3)=cs*pq2(3)+sn*pq2(1)

       zetax=ov4/pmi(1)*trstep* &
         (-pq3(0)*qn2(1)+pq3(1)*qn2(0)+ &
         pq3(2)*qn2(3)-pq3(3)*qn2(2))
        cs=cos(zetax)
        sn=sin(zetax)
        qn3(0)=cs*qn2(0)-sn*qn2(1)
        qn3(1)=cs*qn2(1)+sn*qn2(0)
        qn3(2)=cs*qn2(2)+sn*qn2(3)
        qn3(3)=cs*qn2(3)-sn*qn2(2)
        pq4(0)=cs*pq3(0)-sn*pq3(1)
        pq4(1)=cs*pq3(1)+sn*pq3(0)
        pq4(2)=cs*pq3(2)+sn*pq3(3)
        pq4(3)=cs*pq3(3)-sn*pq3(2)

       zetay=ov8/pmi(2)*trstep* &
         (-pq4(0)*qn3(2)-pq4(1)*qn3(3)+ &
         pq4(2)*qn3(0)+pq4(3)*qn3(1))
        cs=cos(zetay)
        sn=sin(zetay)
        qn2(0)=cs*qn3(0)-sn*qn3(2)
        qn2(1)=cs*qn3(1)-sn*qn3(3)
        qn2(2)=cs*qn3(2)+sn*qn3(0)
        qn2(3)=cs*qn3(3)+sn*qn3(1)
        pq3(0)=cs*pq4(0)-sn*pq4(2)
        pq3(1)=cs*pq4(1)-sn*pq4(3)
        pq3(2)=cs*pq4(2)+sn*pq4(0)
        pq3(3)=cs*pq4(3)+sn*pq4(1)

       zetaz=ov8/pmi(3)*trstep* &
         (-pq3(0)*qn2(3)+pq3(1)*qn2(2)- &
         pq3(2)*qn2(1)+pq3(3)*qn2(0))
        cs=cos(zetaz)
        sn=sin(zetaz)
        qq0=cs*qn2(0)-sn*qn2(3)
        qq1=cs*qn2(1)+sn*qn2(2)
        qq2=cs*qn2(2)-sn*qn2(1)
        qq3=cs*qn2(3)+sn*qn2(0)
        pp0=cs*pq3(0)-sn*pq3(3)
        pp1=cs*pq3(1)+sn*pq3(2)
        pp2=cs*pq3(2)-sn*pq3(1)
        pp3=cs*pq3(3)+sn*pq3(0)

    end do

  end subroutine nosquish   


  subroutine nvtqscl(ttt,qmass,dof,taut,chit,conint)
   implicit none
   real(8) :: ttt,chit,qmass,taut,dof,conint

   integer :: i
   real(8) :: qtke,qrke,qsystmp
   real(8) :: scale

!   calculate kinetic energy and inst temp
   call tkinetic_eng (nm, totm, vx, vy, vz, qtke)
   call rkinetic_eng (nm, jx, jy, jz, qrke)
   call system_temp(nm,qtke,qrke,qsystmp)  

!  update chit to 1/2 step
   chit=chit + ttt*dof*boltz*(qsystmp-temp0)/qmass
   
!  thermostat scale parameter
   scale=exp(-ttt*chit)
      
!  thermostat rigid body velocities

   do i=1,nm
 
      vx(i) = scale*vx(i)
      vy(i) = scale*vy(i)
      vz(i) = scale*vz(i) 

      jx(i) = scale*jx(i)
      jy(i) = scale*jy(i)
      jz(i) = scale*jz(i)

   end do

!  scale kinetic energy

   qtke = qtke*scale**2
   qrke = qrke*scale**2
   call system_temp(nm,qtke,qrke,qsystmp)

!  update chi to full step

   conint=conint+ttt*chit*qmass/taut**2 

!  update chit to full step
   chit=chit + ttt*dof*boltz*(qsystmp-temp0)/qmass

  end subroutine nvtqscl


  subroutine nptqscl_t(ttt,qmass,pmass,dof,taut,chip,chit,conint)
   implicit none
   real(8) :: ttt,chit,chip,pmass,qmass,taut,dof,conint

   integer :: i
   real(8) :: qtke,qrke,qsystmp
   real(8) :: scale

!   calculate kinetic energy and inst temp
   call tkinetic_eng (nm, totm, vx, vy, vz, qtke)
   call rkinetic_eng (nm, jx, jy, jz, qrke)
   call system_temp(nm,qtke,qrke,qsystmp)  

!  update chit to 1/2 step
   chit=chit + 0.5d0*ttt*dof*boltz*(qsystmp-temp0)/qmass + &
               0.5d0*ttt*(pmass*chip**2.d0-boltz*temp0)/qmass
   
!  thermostat scale parameter
   scale=exp(-ttt*chit)
      
!  thermostat rigid body velocities

   do i=1,nm
 
      vx(i) = scale*vx(i)
      vy(i) = scale*vy(i)
      vz(i) = scale*vz(i) 

      jx(i) = scale*jx(i)
      jy(i) = scale*jy(i)
      jz(i) = scale*jz(i)

   end do

!  scale kinetic energy

   qtke = qtke*scale**2.d0
   qrke = qrke*scale**2.d0
   call system_temp(nm,qtke,qrke,qsystmp)

!  update chi to full step

   conint=conint+ttt*chit*(qmass/taut**2+boltz*temp0) 

!  update chit to full step
   chit=chit + 0.5d0*ttt*dof*boltz*(qsystmp-temp0)/qmass + &
               0.5d0*ttt*(pmass*chip**2.d0-boltz*temp0)/qmass

  end subroutine nptqscl_t

  subroutine nptqscl_p(ttt,boxi,pmass,dof,chit,chip,qvir,qpres,conint)
   implicit none
   real(8) :: ttt,chit,chip,pmass,taut,dof,conint

   integer :: i
   real(8) :: qtke,qrke,qsystmp,boxi
   real(8) :: scale,volmi,qpres,qvir

!   calculate kinetic energy and inst temp

   call tkinetic_eng (nm, totm, vx, vy, vz, qtke)

!   call rkinetic_eng (nm, jx, jy, jz, qrke)
!   call system_temp(nm,qtke,qrke,qsystmp)  

   volmi = boxi**3

   qpres = (2.d0*qtke - qvir)/(3.d0*volmi) 

!  update chit to 1/2 step
   chip=chip + 0.5d0*ttt*(3.d0*volmi*(qpres-pres0)/pmass - chip*chit)
   
!  barostat scale parameter
   scale=exp(-ttt*chip)
      
!  thermostat rigid body velocities

   do i=1,nm
 
      vx(i) = scale*vx(i)
      vy(i) = scale*vy(i)
      vz(i) = scale*vz(i) 

!      jx(i) = scale*jx(i)
!      jy(i) = scale*jy(i)
!      jz(i) = scale*jz(i)

   end do

!  scale kinetic energy

   qtke = qtke*scale**2.d0

!   qrke = qrke*scale**2                   ! check

!   call system_temp(nm,qtke,qrke,qsystmp) ! check

!   volmi=boxi**3.d0

   volmi=volmi*exp(3.d0*ttt*chip)  
   
!  update chip to full step
   boxi = volmi**1.d0/3.d0

!!   boxi=boxi*exp(3.d0*ttt*chip)
!!   volmi=boxi**3.d0

   qpres = (2.d0*qtke - qvir)/(3.d0*volmi)

!  update chit to 1/2 step
   chip=chip + 0.5d0*ttt*(3.d0*volmi*(qpres-pres0)/pmass - chip*chit)

  end subroutine nptqscl_p

end module evolve

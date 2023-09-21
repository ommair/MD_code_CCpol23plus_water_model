module utility

  use, intrinsic :: iso_fortran_env
  use common_arrays

  implicit none

  private

  public :: equilibration
  public :: jacobi, gauss, tkinetic_eng, trans_temp, rescale_velocities
  public :: system_temp, rot_temp,rkinetic_eng,site_to_com_vel
  public :: tt,dtt,scalp !,flr,dflr
  public :: matvec1
  public :: equilibrationGRU,system_com,system_mass,system_vcom
  
  ! GRU
  real(8) :: duni_DefValue = 0.0D0;
  ! GRU

contains

! scalar product
  function scalp(a,b)  
        implicit none
        
        real(8) ::  a(3),b(3)
        real(8) :: scalp

        scalp = a(1)*b(1) + a(2)*b(2) + a(3)*b(3)
  end function scalp


! Multiply vetor v by the matrix a. Store in u.

  subroutine matvec1(a,v,u)
    implicit none
    integer :: i,j,nnn
    real(8) ::  a(3,3),v(3),u(3)

    nnn = 3
    do i=1,nnn
       u(i) = 0.d0
       do j=1,nnn
          u(i) = u(i) + a(i,j)*v(j)
       end do
    end do
  
  end subroutine matvec1


! compute the factorial of n
!
  function fact(nmax)
    implicit none
    integer :: nmax
    integer :: i
    real(8) :: fact
    fact = 1.d0
    do i=1,nmax
    fact = fact*i
    end do
    
!    write(200,*)'fact', fact    
  end function fact

  ! tang-toennies damoing function
!  function tt(n, delta, r)
!    implicit none
!    integer :: i
!    integer, intent(in) :: n
!    real(8), intent(in) :: delta, r
!    real(8) :: tt
!    real(8) :: rscaled, sum, ffact, rtoi

!    rscaled = delta * r
!    rtoi = 1.0
!    ffact = 1.0
!    sum = 1.0
!    do i=1, n
!       ffact = ffact * i
!       rtoi = rtoi*rscaled
!       sum = sum + rtoi/ffact
!    end do
!    tt = 1.0 - exp(-delta*r)*sum

!  end function tt

  function tt(n, delta, r)
    implicit none
    integer :: i,ncn
    integer, intent(in) :: n
    real(8), intent(in) :: delta, r
    real(8) :: tt
    real(8) :: term,ssum,deltar

    deltar = delta * r
    term = 1.0d0
    ssum = 1.0d0
    ncn = n
    do i=1, ncn
       term = term*deltar/i
       ssum = ssum + term
    end do
    tt = 1.0 - exp(-deltar)*ssum

  end function tt

  function dtt(n, delta, r)
    implicit none
    integer :: i,ncn
    integer, intent(in) :: n
    real(8), intent(in) :: delta, r
    real(8) :: dtt
    real(8) :: delr, ffact

    ffact = 1.0d0
    ncn = n
    do i=1, ncn
       ffact = ffact * i
    end do
    delr = delta * r
    dtt = delta*exp(-delr)*delr**n/ffact

  end function dtt

  function flr(rab,rct)

    implicit none

    real(8), parameter :: epsprim = 1.d-10
    real(8), parameter :: gamdmp = bohr2a*10.d0
    real(8), parameter :: delta_rb = log(1.d0/epsprim)/log(2.d0)
    real(8) :: r0dmp,flr,rab,rct
    real(8) :: pom,pomij

    r0dmp = rct*bohr2a**(-1.d0)

    pom = exp(gamdmp*(rab-r0dmp))
    pomij = 1.d0/(1.d0+pom)

    flr = pomij**delta_rb
  end function flr

  function dflr(rab,rct)

    implicit none

    real(8), parameter :: epsprim = 1.d-10
    real(8), parameter :: gamdmp = bohr2a*10.d0
    real(8), parameter :: delta_rb = log(1.d0/epsprim)/log(2.d0)
    real(8) :: r0dmp,fflr,dflr,rab,rct
    real(8) :: pom,pomij

    r0dmp = rct*bohr2a**(-1.d0)

    pom = exp(gamdmp*(rab-r0dmp))
    pomij = 1.d0/(1.d0+pom)
    fflr = pomij**delta_rb

    dflr = -gamdmp*delta_rb*pom*pomij*fflr
  end function dflr

  subroutine tkinetic_eng (n, ms, va, vb, vc, ke)
    implicit none
    real(8), dimension(:) :: va, vb, vc  ! dummy variables for velocities
    real(8) :: ms   ! dummy mass  variable
    integer :: n ! dummy n atoms
    real(8) :: ke  ! dummy kinetic energy
    integer :: j

    ke = 0.d0
    do j = 1, n
       ke = ke + (va(j)**2 + vb(j)**2 + vc(j)**2)

!      write(102,*) step,n,ms,va(j),vb(j),vc(j),0.5*ms*(va(j)**2 + vb(j)**2 + vc(j)**2)
    end do 
      ke = 0.5*ms*ke
  end subroutine tkinetic_eng

  subroutine rkinetic_eng (n, ja, jb, jc, ke)
    implicit none
    real(8), dimension(:) :: ja, jb, jc  ! dummy variables for velocities
!    real(8) :: m   ! dummy mass  variable
    integer :: n ! dummy n atoms
    real(8) :: ke  ! dummy kinetic energy
    integer :: j

    ke = 0.d0
    do j = 1, n
       ke = ke + (ja(j)**2/pmi(1) + jb(j)**2/pmi(2) + jc(j)**2/pmi(3))
    end do
      ke = 0.5*ke
  end subroutine rkinetic_eng

  subroutine trans_temp (n, tkke, ttmp)
    implicit none

    real(8) :: tkke  ! dummy kinetic energy
!    real(8) :: kb   ! boltzman const
    real(8) :: ttmp ! dummy temp variable
    integer :: n
    real(8) :: dof

!    boltz = 8.31451115d-1 ! conversion kelvin to eo/k = (10j/(mol*k)
    dof = real(3.d0*n-3.d0)
!    dof = real(3.d0*n)

    ttmp = (2.d0 * tkke) /(dof*boltz)

  end subroutine trans_temp

  subroutine rot_temp (n, rkke, rtmp)
    implicit none

    real(8) :: rkke  ! dummy kinetic energy
!    real(8) :: kb   ! boltzman const
    real(8) :: rtmp ! dummy temp variable
    integer :: n
    real(8) :: dof

!    boltz = 8.31451115d-1 ! conversion kelvin to eo/k = (10j/(mol*k)
    dof = real(3.d0*n)

    rtmp = (2.d0 * rkke) /(dof*boltz)

  end subroutine rot_temp

  subroutine system_temp(n,tkke,rkke,systmp)
    implicit none
    integer :: n
    real(8) :: tkke
    real(8) :: rkke
    real(8) :: systmp 

    systmp = 2.d0*(tkke+rkke)/(3.d0*boltz*(2.d0*n-1.d0))

  end subroutine system_temp

  subroutine rescale_velocities (n, temp0, temp, e1, e2, e3)
    implicit none

    integer :: i
    integer :: n         ! dummy n atoms
    real(8) :: temp0, temp  ! dummy temp variables
    real(8), dimension(:) :: e1, e2, e3  ! dummy vleocity variables
    real(8) :: f

    f = sqrt(temp0/temp)

    do i=1, n

       e1(i) = f * e1(i)
       e2(i) = f * e2(i)
       e3(i) = f * e3(i)

    end do
  end subroutine rescale_velocities

  subroutine equilibration(n,tmp0,ttmp,e1,e2,e3,rtmp,w1,w2,w3)

    implicit none 
  
    real(8) :: tmp0  ! desired temp
    real(8) :: tkke, ttmp  ! transl ke and temp
    real(8) :: rkke, rtmp  ! rot ke and temp
    real(8), dimension(:) :: e1, e2, e3 ! transl velocities
    real(8), dimension(:) :: w1, w2, w3 ! rot velocities
    integer :: n

    call rescale_velocities (n,tmp0,ttmp,e1,e2,e3) ! linear vel rescaling
    call rescale_velocities (n,tmp0,rtmp,w1,w2,w3) ! angular vel rescaling

  end subroutine equilibration

  subroutine system_mass
 
    implicit none

    integer :: i,j,n,is,js   

    sysmass=0.d0

    do i=1,nm
    do is=1,ns
       if (massites(is).ne.0.d0) then      
         sysmass = sysmass + massites(is)    ! total mass of the system
       end if
    end do
    end do

  end subroutine system_mass

  subroutine system_com

    implicit none

    integer :: i,j,n,is,js

    xcom=0.d0
    ycom=0.d0
    zcom=0.d0

    do i=1,nm
    do is=1,ns  ! loop on sites
          xcom = xcom + massites(is)*(x(i)+xs(is,i))/sysmass
          ycom = ycom + massites(is)*(y(i)+ys(is,i))/sysmass
          zcom = zcom + massites(is)*(z(i)+zs(is,i))/sysmass
    end do
    end do 

  end subroutine system_com

  subroutine system_vcom

    implicit none

    integer :: i,j,n,is,js

    vxcom=0.d0
    vycom=0.d0
    vzcom=0.d0

    do i=1,nm
!    do is=1,ns  ! loop on sites
          vxcom = vxcom + totm*vx(i)/sysmass
          vycom = vycom + totm*vy(i)/sysmass
          vzcom = vzcom + totm*vz(i)/sysmass
!    end do
    end do

!    write(*,*) nm*totm,sysmass,vxcom,vycom,vzcom

  end subroutine system_vcom
  
  subroutine equilibrationGRU
    implicit none;
  
!    real(8) :: tmp0  ! desired temp
!    real(8) :: tkke, ttmp  ! transl ke and temp
!    real(8) :: rkke, rtmp  ! rot ke and temp
!    real(8), dimension(:) :: e1, e2, e3 ! transl velocities
!    real(8), dimension(:) :: w1, w2, w3 ! rot velocities
!    integer :: n
!n,tmp0,ttmp,e1,e2,e3,rtmp,w1,w2,w3
!nm,temp0,stpttp,vx,vy,vy,stprtp,jx,jy,jz

    call rescale_velocities (nm,temp0,stpttp,vx,vy,vz) ! linear vel rescaling
    call rescale_velocities (nm,temp0,stprtp,jx,jy,jz) ! angular vel rescaling

  end subroutine equilibrationGRU

  subroutine equilibration_b(n,tmp0,ttmp,e1,e2,e3,rtmp,w1,w2,w3)

    implicit none

    real(8) :: tmp0  ! desired temp
    real(8) :: tkke, ttmp  ! transl ke and temp
    real(8) :: rkke, rtmp  ! rot ke and temp
    real(8), dimension(:) :: e1, e2, e3 ! transl velocities
    real(8), dimension(:) :: w1, w2, w3 ! rot velocities
    integer :: n

    call rescale_velocities (n,tmp0,ttmp,e1,e2,e3) ! linear vel rescaling
    call rescale_velocities (n,tmp0,rtmp,w1,w2,w3) ! angular vel rescaling

  end subroutine equilibration_b

  subroutine site_to_com_vel(gvxx,gvyy,gvzz,omxx,omyy,omzz)

    implicit none
  
    integer :: i,j,is,js
    !real(8), dimension(nsite,nom) :: vxx,vyy,vzz  ! site velocities 
    real(8), dimension(nom) :: gvxx,gvyy,gvzz  ! com velocities
    real(8), dimension(nom) :: omxx,omyy,omzz  ! angular momentum
    
    real(8) :: rotm(1:9) 
    real(8) :: totmass    

   totmass = 0.d0
    do i = 1,nm
       gvxx(i)=0.d0
       gvyy(i)=0.d0
       gvzz(i)=0.d0
       omxx(i)=0.d0
       omyy(i)=0.d0
       omzz(i)=0.d0
    end do

    ! centre of mass momentum 
    do i = 1,nm
       do is=1,ns
          gvxx(i)=gvxx(i)+massites(is)*vxs(is,i)
          gvyy(i)=gvyy(i)+massites(is)*vys(is,i)
          gvzz(i)=gvzz(i)+massites(is)*vzs(is,i)
          totmass = totmass + massites(is)
       end do 
    end do

!    centre of mass velocity
    do i = 1,nm
       gvxx(i)=gvxx(i)/totmass
       gvyy(i)=gvyy(i)/totmass
       gvzz(i)=gvzz(i)/totmass
    end do

    !xs(js,j) 
    do i = 1,nm
       do is=1,ns
        omxx(i)=omxx(i)+massites(is)*(ys(is,i)*vzs(is,i)-zs(is,i)*vys(is,i))     
        omyy(i)=omyy(i)+massites(is)*(zs(is,i)*vxs(is,i)-xs(is,i)*vzs(is,i))
        omzz(i)=omzz(i)+massites(is)*(xs(is,i)*vys(is,i)-ys(is,i)*vxs(is,i)) 
       end do
    end do

    do i = 1,nm

       rotm(1)=q0(i)**2+q1(i)**2-q2(i)**2-q3(i)**2
       rotm(2)=2.d0*(q1(i)*q2(i)-q0(i)*q3(i))
       rotm(3)=2.d0*(q1(i)*q3(i)+q0(i)*q2(i))
       rotm(4)=2.d0*(q1(i)*q2(i)+q0(i)*q3(i))
       rotm(5)=q0(i)**2-q1(i)**2+q2(i)**2-q3(i)**2
       rotm(6)=2.d0*(q2(i)*q3(i)-q0(i)*q1(i))
       rotm(7)=2.d0*(q1(i)*q3(i)-q0(i)*q2(i))
       rotm(8)=2.d0*(q2(i)*q3(i)+q0(i)*q1(i))
       rotm(9)=q0(i)**2-q1(i)**2-q2(i)**2+q3(i)**2

      !  angular momentum in body fixed frame

!       omxx(i)=rotm(1)*omxx(i)+rotm(4)*omyy(i)+rotm(7)*omzz(i)
!       omyy(i)=rotm(2)*omxx(i)+rotm(5)*omyy(i)+rotm(8)*omzz(i)
!       omzz(i)=rotm(3)*omxx(i)+rotm(6)*omyy(i)+rotm(9)*omzz(i)

!       write(*,*) rotm(1),rotm(2),rotm(3),omxx(i),omyy(i),omzz(i)  

    end do
 
  end subroutine site_to_com_vel 

   subroutine jacobi(a,v,n)

!c***********************************************************************
!c
!c     diagonalisation of real symmetric matices by jacobi method
!c
!c     input parameters:
!c
!c     a(n,n) is the matrix to be diagonalised
!c     v(n,n) is the eigenvector matrix
!c     n   is the dimension of the matrices
!c
!c     jacobi processes lower triangle only (upper triangle unchanged)
!c
!c     variable rho sets absolute tolerance on convergence
!c     variable tes is a moving tolerance that diminishes
!c     on each pass until at true convergence tes<rho
!c
!c     author w.smith 1993
!c
!c***********************************************************************

      implicit none

      logical :: pass
      integer :: n,i,j,k
      real(8) :: a,v,rho,tes,scl,v1,v2,v3,u,omg,s,c,tem

      dimension a(n,n),v(n,n)

      rho=1.0d-16
      tes=0.0d0
      scl=0.0d0

!c     initialize eigenvectors

      do i=1,n
        do j=1,n
          v(i,j)=0.0d0
        end do
        v(i,i)=1.0d0
      end do

!c     rescale matrix for optimal accuracy

      do i=1,n
        if(abs(a(i,i)).gt.scl) scl=abs(a(i,i))
      end do
      do i=1,n
        do j=1,i
          a(i,j)=a(i,j)/scl
        end do
      end do

!c     set initial value of moving tolerance

      do i=2,n
        do j=1,i-1
          tes=tes+2.0d0*a(i,j)*a(i,j)
        enddo
      enddo
      tes=sqrt(tes)

!c     recycle until absolute tolerance satisfied

      do while(tes.gt.rho)

        tes=tes/dble(n)
        if(tes.lt.rho)tes=rho

!c     jacobi diagonalisation

        pass=.true.

!c     recycle until moving tolerance satisfied

        do while(pass)

          pass=.false.

          do i=2,n

            do j=1,i-1

              if(abs(a(i,j)).ge.tes)then
                pass=.true.
                v1=a(j,j)
                v2=a(i,j)
                v3=a(i,i)
                u=0.5d0*(v1-v3)
                if(abs(u).lt.rho)then
                  omg=-1.0d0
                else
                  omg=-v2/sqrt(v2*v2+u*u)
                  if(u.lt.0.0d0)omg=-omg
                endif
                s=omg/sqrt(2.0d0*(1.0d0+sqrt(1.0d0-omg*omg)))
                c=sqrt(1.0d0-s*s)
                do k=1,n
                  if(k.ge.i)then
                    tem=a(k,j)*c-a(k,i)*s
                    a(k,i)=a(k,j)*s+a(k,i)*c
                    a(k,j)=tem
                  else if(k.lt.j)then
                    tem=a(j,k)*c-a(i,k)*s
                    a(i,k)=a(j,k)*s+a(i,k)*c
                    a(j,k)=tem
                  else
                    tem=a(k,j)*c-a(i,k)*s
                    a(i,k)=a(k,j)*s+a(i,k)*c
                    a(k,j)=tem
                  endif
                  tem=v(k,j)*c-v(k,i)*s
                  v(k,i)=v(k,j)*s+v(k,i)*c
                  v(k,j)=tem
                enddo
                a(j,j)=v1*c*c+v3*s*s-2.0d0*v2*s*c
                a(i,i)=v1*s*s+v3*c*c+2.0d0*v2*s*c
                a(i,j)=(v1-v3)*s*c+v2*(c*c-s*s)
              endif

            enddo

          enddo

        enddo

      enddo

!c     rescale matrix

      do i=1,n
        do j=1,i
          a(i,j)=scl*a(i,j)
        enddo
      enddo

      return
      end subroutine jacobi

      function duni()

!c*********************************************************************
!c     
!c     dl_poly random number generator based on the universal
!c     random number generator of marsaglia, zaman and tsang
!c     (stats and prob. lett. 8 (1990) 35-39.) it must be
!c     called once to initialise parameters u,c,cd,cm
!c     
!c     copyright daresbury laboratory 1992
!c     author -  w.smith         july 1992
!c     
!c*********************************************************************

      implicit none

      logical new
      integer ir,jr,i,j,k,l,m,ii,jj
      real(4) s,t,u,c,cd,cm,uni
      real(8) duni
      dimension u(97)
      save u,c,cd,cm,uni,ir,jr,new
      data new/.true./

      if(new)then

!c     initial values of i,j,k must be in range 1 to 178 (not all 1)
!c     initial value of l must be in range 0 to 168.

        i=12
        j=34
        k=56
        l=78
!c     
        ir=97
        jr=33
        new=.false.

        do 200 ii=1,97
          s=0.0
          t=0.5
          do 100 jj=1,24
            m=mod(mod(i*j,179)*k,179)
            i=j
            j=k
            k=m
            l=mod(53*l+1,169)
            if(mod(l*m,64).ge.32)s=s+t
            t=0.5*t
  100     continue
          u(ii)=s
  200   continue
        c =  362436.0/16777216.0
        cd= 7654321.0/16777216.0
        cm=16777213.0/16777216.0
      else

!c     calculate random number
        uni=u(ir)-u(jr)
        if(uni.lt.0.0)uni=uni+1.0
        u(ir)=uni
        ir=ir-1
        if(ir.eq.0)ir=97
        jr=jr-1
        if(jr.eq.0)jr=97
        c=c-cd
        if(c.lt.0.0)c=c+cm
        uni=uni-c
        if(uni.lt.0.0)uni=uni+1.0
!		duni=dble(uni)
        duni_DefValue = dble(uni);
      endif
      duni=duni_DefValue;
      
      return
      end function duni


     subroutine gauss(nnn,vxx,vyy,vzz)

!c*********************************************************************
!c     
!c     dl_poly subroutine for constructing velocity arrays
!c     with a gaussian distribution of unit variance.
!c     
!c     based on the box-muller method
!c     
!c     note - this version uses a universal random number 
!c     generator, which generates pseudo-random numbers between
!c     0 and 1. it is based on the algorithm of marsaglia, zaman
!c     and tsang in: stats and prob. lett. 8 (1990) 35-39.
!c     
!c     copyright daresbury laboratory 2007
!c     author - w. smith         nov  2007
!c     
!c*********************************************************************
      
!      use setup_module
      
      implicit none

      integer nnn,ii
      real(8) vxx,vyy,vzz,rrr,rr1,rr2
      
      dimension vxx(nnn),vyy(nnn),vzz(nnn)
      
!c     initialise random number generator
      
      rrr=duni()
      
!c     calculate gaussian random numbers
      
      do ii=1,2*(nnn/2),2
        
        rr1=sqrt(-2.d0*log(duni()))
        rr2=2.d0*pi*duni()
        vxx(ii)=rr1*cos(rr2)
        vxx(ii+1)=rr1*sin(rr2)

        rr1=sqrt(-2.d0*log(duni()))
        rr2=2.d0*pi*duni()
        vyy(ii)=rr1*cos(rr2)
        vyy(ii+1)=rr1*sin(rr2)

        rr1=sqrt(-2.d0*log(duni()))
        rr2=2.d0*pi*duni()
        vzz(ii)=rr1*cos(rr2)
        vzz(ii+1)=rr1*sin(rr2)
        
      enddo
      if(mod(nnn,2).ne.0)then
        
        rr1=sqrt(-2.d0*log(duni()))
        rr2=2.d0*pi*duni()
        vxx(nnn)=rr1*cos(rr2)
        vyy(nnn)=rr1*sin(rr2)
        rr1=sqrt(-2.d0*log(duni()))
        rr2=2.d0*pi*duni()
        vzz(nnn)=rr1*cos(rr2)
        
      endif
      
      return
      end subroutine gauss
 
! not complete !!! please do it later on 
  subroutine rdf

    implicit none

    mxrdf = max(128,int(rcut/drdf))
    ! default bin width
    delrdf = rcut/real(mxrdf) 

  end subroutine rdf

end module utility  

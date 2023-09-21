module molecular_sites
  use, intrinsic :: iso_fortran_env
  use common_arrays

  implicit none
  private

  public :: quat_to_rot
  public :: sites, transform
  

contains

  subroutine sites
    implicit none
  
    integer :: i,j,is,js
    real(8) :: qp(3)
     
!    use quaternions to set site coordinates relative to centre-of-mass
!$omp parallel DEFAULT(SHARED)&
!$omp& private(i,j,is,js,qp)
!$omp do schedule(dynamic)
      do  j=1,nm

!          write(100,*)j, ' molecule'
          do  js=1,ns

              qp(1)=pcoord(1,js)
              qp(2)=pcoord(2,js)
              qp(3)=pcoord(3,js)
              call transform(q0(j),q1(j),q2(j),q3(j),qp,1)
              xs(js,j)=qp(1)
              ys(js,j)=qp(2)
              zs(js,j)=qp(3)
              chgs(js,j)=charge(js)
              apol(js,j)=polarizability(js)*bohr2a**(3.d0)  ! in DLP units

!              apol(js,j)=polarizability(js)*bohr2a**(3.d0)/r4pie0    
!              apol(js,j)=polarizability(js)

          end do
      end do 
!$omp end do
!$omp end parallel
 
  end subroutine sites 

  subroutine transform(qq0,qq1,qq2,qq3,rr,iway)
   
    implicit none
   
    real(8), intent(in) :: qq0,qq1,qq2,qq3
    real(8) :: rott(1:9)
    real(8) :: qrr(3)
    real(8), intent(out) :: rr(3)
    integer :: iway
    
! components of rotation matrix from quaternions
  
    call quat_to_rot(qq0,qq1,qq2,qq3,rott)


!* case iway=1, from  principle axis to box frame 
!   
    if(iway.eq.1) then

      qrr(1) = rr(1)*rott(1) + rr(2)*rott(2) + rr(3)*rott(3)
      qrr(2) = rr(1)*rott(4) + rr(2)*rott(5) + rr(3)*rott(6)
      qrr(3) = rr(1)*rott(7) + rr(2)*rott(8) + rr(3)*rott(9)

      rr(1)=qrr(1)
      rr(2)=qrr(2)
      rr(3)=qrr(3)

    endif

!* case iway=2,  from box frame to  principle axis

    if (iway.eq.2) then

      qrr(1)=rr(1)*rott(1) +rr(2)*rott(4) +rr(3)*rott(7)
      qrr(2)=rr(1)*rott(2) +rr(2)*rott(5) +rr(3)*rott(8)
      qrr(3)=rr(1)*rott(3) +rr(2)*rott(6) +rr(3)*rott(9)

      rr(1)=qrr(1)
      rr(2)=qrr(2)
      rr(3)=qrr(3)

    endif
 

  end subroutine transform

  subroutine quat_to_rot(qq0,qq1,qq2,qq3,rott)

    implicit none
    real(8), intent(in) :: qq0,qq1,qq2,qq3
    real(8), intent(out) :: rott(1:9)

    rott(1) = qq0**2 + qq1**2 - qq2**2 - qq3**2
    rott(2) = 2.0*(qq1*qq2 - qq0*qq3)
    rott(3) = 2.0*(qq1*qq3 + qq0*qq2)
    rott(4) = 2.0*(qq1*qq2 + qq0*qq3)
    rott(5) = qq0**2 - qq1**2 + qq2**2 - qq3**2
    rott(6) = 2.0*(qq2*qq3 - qq0*qq1)
    rott(7) = 2.0*(qq1*qq3 - qq0*qq2)
    rott(8) = 2.0*(qq2*qq3 + qq0*qq1)
    rott(9) = qq0**2 - qq1**2 - qq2**2 + qq3**2

  end subroutine quat_to_rot
 
end module molecular_sites

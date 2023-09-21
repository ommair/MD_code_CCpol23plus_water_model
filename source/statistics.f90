module statistics

  use, intrinsic :: iso_fortran_env
  use common_arrays

  implicit none
  private
  public :: zero,tot_sum,averages_fluctuations

contains

  subroutine zero

    implicit none

    integer :: jr,jd
    integer :: i,j,is,js

!    do j=1,nm
!       ddx(j) = 0.d0
!       ddy(j) = 0.d0
!       ddz(j) = 0.d0
!    end do

    numrdf = 0
      

    sumvdwpe = 0.d0
    sumcpe = 0.d0
    sumpe    = 0.d0
    sumvir   = 0.d0
    sumtke   = 0.d0
    sumrke   = 0.d0
    sumte    = 0.d0
    sumprs   = 0.d0
    sumttp   = 0.d0
    sumrtp   = 0.d0
    sumtmp   = 0.d0   
    sumindpe = 0.d0
    sum3bpe = 0.d0
    sumtotke = 0.d0 
    sumvolm = 0.d0  

    ssqvdwpe = 0.d0
    ssqcpe = 0.d0
    ssqpe    = 0.d0
    ssqvir   = 0.d0
    ssqtke   = 0.d0
    ssqrke   = 0.d0
    ssqte    = 0.d0
    ssqprs   = 0.d0
    ssqttp   = 0.d0
    ssqrtp   = 0.d0
    ssqtmp   = 0.d0
    ssqindpe = 0.d0
    ssq3bpe = 0.d0
    ssqtotke = 0.d0
    ssqvolm = 0.d0

  end subroutine zero

  subroutine tot_sum

    implicit none

! accumulate sums for run averages

    sumvdwpe = sumvdwpe + stpvdwpe
    sumcpe = sumcpe + stpcpe
    sumpe = sumpe + stppe 
    sumvir =  sumvir + stpvir
    sumtke = sumtke + stptke
    sumrke = sumrke + stprke
    sumte = sumte + stpte
    sumprs = sumprs + stpprs
    sumttp = sumttp + stpttp
    sumrtp = sumrtp + stprtp
    sumtmp = sumtmp + stptmp
    sumindpe = sumindpe + stpindpe
    sum3bpe = sum3bpe + stp3bpe
    sumtotke = sumtotke + stptotke
    sumvolm = sumvolm + stpvolm 

    ssqvdwpe = ssqvdwpe + stpvdwpe**2
    ssqcpe =  ssqcpe + stpcpe**2
    ssqpe = ssqpe + stppe**2
    ssqvir = ssqvir + stpvir**2
    ssqtke = ssqtke + stptke**2
    ssqrke = ssqrke + stprke**2
    ssqte = ssqte + stpte**2
    ssqprs = ssqprs + stpprs**2
    ssqttp = ssqttp + stpttp**2
    ssqrtp = ssqrtp + stprtp**2
    ssqtmp = ssqtmp + stptmp**2
    ssqindpe = ssqindpe + stpindpe**2
    ssq3bpe = ssq3bpe + stp3bpe**2
    ssqtotke = ssqtotke + stptotke**2
    ssqvolm = ssqvolm + stpvolm**2
   
 
  end subroutine tot_sum

  subroutine averages_fluctuations
  
    implicit none  
    
    navge = nsteps - nequil - laststp

! calculate thermodynamic averages

    avvdwpe = sumvdwpe/navge
    avcpe = sumcpe/navge
    avpe = sumpe /navge
    avvir = sumvir/navge
    avtke = sumtke/navge
    avrke = sumrke/navge
    avte = sumte /navge
    avprs = sumprs/navge
    avttp = sumttp/navge
    avrtp = sumrtp/navge
    avtmp = sumtmp/navge
    avindpe = sumindpe/navge
    av3bpe = sum3bpe/navge
    avtotke = sumtotke/navge
    avvolm = sumvolm/navge

! calculate root mean square fluctuations    
 
    flcvdwpe = sqrt(abs(ssqvdwpe/navge - avvdwpe**2)) 
    flccpe = sqrt(abs(ssqcpe/navge - avcpe**2))
    flcpe = sqrt(abs(ssqpe/navge - avpe**2))
    flcvir = sqrt(abs(ssqvir/navge - avvir**2))
    flctke = sqrt(abs(ssqtke/navge - avtke**2))
    flcrke = sqrt(abs(ssqrke/navge - avrke**2))
    flcte = sqrt(abs(ssqte/navge - avte**2))
    flcprs = sqrt(abs(ssqprs/navge - avprs**2))
    flcttp = sqrt(abs(ssqttp/navge - avttp**2))
    flcrtp = sqrt(abs(ssqrtp/navge - avrtp**2))
    flctmp = sqrt(abs(ssqtmp/navge - avtmp**2)) 
    flcindpe = sqrt(abs(ssqindpe/navge - avindpe**2)) 
    flc3bpe = sqrt(abs(ssq3bpe/navge - av3bpe**2))
    flctotke = sqrt(abs(ssqtotke/navge - avtotke**2))
    flcvolm = sqrt(abs(ssqvolm/navge - avvolm**2))

  end subroutine averages_fluctuations
 

end module statistics

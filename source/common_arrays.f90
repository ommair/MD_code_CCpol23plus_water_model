module common_arrays
  USE mpi;
  implicit none;

!******************************************************************************
  ! unit of time      (to)    = 1.000000000 x 10**(-12) seconds   (pico-seconds)
  ! unit of length    (lo)    = 1.000000000 x 10**(-10) metres    (angstroms)
  ! unit of mass      (mo)    = 1.660540200 x 10**(-27) kilograms (daltons)
  ! unit of charge    (qo)    = 1.602177330 x 10**(-19) coulombs  (electrons)
  ! unit of energy    (eo)    = 1.660540200 x 10**(-23) joules    (10 j mol^-1)
  ! unit of pressure  (po)    = 1.660540200 x 10**(  7) pascals  (163.882576
  ! atm).
  !                           = 1.638825760 x 10**(  2) atmospheres
!*****************************************************************************

!************************************************************************
!*
!*          index of parameters and common variables
!*          ========================================
!*
  integer, parameter :: nom=5000     ! max number of molecules (32,108,256,500,864,1372,2048,2916)
  integer, parameter :: nsite=50   ! max number of sites 
  integer, parameter :: npimax=50  ! number of pair interactions
  integer, parameter :: ncsmax=3   ! max coreshell pairs 
  integer, parameter :: npolc=3  ! number of polarization sites per molecule
  integer, parameter :: nhist = 1000 !dimension of rdf histogram 
  integer, parameter :: ksqmax = 200 ! max squred vlues of k vectors indices in FOurier sum 
  real(8), parameter :: pi=dacos(-1.d0)  ! pi
  real(8), parameter :: boltz=8.31451115d-1  ! Boltzman Const eo/(k*mol) (dlpoly)
  real(8), parameter :: avsno=6.0222e+23    ! Avogadro number
  real(8), parameter :: r4pie0=138935.4835d0 ! 1/4pi*epsilon in internal units
  real(8), parameter :: ewp = 1.0d-6   ! precision for Ewald sum
  real(8), parameter :: prsunt=0.163882576d0  ! internal units to katm
  real(8), parameter :: bohr2a = 0.529177249d0
  real(8), parameter :: h2kcal = 627.510d0   ! Hartree to kcal/mol conversion
  real(8), parameter :: twosqpi= 2.0d0/sqrt(pi)  

  integer  :: nm                   ! number of molecule in the system
  integer  :: ns                   ! sites per molecules
  integer  :: natom                ! number  of massive sites per molecules  
  integer :: npi                   ! number of pair interactions in field file
  character(len=15) :: system_name ! name of the system in the field file
  character(len=3) :: an(nsite)    ! name of atoms in a molecule
  character(len=3) :: ann(nsite,nom)    ! name of atoms in a molecule

  real(8) :: alpha,alpha3  ! for ewald summation (intraction of permanent charges and point induced dipoles)
  !real(8), parameter :: alphad = 0.35818465241608965 !0.28d0  ! for induction model with GSF 0.28
  real(8) :: alphad
  integer :: kmax1,kmax2,kmax3
  character(len=15) :: initial_config
  real(8) :: hsize  ! size of step for force check
  logical :: check_forces ! (T or F)
  logical :: rand_vcom  ! random velocities
  logical :: rand_quat  ! random quaternios
  logical :: restart

!  RDF varaibles **** starts here
  integer ::  hist_freq ! frequncy of history file data
  integer :: ityp(nsite) ! types of sites (for RDF)
  integer :: numrdf ! number steps used fpr rdf statistics 
  real(8), parameter :: drdf = 0.05d0
  integer :: mxrdf 
  real(8) :: delrdf ! default bin width
! RDF varaiables **** ends here  

  real(8) :: diameter ! molecular diameter
  character(len=3) :: an1(npimax),an2(npimax) ! site pair names in vdw interactions  
  real(8) :: massites(nsite) ! mass of sites 
  real(8) :: coords(3,nsite) ! coordinates of given molecular geometry 
  real(8) :: pcoord(3,nsite) ! principal corrdiates of each sites in molecule
  real(8) :: pmi(3)          ! principal moment of inertia
  real(8) :: charge(nsite)   ! charge on sites
  real(8) :: sig(npimax),eps(npimax) ! lj parameters 
  real(8) :: ssig(nsite,nsite),eeps(nsite,nsite)  ! test
  real(8) :: temp0  ! desired temperature
  real(8) :: pres0  ! desired pressure
  real(8) :: pres   
  real(8) :: ttemp  ! temp related to translation KE
  real(8) :: rtemp  ! temp rlated to rotational KE
  real(8) :: deltmp ! tolerence in temp
  real(8) :: volm,volm0,volmt   ! volume of simulation box
  real(8) :: totm   ! mass of single molecule
  real(8) :: dt     ! time step size
  real(8) :: rcut,rcutsd ! real space cutoff radius / squared
  real(8) :: rcut3,rcutsd3
  real(8) :: box,box0,boxt,hbox    ! box length in each direction / half of box 
  integer :: nsteps,nequil ! total steps, no of equilibration steps
  integer :: step,stepr,navge,laststp ! step number , no of steps over which avg is taken

 !  arrays for molecular centers  
  real(8) :: x(nom),y(nom),z(nom) ! position of com of molecule
  real(8) :: q0(nom),q1(nom),q2(nom),q3(nom) ! quaternions   
  real(8) :: vx(nom),vy(nom),vz(nom)  !  velocity of com of molecule 
  real(8) :: jx(nom),jy(nom),jz(nom)  ! angular momentum of com of molecule
  real(8) :: fx(nom),fy(nom),fz(nom)  ! force on com of molecule
  real(8) :: tx(nom),ty(nom),tz(nom)  ! torque on com of molecule

  real(8) :: xs(nsite,nom), ys(nsite,nom), zs(nsite,nom) ! sites corrds
  real(8) :: vxs(nsite,nom), vys(nsite,nom), vzs(nsite,nom) ! sites velocites
  real(8) :: fxs(nsite,nom),fys(nsite,nom),fzs(nsite,nom) ! forces on sites
  real(8) :: chgs(nsite,nom)                              ! charges on sites

  real(8) :: xcom,ycom,zcom,vxcom,vycom,vzcom,sysmass

  integer :: integrator ! integer to select integrator
  character(len=15) :: pot_name   ! water potential
  character(len=10) :: ensemble   ! ensemble for simulation

  real(8) :: vir_test,vircom
  real(8) :: vdwpe,vir,tke,rke ! PE, virial, trans. KE, rot. KE
  real(8) :: rcpe,kcpe,scpe ! real, kspace, self coulumb pe
  real(8) :: vir_vdw,vir_rc,vir_kc,vir_self,vir_rind,vir_kind,vir_selfind,vir_shtind
  real(8) :: indpe          ! induction energy (point dipole polarization model)
  real(8) :: fs2pe,fs3pe          ! 3b energy  
  real(8) :: vdwlrc,virlrc  ! long corrections for VDW and virial
  real(8) :: sqf,sqt        
  real(8) :: stpvdwpe ! step vdw pe
  real(8) :: stpcpe   ! sum of real, kspace, self coulumb pe
  real(8) :: stpindpe ! step induction energy
  real(8) :: stp3bpe ! step 3b energy
  real(8) :: stppe          ! step total pe
  real(8) :: stptke ! step inst trans. ke
  real(8) :: stpttp ! step inst trans. temp
  real(8) :: stprke ! inst rot. ke
  real(8) :: stprtp ! inst rot. temp
  real(8) :: stpvir ! step vir
  real(8) :: stpte  ! step total energy
  real(8) :: stpprs ! step presurre 
  real(8) :: stptotke ! step tot kinetic energy
  real(8) :: stptmp ! step total temperature
  real(8) :: stpvolm ! step volume

  real(8) :: sumttp,sumrtp,sumprs,sumvolm
  real(8) :: sumpe,sumvdwpe,sumcpe,sumvir,sumtke,sumrke,sumte,sumtmp
  real(8) :: sumindpe, sum3bpe,sumtotke

  real(8) :: ssqvdwpe,ssqcpe,ssqpe,ssqvir,ssqtke,ssqrke,ssqte,ssqprs
  real(8) :: ssqttp,ssqrtp,ssqtmp,ssqindpe,ssq3bpe,ssqtotke,ssqvolm

  real(8) :: avvdwpe,avcpe,avindpe,av3bpe,avtotke,avvolm
  real(8) :: avpe,avvir,avtke,avrke,avte,avprs,avttp,avrtp,avtmp  
  real(8) :: flcvdwpe,flccpe,flcindpe,flc3bpe,flctotke,flcvolm 
  real(8) :: flcpe,flcvir,flctke,flcrke,flcte,flcprs,flcttp,flcrtp,flctmp

  real(8) :: t, tr, taut, taup ! step time

  ! sapt5s,  ccpol5 or ccpol8s parameters variables

  real(8) :: a0(npimax),a1(npimax),a2(npimax),a3(npimax)
  real(8) :: expalp(npimax), beta(npimax),aalp(npimax)
  real(8) :: d1(npimax), d6(npimax),d8(npimax),d10(npimax)
  real(8) :: c0(npimax),c1(npimax),c2(npimax),c3(npimax)
  real(8) :: c6(npimax),c8(npimax),c10(npimax)
  real(8) :: domo(npimax),qa(npimax),qb(npimax)

  ! induction model parameter

  logical :: induction      ! true or false (wants to include induction model) 
  character(len=10) :: induction_type         ! (damped or undamped)
  character(len=18) :: ind_model              ! physical dipole or point dipole  
  character(len=10) :: method       ! method for long range correctin to avoid total energy drift 

  real(8) :: E0x(nsite,nom),E0y(nsite,nom),E0z(nsite,nom)    ! efield of static charges
  real(8) :: Eidmx(nsite,nom),Eidmy(nsite,nom),Eidmz(nsite,nom) ! efld of ind dipoles

  real(8) :: rE0x(nsite,nom),rE0y(nsite,nom),rE0z(nsite,nom)    ! efield of static charges (real space)
  real(8) :: kE0x(nsite,nom),kE0y(nsite,nom),kE0z(nsite,nom)    ! efield of static charges (recp space)
  real(8) :: slfE0x(nsite,nom),slfE0y(nsite,nom),slfE0z(nsite,nom)    ! efield of static charges (self correction)
  real(8) :: srfE0x(nsite,nom),srfE0y(nsite,nom),srfE0z(nsite,nom) 
  real(8) :: shtE0x(nsite,nom),shtE0y(nsite,nom),shtE0z(nsite,nom)

  real(8) :: rEidmx(nsite,nom),rEidmy(nsite,nom),rEidmz(nsite,nom) ! efld of ind dipoles (real space)
  real(8) :: kEidmx(nsite,nom),kEidmy(nsite,nom),kEidmz(nsite,nom) ! efld of ind dipoles (recp space)
  real(8) :: slfEidmx(nsite,nom),slfEidmy(nsite,nom),slfEidmz(nsite,nom) ! efld of ind dipoles (self correction)
  real(8) :: srfEidmx(nsite,nom),srfEidmy(nsite,nom),srfEidmz(nsite,nom)
  real(8) :: shtEidmx(nsite,nom),shtEidmy(nsite,nom),shtEidmz(nsite,nom)
  real(8) :: Etotx(nsite,nom),Etoty(nsite,nom),Etotz(nsite,nom)

  real(8) :: idmx(nsite,nom),idmy(nsite,nom),idmz(nsite,nom) ! induced dipoles

  real(8) :: polarizability(nsite)
  real(8) :: apol(nsite,nom)
  integer :: npol

    real(8), parameter :: delta2 = 1.65d0/bohr2a  ! 1/Bohr to 1/Ang
    real(8), parameter :: delta3 = 1.55d0/bohr2a

!  real(8), parameter :: delta2 = 1.756d0/bohr2a  ! 1/Bohr to 1/Ang
!  real(8), parameter :: delta3 = 1.65d0/bohr2a

  ! core-shell induction model params
  integer  :: ncs                  ! number of core shell pairs
  integer :: icore(ncsmax),ishell(ncsmax)    ! site number for core and shell
  real(8) :: spring_k(ncsmax)         ! spring constant k in internal units
  real(8) :: spge ! spring energy (core-shell)

  ! 3B interaction arrays
  character(len=3) :: a3b(npimax),b3b(npimax) ! site pair names in 3b interactions
  logical :: nonadditive_3B      ! true or false (wants to include 3b potential)
  integer :: ipar
  real(8) :: params(10)
  real(8) :: eta3, beta3, g3 ! FS^2 parameters
  real(8) :: bS3(nsite,nsite),gS3(nsite,nsite),r0S3(nsite,nsite)
  character(len=3) ::  buf1(nsite), buf2(nsite)
  real(8) :: cs3(5000)
  real(8) :: ind_beta(nsite,nsite),ind_gamma(nsite,nsite),ind_r0(nsite,nsite) 
  real(8) :: fxfs2(nom),fyfs2(nom),fzfs2(nom)  
  real(8) :: txfs2(nom),tyfs2(nom),tzfs2(nom)  

  real(8) :: fxfs3(nom),fyfs3(nom),fzfs3(nom)
  real(8) :: txfs3(nom),tyfs3(nom),tzfs3(nom)
  
  ! OMP
  include 'omp_lib.h'
  integer(4) :: OMP_NumOfThreads = 2;
  integer(4) :: id;
  ! OMP
  ! TMP
  real(8), allocatable :: rndArr(:);
  logical :: rndInit = .FALSE.;
  integer(4) :: rndIdx = 1;
  ! TMP
  ! MPI
  INTEGER(4) :: MPIstat(MPI_STATUS_SIZE);
  INTEGER(4) :: MPIierr = 0, MPIcomm_size = 1, MPIrank = 0, MPIchunk_size = 1, MPIaddr = 0, MPIiArr(0:1);
  INTEGER(4), PARAMETER :: MPI_TAG_CommReady = 10, MPI_TAG_SendDataIdx = 11, MPI_TAG_SendData = 12;
  INTEGER(4), ALLOCATABLE :: MPItskIdx(:);
  ! DBG
  INTEGER(4) :: dbgUnit = 0;
  ! DBG
  ! MPI

contains

subroutine InitRandom
	integer(4) i, N, ierr, seed;
	open(101, file = 'Random.num', STATUS = 'OLD', IOSTAT = ierr);
	if (ierr .EQ. 0) then
		read(101, *) N;
		allocate(rndArr(1:N)); 
		do i = 1, N
			read(101, *) rndArr(i);
		enddo
	else
		open(101, file = 'Random.num');
		N = 10000;
		allocate(rndArr(1:N)); 
		seed = 96873683;
		call random_seed(seed);
		write(101, *) N;
		do i = 1, N
			call random_number(rndArr(i));
			write(101, '(1pE26.17E3)') rndArr(i);
		enddo
	endif
	close(101);
	rndInit = .TRUE.;
	return;
end subroutine InitRandom

subroutine random_number0(res)
	real(8) res;

	if (.NOT. rndInit) call InitRandom;
	if (rndIdx .GT. SIZE(rndArr)) rndIdx = 1;
	res = rndArr(rndIdx);
	rndIdx = rndIdx + 1;
	return;
end subroutine random_number0

subroutine random_number1(res)
	real(8) res(:);
	integer(4) i;

	if (.NOT. rndInit) call InitRandom;
	do i = 1, size(res)
		if (rndIdx .GT. SIZE(rndArr)) rndIdx = 1;
		res(i) = rndArr(rndIdx);
		rndIdx = rndIdx + 1;
	enddo
	return;
end subroutine random_number1

subroutine ShuffleArrayI(N, arr)
	integer(4) N, arr(0:N - 1);
	real(8), allocatable :: arrRand(:);

	allocate(arrRand(0:N - 1));
	call random_number(arrRand);
	! DBG
!	write(1000, '(I6)') arr;
!	write(1001, '(1pE16.6E3)') arrRand;
	! DBG
	call QuickRealL(N, arrRand, arr);
	! DBG
!	write(1002, '(I6)') arr;
!	write(1003, '(1pE16.6E3)') arrRand;
!	STOP;
	! DBG
	deallocate(arrRand);
	return;
end subroutine ShuffleArrayI

SUBROUTINE QuickRealL(N, D, ib, INFO, ID)
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
	INTEGER(4) N;
	INTEGER(4), optional :: INFO;
	CHARACTER, optional :: ID;
!     ..
!     .. Array Arguments ..
	REAL(8) D(*)
	INTEGER(4) ib(*);
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
	INTEGER(4), PARAMETER :: SELECT = 20;
!     ..
!     .. Local Scalars ..
	INTEGER(4) DIR, ENDD, I, J, START, STKPNT;
	integer(4) tmpI, nfo;
	real(8) D1, D2, D3, DMNMX, TMP, tmpR1;
!     ..
!     .. Local Arrays ..
	INTEGER(4) STACK(2, 32);
!     ..
!     .. External Functions ..
!      LOGICAL            LSAME
!      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
!      EXTERNAL           XERBLA
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
	nfo = 0;
	DIR = -1;
	if (present(ID)) then
		IF (ID .EQ. 'D') THEN
			DIR = 0;
		ELSEIF (ID .EQ. 'I') THEN
			DIR = 1;
		ENDIF
	else
		DIR = 1;
	endif
	IF (DIR .EQ. -1) THEN
		nfo = -1;
	ELSEIF (N .LT. 0) THEN
		nfo = -2;
	ENDIF
	if (present(INFO)) INFO = nfo;
	IF (nfo .NE. 0) THEN
!		CALL THSTOP('DLASRT Error');
		STOP 'DLASRT Error';
!		RETURN;
	ENDIF
!
!     Quick return if possible
!
      IF (N .LE. 1) RETURN;

      STKPNT = 1
      STACK( 1, 1 ) = 1
      STACK( 2, 1 ) = N
   10 CONTINUE
      START = STACK( 1, STKPNT )
      ENDD = STACK( 2, STKPNT )
      STKPNT = STKPNT - 1
      IF (((ENDD - START) .LE. SELECT) .AND. ((ENDD - START) .GT. 0)) THEN
!
!        Do Insertion sort on D( START:ENDD )
!
         IF (DIR .EQ. 0) THEN
!
!           Sort into decreasing order
!
            DO 30 I = START + 1, ENDD
               DO 20 J = I, START + 1, -1
                  IF( D( J ).GT.D( J-1 ) ) THEN
                     DMNMX = D( J )
                     D( J ) = D( J-1 )
                     D( J-1 ) = DMNMX
                     tmpI = ib(j);
                     ib(j) = ib(j-1);
                     ib(j-1) = tmpI;
                  ELSE
                     GO TO 30
                  END IF
   20          CONTINUE
   30       CONTINUE

         ELSE
!
!           Sort into increasing order
!
            DO 50 I = START + 1, ENDD
               DO 40 J = I, START + 1, -1
                  IF( D( J ).LT.D( J-1 ) ) THEN
                     DMNMX = D( J )
                     D( J ) = D( J-1 )
                     D( J-1 ) = DMNMX
                     tmpI = ib(j);
                     ib(j) = ib(j-1);
                     ib(j-1) = tmpI;
                  ELSE
                     GO TO 50
                  END IF
   40          CONTINUE
   50       CONTINUE

         END IF

      ELSEIF ((ENDD - START) .GT. SELECT) THEN
!
!        Partition D( START:ENDD ) and stack parts, largest one first
!
!        Choose partition entry as median of 3
!
         D1 = D( START )
         D2 = D( ENDD )
         I = ( START+ENDD ) / 2
         D3 = D( I )
         IF( D1.LT.D2 ) THEN
            IF( D3.LT.D1 ) THEN
               DMNMX = D1
            ELSE IF( D3.LT.D2 ) THEN
               DMNMX = D3
            ELSE
               DMNMX = D2
            END IF
         ELSE
            IF( D3.LT.D2 ) THEN
               DMNMX = D2
            ELSE IF( D3.LT.D1 ) THEN
               DMNMX = D3
            ELSE
               DMNMX = D1
            END IF
         END IF

         IF (DIR .EQ. 0) THEN
!
!           Sort into decreasing order
!
            I = START - 1
            J = ENDD + 1
   60       CONTINUE
   70       CONTINUE
            J = J - 1
            IF( D( J ).LT.DMNMX ) GO TO 70
   80       CONTINUE
            I = I + 1
            IF( D( I ).GT.DMNMX ) GO TO 80
            IF( I.LT.J ) THEN
               TMP = D( I )
               D( I ) = D( J )
               D( J ) = TMP
               tmpI = ib(i);
               ib(i) = ib(j);
               ib(j) = tmpI;
               GO TO 60
            END IF
            IF( J-START.GT.ENDD-J-1 ) THEN
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = START
               STACK( 2, STKPNT ) = J
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = J + 1
               STACK( 2, STKPNT ) = ENDD
            ELSE
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = J + 1
               STACK( 2, STKPNT ) = ENDD
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = START
               STACK( 2, STKPNT ) = J
            END IF
         ELSE
!
!           Sort into increasing order
!
            I = START - 1
            J = ENDD + 1
   90       CONTINUE
  100       CONTINUE
            J = J - 1
            IF( D( J ).GT.DMNMX ) GO TO 100
  110       CONTINUE
            I = I + 1
            IF( D( I ).LT.DMNMX ) GO TO 110
            IF( I.LT.J ) THEN
               TMP = D( I )
               D( I ) = D( J )
               D( J ) = TMP
               tmpI = ib(i);
               ib(i) = ib(j);
               ib(j) = tmpI;
               GO TO 90
            END IF
            IF( J-START.GT.ENDD-J-1 ) THEN
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = START
               STACK( 2, STKPNT ) = J
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = J + 1
               STACK( 2, STKPNT ) = ENDD
            ELSE
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = J + 1
               STACK( 2, STKPNT ) = ENDD
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = START
               STACK( 2, STKPNT ) = J
            END IF
         END IF
      END IF
      IF( STKPNT.GT.0 ) GO TO 10
      RETURN
!
!     End of DLASRT
!
!END SUBROUTINE DLASRT
END SUBROUTINE QuickRealL

end module common_arrays

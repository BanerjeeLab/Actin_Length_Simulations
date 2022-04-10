	program code

	implicit none

	real*8, parameter :: time=600.d0, t0=0			!time-60
	real*8, parameter :: dt=1.d0,dtmax=0.001d0
	real*8, parameter :: sdur=2.d2
	real*8, parameter :: er=1.d-6, dl=1.d0/370.d0
	real*8, parameter :: vton=600.d0

	real*8,parameter :: rho0 = 10.d0			!! actin density in micro-Molar
	real*8,parameter :: rhop = 0.5d0				!! profilin density in micro-Molar
!	real*8,parameter :: rhoc = 5.d0				!! capping density in nano-Molar
!	real*8,parameter :: rhof = 5.d0				!! formin density in nano-Molar
	real*8,parameter :: V = 200.d0				!! volume in micro-m^3

	integer, parameter :: N=nint(rho0*V*vton)		!! actin monomer pool

	real*8,parameter :: kp0 = 11.6d0/vton, km0 = 1.4d0		!! polymerization in barbed end
	real*8,parameter :: k1p = 35.7d0/vton, k1m = 1.63d8		!! dimer formation : reaction-b@Sept
	real*8,parameter :: k2p = 2.18d0/vton, k2m = 1.3d3		!! trimer formation : reaction-f@Sept
	real*8,parameter :: kpp = 40.0d0/vton, kpm = 0.75d0		!! profilin-monomer reaction : reaction4@Bieling
	real*8,parameter :: kpp0 = 11.d0/vton, kpm0 = 1450.d0		!! profilin-actin -filament reaction : reaction2@Bieling
	real*8,parameter :: kpam0 = 50.d0
	real*8,parameter :: kfp0=29.1d0/vton, kfm0=8.1d-5			!! formin filament reactions : reaction-Table1@Guillaume
	real*8,parameter :: kcp0=12.8d0/vton, kcm0=2.d-4			!! capping reactions : reaction-Table1@Guillaume
	real*8,parameter :: kfcp0=0.21d0/vton, kfcm0=1.8d-3, kfcm20=4.2d-3	!! cap-formin interaction : reaction-Table1@Guillaume
	real*8,parameter :: kcfp0=1.6d0/vton, kcfm0=6.2d-3, kcfm20=6.2d-3	!! cap-formin interaction : reaction-Table1@Guillaume
	real*8,parameter :: factor1=0.2d0					!! formin: barbed end growth by only actin (mdia1 - 1) (mdia2 - 0.2)
	real*8,parameter :: factor2=1.5d0					!! formin: barbed end growth by profilin-actin (mdia1 - 6) (mdia2 - 1.5)
	real*8,parameter :: factor3=3.d0
        
	real*8,parameter :: kfnuc=1.d-5/(vton*vton)					!! formin mediated nucleation from dimer : Courtemanche (2.d-4/(vton*vton))
        
	integer,parameter :: nf_max=5000			!nint(1.d0*N/3.d0)

!	integer,parameter :: Nfr=nint(rhof*V*vton/1000.d0)				!! formin pool
!	integer, parameter :: Nc=nint(rhoc*V*vton/1000.d0)	!! capping pool
	integer, parameter :: Np0=nint(rhop*V*vton)		!! profilin pool
	integer, parameter :: nstep=nint(time)

	real*8, parameter :: delt=time/nstep

!--------------- input parameters ------------------

	integer :: Nfr, Nc

!---------------------------------------------------

	real*8 :: t,r1,r2,pl,ph
	real*8 :: tau
	real*8 :: nm, nd, np, npa
	real*8 :: beta_a(1:nf_max),beta_d(1:nf_max)
	real*8 :: beta_cp(1:nf_max),beta_cm(1:nf_max)
	real*8 :: beta_fp(1:nf_max),beta_fm(1:nf_max)
	real*8 :: beta_fcp(1:nf_max),beta_fcm(1:nf_max),beta_fcm2(1:nf_max)
	real*8 :: beta_cfp(1:nf_max),beta_cfm(1:nf_max),beta_cfm2(1:nf_max), fc_sum
	real*8 :: beta_pp(1:nf_max),beta_pm(1:nf_max),beta_pam(1:nf_max)
	real*8 :: beta_fma, beta_fmd, beta_tma, beta_tmd
	real*8 :: beta_na, beta_nd
	real*8 :: beta_npa, beta_npd
	real*8 :: beta_fnuc1, beta_fnuc2, beta_fp3, beta_fp4, beta_fp5, beta_fm3, beta_fm4, beta_fm5

	real*8 :: kp(1:nf_max),km(1:nf_max), kcp(1:nf_max),kcm(1:nf_max)
	real*8 :: kfp(1:nf_max),kfm(1:nf_max)
	real*8 :: kfcp(1:nf_max),kfcm(1:nf_max),kfcm2(1:nf_max),kcfp(1:nf_max),kcfm(1:nf_max),kcfm2(1:nf_max)
	real*8 :: kfpp(1:nf_max),kfpm(1:nf_max),kfpam(1:nf_max)
	real*8 :: alpha0, beta_sum0, beta_sum1, beta_sum
	real*8 :: dum1, dum2, dum_sum, dumran
	real*8 :: pdum1,pdum2

	integer :: m, dm, tm, fm			! m=A1, dm=A2, tm=A3, fm=A4
	integer :: p, pa				! p=P, pa=PA
	integer :: df, tf, ff				! df=A2F, tf=A3F, ff=A4F
	integer :: cap, frm, cpfr
	integer :: l5, lp5
	integer :: nf
	integer :: nnt, Rid

	integer :: i,j,k,nt,nk
	integer :: l(1:nf_max)
	integer :: ftrimer, filcount
	integer :: Bstate(1:nf_max)
	integer :: dbug

	integer :: kk, cnt, scount, rcount
	INTEGER :: seed(12)

	character(len=70) :: fname

!------------------------------------------------------

	open(unit=100,file='data/m.txt',status='unknown')
	open(unit=101,file='data/ml.txt',status='unknown')
	open(unit=200,file='data/l.txt',status='unknown')
	open(unit=201,file='data/pl.txt',status='unknown')
	open(unit=202,file='data/state.txt',status='unknown')
	open(unit=203,file='data/fil.txt',status='unknown')
	open(unit=500,file='ran.txt',status='unknown')
	open(unit=900,file='xx.txt',status='unknown')

!--------------------- read parameter -----------------

	read(*,*) dum1,dum2		!! concentration in nano-molar
	Nc=nint(dum1*vton*V*1.d-3)
	Nfr=nint(dum2*vton*V*1.d-3)

!---------------------- seed random -------------------

  	read(500,*) dumran

  	seed= 987877+nint(dumran*1091)

        CALL RANDOM_SEED (PUT=seed)

	write(*,*) Nc, Nfr
!  	----------Barbed end state:-----------
!	Bstate(i) = 1 FA 
!	Bstate(i) = 2 FA+FR
!	Bstate(i) = 3 FA+CP 
!	Bstate(i) = 4 FA+CP-FR or FA+FR-CP 
!	Bstate(i) = 5 FA+FR+PA // Bstate(i) = 6 FA+PA

	t=0
	nk=0
	scount = 1		! count random number gen seed
	rcount = 0
	Bstate = 0		! flag for filament barbed end state : 0 = not nucleated
	filcount = 0		! count number of nucleation i.e., when A1+A3 -> A4
	nf = 0			! count active(growing>4) filaments
	dbug=98
!----------- reaction rates -----------

!------- no f-actin present at t=0 -----------
	kp=0.d0
	km=0.d0
!------- hence, no capping reaction ----------
	kcp=0.d0
	kcm=0.d0
!------- and no formin reaction ----------
	kfp=0.d0
	kfm=0.d0
!------- no cap formin interaction -------
	kfcp=0.d0
	kfcm=0.d0
	kfcm2=0.d0
	kcfp=0.d0
	kcfm=0.d0
	kcfm2=0.d0
!------- and no profilin-filament reaction -------
	kfpp=0.d0
	kfpm=0.d0
	kfpam=0.d0


!----------- initialisation -----------

	l=0
	dm=20
	tm=2
	fm=0
	df=0
	tf=0
	ff=0
	cap=Nc
	pa=0
	p=Np0-pa
	frm=Nfr-df-tf-ff
	cpfr=0
	m=N-sum(l)-2*dm-3*tm-4*fm -2*df-3*tf-4*ff -pa

!	write(200,*) t,(l(i), i=1,40),sum(l)+m+2*dm+3*tm+4*fm

!***************************************************************
!------------------- time loop ------------------------
!***************************************************************

	do while (t.le.time)
	ftrimer = 0		! flag for trimer formation, 1 for new trimer
	l5 = 0			! count of filament of length 5
	lp5 = 0			! count of filament of length 5 + PA

!------------------ propensity sum ---------------------

	if (filcount.gt.0) then
	do i=1,filcount

	if (l(i).eq.5) then
	if (Bstate(i).eq.6) then
	lp5=lp5+1
	else
	l5=l5+1
	endif
	endif

	enddo
	endif

	beta_tma = k2p*m*dm/(V)
	beta_tmd = k2m*tm

	beta_fma = kp0*m*tm/(V)
	beta_fmd = km0*fm

	beta_na = kp0*m*fm/(V)
	beta_nd = km0*l5
	beta_npa = kpp0*pa*fm/(V)
	beta_npd = kpam0*lp5

	beta_fnuc1 = kfnuc*frm*(m/V)*(m/V)
	beta_fnuc2 = kfnuc*frm*(pa/V)*(m/V)

	beta_fp3 = factor2*kp0*df*pa/V
	beta_fp4 = factor2*kp0*tf*pa/V
	beta_fp5 = factor2*kp0*ff*pa/V

	do i=1,nf_max
	beta_a(i) = kp(i)*m/V
	if (l(i).gt.5) then
	beta_d(i) = km(i)
	else
	beta_d(i) = 0.d0
	endif

	beta_cp(i) = kcp(i)*cap/V
	beta_cm(i) = kcm(i)
	beta_fp(i) = kfp(i)*frm/V
	beta_fm(i) = kfm(i)

	beta_fcp(i) = kfcp(i)*cap/V
   	beta_fcm(i) = kfcm(i)
	beta_fcm2(i) = kfcm2(i)
	beta_cfp(i) = kcfp(i)*frm/V
	beta_cfm(i) = kcfm(i)
	beta_cfm2(i) = kcfm2(i)

	beta_pp(i) = kfpp(i)*pa/V
	beta_pm(i) = kfpm(i)
	beta_pam(i) = kfpam(i)
	enddo

	fc_sum = sum(beta_fcp) + sum(beta_fcm) + sum(beta_fcm2) + sum(beta_cfp) + sum(beta_cfm) + sum(beta_cfm2)
	beta_sum0 = sum(beta_a) + sum(beta_d) + sum(beta_cp) + sum(beta_cm) + sum(beta_fp) + &
                    sum(beta_fm) + sum(beta_pp) + sum(beta_pm) + sum(beta_pam)
	beta_sum1 = fc_sum + beta_sum0 + beta_fma + beta_fmd + beta_tma + beta_tmd + beta_na + beta_npa + beta_nd + beta_npd
	alpha0 = beta_sum1 + beta_fnuc1 + beta_fnuc2 + beta_fp3 + beta_fp4 + beta_fp5

	if(alpha0.lt.0.d0) then
	write(*,*) '*** ERROR *** PROPENSITY ZERO', alpha0, Rid
	write(*,*) 'propensities:', fc_sum, beta_sum0, beta_sum1
	GOTO 100
	endif

!----------------- next reaction time: tau --------------

	CALL RANDOM_NUMBER(r1)
	if(r1.lt.1.d-8)then		!and hope the next one is not zero too :)
	CALL RANDOM_NUMBER(r1)
	endif

	tau = (1.d0/alpha0)*log(1.d0/r1)

	if(tau.lt.0.d0) then
	write(*,*) '*** ERROR *** TIME STEP ZERO', tau, Rid
	GOTO 100
	endif

!----------------------- Progress dimerization as ODEs for time t -> t+tau ------------------------

	nm=m
	nd=dm
	npa=pa
	np=p
	do k=1,1000
	nd = nd + ( (k1p/k1m)*(nm*(nm-1))/V - nd )*dt
	npa = npa + ( (kpp/k1m)*(nm*np)/V - (kpm/k1m)*npa )*dt
	np = np - ( (kpp/k1m)*(nm*np)/V - (kpm/k1m)*npa )*dt
	nm = nm - 2.d0*( (k1p/k1m)*(nm*(nm-1))/V - nd )*dt - ( (kpp/k1m)*(nm*np)/V - (kpm/k1m)*npa )*dt
	enddo
	dm=nint(nd)
	p=nint(np)
	pa=nint(npa)
	m=N-sum(l)-2*dm-3*tm-4*fm-2*df-3*tf-4*ff-pa

	if(m.lt.er) then
	write(*,*) '*** ERROR *** MONOMER ZERO', m
	GOTO 100
	endif
	
!---------------------------------------------------------------------

!----------------------- Progress with stochastic reactions ---------

	CALL RANDOM_NUMBER(r2)

	beta_sum = 0.d0

!---------------- Trimerization dm+m -> tm ------------------	

	pl=beta_sum/alpha0
	ph=(beta_sum+beta_tma)/alpha0
	
	if (r2.ge.pl.and.r2.lt.ph) then	
	tm = tm+1			
	dm = dm-1
	m = m-1
	Rid=1
	!write(*,*)t,'A3+1', pl, r2, ph
	endif

	pl=ph
	ph=(beta_sum+beta_tma+beta_tmd)/alpha0
	
	if (r2.ge.pl.and.r2.lt.ph) then	
	tm = tm-1			
	dm = dm+1
	m = m+1
	Rid=2
	endif

!---------------- Formin mediated dimerization frm + m + m -> df ------------------	

	pl=ph
	ph=(beta_sum+beta_tma+beta_tmd+beta_fnuc1)/alpha0
	
	if (r2.ge.pl.and.r2.lt.ph) then	
	m = m-2	
	frm = frm-1		
	df = df+1
	Rid=3
	endif

!---------------- Formin mediated dimerization frm + m + pa -> df ------------------	

	pl=ph
	ph=(beta_sum+beta_tma+beta_tmd+beta_fnuc1+beta_fnuc2)/alpha0
	
	if (r2.ge.pl.and.r2.lt.ph) then	
	m = m-1	
	pa = pa-1
	frm = frm-1		
	df = df+1
	p = p+1
	Rid=4
	endif

!---------------- Formin reaction 3 df + pa -> tf ------------------	

	pl=ph
	ph=(beta_sum+beta_tma+beta_tmd+beta_fnuc1+beta_fnuc2+beta_fp3)/alpha0
	
	if (r2.ge.pl.and.r2.lt.ph) then	
	pa = pa-1		
	df = df-1
	tf = tf+1
	p = p+1
	Rid=5
	endif

!---------------- Formin reaction 4 tf + pa -> ff ------------------	

	pl=ph
	ph=(beta_sum+beta_tma+beta_tmd+beta_fnuc1+beta_fnuc2+beta_fp3+beta_fp4)/alpha0
	
	if (r2.ge.pl.and.r2.lt.ph) then	
	pa = pa-1		
	tf = tf-1
	ff = ff+1
	p = p+1
	Rid=6
	endif

!---------------- Formin reaction 5 ff + pa -> new filament ------------------	

	pl=ph
	ph=(beta_sum+beta_tma+beta_tmd+beta_fnuc1+beta_fnuc2+beta_fp3+beta_fp4+beta_fp5)/alpha0

	!********* add a new filament ********
	if (r2.ge.pl.and.r2.lt.ph) then	
	filcount = filcount + 1
	kp(filcount)=factor1*kp0
	kfpp(filcount)=factor2*kp0
	km(filcount)=km0
	kfm(filcount)=kfm0
	kfcp(filcount)=kfcp0
	l(filcount)=5			
	ff = ff-1
	pa = pa-1
	p = p+1
	nf = nf + 1
	Bstate(filcount) = 2
	Rid=7
	write(900,*) t, 'formin-bound filament created- FA4+PA -> FA-FR + P','index=',filcount
	endif

	beta_sum = beta_sum+beta_tma+beta_tmd+beta_fnuc1+beta_fnuc2+beta_fp3+beta_fp4+beta_fp5


!------------------------ Nucleation --------------------------

!---------------- tetramer formation (A1+A3 -> A4) ---------------
	pl=beta_sum/alpha0
	ph=(beta_sum+beta_fma)/alpha0
	
	if (r2.ge.pl.and.r2.lt.ph) then	
	fm = fm+1			
	tm = tm-1
	m = m-1
	Rid=8
	endif

	pl=ph
	ph=(beta_sum+beta_fma+beta_fmd)/alpha0
	
	if (r2.ge.pl.and.r2.lt.ph) then	
	fm = fm-1			
	tm = tm+1
	m = m+1
	Rid=9
	endif

!----------------- filament formation from A1: A4+A1 -> FA ---------------
	pl=ph
	ph=(beta_sum+beta_fma+beta_fmd+beta_na)/alpha0
	
	!********* add a new filament ********
	if (r2.ge.pl.and.r2.lt.ph) then	
	filcount = filcount + 1
	kp(filcount)=kp0
	km(filcount)=km0
	kcp(filcount)=kcp0
	kfp(filcount)=kfp0
	kfpp(filcount)=kpp0
	l(filcount)=5			
	fm = fm-1
	m = m-1
	nf = nf + 1
	Bstate(filcount) = 1
	Rid=10
	write(900,*) t, 'new filament created- A4+A1 -> FA','index=',filcount
	endif



!----------------- filament formation from PA: A4+PA -> FA-P ---------------
	pl=ph
	ph=(beta_sum+beta_fma+beta_fmd+beta_na+beta_npa)/alpha0
	
	!********* add a new filament ********
	if (r2.ge.pl.and.r2.lt.ph) then	
	filcount = filcount + 1
	kfpm(filcount) = kpm0
	kfpam(filcount) = kpam0
	l(filcount)=5			
	fm = fm-1
	pa = pa-1
	nf = nf + 1
	Bstate(filcount) = 6
	Rid=11
	write(900,*) t, 'new filament created- A4+PA -> FA-P','index=',filcount
	endif



	!********* destroy FA filament ********
	pl=ph
	ph=(beta_sum+beta_fma+beta_fmd+beta_na+beta_npa+beta_nd)/alpha0
	
	if (r2.ge.pl.and.r2.lt.ph) then	
	!------- find the first filament of l==5 -------
	do i=1,filcount
	if(l(i).eq.5.and.Bstate(i).le.2)then	!! select the first FA or FA-FR of length 5 to break

	kp(i)=0.d0
	km(i)=0.d0
	kcp(i)=0.d0
	kcm(i)=0.d0
	kfp(i)=0.d0
	kfm(i)=0.d0
	kfpp(i)=0.d0
	kfpm(i)=0.d0
	kfpam(i)=0.d0
	l(i)=0
	fm = fm+1
	nf = nf-1
	Rid=12
	write(900,*) 'fil destroy','index=',i,'Bstate=',Bstate(i)
	
	if (Bstate(i).eq.1) then		!! destroy a filament with free-barbed end
	m = m+1
	Bstate(i) = 0
	write(900,*) t, 'A5 -> A4 + A1','index=',i
	exit
	endif

	if (Bstate(i).eq.2) then		!! destroy a filament with formin-bound-barbed end
	m = m+1
	frm=frm+1
	Bstate(i) = 0
	write(900,*) t, 'A5F-> A4 + A1 + F','index=',i
	exit
	endif

	endif 	!! l(i)=5 if loop
	enddo	!! l(i) do loop
	
	endif


	!********* destroy FA-P filament ********
	pl=ph
	ph=(beta_sum+beta_fma+beta_fmd+beta_na+beta_npa+beta_nd+beta_npd)/alpha0
	
	if (r2.ge.pl.and.r2.lt.ph) then	
	!------- find the first filament of l==5 -------
	do i=1,filcount
	if(l(i).eq.5.and.Bstate(i).eq.6)then	!! select the first FA-P of length 5 to break

	kp(i)=0.d0
	km(i)=0.d0
	kcp(i)=0.d0
	kcm(i)=0.d0
	kfp(i)=0.d0
	kfm(i)=0.d0
	kfpp(i)=0.d0
	kfpm(i)=0.d0
	kfpam(i)=0.d0
	l(i)=0
	fm = fm+1
	nf = nf -1
	pa=pa+1
	Bstate(i) = 0
	Rid=13
	write(900,*) t, 'A5-P -> A4 + PA','index=',i
	exit

	endif 	!! l(i)=5 if loop
	enddo	!! l(i) do loop
	
	endif

!---------------- filament growth and decay ------------------

	beta_sum = beta_sum+beta_fma+beta_fmd+beta_na+beta_npa+beta_nd+beta_npd

	if (filcount.ge.1) then
	do i=1,filcount		!!---- nf-loop ----!!

!---------------- G-actin unbinding (B->B-1) ------------------

	pl=beta_sum/alpha0
	ph=(beta_sum+beta_d(i))/alpha0
	
	if (r2.ge.pl.and.r2.lt.ph) then				
	l(i) = l(i)-1
	m = m+1
	rcount = rcount + 1
	Rid=14
	endif


!---------------- G-actin binding (B->B+1) ------------------

	pl=ph
	ph=(beta_sum+beta_a(i)+beta_d(i))/alpha0

	if (r2.ge.pl.and.r2.lt.ph) then
	l(i) = l(i)+1
	m = m-1
	rcount = rcount + 1
	Rid=15
	endif

!---------------- profilin-actin binding (B->B+1) ------------------

	pl=ph
	ph=(beta_sum+beta_a(i)+beta_d(i)+beta_pp(i))/alpha0

	if (r2.ge.pl.and.r2.lt.ph) then
	l(i) = l(i)+1
	pa = pa-1
	if (Bstate(i).eq.1) then
	kfpm(i) = kpm0
	kfpam(i) = kpam0
	kfpp(i)=0.d0
	kp(i)=0.d0
	km(i)=0.d0
	kcp(i)=0.d0
	kfp(i)=0.d0
	Bstate(i) = 6
	elseif (Bstate(i).eq.2) then
	kfpm(i) = factor3*kpm0
	kfpam(i) = kpam0
	kfpp(i)=0.d0
	kp(i)=0.d0
	km(i)=0.d0
	kcp(i)=0.d0
	kfm(i)=0.d0
	kfcp(i)=0.d0
	Bstate(i) = 5
	endif
	rcount = rcount + 1
	Rid=16
	endif

!---------------- profilin unbinding (B,P->B+P) ------------------ **

	pl=ph
	ph=(beta_sum+beta_a(i)+beta_d(i)+beta_pp(i)+beta_pm(i))/alpha0
	
	if (r2.ge.pl.and.r2.lt.ph) then				
	p = p+1
	if (Bstate(i).eq.6) then
	kfpm(i) = 0.d0
	kfpam(i) = 0.d0
	kp(i)=kp0
	km(i)=km0
	kcp(i)=kcp0
	kfp(i)=kfp0
	kfpp(i)=kpp0
	Bstate(i) = 1
	elseif (Bstate(i).eq.5) then
	kfpm(i) = 0.d0
	kfpam(i) = 0.d0
	kp(i)=factor1*kp0
	kfpp(i)=factor2*kpp0
	kfp(i)=0.d0
	kcp(i)=0.d0
	kfcp(i)=kfcp0
	km(i)=km0
	Bstate(i) = 2
	endif
	Rid=17
	endif

!---------------- profilin-actin unbinding (B->B-1) ------------------

	pl=ph
	ph=(beta_sum+beta_a(i)+beta_d(i)+beta_pp(i)+beta_pm(i)+beta_pam(i))/alpha0
	
	if (r2.ge.pl.and.r2.lt.ph) then				
	l(i) = l(i)-1
	pa = pa+1
	if (Bstate(i).eq.6) then
	kfpm(i) = 0.d0
	kfpam(i) = 0.d0
	kp(i)=kp0
	km(i)=km0
	kcp(i)=kcp0
	kfp(i)=kfp0
	kfpp(i)=kpp0
	Bstate(i) = 1
	elseif (Bstate(i).eq.5) then
	kfpm(i) = 0.d0
	kfpam(i) = 0.d0
	kp(i)=factor1*kp0
	kfpp(i)=factor2*kpp0
	kfp(i)=0.d0
	kcp(i)=0.d0
	kfcp(i)=kfcp0
	km(i)=km0
	Bstate(i) = 2
	endif
	rcount = rcount + 1
	Rid=18
	endif



!---------------- capping ------------------

	pl=ph
	ph=(beta_sum+beta_a(i)+beta_d(i)+beta_pp(i)+beta_pm(i)+beta_pam(i)+beta_cp(i))/alpha0

	if (r2.ge.pl.and.r2.lt.ph) then
	cap=cap-1
	kcm(i)=kcm0
	kcfp(i)=kcfp0
	kfpp(i)=0.d0
	kp(i)=0.d0
	km(i)=0.d0
	kcp(i)=0.d0
	kfp(i)=0.d0
	Bstate(i) = 3
	Rid=19
	write(900,*)t,'CP+1', Nc-cap, 'index=', i

	if((Nc-cap + Nfr-frm-df-tf-ff).gt.nf + cpfr) then
	write(*,*) '*** ERROR *** UNPHYSICAL BINDING',i,Bstate(i),Rid
	GOTO 100
	endif

	endif


!---------------- de-capping ------------------

	pl=ph
	ph=(beta_sum+beta_a(i)+beta_d(i)+beta_pp(i)+beta_pm(i)+beta_pam(i)+beta_cp(i)+beta_cm(i))/alpha0

	if (r2.ge.pl.and.r2.lt.ph) then
	cap=cap+1
	kcm(i)=0.d0
	kcfp(i)=0.d0
	kp(i)=kp0
	km(i)=km0
	kcp(i)=kcp0
	kfp(i)=kfp0
	kfpp(i)=kpp0
	Bstate(i) = 1
 	Rid=20
	write(900,*)t,'CP-1', Nc-cap, 'index=', i
	endif

!---------------- formin binding ------------------

	pl=ph
	ph=(beta_sum+beta_a(i)+beta_d(i)+beta_pp(i)+beta_pm(i)+beta_pam(i)+beta_cp(i)+beta_cm(i)+beta_fp(i))/alpha0

	if (r2.ge.pl.and.r2.lt.ph) then
	kp(i)=factor1*kp0
	kfpp(i)=factor2*kpp0
	kfp(i)=0.d0
	kcp(i)=0.d0
	kfm(i)=kfm0
	kfcp(i)=kfcp0
	frm=frm-1
	Bstate(i) = 2
	Rid=21
	write(900,*)t,'FR+1', Nfr-frm, 'index=', i

	if((Nc-cap + Nfr-frm-df-tf-ff).gt.nf + cpfr) then
	write(*,*) '*** ERROR *** UNPHYSICAL BINDING',i,Bstate(i),Rid
	GOTO 100
	endif

	endif


!---------------- formin unbinding ------------------

	pl=ph
	ph=(beta_sum+beta_a(i)+beta_d(i)+beta_pp(i)+beta_pm(i)+beta_pam(i)+beta_cp(i)+beta_cm(i)+beta_fp(i)+beta_fm(i))/alpha0

	if (r2.ge.pl.and.r2.lt.ph) then
	kp(i)=kp0
	kfpp(i)=kpp0
	kfp(i)=kfp0
	kcp(i)=kcp0
	kfp(i)=kfp0
	kfm(i)=0.d0
	kfcp(i)=0.d0
	frm=frm+1
	Bstate(i) = 1
	Rid=22
	write(900,*)t,'FR-1', Nfr-frm, 'index=', i
	endif

!============= formin capping interaction ==============


!---------------- BFC formation (BF+C -> BFC) ---------------

	pl=ph
	ph=(beta_sum+beta_a(i)+beta_d(i)+beta_pp(i)+beta_pm(i)+beta_pam(i)+beta_cp(i)+beta_cm(i)+ &
           beta_fp(i)+beta_fm(i)+beta_fcp(i))/alpha0

!	write(100,*)'(BF+C -> BFC)',t,i,pl,ph,r2

	if (r2.ge.pl.and.r2.lt.ph) then
	kp(i)=0.d0
	km(i)=0.d0
	kfpp(i)=0.d0
	kcp(i)=0.d0
	kcm(i)=0.d0
	kfp(i)=0.d0
        kfm(i)=0.d0
	kfcp(i)=0.d0
	kfcm(i)=kfcm0
	kfcm2(i)=kfcm20
	cap=cap-1
	Bstate(i) = 4
	cpfr=cpfr+1
	Rid=23
	endif

!---------------- BFC dissociation ------------
!---------------- BFC -> BF + C ---------------
	pl=ph
	ph=(beta_sum+beta_a(i)+beta_d(i)+beta_pp(i)+beta_pm(i)+beta_pam(i)+beta_cp(i)+beta_cm(i)+ &
           beta_fp(i)+beta_fm(i)+beta_fcp(i)+beta_fcm(i))/alpha0
	
!	write(100,*)'(BFC -> BF + C)',t,i,pl,ph,r2

	if (r2.ge.pl.and.r2.lt.ph) then
	kp(i)=factor1*kp0
	kfpp(i)=factor2*kpp0
	km(i)=km0
       	kfp(i)=0.d0
   	kfm(i)=kfm0
	kcp(i)=0.d0
	kcm(i)=0.d0
	kfcp(i)=kfcp0
	kfcm(i)=0.d0
	kfcm2(i)=0.d0
	cap=cap+1
	Bstate(i) = 2
	cpfr=cpfr-1
	Rid=24
	endif

!---------------- BFC -> BC + F ---------------
	pl=ph
	ph=(beta_sum+beta_a(i)+beta_d(i)+beta_pp(i)+beta_pm(i)+beta_pam(i)+beta_cp(i)+beta_cm(i)+ &
           beta_fp(i)+beta_fm(i)+beta_fcp(i)+beta_fcm(i)+beta_fcm2(i))/alpha0

!	write(100,*)'(BFC -> BC + F)',t,i,pl,ph,r2

	if (r2.ge.pl.and.r2.lt.ph) then
	kp(i)=0.d0
	km(i)=0.d0
	kcp(i)=0.d0
	kcm(i)=kcm0
	kfp(i)=0.d0
        kfm(i)=0.d0
	kcfp(i)=kcfp0
	kfcm(i)=0.d0
	kfcm2(i)=0.d0
	frm=frm+1
	Bstate(i) = 3
	cpfr=cpfr-1
	Rid=25
	endif

!---------------- BCF formation (BC + F -> BCF) ---------------
	pl=ph
	ph=(beta_sum+beta_a(i)+beta_d(i)+beta_pp(i)+beta_pm(i)+beta_pam(i)+beta_cp(i)+beta_cm(i)+ &
           beta_fp(i)+beta_fm(i)+beta_fcp(i)+beta_fcm(i)+beta_fcm2(i)+beta_cfp(i))/alpha0

!	write(100,*)'(BC + F -> BCF)',t,i,pl,ph,r2

	if (r2.ge.pl.and.r2.lt.ph) then
	kcm(i)=0.d0
	kcfp(i)=0.d0
	kcfm(i)=kcfm0
	kcfm2(i)=kcfm20
	frm=frm-1
	Bstate(i) = 4
	cpfr=cpfr+1
	Rid=26
	endif

!---------------- BCF dissociation ------------
!---------------- BCF -> BC + F ---------------
	pl=ph
	ph=(beta_sum+beta_a(i)+beta_d(i)+beta_pp(i)+beta_pm(i)+beta_pam(i)+beta_cp(i)+beta_cm(i)+ &
           beta_fp(i)+beta_fm(i)+beta_fcp(i)+beta_fcm(i)+beta_fcm2(i)+beta_cfp(i)+beta_cfm(i))/alpha0

!	write(100,*)'(BCF -> BC + F)',t,i,pl,ph,r2

	if (r2.ge.pl.and.r2.lt.ph) then
	kcm(i)=kcm0
	kcfp(i)=kcfp0
	kcfm(i)=0.d0
	kcfm2(i)=0.d0
	frm=frm+1
	Bstate(i) = 3
	cpfr=cpfr-1
	Rid=27
	endif

!---------------- BCF -> BF + C ---------------
	pl=ph
	ph=(beta_sum+beta_a(i)+beta_d(i)+beta_pp(i)+beta_pm(i)+beta_pam(i)+beta_cp(i)+beta_cm(i)+ &
           beta_fp(i)+beta_fm(i)+beta_fcp(i)+beta_fcm(i)+beta_fcm2(i)+beta_cfp(i)+beta_cfm(i)+ &
           beta_cfm2(i))/alpha0

!	write(100,*)'(BCF -> BF + C)',t,i,pl,ph,r2

	if (r2.ge.pl.and.r2.lt.ph) then
	kp(i)=factor1*kp0
	kfpp(i)=factor2*kpp0
	km(i)=km0
        kfm(i)=kfm0
	kfcp(i)=kfcp0
	kcfm(i)=0.d0
	kcfm2(i)=0.d0
	cap=cap+1
	Bstate(i) = 2
	cpfr=cpfr-1
	Rid=28
	endif


!----------------- pass ass propensity sum as beta_sum -------------------

	beta_sum = beta_sum+beta_a(i)+beta_d(i)+beta_pp(i)+beta_pm(i)+beta_pam(i)+beta_cp(i)+beta_cm(i)+ &
                   beta_fp(i)+beta_fm(i)+beta_fcp(i)+beta_fcm(i)+beta_fcm2(i)+beta_cfp(i)+beta_cfm(i)+ &
                   beta_cfm2(i)

	enddo	!! nf-loop ends
	endif

!	write(200,*) t,m,dm,tm,fm, l(1), l(2)

!------------------ t -> t+tau ---------------

	t = t + tau

!----------------- write ---------------------

	if(t/(delt+nk*delt).ge.1.d0.and.t.ge.t0) then
!	write(202,*) t, rcount, float(sum(l)+4*fm)/V, float(nf)/V
	do i=1,nf_max
	if(abs(l(i)).gt.er) then
	write(201,*) l(i)
	endif
	enddo
	endif

	if(t/(delt+nk*delt).ge.1.d0) then
	write(200,*) t,float(sum(l))/nf, float(sum(l))/N				! t,(l(i), i=1,nf),m,dm,tm
	write(100,*) t, nf, filcount , Nc-cap, p, pa, Nfr-frm, rcount,m,dm,tm
!	if (nf.gt.0) write(101,*) t, float(sum(l)+4*fm)/nf
	nk=nk+1
	endif


        if(t.ge.scount*sdur) then
        seed = nint(seed/(1.001d0))
	scount = scount+1
	write(*,*) 'seed=',seed(1)
        endif


!	--------------- ERROR message --------------

	do i=1,nf_max
	if(l(i).lt.-er) then
	write(*,*) '*** ERROR *** NEGATIVE LENGTH',i,l(i),Rid
	GOTO 100
	endif
	enddo

	if(sum(l)+pa+m+2*dm+3*tm+4*fm+2*df+3*tf+4*ff.ne.N) then
	write(*,*) '*** ERROR *** POOL NOT CONSERVED', N-(sum(l)+pa+m+2*dm+3*tm+4*fm+2*df+3*tf+4*ff), Rid
	GOTO 100
	endif

	enddo

!------------------ time loop ended -----------------------

	write(*,*)'nf=',nf, 'maxsize=',maxval(l/370.d0)


100	WRITE(*,*) "CODE ENDS"

	stop
	end program code



	!open(unit=100,file='data/m.txt',status='unknown')
	!open(unit=101,file='data/ml.txt',status='unknown')
	!open(unit=200,file='data/l.txt',status='unknown')
	!open(unit=201,file='data/pl.txt',status='unknown')
	!open(unit=202,file='data/state.txt',status='unknown')
	!open(unit=500,file='ran.txt',status='unknown')





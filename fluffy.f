      module interp
        implicit none
        include 'physconst.inc' !Rgas is in erg K-1 mol-1 (convert to J by x10^-7)
        include 'cloudvars.inc'
        include 'declaration.inc'
        contains
!----------------------------------------------------------------------!
!       Find temp at given alt by interpolation                        !
!----------------------------------------------------------------------!
      function TatZ(Pz,Pin,Tin) result(Tz)
      real*8, intent(in) :: Pz ! input
      real*8,intent(in), dimension (1 : 21) :: Pin, Tin ! input
      real*8 :: Tz ! output
      integer :: l ! only a location
      l = minloc(abs(Pin-Pz),1)
      if ( abs(Pin(l+1)-Pz) <= abs(Pin(l-1)-Pz) ) then
        Tz = Tin(l) + (Pz-Pin(l))*(Tin(l+1)-Tin(l))/(Pin(l+1)-Pin(l))
      elseif ( abs(Pin(l+1)-Pz) > abs(Pin(l-1)-Pz) ) then
        Tz = Tin(l-1)+(Pz-Pin(l-1))*(Tin(l)-Tin(l-1))/(Pin(l)-Pin(l-1))
      end if
      end function TatZ
!----------------------------------------------------------------------!
!       Find local temperature gradient                                !
!----------------------------------------------------------------------!
      function GRADI(Pz,Pin,Tin) result(tau)
      real*8, intent(in) :: Pz ! input
      real*8,intent(in), dimension (1 : 21) :: Pin, Tin ! input
      real*8 :: tau ! output
      integer :: l ! only a location
      l = minloc(abs(Pin-Pz),1)
!      temp = Tmes(l)
      if ( abs(Pin(l+1)-Pz) <= abs(Pin(l-1)-Pz) ) then
        tau = (Tin(l+1)-Tin(l))/(Pin(l+1)-Pin(l))
      elseif ( abs(Pin(l+1)-Pz) > abs(Pin(l-1)-Pz) ) then
        tau = (Tin(l)-Tin(l-1))/(Pin(l)-Pin(l-1))
      end if
      end function GRADI
!----------------------------------------------------------------------!
!    Find eddy diffusion coefficient at given alt                      !
!    This is Equation 5 from Ackerman & Marley 2001                    !
!----------------------------------------------------------------------!
      function EDDY(H, L, ro) 	  result(K)
      real*8, intent(in) :: H, L, ro
      real*8, parameter :: F = (stefco*1.D-3)*(Tf)**4.        ! Flux [J/m2 s]
      real*8, parameter :: rg = Rgas*1.D-7
      real :: K
      K = H*1.D5/3.*((L/H)**(4./3.))*(1.D3*rg*F/(mu*ro*cp))**(1./3.)
      end function EDDY
!----------------------------------------------------------------------!
!       A function to calculate the saturation vapour mixing ratio for !
!       a given temperature and ambient pressure                       !
!----------------------------------------------------------------------!
      function SatVap(T,P)  result(Qs)
        real*8, intent(in) :: T, P  ! input temperature, pressure    !
        real*8 :: Qs    ! output is saturation vapour MR !
        qs = (exp(10.53 - 2161./T - 86596./T**2.))/P
      end function SatVap
!----------------------------------------------------------------------!
!	find r_w [m]
!----------------------------------------------------------------------!
      function AR_W(b, c) result(rw)
         real*8, intent(in) :: b, c
         real*8 :: rw
      	rw = (1./2.)*(sqrt(b**2-4*c)-b)
      end function AR_W
!----------------------------------------------------------------------!
      end module interp

      program fluffy
        use interp
        implicit none

      open (unit=20,file='pressure.txt',action="read",status="old")
      read(20, *) P;
      close(UNIT=20)
      open (unit=20,file='temp.txt',action="read",status="old")
      read(20, *) T;
      close(UNIT=20)

!     SET SOME INITIAL VALUES
      z(1) = 0.
      Tz(1) = 165.
      Pz(1) = 1.
      Hs = (Rgas*1.D-7)*Tf/(mu*g)!*Tz(1)/(mu*g)
      Qv(1) = qbelow
      Qc(1) = 0.d0
      Qt(1) = Qv(1) + Qc(1)
!     NOW CALCULATE THE LOT
      open (unit=20,file="output.txt",action="write",status="replace")
      do i=1, layers-1
        z(i+1) = i*Dz
!       ESTIMATING PRESSURE & DENSITY USING HYDROSTATIC EQUILIBRIUM
        Pz(i+1) = exp(-z(i+1)/Hs)
        rho_a = rho_0*exp(-z(i+1)/HS)
!       GETTING TEMPERATURE FROM MEASURED PROFILE ar_wBY INTERPOLATION
        Tz(i+1) = TatZ(Pz(i+1),P,T)
!       FIND LOCAL SCALE HEIGHT ... or don't
!        Hs = (Rgas*1.D-7)*Tz(i+1)/(mu*g)
!       CALCULATE LOCAL TEMPERATURE LAPSE RATE dT/dz
!       USING
!       dT/dz = [dT/dP]*[dP/dz] = [dT/dP]*[-rho_a * g]
        Tau = GRADI(Pz(i+1),P,T)*(-rho_a*g)*1.d1
!       CALCULATE LOCAL MIXING LENGTH
        L = Hs*max(A,abs(Tau/Tau_a))
!       CALCULATE EDDY DIFFUSION COEFFICIENT
        K = EDDY(Hs,L,rho_a)
         if ( K < Kmin ) then
	        K = Kmin ! have prescribed a minimum value for eddy diff. coeff.
       	 endif
!       FIND THE CONVECTIVE VELOCITY FROM MIXING LENGTH THEORY
        w = (K/L)*(1.D-7)
!
!       AT LAST, THE GOOD BIT
!       SATURATION MIXING RATIO
        Qs(i+1) = SatVap(Tz(i+1),Pz(i+1))
!       ACTUAL MIXING RATIO OF VAPOUR => EQ 2 OF ACKERMAN & MARLEY 2001
        Qv(i+1) = min( Qv(i), (1.D0+Sc)*Qs(i+1) )
!       FIND AMOUNT OF "NEW" CONDENSATE BORN IN THIS LAYER
        cond = max( 0.d0, Qv(i)-(1.D0+Sc)*Qs(i+1) )
!       FIND RATE OF CHANGE OF MIXING RATIO WITH HEIGHT
        Dqt = -(frain*w*(Qc(i)+cond)/(K))*1.d4
!       UPDATE THE TOTAL MIXING RATIO
        Qt(i+1) = Qt(i) + Dqt*Dz*1.d3
!       FINALLY, SET THE NEW MIXING RATIO OF CONDENSATE
        Qc(i+1) = Qt(i+1)-Qv(i+1)
!        write (20,*) Pz(i+1),Tz(i+1),Qc(i+1),Qv(i+1),Qt(i+1),Qs(i+1)
!       NOW THE SECOND PART
        eta = ((kb*Tz(i+1)/eps)**0.16)/(pi*(d**2)*1.22)
        eta = sqrt(pi*em*kb*Tz(i+1)/1.d3)*eta
        eta = ((1.D2**2)*5./16.)*eta
!       MEAN FREE PATH
        MFP = (Rgas*1.D-7)*Tz(i+1)/(sqrt(2.)*pi*(d**2.)*Navog*Pz(i+1))
        MFP = 0.1*MFP
        rw = AR_W(1.26*MFP, -9.*w*eta/(2.*g*((rho_amm-rho_a)*1.d3)))
        Nre = 2.*rw*rho_a*w/eta
        write (20,*) Nre
      enddo
!
      close (20)
!      print*, 1.d1,1d1
      end program fluffy

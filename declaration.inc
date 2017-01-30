!---------------------------------------------------------------------!
!     CALCULATION STUFF
!---------------------------------------------------------------------!
      integer :: i
      real*8 :: Dqt ! rate of change of total mixing ratio with height [1/km]
!---------------------------------------------------------------------!
!     REAL, ACTUAL THINGS IN THE WORLD
!---------------------------------------------------------------------!
      real*8 :: rho_a ! atmospheric density [g/cm3]
      real*8 :: Hs ! atmospheric scale height [km]
      real*8, dimension (1 : 21) :: P, T ! input pressure/temperature profile [bars, K]
      real*8, dimension(1 : layers) :: z ! elevation [km]
      real*8, dimension(1 : layers) :: Pz, Tz ! interpolated pressure/temperature profile [bars, K]
      real*8 :: Tau ! local temperature lapse rate dT/dz ![K/km]
      real*8 :: L ! mixing length [km]
      real*8 :: K ! eddy diffusion coefficient [cm2/s]
      real*8 :: w ! convective velocity [m/s]
!---------------------------------------------------------------------!
!     THE MIXING RATIOS
!---------------------------------------------------------------------!
      real*8, dimension(1 : layers) :: Qs ! saturation vapour mixing ratio
      real*8, dimension(1 : layers) :: Qv ! vapour mixing ratio
      real*8, dimension(1 : layers) :: Qc ! condensate mixing ratio
      real*8, dimension(1 : layers) :: Qt ! total mixing ratio {Qt = Qv + Qc}
      real*8 :: cond
!---------------------------------------------------------------------!
!     NOT SURE WHERE THESE BELONG
!---------------------------------------------------------------------!
      real*8 :: eta
      real*8 :: MFP ! molecular mean free path
      real*8 :: rw
      real*8 :: Nre ! Reynold's number

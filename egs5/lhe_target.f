!-----------------------------------------------------------------------
! Step 1: Initialization
!-----------------------------------------------------------------------

      implicit none

!     Main EGS header file      
      include 'include/egs5_h.f'

!     bounds contains the COMMON block BOUNDS with the variables
!       ECUT: An array of region-dependent charged particle cutoff
!       energies in MeV.
!       PCUT: An array of region-dependent photon cutoff energies in
!       MeV.
      include 'include/egs5_bounds.f'
!     media contains the COMMON block MEDIA with the variables
!       NMED number of media being used
!       MEDIA Array containing the density of the media in g/cm3
!       IRAYLM Array of flags for turning on Rayleigh scattering in the
!       various regions.
!       CHARD Array used to define the dimensions in cm (default 1 cm)
!       MED Array containing the medium index for each region
      include 'include/egs5_media.f'
!     usersrc contains the COMMON block USERSC with the variable 
!       EMAXE: The maximum total energy of any electron in the
!       simulation
      include 'include/egs5_usersc.f'
!     edge contains the COMMON block EDGE2 with the variable IEDGFL     
      include 'include/egs5_edge.f'
!     misc contains the COMMON block MISC with the variable IMPACR
      include 'include/egs5_misc.f'
!     randomm contains the COMMON block RLUXDAT with the variables
!       INSEED: The initial seed
!       LUXLEV: Luxury level used by RANLUX
      include 'include/randomm.f'
!     useful contains the COMMON block USEFUL with the variable
!       RM: Rest mass of the electron
      include 'include/egs5_useful.f'
!
      include 'include/egs5_epcont.f'

!     Set the target thickness in cm. This value will be used by the
!     subroutine howfar
      common/geo/tgtdz
      real*8 tgtdz

      real*8 ei, ekin, xi, yi, zi, ui, vi, wi, wti, p
      integer iqi, iri

!     Variables needed for parsing of LHE file  
      character*80 line
      character*6 eventheader / '<event' /
      integer nup, idup, istup, mothup(2), icolup(2)
      real*8 pup(5)

!     Variables needed for writing an LHE file
      common/event/idp, xwg, scal, aqed, aqcd
      integer idp
      real*8 xwg, scal, aqed, aqcd


!     Initialize the array that will contain the names of the media. In
!     this case, the array is of length 1 since only the target material
!     needs to be defined.
!     NOTE: The names of the media need to be exactly 24 characters
!     long.
      integer i, j
      character*24 medarr(1)

      open(UNIT=6, FILE='egs5job.out', STATUS='UNKNOWN')
      open(UNIT=100, FILE='out.lhe', STATUS='UNKNOWN')
      write(100, 500) '<LesHouchesEvents version="1.0">'
500   format(A)
 
!-----------------------------------------------------------------------
! Step 2: pegs5-call
!-----------------------------------------------------------------------
!     Initialize some general variables
      call block_set                 

!     Define the media before calling PEGS5
      nmed=1
      medarr(1)='W-RAY                   '

      do j=1,nmed
        do i=1,24
          media(i,j)=medarr(j)(i:i)
        end do
      end do

      chard(1) = 0.0001

      call pegs5

!-----------------------------------------------------------------------
! Step 3: Pre-hatch-call-initialization
!-----------------------------------------------------------------------

!     The number of regions in the geometry. In this case, there will be
!     three regions: | Vacuum | W | Vacuum | 
      nreg = 3

!     Specify what type of media fills each region
      med(1) = 0
      med(2) = 1
      med(3) = 0

!     Terminate electron histories at 0.521 MeV in the W target
      ecut(2) = 0.521
!     Terminate photon histories at 0.001 MeV in the W target
      pcut(2) = 0.001
!     This turns on explicit modeling of K and L-edge fluorescent
!     photons.
      iedgfl(2) = 1
!     This turns on Rayleigh scattering
      iraylr(2) = 1
!     This turns on electron impact ionization
      impacr(2) = 1

!     --------------------------------------------------------
!     Random number seeds.  Must be defined before call hatch
!     or defaults will be used.  inseed (1- 2^31)
!     --------------------------------------------------------
      inseed = 1
!     The "luxury level" sets the level of tests for randomness.
      luxlev = 1

!     Initialize the ranlux random-number generator
      call rluxinit

!-----------------------------------------------------------------------
! Step 4:  Determination-of-incident-particle-parameters
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
! Step 5:   hatch-call
!-----------------------------------------------------------------------

!     Maximum total energy
      if (iqi.ne.0) then
        emaxe = ei
      else
        emaxe = ei + RM
      end if

      open(UNIT=KMPI,FILE='pgs5job.pegs5dat',STATUS='old')
      open(UNIT=KMPO,FILE='egs5job.dummy',STATUS='unknown')

      call hatch

      close(UNIT=KMPI)
      close(UNIT=KMPO)
!-----------------------------------------------------------------------
! Step 6:  Initialization-for-howfar
!-----------------------------------------------------------------------
      tgtdz = 0.002

!-----------------------------------------------------------------------
! Step 7:  Initialization-for-ausgab
!-----------------------------------------------------------------------

      open(UNIT=15, FILE='input.lhe', STATUS='OLD')
100   read(UNIT=15,FMT='(A)',END=250) line
      if (line(1:6).eq.eventheader) go to 150
      go to 100
150   continue
      write(100, 501) '<event>'
501   format(A)
      write(100, *) 0, idp, xwg, scal, aqed, aqcd

      read(UNIT=15,FMT=*,END=250) nup, idp, xwg, scal, aqed, aqcd
      print *, nup, idp, xwg, scal, aqed, aqcd
200   if (nup.eq.0) then
        write(100, 502) '</event>'
502   format(A)
        go to 100
      end if
!      print *, 'Total numer of particles left', nup
      read(UNIT=15,FMT=*,END=250) idup, istup, mothup, icolup, pup
      print *, idup, istup, mothup, icolup, pup
      nup = nup - 1

      if (istup.ne.1) go to 200
      if (idup.eq.11) then
!     Specify that the incident particle is an electron
        iqi = -1
      else if (idup.eq.22) then
        iqi = 0
      else
        go to 200
      end if

!     Set the coordinates of the incident particle
      xi = 0.0
      yi = 0.0
      zi = 0.0

!     Calculate the momentum of the particle
      p = sqrt(pup(4)**2 - pup(5)**2)
!     Total energy of the incident particle in MeV
      ei = pup(4)*1000.
      ekin = ei+iqi*RM
!     Set the direction cosines
      ui = pup(1)/p
      vi = pup(2)/p
      wi = pup(3)/p
!     Set the incident region to be 2 
      iri = 2
!     Weight factor in importance sampling
      wti = 1.0

      call shower(iqi, ei, xi, yi, zi, ui, vi, wi, iri, wti)
      go to 200

      close(UNIT=7) 

250   continue
      write(100, *) '</LesHouchesEvents>'
      stop
      end

!
!
!
      subroutine ausgab(iarg)

      implicit none
!     Main EGS header file
      include 'include/egs5_h.f'
!     The COMMON block STACK contains the following useful variables
!       x, y, z: Position of a particle 
!       u, v, w: Direction cosines of a particle
!       ir: Index of particle's current region
!       np: The particle currently being pointed to
      include 'include/egs5_stack.f'
      include 'include/egs5_useful.f'
      include 'include/egs5_epcont.f'

!     Arguments
      integer iarg

      real*8 kine, thetaz
      
      integer id, ist, spin
      integer i
  
      real p, px, py, pz, energy, mass

      if (iarg.ne.3) return

      mass = RM
      p = sqrt(e(np)**2-0.511**2)*0.001
      if (ir(np).eq.3) then
        if (iq(np).eq.0) then
          id = 22
          mass = 0
          p = e(np)*0.001
        else if (iq(np).eq.-1) then
          id = 11
        else if (iq(np).eq.1) then
          id = -11
        else
          stop
        end if

        ist = 1
        px = p * u(np)
        py = p * v(np)
        pz = p * w(np)
        energy = e(np) * 0.001
       
        write(100, *) id,ist,0,0,0, 0, px, py, pz, energy, mass,  0, 0
        print *, id, ist, 0, 0, 0, 0, px, py, pz, energy, mass,  0, 0
      end if

      return
      end
!
!
!

      subroutine howfar

      implicit none

!     Main EGS header file
      include 'include/egs5_h.f'
!     COMMONs required by EGS5 code
      include 'include/egs5_epcont.f'
!     The COMMON block STACK contains the following useful variables
!       x, y, z: Position of a particle 
!       u, v, w: Direction cosines of a particle
!       ir: Index of particle's current region
!       np: The particle currently being pointed to
      include 'include/egs5_stack.f'

!     The target thickness. 
      common/geo/tgtdz
      real*8 tgtdz
    
!     The distance from the particle to the target boundary
      real*8 d

!     First, check if the particle has exited the target i.e. not in
!     region 2. If it's outside the target, discard the particle.
      if ( ir(np).eq.3 ) then
!     Setting idisc to a non-zero value tells the simulation to
!     discard the particle.       
         
        idisc = 1
        return
!     Now consider the case where the particle is in the target and
!     traversing forward
      else if ( ir(np).eq.2 ) then
        if ( w(np).gt.0.0) then
!     Calculate the distance from the particle to the downstream face of
!     the target (d).  If the  distance is greater than the step size
!     (ustep) then continue stepping. Otherwise, set ustep to d and set
!     the new region (irnew) to 3 to indicate the particle will exit.
          d = (tgtdz - z(np))/w(np)
          if (d.gt.ustep) then
            return
          else
            ustep = d
            irnew = 3
            return
          end if
!     In the case where the particle is in the target and traversing
!     backwards, let the particle continue stepping if the step size is
!     smaller than the distance to the upstream face of the target.
        else if ( w(np).lt.0.0 ) then
          d = -z(np)/w(np)
          if ( d.gt.ustep ) then
            return
          else 
            ustep = d
            irnew = 1
            return
          end if
        else if ( w(np).eq.0.0 ) then
          return
        end if
!     If the particle is upstream of the target, then it's the incident
!     particle.  In this case, as long as the particle is pointed
!     towards the target, set the position to the upstream face of the
!     target and let it propagate.
      else if ( ir(np).eq.1 ) then
        if ( w(np).gt.0.0 ) then
          ustep = 0.0
          irnew = 2
          return
!     The particle must have been reflected so discard it
        else
          idisc = 1
          return
        end if
      end if 
      end

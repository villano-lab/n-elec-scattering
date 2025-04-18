c***********************************************************************
c*                                                                     *
c*  SOURCES: Version 4A       (Last Edited: 11/26/21)                  *
c*                                                                     *
c***********************************************************************


c***********************************************************************
c                                                                      *
c Written By:                                                          *
c                                                                      *
c   W.B.Wilson   (LANL, T-2, Nuclear Theory & Applications)            *
c   R.T.Perry    (LANL, NIS-7, Safeguard Systems)                      *
c   J.E.Stewart  (LANL, NIS-5, Safeguards, Science, & Technology)      *
c   W.S.Charlton (LANL, NIS-7, Safeguard Systems)                      *
c   T.H.Brown    (LANL, ESH-12, Policy and Program Analysis)           *
c   T.A.Parish   (Texas A&M University)                                *
c                                                                      *
c Nuclear Data Contributors:                                           *
c                                                                      *
c   T.R.England  (LANL, T-2, Nuclear Theory & Applications)            *
c   E.D.Arthur   (LANL, NMSM-PO, Nuclear Component Readiness)          *
c   D.G.Madland  (LANL, T-2, Nuclear theory & Applications)            *
c                                                                      *
c***********************************************************************
c***********************************************************************
c                                                                      *
c Overview:                                                            *
c                                                                      *
c   The sources code calculates (alpha,n), spontaneous fission,        *
c and delayed neutron source magnitudes and spectra using available    *
c library data.  Data sources are identified on leading comment        *
c records in the four input data library files.                        *
c                                                                      *
c   (alpha,n) neutron sources are calculated with 4th.-order fits      *
c to stopping cross sections and linear-linear cross section interpo-  *
c lation points. Spectra are calculated with the approximation of      *
c whitmore (phys.rev.78,799(1950), assuming isotropic neutron emission *
c in the center-of-mass system.  Branchings to product levels were     *
c produced with model code calculations.                               *
c                                                                      *
c   Spontaneous fission sources are calculated with s.f. branching and *
c s.f. nu-bar values gathered from available sources.  Spectra are     *
c calculated with watt spectrum parameters.                            *
c                                                                      *
c   Delayed neutron sources and spectra are calculated with pn (i.e.,  *
c branching) and independent d.n. spectra of england,et al.,           *
c (la-ur-82-841)                                                       *
c                                                                      *
c***********************************************************************


c***********************************************************************
c
c Updates by Piotr Krawczun
c
c 26 Nov 2021:
c 
c (645) The Big Loop iterates over every target nuclide. A write statement was
c added to show the loop count so that we know the program didn't freeze.
c
c Dimensions of arrays totlev, sl, etopl, ebotl were increased from 20 to 100
c in order to avoid overindexing.
c
c (894) SOURCES linearly interpolates cross sections from tape3 to estimate
c the value at a specific energy that is not in the library. A warning will
c be now printed out on the screen whenever calculating the slope involves
c division by zero. No fix for this problem has been implemented yet but
c we want at least to know if and when it happens. Values such as NaN or
c Infinity in the output files are symptomatic of this problem.
c
c Another modification allows to work with multiple input files. They have to
c be listed in the file called 'filelist'. First line is the number of files in the
c list and then each file name must be in a new line. No spaces are allowed in the
c file names. Max file name length is 50. If the filelist is not found the code goes
c into default and looks for 'tape1'. This upgrade is therefore optional. Max number of
c files is 1000.
c***********************************************************************

c***********************************************************************
c Description of File Structure
c***********************************************************************
c tape1 = user input
c tape2 = stopping cross section expansion coefficients library.
c tape3 = targets (alpha,n) cross section library.
c tape4 = targets (alpha,n) product level branching library.
c tape5 = sources decay data library.
c tape6 = neutron source magnitudes output file.
c tape7 = absolute neutron spectra output file.
c tape8 = normalized neutron spectra output file.
c tape9 = neutron source spectra by product level output file.
c outp  = summary of input and output.  

c***********************************************************************
c Definition of Variables (List May be Incomplete)
c***********************************************************************
c
c a          = watt spec. parameter: watt(e)=w*exp(-e/a)*sinh(sqrt(be))
c adef       = default value of watt spectrum parameter a
c alam       = source nuclide decay constant (/sec).
c alph       = mass of alpha particle (4).
c aneut      = mass of neutron (1).
c apro       = mass of product nuclide.
c aps        = source rate of alphas at eal(l) from source nuclide k
c aq(k)      = density (atoms/cc) of source nuclide k
c at         = input fraction of target atoms that are nuclide i
c atar       = mass of target nuclide.
c azm(j)     = fraction of target material atoms that are element j
c b          = watt spec. parameter: watt(e)=w*exp(-e/a)*sinh(sqrt(be))
c barnu      = source nuclide spontaneous fission nu-bar value
c bdef       = default value of watt spectrum parameter b
c beamn      = neutrons/sec per microamp of beam on target
c bfsf       = source nuclide decay branching fraction for s.f.
c bx         = branching fraction of alphas at ea reacting with target i
c              and producing product level il
c c(jt)      = solid stopping cross section coef. for terms jt=1,5;
c              gas stopping cross section coef. for terms jt=6,10.
c c1         =          erf() arguments in
c c2         =          watt spectrum integral
c c3         =          calculation.
c c4         =
c c1s        =          exp() arguments in
c c2s        =          watt spectrum integral
c c3s        =          calculation.
c c4s        =
c cx(m)      = (alpha,n) cross section (mb.) at alpha energy pt. m
c czm(jt,j)  = scx(j) coef. jt of elemental constituent
c dcx        = ln of element j stopping cross section.
c de         = that part (mev) of dele in neutron energy group n
c dea        = alpha energy grid width (mev) used for beam or source i
c dele       = width (mev) of neutron energy range possible from alphas
c              reacting at ea with target i populating product nuclide l
c e(ip)      = energy at point ip of library cross-section evaluation
c e90        = minimum alpha energy (mev) at which neutrons may be produ
c              at 90 degrees in the center-of-mass system in reactions
c              with target i populating product level il.
c ea         = midpoint of alpha energy group
c eal(l)     = energy (mev) of source alpha l
c eamax      = maximum alpha energy (mev) of source nuclide k
c eamin      = 0.001 mev or, if greater, target (alpha,n) threshold
c ebeam      = alpha beam energy (mev), used only if nq=0
c ee(m)      = lower energy bound (mev) of alpha energy group m,
c              groups numbered in ascending alpha energy.
c eft        = collection of erf() terms used in watt spec. integral
c el(il)     = excitation energy (mev) of product nuclide level il
c elg        = ln. of alpha energy
c en(n)      = upper energy bound (mev) of neutron energy group n,
c              groups numbered in descending neutron energy
c enmax      = maximum energy of neutrons produced by alphas at ea
c              reacting with target i and populating product level il.
c enmin      = minimum energy (mev) of neutrons produced by alphas at ea
c              reacting with target i and populating level il.
c enrsep     = cross section linear interpolation intercept
c ep(ip)     = energy value ip at which branching fractions
c              f(il,ip) are evaluated for target i
c f(il,ip)   = fraction of (alpha,n) reactions with target i at energy
c              e(ip) resulting in the production of product level il.
c fact       = factor used in the calculation of p(m) to account
c              for units of ee,cx,and scx.
c fal(l)     = fraction of source nuclide decays by any mode resulting
c              in the emission of an alpha at eal(l).
c fdng(nd)   = fractional delayed neutron spectrum contribution in
c              input 100kev-wide bin structure.
c fwatt      = fraction of source k s.f. watt spectrum included in
c              neutron energy group structure
c gtmg       = sum of (alpha,n) + s.f. + d.n. multigroup neutron
c              spectra values for all sources and targets
c gtmga      = sum of (alpha,n) multigroup neutron spectra values for al
c              sources and targets
c gttqan     = total (alpha,n) neutron source from all sources and targe
c gtsan(n)   = group n spectrum value of total (alpha,n) spectrum
c i          = target nuclide index
c icon       = 0, indicates last library title card
c            > 0, indicates more library title card(s) to follow.
c id         = 1, calculate source magnitudes only
c            = 2, calculate source spectra also
c idnnq      = number of source nuclides that are d.n. precursors.
c idq        = source nuclide ident from library(lsq+10*laq+10000*lzq)
c idt        = target nuclide ident from input(lst+10*lat+10000*lzt)
c il         = product nuclide level index(1=ground,2=1st excited,etc.)
c ip         = index for input cross section or level branching data
c isfnq      = number of source nuclides decaying by s.f.
c isg        = 0, use solid stopping cross sections,
c            = 1, use gas stopping cross sections if available,
c              otherwise use solid stopping cross sections.
c iz         = element z value
c j          = elemental stopping cross section constituent index.
c jdt        = target ident read from library
c jl         = number of product nuclide levels
c jn         = temp element index used in ordering
c jp         = number of product level branching data points
c jps        = number of points in library (alpha,n) cross-section
c              evaluation for nuclide i.
c jq(k)      = source nuclide k ident read from user input;
c            = state+10*mass+10000*z
c jsm(iz)    = chemical symbol of element with z=iz
c jt         = scx coefficient index
c ajwd1      =       table i column heading words for incident
c ajwd2      =       alpha beam or decay alpha source
c jzm(j)     = z value of stopping cross section elemental constituent
c k          = source nuclide index
c kn         = temp source nuclide index used in ordering
c l          = source nuclide alpha spectrum index
c laq        = source nuclide a value
c lat        = target nuclide a value
c lp         = index of level branching data frame containing ea
c lsq        = source nuclide state (1=gnd.,1=1st.iso.,etc.)
c lst        = target nuclide state (0=gnd.,1=1st.iso.,etc.)
c lzq        = source nuclide z value
c lzt        = target nuclide z value
c m          = alpha energy group index
c mm         = alpha energy group containing energy of source nuclide k
c              alpha eal(l).
c amq        = source nuclide isomer tag (blank or m)
c amt        = target nuclide isomer tag (blank or m)
c n          = neutron energy group index
c nag        = number of alpha energy groups
c nal        = number of alpha energies in source nuclide k library spec
c nd         = delayed neutron energy structure group index.
c ndn        = number of 100kev-wide delayed neutron energy groups
c              required to describe the delayed neutron spectrum
c              associated with precursor k.
c nh         = nng/2
c nng        = number of neutron spectrum energy groups
c nq         = number of source nuclides to be evaluated;if=0,source is
c              incident alpha particle beam at ebeam.
c nt         = number of target nuclides specified in user input;
c              if=0, neutrons from s.f. only.
c nz         = number of stopping cross section elemental constituents
c p(m)       = probability of an alpha at ee(m) reacting with target i
c              and populating product level il rather than slowing
c              below eamin.
c pi         = 3.14159
c pval       = interpolated value of p at eal(l) or at ebeam.
c q          = target nuclide i (alpha,n) reaction q value (mev).
c qan        = neutron source rate (#/cc/sec) due to source nuclide k
c              alphas at eal(l) reacting at some energy above eamin with
c              target i and populating product level il.
c qlev       = target nuclide i (alpha,n) reaction q value (mev) for
c              production of product nuclide level il.
c qsf        = neutron source k rate due to s.f. of source k
c r(m)       = ratio of (alpha,n) cross section to alpha stopping cross
c              at ee(m); this ratio is the integrand of p at ee(m).
c rr(m)      = calculated fraction of target i product level il
c              reactions of source k alphas l occuring in alpha energy g
c s(n)       = accumulated group n neutron spectrum contribution from
c              source k alpha l on target i
c sa         = sqrt(a), used in watt spectrum integration
c san(n)     = accumulated neutron group n (alpha,n) spec. contribution
c              for source k alphas on target i.
c sbaso4     = sqrt(b * a**2 / 4)
c se1        = sqrt of maximum energy of neutron structure
c seh        = sqrt of neutron group upper energy bound
c sel        = sqrt of neutron group lower energy bound
c sen1       = sqrt of minimum energy of neutron structure
c senmax     = sqrt of maximum neutron energy from reactions at ea
c senmin     = sqrt of minimum neutron energy from reactions at ea
c scx(m)     = stopping cross section (ev/10**15 atoms/cm**2) at ee(m).
c sdn(n)     = group n contribution to delayed neutron spectrum of
c              nuclide k.
c sl(n,il)   = accumulated group n neutron spectrum from production of
c              target i product nuclide level il by source k alphas.
c slope      = cross section linear interpolation slope
c smga       = sum of source k on target i (alpha,n) multigroup
c              neutron spectrum values
c smgs       = check sum of source k multigroup neutron spectrum values
c spba       = sqrt(pi * b * a)
c ssf(n)     = group n source k s.f. neutron spectrum contribution.
c sbtqan     = total (alpha,n) source for target i / source k
c term1      =       combinations of mass and energy
c term2      =       quantities used in the calculation of
c term3      =       values of enmin, and enmax.
c thre       = threshold energy (mev) of (alpha,n) reactions on target i
c              populating product nuclide energy level il.
c tmga       = sum of target i (alpha,n) multigroup neutron spectra
c              values for all sources
c tmgs       = check sum of multigroup spectra values for all sources
c totqan     = accumulated total of qan values all sources on target i.
c totqsf     = accumulated total of qsf values.
c ts(n)      = total (alpha,n)+s.f. neutron spectrum in group n
c tsan(n)    = total group n target i (alpha,n) neutron spectrum
c tsdn(n)    = group n contribution to total delayed neutron spectrum
c tssf(n)    = total s.f. neutron spectrum in group n
c w          = watt spec. parameter: watt(e)=w*exp(-e/a)*sinh(sqrt(be))
c x(ip)      = cross section value at library evaluation point ip.
c **********************************************************************

c=======================================================================
c  SOURCES Main Program (6/97)
c=======================================================================

      program source

c----------------------
c  Common Data Storage
c----------------------
      character*7 title
      character*2 jsm(105)
      dimension title(11), amass(105)
      common /names/ jsm
      common /masses/ amass
c The following variables will be needed to read filelist
      integer nof ! number of files !pk
      integer count ! generic loop counter !pk
      integer leng ! length of a filename !pk
      logical exists ! does filelist exist? !pk
      character*50 fnames !names of files !pk
      character letter !pk
      dimension fnames(1000) !up to 1000 input files !pk

c------------------
c  Read file names !pk
c------------------
  90  format(i3)
 100  format(a50)

      inquire(file='filelist', EXIST=exists)
      if (exists) then
       open(unit=99,file='filelist',status='old')
       read(99,90) nof
       if (nof.eq.0) then
        goto 213
       endif
 150   read(99,100,end=200) (fnames(n),n=1,nof)
       goto 210
 200   stop 'Unexpected end of file reached.'
 210   continue
       close(unit=99)
      else
 213   write(0,*) 'No filelist found. tape1 will be used as input.'
       nof = 1;
       fnames(1)='tape1'
      endif
      close(unit=99)

      do 300 count=1,nof
       do 290 jj=1,50 ! 50 is the length of that array
        letter = fnames(count)(jj:jj)
        if (letter.eq.' ') then
         leng=jj-1
         goto 295
        endif
 290   continue
 295  continue

      inquire(file=fnames(count)(1:leng), EXIST=exists)
      if (exists) then
       write(0,*) 'Input file: ', fnames(count)(1:leng)
      else
       write(0,*) 'Input file ', fnames(count)(1:leng),
     1' not found.'
       goto 300 ! skip non-existing input files
      endif

c-----------------
c  File Structure
c-----------------
      open(unit=1,file=fnames(count)(1:leng),status='old')
      open(unit=2,file='tape2',status='old')
      open(unit=3,file='tape3',status='old')
      open(unit=4,file='tape4',status='old')
      open(unit=5,file='tape5',status='old')
      open(unit=6,file=fnames(count)(1:leng)//'_tape6',
     1status='unknown')
      open(unit=7,file=fnames(count)(1:leng)//'_tape7',
     1status='unknown')
      open(unit=8,file=fnames(count)(1:leng)//'_tape8',
     1status='unknown')
      open(unit=10,file=fnames(count)(1:leng)//'_tape9',
     1status='unknown')
      open(unit=11,file=fnames(count)(1:leng)//'_outp',
     1status='unknown')

c------------------
c  Read From Tape1
c------------------
      read (1,17,end=99) title
      read(1,*)idd,id
      if (idd.eq.1) then
       call homo(title,idd,id)
      elseif (idd.eq.2) then
       call interf(title,idd,id)
      elseif (idd.eq.3) then
       call homo(title,idd,id)
      elseif (idd.eq.4) then
       call three(title,idd,id)
      else
       stop 'Error in idd input!'
      endif
   99 continue

c-------------------------
c  Close All Output Files
c-------------------------
      close(unit=1)
      close(unit=2)
      close(unit=3)
      close(unit=4)
      close(unit=5)
      close(unit=6)
      if (id.eq.1) then
       close(unit=7, status='delete')
       close(unit=8, status='delete')
       close(unit=10, status='delete')
      elseif (idd.eq.2) then
       close(unit=10, status='delete')
      else
       close(unit=7)
       close(unit=8)
       close(unit=10)
      endif

 300  continue

c--------------------
c  Format Statement
c--------------------
   17 format(11a7)

c------------------------------
c  End of Execution of SOURCES
c------------------------------
      Stop 'Normal Execution'
      End



c=======================================================================
c  Homogeneous Problem Subroutine (6/97)
c=======================================================================

      subroutine homo(title,idd,id)


c---------
c Storage
c---------
      character*7 title
      character*2 jsm(105)
      character*1 amt,amq,ajwd1
      character*4 ajwd2,ajwd3
      dimension title(11), jzm(20), azm(20), czm(9,20), c(14), jq(300),
     1aq(300),el(100), ep(200),f(100,200), eal(30), fal(30),e(1100),x(11
     2 00), ee(4001), scx(4001), en(751), san(750), ssf(750),
     3 ts(750), cx(4001), r(4001), rr(4001), p(4001), tsan(750
     4 ), fdng(1000), tsdn(750), sdn(750), gtsan(750), s(750), tssf(750)
     5 ,totlev(100),sl(750,100),etopl(100),ebotl(100) !pk
      common /names/ jsm
      common /masses/ amass(105)
c     dimension frclev(20) 
c el() and f(,) modified by v.t.

c-------------------
c Format Statements
c-------------------
   10 format(i3,11a7)
   15 format(3x,11a7)
   17 format(11a7)
   20 format(2i3)
   30 format(i8,e12.5)
   40 format(i3,8e8.3,f10.3,/,3x,5f8.3)
   50 format(' failed to find stopping cross section coefficients on',
     1' tape 2 for iz=',i3,' stop')
   60 format(6e12.5)
   70 format(/,1x,'Neutron Multigroup Structure (MeV)')
   80 format(1p8e10.3)
   81 format(35x,'Total (all groups): ',1pE10.3,' neutrons',a1,'sec-',
     1a4,a4,'.',/,31x,'Average Neutron Energy: ',1pE10.3,' MeV.')
   83 format(' Title    gttqan  ebaran   totqsf  ebarsf   totqdn  ',
     1 'ebardn   qtotal  ebrall')
   84 format(1x,a7,4(1pe10.3,7x))
   85 format(1x,a7,4(1pe10.3,0pf7.3))
   90 format (i3,5e12.5)
   95 format(2i6)
   96 format(/,26x,'Neutron Source Magnitudes',/,26x,25(1h~),/)
   97 format(/,23x,'Absolute Neutron Source Spectra',/,23x,31(1h~),/)
   98 format(/,22x,'Normalized Neutron Source Spectra',/,22x,
     133(1h~),/)
   99 format(/,16x,'Neutron Source Spectra by Nuclide Energy Level',
     1/,16x,46(1h~),/)
  100 format(//,1h1,30x,'Table I',/,31x,7(1h=),//,12x,'(alpha,n) ',
     1'Neutron Production by Target Per Source Alpha',///,17x,
     2'target    alpha   alpha   alphas/sec     p(e)     neuts/sec',
     3/,7x,'target  atom frac.  source  energy  ',a1,a4,a4,2x,
     4' neut/alpha  ',a1,a4,a4,/,1h+,6x,6(1h_),2x,10(1h_),2x,6(1h_),2x,
     5 6(1h_),3(12h  __________))
  102 format(1h+,77x,37hfractional ground/excited level split,/,
     1 1h+,77x,49(1h_))
  105 format(i8,i4)
  110 format (i8,2i4)
  120 format (8e10.3)
  130 format (52h failed to find alpha-n cross section on tape 3 for ,i8
     1 ,5h stop)
  140 format (i8,2i3)
  150 format (66h failed to find comp.nucl.neut.emiss.lev. branching on
     1tape 4 for ,i8,5h stop)
  160 format (56h failed to find alpha/s.f./dn source data on tape 5 for
     1 ,i8,5h stop)
  170 format (17h alpha energy of ,1pe12.5,62hexceeds current 6.5 mev.
     1 limit of most cross sections, stop.)
  180 format (7h target,i8,40h cross-section data not available at ee(
     1 ,i4,3h)= ,1pe12.5,6h stop.)
  190 format(7x,a2,i3,a1,1pe12.4,3x,5hbeam ,0pf8.3,1p3e12.4)
  200 format(33x,f8.3,1p3e12.4)
  205 format(1h+,(76x,7f7.5,/,1x))
  210 format(7x,a2,i3,a1,1pe12.4,2x,a2,i3,a1,0pf8.3,1p3e12.4)
  220 format (14h alpha energy ,1pe12.5,45h exceeds range of level branc
     1hing data, stop.)
  230 format (//,' (a,n) neutrons/sec-microamp ',f6.3,' mev a on ',a2,
     1i3,a1,' in target')
  240 format (//,' (a,n) neutrons/sec-cc from ',1pe12.5,' at/cc ',a2,
     1i3,a1,' alphas on ',a2,i3,a1,' in target')
  250 format (1h+,66x,10(1h_),/,55x,'Total:',4x,1pe12.4,/)
  251 format (1h+,66x,10(1h_),/,41x,'Total (this target):',
     14x,1pe12.4,/)
  252 format (1h+,66x,10(1h_),/,41x,'Total (all targets):',
     14x,1pe12.4,/)
  255 format(1h )
  260 format (///,' Total (alpha,n) neutron spectrum this target')
  261 format(//,' Neutron spectrum from ',a2,i3,a1,' alphas on ',a2,i3,
     1a1,' via the ',f6.2,'-MeV product level')
  262 format(//,' Neutron spectrum from ',F5.2,' MeV alphas on ',a2,i3,
     1a1,' via the ',f6.2,'-MeV product level')
  270 format(///,' Grand total (alpha,n) neutron spectrum, all targets,'
     1,1x,'all sources')
  280 format (////,1h2,35x,'Table II',/,36x,8(1h=),//,21x,
     1'Spontaneous Fission Neutron Production',///,8x,2(8h source ),
     2'atoms  dk constant  sf decay    nu    neutrons',/,8x,'nuclide
     3per cm**3    (/second)   branching   bar   sec/cm**3',
     4/,1h+,7x,7(1h_),2x,12(1h_),2x,11(1h_),2x,9(1h_),2x,5(1h_),2x,
     59(1h_))
  290 format (9x,a2,i3,a1,1p2e13.4,e12.3,0pf7.3,1pe11.3)
  300 format (39h no watt spectrum parameters found for ,a2,i3,a1,22h, d
     1efault values used:,/,3x,3ha= ,1pe12.5,5h, b= ,e12.5)
  310 format (//,' S.F. neutrons/cc from ',1pe12.5,' at/cc ',a2,i3,a1)
  320 format (1h+,61x,9(1h_),/,52x,'Total:',4x,1pe9.3)
  330 format (///,' Total S.F. neutron spectrum')
  340 format (////,1h3,35x,'Table III',/,36x,9(1h=),//,22x,
     1'Delayed Neutron Production',///,8x,2(8h source ),
     2'atoms  dk constant  branching neutrons'/,/,8x,
     3'nuclide   per cm**3     (/second)     (pn)sec/cm**3',
     4/,1h+,7x,7(1h_),2x,12(1h_),2x,11(1h_),2x,9(1h_),2x,9(1h_))
  350 format (9x,a2,i3,a1,1p2e13.4,e12.3,e11.3)
  360 format (//,' Delayed neutrons/sec-cc from',1pe12.5,' at/cc ',
     1a2,i3,a1)
  370 format (1h+,54x,9(1h_),/,53x,1pe11.3)
  380 format (///,' Total delayed neutron spectrum')
  390 format (////,' Total Neutron Spectrum')

 1900 format('SOURCES 4A Calculation',/,'<<<<<<<<<<<>>>>>>>>>>>',//)
 2000 format(/,22x,'Summary of Input')
 2010 format(22x,'================')
 2020 format(' ')
 2021 format(/,80(1h-))
 2030 format('Title: ', 11a7)
 2040 format('Homogeneous problem input (idd=',i2,')')
 2043 format('Interface problem input (idd=',i2,')')
 2047 format('Beam problem input (idd=',i2,')')
 2048 format('Magnitudes and spectra computed.')
 2049 format('Magnitudes only computed.')
 2050 format('Number of elemental constituents:',i3)
 2060 format('Solid stopping cross-sections used (isg=',i2,')')
 2070 format('Gas stopping cross-sections used (isg=',i2,')')
 2080 format('Elemental Constituents:')
 2090 format(' ')
 2100 format('  Z-value  Atom Fraction')
 2110 format('  -------  -------------')
 2120 format(4x,i3,5x,f12.10)
 2130 format('Number of neutron spectrum energy groups:', i4)
 2140 format('Maximum neutron energy is ',1pE10.3,' MeV.')
 2150 format('Minimum neutron energy is ',1pE10.3,' MeV.')
 2160 format('Energy Group Structure:')
 2170 format('  Group  Upper-Bound  Lower-Bound')
 2180 format('  -----  -----------  -----------')
 2190 format(2x,i4,4x,1pE10.3,3x,1pE10.3)
 2193 format(/,'Warning: First Group Upper-Bound Does Not Coincide ',
     1'With Maximum Neutron Energy.')
 2195	format(/,'Warning: Energy Group Structure Not Entered in ',
     1'Decreasing Order.',/,5x,'Group ',i3,' Upper-Bound is Greater ',
     2'Than Group ',i3,' Upper-Bound.') 
 2200 format('Number of source nuclides to be evaluated:',i3)
 2210 format('Alpha beam energy is ',1pE10.3,' Mev.')
 2220 format('Source Nuclides:')
 2230 format('    ZAID     Atom Density (g/cc)')
 2240 format('    ----     -------------------')
 2250 format(3x,i6,4x,1pE10.3)
 2260 format('Number of target nuclides to be used:',i3)
 2270 format(i5,' Alpha energy groups used.')
 2280 format('Target Nuclides:')
 2290 format('    ZAID     Atom Fraction')
 2300 format('    ----     -------------')
 2310 format(3x,i6,4x,1pE10.3)

 3000 format(//,22x,'Summary of Output')
 3010 format(22x,'=================')
 3020 format('Total (alpha,n) neutron source from all sources and ',
     1'targets: ',1pE10.3,' n/sec-cm^3.')
 3021 format('Total (alpha,n) neutron source from all sources and ',
     1'targets: ',1pE10.3,' n/sec-microamp.')
 3022 format(/,'Total spontaneous fission neutron source from all ',
     1'sources and targets: ',1pE10.3,' n/sec-cm^3.')
 3024 format(/,'Total delayed neutron source from all sources and ',
     1'targets: ',1pE10.3,' n/sec-cm^3.')
 3026 format(/,'Total neutron source from all sources and ',
     1'targets: ',1pE10.3,' n/sec-cm^3.')
 3030 format(/,'Average (alpha,n) neutron energy: ',1pE10.3,' MeV.')
 3032 format(/,'Average spontaneous fission neutron energy: ',
     11pE12.5,' MeV.')
 3034 format(/,'Average delayed neutron energy: ',1pE10.3,' MeV.')
 3036 format(/,'Average neutron energy from all sources: ',1pE10.3,
     1' MeV.')
 3050 format(/,'Normalized Neutron Energy Spectrum by Energy Group ',
     1'for All Sources:',//,10x,'Group',5x,'Contribution',/,10x,
     2'-----',5x,'------------')
 3051 format(11x,i3,7x,1pE10.3)
 3500 format(/,'Portion of Total Neutron Source Rate Accounted for ',
     1'in the Total Energy Spectrum: ',F5.1,'%.')
 3501 format(/,'Portion of (a,n) Neutron Source Rate Accounted for ',
     1'in the (a,n) Energy Spectrum: ',F5.1,'%.')
 3502 format(/,'Portion of Spontaneous Fission Neutron Source Rate ',
     1'Accounted For in the Spontaneous Fission Energy Spectrum: ',
     2F5.1,'%.')
 3503 format(/,'Portion of Delayed Neutron Neutron Source Rate ',
     1'Accounted for in the Delayed Neutron Energy Spectrum: ',
     2F5.1,'%.')
 3510 format(/,'Warning: Energy Group Structure May Be Poorly Chosen.')
c***********************************************************************

c-------------
c  Parameters
c-------------
      rewind 1
      aneut=1.
      alph=4.
      zalp=2.
      pi=3.14159
      adef=0.82
      bdef=4.6
  400 totqsf=0.
      isfnq=0
      idnnq=0
      gttqan=0.

c-----------------------------------
c  Read Input Parameters from Tape1
c-----------------------------------
      read (1,17,end=1440) title
      read(1,*)idd,id
      write(6,1900)
      write(7,1900)
      write(8,1900)
      write(10,1900)
      write(11,1900)
      write(6,96)
      write(6,2030)title
      write(11,2000)
      write(11,2010)
      write(11,2020)
      write(11,2030)title
      if (idd.eq.1) then
       write(11,2040)idd
      elseif (idd.eq.3) then
       write(11,2047)idd
      else
       stop 'Error in idd input!'
      endif
      if (id.eq.1) then
       write(11,2049)
      elseif (id.eq.2) then
       write(11,2048)
      else
       stop 'Error in id input!'
      endif
  410 if(id.eq.2) then
       write(7,97)
       write(7,2030)title
       write(8,98)
       write(8,2030) title
       write(10,99)
       write(10,2030)title
      endif

      read (1,*) nz,isg
      write(11,2050)nz
      if (isg.eq.0) then
       write(11,2060)isg
      elseif (isg.eq.1) then
       write(11,2070)isg
      else
       stop 'Error in isg input!'
      endif
      if (nz.le.0) go to 525

c---------------------------------
c read input material constituents
c---------------------------------
      write(11,2090)
      write(11,2080)
      write(11,2090)
      write(11,2100)
      write(11,2110)
      do 420 j=1,nz
       read (1,*) jzm(j),azm(j)
       write(11,2120) jzm(j),azm(j)
  420 continue
      write(11,2090)

c----------------------------------------
c order the nz material constituents by z
c----------------------------------------
      do 440 j=1,nz
       isw=0
       do 430 jn=2,nz
	if (jzm(jn).gt.jzm(jn-1)) go to 430
	jtem=jzm(jn-1)
	atem=azm(jn-1)
	jzm(jn-1)=jzm(jn)
	azm(jn-1)=azm(jn)
	jzm(jn)=jtem
	azm(jn)=atem
	isw=1
  430  continue
       if (isw.eq.0) go to 450
  440 continue

c----------------------------------------------------
c read tape 2 for stopping cross-section coefficients
c----------------------------------------------------
  450 rewind 2
  460 read (2,10) icont
      if (icont.ne.0) go to 460
      do 520 j=1,nz
  470  read (2,40,end=510) iz,c
  480  if (jzm(j)-iz) 510,490,470
  490  continue
       do 500 jt=1,9
  500  czm(jt,j)=c(jt)
       if (isg.eq.1.and.c(10).gt.0.) then
	do 505 jt=1,5
  505   czm(jt,j)=c(jt+9)
       endif
       go to 520
  510  write (6,50) jzm(j)
       stop 'Element not found in tape2'
  520 continue
  525 if (id.eq.1) go to 580

c------------------------------------------------------------
c construct neutron group structure desired from tape 1 input
c neutron groups ordered in decreasing energy
c------------------------------------------------------------
      read (1,*) nng,enmax,enmin
      if (nng.gt.0) then
       write(11,2130)nng
      else
       idummy=-nng
       write(11,2130)idummy
      endif
      write(11,2140)enmax
      write(11,2150)enmin
      write(11,2090)
      write(11,2160)
      write(11,2090)
      write(11,2170)
      write(11,2180)
      nngp1=nng+1
      if (nng.gt.0) go to 540
      nng=-nng
      nngp1=nng+1
      en(nngp1)=enmin
      read (1,*) (en(n),n=1,nng)
      do 527 n=1,nng
	write(11,2190) n,en(n),en(n+1)
  527 continue
      if (en(1).ne.enmax) write(11,2193)
	do 528 n=2,nng
	 nm1=n-1
	 if (en(n).gt.en(n-1)) then
	  write(11,2195)n,nm1
	  stop 'Energy Structure Incorrectly Entered!'
	 endif
  528 continue
      if (en(1).gt.en(2)) go to 560
      nh=nngp1/2
      do 530 n=1,nh
      etem=en(nngp1-n+1)
      en(nngp1-n+1)=en(n)
  530 en(n)=etem
      go to 560
  540 fnng=nng
      den=(enmax-enmin)/fnng
      en(nngp1)=enmin
      en(1)=enmax
      do 550 n=2,nng
       nnm1=n-1
       fnm1=n-1
       en(n)=en(1)-fnm1*den
       write(11,2190)nnm1,en(nnm1),en(n)
  550 continue
      write(11,2190)nnm1,en(nng),en(nngp1)
  560 write(7,70)
      write(7,80) (en(n),n=1,nngp1)
      write(7,2021)
      write(8,70)
      write(8,80) (en(n),n=1,nngp1)
      write(8,2021)
      write(10,70)
      write(10,80) (en(n),n=1,nngp1)
      write(10,2021)

c----------------------------
c zero total spectrum storage
c----------------------------
      do 570 n=1,nng
       ts(n)=0.
       gtsan(n)=0.
  570 tssf(n)=0.
      etpall=0.
      ebtall=0.
      etopan=0.
      ebotan=0.
      etopsf=0.
      ebotsf=0.
      etopdn=0.
      ebotdn=0.
      gttqan=0.
      totqsf=0.
      totqdn=0.
      qtotal=0.
      ebarsf=0.
      ebardn=0.
      ebrall=0.

c---------------------------------------------
c read alpha/s.f./d.n. sources from user input
c---------------------------------------------
  580 if (idd.eq.1) then
       read(1,*)nq
      elseif (idd.eq.3) then
       read(1,*)ebeam
       nq=0
      else
       stop 'Error in idd input!'
      endif
      write(11,2090)
      if (idd.eq.1) then
       write(11,2200)nq
      else
       write(11,2210)ebeam
       go to 620
      endif
      if (nq.eq.0) then
       stop 'Error in nq input!'
      endif
      write(11,2090)
      write(11,2220)
      write(11,2090)
      write(11,2230)
      write(11,2240)
      do 590 k=1,nq
       read (1,*) jq(k),aq(k)
       write(11,2250) jq(k),aq(k)
  590 continue

c-----------------------------------------------
c order sources by z-a-state id = s+10*a+10000*z
c-----------------------------------------------
      if (nq.eq.1) go to 620
      do 610 k=1,nq
       isw=0
       do 600 kn=2,nq
	if (jq(kn).gt.jq(kn-1)) go to 600
	jtem=jq(kn-1)
	atem=aq(kn-1)
	jq(kn-1)=jq(kn)
	aq(kn-1)=aq(kn)
	jq(kn)=jtem
	aq(kn)=atem
	isw=1
  600  continue
       if (isw.eq.0) go to 620
  610 continue

c----------------------------------------
c read target information from user input
c----------------------------------------
  620 if (nz.eq.0) go to 1140
      read (1,*) nt,nag
      write(11,2090)
      write(11,2260)nt
      write(11,2270)nag
      write(11,2090)
      write(11,2280)
      write(11,2090)
      write(11,2290)
      write(11,2300)
      if (nt.eq.0) go to 1140

c----------------------------------------------------
c beginning of big loop over target nuclides
c         for nt=0 calculate spontaneous fission only
c----------------------------------------------------
      rewind 3
      if (id.eq.2) rewind 4
      if (idd.eq.3) then
	 ajwd1='/'
	 ajwd2='micr'
	 ajwd3='oamp'
	else
	 ajwd1='/'
	 ajwd2='cm**'
	 ajwd3='3   '
      endif
      write (6,100) ajwd1,ajwd2,ajwd3,ajwd1,ajwd2,ajwd3

c--------------------------
c loop on target nuclides i
c--------------------------
  630 read (3,10) icont
      if (icont.ne.0) go to 630
  640 if (id.eq.2) read (4,10) icont
      if (icont.ne.0) go to 640
      do 1120 i=1,nt
  644 format('Target nuclides complete: ', i3, ' out of ', i3) !pk
  645 write(0,644) i,nt !pk
c The last line counts iterations of the big loop so we know
c the program didn't freeze.
      etpant=0.
      ebtant=0.
      totqan=0.
      if(id.eq.1) go to 660
      do 650 n=1,nng
  650 tsan(n)=0.

c------------------------------------------------------
c find alpha-n cross sections for this target on tape 3
c------------------------------------------------------
  660 read (1,*) idt,at
      write(11,2310)idt,at
  665 read (3,105,end=680) jdt,jps
c      print *,jdt,jps 
  670 read (3,120) (e(ip),x(ip),ip=1,jps)
c      print *,(e(ip),x(ip),ip=1,jps)
      if (idt-jdt) 680,690,665
  680 write (6,130) idt
      stop 'Target nuclide not found in tape3'
  690 lzt=idt/10000
      lat=idt/10-lzt*1000
      atar=lat
      lst=idt-lzt*10000-lat*10
      amt=' '
      if (lst.ne.0) amt='m'
      apro=atar+3.
      do 700 ip=2,jps
       eamin=e(ip-1)
       if (x(ip).ne.0.) go to 710
  700 continue
  710 if (eamin.lt.0.001) eamin=0.001
      if (id.eq.1) go to 760

c----------------------------------------------------------------
c find product nuclide level branchings for this target of tape 4
c----------------------------------------------------------------
  720 read (4,140,end=750) jdt,jl,jp
  730 read (4,120) q,(el(il),il=1,jl)
      do 740 ip=1,jp
  740 read (4,120) ep(ip),(f(il,ip),il=1,jl)
      if (idt-jdt) 750,760,720
  750 write (6,150) idt
      stop 'Target nuclide not found in tape4'
  760 rewind 5
      if (idd.eq.3) then 
       sbtqan=0.
       if(id.ne.1) then
	do 765 il=1,jl
  765   totlev(il)=0.
       endif
	 nq=1
	 go to 771
      endif
  770 read (5,10) icont
      if (icont.ne.0) go to 770
  771 continue
c--------------------------------------------
c loop on source nuclides k for this target i
c contribution to (alpha,n) neutrons
c--------------------------------------------
      do 1100 k=1,nq
	if (idd.eq.3) goto 850
      etop=0.
      ebot=0.
      do 775 il=1,jl
       etopl(il)=0.
       ebotl(il)=0.
       do 774 n=1,nng
  774  sl(n,il)=0.
  775 totlev(il)=0.
      sbtqan=0.
      if(id.eq.1) go to 790
      do 780 n=1,nng
  780 san(n)=0.
      if(nq.le.0) go to 850

c-------------------------------------------
c read source nuclide decay data from tape 5
c-------------------------------------------
  790 read (5,110,end=830) idq,nal,ndn
  800 read (5,60) alam,bfsf,barnu,a,b,bfdn
      if (nal.eq.0) go to 810
      read (5,120) (eal(l),fal(l),l=1,nal)
  810 if (ndn.eq.0) go to 820
      read (5,120) (fdng(nd),nd=1,ndn)
  820 if (idq-jq(k)) 790,840,830
  830 write (6,160) jq(k)
      stop 'Source nuclide not found in tape5'
  840 if (i.eq.1.and.bfsf.gt.0.) isfnq=isfnq+1
      if (i.eq.1.and.bfdn.gt.0.) idnnq=idnnq+1
      if (nal.eq.0) go to 1100
      if(eal(nal).gt.e(1)) go to 845
      write(6,842) idq,idt
  842 format(i8,13h alphas below,i8,20h (alpha,n) threshold)
      go to 1100
  845 lzq=idq/10000
      laq=idq/10-lzq*1000
      lsq=idq-lzq*10000-laq*10
      amq=' '
      if (lsq.ne.0) amq='m'

c---------------------------------------------------------
c set alpha energy grid ee(m) for this target/source combo
c---------------------------------------------------------
      eamax=eal(nal)
  850 continue
      if (idd.eq.3) eamax=ebeam
      if(eamax.le.eamin) then
       qan=0.
       mm=0
       p(nagp1)=0.
       pval=0.
       go to 945
      endif
      if (eamax.le.9.8) go to 860                      !modified - vk
      write (6,170) eamax
      stop 'Maximum alpha energy above 9.8 MeV'
  860 fnag=nag
      nagp1=nag+1
      dea=(eamax-eamin)/fnag
      ee(1)=eamin
      do 870 m=2,nag
       fm=m-1
  870 ee(m)=ee(1)+dea*fm
      ee(nagp1)=eamax

c----------------------------------------------------
c for each energy ee(m) calc. and store cross section
c----------------------------------------------------
      ip=1
      do 900 m=1,nagp1
  880  if (ee(m).ge.e(ip).and.ee(m).le.e(ip+1)) go to 890
       ip=ip+1
       if (ip.lt.jps) go to 880
       write (6,180) idt,m,ee(m)
       stop 'Error in energy structure!' 
  890  slope=(x(ip+1)-x(ip))/(e(ip+1)-e(ip))
c If the denominator is zero we'll get an Infinity or a NaN !pk
c I'm not fixing it but I want to be able to know when it happens. !pk
  894 format('Warning: Division by zero in the ',i3,'th target a',
     1'nd',i3,'th source around the',i5,'th alpha energy group.') !pk
       denomi = e(ip+1)-e(ip) !pk
       if (denomi.eq.0.0) then !pk
        write(0,894)i,k,m !pk
       endif !pk
       enrsep=(x(ip)*e(ip+1)-x(ip+1)*e(ip))/(e(ip+1)-e(ip))
  900 cx(m)=enrsep+ee(m)*slope

c-------------------------------------------------------------
c for each energy ee(m) calc. and store stopping cross section
c-------------------------------------------------------------
      do 930 m=1,nagp1
       scx(m)=0.
       do 920 j=1,nz
	zmat=jzm(j)
	amat=amass(jzm(j))
c-----  nuclear stopping---
	term=zalp**0.66667+zmat**0.66667
	sterm=sqrt(term)
	bot=zalp*zmat*(alph+amat)*sterm
	rep=32530.*amat*ee(m)/bot
	if(rep.lt.0.001) then
	 dcx=1.593*sqrt(rep)
	 go to 905
	endif
	if(rep.lt.10.) then
	 bot=1.+6.8*rep+3.4*rep**1.5
	 dcx=1.7*sqrt(rep)*alog(rep+2.71828)/bot
	 go to 905
	endif
	dcx=alog(0.47*rep)/(2.*rep)
c-----  electronic stopping---
  905   if(ee(m).gt.30.) go to 906
	slow=czm(1,j)*(1000.*ee(m))**czm(2,j)
	shigh=(czm(3,j)/ee(m))*alog(1.+czm(4,j)/ee(m)+czm(5,j)*ee(m))
	dcx=dcx+slow*shigh/(slow+shigh)
	go to 907
  906   eil=alog(1/ee(m))
	arg=czm(6,j)
	arg=arg+eil*(czm(7,j)+czm(8,j)*eil+czm(9,j)*eil*eil)
	dcx=dcx+exp(arg)
  907   continue
  920  scx(m)=scx(m)+azm(j)*dcx
  930 continue

c-----------------------------------------------------------
c for each energy ee(m) calc. ratio r(m), then integral p(m)
c-----------------------------------------------------------
      r(1)=cx(1)/scx(1)
      p(1)=0.
      fact=1.0e-06*at
      do 940 m=2,nagp1
       r(m)=cx(m)/scx(m)
  940 p(m)=p(m-1)+fact*(r(m-1)+r(m))*dea/2.
  945 if (idd.ne.3) go to 950

c---------------------------------------------------------------------
c if nq=0 calculate neutrons from alpha beam: 1 microamp=3.1209+12 a/s
c---------------------------------------------------------------------
      aps=3.1209e+12
      beamn=aps*p(nagp1)
      pval=p(nagp1)
      write (6,190) jsm(lzt),lat,amt,at,ebeam,aps,pval,beamn
      nal=1
      eal(1)=ebeam
      qan=beamn
      sbtqan=sbtqan+qan
      totqan=totqan+qan
      if (id.eq.1) go to 1101
      mm=nag

c----------------------------------------------
c if nq>0 calculate neutrons from source alphas
c----------------------------------------------
  950 do 1080 l=1,nal
      if (idd.eq.3) go to 1000
      do 960 m=1,nag
       mm=m
  960 if (eal(l).ge.ee(m).and.eal(l).le.ee(m+1)) go to 970
      mm=0
      pval=0.
      go to 975
  970 pval=p(mm)+(eal(l)-ee(mm))*(p(mm+1)-p(mm))/(ee(mm+1)-ee(mm))
  975 aps=aq(k)*alam*fal(l)
      qan=aps*pval
      sbtqan=sbtqan+qan
      totqan=totqan+qan
      if (l.eq.1) go to 980
      write (6,200) eal(l),aps,pval,qan
      if(nal.ne.1.and.l.eq.nal) write(6,250) sbtqan
      go to 990
  980 write(6,210)jsm(lzt),lat,amt,at,jsm(lzq),laq,amq,eal(l),aps,
     1pval,qan
      if(nal.eq.1) write(6,255)
  990 if (id.eq.1.or.mm.eq.0) go to 1080

c-----------------------------------------------------------------
c calculate (alpha,n) neutron spectrum contribution in multigroups
c-----------------------------------------------------------------
 1000 mmm1=mm-1
      do 1010 m=1,mmm1
 1010 rr(m)=(p(m+1)-p(m))/pval
      rr(mm)=(pval-p(mm))/pval
      do 1020 n=1,nng
 1020 s(n)=0.
      do 1060 il=1,jl
       qlev=q-el(il)
       e90=-qlev*apro/(apro-alph)
       thre=-qlev*(aneut+apro)/(aneut+apro-alph)
       if (qlev.gt.0.) thre=0.
       do 1055 m=1,mm
	ea=(ee(m)+ee(m+1))/2.
	if (m.eq.mm) ea=(eal(l)+ee(mm))/2.
	if (ea.lt.thre) go to 1055
	do 1030 ip=2,jp
	 lp=ip
 1030   if (ep(ip-1).le.ea.and.ep(ip).ge.ea) go to 1040
	write (6,220)
	stop 'Error in alpha energy structure'
 1040   lp1=lp-1
	bx=f(il,lp1)
	bx=bx+(ea-ep(lp1))*(f(il,lp)-f(il,lp1))/(ep(lp)-ep(lp1))
	term1=sqrt(4.*alph*aneut*ea)/(2.*(aneut+apro))
	term2=alph*aneut*ea/((aneut+apro)*(aneut+apro))
	term3=(apro*ea+apro*qlev-alph*ea)/(aneut+apro)
	senmax=term1+sqrt(term2+term3)
	enmax=senmax*senmax
	if (ea.le.e90) senmin=term1-sqrt(term2+term3)
	if (ea.gt.e90) senmin=-term1+sqrt(term2+term3)
	enmin=senmin*senmin
	ebar=(enmin+enmax)/2.
	val=qan*rr(m)*bx
	etop=etop+val*ebar
	etopl(il)=etopl(il)+val*ebar
	ebotl(il)=ebotl(il)+val
	ebot=ebot+val
	dele=enmax-enmin
	do 1050 n=1,nng
	 if (en(n).lt.enmin) go to 1050
	 if (en(n+1).gt.enmax) go to 1050
	 de=en(n)-en(n+1)
	 if (en(n+1).lt.enmin) de=de-(enmin-en(n+1))
	 if (en(n).gt.enmax) de=de-(en(n)-enmax)
	 gpadd=qan*rr(m)*bx*de/dele
	 s(n)=s(n)+gpadd
	 sl(n,il)=sl(n,il)+gpadd
	 totlev(il)=totlev(il)+gpadd
 1050   continue
 1055  continue
 1060 continue
      do 1070 n=1,nng
       san(n)=san(n)+s(n)
       gtsan(n)=gtsan(n)+s(n)
       ts(n)=ts(n)+s(n)
 1070 tsan(n)=tsan(n)+s(n)
 1075 continue
 1080 continue
      if(id.eq.1) go to 1100

c----------------------------------------------------
c output this (alpha,n) neutron spectrum contribution
c for this source/target combination
c----------------------------------------------------
      ebarqt=0.
      if(ebot.le.0.) go to 1085
      ebarqt=etop/ebot
      etpant=etpant+etop
      ebtant=ebtant+ebot
 1085 continue
      smga=0.
      do 1090 n=1,nng
 1090 smga=smga+san(n)
c     fracgp=smga/sbtqan
c     do 1095 il=1,jl
c1095 frclev(il)=totlev(il)/sbtqan
c     write(6,205) (frclev(il),il=1,jl)
      if (idd.eq.3) write(7,230)ebeam,jsm(lzt),lat,amt
      if (idd.ne.3) write(7,240)aq(k),jsm(lzq),laq,amq,jsm(lzt),lat,amt
      write(7,80) (san(n),n=1,nng)
      write(7,81) smga,ajwd1,ajwd2,ajwd3,ebarqt
      smga=0.
      if (idd.eq.3) write(8,230)ebeam,jsm(lzt),lat,amt
      if (idd.ne.3) write(8,240)aq(k),jsm(lzq),laq,amq,jsm(lzt),lat,amt
      do 1097 n=1,nng
       san(n)=san(n)/sbtqan
 1097 smga=smga+san(n)
      write(8,80) (san(n),n=1,nng)
      write(8,81) smga,ajwd1,ajwd2,ajwd3,ebarqt
      do 1098 il=1,jl
       if(ebotl(il).le.0.) go to 1098
       ebarl=etopl(il)/ebotl(il)
       if (idd.eq.3) then
	write(10,262) ebeam,jsm(lzt),lat,amt,el(il)
       else
	write(10,261) jsm(lzq),laq,amq,jsm(lzt),lat,amt,el(il)
       endif
       write(10,80) (sl(n,il),n=1,nng)
       write(10,81) totlev(il),ajwd1,ajwd2,ajwd3,ebarl
 1098 continue
      if(idd.eq.3) go to 1101
 1100 continue
 1101 gttqan=gttqan+totqan

c--------------------------------------
c output total (alpha,n) source for all
c alpha sources on this target nuclide
c--------------------------------------
      if (nq.gt.1) write (6,251) totqan
      if(id.eq.1) go to 1120
      tmga=0.
      do 1110 n=1,nng
 1110 tmga=tmga+tsan(n)
c     fracgp=tmga/totqan
      ebart=etpant/ebtant
      etopan=etopan+etpant
      ebotan=ebotan+ebtant
      if(nq.le.1) go to 1120
      write(7,260)
      write(7,80) (tsan(n),n=1,nng)
      write(7,81) tmga,ajwd1,ajwd2,ajwd3,ebart
      write(7,2090)
      write(7,2021)
      write(7,2090)
      tmga=0.
      do 1115 n=1,nng
       tsan(n)=tsan(n)/totqan
 1115 tmga=tmga+tsan(n)
      write(8,260)
      write(8,80) (tsan(n),n=1,nng)
      write(8,81) tmga,ajwd1,ajwd2,ajwd3,ebart
      write(8,2090)
      write(8,2021)
      write(8,2090)
      write(10,2090)
      write(10,2021)
      write(10,2090)
 1120 continue

c----------------------------------------------------------------------
c output grand total (alpha,n) neutron source: all targets, all sources
c----------------------------------------------------------------------
      if (nt.gt.1) write (6,252) gttqan
      if(id.eq.1) go to 1140
      gtmga=0.
      do 1130 n=1,nng
 1130 gtmga=gtmga+gtsan(n)
      ebaran=etopan/ebotan
c     fracgp=gtmga/gttqan
      write(7,2021)
      write(7,270)
      write(7,80) (gtsan(n),n=1,nng)
      write(7,81) gtmga,ajwd1,ajwd2,ajwd3,ebaran
      write(7,2090)
      write(7,2021)
      write(7,2021)
      write(7,2090)
      gtmga=0.
      do 1135 n=1,nng
       dummy1=dummy1+gtsan(n)
       gtsan(n)=gtsan(n)/gttqan
 1135 gtmga=gtmga+gtsan(n)
      write(8,2021)
      write(8,270)
      write(8,80) (gtsan(n),n=1,nng)
      write(8,81) gtmga,ajwd1,ajwd2,ajwd3,ebaran
      write(8,2090)
      write(8,2021)
      write(8,2021)
      write(8,2090)
      write(8,2090)

c-----------------------------
c calculate s.f.neutron source
c-----------------------------
 1140 if (idd.eq.3) go to 1400
      if (isfnq.eq.0.and.nt.gt.0) go to 1260
      rewind 5
 1150 read (5,10) icont
      if (icont.ne.0) go to 1150
      write (6,280)
      tmgs=0.
      etopsf=0.
      ebotsf=0.
      do 1250 k=1,nq
 1160 read (5,110,end=1200) idq,nal,ndn
 1170 read (5,60) alam,bfsf,barnu,a,b,bfdn
      if (nal.eq.0) go to 1180
      read (5,120) (eal(l),fal(l),l=1,nal)
 1180 if (ndn.eq.0) go to 1190
      read (5,120) (fdng(nd),nd=1,ndn)
 1190 if (idq-jq(k)) 1160,1210,1200
 1200 write (6,160) jq(k)
      stop 'S.F. source nuclide not found on tape5'
 1210 qsf=aq(k)*alam*bfsf*barnu
      if(qsf.le.0.) go to 1250
      ebar=0.25*a*a*b + 1.5*a
      etopsf=etopsf+ebar*qsf
      ebotsf=ebotsf+qsf
      totqsf=totqsf+qsf
      lzq=idq/10000
      laq=idq/10-lzq*1000
      lsq=idq-lzq*10000-laq*10
      amq=' '
      if (lsq.ne.0) amq='m'

c---------------------------------------------
c output this s.f. neutron source contribution
c---------------------------------------------
      write (6,290) jsm(lzq),laq,amq,aq(k),alam,bfsf,barnu,qsf
      if (id.eq.1) go to 1250
      do 1220 n=1,nng
 1220 ssf(n)=0.

c----------------------------------------------------------------
c calc. frac. watt spectrum contained in neutron energy structure
c----------------------------------------------------------------
      if (a.gt.0.0.and.b.gt.0.0) go to 1230
      a=adef
      b=bdef
      write (7,300) jsm(lzq),laq,amq,a,b
 1230 se1=sqrt(en(1))
      sen1=sqrt(en(nng+1))
      sa=sqrt(a)
      sbaso4=sqrt(b*a*a/4.)
      c1=(sen1-sbaso4)/sa
      c2=(se1-sbaso4)/sa
      c3=(sen1+sbaso4)/sa
      c4=(se1+sbaso4)/sa
      if (c1.lt.0..and.c2.lt.0.) eft=0.5*(ERF1(-c1)-ERF1(-c2))
      if (c1.lt.0..and.c2.ge.0.) eft=0.5*(ERF1(-c1)+ERF1(c2))
      if (c1.ge.0..and.c2.ge.0.) eft=0.5*(ERF1(c2)-ERF1(c1))
      eft=eft+0.5*(ERF1(c4)-ERF1(c3))
      c1s=c1*c1
      c2s=c2*c2
      c3s=c3*c3
      c4s=c4*c4
      spba=sqrt(pi*b*a)
c     fwatt=eft+(1./spba)*(exp(-c1s)-exp(-c2s)-exp(-c3s)+exp(-c4s))

c-------------------------------------------
c calculate multigroup s.f. neutron spectrum
c-------------------------------------------
      w=qsf
      smgs=0.
      do 1240 n=1,nng
       seh=sqrt(en(n))
       sel=sqrt(en(n+1))
       c1=(sel-sbaso4)/sa
       c2=(seh-sbaso4)/sa
       c3=(sel+sbaso4)/sa
       c4=(seh+sbaso4)/sa
       c1s=c1*c1
       c2s=c2*c2
       c3s=c3*c3
       c4s=c4*c4
       if (c1.lt.0..and.c2.lt.0.) eft=0.5*(ERF1(-c1)-ERF1(-c2))
       if (c1.lt.0..and.c2.ge.0.) eft=0.5*(ERF1(-c1)+ERF1(c2))
       if (c1.ge.0..and.c2.ge.0.) eft=0.5*(ERF1(c2)-ERF1(c1))
       eft=eft+0.5*(ERF1(c4)-ERF1(c3))
       ssf(n)=(1./spba)*(exp(-c1s)-exp(-c2s)-exp(-c3s)+exp(-c4s))
       ssf(n)=w*(eft+ssf(n))
       ts(n)=ts(n)+ssf(n)
       smgs=smgs+ssf(n)
 1240 tssf(n)=tssf(n)+ssf(n)
      tmgs=tmgs+smgs

c-----------------------------------------------
c output this s.f. neutron spectrum contribution
c-----------------------------------------------
      write (7,310) aq(k),jsm(lzq),laq,amq
      write (7,80) (ssf(n),n=1,nng)
      write (7,81) smgs,ajwd1,ajwd2,ajwd3,ebar
      smgs=0.
      do 1245 n=1,nng
       ssf(n)=ssf(n)/qsf
 1245 smgs=smgs+ssf(n)
      write (8,310) aq(k),jsm(lzq),laq,amq
      write (8,80) (ssf(n),n=1,nng)
      write (8,81) smgs,ajwd1,ajwd2,ajwd3,ebar
 1250 continue

c----------------------------------------------
c output total s.f. neutron source and spectrum
c----------------------------------------------
      if(isfnq.gt.1) write(6,320) totqsf
      if(id.eq.1.or.isfnq.lt.2) go to 1260
      ebarsf=etopsf/ebotsf
c     fracgp=tmgs/totqsf
      write(7,2021)
      write(7,330)
      write(7,80) (tssf(n),n=1,nng)
      write(7,81) tmgs,ajwd1,ajwd2,ajwd3,ebarsf
      write(7,2090)
      write(7,2021)
      write(7,2090)
      tmgs=0.
      do 1255 n=1,nng
	dummy2=dummy2+tssf(n)
        tssf(n)=tssf(n)/totqsf
 1255 tmgs=tmgs+tssf(n)
      write(8,2021)
      write(8,330)
      write(8,80) (tssf(n),n=1,nng)
      write(8,81) tmgs,ajwd1,ajwd2,ajwd3,ebarsf
      write(8,2090)
      write(8,2021)
      write(8,2090)

c----------------------------------
c calculate delayed neutron sources
c----------------------------------
 1260 if (idnnq.eq.0) go to 1400
      rewind 5
      totqdn=0.
      if (id.eq.1) go to 1280
      do 1270 n=1,nng
 1270 tsdn(n)=0.
 1280 read (5,10) icont
      if (icont.ne.0) go to 1280
      write (6,340)
      tmgs=0.
      etopdn=0.
      ebotdn=0.
      do 1390 k=1,nq
 1290  read (5,110,end=1330) idq,nal,ndn
 1300 read (5,60) alam,bfsf,barnu,a,b,bfdn
      if (nal.eq.0) go to 1310
      read (5,120) (eal(l),fal(l),l=1,nal)
 1310 if (ndn.eq.0) go to 1320
      read (5,120) (fdng(nd),nd=1,ndn)
 1320 if (idq-jq(k)) 1290,1340,1330
 1330 write (6,160)
      stop 'D.N. source nuclide not found on tape5'
 1340 qdn=aq(k)*alam*bfdn
      etpdnq=0.
      ebtdnq=0.
      if(qdn.le.0.) go to 1390
      totqdn=totqdn+qdn
      lzq=idq/10000
      laq=idq/10-lzq*1000
      lsq=idq-laq*10-lzq*10000
      amq=' '
      if (lsq.eq.1) amq='m'
      if (lsq.eq.2) amq='n'

c-----------------------------------------------
c output the delayed neutron source contribution
c-----------------------------------------------
      write (6,350) jsm(lzq),laq,amq,aq(k),alam,bfdn,qdn
      if (id.eq.1) go to 1390

c---------------------------------------------------
c calculate fraction of this d.n. spectrum contained
c in the neutron energy structure.
c---------------------------------------------------
      fdn=1.
      fndn=ndn
      top=fndn*.01
      if (en(nngp1).eq.0..and.en(1).ge.top) go to 1355
      fdn=0.
      do 1350 nd=1,ndn
       fnd=nd
       top=fnd*.01
       bot=top-.01
       ebar=(top+bot)/2.
       val=qdn*fdng(nd)
       etpdnq=etpdnq+val*ebar
       ebtdnq=ebtdnq+val
       if (top.le.en(nngp1)) go to 1350
       if (bot.ge.en(1)) go to 1350
       frac=1.
       botlim=bot
       if (bot.lt.en(nngp1)) botlim=en(nngp1)
       toplim=top
       if (top.gt.en(1))  toplim=en(1)
       frac=(toplim-botlim)/.01
       fdn=fdn+frac*fdng(nd)
 1350 continue

c----------------------------------------------
c calculate multigroup delayed neutron spectrum
c----------------------------------------------
 1355 smgs=0.
      do 1380 n=1,nng
       sdn(n)=0.
       do 1370 nd=1,ndn
        fnd=nd
        top=fnd*.01
        bot=top-.01
        if(bot.ge.en(n)) go to 1375
        if(top.lt.en(n+1)) go to 1370
        frac=1.
        if(bot.ge.en(n+1).and.top.le.en(n)) go to 1360
        if(top.gt.en(n)) top=en(n)
        if(bot.lt.en(n+1)) bot=en(n+1)
        frac=(top-bot)/.01
 1360   sdn(n)=sdn(n)+frac*fdng(nd)*qdn
 1370  continue
 1375  smgs=smgs+sdn(n)
       tsdn(n)=tsdn(n)+sdn(n)
       ts(n)=ts(n)+sdn(n)
 1380 continue
      tmgs=tmgs+smgs

c--------------------------------------------------
c output this delayed neutron spectrum contribution
c--------------------------------------------------
      etopdn=etopdn+etpdnq
      ebotdn=ebotdn+ebtdnq
      ebrdnq=etpdnq/ebtdnq
      write (7,360) aq(k),jsm(lzq),laq,amq
      write (7,80) (sdn(n),n=1,nng)
      write (7,81) smgs,ajwd1,ajwd2,ajwd3,ebrdnq
      smgs=0.
      do 1385 n=1,nng
       sdn(n)=sdn(n)/qdn
 1385 smgs=smgs+sdn(n)
      write (8,360) aq(k),jsm(lzq),laq,amq
      write (8,80) (sdn(n),n=1,nng)
      write (8,81) smgs,ajwd1,ajwd2,ajwd3,ebrdnq
 1390 continue

c-------------------------------------------------
c output total delayed neutron source and spectrum
c-------------------------------------------------
      if(idnnq.gt.1) write(6,370) totqdn
      if(id.eq.1.or.idnnq.lt.2) go to 1400
      ebardn=etopdn/ebotdn
c     fracgp=tmgs/totqdn
      write(7,2021)
      write(7,380)
      write(7,80) (tsdn(n),n=1,nng)
      write(7,81) tmgs,ajwd1,ajwd2,ajwd3,ebardn
      write(7,2090)
      write(7,2021)
      write(7,2090)
      tmgs=0.
      do 1395 n=1,nng
       dummy3=dummy3+tsdn(n)
       tsdn(n)=tsdn(n)/totqdn
 1395 tmgs=tmgs+tsdn(n)
      write(8,2021)
      write(8,380)
      write(8,80) (tsdn(n),n=1,nng)
      write(8,81) tmgs,ajwd1,ajwd2,ajwd3,ebardn
      write(8,2090)
      write(8,2021)
      write(8,2090)
 1400 continue
      if (id.eq.1) go to 1430
      itest=0
      if(nt.gt.0) itest=1
      if(isfnq.gt.0) itest=itest+1
      if(idnnq.gt.0) itest=itest+1
      if(itest.lt.2) then
       do 1405 n=1,nng
        dummy=dummy+ts(n)
        ts(n)=ts(n)/gttqan
 1405  gtmg=gtmg+ts(n)
       go to 1430
      endif

c------------------------------------------------------------
c output grand total (alpha,n) + s.f. + d.n. neutron spectrum
c------------------------------------------------------------
      gtq=gttqan+totqsf+totqdn
      dummy=0.0
      gtmg=0.
      do 1420 n=1,nng
 1420 gtmg=gtmg+ts(n)
c     fracgp=gtmg/gtq
      etpall=etopan+etopsf+etopdn
      ebtall=ebotan+ebotsf+ebotdn
      ebrall=etpall/ebtall
      write(7,2090)
      write(7,2021)
      write(7,390)
      write(7,80) (ts(n),n=1,nng)
      write(7,81) gtmg,ajwd1,ajwd2,ajwd3,ebrall
      gtmg=0.
      do 1425 n=1,nng
       dummy=dummy+ts(n)
       ts(n)=ts(n)/gtq
 1425 gtmg=gtmg+ts(n)
      write(8,2090)
      write(8,2021)
      write(8,390)
      write(8,80) (ts(n),n=1,nng)
      write(8,81) gtmg,ajwd1,ajwd2,ajwd3,ebrall
 1430 continue

c-----------------
c  Output Summary
c-----------------
 1440 continue
      qtotal=gttqan+totqsf+totqdn
      ebrall=ebaran*gttqan
      if (isfnq.ne.0) ebrall=ebrall+(ebarsf*totqsf)
      if (idnnq.ne.0) ebrall=ebrall+(ebardn*totqdn)
      ebrall=ebrall/qtotal
      if (nq.eq.0) ebrall=ebaran
      write(11,2090)
      write(11,2090)
      write(11,3000)
      write(11,3010)
      write(11,2090)
      if (idd.eq.1) then
       write(11,3020)gttqan
       write(11,3022)totqsf
       write(11,3024)totqdn
       write(11,3026)qtotal
      else
       write(11,3021)gttqan
      endif
      if (id.gt.1) then
       write(11,3030)ebaran
       if (idd.eq.1) then
	    write(11,3032)ebarsf
	    write(11,3034)ebardn
	    write(11,3036)ebrall
       endif
	   write(11,3050)
	   do 1450 n=1,nng
	    write(11,3051)n,ts(n)
1450   continue
       frac=100.0*(dummy/qtotal)
	 write(11,3500)frac
	 if (frac.le.90.0) write(11,3510)
	 if (itest.gt.2) then
	  frac=100.0*(dummy1/gttqan)
	  write(11,3501)frac
	 endif
	 if (isfnq.ne.0) then
	  frac=100.0*(dummy2/totqsf)
	  write(11,3502)frac
	  if (frac.le.90.0) write(11,3510)
	 endif
	 if (idnnq.ne.0) then
	  frac=100.0*(dummy3/totqdn)
	  write(11,3503)frac
	  if (frac.le.90.0) write(11,3510)
	 endif	  
      endif

c------------------------
c  End Execution of HOMO
c------------------------
      return
      end



c=======================================================================
c  Interface Problem Main Subroutine (8/96)
c=======================================================================

      subroutine interf(title,idd,id)

c----------
c  Storage
c----------
      
      character*7 title
      dimension title(11), ast(4001), eal(4001)

c------------------
c  Call Statements
c------------------

      call sainf(nag,title,idd,id,ast,eal)
      call sninf(nag,title,id,ast,eal)

c--------------------------------
c  End of Execution of INTERFACE
c--------------------------------

      return
      end


c=======================================================================
c  Interface Problem Alpha Source Subroutine (8/96)
c=======================================================================

      subroutine sainf(nag,title,idd,id,ast,eal)

c---------
c storage
c---------

      character*7 title
      character*2 jsm(105)
      dimension title(11), jzm(20), azm(20), czm(9,20), c(14), jq(300),
     1 aq(300), eala(30), eal(4001), ee(4001), scx(4001), r(4001), 
     2 p(4001), fdng(1000), ast(4001),fal(30),ea(4001),ps(4001),
     3 faq(300)
       common /names/ jsm
       common /masses/ amass(105)

c--------------------
c  Format Statements
c--------------------
   10 format (i3,11a7)
   15 format(11a7)
   16 format(11a7,//)
   20 format (2i3,f7.3)
   30 format (i8,e12.5)
   40 format(i3,8e8.3,f10.3,/,3x,5f8.3)
   50 format (69h failed to find stopping cross section coefficients on
     1tape 2 for iz=,i3,5h stop)
   60 format (6e12.5)
   90 format (i3,5e12.5)
   95 format(2i6)
   96 format(/,26x,'Neutron Source Magnitudes',/,26x,25(1h~),/)
   97 format(/,23x,'Absolute Neutron Source Spectra',/,23x,31(1h~),/)
   98 format(/,22x,'Normalized Neutron Source Spectra',/,22x,
     133(1h~),/)
   99 format(/,16x,'Neutron Source Spectra by Nuclide Energy Level',
     1/,16x,46(1h~),/)
  110 format (i8,2i4)
  120 format (8e10.3)
  130 format (1x,'SOURCE DOES NOT HAVE ALPHA EMITTER')
  131 format(' WARNING-AN ALPHA ENERGY THAT EXCEEDED ',f7.3,' MEV WAS',
     1' SET TO:',f7.3,' MeV')
  132 format(' NUMBER OF STOPPING ELEMENTS ON SOURCE SIDE:',i3,/,
     1' SOLID OR GAS STOPPING POWER TRIGER:',i3,/,
     2' MAXIMUM ENERGY FOR ALPHA SPECTRA:',f7.3,/)
  133 format(' ISOTOPE IDENTIFICATION:',i6,/
     1' ISOTOPE FRACTION:',e15.6,/)
  134 format(' NUMBER OF SOURCES:',i3)
  135 format(/,' SOURCE IDENTIFICATION: ',i6,/
     1' SOURCE ISOTOPIC FRACTION:',e15.6)
  136  format(//,' ALPHA SPECTRA ENERGY BOUNDARIES',/)
  137  format(5e14.6)
  138 format(' MAXIMUM ALPHA DECAY ENERGIES')
  160 format (56h failed to find alpha/s.f./dn source data on tape 5 for
     1 ,i8,5h stop)

 1900 format('SOURCES 4A Calculation',/,'<<<<<<<<<<<>>>>>>>>>>>',//)
 2000 format(/,22x,'Summary of Input')
 2010 format(22x,'================')
 2020 format(' ')
 2021 format(/,80(1h-))
 2030 format('Title: ', 11a7,/)
 2043 format('Interface problem input (idd=',i2,')')
 2048 format('Magnitudes and spectra computed.')
 2049 format('Magnitudes only computed.')
 2050 format('Number of elemental constituents on source side:',i3)
 2060 format('Solid stopping cross-sections used (isg=',i2,') on',
     1 1x,'source side.')
 2070 format('Gas stopping cross-sections used (isg=',i2,')',
     1 1x,'on sources side.')
 2072 format('Maximum energy for alpha spectra: ',1pE10.3,' MeV.')
 2074 format('Minimum energy for alpha spectra: ',1pE10.3,' MeV.')
 2080 format('Elemental Constituents on Source Side:')
 2090 format(' ')
 2100 format('  Z-value  Atom Fraction')
 2110 format('  -------  -------------')
 2120 format(4x,i3,5x,f12.10)
 2200 format('Number of source nuclides to be evaluated:',i3)
 2220 format('Source Nulcides:')
 2230 format('    ZAID     Atom Fraction')
 2240 format('    ----     -------------')
 2250 format(3x,i6,4x,1pE10.3)
 2270 format(i5,' alpha energy groups used at interface.')

c-----------------
c  Set Parameters
c-----------------
      rewind 2
      rewind 3
      rewind 4
      rewind 5
c     aneut=1.
      alph=4.
      zalp=2.
c     pi=3.14159
  400 continue

c------------------
c  Read from Tape1
c------------------
      write(6,1900)
      write(7,1900)
      write(8,1900)
      write(10,1900)
      write(11,1900)
      write(11,2000)
      write(11,2010)
      write(11,2020)
      write(11,2030)title
      write(11,2043)idd
      if (id.eq.1) then
       write(11,2049)
      elseif (id.eq.2) then
       write(11,2048)
       write(7,97)
       write(8,98)
      else
       stop 'Error in id input!'
      endif
      write(6,96)
      write(6,2030)title
      write(7,2030)title
      write(8,2030)title
      read (1,*) nz,isg,em,eamin
      write(11,2050)nz
      if (isg.eq.0) then
       write(11,2060)isg
      elseif (isg.eq.1) then
       write(11,2070)isg
      else
       stop 'Error in isg input!'
      endif
      write(11,2072)em
      write(11,2074)eamin
      write(11,2090)
      write(11,2080)
      write(11,2090)
      write(11,2100)
      write(11,2110)
      if (eamin.lt.0.001) eamin=0.001

c---------------------------------
c read input material constituents
c---------------------------------
      do 420 j=1,nz
       read (1,*) jzm(j),azm(j)
       write(11,2120)jzm(j),azm(j)
  420 continue

c----------------------------------------
c order the nz material constituents by z
c----------------------------------------
      do 440 j=1,nz
       isw=0
       do 430 jn=2,nz
        if (jzm(jn).gt.jzm(jn-1)) go to 430
        jtem=jzm(jn-1)
        atem=azm(jn-1)
        jzm(jn-1)=jzm(jn)
        azm(jn-1)=azm(jn)
        jzm(jn)=jtem
        azm(jn)=atem
        isw=1
  430  continue
       if (isw.eq.0) go to 450
  440 continue

c----------------------------------------------------
c read tape 2 for stopping cross-section coefficients
c----------------------------------------------------
  450 rewind 2
  460 read (2,10) icont
      if (icont.ne.0) go to 460
      do 520 j=1,nz
  470 read (2,40,end=510) iz,c
  480 if (jzm(j)-iz) 510,490,470
  490 continue
      do 500 jt=1,9
  500 czm(jt,j)=c(jt)
      if (isg.eq.1.and.c(10).gt.0.) then
      do 505 jt=1,5
  505 czm(jt,j)=c(jt+9)
      endif
      go to 520
  510 write (6,50) jzm(j)
      stop 1
  520 continue
      read (1,*)nag
      write(11,2090)
      write(11,2270)nag
      do 948 m=1,nag
  948 ast(m)=0.0

c---------------------------------------------------------
c set alpha energy grid ee(m) for this target/source combo
c---------------------------------------------------------
      eamax=em
  860 fnag=nag
      nagp1=nag+1
      dea=(eamax-eamin)/fnag
      ee(1)=eamin
      do 870 m=2,nag
       fm=m-1
  870 ee(m)=ee(1)+dea*fm
      ee(nagp1)=eamax
      do 871 m=1,nag
  871 ea(m)=(ee(m+1) + ee(m))/2.

c-------------------------------------------------------------
c for each energy ee(m) calc. and store stopping cross section
c-------------------------------------------------------------
      do 930 m=1,nagp1
       scx(m)=0.
       do 920 j=1,nz
        zmat=jzm(j)
        amat=amass(jzm(j))
c-----  nuclear stopping---
        term=zalp**0.66667 + zmat**0.66667
        sterm=sqrt(term)
        bot=zalp*zmat*(alph+amat)*sterm
        rep=32530.*amat*ee(m)/bot
        if(rep.lt.0.001) then
         dcx=1.593*sqrt(rep)
         go to 905
         endif
        if(rep.lt.10.) then
         bot=1.+6.8*rep+3.4*rep**1.5
         dcx=1.7*sqrt(rep)*alog(rep+2.71828)/bot
         go to 905
        endif
        dcx=alog(0.47*rep)/(2.*rep)
c-----  electronic stopping---
  905   if(ee(m).gt.30.) go to 906
        slow=czm(1,j)*(1000.*ee(m))**czm(2,j)
        shigh=(czm(3,j)/ee(m))
        shigh=shigh*alog(1.+czm(4,j)/ee(m)+czm(5,j)*ee(m))
        dcx=dcx+slow*shigh/(slow+shigh)
        go to 907
  906   eil=alog(1/ee(m))
        arg=czm(6,j)+eil*(czm(7,j)+czm(8,j)*eil+czm(9,j)*eil*eil)
        dcx=dcx+exp(arg)
  907  continue
  920  scx(m)=scx(m)+azm(j)*dcx
  930 continue

c-------------------------------------------
c read alpha sources from user input
c---------------------------------------------
  580 read (1,*) nq
      write(11,2200)nq
      write(11,2090)
      write(11,2220)
      write(11,2090)
      write(11,2230)
      write(11,2240)
      do 590 k=1,nq
       read (1,*) jq(k),faq(k)
       write(11,2250)jq(k),faq(k)
  590 continue

c-----------------------------------------------
c order sources by z-a-state id = s+10*a+10000*z
c-----------------------------------------------
      if (nq.eq.1) go to 620
      do 610 k=1,nq
       isw=0
       do 600 kn=2,nq
        if (jq(kn).gt.jq(kn-1)) go to 600
        jtem=jq(kn-1)
        atem=aq(kn-1)
        jq(kn-1)=jq(kn)
        aq(kn-1)=aq(kn)
        jq(kn)=jtem
        aq(kn)=atem
        isw=1
  600  continue
       if (isw.eq.0) go to 620
  610 continue

c----------------------------------------
c read target information from user input
c----------------------------------------
  620 continue

c-------------------------------------------
c read source nuclide decay data from tape 5
c-------------------------------------------
  770 read (5,10) icont
      if (icont.ne.0) go to 770
      do 1100 k=1,nq
  790 read (5,110,end=830) idq,nal,ndn
  800 read (5,60) alam,bfsf,barnu,a,b,bfdn
      if (nal.eq.0) go to 810
      read (5,120) (eala(l),fal(l),l=1,nal)
  810 if (ndn.eq.0) go to 820
      read (5,120) (fdng(nd),nd=1,ndn)
  820 if (idq-jq(k)) 790,840,830
  830 write (6,160) jq(k)
      stop 4
  840 continue
      if(nal.eq.0)write(6,130)
      if(nal.eq.0)go to 841
  841 continue

c---------------------------------------------------------
c for each energy group cal intergal 1/scx
c---------------------------------------------------------
      r(1)=1./scx(1)
      fact=alam*faq(k)*2.5e+20
      do 940 m=1,nag
      ps(m)= 0.0
      r(m+1)= 1.0/scx(m+1)
      p(m)=fact*(r(m+1)+r(m))*dea/2.
  940 continue
c---------------------------------------------------------
      do 946 mm = 1,nal
       if(eala(mm).gt.em)write(6,131)em,em
       if(eala(mm).gt.9.8)eala(mm)=9.8               !modified - vk
      do 945  m = 1,nag
      f=1.
      if(eala(mm).ge.ee(m+1))go to 944
      if(eala(mm).gt.ee(m))f=(eala(mm)-ee(m))/dea
      if(eala(mm).lt.ee(m))f=0.0
  944 ps(m) = ps(m) +fal(mm)*p(m)*f
  945 continue
  946 continue
      do 947 m=1,nag
      ast(m) = ast(m) + ps(m)
  947 continue
c----------------------------------------------------------
 1100 continue
      knt=nag
      do 2 i=1,nag
       if(ast(knt).gt.0.0)go to 3
       knt=knt-1
 2    continue
 3    nag=knt
      if(nag.eq.0)write(6,130)
      if(nag.eq.0)stop 20
      do 1 i=1,nag
       eal(i)=ea(i)
 1    continue

c-------------------------
c  End Execution of SAINF
c-------------------------
      return
 1440 continue
      stop
      end



c=======================================================================
c  Interface Problem Neutron Source Subroutine (8/96)
c=======================================================================

      subroutine sninf(nal,title,id,ast,eal)

c---------
c storage
c---------
      character*7 title
      character*2 jsm(105)
      character*1 amt
      dimension title(11), jzm(20), azm(20), czm(9,20), c(14),
     1 el(20), ep(200), f(20,200), e(1100), x(1100),
     2 ee(4001), scx(4001), en(751), san(750), eal(4001),
     3 ts(750), cx(4001), r(4001), rr(4001), p(4001), tsan(750),
     4 gtsan(750), s(750), totlev(20), ast(4001),
     5 sl(750,20), etopl(20), ebotl(20)
      common /names/ jsm
      common /masses/ amass(105)
c     dimension frclev(20)

c--------------------
c formats statements
c--------------------
   10 format (11a7)
   15 format(//,11a7,/)
   20 format (3i3)
   30 format (i8,e12.5)
   40 format(i3,8e8.3,f10.3,/,3x,5f8.3)
   50 format (69h failed to find stopping cross section coefficients on
     1tape 2 for iz=,i3,5h stop)
   60 format (6e12.5)
   70 format(/,1x,'Neutron Multigroup Structure (MeV)')
   80 format (1p8e10.3)
   81 format(31x,'Total (all groups): ',1pE10.3,' neutrons/sec-cm^2.',
     1/,27x,'Average Neutron Energy: ',1pE10.3,' MeV.')
   90 format (i3,5e12.5)
   95 format(2i6)
  100 format(//,1h1,30x,'Table I',/,31x,7(1h=),//,12x,'(alpha,n) ',
     1'Neutron Production by Target and Alpha Energy',///,17x,
     2'target            alpha   alphas/sec     p(e)     neuts/sec',
     3/,7x,'target  atom frac.          energy     /cm^2',4x,
     4'neut/alpha',4x,'/cm^2',/,1h+,6x,6(1h_),2x,10(1h_),2x,6(1h_),2x,
     5 6(1h_),3(12h  __________))
c 102 format(1h+,77x,37hfractional ground/excited level split,/,
c    1 1h+,77x,49(1h_))
  105 format(i8,i4)
  120 format (8e10.3)
  130 format (52h failed to find alpha-n cross section on tape 3 for ,i8
     1 ,5h stop)
  140 format (i8,2i3)
  150 format (66h failed to find comp.nucl.neut.emiss.lev. branching on
     1tape 4 for ,i8,5h stop)
  170 format (/,17h alpha energy of ,1pe12.5,62h exceeds current 6.5 MeV
     1 limit of most cross sections, stop.  )
  180 format (7h target,i8,40h cross-section data not available at ee(
     1 ,i4,3h)= ,1pe12.5,6h stop.)
  200 format(33x,f8.3,1p3e12.4)
c 205 format(1h+,(76x,7f7.5,/,1x))
  210 format(7x,a2,i3,a1,1pe12.4,8x,         0pf8.3,1p3e12.4)
  220 format (14h alpha energy ,1pe12.5,45h exceeds range of level branc
     1hing data, stop.)
  240 format(/,' (a,n) neutrons from alphas on ',a2,i3,a1,' in target.')
  250 format(1h+,66x,10(1h_),/,55x,'Total:',4x,1pe12.4,/)
  251 format(1h+,66x,10(1h_),/,41x,'Total (this target):',
     14x,1pe12.4,/)
  252 format (1h+,66x,10(1h_),/,41x,'Total (all targets):',
     14x,1pe12.4,/)
  255 format(1h )
  260 format (///,' Total (alpha,n) neutron spectrum this target')
  270 format(///,' Grand total (alpha,n) neutron spectrum, all targets,'
     1 ,1x,'all sources')
  390 format (////,' Total Neutron Spectrum')

 2021 format(/,80(1h-))
 2030 format('Target title: ', 11a7)
 2050 format('Number of elemental constituents on target side:',i3)
 2060 format('Solid stopping cross-sections used (isg=',i2,') on',
     1 1x,'target side.')
 2070 format('Gas stopping cross-sections used (isg=',i2,')',
     1 1x,'on target side.')
 2082 format('Elemental Constituents on Target Side:')
 2090 format(' ')
 2100 format('  Z-value  Atom Fraction')
 2110 format('  -------  -------------')
 2120 format(4x,i3,5x,f12.10)
 2130 format('Number of neutron spectrum energy groups:', i4)
 2140 format('Maximum neutron energy is ',1pE12.5,' MeV.')
 2150 format('Minimum neutron energy is ',1pE12.5,' MeV.')
 2160 format('Energy Group Structure:')
 2170 format('  Group  Upper-Bound  Lower-Bound')
 2180 format('  -----  -----------  -----------')
 2190 format(2x,i4,4x,1pE12.5,3x,1pE12.5)
 2193 format(/,'Warning: First Group Upper-Bound Does Not Coincide ',
     1'With Maximum Neutron Energy.')
 2195	format(/,'Warning: Energy Group Structure Not Entered in ',
     1'Decreasing Order.',/,5x,'Group ',i3,' Upper-Bound is Greater ',
     2'Than Group ',i3,' Upper-Bound.') 
 2260 format(/,'Number of target nuclides to be used:',i3)
 2270 format(i5,' alpha energy groups used in target calculation.')
 2280 format('Target Nuclides:')
 2290 format('    ZAID     Atom Fraction')
 2300 format('    ----     -------------')
 2310 format(3x,i6,4x,1pE12.5)

 3000 format(//,22x,'Summary of Output')
 3010 format(22x,'=================')
 3020 format('Total (alpha,n) neutron source from all sources and ',
     1'targets: ',1pE12.5,' n/sec-cm^2.')
 3030 format(/,'Average (alpha,n) neutron energy: ',1pE12.5,' MeV.')
 3050 format(/,'Normalized Neutron Energy Spectrum by Energy Group ',
     1'for All Sources:',//,10x,'Group',5x,'Contribution',/,10x,
     2'-----',5x,'------------')
 3051 format(11x,i3,7x,1pE12.5)
 3500 format(/,'Portion of Total Neutron Source Rate Accounted for ',
     1'in the Energy Spectrum: ',F5.1,'%.')
 3510 format(/,'Warning: Energy Group Structure May Be Poorly Chosen.')

c-----------------
c  Set Parameters
c-----------------
      rewind 2
      rewind 3
      rewind 4
      rewind 5
      aneut=1.
      alph=4.
      zalp=2.
  400 gttqan=0.

c------------------
c  Read from Tape1
c------------------
      read (1,10)title
      write(11,2090)
      write(11,2030)title
      read (1,*) nz,isg
      write(11,2050)nz
      if (isg.eq.0) then
      write(11,2060)isg
      elseif (isg.eq.1) then
       write(11,2070)
      else
       stop 'Error in isg input!'
      endif
      write(11,2090)
      write(11,2082)
      write(11,2090)
      write(11,2100)
      write(11,2110)
  410 continue

c---------------------------------
c read input material constituents
c---------------------------------
      do 420 j=1,nz
       read (1,*) jzm(j),azm(j)
  420 write(11,2120)jzm(j),azm(j)

c----------------------------------------
c order the nz material constituents by z
c----------------------------------------
      do 440 j=1,nz
       isw=0
       do 430 jn=2,nz
        if (jzm(jn).gt.jzm(jn-1)) go to 430
        jtem=jzm(jn-1)
        atem=azm(jn-1)
        jzm(jn-1)=jzm(jn)
        azm(jn-1)=azm(jn)
        jzm(jn)=jtem
        azm(jn)=atem
        isw=1
  430  continue
       if (isw.eq.0) go to 450
  440 continue

c----------------------------------------------------
c read tape 2 for stopping cross-section coefficients
c----------------------------------------------------
  450 rewind 2
  460 read (2,20) icont
      if (icont.ne.0) go to 460
      do 520 j=1,nz
  470 read (2,40,end=510) iz,c
  480 if (jzm(j)-iz) 510,490,470
  490 continue
      do 500 jt=1,9
  500 czm(jt,j)=c(jt)
      if (isg.eq.1.and.c(10).gt.0.) then
      do 505 jt=1,5
  505 czm(jt,j)=c(jt+9)
      endif
      go to 520
  510 write (6,50) jzm(j)
      stop 1
  520 continue
  525 if (id.eq.1) go to 580

c------------------------------------------------------------
c construct neutron group structure desired from tape 1 input
c neutron groups ordered in decreasing energy
c------------------------------------------------------------
      read (1,*) nng,enmax,enmin
      write(11,2090)
      write(11,2130)nng
      write(11,2140)enmax
      write(11,2150)enmin
      write(11,2090)
      write(11,2160)
      write(11,2090)
      write(11,2170)
      write(11,2180)
      nngp1=nng+1
      if (nng.gt.0) go to 540
      nng=-nng
      nngp1=nng+1
      read (1,*) (en(n),n=1,nng)
      en(nngp1)=enmin
      if (en(1).gt.en(2)) go to 560
      nh=nngp1/2
      do 530 n=1,nh
       etem=en(nngp1-n+1)
       en(nngp1-n+1)=en(n)
  530 en(n)=etem
      go to 560
  540 en(nngp1)=enmin
      fnng=nng
      den=(enmax-enmin)/fnng
      en(1)=enmax
      do 550 n=2,nng
       fnm1=n-1
  550 en(n)=en(1)-fnm1*den
  560 do 562 n=1,nng
       write(11,2190)n,en(n),en(n+1)
  562 continue
      if (en(1).ne.enmax) write(11,2193)
	do 528 n=2,nng
	 nm1=n-1
	 if (en(n).gt.en(n-1)) then
	  write(11,2195)n,nm1
	  stop 'Energy Structure Incorrectly Entered!'
	 endif
  528 continue
      write (7,70)
      write (7,80) (en(n),n=1,nngp1)
      write (8,70)
      write (8,80) (en(n),n=1,nngp1)

c----------------------------
c zero total spectrum storage
c----------------------------
      do 570 n=1,nng
       ts(n)=0.
       gtsan(n)=0.
  570 etopan=0.
      ebotan=0.
      gttqan=0.

c----------------------------------------------------
c    alpha sources treated as if one isotope ie nq=1
c----------------------------------------------------
  580 continue
      nq=1

c----------------------------------------
c read target information from user input
c----------------------------------------
  620 continue
      read (1,*) nt, nag
      write(11,2260)nt
      write(11,2270)nag
      write(11,2090)
      write(11,2280)
      write(11,2090)
      write(11,2290)
      write(11,2300)
      if (nt.eq.0) go to 1440

c----------------------------------------------------
c beginning of big loop over target nuclides
c         for nt=0 calculate spontaneous fission only
c----------------------------------------------------
      rewind 3
      if (id.eq.2) rewind 4
      write (6,100)
c     if(id.eq.2) write(6,102)

c--------------------------
c loop on target nuclides i
c--------------------------
  630 read (3,20) icont
      if (icont.ne.0) go to 630
  640 if (id.eq.2) read (4,20) icont
      if (icont.ne.0) go to 640
      do 1120 i=1,nt
      etpant=0.
      ebtant=0.
      totqan=0.
      if(id.eq.1) go to 660
      do 650 n=1,nng
  650 tsan(n)=0.

c------------------------------------------------------
c find alpha-n cross sections for this target on tape 3
c------------------------------------------------------
  660 read (1,*) idt,at
      write(11,2310)idt,at
  665 read (3,105,end=680) jdt,jps
  670 read (3,120) (e(ip),x(ip),ip=1,jps)
      if (idt-jdt) 680,690,665
  680 write (6,130) idt
      stop 2
  690 lzt=idt/10000
      lat=idt/10-lzt*1000
      atar=lat
      lst=idt-lzt*10000-lat*10
      amt=' '
      if (lst.ne.0) amt='m'
      apro=atar+3.
      do 700 ip=2,jps
      eamin=e(ip-1)
      if (x(ip).ne.0.) go to 710
  700 continue
  710 if (eamin.lt.0.001) eamin=0.001
      if (id.eq.1) go to 760

c----------------------------------------------------------------
c find product nuclide level branchings for this target of tape 4
c----------------------------------------------------------------
  720 read (4,140,end=750) jdt,jl,jp
  730 read (4,120) q,(el(il),il=1,jl)
      do 740 ip=1,jp
  740 read (4,120) ep(ip),(f(il,ip),il=1,jl)
      if (idt-jdt) 750,760,720
  750 write (6,150) idt
      stop 3
  760 rewind 5
      if (nq.ne.0) go to 770
      sbtqan=0.
  770 read (5,20) icont
      if (icont.ne.0) go to 770

c--------------------------------------------
c loop on source nuclides k for this target i
c contribution to (alpha,n) neutrons
c--------------------------------------------
      do 1100 k=1,nq
      etop=0.
      ebot=0.
      do 775 il=1,jl
      etopl(il)=0.
      ebotl(il)=0.
      do 774 n=1,nng
  774 sl(n,il)=0.
  775 totlev(il)=0.
      sbtqan=0.
      if(id.eq.1) go to 790
      do 780 n=1,nng
  780 san(n)=0.
  790 continue

c---------------------------------------------------------
c set alpha energy grid ee(m) for this target/source combo
c---------------------------------------------------------
      eamax=eal(nal)
  850 continue
      if(eamax.le.eamin) then
      qan=0.
      mm=0
      p(nagp1)=0.
      pval=0.
      go to 945
      endif
      if (eamax.le.9.8) go to 860                     !modified - vk
      write (6,170) eamax
      stop 5
  860 fnag=nag
      nagp1=nag+1
      dea=(eamax-eamin)/fnag
      ee(1)=eamin
      do 870 m=2,nag
      fm=m-1
  870 ee(m)=ee(1)+dea*fm
      ee(nagp1)=eamax

c----------------------------------------------------
c for each energy ee(m) calc. and store cross section
c----------------------------------------------------
      ip=1
      do 900 m=1,nagp1
  880 if (ee(m).ge.e(ip).and.ee(m).le.e(ip+1)) go to 890
      ip=ip+1
      if (ip.lt.jps) go to 880
      write (6,180) idt,m,ee(m)
      stop 6
  890 slope=(x(ip+1)-x(ip))/(e(ip+1)-e(ip))
      enrsep=(x(ip)*e(ip+1)-x(ip+1)*e(ip))/(e(ip+1)-e(ip))
  900 cx(m)=enrsep+ee(m)*slope

c-------------------------------------------------------------
c for each energy ee(m) calc. and store stopping cross section
c-------------------------------------------------------------
      do 930 m=1,nagp1
      scx(m)=0.
      do 920 j=1,nz
      zmat=jzm(j)
      amat=amass(jzm(j))
c-----nuclear stopping---
      term=zalp**0.66667 + zmat**0.66667
      sterm=sqrt(term)
      bot=zalp*zmat*(alph+amat)*sterm
      rep=32530.*amat*ee(m)/bot
      if(rep.lt.0.001) then
      dcx=1.593*sqrt(rep)
      go to 905
      endif
      if(rep.lt.10.) then
      bot=1.+6.8*rep+3.4*rep**1.5
      dcx=1.7*sqrt(rep)*alog(rep+2.71828)/bot
      go to 905
      endif
      dcx=alog(0.47*rep)/(2.*rep)
c-----electronic stopping---
  905 if(ee(m).gt.30.) go to 906
      slow=czm(1,j)*(1000.*ee(m))**czm(2,j)
      shigh=(czm(3,j)/ee(m))*alog(1.+czm(4,j)/ee(m) + czm(5,j)*ee(m))
      dcx=dcx+slow*shigh/(slow+shigh)
      go to 907
  906 eil=alog(1/ee(m))
      arg=czm(6,j)+czm(7,j)*eil+czm(8,j)*eil*eil+czm(9,j)*eil*eil*eil
      dcx=dcx+exp(arg)
  907 continue
  920 scx(m)=scx(m)+azm(j)*dcx
  930 continue

c-----------------------------------------------------------
c for each energy ee(m) calc. ratio r(m), then integral p(m)
c-----------------------------------------------------------
      r(1)=cx(1)/scx(1)
      p(1)=0.
      fact=1.0e-06*at
      do 940 m=2,nagp1
      r(m)=cx(m)/scx(m)
  940 p(m)=p(m-1)+fact*(r(m-1)+r(m))*dea/2.
  945 continue

c----------------------------------------------
c  calculate neutrons from alphas
c----------------------------------------------
  950 do 1080 l=1,nal
      if (nq.eq.0) go to 1000
      do 960 m=1,nag
      mm=m
  960 if (eal(l).ge.ee(m).and.eal(l).le.ee(m+1)) go to 970
      mm=0
      pval=0.
      go to 975
  970 pval=p(mm)+(eal(l)-ee(mm))*(p(mm+1)-p(mm))/(ee(mm+1)-ee(mm))
  975 aps=ast(l)
      qan=aps*pval
      sbtqan=sbtqan+qan
      totqan=totqan+qan
      if (l.eq.1) go to 980
      write (6,200) eal(l),aps,pval,qan
      if(nal.ne.1.and.l.eq.nal) write(6,250) sbtqan
      go to 990
  980 write(6,210)jsm(lzt),lat,amt,at,eal(l),aps,pval,qan
      if(nal.eq.1) write(6,255)
  990 if (id.eq.1.or.mm.eq.0) go to 1080

c-----------------------------------------------------------------
c calculate (alpha,n) neutron spectrum contribution in multigroups
c-----------------------------------------------------------------
 1000 mmm1=mm-1
      do 1010 m=1,mmm1
 1010 rr(m)=(p(m+1)-p(m))/pval
      rr(mm)=(pval-p(mm))/pval
      do 1020 n=1,nng
 1020 s(n)=0.
      do 1060 il=1,jl
      qlev=q-el(il)
      e90=-qlev*apro/(apro-alph)
      thre=-qlev*(aneut+apro)/(aneut+apro-alph)
      if (qlev.gt.0.) thre=0.
      do 1060 m=1,mm
      ea=(ee(m)+ee(m+1))/2.
      if (m.eq.mm) ea=(eal(l)+ee(mm))/2.
      if (ea.lt.thre) go to 1055
      do 1030 ip=2,jp
      lp=ip
 1030 if (ep(ip-1).le.ea.and.ep(ip).ge.ea) go to 1040
      write (6,220)ea
      stop 7
 1040 lp1=lp-1
      bx=f(il,lp1)+(ea-ep(lp1))*(f(il,lp)-f(il,lp1))/(ep(lp)-ep(lp1))
      term1=sqrt(4.*alph*aneut*ea)/(2.*(aneut+apro))
      term2=alph*aneut*ea/((aneut+apro)*(aneut+apro))
      term3=(apro*ea+apro*qlev-alph*ea)/(aneut+apro)
      senmax=term1+sqrt(term2+term3)
      enmax=senmax*senmax
      if (ea.le.e90) senmin=term1-sqrt(term2+term3)
      if (ea.gt.e90) senmin=-term1+sqrt(term2+term3)
      enmin=senmin*senmin
      ebar=(enmin+enmax)/2.
      val=qan*rr(m)*bx
      etop=etop+val*ebar
      etopl(il)=etopl(il)+val*ebar
      ebotl(il)=ebotl(il)+val
      ebot=ebot+val
      dele=enmax-enmin
      do 1050 n=1,nng
      if (en(n).lt.enmin) go to 1050
      if (en(n+1).gt.enmax) go to 1050
      de=en(n)-en(n+1)
      if (en(n+1).lt.enmin) de=de-(enmin-en(n+1))
      if (en(n).gt.enmax) de=de-(en(n)-enmax)
      gpadd=qan*rr(m)*bx*de/dele
      s(n)=s(n)+gpadd
      sl(n,il)=sl(n,il)+gpadd
      totlev(il)=totlev(il)+gpadd
 1050 continue
 1055 continue
 1060 continue
      do 1070 n=1,nng
      san(n)=san(n)+s(n)
      gtsan(n)=gtsan(n)+s(n)
      ts(n)=ts(n)+s(n)
 1070 tsan(n)=tsan(n)+s(n)
 1075 continue
 1080 continue
      if(id.eq.1) go to 1100

c----------------------------------------------------
c output this (alpha,n) neutron spectrum contribution
c for this source/target combination
c----------------------------------------------------
      ebarqt=0.
      if(ebot.le.0.) go to 1085
      ebarqt=etop/ebot
      etpant=etpant+etop
      ebtant=ebtant+ebot
 1085 continue
      smga=0.
      do 1090 n=1,nng
 1090 smga=smga+san(n)
c     fracgp=smga/sbtqan
c     do 1095 il=1,jl
c1095 frclev(il)=totlev(il)/sbtqan
c     write(6,205) (frclev(il),il=1,jl)
      if (nq.ne.0) write (7,240) jsm(lzt),lat,amt
      write (7,80) (san(n),n=1,nng)
      write (7,81) smga,ebarqt
      smga=0.
      if (nq.ne.0) write (8,240) jsm(lzt),lat,amt
      do 1097 n=1,nng
      san(n)=san(n)/sbtqan
 1097 smga=smga+san(n)
      write (8,80) (san(n),n=1,nng)
      write (8,81) smga, ebarqt
      do 1098 il=1,jl
      if(ebotl(il).le.0.) go to 1098
c     ebarl=etopl(il)/ebotl(il)
 1098 continue
 1100 continue
 1101 gttqan=gttqan+totqan

c--------------------------------------
c output total (alpha,n) source for all
c alpha sources on this target nuclide
c--------------------------------------
      if(id.eq.1) go to 1120
      tmga=0.
      do 1110 n=1,nng
 1110 tmga=tmga+tsan(n)
c     fracgp=tmga/totqan
      ebart=etpant/ebtant
      etopan=etopan+etpant
      ebotan=ebotan+ebtant
      if(nq.le.1) go to 1120
      write(7,260)
      write(7,80) (tsan(n),n=1,nng)
      write(7,81) tmga,ebart
      write(7,2090)
      write(7,2021)
      write(7,2090)
      tmga=0.
      do 1115 n=1,nng
       tsan(n)=tsan(n)/totqan
 1115 tmga=tmga+tsan(n)
      write(8,260)
      write(8,80) (tsan(n),n=1,nng)
      write(8,81) tmga,ebart
      write(8,2090)
      write(8,2021)
      write(8,2090)
      write(10,2090)
      write(10,2021)
      write(10,2090)
 1120 continue

c----------------------------------------------------------------------
c output grand total (alpha,n) neutron source: all targets, all sources
c----------------------------------------------------------------------
      if (nt.gt.1) write (6,252) gttqan
      if(id.eq.1) go to 1440
      gtmga=0.
      do 1130 n=1,nng
 1130 gtmga=gtmga+gtsan(n)
      ebaran=etopan/ebotan
c     fracgp=gtmga/gttqan
      write(7,2090)
      write(7,2021)
      write(7,390)
      write (7,80) (gtsan(n),n=1,nng)
      write (7,81) gtmga,ebaran
      gtmga=0.
	dummy=0.0
      do 1135 n=1,nng
	 dummy=dummy+gtsan(n)
       gtsan(n)=gtsan(n)/gttqan
 1135 gtmga=gtmga+gtsan(n)
      write(8,2090)
      write(8,2021)
      write(8,390)
      write(8,80) (gtsan(n),n=1,nng)
      write(8,81) gtmga,ebaran

c-----------------
c  Output Summary
c-----------------
 1440 continue
      write(11,2090)
      write(11,2090)
      write(11,3000)
      write(11,3010)
      write(11,2090)
      write(11,3020)gttqan
      if (id.gt.1) then
       write(11,3030)ebaran
	 write(11,3050)
	 do 1450 n=1,nng
	  write(11,3051)n,gtsan(n)
 1450  continue
       frac=100.0*(1.0-(abs(dummy-gttqan)/gttqan))
	 write(11,3500)frac
	 if (frac.le.85.0) write(11,3510)
      endif

c-------------------------
c  End Execution of SNINF
c-------------------------
      return
      end


	
c=======================================================================
c  Error Function Routines (7/97)
c=======================================================================

      FUNCTION erf1(x)
	  REAL erf1,x
	  REAL gammp
	  if(x.lt.0.)then
	   erf1=-gammp(.5,x**2)
	  else
	   erf1=gammp(.5,x**2)
	  endif
      return
      END


      FUNCTION gammp(a,x)
	  REAL a,gammp,x
	  REAL gammcf,gamser,gln
	  if(x.lt.0..or.a.le.0.)pause 'bad arguments in gammp'
	  if(x.lt.a+1.)then
	  call gser(gamser,a,x,gln)
	  gammp=gamser
	  else
	  call gcf(gammcf,a,x,gln)
	  gammp=1.-gammcf
	  endif
      return
      END

      SUBROUTINE gcf(gammcf,a,x,gln)
	  INTEGER ITMAX
	  REAL a,gammcf,gln,x,EPS,FPMIN
	  PARAMETER (ITMAX=100,EPS=3.e-7,FPMIN=1.e-30)
	  INTEGER i
	  REAL an,b,c,d,del,h,gammln
	  gln=gammln(a)
	  b=x+1.-a
	  c=1./FPMIN
	  d=1./b
	  h=d
	  do 11 i=1,ITMAX
	   an=-i*(i-a)
	   b=b+2.
	   d=an*d+b
	   if(abs(d).lt.FPMIN)d=FPMIN
	   c=b+an/c
	   if(abs(c).lt.FPMIN)c=FPMIN
	   d=1./d
	   del=d*c
	   h=h*del
	   if(abs(del-1.).lt.EPS)goto 1
11    continue
	  pause 'a too large, ITMAX too small in gcf'

1     gammcf=exp(-x+a*log(x)-gln)*h
      return
      END


      SUBROUTINE gser(gamser,a,x,gln)
	  INTEGER ITMAX
	  REAL a,gamser,gln,x,EPS
	  PARAMETER (ITMAX=100,EPS=3.e-7)
	  INTEGER n
	  REAL ap,del,sum,gammln
	  gln=gammln(a)
	  if(x.le.0.)then
	   if(x.lt.0.)pause 'x < 0 in gser'
	   gamser=0.
	   return
	  endif
	  ap=a
	  sum=1./a
	  del=sum
	  do 11 n=1,ITMAX
	   ap=ap+1.
	   del=del*x/ap
	   sum=sum+del
	   if(abs(del).lt.abs(sum)*EPS)goto 1
11    continue
	  pause 'a too large, ITMAX too small in gser'
1     gamser=sum*exp(-x+a*log(x)-gln)

      return
      END

      FUNCTION gammln(xx)
	  REAL gammln,xx
	  INTEGER j
	  DOUBLE PRECISION ser,stp,tmp,x,y,cof(6)
	  SAVE cof,stp
	  DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,
     *24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,
     *-.5395239384953d-5,2.5066282746310005d0/
	  x=xx
	  y=x
	  tmp=x+5.5d0
	  tmp=(x+0.5d0)*log(tmp)-tmp
	  ser=1.000000000190015d0
	  do 11 j=1,6
	   y=y+1.d0
	   ser=ser+cof(j)/y
11    continue
	  gammln=tmp+log(stp*ser/x)

      return
      END


c=======================================================================
c  Three Region Problem Main Subroutine (1/99)
c=======================================================================

      subroutine three(title,idd,id)

c----------
c  Storage
c----------
      
      character*7 title
      dimension title(11)

c------------------
c  Call Statements
c------------------

      call input(title,idd,id)
      call processor(title,idd,id)

c----------------------------
c  End of Execution of THREE
c----------------------------

      return
      end


c=======================================================================
c  Three Region Problem Input Subroutine (1/99)
c=======================================================================

      subroutine input(title,idd,id)

c----------
c  Storage
c----------
      character*7 title, title1, title2, title3
      character*2 jsm(105)
      character*1 amt(20), amtc(20)
      dimension title(11), title1(11), title2(11), title3(11), 
     1 ee(4001), ea(4001), en(751), cg(4001), jza(20), aza(20),
     2 jq(300), faq(300), c(14), cza(9,20), fal(30,20), eala(30,20), 
     3 fdng(1000), jzb(20), azb(20), czb(9,20), idb(20), atb(20),
     4 jps(20), e(1100,20), x(1100,20), lzt(20), lat(20), lst(20),
     5 apro(20), jl(20), jp(20), q(20), el(20,20), ep(200,20),
     6 f(20,200,20), jzc(20), azc(20), czc(9,20), idc(20), atc(20),
     7 jpsc(20), ec(1100,20), xc(1100,20), lztc(20), latc(20), 
     8 lstc(20), aproc(20), jlc(20), jpc(20), qc(20), elc(20,20), 
     9 epc(200,20), fc(20,200,20), nal(20), alam(20), eal(4001)

       common /names/ jsm
       common /masses/ amass(105)
       common /vars/ jza, aza, jq, faq, cza, fal, eala, fdng, jzb, 
     1 azb, czb, idb, atb, jps, e, x, lzt, lat, lst, apro, jl, jp, 
     2 q, el, ep, f, jzc, azc, czc, idc, atc, jpsc, ec, xc, lztc, 
     3 latc, lstc, aproc, jlc, jpc, qc, elc, epc, fc, nza, isga, nq,
     5 nzb, isgb, anumb, thick, ntb, nzc, isgc, ntc, amt, amtc, nal,
     6 alam, title1, title2, title3
       common /grids/ ee, en, cg, nag, eamax, eamin, nng, 
     1 enmax, enmin, ncg, dea, eal

c--------------------
c  Format Statements
c--------------------
   10 format (i3,11a7)
   15 format(11a7)
   16 format(11a7,//)
   20 format (2i3,f7.3)
   30 format (i8,e12.5)
   40 format(i3,8e8.3,f10.3,/,3x,5f8.3)
   50 format (69h failed to find stopping cross section coefficients on
     1tape 2 for iz=,i3,5h stop)
   60 format (6e12.5)
   70 format(/,1x,'Neutron Multigroup Structure (MeV)')
   80 format (1p8e10.3)
   90 format (i3,5e12.5)
   95 format(2i6)
   96 format(/,26x,'Neutron Source Magnitudes',/,26x,25(1h~),/)
   97 format(/,23x,'Absolute Neutron Source Spectra',/,23x,31(1h~),/)
   98 format(/,22x,'Normalized Neutron Source Spectra',/,22x,
     133(1h~),/)
   99 format(/,16x,'Neutron Source Spectra by Nuclide Energy Level',
     1/,16x,46(1h~),/)
  100 format(//,1h1,30x,'Table I',/,31x,7(1h=),//,12x,'(alpha,n) ',
     1'Neutron Production by Target and Alpha Energy',///,17x,
     2'target            alpha   alphas/sec     p(e)     neuts/sec',
     3/,7x,'target  atom frac.          energy     /cm^2',4x,
     4'neut/alpha',4x,'/cm^2',/,1h+,6x,6(1h_),2x,10(1h_),2x,6(1h_),2x,
     5 6(1h_),3(12h  __________))
  105 format(i8,i4)
  110 format (i8,2i4)
  120 format (8e10.3)
  130 format (1x,'SOURCE DOES NOT HAVE ALPHA EMITTER ', i6)
  131 format(' WARNING-AN ALPHA ENERGY THAT EXCEEDED ',f7.3,' MEV WAS',
     1' SET TO:',f7.3,' MeV')
  132 format(' NUMBER OF STOPPING ELEMENTS ON SOURCE SIDE:',i3,/,
     1' SOLID OR GAS STOPPING POWER TRIGER:',i3,/,
     2' MAXIMUM ENERGY FOR ALPHA SPECTRA:',f7.3,/)
  133 format(' ISOTOPE IDENTIFICATION:',i6,/
     1' ISOTOPE FRACTION:',e15.6,/)
  134 format(' NUMBER OF SOURCES:',i3)
  135 format(/,' SOURCE IDENTIFICATION: ',i6,/
     1' SOURCE ISOTOPIC FRACTION:',e15.6)
  136  format(//,' ALPHA SPECTRA ENERGY BOUNDARIES',/)
  137  format(5e14.6)
  138 format(' MAXIMUM ALPHA DECAY ENERGIES')
  140 format (i8,2i3)
  150 format (66h failed to find comp.nucl.neut.emiss.lev. branching on
     1tape 4 for ,i8,5h stop)
  160 format (56h failed to find alpha/s.f./dn source data on tape 5 for
     1 ,i8,5h stop)

 1900 format('SOURCES 4A Calculation',/,'<<<<<<<<<<<>>>>>>>>>>>',//)
 2000 format(/,22x,'Summary of Input')
 2010 format(22x,'================')
 2020 format(' ')
 2021 format(/,80(1h-))
 2030 format('Title: ', 11a7,/)
 2043 format('Three region problem input (idd=',i2,')')
 2048 format('Magnitudes and spectra computed.')
 2049 format('Magnitudes only computed.')
 2072 format('Maximum energy for alpha spectra: ',1pE10.3,' MeV.')
 2074 format('Minimum energy for alpha spectra: ',1pE10.3,' MeV.')
 2080 format('Elemental Constituents in Region A:')
 2081 format('Elemental Constituents in Region B:')
 2082 format('Elemental Constituents in Region C:')
 2090 format(' ')
 2100 format('  Z-value  Atom Fraction')
 2110 format('  -------  -------------')
 2120 format(4x,i3,5x,f12.10)
 2130 format('Number of neutron spectrum energy groups:', i4)
 2140 format('Maximum neutron energy is ',1pE12.5,' MeV.')
 2150 format('Minimum neutron energy is ',1pE12.5,' MeV.')
 2160 format('Neutron Energy Group Structure:')
 2170 format('  Group  Upper-Bound  Lower-Bound')
 2180 format('  -----  -----------  -----------')
 2190 format(2x,i4,4x,1pE12.5,3x,1pE12.5)
 2193 format(/,'Warning: First Group Upper-Bound Does Not Coincide ',
     1'With Maximum Neutron Energy.')
 2195	format(/,'Warning: Energy Group Structure Not Entered in ',
     1'Decreasing Order.',/,5x,'Group ',i3,' Upper-Bound is Greater ',
     2'Than Group ',i3,' Upper-Bound.') 
 2200 format('Number of source nuclides to be evaluated:',i3)
 2220 format('Source Nuclides in Region A:')
 2230 format('    ZAID     Atom Fraction')
 2240 format('    ----     -------------')
 2250 format(3x,i6,4x,1pE10.3)
 2260 format(/,'Number of target nuclides in region B:',i3)
 2261 format(/,'Number of target nuclides in region C:',i3)
 2270 format(i5,' alpha energy groups used at each interface.')
 2280 format('Target Nuclides in Region B:')
 2281 format('Target Nuclides in Region C:')

 3030 format(/,'Region A Title: ', 11a7,/)
 3040 format(/,'Region B Title: ', 11a7,/)
 3050 format(/,'Region C Title: ', 11a7,/)
 3130 format('Number of angular groups:', i4)
 3160 format('Angular Group Structure:',/)
 3170 format('  Group  Upper-Bound  Lower-Bound')
 3180 format('  -----  -----------  -----------')
 3190 format(2x,i4,4x,1pE12.5,3x,1pE12.5)
 3200 format('Number of elemental constituents in region A:',i3)
 3201 format('Number of elemental constituents in region B:',i3)
 3202 format('Number of elemental constituents in region C:',i3)
 3210 format('Solid stopping cross-sections used (isga=',i2,') in',
     1 1x,'region A.')
 3211 format('Solid stopping cross-sections used (isgb=',i2,') in',
     1 1x,'region B.')
 3212 format('Solid stopping cross-sections used (isgc=',i2,') in',
     1 1x,'region C.')
 3220 format('Gas stopping cross-sections used (isga=',i2,')',
     1 1x,'in region A.')
 3221 format('Gas stopping cross-sections used (isgb=',i2,')',
     1 1x,'in region B.')
 3222 format('Gas stopping cross-sections used (isgc=',i2,')',
     1 1x,'in region C.')
 3231 format('Material B atom density: ',1pE12.5,' atoms/cc.')
 3240 format('Interface region thickness: ',1pE12.5,' cm.')

c-----------------
c  Set parameters
c-----------------
      pi=4.0*atan(1.0)

c----------------
c  Write headers
c----------------
      write(6,1900)
      write(7,1900)
      write(8,1900)
      write(10,1900)
      write(11,1900)
      write(11,2000)
      write(11,2010)
      write(11,2020)
      write(11,2030)title
      write(11,2043)idd
      if (id.eq.1) then
       write(11,2049)
      elseif (id.eq.2) then
       write(11,2048)
       write(7,97)
       write(8,98)
      else
       stop 'Error in id input!'
      endif
      write(6,96)
      write(6,2030)title
      write(7,2030)title
      write(8,2030)title

c--------------------------------------
c  Create alpha energy group structure
c--------------------------------------
      read (1,*)nag,eamax,eamin
      write(11,2090)
      write(11,2270)nag
      write(11,2072)eamax
      write(11,2074)eamin
      if (eamin.lt.0.001) eamin=0.001
  200 fnag=nag
      nagp1=nag+1
      dea=(eamax-eamin)/fnag
      ee(1)=eamin
      do 210 m=2,nag
       fm=m-1
  210 ee(m)=ee(1)+dea*fm
      ee(nagp1)=eamax
      do 211 m=1,nag
  211 ea(m)=(ee(m+1) + ee(m))/2.

c----------------------------------------
c  Create neutron energy group structure
c----------------------------------------
      read (1,*) nng,enmax,enmin
      write(11,2090)
      write(11,2130)nng
      write(11,2140)enmax
      write(11,2150)enmin
      write(11,2090)
      write(11,2160)
      write(11,2090)
      write(11,2170)
      write(11,2180)
      nngp1=nng+1
      if (nng.gt.0) go to 240
      nng=-nng
      nngp1=nng+1
      read (1,*) (en(n),n=1,nng)
      en(nngp1)=enmin
      if (en(1).gt.en(2)) go to 260
      nh=nngp1/2
      do 231 n=1,nh
       etem=en(nngp1-n+1)
       en(nngp1-n+1)=en(n)
  231 en(n)=etem
      go to 260
  240 en(nngp1)=enmin
      fnng=nng
      den=(enmax-enmin)/fnng
      en(1)=enmax
      do 250 n=2,nng
       fnm1=n-1
  250 en(n)=en(1)-fnm1*den
  260 do 262 n=1,nng
       write(11,2190)n,en(n),en(n+1)
  262 continue
      if (en(1).ne.enmax) write(11,2193)
      do 268 n=2,nng
       nm1=n-1
       if (en(n).gt.en(n-1)) then
	write(11,2195)n,nm1
	stop 'Energy Structure Incorrectly Entered!'
       endif
  268 continue
      write (7,70)
      write (7,80) (en(n),n=1,nngp1)
      write (8,70)
      write (8,80) (en(n),n=1,nngp1)

c---------------------------------
c  Create angular group structure
c---------------------------------
      read (1,*)ncg
      write(11,2090)
      write(11,3130)ncg
      write(11,2090)
      write(11,3160)
      write(11,3170)
      write(11,3180)
      fncg=ncg
      cg(1)=0
      do 270 i=1,ncg
        cg(i+1)=cg(i)+(pi/2.0/fncg)
        write(11,3190)i,cg(i+1),cg(i)
  270 continue

c------------------------------------------
c  Read elmental constituents for region A
c------------------------------------------
      read(1,15)title1
      write(11,2090)
      write(11,3030)title1
      read(1,*)nza,isga
      write(11,3200)nza
      if (isga.eq.0) then
       write(11,3210)isga
      elseif (isga.eq.1) then
       write(11,3220)isga
      else
       stop 'Error in isga input!'
      endif
      write(11,2080)
      write(11,2090)
      write(11,2100)
      write(11,2110)
      do 280 j=1,nza
       read (1,*) jza(j),aza(j)
       write(11,2120)jza(j),aza(j)
  280 continue

c----------------------------------------
c order the nz material constituents by z
c----------------------------------------
      do 300 j=1,nza
       isw=0
       do 290 jn=2,nza
        if (jza(jn).gt.jza(jn-1)) go to 290
        jtem=jza(jn-1)
        atem=aza(jn-1)
        jza(jn-1)=jza(jn)
        aza(jn-1)=aza(jn)
        jza(jn)=jtem
        aza(jn)=atem
        isw=1
  290  continue
       if (isw.eq.0) go to 310
  300 continue

c----------------------------------------------------
c read tape 2 for stopping cross-section coefficients
c----------------------------------------------------
  310 rewind 2
  320 read (2,10) icont
      if (icont.ne.0) go to 320
      do 390 j=1,nza
  330  read (2,40,end=380) iz,c
  340  if (jza(j)-iz) 380,350,330
  350  continue
       do 360 jt=1,9
  360  cza(jt,j)=c(jt)
       if (isga.eq.1.and.c(10).gt.0.) then
       do 370 jt=1,5
  370  cza(jt,j)=c(jt+9)
       endif
       go to 390
  380  write (6,50) jza(j)
       stop 1
  390 continue

c-------------------------------------------
c read alpha sources from user input
c---------------------------------------------
      read (1,*) nq
      write(11,2090)
      write(11,2200)nq
      write(11,2090)
      write(11,2220)
      write(11,2090)
      write(11,2230)
      write(11,2240)
      do 400 k=1,nq
       read (1,*) jq(k),faq(k)
       write(11,2250)jq(k),faq(k)
  400 continue

c-----------------------------------------------
c order sources by z-a-state id = s+10*a+10000*z
c-----------------------------------------------
      if (nq.eq.1) go to 430
      do 420 k=1,nq
       isw=0
       do 410 kn=2,nq
        if (jq(kn).gt.jq(kn-1)) go to 410
        jtem=jq(kn-1)
        atem=faq(kn-1)
        jq(kn-1)=jq(kn)
        faq(kn-1)=faq(kn)
        jq(kn)=jtem
        faq(kn)=atem
        isw=1
  410  continue
       if (isw.eq.0) go to 430
  420 continue

c-------------------------------------------
c read source nuclide decay data from tape 5
c-------------------------------------------
  430 continue
  440 read (5,10) icont
      if (icont.ne.0) go to 440
      do 510 k=1,nq
  450  read (5,110,end=480) idq,nal(k),ndn
       read (5,60) alam(k),bfsf,barnu,a,b,bfdn
       if (nal(k).eq.0) go to 460
       read (5,120) (eala(l,k),fal(l,k),l=1,nal(k))
  460  if (ndn.eq.0) go to 470
       read (5,120) (fdng(nd),nd=1,ndn)
  470  if (idq-jq(k)) 450,490,480
  480  write (6,160) jq(k)
       stop 'Source nuclide not in library'
  490  continue
       if(nal(k).eq.0)write(6,130)
       if(nal(k).eq.0)go to 500
  500  continue
  510 continue

c------------------------------------------
c  Read material constituents for region B
c------------------------------------------
      read (1,15)title2
      write(11,2090)
      write(11,3040)title2
      read (1,*) nzb,isgb,anumb,thick
      write(11,3201)nzb
      if (isgb.eq.0) then
      write(11,3211)isgb
      elseif (isgb.eq.1) then
       write(11,3221)isgb
      else
       stop 'Error in isgb input!'
      endif
      write(11,3231)anumb
      write(11,3240)thick
      write(11,2090)
      write(11,2081)
      write(11,2090)
      write(11,2100)
      write(11,2110)

c---------------------------------
c read input material constituents
c---------------------------------
      do 520 j=1,nzb
       read (1,*) jzb(j),azb(j)
  520 write(11,2120)jzb(j),azb(j)

c----------------------------------------
c order the nz material constituents by z
c----------------------------------------
      do 540 j=1,nzb
       isw=0
       do 530 jn=2,nzb
        if (jzb(jn).gt.jzb(jn-1)) go to 530
        jtem=jzb(jn-1)
        atem=azb(jn-1)
        jzb(jn-1)=jzb(jn)
        azb(jn-1)=azb(jn)
        jzb(jn)=jtem
        azb(jn)=atem
        isw=1
  530  continue
       if (isw.eq.0) go to 550
  540 continue

c----------------------------------------------------
c read tape 2 for stopping cross-section coefficients
c----------------------------------------------------
  550 rewind 2
  560 read (2,20) icont
      if (icont.ne.0) go to 560
      do 620 j=1,nzb
  570  read (2,40,end=610) iz,c
  580  if (jzb(j)-iz) 610,590,570
  590  continue
       do 600 jt=1,9
  600  czb(jt,j)=c(jt)
       if (isgb.eq.1.and.c(10).gt.0.) then
        do 605 jt=1,5
  605   czb(jt,j)=c(jt+9)
       endif
       go to 620
  610  write (6,50) jzb(j)
       stop 1
  620 continue
  625 if (id.eq.1) go to 720

c----------------------------------------
c read target information from user input
c----------------------------------------
  720 continue
      read (1,*) ntb
      write(11,2260)ntb
      write(11,2090)
      write(11,2280)
      write(11,2090)
      write(11,2230)
      write(11,2240)

c--------------------------
c loop on target nuclides i
c--------------------------
      rewind 3
      if (id.eq.2) rewind 4
      write (6,100)
  730 read (3,20) icont
      if (icont.ne.0) go to 730
  740 if (id.eq.2) read (4,20) icont
      if (icont.ne.0) go to 740
      do 880 i=1,ntb

c------------------------------------------------------
c find alpha-n cross sections for this target on tape 3
c------------------------------------------------------
  760  read (1,*) idb(i),atb(i)
       write(11,2250)idb(i),atb(i)
  765  read (3,105,end=780) jdt,jps(i)
  770  read (3,120) (e(ip,i),x(ip,i),ip=1,jps(i))
       if (idb(i)-jdt) 780,790,765
  780  write (6,130) idb(i)
       stop 2
  790  lzt(i)=idb(i)/10000
       lat(i)=idb(i)/10-lzt(i)*1000
       atar=lat(i)
       lst(i)=idb(i)-lzt(i)*10000-lat(i)*10
       amt(i)=' '
       if (lst(i).ne.0) amt(i)='m'
       apro(i)=atar+3.
c       do 800 ip=2,jps(i)
c        eamin=e(ip-1,i)
c        if (x(ip,i).ne.0.) go to 810
c  800  continue
  810  if (eamin.lt.0.001) eamin=0.001
       if (id.eq.1) go to 860

c----------------------------------------------------------------
c find product nuclide level branchings for this target of tape 4
c----------------------------------------------------------------
  820  read (4,140,end=850) jdt,jl(i),jp(i)
  830  read (4,120) q(i),(el(il,i),il=1,jl(i))
       do 840 ip=1,jp(i)
  840  read (4,120) ep(ip,i),(f(il,ip,i),il=1,jl(i))
       if (idb(i)-jdt) 850,860,820
  850  write (6,150) idb(i)
       stop 3
  860  rewind 5
       if (nq.ne.0) go to 870
       sbtqan=0.
  870  read (5,20) icont
       if (icont.ne.0) go to 870
  880 continue

c------------------------------------------
c  Read material constituents for region C
c------------------------------------------
      read (1,15)title3
      write(11,2090)
      write(11,3050)title3
      read (1,*) nzc,isgc
      write(11,3202)nzc
      if (isgc.eq.0) then
      write(11,3212)isgc
      elseif (isgc.eq.1) then
       write(11,3222)isgc
      else
       stop 'Error in isgc input!'
      endif
      write(11,2090)
      write(11,2082)
      write(11,2090)
      write(11,2100)
      write(11,2110)

c---------------------------------
c read input material constituents
c---------------------------------
      do 1520 j=1,nzc
       read (1,*) jzc(j),azc(j)
 1520 write(11,2120)jzc(j),azc(j)

c----------------------------------------
c order the nz material constituents by z
c----------------------------------------
      do 1540 j=1,nzc
       isw=0
       do 1530 jn=2,nzc
        if (jzc(jn).gt.jzc(jn-1)) go to 1530
        jtem=jzc(jn-1)
        atem=azc(jn-1)
        jzc(jn-1)=jzc(jn)
        azc(jn-1)=azc(jn)
        jzc(jn)=jtem
        azc(jn)=atem
        isw=1
 1530  continue
      if (isw.eq.0) go to 1550
 1540 continue

c----------------------------------------------------
c read tape 2 for stopping cross-section coefficients
c----------------------------------------------------
 1550 rewind 2
 1560 read (2,20) icont
      if (icont.ne.0) go to 1560
      do 1620 j=1,nzc
 1570  read (2,40,end=1610) iz,c
 1580  if (jzc(j)-iz) 1610,1590,1570
 1590  continue
       do 1600 jt=1,9
 1600  czc(jt,j)=c(jt)
       if (isgc.eq.1.and.c(10).gt.0.) then
        do 1605 jt=1,5
 1605   czc(jt,j)=c(jt+9)
       endif
       go to 1620
 1610  write (6,50) jzc(j)
       stop 1
 1620 continue
 1625 if (id.eq.1) go to 1720

c----------------------------------------
c read target information from user input
c----------------------------------------
 1720 continue
      read (1,*) ntc
      write(11,2261)ntc
      write(11,2090)
      write(11,2281)
      write(11,2090)
      write(11,2230)
      write(11,2240)

c--------------------------
c loop on target nuclides i
c--------------------------
      rewind 3
      if (id.eq.2) rewind 4
 1730 read (3,20) icont
      if (icont.ne.0) go to 1730
 1740 if (id.eq.2) read (4,20) icont
      if (icont.ne.0) go to 1740
      do 1880 i=1,ntc

c------------------------------------------------------
c find alpha-n cross sections for this target on tape 3
c------------------------------------------------------
 1760  read (1,*) idc(i),atc(i)
       write(11,2250)idc(i),atc(i)
 1765  read (3,105,end=1780) jdt,jpsc(i)
 1770  read (3,120) (ec(ip,i),xc(ip,i),ip=1,jpsc(i))
       if (idc(i)-jdt) 1780,1790,1765
 1780  write (6,130) idc(i)
       stop 2
 1790  lztc(i)=idc(i)/10000
       latc(i)=idc(i)/10-lztc(i)*1000
       atar=latc(i)
       lstc(i)=idc(i)-lztc(i)*10000-latc(i)*10
       amtc(i)=' '
       if (lstc(i).ne.0) amtc(i)='m'
       aproc(i)=atar+3.
c       do 1800 ip=2,jpsc(i)
c        eamin=ec(ip-1,i)
c        if (xc(ip,i).ne.0.) go to 1810
c 1800  continue
 1810  if (eamin.lt.0.001) eamin=0.001
       if (id.eq.1) go to 1860

c----------------------------------------------------------------
c find product nuclide level branchings for this target of tape 4
c----------------------------------------------------------------
 1820  read (4,140,end=1850) jdt,jlc(i),jpc(i)
 1830  read (4,120) qc(i),(elc(il,i),il=1,jlc(i))
       do 1840 ip=1,jpc(i)
 1840  read (4,120) epc(ip,i),(fc(il,ip,i),il=1,jlc(i))
       if (idc(i)-jdt) 1850,1860,1820
 1850  write (6,150) idc(i)
       stop 3
 1860  rewind 5
       if (nq.ne.0) go to 1870
       sbtqan=0.
 1870  read (5,20) icont
       if (icont.ne.0) go to 1870
 1880 continue

c-------------------------
c  End Execution of INPUT
c-------------------------
      return
      end


c=======================================================================
c  Three Region Problem Processor Subroutine (1/99)
c=======================================================================

      subroutine processor(title,idd,id)

c----------
c  Storage
c----------
      character*7 title, title1, title2, title3
      character*2 jsm(105)
      character*1 amt(20), amtc(20)
      dimension title(11), r(4001), ps(4001), p(4001), nal(20),
     1 ee(4001), en(751), cg(4001), jza(20), aza(20),
     2 jq(300), faq(300), cza(9,20), fal(30,20), eala(30,20), 
     3 fdng(1000), jzb(20), azb(20), czb(9,20), idb(20), atb(20),
     4 jps(20), e(1100,20), x(1100,20), lzt(20), lat(20), lst(20),
     5 apro(20), jl(20), jp(20), q(20), el(20,20), ep(200,20),
     6 f(20,200,20), jzc(20), azc(20), czc(9,20), idc(20), atc(20),
     7 jpsc(20), ec(1100,20), xc(1100,20), lztc(20), latc(20), 
     8 lstc(20), aproc(20), jlc(20), jpc(20), qc(20), elc(20,20), 
     9 epc(200,20), fc(20,200,20), scxa(4001), scxb(4001), 
     1 astab(4001), eal(4001), alam(20), D(4001),
     2 itrans(4001,4001), astbc(4001), title1(11), title2(11), 
     3 title3(11)
      real nstabb(750), nstbcb(750), nstbcc(750), nst(750), 
     1 nstnorm(750)

       common /names/ jsm
       common /masses/ amass(105)
       common /vars/ jza, aza, jq, faq, cza, fal, eala, fdng, jzb, 
     1 azb, czb, idb, atb, jps, e, x, lzt, lat, lst, apro, jl, jp, 
     2 q, el, ep, f, jzc, azc, czc, idc, atc, jpsc, ec, xc, lztc, 
     3 latc, lstc, aproc, jlc, jpc, qc, elc, epc, fc, nza, isga, nq,
     5 nzb, isgb, anumb, thick, ntb, nzc, isgc, ntc, amt, amtc, nal,
     6 alam, title1, title2, title3
       common /grids/ ee, en, cg, nag, eamax, eamin, nng, 
     1 enmax, enmin, ncg, dea, eal

c--------------------
c  Format Statements
c--------------------
   16 format(//,10(8h========),//,'Title: ',11a7)
   80 format (1p8e10.3)
   81 format(31x,'Total (all groups): ',1pE10.3,' neutrons/sec-cm^2.',
     1/,27x,'Average Neutron Energy: ',1pE10.3,' MeV.')
  130 format (1x,'SOURCE DOES NOT HAVE ALPHA EMITTER ', i6)
  131 format(' WARNING-AN ALPHA ENERGY THAT EXCEEDED ',f7.3,' MEV WAS',
     1' SET TO:',f7.3,' MeV')
  390 format (////,' Total Neutron Spectrum')
 2021 format(/,80(1h-))
 3000 format(//,22x,'Summary of Output')
 3010 format(22x,'=================')
 3020 format('Total (alpha,n) neutron source from all sources and ',
     1'targets: ',1pE12.5,' n/sec-cm^2.')
 3030 format(/,'Average (alpha,n) neutron energy: ',1pE12.5,' MeV.')
 3050 format(/,'Normalized Neutron Energy Spectrum by Energy Group ',
     1'for All Sources:',//,10x,'Group',5x,'Contribution',/,10x,
     2'-----',5x,'------------')
 3051 format(11x,i3,7x,1pE12.5)
 3500 format(/,'Portion of Total Neutron Source Rate Accounted for ',
     1'in the Energy Spectrum: ',F5.1,'%.')
 3510 format(/,'Warning: Energy Group Structure May Be Poorly Chosen.')

c-------------------------------------------------------------------
c   Calculate stopping-power cross sections for regions A, B, and C
c-------------------------------------------------------------------
      call stopxs(nag, nza, jza, aza, ee, cza, scxa)
      call stopxs(nag, nzb, jzb, azb, ee, czb, scxb)
     
c-------------------------------------------------------
c   Calculate alpha source term at the interface for ab
c-------------------------------------------------------
      do 230 m=1,nag
       astab(m)=0.0
  230 continue

      fnag=nag
      do 250 k=1,nq
       r(1)=1./scxa(1)
       fact=alam(k)*faq(k)*2.5e+20
       do 240 m=1,nag
        ps(m)= 0.0
        r(m+1)= 1.0/scxa(m+1)
        p(m)=fact*(r(m+1)+r(m))*dea/2.
  240  continue
       do 246 mm = 1,nal(k)
        if(eala(mm,k).gt.eamax)write(6,131)eamax,eamax
        if(eala(mm,k).gt.9.8)eala(mm,k)=9.8                 !modified - vk
        do 245 m = 1,nag
         fm=1.
         if(eala(mm,k).ge.ee(m+1))go to 244
         if(eala(mm,k).gt.ee(m))fm=(eala(mm,k)-ee(m))/dea
         if(eala(mm,k).lt.ee(m))fm=0.0
  244    ps(m) = ps(m)+fal(mm,k)*p(m)*fm
  245   continue
  246  continue
       do 247 m=1,nag
        astab(m) = astab(m) + ps(m)
  247  continue
  250 continue

      knt=nag
      do 260 i=1,nag
       if(astab(knt).gt.0.0)go to 270
       knt=knt-1
  260 continue
  270 naga=knt
      if(naga.eq.0)write(6,130)
      if(naga.eq.0)stop 'All astab are zero'

c--------------------------------------
c   Calculate the ab transition energy
c--------------------------------------
      do 305 i=1,ncg
       if (cos(cg(i)).le.0.0) then
        stop 'Error in angular structure'
       endif
       D(i)=thick/cos(cg(i))
       do 304 ig=1,nag
        itrans(i,ig)=ig
        sum=(dea/2.)*((1./scxb(ig))+(1./scxb(ig+1)))/anumb
        if (sum.gt.D(i)) goto 304
  302   sum=0
        itrans(i,ig)=itrans(i,ig)+1
        if (itrans(i,ig).gt.nag) goto 304
        do 303 m=ig,itrans(i,ig)
         sum=sum+(dea/2.)*((1./scxb(m))+(1./scxb(m+1)))/anumb
  303   continue
        if (sum.lt.D(i)) goto 302
  304  continue
  305 continue

c-------------------------------------------------------
c   Calculate alpha source term at the interface for bc
c-------------------------------------------------------
      do 330 m=1,nag
       astbc(m)=0.0
  330 continue

      fnag=nag
      do 350 k=1,nq
       fact=alam(k)*faq(k)*1.25e+20
       do 348 i=1,ncg
        r(1)=1./scxa(1)
        do 340 m=1,nag
         ps(m)= 0.0
         r(m+1)= 1.0/scxa(m+1)
         p(m)=fact*(r(m+1)+r(m))*dea/2.
  340   continue
        do 346 mm = 1,nal(k)
         if(eala(mm,k).gt.eamax)write(6,131)eamax,eamax
         if(eala(mm,k).gt.9.8)eala(mm,k)=9.8             ! modified - vk
         do 345 m = 1,nag
          fm=1.
          if(eala(mm,k).ge.ee(m+1))go to 344
          if(eala(mm,k).gt.ee(m))fm=(eala(mm,k)-ee(m))/dea
          if(eala(mm,k).lt.ee(m))fm=0.0
  344     pm=1.
          if(eala(mm,k).lt.ee(itrans(i,m)))pm=0.0
          ps(m) = ps(m)+fal(mm,k)*p(m)*fm*pm
  345    continue
  346   continue
        do 347 m=1,nag
         astbc(m)=astbc(m)+ps(m)*(cos(2.*cg(i))-cos(2.*cg(i+1)))
  347   continue
  348  continue
  350 continue

      knt=nag
      do 360 i=1,nag
       if(astab(knt).gt.0.0)go to 370
       knt=knt-1
  360 continue
  370 nagb=knt
      if(nagb.eq.0)write(6,130)
      if(nagb.eq.0)stop 'All astbc are zero'
      do 380 i=1,nag
       eal(i)=(ee(i+1) + ee(i))/2.
  380 continue

c--------------------------------------------------------------
c   Calculate neutron source from ab interface due to region B
c--------------------------------------------------------------
      if (ntb.ne.0) then
       title1(1)='Alphas '
       title1(2)='at inte' 
       title1(3)='rface a'
       title1(4)='b using'
       title1(5)=' region'
       title1(6)=' B mate'
       title1(7)='rials f'
       title1(8)='or neut'
       title1(9)='ron pro'
       title1(10)='duction'
       call neutron(title1,ntb,idb,atb,jp,ep,f,lzt,lat,lst,amt,el,
     1 apro,q,e,x,nzb,jzb,azb,czb,jps,jl,astab,nstabb)

c--------------------------------------------------------------
c   Calculate neutron source from bc interface due to region B
c--------------------------------------------------------------
       title1(1)='Alphas '
       title1(2)='at inte' 
       title1(3)='rface b'
       title1(4)='c using'
       title1(5)=' region'
       title1(6)=' B mate'
       title1(7)='rials f'
       title1(8)='or neut'
       title1(9)='ron pro'
       title1(10)='duction'
       call neutron(title1,ntb,idb,atb,jp,ep,f,lzt,lat,lst,amt,el,
     1 apro,q,e,x,nzb,jzb,azb,czb,jps,jl,astbc,nstbcb)
      endif

c--------------------------------------------------------------
c   Calculate neutron source from bc interface due to region C
c--------------------------------------------------------------
      title1(1)='Alphas '
      title1(2)='at inte' 
      title1(3)='rface b'
      title1(4)='c using'
      title1(5)=' region'
      title1(6)=' C mate'
      title1(7)='rials f'
      title1(8)='or neut'
      title1(9)='ron pro'
      title1(10)='duction'
      call neutron(title1,ntc,idc,atc,jpc,epc,fc,lztc,latc,lstc,amtc,
     1elc,aproc,qc,ec,xc,nzc,jzc,azc,czc,jpsc,jlc,astbc,nstbcc)

c----------------------------------
c   Calculate total neutron source
c----------------------------------
      tot=0
      title1(1)='Total n'
      title1(2)='eutron ' 
      title1(3)='product'
      title1(4)='ion fro'
      title1(5)='m all i'
      title1(6)='nterfac'
      title1(7)='es     '
      title1(8)='       '
      title1(9)='       '
      title1(10)='       '
      write(7,16)title1
      write(8,16)title1
      ebaran=0.0
      do 500 m=1,nng
       nst(m)=nstbcc(m)+nstabb(m)-nstbcb(m)
       tot=tot+nst(m)
       ebaran=ebaran+((en(m)+en(m+1))*nst(m)/2.0)
  500 continue
      ebaran=ebaran/tot
      write(7,*)
      write(7,2021)
      write(7,390)
      write(7,80) (nst(n),n=1,nng)
      write(7,81) tot,ebaran
      gtmga=0.0
      do 600 n=1,nng
       nstnorm(n)=nst(n)/tot
       gtmga=gtmga+nstnorm(n)
  600 continue
      write(8,*)
      write(8,2021)
      write(8,390)
      write(8,80) (nstnorm(n),n=1,nng)
      write(8,81) gtmga,ebaran

c------------------
c   Output Summary
c------------------
      write(11,*)
      write(11,*)
      write(11,3000)
      write(11,3010)
      write(11,*)
      write(11,3020)tot
      if (id.gt.1) then
       write(11,3030)ebaran
       write(11,3050)
       do 1450 n=1,nng
        write(11,3051)n,nstnorm(n)
 1450  continue
      endif

c------------------------------------
c   Terminate execution of PROCESSOR
c------------------------------------
      return
      end

c=======================================================================
c  Three Region Problem Stopping Cross Section Calculator (1/99)
c=======================================================================

      subroutine stopxs(nag,nz,jzm,azm,ee,czm,scx)

c----------
c  Storage
c----------
      character*2 jsm(105)
      dimension ee(4001), jzm(20), azm(20), czm(9,20), scx(4001)

      common /names/ jsm
      common /masses/ amass(105)

c-----------------
c  Set parameters
c-----------------
      zalp=2.
      alph=4.

c---------------------------------------------------------------
c  For each energy ee(m) calc. and store stopping cross section
c---------------------------------------------------------------
      do 930 m=1,nag+1
       scx(m)=0.
       do 920 j=1,nz
        zmat=jzm(j)
        amat=amass(jzm(j))
c-----  nuclear stopping---
        term=zalp**0.66667 + zmat**0.66667
        sterm=sqrt(term)
        bot=zalp*zmat*(alph+amat)*sterm
        rep=32530.*amat*ee(m)/bot
        if(rep.lt.0.001) then
         dcx=1.593*sqrt(rep)
         go to 905
        endif
        if(rep.lt.10.) then
         bot=1.+6.8*rep+3.4*rep**1.5
         dcx=1.7*sqrt(rep)*alog(rep+2.71828)/bot
         go to 905
        endif
        dcx=alog(0.47*rep)/(2.*rep)
c-----  electronic stopping---
  905   if(ee(m).gt.30.) go to 906
        slow=czm(1,j)*(1000.*ee(m))**czm(2,j)
        shigh=(czm(3,j)/ee(m))
        shigh=shigh*alog(1.+czm(4,j)/ee(m)+czm(5,j)*ee(m))
        dcx=dcx+slow*shigh/(slow+shigh)
        go to 907
  906   eil=alog(1/ee(m))
        arg=czm(6,j)+eil*(czm(7,j)+czm(8,j)*eil+czm(9,j)*eil*eil)
        dcx=dcx+exp(arg)
  907  continue
  920  scx(m)=scx(m)+azm(j)*dcx
  930 continue

c--------------------------------
c  Terminate execution of STOPXS
c--------------------------------
      return
      end



c=======================================================================
c  Three Region Problem Neutron Source Subroutine (1/99)
c=======================================================================

      subroutine neutron(title,nt,idt,at,jp,ep,f,lzt,lat,lst,amt,el,
     1 apro,q,e,x,nz,jz,az,cz,jps,jl,ast,gtsan)

c----------
c  Storage
c----------
      character*7 title
      character*2 jsm(105)
      character*1 amt(20)
      dimension ee(4001), en(751), cg(4001), e(1100,20), el(20,20),
     1 x(1100,20), idt(20), at(20), scx(4001), jps(20), jl(20), 
     2 cx(4001), eal(4001), apro(20), q(20), jp(20), ep(200,20),
     3 f(20,200,20), lat(20), ast(4001), san(750), ts(750), lzt(20),
     4 r(4001), rr(4001), p(4001), tsan(750), gtsan(750), s(750), 
     5 totlev(20), sl(750,20), etopl(20), ebotl(20), jz(20), az(20),
     6 cz(9,20), gtnorm(750), title(11)

       common /names/ jsm
       common /masses/ amass(105)
       common /grids/ ee, en, cg, nag, eamax, eamin, nng, 
     1 enmax, enmin, ncg, dea, eal

c--------------------
c formats statements
c--------------------
   10 format (11a7)
   15 format(//,11a7,/)
   16 format(//,10(8h========),//,'Title: ',11a7)
   20 format (3i3)
   30 format (i8,e12.5)
   40 format(i3,8e8.3,f10.3,/,3x,5f8.3)
   50 format (69h failed to find stopping cross section coefficients on
     1tape 2 for iz=,i3,5h stop)
   60 format (6e12.5)
   70 format(/,1x,'Neutron Multigroup Structure (MeV)')
   79 format(/,1x,'Neutron Spectrum (neuts/cm^2-sec)')
   80 format (1p8e10.3)
   81 format(31x,'Total (all groups): ',1pE10.3,' neutrons/sec-cm^2.',
     1/,27x,'Average Neutron Energy: ',1pE10.3,' MeV.')
   90 format (i3,5e12.5)
   95 format(2i6)
  100 format(//,1h1,30x,'Table I',/,31x,7(1h=),//,12x,'(alpha,n) ',
     1'Neutron Production by Target and Alpha Energy',///,17x,
     2'target            alpha   alphas/sec     p(e)     neuts/sec',
     3/,7x,'target  atom frac.          energy     /cm^2',4x,
     4'neut/alpha',4x,'/cm^2',/,1h+,6x,6(1h_),2x,10(1h_),2x,6(1h_),2x,
     5 6(1h_),3(12h  __________))
c 102 format(1h+,77x,37hfractional ground/excited level split,/,
c    1 1h+,77x,49(1h_))
  105 format(i8,i4)
  120 format (8e10.3)
  130 format (52h failed to find alpha-n cross section on tape 3 for ,i8
     1 ,5h stop)
  140 format (i8,2i3)
  150 format (66h failed to find comp.nucl.neut.emiss.lev. branching on
     1tape 4 for ,i8,5h stop)
  170 format (/,17h alpha energy of ,1pe12.5,62h exceeds current 6.5 MeV
     1 limit of most cross sections, stop.  )
  180 format (7h target,i8,40h cross-section data not available at ee(
     1 ,i4,3h)= ,1pe12.5,6h stop.)
  200 format(33x,f8.3,1p3e12.4)
c 205 format(1h+,(76x,7f7.5,/,1x))
  210 format(7x,a2,i3,a1,1pe12.4,8x,         0pf8.3,1p3e12.4)
  220 format (14h alpha energy ,1pe12.5,45h exceeds range of level branc
     1hing data, stop.)
  240 format(/,' (a,n) neutrons from alphas on ',a2,i3,a1,' in target.')
  250 format(1h+,66x,10(1h_),/,55x,'Total:',4x,1pe12.4,/)
  251 format(1h+,66x,10(1h_),/,41x,'Total (this target):',
     14x,1pe12.4,/)
  252 format (1h+,66x,10(1h_),/,41x,'Total (all targets):',
     14x,1pe12.4,/)
  255 format(1h )
  260 format (///,' Total (alpha,n) neutron spectrum this target')
  270 format(///,' Grand total (alpha,n) neutron spectrum, all targets,'
     1 ,1x,'all sources')
  390 format (////,' Total Neutron Spectrum')

 2021 format(/,80(1h-))
 2030 format('Target title: ', 11a7)
 2050 format('Number of elemental constituents on target side:',i3)
 2060 format('Solid stopping cross-sections used (isg=',i2,') on',
     1 1x,'target side.')
 2070 format('Gas stopping cross-sections used (isg=',i2,')',
     1 1x,'on target side.')
 2082 format('Elemental Constituents on Target Side:')
 2090 format(' ')
 2100 format('  Z-value  Atom Fraction')
 2110 format('  -------  -------------')
 2120 format(4x,i3,5x,f12.10)
 2130 format('Number of neutron spectrum energy groups:', i4)
 2140 format('Maximum neutron energy is ',1pE12.5,' MeV.')
 2150 format('Minimum neutron energy is ',1pE12.5,' MeV.')
 2160 format('Energy Group Structure:')
 2170 format('  Group  Upper-Bound  Lower-Bound')
 2180 format('  -----  -----------  -----------')
 2190 format(2x,i4,4x,1pE12.5,3x,1pE12.5)
 2193 format(/,'Warning: First Group Upper-Bound Does Not Coincide ',
     1'With Maximum Neutron Energy.')
 2195	format(/,'Warning: Energy Group Structure Not Entered in ',
     1'Decreasing Order.',/,5x,'Group ',i3,' Upper-Bound is Greater ',
     2'Than Group ',i3,' Upper-Bound.') 
 2260 format(/,'Number of target nuclides to be used:',i3)
 2270 format(i5,' alpha energy groups used in target calculation.')
 2280 format('Target Nuclides:')
 2290 format('    ZAID     Atom Fraction')
 2300 format('    ----     -------------')
 2310 format(3x,i6,4x,1pE12.5)

 3000 format(//,22x,'Summary of Output')
 3010 format(22x,'=================')
 3020 format('Total (alpha,n) neutron source from all sources and ',
     1'targets: ',1pE12.5,' n/sec-cm^2.')
 3030 format(/,'Average (alpha,n) neutron energy: ',1pE12.5,' MeV.')
 3050 format(/,'Normalized Neutron Energy Spectrum by Energy Group ',
     1'for All Sources:',//,10x,'Group',5x,'Contribution',/,10x,
     2'-----',5x,'------------')
 3051 format(11x,i3,7x,1pE12.5)
 3500 format(/,'Portion of Total Neutron Source Rate Accounted for ',
     1'in the Energy Spectrum: ',F5.1,'%.')
 3510 format(/,'Warning: Energy Group Structure May Be Poorly Chosen.')

c---------------
c write headers
c---------------
      write(6,16)title
      write(7,16)title
      write(8,16)title

c--------------------------
c loop on target nuclides i
c--------------------------
      alph=4.
      aneut=1.
      do 1120 i=1,nt
       etpant=0.
       ebtant=0.
       totqan=0.
       if(id.eq.1) go to 790
       do 650 n=1,nng
        tsan(n)=0.
  650  continue
       do 700 ip=2,jps(i)
        eamin=e(ip-1,i)
        if (x(ip,i).ne.0.) go to 710
  700  continue
  710  if (eamin.lt.0.001) eamin=0.001
       etop=0.
       ebot=0.
       do 775 il=1,jl(i)
        etopl(il)=0.
        ebotl(il)=0.
        do 774 n=1,nng
         sl(n,il)=0.
  774   continue
        totlev(il)=0.
  775  continue
       sbtqan=0.
       if(id.eq.1) go to 790
       do 780 n=1,nng
        san(n)=0.
  780  continue
  790  continue

c---------------------------------------------------------
c set alpha energy grid ee(m) for this target/source combo
c---------------------------------------------------------
      eamax=eal(nag)
  850 continue
      if(eamax.le.eamin) then
      qan=0.
      mm=0
      p(nag+1)=0.
      pval=0.
      go to 945
      endif
      if (eamax.le.9.8) go to 860                    !modified - vk
      write (6,170) eamax
      stop 5
  860 fnag=nag
      nagp1=nag+1
      dea=(eamax-eamin)/fnag
      ee(1)=eamin
      do 870 m=2,nag
       fm=m-1
  870 ee(m)=ee(1)+dea*fm
      ee(nagp1)=eamax

c----------------------------------------------------
c for each energy ee(m) calc. and store cross section
c----------------------------------------------------
       ip=1
       do 900 m=1,nag+1
  880   if (ee(m).ge.e(ip,i).and.ee(m).le.e(ip+1,i)) go to 890
        ip=ip+1
        if (ip.lt.jps(i)) go to 880
        write (6,180) idt(i),m,ee(m)
        stop 6
  890   slope=(x(ip+1,i)-x(ip,i))/(e(ip+1,i)-e(ip,i))
        enrsep=(x(ip,i)*e(ip+1,i)-x(ip+1,i)*e(ip,i))/
     1   (e(ip+1,i)-e(ip,i))
        cx(m)=enrsep+ee(m)*slope
  900  continue

c-------------------------------------------------------------
c for each energy ee(m) calc. and store stopping cross section
c-------------------------------------------------------------
      call stopxs(nag, nz, jz, az, ee, cz, scx)

c-----------------------------------------------------------
c for each energy ee(m) calc. ratio r(m), then integral p(m)
c-----------------------------------------------------------
       r(1)=cx(1)/scx(1)
       p(1)=0.
       fact=1.0e-06*at(i)
       do 940 m=2,nag+1
        r(m)=cx(m)/scx(m)
        p(m)=p(m-1)+fact*(r(m-1)+r(m))*dea/2.
  940  continue
  945  continue

c----------------------------------------------
c  calculate neutrons from alphas
c----------------------------------------------
  950 do 1080 l=1,nag
      do 960 m=1,nag
       mm=m
  960 if (eal(l).ge.ee(m).and.eal(l).le.ee(m+1)) go to 970
      mm=0
      pval=0.
      go to 975
  970 pval=p(mm)+(eal(l)-ee(mm))*(p(mm+1)-p(mm))/(ee(mm+1)-ee(mm))
  975 aps=ast(l)
      qan=aps*pval
      sbtqan=sbtqan+qan
      totqan=totqan+qan
      if (l.eq.1) go to 980
      write (6,200) eal(l),aps,pval,qan
      if(nag.ne.1.and.l.eq.nag) write(6,250) sbtqan
      go to 990
  980 write(6,210)jsm(lzt(i)),lat(i),amt(i),at(i),eal(l),aps,pval,qan
      if(nag.eq.1) write(6,255)
  990 if (id.eq.1.or.mm.eq.0) go to 1080

c-----------------------------------------------------------------
c calculate (alpha,n) neutron spectrum contribution in multigroups
c-----------------------------------------------------------------
 1000 mmm1=mm-1
      do 1010 m=1,mmm1
 1010 rr(m)=(p(m+1)-p(m))/pval
      rr(mm)=(pval-p(mm))/pval
      do 1020 n=1,nng
 1020 s(n)=0.
      do 1060 il=1,jl(i)
      qlev=q(i)-el(il,i)
      e90=-qlev*apro(i)/(apro(i)-alph)
      thre=-qlev*(aneut+apro(i))/(aneut+apro(i)-alph)
      if (qlev.gt.0.) thre=0.
      do 1060 m=1,mm
      ea=(ee(m)+ee(m+1))/2.
      if (m.eq.mm) ea=(eal(l)+ee(mm))/2.
      if (ea.lt.thre) go to 1055
      do 1030 ip=2,jp(i)
      lp=ip
 1030 if (ep(ip-1,i).le.ea.and.ep(ip,i).ge.ea) go to 1040
      write (6,220)ea
      stop 7
 1040 lp1=lp-1
      bx=f(il,lp1,i)+(ea-ep(lp1,i))*(f(il,lp,i)-f(il,lp1,i))/
     1 (ep(lp,i)-ep(lp1,i))
      term1=sqrt(4.*alph*aneut*ea)/(2.*(aneut+apro(i)))
      term2=alph*aneut*ea/((aneut+apro(i))*(aneut+apro(i)))
      term3=(apro(i)*ea+apro(i)*qlev-alph*ea)/(aneut+apro(i))
      senmax=term1+sqrt(term2+term3)
      enmax=senmax*senmax
      if (ea.le.e90) senmin=term1-sqrt(term2+term3)
      if (ea.gt.e90) senmin=-term1+sqrt(term2+term3)
      enmin=senmin*senmin
      ebar=(enmin+enmax)/2.
      val=qan*rr(m)*bx
      etop=etop+val*ebar
      etopl(il)=etopl(il)+val*ebar
      ebotl(il)=ebotl(il)+val
      ebot=ebot+val
      dele=enmax-enmin
      do 1050 n=1,nng
      if (en(n).lt.enmin) go to 1050
      if (en(n+1).gt.enmax) go to 1050
      de=en(n)-en(n+1)
      if (en(n+1).lt.enmin) de=de-(enmin-en(n+1))
      if (en(n).gt.enmax) de=de-(en(n)-enmax)
      gpadd=qan*rr(m)*bx*de/dele
      s(n)=s(n)+gpadd
      sl(n,il)=sl(n,il)+gpadd
      totlev(il)=totlev(il)+gpadd
 1050 continue
 1055 continue
 1060 continue
      do 1070 n=1,nng
      san(n)=san(n)+s(n)
      gtsan(n)=gtsan(n)+s(n)
      ts(n)=ts(n)+s(n)
 1070 tsan(n)=tsan(n)+s(n)
 1075 continue
 1080 continue
      if(id.eq.1) go to 1200

c----------------------------------------------------
c output this (alpha,n) neutron spectrum contribution
c for this source/target combination
c----------------------------------------------------
      ebarqt=0.
      if(ebot.le.0.) go to 1085
      ebarqt=etop/ebot
      etpant=etpant+etop
      ebtant=ebtant+ebot
 1085 continue
      smga=0.
      do 1090 n=1,nng
 1090 smga=smga+san(n)
c     fracgp=smga/sbtqan
c     do 1095 il=1,jl(i)
c1095 frclev(il)=totlev(il)/sbtqan
c     write(6,205) (frclev(il),il=1,jl(i))
      if (nq.ne.0) write (7,240) jsm(lzt(i)),lat(i),amt(i)
      write (7,79)
      write (7,80) (san(n),n=1,nng)
      write (7,81) smga,ebarqt
      smga=0.
      if (nq.ne.0) write (8,240) jsm(lzt(i)),lat(i),amt(i)
      do 1097 n=1,nng
      san(n)=san(n)/sbtqan
 1097 smga=smga+san(n)
      write (8,80) (san(n),n=1,nng)
      write (8,81) smga, ebarqt
      do 1098 il=1,jl(i)
      if(ebotl(il).le.0.) go to 1098
c     ebarl=etopl(il)/ebotl(il)
 1098 continue
 1101 gttqan=gttqan+totqan

c--------------------------------------
c output total (alpha,n) source for all
c alpha sources on this target nuclide
c--------------------------------------
      if(id.eq.1) go to 1120
      tmga=0.
      do 1110 n=1,nng
 1110 tmga=tmga+tsan(n)
c     fracgp=tmga/totqan
      ebart=etpant/ebtant
      etopan=etopan+etpant
      ebotan=ebotan+ebtant
      if(nq.le.1) go to 1120
      write(7,260)
      write(7,80) (tsan(n),n=1,nng)
      write(7,81) tmga,ebart
      write(7,2090)
      write(7,2021)
      write(7,2090)
      tmga=0.
      do 1115 n=1,nng
       tsan(n)=tsan(n)/totqan
 1115 tmga=tmga+tsan(n)
      write(8,260)
      write(8,80) (tsan(n),n=1,nng)
      write(8,81) tmga,ebart
      write(8,2090)
      write(8,2021)
      write(8,2090)
      write(10,2090)
      write(10,2021)
      write(10,2090)
 1120 continue

c----------------------------------------------------------------------
c output grand total (alpha,n) neutron source: all targets, all sources
c----------------------------------------------------------------------
      if (nt.gt.1) write (6,252) gttqan
      if(id.eq.1) go to 1200
      gtmga=0.
      do 1130 n=1,nng
 1130 gtmga=gtmga+gtsan(n)
      ebaran=etopan/ebotan
c     fracgp=gtmga/gttqan
      write(7,2090)
      write(7,2021)
      write(7,390)
      write(7,80) (gtsan(n),n=1,nng)
      write(7,81) gtmga,ebaran
      gtmga=0.
      dummy=0.0
      do 1135 n=1,nng
       dummy=dummy+gtsan(n)
       gtnorm(n)=gtsan(n)/gttqan
 1135 gtmga=gtmga+gtnorm(n)
      write(8,2090)
      write(8,2021)
      write(8,390)
      write(8,80) (gtnorm(n),n=1,nng)
      write(8,81) gtmga,ebaran
 1200 continue

c------------------------------
c  End of Execution of NEUTRON
c------------------------------
      return
      end

c=======================================================================
c  Block Data Statement (7/97)
c=======================================================================

      block data initvl
      character*2 jsm(105)
      dimension amass(105)
      common /names/ jsm
      common /masses/ amass

      data jsm /' h','he','li','be',' b',' c',' n',' o',' f','ne','na',
     1 'mg','al','si',' p',' s','cl','ar',' k','ca','sc','ti',' v','cr',
     2 'mn','fe','co','ni','cu','zn','ga','ge','as','se','br','kr','rb',
     3 'sr',' y','zr','nb','mo','tc','ru','rh','pd','ag','cd','in','sn',
     4 'sb','te',' i','xe','cs','ba','la','ce','pr','nd','pm','sm','eu',
     5 'gd','tb','dy','ho','er','tm','yb','lu','hf','ta',' w','re','os',
     6 'ir','pt','au','hg','tl','pb','bi','po','at','rn','fr','ra','ac',
     7 'th','pa',' u','np','pu','am','cm','bk','cf','es','fm','md','no',
     8 'lr','ru','ha'/

      data amass /1.0079,4.0026,6.941,9.012,10.81,12.011,14.007,16.,19.,
     120.18,22.99,24.31,26.98,28.09,30.97,32.06,35.45,39.95,39.10,40.08,
     244.96,47.9,50.94,52.,54.94,55.85,58.93,58.7,63.55,65.38,69.72,
     372.59,74.92,78.96,79.9,83.8,85.47,87.62,88.91,91.22,92.91,95.94,
     499.,101.07,102.91,106.4,107.87,112.41,114.82,118.69,121.75,127.6,
     5126.90,131.3,132.9,137.33,138.91,140.1,140.9,144.2,
     6147.,150.4,151.96,157.25,158.9,162.5,164.9,167.3,168.9,173.,175.,
     7178.,181.,184.,186.,190.,192.,195.,196.,201.,204.,207.,209.,209.,
     8217.,222.,223.,226.,227.,232.,231.,238.,238.,239.,241.,243.,249.,
     9249.,252.,257.,258.,259.,260.,261.,262./

      end


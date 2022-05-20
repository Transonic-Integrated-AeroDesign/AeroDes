      program classiceq
      implicit none
      integer lxx,ndatx,nal,itx,i,ndat,ipolar,km,k,kp,m,it
      integer image,ittd,nalmin,nalmax,n,incidence,imarg,kx,inewton
      parameter(lxx=101,ndatx=21)
      real pi,eps,degrad,radeg,us3,omega,epser,rho,amu,xcg
      real statmarg,amlb,bm,cxm,cam,am,xacm,dm,em,atw,arm,rf,bt
      real cxt,cat,at,xact,flapcent,dCltda0,arteff,eteff
      real awt,art,dt,axflap,tetaflap,dCltdtt,dCmqdtt,rav,ruh
      real Cdbrake,aref,Ueq,Ref,hpower,Pcent,rcor,pcor,hpowatt
      real T,dTdv,tr0,dtrdv,prod,dcz,dcxm,dcxp,dCldam0
      real Clm00,dCmacdam0,Cmacm00,dCmdam0,Cmm00,dCldat0,Clt00
      real dCmacdat0,Cmact00,dCmdat0,Cmt00,dClda,Cl0,dCmacda,Cmac0
      real dCmda,Cm0,xac,dCldat,Cltt0,dCmdat,Cmact0
      real amarg,aleq,Cleq,Cweq,reyeq,dCmacdam,Cmacm0
      real dynaref,Cteq,Cdeq,Cdmeq,Clmeq,Cdim,Cdm0,reyf,Cdf0
      real Cdfeq,Clteq,Cdit,reyt,Cdt0,Cdteq,beteq,Cmac,Cmeq
      real dCldam,Clm0,dCmdam,Cmm0,Cmmeq,dcxdcz,af,aleqd,alph
      real d2cxdcz2,daleq,dbeteq,beteqd,theq,theqd,alphd,dum
      real ClCdeq,reseq,weq,winglift,dyn,taillift,tt,ttd,Clmnal
      real Cmmnal,Cdmnal,Clm,Cdm,Cmm,Clt,Cdt,Cmt,Cl,Cm,Cd,Cl3,Cd2
      real ratio1,ratio2,Cmcg,Cmcgeq,Clt0,Cmt0,dalph,dUeq
      real reym,aleq0,beteq0,Cleq0,Cdeq0,Ueq0,reseq0,xcpm,xcpt
      real zeng,dna,lna,tfd,dClmdtf,dCmmdtf,dCdtf0,dCdtf1,dCdtf2
      real ar,an,tf,reyr,Cdr0,Cdreq,reyn,Cdn0,Cdneq
      real vr(ndatx),tr(ndatx)
      real cx(lxx),cz(lxx),cq(lxx),inc(lxx)
      integer kxtrm(lxx)
      real incd,lref,lf,mass,mg
      doubleprecision det,usdet,b1,b2,b3
      doubleprecision aa(3,3),bb(3)
      character*4 title(18),prop(18),nc(5),yw(5)
      open(unit=12,file='classicpolar.dat',form='formatted')
      open(unit=13,file='classicpolar.prandtline',form='formatted')
      open(unit=14,file='classicpolar.cdclcq',form='formatted')
      open(unit=15,file='classiceq.data',form='formatted')
      open(unit=16,file='classiceq.cdclcqm',form='formatted')
      open(unit=17,file='classiceq.cdclcqt',form='formatted')
      open(unit=18,file='classiceq.cdclcq',form='formatted')
      open(unit=19,file='classiceq.cdmcln',form='formatted')
      open(unit=20,file='classiceq.cmcg',form='formatted')
      open(unit=21,file='classiceq.xtabvcl',form='formatted')
      open(unit=22,file='classiceq.cdclcqmeq',form='formatted')
      open(unit=23,file='classiceq.cdclcqcmcgeq',form='formatted')
      open(unit=24,file='classiceq.list',form='formatted')
      open(unit=25,file='classiceq.stab',form='formatted')
      open(unit=26,file='classiceq.itres',form='formatted')
      data nc/'!not',' con','verg','ed!!','!!!!'/
      data yw/'!ale','q ou','t of',' ran','ge!!'/
c*****constants
      pi=2.*asin(1.)
      eps=1.e-5
      degrad=pi/180.
      radeg=1./degrad
      nal=30
      dalph=degrad
      us3=1./3.
      write(6,*)
c*****Newton method parameters
      read(15,*)itx
      read(15,*)omega
      read(15,*)epser
c*****airflow and weight parameters
      read(15,*)rho
      read(15,*)amu
      read(15,*)mass
      read(15,*)xcg
      read(15,*)statmarg
      amlb=2.205*mass
      mg=mass*9.81
c*****data for main wing
      read(15,*)bm
      read(15,*)cxm
      read(15,*)cam
      read(15,*)am
      read(15,*)xacm
      read(15,*)dm
      read(15,*)em
      read(15,*)atw
      read(15,*)tfd
      read(15,*)dClmdtf
      read(15,*)dCmmdtf
      read(15,*)dCdtf0
      read(15,*)dCdtf1
      read(15,*)dCdtf2
      arm=bm*bm/am
      tf=degrad*tfd
c*****data for fuselage
      read(15,*)lf
      read(15,*)rf
      af=2.*pi*rf*lf
c*****data for tail
      read(15,*)bt
      read(15,*)cxt
      read(15,*)cat
      read(15,*)at
      read(15,*)xact
      read(15,*)ttd
      read(15,*)flapcent
      read(15,*)dCltda0
      read(15,*)arteff
      read(15,*)eteff
      read(15,*)dt
      read(15,*)awt
      art=bt*bt/at
      tt=degrad*ttd
      dt=dt*cxt
      axflap=1.-flapcent/100.
      tetaflap=acos(1.-2.*axflap)
      dCltdtt=(pi-tetaflap)/pi
      dCmqdtt=-sin(tetaflap)*axflap
c*****rudder and airbrake data
      read(15,*)Cdbrake
      read(15,*)rav
      read(15,*)ruh
      ar=rav*ruh
c*****reference parameters
      aref=am+at
      lref=lf
      Ueq=100.0
      Ref=rho*Ueq*cam/amu
c*****nacelle data
      read(15,*)zeng
      read(15,*)dna
      read(15,*)lna
c*****engine data
      read(15,*)hpower
      read(15,*)Pcent
      read(15,*)ndat
      read(15,1000)prop
      read(15,1000)title
      do 1 i=1,ndat
         read(15,*)vr(i),tr(i)
 1    continue
      an=pi*dna*lna
      hpowatt=745.7*hpower
      Ueq=0.
      call thrust(ndatx,ndat,Ueq,vr,tr,T,dTdv)
c*****thrust and thrust slope correction with density
      rcor=(rho/1.225)**us3
      pcor=(Pcent/100.)**us3
      tr0=rcor*pcor**2*(T+dTdv*Ueq)
      dtrdv=rcor*dTdv
c*****write data
      write(6,1002)
      write(6,1002)
      write(6,*)'*************************'
      write(6,*)'Newton method parameters:'
      write(6,*)'                         itx=',itx
      write(6,*)'                       omega=',omega
      write(6,*)'                       epser=',epser
      write(6,*)'***************'
      write(6,*)'air parameters:      density=',rho,'     (kg/m**3)'
      write(6,*)'               dynamic visc.=',amu,' (N.s/m**2)'
      write(6,*)'*************'
      write(6,*)'gravity data:           mass=',mass,' (kg) = '
     &     ,amlb,' (lb)'
      write(6,*)'              location of CG=',xcg,' (m)'
      write(6,*)'               static margin=',statmarg,' (if <0,'
     &     ,' not used to place xcg)'
      write(6,*)'***************'
      write(6,*)'main wing data:         span=',bm,' (m)'
      write(6,*)'  maximum wing or fuse chord=',cxm,' (m)'
      write(6,*)'  average wing or fuse chord=',cam,' (m)'
      write(6,*)'    wing+fuse projected area=',am,' (m**2)'
      write(6,*)'   location of w+f a.c. xacm=',xacm,' (m)'
      write(6,*)'                 wing AR arm=',arm
      write(6,*)'        wing relative camber=',dm
      write(6,*)'        wing+fuse efficiency=',em
      write(6,*)'downwash of tail on wing atw=',atw
      write(6,*)'          flap setting angle=',tfd,' (deg)'
      write(6,*)'flap setting influence on CL=',dClmdtf
      write(6,*)'flap setting influence on CM=',dCmmdtf
      write(6,*)'flap setting influence on CD=',dCdtf0
      write(6,*)'flap setting influence on CD=',dCdtf1
      write(6,*)'flap setting influence on CD=',dCdtf2
      write(6,*)'**************'
      write(6,*)'fuselage data:        length=',lf,' (m)'
      write(6,*)'                      radius=',rf,' (m)'
      write(6,*)'               fuselage area=',af,' (m**2)'
      write(6,*)'*********************'
      write(6,*)'tail/stabilizer data:   span=',bt,' (m)'
      write(6,*)'               maximum chord=',cxt,' (m)'
      write(6,*)'               average chord=',cat,' (m)'
      write(6,*)'                tail area at=',at,' (m**2)'
      write(6,*)'        tail relative camber=',dt
      write(6,*)' location of tail a.c. xact==',xact,' (m)'
      write(6,*)'          tail setting angle=',tt,' (rd) ='
     &     ,ttd,' (deg)'
      write(6,*)'     flap percentage of tail=',flapcent,' (%)'
      write(6,*)'given taillift slope dCltda0=',dCltda0
      write(6,*)' effective tail aspect ratio=',arteff
      write(6,*)'             tail efficiency=',eteff
      write(6,*)'    downwash of wing on tail=',awt
      write(6,*)'*****************************'
      write(6,*)'           airbrake: CDbrake=',Cdbrake
      write(6,*)'***********'
      write(6,*)'rudder data:'
      write(6,*)'    rudder average chord rav=',rav,' (m)'
      write(6,*)'           rudder height ruh=',ruh,' (m)'
      write(6,*)'              rudder area ar=',ar,' (m**2)'
      write(6,*)'***************'
      write(6,*)'reference data:  ref. length=',lref,' (m)'
      write(6,*)'                   ref. area=',aref,' (m**2)'
      write(6,*)'   reference Reynolds number=',Ref
      write(6,*)'************'
      write(6,*)' propulsion system reference=',prop
      write(6,*)'      z-coordinate of engine=',zeng,' (m)'
      write(6,*)'         diameter of nacelle=',dna,' (m)'
      write(6,*)'           length of nacelle=',lna,' (m)'
      write(6,*)'             nacelle area an=',an,' (m**2)'
      write(6,*)'             static traction=',tr0,' (N)'
      write(6,*)'              traction slope=',dtrdv,' (Ns/m)'
      write(6,*)'                engine power=',hpowatt,' (W) ='
     &     ,hpower,' (hp)'
      write(6,*)'percent of engine power avai=',Pcent
      write(6,*)'                        ndat=',ndat
      write(6,*)'    i=     vr(i)=      tr(i)='
      write(6,1004)(i,vr(i),tr(i),i=1,ndat)
      write(6,*)
      write(24,*)'*************************'
      write(24,*)'Newton method parameters:'
      write(24,*)'                         itx=',itx
      write(24,*)'                       omega=',omega
      write(24,*)'                       epser=',epser
      write(24,*)'***************'
      write(24,*)'air parameters:      density=',rho,'     (kg/m**3)'
      write(24,*)'                   kin.visc.=',amu,' (N.s/m**2)'
      write(24,*)'*************'
      write(24,*)'gravity data:           mass=',mass,' (kg) = '
     &     ,amlb,' (lb)'
      write(24,*)'              location of CG=',xcg,' (m)'
      write(24,*)'               static margin=',statmarg,' (if <0,'
     &     ,' not used to place xcg)'
      write(24,*)'***************'
      write(24,*)'main wing data:         span=',bm,' (m)'
      write(24,*)'  maximum wing or fuse chord=',cxm,' (m)'
      write(24,*)'  average wing or fuse chord=',cam,' (m)'
      write(24,*)'    wing+fuse projected area=',am,' (m**2)'
      write(24,*)'   location of w+f a.c. xacm=',xacm,' (m)'
      write(24,*)'                 wing AR arm=',arm
      write(24,*)'        wing relative camber=',dm
      write(24,*)'        wing+fuse efficiency=',em
      write(24,*)'downwash of tail on wing atw=',atw
      write(24,*)'          flap setting angle=',tfd,' (deg)'
      write(24,*)'flap setting influence on CL=',dClmdtf
      write(24,*)'flap setting influence on CM=',dCmmdtf
      write(24,*)'flap setting influence on CD=',dCdtf0
      write(24,*)'flap setting influence on CD=',dCdtf1
      write(24,*)'flap setting influence on CD=',dCdtf2
      write(24,*)'**************'
      write(24,*)'fuselage data:        length=',lf,' (m)'
      write(24,*)'                      radius=',rf,' (m)'
      write(24,*)'               fuselage area=',af,' (m**2)'
      write(24,*)'*********************'
      write(24,*)'tail/stabilizer data:   span=',bt,' (m)'
      write(24,*)'               maximum chord=',cxt,' (m)'
      write(24,*)'               average chord=',cat,' (m)'
      write(24,*)'                   tail area=',at,' (m**2)'
      write(24,*)'        tail relative camber=',dt
      write(24,*)'  location of tail a.c. xact=',xact,' (m)'
      write(24,*)'           tail aspect ratio=',art
      write(24,*)'             tail efficiency=',eteff
      write(24,*)'          tail setting angle=',tt,' (rd) ='
     &     ,ttd,' (deg)'
      write(24,*)'     flap percentage of tail=',flapcent,' (%)'
      write(24,*)'given taillift slope dCltda0=',dCltda0
      write(24,*)' effective tail aspect ratio=',arteff
      write(24,*)'             tail efficiency=',eteff
      write(24,*)'     flap percentage of tail=',flapcent,' (%)'
      write(24,*)'          downwash parameter=',awt
      write(24,*)'*****************************'
      write(24,*)'           airbrake: CDbrake=',Cdbrake
      write(24,*)'***********'
      write(24,*)'rudder data:'
      write(24,*)'    rudder average chord rav=',rav,' (m)'
      write(24,*)'           rudder height ruh=',ruh,' (m)'
      write(24,*)'              rudder area ar=',ar,' (m**2)'
      write(24,*)'***************'
      write(24,*)'reference data:  ref. length=',lref,' (m)'
      write(24,*)'                   ref. area=',aref,' (m**2)'
      write(24,*)'************'
      write(24,1000)prop
      write(24,*)' propulsion system reference=',prop
      write(24,*)'      z-coordinate of engine=',zeng,' (m)'
      write(24,*)'         diameter of nacelle=',dna,' (m)'
      write(24,*)'           length of nacelle=',lna,' (m)'
      write(24,*)'             nacelle area an=',an,' (m**2)'
      write(24,*)'             static traction=',tr0,' (N)'
      write(24,*)'              traction slope=',dtrdv,' (Ns/m)'
      write(24,*)'                engine power=',hpowatt,' (W) ='
     &     ,hpower,' (hp)'
      write(24,*)'percent of engine power avai=',Pcent
      write(24,*)'                        ndat=',ndat
      write(24,*)'    i=     vr(i)=      tr(i)='
      write(24,1004)(i,vr(i),tr(i),i=1,ndat)
      write(24,*)
      write(25,*)'******Pcent=',Pcent,'******'
c*****use polar data?
      write(6,*)
      write(6,*)'use polar data? Y/N=1/0'
      read(5,*)ipolar
c################### if no polar available
      if(ipolar.ne.1)goto 5
c*****polar data
      write(6,*)
      read(12,1000)title
      write(6,1000)title
      write(24,1000)title
      write(6,*)
      write(6,*)'*****extrema of the Cl(alpha) function:'
      prod=1.
      km=1
      do 2 k=1,lxx
         read(12,*,end=3)inc(k),cz(k),cx(k),dum,cq(k)
         kxtrm(k)=0
         dcz=cz(k)-cz(km)
         prod=prod*dcz
         if(prod.lt.0.)then
            write(6,*)'kxtrm=',km,' cz(kxtrm)=',cz(km)
            kxtrm(km)=km
         endif
         prod=dcz
         km=k
 2    continue
 3    continue
      kx=k-1
      write(6,*)
      write(6,*)'*************wing polar at zero setting angle 
     &     (profile data from Xfoil):'
      write(6,*)'  k=      inc=        cz=         cx=         cm=   '
      do 4 k=1,kx
         kp=k+1
         if(kp.gt.kx)kp=kx
         km=k-1
         if(km.lt.1)km=1
         prod=1.
         if(k.gt.1.and.k.lt.kx)then
         dcxm=((cx(k)-cx(km))*(cz(kp)-cz(k))*(cz(k)+cz(kp)-2.*cz(km))
     &        -(cx(kp)-cx(k))*(cz(k)-cz(km))**2)
     &        /((cz(kp)-cz(k))*(cz(kp)-cz(km))*(cz(k)-cz(km)))
         dcxp=((cx(k)-cx(kp))*(cz(km)-cz(k))*(cz(k)+cz(km)-2.*cz(kp))
     &        -(cx(km)-cx(k))*(cz(k)-cz(kp))**2)
     &        /((cz(km)-cz(k))*(cz(km)-cz(kp))*(cz(k)-cz(kp)))
         prod=dcxm*dcxp
         endif
         if(prod.lt.-eps.and.(kxtrm(km).ne.0.or.kxtrm(kp).ne.0))then
            write(6,*)'bad data distribution: interpolate a new data'
         endif
         incd=inc(k)
         write(13,*)cx(k),cz(k),cq(k),incd
         incd=inc(k)
         inc(k)=degrad*inc(k)
         inc(k)=inc(k)
         if(inc(km)*inc(k).le.0.and.k.ne.1)then
c*****zero incidence coefficients with assumption of small angles
c     lift coefficients of main wing
            dCldam0=(cz(k)-cz(km))/(inc(k)-inc(km))
            Clm00=cz(km)+dCldam0*(-inc(km))
c     moment coefficients of the main wing Cmacm and Cmmo
            dCmacdam0=(cq(k)-cq(km))/(inc(k)-inc(km))
            Cmacm00=cq(km)+dCmacdam0*(-inc(km))
            dCmdam0=dCmacdam0-xacm*dCldam0/cam
            Cmm00=Cmacm00-xacm*Clm00/cam
c     lift coefficients of tail from data file
            dCldat0=dCltda0
            Clt00=dCldat0*(dCltdtt*tt+2.*dt)
c     moment coefficients of tail Cmact and Cmto
            dCmacdat0=0.
            Cmact00=-pi*dt+dCmqdtt*tt
            dCmdat0=dCmacdat0-xact*dCldat0/cat
            Cmt00=Cmact00-xact*Clt00/cat
         endif
         if(k.eq.1.and.inc(1).ge.0.0)then
c     lift coefficients of main wing
            dCldam0=(cz(2)-cz(1))/(inc(2)-inc(1))
            Clm00=cz(1)+dCldam0*(-inc(1))
c     moment coefficients of the main wing Cmacm and Cmmo
            dCmacdam0=(cq(2)-cq(1))/(inc(2)-inc(1))
            Cmacm00=cq(1)+dCmacdam0*(-inc(1))
            dCmdam0=dCmacdam0-xacm*dCldam0/cam
            Cmm00=Cmacm00-xacm*Clm00/cam
c     lift coefficient of tail from data file
            dCldat0=dCltda0
            Clt00=dCldat0*(2.*dt+dCltdtt*tt)
c     moment coefficients of tail Cmact and Cmto
            dCmacdat0=0.
            Cmact00=-pi*dt+dCmqdtt*tt
            dCmdat0=dCmacdat0-xact*dCldat0/cat
            Cmt00=Cmact00-xact*Clt00/cat
         endif
c     aerodynamic center of airplane
         xac=(am*dCldam0*xacm+at*dCltda0*xact)
     &        /(am*dCldam0+at*dCltda0)
         write(14,*)cx(k),cz(k),cq(k),incd
         write(6,1001)k,incd,cz(k),cx(k),cq(k)
 4    continue
      write(6,*)
      write(6,*)'******extrema pointer:'
      write(6,*)'kxtrm(k)=',(kxtrm(k),k=1,kx)
 5    continue
c################### if no polar available
      if(ipolar.ne.1)then
         dCldam0=2.*pi/(1.+2./arm)
         Clm00=dCldam0*2.*dm
         dCmacdam0=0.
         dCmacdat0=0.
         Cmact0=-pi*dt+dCmqdtt*tt
         dCmdam0=-xacm*dCldam0/cam
         Cmm00=-pi*dm-xacm*Clm00/cam
         dCldat0=dCltda0
         Clt00=(2.*pi*2.*dt+dCltdtt*tt)/(1.+2./arteff)
         dCmdat0=-xact*dCldat0/cat
         Cmt00=Cmact0-xact*Clt00/cat
      endif
c     add downwash of wing on tail
      dCldat=dCltda0*(1.+awt*dCldam0/(pi*arm))
      Clt0=Clt00+dCltda0*awt*Clm00/(pi*arm)
c     add downwash of tail on wing
      dCldam=dCldam0*(1.+atw*dCltda0/(pi*arteff))
      Clm0=Clm00+dCldam*atw*Clt00/(pi*arteff)
c     global lift
      dClda=(am*dCldam+at*dCldat)/aref
      Cl0=(am*Clm0+at*Clt0)/aref
c     global moments
      dCmacda=(am*cam*dCmacdam0+at*cat*dCmacdat0)/(lref*aref)
      Cmac0=(am*cam*Cmacm00+at*cat*Cmact00)/(lref*aref)
      dCmda=(am*cam*dCmdam0+at*cat*dCmdat0)/(lref*aref)
      Cm0=(am*cam*Cmm00+at*cat*Cmt00)/(lref*aref)
      dCmdat=dCmdat0
      write(6,*)'****************************************'
      write(6,*)'zero incidence aerodynamic coefficients:'
      write(6,*)'**************lift:'
      write(6,*)' wing:     dCldam0=',dCldam0
     &     ,'   Clm00=',Clm00
      write(6,*)' tail:     dCldat0=',dCldat0
     &     ,'   Clt00=',Clt00
      write(6,*)'total:       dClda=',dClda
     &     ,'     Cl0=',Cl0
      write(6,*)'*****moment at xac:'
      write(6,*)' wing:   dCmacdam0=',dCmacdam0
     &     ,' Cmacm00=',Cmacm00
      write(6,*)' tail:   dCmacdat0=',dCmacdat0
     &     ,' Cmact00=',Cmact00
      write(6,*)'total:     dCmacda=',dCmacda
     &     ,'   Cmac0=',Cmac0
      write(6,*)'*****moment at x=0:'
      write(6,*)' wing:     dCmdam0=',dCmdam0
     &     ,'   Cmm00=',Cmm00
      write(6,*)' tail:     dCmdat0=',dCmdat0
     &     ,'   Cmt00=',Cmt00
      write(6,*)'total:       dCmda=',dCmda
     &     ,'     Cm0=',Cm0
      write(24,*)'****************************************'
      write(24,*)'zero incidence aerodynamic coefficients:'
      write(24,*)'**************lift:'
      write(24,*)' wing:     dCldam0=',dCldam0
     &     ,'   Clm00=',Clm00
      write(24,*)' tail:     dCldat0=',dCldat0
     &     ,'   Clt00=',Clt00
      write(24,*)'total:       dClda=',dClda
     &     ,'     Cl0=',Cl0
      write(24,*)'*****moment at xac:'
      write(24,*)' wing:   dCmacdam0=',dCmacdam0
     &     ,' Cmacm00=',Cmacm00
      write(24,*)' tail:   dCmacdat0=',dCmacdat0
     &     ,' Cmact00=',Cmact00
      write(24,*)'total:     dCmacda=',dCmacda
     &     ,'   Cmac0=',Cmac0
      write(24,*)'*****moment at x=0:'
      write(24,*)' wing:     dCmdam0=',dCmdam0
     &     ,'   Cmm00=',Cmm00
      write(24,*)' tail:     dCmdat0=',dCmdat0
     &     ,'   Cmt00=',Cmt00
      write(24,*)'total:       dCmda=',dCmda
     &     ,'     Cm0=',Cm0
c     global aerodynamic center
      xac=(am*dCldam0*xacm+at*dCltda0*xact)
     &     /(am*dCldam0+at*dCltda0)
      if(statmarg.gt.0.)then
         xcg=xac-statmarg*lref
      endif
      amarg=100.*(xac-xcg-eps)/lref
c*****linear model results
      aleq=-(Cm0+xcg*Cl0/lref)
     &     /(dCmda+xcg*dClda/lref)
      Cleq=Cl0+dClda*aleq
      Cweq=Cleq
      write(6,1002)
      write(24,1002)
      if(Cleq.le.0.)then
         write(6,*)'################# No linear solution with Cleq<0'
         Ueq=sqrt(2.*mg/(rho*aref*abs(Cleq)))
      else
         Ueq=sqrt(2.*mg/(rho*aref*Cleq))
      endif
      reyeq=Ref*Ueq/100.
      call thrust(ndatx,ndat,Ueq,vr,tr,T,dTdv)
      dynaref=0.5*rho*Ueq**2*aref
      Cteq=rcor*pcor**2*(T+dTdv*Ueq)/dynaref
c################### if no polar available
      if(ipolar.ne.1)then
         reym=reyeq*cam/lref
         Cdm0=1.328/sqrt(reym)
         if(reym.gt.1.0E5)then
            Cdm0=.072/reym**.2
         endif
         Cdeq=2.*Cdm0
         go to 8
      endif
c     search for point on wing polar and interpolate Cd
         m=1
         prod=aleq+0.5*pi
         do 6 k=2,kx-1
            prod=prod*(aleq-inc(k))
            if(prod.le.-eps)goto 7
            prod=1.
            m=k
 6          continue
 7       continue
         Cdeq=cx(m)+(aleq-inc(m))*(cx(m+1)-cx(m))
     &        /(inc(m+1)-inc(m))
 8       continue
c     wing induced drag estimate
         Cdmeq=Cdeq
         Clmeq=Clm0+dCldam0*aleq
         Cdim=Clmeq**2/(pi*em*arm)
         Cdm0=Cdmeq-Cdim
c     calculate fuselage friction drag
         reyf=reyeq*lf/cam
         Cdf0=1.328/sqrt(reyf)
         if(reyf.gt.1.0E5)then
            Cdf0=.072/reyf**.2
         endif
         Cdfeq=Cdf0
c     calculate tail induced drag and friction drag
         Clteq=Clt0+dCldat*aleq
         Cdit=2.0*Clteq**2/(pi*eteff*arteff)
         reyt=reyeq*cat/cam
         Cdt0=1.328/sqrt(reyt)
         if(reyt.gt.1.0E5)then
            Cdt0=.072/reyt**.2
         endif
         Cdteq=Cdit+2.0*Cdt0
c     calculate rudder drag for 2 sides
         reyr=reyeq*rav/cam
         Cdr0=1.328/sqrt(reyr)
         if(reyr.gt.1.0E5)then
            Cdr0=.072/reyr**.2
         endif
         Cdreq=2.0*Cdr0
c     if winglets are used as rudder, double contribution
c         Cdreq=4.0*Cdr0
c     calculate nacelles drag for 2 engines
         reyn=reyeq*lna/cam
         Cdn0=1.328/sqrt(reyn)
         if(reyn.gt.1.0E5)then
            Cdn0=.072/reyn**.2
         endif
         Cdneq=2.0*Cdn0
c     total drag
         Cdeq=(am*Cdmeq+af*Cdfeq+at*Cdteq+ar*Cdreq+an*Cdneq)/aref
         Cdeq=Cdeq+Cdbrake
         beteq=(Cteq-Cdeq)/Cweq
         Cmac=Cmac0+dCmacda*aleq
         Cmeq=Cm0+dCmda*aleq
         bb(1)=-(Cmeq+xcg*Cweq*cos(aleq+beteq)/lref)
         bb(2)=-(Cleq-Cweq*cos(beteq)+Cteq*sin(aleq))
         bb(3)=-(Cdeq+Cweq*sin(beteq)-Cteq*cos(aleq))
         write(6,*)'****linear model results: m=',m
         write(6,*)'aleq=',aleq,' Ueq=',Ueq,'beteq=',beteq
         write(6,*)'Cleq=',Cleq,'Cdeq=',Cdeq,' Cteq=',Cteq
     &        ,'CM,ac=',Cmac
         write(24,*)'aleq=',aleq,' Ueq=',Ueq,'beteq=',beteq
         write(24,*)'Cleq=',Cleq,'Cdeq=',Cdeq,' Cteq=',Cteq
     &        ,'CM,ac=',Cmac
c*****non-linear equations residuals with linear solution
         write(6,*)'*****************non-linear equations residuals'
     &        ,' from linear solution:'
         write(6,*)'equ1=',bb(1),'equ2=',bb(2),'equ3=',bb(3)
         write(24,*)'*****************non-linear equations residuals'
     &        ,' from linear solution:'
         write(24,*)'equ1=',bb(1),'equ2=',bb(2),'equ3=',bb(3)
c*****best design input data
      read(15,*)aleq
      read(15,*)Ueq
      read(15,*)beteq
      read(15,*)Cleq
      read(15,*)Cdeq
      read(15,*)Cmac
      aleq0=aleq
      Ueq0=Ueq
      beteq0=beteq
      Cleq0=Cleq
      Cdeq0=Cdeq
      Cmac0=Cmac
      write(6,*)'*****best design parameters:'
      write(6,*)'aleq=',aleq,' Ueq=',Ueq,'beteq=',beteq
      write(6,*)'Cleq=',Cleq,'Cdeq=',Cdeq,' CM,ac=',Cmac
      write(24,*)'*****best design parameters:'
      write(24,*)'aleq=',aleq,' Ueq=',Ueq,'beteq=',beteq
      write(24,*)'Cleq=',Cleq,'Cdeq=',Cdeq,' CM,ac=',Cmac
      Cteq=Cdeq
      Cweq=Cleq
      Cmeq=Cmac-xac*Cleq*cos(aleq)/lref
      Cmcg=Cmeq+xcg*Cweq*cos(aleq+beteq)/lref
      amarg=100.*(xac-xcg-eps)/lref
 100  continue
      write(6,*)'begin non-linear equilibrium algorithm'
      write(24,*)'begin non-linear equilibrium algorithm'
c*****read tail setting angle ttd (deg)
      write(6,*)'ttd=? (-100 exit)'
      write(24,*)'ttdd=? (-100 exit)'
      read(5,*)ttd
      write(6,*)' tail ttd=',ttd,'flap tfd=',tfd
      write(24,*)' tail ttd=',ttd,'flap tfd=',tfd
      if(ttd.lt.-99.)goto 300
      tt=degrad*ttd
c*****update with best design input data
      aleq=aleq0
      Ueq=Ueq0
      beteq=beteq0
      Cleq=Cleq0
      Cdeq=Cdeq0
      Cmac=Cmac0
      Cteq=Cdeq0
      Cweq=Cleq0
      Cmeq=Cmac-xac*Cleq*cos(aleq)/lref
      Cmcg=Cmeq+xcg*Cweq*cos(aleq+beteq)/lref
c     aerodynamic center moment coefficient from design data:
      write(6,*)'********************************************'
      write(6,*)'xcg, xac (xcg calculated or given):'
      write(6,*)'center of gravity:       xcg=',xcg,' (m)'
      write(6,*)'aerodynamic center:      xac=',xac,' (m)'
      write(6,*)'              moment at a.c.=',Cmac
      write(6,1003)amarg,lref
      write(24,*)'*******************************************'
      write(24,*)'xcg, xac (xcg calculate or given):'
      write(24,1002)
      write(24,*)'center of gravity:       xcg=',xcg,' (m)'
      write(24,*)'aerodynamic center:      xac=',xac,' (m)'
      write(24,*)'              moment at a.c.=',Cmac
      write(24,1003)amarg,lref
c*****iterations
      do 200 it=1,itx
c################### if no polar available
         if(ipolar.ne.1)goto 11
c     search for point on polar of main wing
         m=1
         prod=aleq+.5*pi
         do 9 k=2,kx-1
            prod=prod*(aleq-inc(k))
            if(prod.le.eps)goto 10
            prod=1.
            m=k
 9          continue
 10      continue
c     lift coefficients of the main wing
         dCldam=(cz(m+1)-cz(m))/(inc(m+1)-inc(m))
         Clm0=cz(m)+dCldam*(-inc(m))
c     moment coefficients of main wing and tail
         dCmacdam=(cq(m+1)-cq(m))/(inc(m+1)-inc(m))
         Cmacm0=cq(m)+dCmacdam*(-inc(m))
 11      continue
c     lift coefficient of tail with downwash of wing on tail
         dCldat=(2.*pi+awt*dCldam/(pi*arm))/(1.+2./arteff)
         Clt0=(2.*pi*(2.*dt+awt*Clm0/(pi*arm))+dCltdtt*tt)
     &        /(1.+2./arteff)
c     moment coefficients of tail
         dCmacdat0=0.
         Cmact0=-pi*dt+dCmqdtt*tt
         dCmdat=dCmacdat0-xact*dCldat/cat
         Cmt0=Cmact0-xact*Clt0/cat
c     add downwash of tail on wing
         dCldam=dCldam*(1.+atw*dCldat/(pi*arteff))
         Clm0=Clm0+dCldam*atw*Clt0/(pi*arteff)
         Clmeq=dCldam*aleq+Clm0
c     add flap influence on CLm
         Clm0=Clm0+dClmdtf*tf
         dCmdam=dCmacdam-xacm*dCldam/cam
         Cmm0=Cmacm0-xacm*Clm0/cam
c     add flap influence on CMm
         Cmm0=Cmm0+dCmmdtf*tf
         Cmmeq=dCmdam*aleq+Cmm0
c     global moment coefficient
         dCmda=(am*cam*dCmdam+at*cat*dCmdat)/(aref*lref)
         Cm0=(am*cam*Cmm0+at*cat*Cmt0)/(aref*lref)
         Cmeq=dCmda*aleq+Cm0
c*****relaxation of equ1
         daleq=-(Cmeq+xcg*Cweq*cos(aleq+beteq)/lref)
     &        /(dCmda-xcg*Cweq*sin(aleq+beteq)/lref)
c         daleq=0.
         aleq=aleq+omega*daleq
         aleqd=radeg*aleq
         if(ipolar.ne.1)then
c     -main wing
            Clmeq=dCldam0*aleq+Clm00
            Cmmeq=dCmdam0*aleq+Cmm00
            Cleq=dClda*aleq+Cl0
            Cmeq=dCmda*aleq+Cm0
         endif
c     global coefficients
         dClda=(am*dCldam+at*dCldat)/aref
         Cl0=(am*Clm0+at*Clt0)/aref
         Cleq=dClda*aleq+Cl0
c*****relaxation of equ2
         dUeq=-(Cleq-Cweq*cos(beteq)+Cteq*sin(aleq))
     &        /(2.*Cweq*cos(beteq)/Ueq)
c         dUeq=0.
         Ueq=Ueq+omega*dUeq
c     engine operating point
         call thrust(ndatx,ndat,Ueq,vr,tr,T,dTdv)
         dynaref=0.5*rho*aref*Ueq**2
         Cteq=rcor*pcor**2*(T+dTdv*Ueq)/dynaref
         Cweq=mg/dynaref
         reyeq=Ref*Ueq/100.
c     velocity and Reynolds number at equilibrium
         reym=reyeq
c     dynamic pressure
         dyn=.5*rho*Ueq*Ueq
c################### if no polar available
         if(ipolar.ne.1)then
         else
         if(kxtrm(m).ne.0)m=m+1
         if(m.lt.2)m=2
         if(m.gt.kx-1)m=kx-1
         dcxdcz=(cx(m)-cx(m-1))/(cz(m)-cz(m-1))
         d2cxdcz2=((cx(m+1)-cx(m))*(cz(m)-cz(m-1))
     &        -(cx(m)-cx(m-1))*(cz(m+1)-cz(m)))
     &     /((cz(m+1)-cz(m))*(cz(m+1)-cz(m-1))*(cz(m)-cz(m-1)))
c     main wing drag
         Cdmeq=cx(m-1)+dcxdcz*(Clmeq-cz(m-1))
     &        +d2cxdcz2*(Clmeq-cz(m-1))*(Clmeq-cz(m))
         endif
c     wing induced drag estimation
         Cdim=Clmeq**2/(pi*em*arm)
         Cdm0=Cdmeq-Cdim
c     add drag of flap
         Cdmeq=Cdmeq+(dCdtf0+dCdtf1*aleq+dCdtf2*aleq**2)*tf
c     fuselage friction drag
         reyf=reyeq*lf/cam
         Cdf0=1.328/sqrt(reyf)
         if(reyf.gt.1.0E5)then
            Cdf0=.072/reyf**.2
         endif
         Cdfeq=Cdf0
c     tail friction drag
         reyt=reyeq*cat/cam
         Cdt0=1.328/sqrt(reyt)
         if(reyt.gt.1.0E5)then
            Cdt0=.072/reyt**.2
         endif
c     tail induced drag estimation
         Clteq=dCldat*aleq+Clt0
         Cdit=Clteq**2/(pi*eteff*arteff)
         Cdt0=2.*Cdt0
c     tail total drag
         Cdteq=Cdt0+Cdit
c     calculate rudder drag for 2 sides
         reyr=reyeq*rav/cam
         Cdr0=1.328/sqrt(reyr)
         if(reyr.gt.1.0E5)then
            Cdr0=.072/reyr**.2
         endif
         Cdreq=2.0*Cdr0
c     if winglets are used as rudder, double contribution
c         Cdreq=4.0*Cdr0
c     calculate nacelles drag for 2 engines
         reyn=reyeq*lna/cam
         Cdn0=1.328/sqrt(reyn)
         if(reyn.gt.1.0E5)then
            Cdn0=.072/reyn**.2
         endif
         Cdneq=2.0*Cdn0
c################### if no polar available
         if(ipolar.ne.1)then
c     viscous and induced drag of main wing
            reym=reyeq
            Cdim=Clmeq**2/(pi*em*arm)
            Cdm0=1.328/sqrt(reym)
            if(reym.gt.1.0E5)then
               Cdm0=.072/reym**.2
            endif
            Cdm0=2.*Cdm0
            Cdmeq=Cdm0+Cdim
         endif
c     total drag
         Cdeq=(am*Cdmeq+af*Cdfeq+at*Cdteq+ar*Cdreq+an*Cdneq)/aref
         Cdeq=Cdeq+Cdbrake
c*****relaxation of equ3
         dbeteq=-(Cdeq+Cweq*sin(beteq)-Cteq*cos(aleq))
     &        /(Cweq*cos(beteq))
c         dbeteq=0.
         beteq=beteq+omega*dbeteq
c     angle of descent
         if(abs(beteq).lt.1.)then
            beteq=asin(beteq)
         else
            beteq=0.5*pi
         endif
         beteqd=radeg*beteq
         theq=aleq+beteq
         theqd=aleqd+beteqd
c     finess coefficient
         ClCdeq=Cleq/Cdeq
c     residuals
         bb(1)=-(Cmeq+xcg*Cweq*cos(aleq+beteq)/lref)
         bb(2)=-(Cleq-Cweq*cos(beteq)+Cteq*sin(aleq))
         bb(3)=-(Cdeq+Cweq*sin(beteq)-Cteq*cos(aleq))
         reseq=bb(1)**2+bb(2)**2+bb(3)**2
         reseq=sqrt(reseq)
         reseq0=reseq
         inewton=0
         det=1.0
         usdet=1./det
         if(reseq0.lt.10.*eps)then
c     Newton's Method
            inewton=1
            b1=bb(1)
            b2=bb(2)
            b3=bb(3)
            aa(1,1)=dCmda-xcg*Cweq*sin(aleq+beteq)/lref
            aa(1,2)=2.*bb(1)/Ueq
            aa(1,3)=-xcg*Cweq*sin(aleq+beteq)/lref
            aa(2,1)=dClda+Cteq*cos(aleq)
            aa(2,2)=2.*Cweq*cos(beteq)/Ueq
            aa(2,3)=Cweq*sin(beteq)
            aa(3,1)=0.*Cleq/(pi*em*arm)*dClda+Cteq*sin(aleq)
            aa(3,2)=2.*bb(3)/Ueq-dTdv*cos(aleq)/dynaref
            aa(3,3)=Cweq*cos(beteq)
c
            call mat3(usdet,aa,bb)
c
            aleq=aleq+bb(1)
            Ueq=Ueq+bb(2)
            beteq=beteq+bb(3)
            reseq=b1**2+b2**2+b3**2
            reseq=sqrt(reseq)
            det=1./usdet
         endif
c     engine operating point
         call thrust(ndatx,ndat,Ueq,vr,tr,T,dTdv)
         dynaref=0.5*rho*Ueq**2*aref
         Cteq=rcor*pcor**2*(T+dTdv*Ueq)/dynaref
         write(26,*)it,inewton,reseq
         if(reseq.lt.epser)goto 12
 200  continue
 12   continue
c*****non-linear equations residuals
      if(inewton.eq.0)then
         write(6,*)'relaxation:'
         write(6,*)'daleq=',daleq
     &        ,' dUeq=',dUeq,'dbeteq=',dbeteq
         write(6,*)' equ1=',bb(1)
     &        ,' equ2=',bb(2),'  equ3=',bb(3)
      elseif(inewton.eq.1)then
            write(6,*)'Newton:'
            write(6,*)'  det=',det
            write(6,*)'daleq=',bb(1)
     &           ,' dUeq=',bb(2),'dbeteq=',bb(3)
            write(6,*)' equ1=',b1
     &           ,' equ2=',b2,'  equ3=',b3
      endif
c     reults
      weq=Ueq*sin(beteq)
      winglift=dyn*am*clmeq
      taillift=dyn*at*clteq
      tr0=rcor*pcor**2*(T+dTdv*Ueq)
      xcpm=xacm-cam*Cmacm0/Clmeq
      xcpt=xact-cat*Cmact0/Clteq
      write(6,*)'*****************************'
      write(6,*)'          tail setting angle=',tt,' (rd) ='
     &     ,ttd,' (deg)'
      if(it.ge.itx)then
         write(6,*)'equilibrium solution:   it=',it
     &     ,' error=',reseq,' ',nc
      else
         if(aleq.gt.inc(1))then
            write(6,*)'equilibrium solution:   it=',it
     &     ,' error=',reseq
         else
            write(6,*)'equilibrium solution:   it=',it
     &     ,' error=',reseq,' ',yw
         endif
      endif
      write(6,*)
      write(6,*)'global results:              '
     &     ,'        (reference area=',aref,' (m**2))'
      write(6,*)'                   incidence=',aleq,' (rd) ='
     &     ,aleqd,' (deg)'
      write(6,*)'                 climb slope=',beteq,' (rd) ='
     &     ,beteqd,' (deg)'
      write(6,*)'            airplane setting=',theq,' (rd) ='
     &     ,theqd,' (deg)'
      write(6,*)'                    velocity=',Ueq,' (m/s)'
      write(6,*)'                  climb rate=',weq,' (m/s)'
      write(6,*)'            dynamic pressure=',dyn,' (Pa)'
      write(6,*)'       dynamic pressure*aref=',dynaref,' (n)'
      write(6,*)'             reynolds number=',reyeq
      write(6,*)'                     lift CL=',Cleq
      write(6,*)'       weight coefficient CW=',Cweq
      write(6,*)'                moment CM,ac=',Cmac
      write(6,*)'                     drag CD=',Cdeq
      write(6,*)'       thrust coefficient CT=',Cteq
      write(6,*)'                     Cdbrake=',Cdbrake
      write(6,*)'                       CL/CD=',ClCdeq
      write(6,*)'                 moment CM,O=',Cmeq
      write(6,*)
      write(6,*)'main wing:                   '
     &     ,'          (ref/wing area=',am,' (m**2))'
      write(6,*)'                    lift CLm=',Clmeq
     &     ,'    Lm=',winglift,' (N)'
      write(6,*)'  wing pitching moment CMacm=',Cmacm0
      write(6,*)'               wing Reynolds=',reym
      write(6,*)' estimated induced drag CDim=',Cdim
      write(6,*)'                    drag CDm=',Cdmeq
      write(6,*)'                        xcpm=',xcpm,' (m)'
      write(6,*)'tail:                        '
     &     ,'          (ref/tail area=',at,' (m**2))'
      write(6,*)'                    lift CLt=',Clteq
     &     ,'    Lt=',taillift,' (N)'
      write(6,*)' tail pitching moment CM,act=',Cmact0
      write(6,*)'               tail Reynolds=',reyt
      write(6,*)' estimated induced drag CDit=',Cdit
      write(6,*)'                    drag CDt=',Cdteq
      write(6,*)'                        xcpt=',xcpt,' (m)'
      write(6,*)'fuselage:                    '
     &     ,'          (ref/fus. area=',af,' (m**2))'
      write(6,*)'           fuselage Reynolds=',reyf
      write(6,*)'                    drag CDf=',Cdfeq
      write(6,*)'rudder:                      '
     &     ,'          (ref/rud. area=',ar,' (m**2)'
      write(6,*)'             rudder Reynolds=',reyr
      write(6,*)'                    drag Cdr=',Cdreq
      write(6,*)'nacelle:                     '
     &     ,'          (ref/wet. area=',an,' (m**2)'
      write(6,*)'            nacelle Reynolds=',reyn
      write(6,*)'                    drag Cdn=',Cdneq
      write(6,1005)Pcent,tr0
      if(it.lt.itx)then
         imarg=amarg+0.1
         ittd=ttd+0.1*sign(1.0,ttd)
         write(25,*)imarg,ittd,aleq,beteq,Ueq,tr0
      endif
      write(24,*)'*****************************'
      write(24,*)'          tail setting angle=',tt,' (rd) ='
     &     ,ttd,' (deg)'
      if(it.ge.itx)then
         write(24,*)'equilibrium solution:   it=',it
     &     ,' error=',daleq,' ',nc
      else
         write(24,*)'equilibrium solution:   it=',it
     &     ,' error=',daleq
      endif
      write(24,*)
      write(24,*)'global results:              '
     &     ,'        (reference area=',aref,' (m**2))'
      write(24,*)'                   incidence=',aleq,' (rd) ='
     &     ,aleqd,' (deg)'
      write(24,*)'                 climb slope=',beteq,' (rd) ='
     &     ,beteqd,' (deg)'
      write(24,*)'            airplane setting=',theq,' (rd) ='
     &     ,theqd,' (deg)'
      write(24,*)'                    velocity=',Ueq,' (m/s)'
      write(24,*)'                  climb rate=',weq,' (m/s)'
      write(24,*)'            dynamic pressure=',dyn,' (Pa)'
      write(24,*)'       dynamic pressure*aref=',dynaref,' (n)'
      write(24,*)'             reynolds number=',reyeq
      write(24,*)'                     lift CL=',Cleq
      write(24,*)'       weight coefficient CW=',Cweq
      write(24,*)'                moment CM,ac=',Cmac
      write(24,*)'                     drag CD=',Cdeq
      write(24,*)'       thrust coefficient CT=',Cteq
      write(24,*)'                     Cdbrake=',Cdbrake
      write(24,*)'                       CL/CD=',ClCdeq
      write(24,*)'                 moment CM,O=',Cmeq
      write(24,*)
      write(24,*)'main wing:                   '
     &     ,'          (ref/wing area=',am,' (m**2))'
      write(24,*)'                    lift CLm=',Clmeq
      write(24,*)'  wing pitching moment CMacm=',Cmacm0
      write(24,*)'               wing Reynolds=',reym
      write(24,*)' estimated induced drag CDim=',Cdim
      write(24,*)'                    drag CDm=',Cdmeq
      write(24,*)'                        xcpm=',xcpm,' (m)'
      write(24,*)'tail:                        '
     &     ,'          (ref/tail area=',at,' (m**2))'
      write(24,*)'                    lift CLt=',Clteq
     &     ,'    Lt=',taillift,' (N)'
      write(24,*)' tail pitching moment CM,act=',Cmact0
      write(24,*)'               tail Reynolds=',reyt
      write(24,*)' estimated induced drag CDit=',Cdit
      write(24,*)'                    drag CDt=',Cdteq
      write(24,*)'                        xcpt=',xcpt,' (m)'
      write(24,*)'fuselage:                    '
     &     ,'          (ref/fus. area=',af,' (m**2))'
      write(24,*)'           fuselage Reynolds=',reyf
      write(24,*)'                    drag CDf=',Cdfeq
      write(24,*)'rudder:                      '
     &     ,'          (ref/rud. area=',ar,' (m**2)'
      write(24,*)'             rudder Reynolds=',reyr
      write(24,*)'                    drag Cdr=',Cdreq
      write(24,*)'nacelle:                     '
     &     ,'          (ref/wet. area=',an,' (m**2)'
      write(24,*)'            nacelle Reynolds=',reyn
      write(24,*)'                    drag Cdn=',Cdneq
      write(24,1005)Pcent,tr0
      goto 100
 300  continue
      ttd=tt*radeg
      aleqd=radeg*aleq
      write(21,*)xcg,ttd,aleqd,beteqd,Ueq,Clmeq
c*****loop on incidence
      nalmin=radeg*inc(1)
      nalmax=nalmin+nal
      do 16 n=nalmin,nalmax
         alphd=n
         alph=degrad*alphd
c*****lift drag and moment of main wing
c################### if no polar available
         if(ipolar.ne.1)goto 15
c     search for point on polar of main wing
         m=1
         prod=alph+.5*pi
         do 13 k=2,kx-1
            prod=prod*(alph-inc(k))
            if(prod.le.eps)goto 14
            prod=1.
            m=k
 13      continue
 14      continue
         dCldam=(cz(m+1)-cz(m))/(inc(m+1)-inc(m))
         Clm0=cz(m)+dCldam*(-inc(m))
         Clmnal=dCldam*alph+Clm0
         dCmdam=(cq(m+1)-cq(m))/(inc(m+1)-inc(m))
         Cmm0=cq(m)+dCmdam*(-inc(m))
         dCmdam=dCmdam-xacm*dCldam/cam
         Cmm0=Cmm0-xacm*Clm0/cam
         Cmmnal=dCmdam*alph+Cmm0
         if(kxtrm(m).ne.0)m=m+1
         if(m.lt.2)m=2
         if(m.gt.kx-1)m=kx-1
         dcxdcz=(cx(m)-cx(m-1))/(cz(m)-cz(m-1))
         d2cxdcz2=((cx(m+1)-cx(m))*(cz(m)-cz(m-1))
     &        -(cx(m)-cx(m-1))*(cz(m+1)-cz(m)))
     &     /((cz(m+1)-cz(m))*(cz(m+1)-cz(m-1))*(cz(m)-cz(m-1)))
         Cdmnal=cx(m-1)+dcxdcz*(Clmnal-cz(m-1))
     &        +d2cxdcz2*(Clmnal-cz(m-1))*(Clmnal-cz(m))
         dClda=(am*dCldam+at*dCldat)/aref
         Cl0=(am*Clm0+at*Clt0)/aref
         dCmda=(am*cam*dCmdam+at*cat*dCmdat)/(aref*lref)
         Cm0=(am*cam*Cmm0+at*cat*Cmt0)/(aref*lref)
 15      continue
c################### if no polar available
         if(ipolar.ne.1)then
            Clm=dCldam0*alph+Clm00
            Cdm=Cdm0+Clm**2/(pi*em*arm)
            Cmm=dCmdam0*alph+Cmm00
         else
            Clm=Clmnal
            Cdm=Cdmnal
            Cmm=Cmmnal
         endif
c*****lift drag and moment of tail
         Clt=dCldat*alph+Clt0
         Cdt=Cdt0+Clt**2/(pi*eteff*arteff)
         Cmt=dCmdat*alph+Cmt0
c*****global lift drag and moment
         Cl=dClda*alph+Cl0
         Cm=dCmda*alph+Cm0
         Cd=(am*Cdm+at*Cdt+af*Cdf0)/aref
         Cl3=Cl**3
         Cd2=Cd**2
         ratio1=Cl/Cd
         ratio2=Cl3/Cd2
         incidence=n
         write(16,*)Cdm,Clm,Cmm,incidence
         write(17,*)Cdt,Clt,Cmt,incidence
         write(18,*)Cd,Cl,Cm,incidence
         write(19,*)Cl,Cd,ratio1,Cd2,Cl3,ratio2,incidence
         Cmcg=Cm+xcg/lref*Cl
         write(20,*)incidence,Cmcg
 16   continue
      Cmmeq=Cmmeq+xacm*Clmeq/cam
      Cmcgeq=Cmeq+xcg*Cleq/lref
      write(22,*)Cdmeq,Clmeq,Cmmeq,aleqd
      write(23,*)Cdeq,Cleq,Cmeq,aleqd,Cmcgeq
c*****files
      write(6,*)
      write(6,*)'******input file'
      write(6,*)'polaruav.dat           :polar from prandtline & XFOIL'
     &     ,' (optional)'
      write(6,*)'classiceq.data         :parameters file'
      write(6,*)'******output files'
      write(6,*)'polaruav.prandtline    :CDm,CLm,CMm,alpha prandtline'
      write(6,*)'polarclassic.cdclcq    :CDm,CLm,CMm,alpha'
      write(6,*)'classiceq.cdclcqm      :CDm,CLm,CMm,alpha'
      write(6,*)'classiceq.cdclcqt      :CDt,CLt,CMt,alpha'
      write(6,*)'classiceq.cdclcq       :CD,CL,CM,alpha'
      write(6,*)'classieq.cdmcln        :CD,CL,CL/CD,CD**2,CL**3'
     &     ,'CL3/CD2',',alpha'
      write(6,*)'classiceq.cmcg         :alpha,CMcg'
      write(6,*)'classiceq.xtabvcl      :xcg,tt,aleq,beteq,Ueq,Clmeq'
      write(6,*)'classiceq.cdclcqmeq    :CDmeq,CLmeq,CMmeq,aleq'
      write(6,*)'classiceq.cdclcqcmcgeq :CDeq,CLeq,CMeq,aleq,CMcgeq'
      write(6,*)'classiceq.list         :listing of results'
      write(6,*)'classiceq.stab         :imarg,ittd,aleq,beteq,Ueq,tr0'
      write(6,*)'classiceq.itres        :it,residual'
 1000 format(18a4)
 1001 format(1x,i4,4x,f8.4,4x,f8.4,4x,f8.4,4x,f8.4)
 1002 format(' #####################################################',
     &     '###############')
 1003 format('static margin=',f5.1,' %'
     &     ,'                  (reference lref=',f12.7,'       (m)')
 1004 format(3x,i3,4x,f8.4,4x,e12.4)
 1005 format(' Pcent=',f8.0,'        thrust=',f12.4,'       (N)')
      end

      subroutine thrust(ndatx,ndat,v,vr,tr,T,dTdv)
      implicit none
      integer ndatx,ndat,i
      real v,T,dTdv,prod
      real vr(ndatx),tr(ndatx)
      prod=v-vr(1)
      do 1 i=2,ndat
         prod=prod*(v-vr(i))
         if(prod.le.0.)then
            dTdv=(tr(i)-tr(i-1))/(vr(i)-vr(i-1))
            T=tr(i-1)+dTdv*(-vr(i-1))
            goto 2
         endif
         prod=1.
 1    continue
      dTdv=(tr(ndat)-tr(ndat-1))/(vr(ndat)-vr(ndat-1))
      T=tr(ndat-1)+dTdv*(-vr(ndat-1))
 2    continue
      return
      end

      subroutine mat3(usdet,a,b)
      implicit none
      doubleprecision usdet,b1,b2,b3
      doubleprecision a(3,3),b(3)
      usdet=1./(a(1,1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2))
     &	       -a(2,1)*(a(1,2)*a(3,3)-a(1,3)*a(3,2))
     &	       +a(3,1)*(a(1,2)*a(2,3)-a(1,3)*a(2,2)))
      b1=(b(1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2))
     &   -b(2)*(a(1,2)*a(3,3)-a(1,3)*a(3,2))
     &   +b(3)*(a(1,2)*a(2,3)-a(1,3)*a(2,2)))*usdet
      b2=(-b(1)*(a(2,1)*a(3,3)-a(2,3)*a(3,1))
     &    +b(2)*(a(1,1)*a(3,3)-a(1,3)*a(3,1))
     &    -b(3)*(a(1,1)*a(2,3)-a(1,3)*a(2,1)))*usdet
      b3=(b(1)*(a(2,1)*a(3,2)-a(2,2)*a(3,1))
     &   -b(2)*(a(1,1)*a(3,2)-a(1,2)*a(3,1))
     &   +b(3)*(a(1,1)*a(2,2)-a(1,2)*a(2,1)))*usdet
      b(1)=b1
      b(2)=b2
      b(3)=b3
      return
      end


      program canareq
      implicit none
      integer lxx,ndatx,ndat,nal,itx,i,km,k,kx,kp,k0,m,iter,it
      integer inewton,imarg,itfd,nalmin,nalmax,n,incidence
      parameter(lxx=101,ndatx=21)
      integer kxtrm(lxx)
      real pi,eps,degrad,radeg,us3,omega,epser,rho,amu,mass
      real xcg,statmarg,amlb,mg,bm,cxm,cam,am,xac,tfd,dClmdtf
      real dCmmdtf,dCdtf0,dCdtf1,dCdtf2,dm,em,tf,rf,hf,af
      real Cdbrake,aref,lref,hpower,Pcent,Ueq,Ref
      real T,dTdv,tr0,dtrdv,arm,prod,dum,dcz,dcxm,dcxp,incd
      real dCldam0,Clm00,dCmacdam0,dCmdam0,Cmacm00,Cmm00,Cdm,dClda
      real dCmda,aleq,aleqd,Cleq,dynaref,Cdeq,Cteq,Cweq,Clm0
      real beteq,Clmeq,Cmac0,Cmcg,dCldam,dCmacda,dCmdam
      real dCmeq,Cmac,amarg,dcxdcz,d2cxdcz2,Cdmeq,Cl0,Cm0,Cdm0,Cdim
      real daleq,dUeq,reyeq,reym,Cdf0,Cdfeq,dbeteq,beteqd,theq,theqd
      real ClCdeq,reseq,winglift,alphd,alph,Cl,Cm,Cd,Cl3,Cd2,ratio1
      real ratio2,Cmmeq,Cmeq,Cmm0,dyn,hpokwatt,lf,reyf,weq
      real reseq0,bc,cxc,cac,ac,xacc,dc,ec,tcd,tc,awc,acw,xacm
      real arc,dCmacdac0,Cmacc00,Clc00,dCldac0,Cmac00,Cmc00,dCmdac0
      real Clceq,canarlift,reyc,Cdic,Cdceq,Clc0,Cmacc0,Cmacm0,Cmacm
      real Cmc0,dCldac,dCmacdac,dCmacdam,dCmdac,Cdc0,Cmacc
      real aleq0,Ueq0,beteq0,Cleq0,Cdeq0,rcor,xcpm,xcpc,dCldalind
      real zeng,rav,ruh,ar,reyr,Cdr0,Cdreq,dna,lna,reyn,an
      real Cdn0,Cdneq,dClmda0,dClcda0,armeff,arceff
      real vr(ndatx),tr(ndatx)
      real cx(lxx),cz(lxx),cq(lxx),inc(lxx)
      doubleprecision det,usdet,b1,b2,b3
      doubleprecision aa(3,3),bb(3)
      character*4 title(18),prop(10),nc(5),yw(5)
      open(unit=12,file='canarpolar.dat',form='formatted')
      open(unit=13,file='canarpolar.prandtline',form='formatted')
      open(unit=14,file='canarpolar.cdclcq',form='formatted')
      open(unit=15,file='canareq.data',form='formatted')
      open(unit=16,file='canareq.cdclcqm',form='formatted')
      open(unit=17,file='canareq.itres',form='formatted')
      open(unit=18,file='canareq.cdclcq',form='formatted')
      open(unit=19,file='canareq.cdmcln',form='formatted')
      open(unit=20,file='canareq.cmcg',form='formatted')
      open(unit=21,file='canareq.xtabvcl',form='formatted')
      open(unit=22,file='canareq.cdclcqmeq',form='formatted')
      open(unit=23,file='canareq.cdclcqcmcgeq',form='formatted')
      open(unit=24,file='canareq.list',form='formatted')
      open(unit=25,file='canareq.stab',form='formatted')
      data nc/'!not',' con','verg','ed!!','!!!!'/
      data yw/'!ale','q ou','t of',' ran','ge!!'/
c*****constants
      pi=2.*asin(1.)
      eps=1.e-6
      us3=1.0/3.0
      degrad=pi/180.
      radeg=1./degrad
      nal=30
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
c*****data for wing
      read(15,*)bm
      read(15,*)cxm
      read(15,*)cam
      read(15,*)am
      read(15,*)xacm
      read(15,*)dClmda0
      read(15,*)armeff
      read(15,*)dm
      read(15,*)em
      read(15,*)acw
      read(15,*)dCldalind
      read(15,*)tfd
      read(15,*)dClmdtf
      read(15,*)dCmmdtf
      read(15,*)dCdtf0
      read(15,*)dCdtf1
      read(15,*)dCdtf2
      arm=armeff
      tf=degrad*tfd
c*****data for fuselage
      read(15,*)lf
      read(15,*)rf
      hf=2.*pi*rf
      af=lf*hf
c*****data for canard
      read(15,*)bc
      read(15,*)cxc
      read(15,*)cac
      read(15,*)ac
      read(15,*)xacc
      read(15,*)dClcda0
      read(15,*)arceff
      read(15,*)dc
      read(15,*)ec
      read(15,*)awc
      arc=arceff
c*****airbrake data
      read(15,*)Cdbrake
c*****data for rudder
      read(15,*)rav
      read(15,*)ruh
      ar=rav*ruh
c*****reference parameters
      aref=am+ac
      lref=lf
      Ueq=100.
      Ref=rho*Ueq*cam/amu
c*****engine data
      read(15,*)zeng
      read(15,*)dna
      read(15,*)lna
      read(15,*)hpower
      read(15,*)Pcent
      read(15,*)ndat
      read(15,1000)prop
      read(15,1000)title
      do 1 i=1,ndat
         read(15,*)vr(i),tr(i)
 1    continue
      an=pi*dna*lna
      hpokwatt=745.7*hpower/1000.
      Ueq=0.
      call thrust(ndatx,ndat,Ueq,vr,tr,T,dTdv)
c*****thrust and thrust slope correction with density
      rcor=(rho/1.225)**us3
      rcor=rcor*(Pcent*Pcent/10000.0)**us3
      tr0=rcor*(T+dTdv*Ueq)
      dtrdv=rcor*dTdv
c*****write data
      write(6,1002)
      write(6,1002)
      write(6,*)'***********************'
      write(6,*)'convergence parameters:'
      write(6,*)'                         itx=',itx
      write(6,*)'                       omega=',omega
      write(6,*)' convergence tolerance epser=',epser
      write(6,*)'***************'
      write(6,*)'air parameters:      density=',rho,' (kg/m**3)'
      write(6,*)'               dynamic visc.=',amu,' (m**2/s)'
      write(6,*)'*************'
      write(6,*)'gravity data:           mass=',mass,' (kg) = '
     &     ,amlb,' (lb)'
      write(6,*)'              location of CG=',xcg,' (m)'
      write(6,*)'               static margin=',statmarg,' (if <0,'
     &     ,' not used to place xcg)'
      write(6,*)'**********'
      write(6,*)'wing data:              span=',bm,' (m)'
      write(6,*)'          root chord of wing=',cxm,' (m)'
      write(6,*)'               average chord=',cam,' (m)'
      write(6,*)'                   wing area=',am,' (m**2)'
      write(6,*)'  aerodynamic center of wing=',xacm,' (m)'
      write(6,*)'     wing lift slope dClmda0=',dClmda0
      write(6,*)'    wing effective AR armeff=',armeff
      write(6,*)'     relative camber of wing=',dm
      write(6,*)'   wing effective efficiency=',em
      write(6,*)' influence fo canard on wing=',acw
      write(6,*)'influence of canar dCldalind=',dCldalind
      write(6,*)'          flap setting angle=',tfd,' (deg)'
      write(6,*)'flap setting influence on Cl=',dClmdtf
      write(6,*)'flap setting influence on Cm=',dCmmdtf
      write(6,*)'flap setting influence on Cd=',dCdtf0
      write(6,*)'flap setting influence on Cd=',dCdtf1
      write(6,*)'flap setting influence on Cd=',dCdtf2
      write(6,*)'**************'
      write(6,*)'fuselage data:        length=',lf,' (m)'
      write(6,*)'                      radius=',rf,' (m)'
      write(6,*)'   fuselage circumference hf=',hf,' (m)'
      write(6,*)'               fuselage area=',af,' (m**2)'
      write(6,*)'************'
      write(6,*)'canard data:      canard span=',bc,' (m)'
      write(6,*)'           canard root chord=',cxc,' (m)'
      write(6,*)'canar mean aerodynamic chord=',cac,' (m)'
      write(6,*)'        canard planform area=',ac,' (m**2)'
      write(6,*)'aerodynamic center of canard=',xacc,' (m)'
      write(6,*)'   canard lift slope dClcda0=',dClcda0
      write(6,*)'  canard effective AR arceff=',arceff
      write(6,*)'   relative camber of canard=',dc
      write(6,*)'           canard efficiency=',ec
      write(6,*)' influence of wing on canard=',awc
      write(6,*)'********************'
      write(6,*)'           airbrake: CDbrake=',Cdbrake
      write(6,*)'************'
      write(6,*)'rudder data:'
      write(6,*)'        rudder average chord=',rav,' (m)'
      write(6,*)'               rudder height=',ruh,' (m)'
      write(6,*)'                 rudder area=',ar,' (m**2)'
      write(6,*)'***************'
      write(6,*)'reference data:  ref. length=',lref,' (m)'
      write(6,*)'                   ref. area=',aref,' (m**2)'
      write(6,*)'   reference Reynolds number=',Ref
      write(6,*)'*****************************'
      write(6,*)' propulsion system reference=   ',prop
      write(6,*)'      z-coordinate of engine=',zeng,' (m)'
      write(6,*)'         diameter of nacelle=',dna,' (m)'
      write(6,*)'           length of nacelle=',lna,' (m)'
      write(6,*)'             static traction=',tr0,' (N)'
      write(6,*)'              traction slope=',dtrdv,' (Ns/m)'
      write(6,*)'                engine power=',hpokwatt,' (kW) ='
     &     ,hpower,' (hp)'
      write(6,*)'percent of engine power avai=',Pcent
      write(6,*)'                        ndat=',ndat
      write(6,*)'    i=     vr(i)=      tr(i)='
      write(6,1004)(i,vr(i),tr(i),i=1,ndat)
      write(6,*)
      write(24,1002)
      write(24,1002)
      write(24,*)'***********************'
      write(24,*)'convergence parameters:'
      write(24,*)'                         itx=',itx
      write(24,*)'                       omega=',omega
      write(24,*)' convergence tolerance epser=',epser
      write(24,*)'***************'
      write(24,*)'air parameters:      density=',rho,' (kg/m**3)'
      write(24,*)'               dynamic visc.=',amu,' (m**2/s)'
      write(24,*)'*************'
      write(24,*)'gravity data:           mass=',mass,' (kg) = '
     &     ,amlb,' (lb)'
      write(24,*)'              location of CG=',xcg,' (m)'
      write(24,*)'               static margin=',statmarg,' (if <0,'
     &     ,' not used to place xcg)'
      write(24,*)'***************'
      write(24,*)'wing data:              span=',bm,' (m)'
      write(24,*)'          root chord of wing=',cxm,' (m)'
      write(24,*)'               average chord=',cam,' (m)'
      write(24,*)'                   wing area=',am,' (m**2)'
      write(24,*)'  aerodynamic center of wing=',xacm,' (m)'
      write(24,*)'     wing lift slope dClmda0=',dClmda0
      write(24,*)'    wing effective AR armeff=',armeff
      write(24,*)'     relative camber of wing=',dm
      write(24,*)'   wing effective efficiency=',em
      write(24,*)' influence fo canard on wing=',acw
      write(24,*)'influence of canar dCldalind=',dCldalind
      write(24,*)'          flap setting angle=',tfd,' (deg)'
      write(24,*)'flap setting influence on Cl=',dClmdtf
      write(24,*)'flap setting influence on Cm=',dCmmdtf
      write(24,*)'flap setting influence on Cd=',dCdtf0
      write(24,*)'flap setting influence on Cd=',dCdtf1
      write(24,*)'flap setting influence on Cd=',dCdtf2
      write(24,*)'**************'
      write(24,*)'fuselage data:        length=',lf,' (m)'
      write(24,*)'                      radius=',rf,' (m)'
      write(24,*)'   fuselage circumference hf=',hf,' (m)'
      write(24,*)'               fuselage area=',af,' (m**2)'
      write(24,*)'***********'
      write(24,*)'canard data:      canard span=',bc,' (m)'
      write(24,*)'           canard root chord=',cxc,' (m)'
      write(24,*)'canar mean aerodynamic chord=',cac,' (m)'
      write(24,*)'        canard planform area=',ac,' (m**2)'
      write(24,*)'aerodynamic center of canard=',xacc,' (m)'
      write(24,*)'   canard lift slope dClcda0=',dClcda0
      write(24,*)'  canard effective AR arceff=',arceff
      write(24,*)'   relative camber of canard=',dc
      write(24,*)'           canard efficiency=',ec
      write(24,*)' influence of wing on canard=',awc
      write(24,*)'************'
      write(24,*)'rudder data:'
      write(24,*)'        rudder average chord=',rav,' (m)'
      write(24,*)'               rudder height=',ruh,' (m)'
      write(24,*)'                 rudder area=',ar,' (m**2)'
      write(24,*)'*****************************'
      write(24,*)' propulsion system reference=   ',prop
      write(24,*)'      z-coordinate of engine=',zeng,' (m)'
      write(24,*)'         diameter of nacelle=',dna,' (m)'
      write(24,*)'           length of nacelle=',lna,' (m)'
      write(24,*)'             static traction=',tr0,' (N)'
      write(24,*)'              traction slope=',dtrdv,' (Ns/m)'
      write(24,*)'                engine power=',hpokwatt,' (kW) ='
     &     ,hpower,' (hp)'
      write(24,*)'percent of engine power avai=',Pcent
      write(24,*)'                        ndat=',ndat
      write(24,*)'    i=     vr(i)=      tr(i)='
      write(24,1004)(i,vr(i),tr(i),i=1,ndat)
      write(24,*)
      write(25,*)'******Pcent=',Pcent,'******'
c*****use polar data
      write(6,*)
      read(12,1000)title
      write(6,1000)title
      write(24,1000)title
      write(6,*)
      write(6,*)'*********extrema of the Cl(alpha) function:'
      km=1
      prod=1.
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
      write(6,*)'*************wing polar (profile data from Xfoil):'
      write(6,*)'  k=      inc=        cz=         cx=         cm=   '
      write(24,*)
      write(24,*)'*************wing polar (profile data from Xfoil):'
      write(24,*)'  k=      inc=        cz=         cx=         cm=   '
      do 4 k=1,kx
         kp=k+1
         if(kp.gt.kx)kp=kx
         km=k-1
         if(km.lt.1)km=1
         if(k.gt.1.and.k.lt.kx)then
         dcxm=((cx(k)-cx(km))*(cz(kp)-cz(k))*(cz(k)+cz(kp)-2.*cz(km))
     &        -(cx(kp)-cx(k))*(cz(k)-cz(km))**2)
     &        /((cz(kp)-cz(k))*(cz(kp)-cz(km))*(cz(k)-cz(km)))
         dcxp=((cx(k)-cx(kp))*(cz(km)-cz(k))*(cz(k)+cz(km)-2.*cz(kp))
     &        -(cx(km)-cx(k))*(cz(k)-cz(kp))**2)
     &        /((cz(km)-cz(k))*(cz(km)-cz(kp))*(cz(k)-cz(kp)))
         endif
         prod=dcxm*dcxp
         if(prod.lt.-eps.and.(kxtrm(km).ne.0.or.kxtrm(kp).ne.0))then
            write(6,*)'bad data distribution: interpolate a new data'
         endif
         write(13,*)cx(k),cz(k),cq(k),inc(k)
         incd=inc(k)
         inc(k)=degrad*inc(k)
         if(inc(km)*inc(k).le.0.0.and.k.ne.1)then
            k0=km
c*****zero incidence coefficients with assumption of small angles
c     lift coefficients of main wing
            dCldam0=(cz(k)-cz(km))/(inc(k)-inc(km))
            Clm00=cz(km)+dCldam0*(-inc(km))
c     add flap influence on Clm
            Clm00=Clm00+dClmdtf*tf
c     moment coefficients of the main wing Cmacm and Cmmo
            dCmacdam0=(cq(k)-cq(km))/(inc(k)-inc(km))
            Cmacm00=cq(km)+dCmacdam0*(-inc(km))
            dCmdam0=dCmacdam0-xacm*dCldam0/cam
            Cmm00=Cmacm00-xacm*Clm00/cam
c     add flap influence on Cmm
            Cmacm00=Cmacm00+dCmmdtf*tf
            Cmm00=Cmm00+dCmmdtf*tf
c     lift coefficients of single canard
            dCldac0=2.0*pi/(1.0+2.0/arc)
            Clc00=dCldac0*tc
c     moment coefficients of single canard Cmacc and Cmco
            dCmacdac0=0.
            Cmacc00=0.
            dCmdac0=dCmacdac0-xacc*dCldac0/cac
            Cmc00=Cmacc00-xacc*Clc00/cac
c     global coefficients wing+canards (ac is area of 2 canards)
c     global lift
            dClda=(am*dCldam0+ac*dCldac0)/aref
            Cl0=(am*Clm00+ac*Clc00)/aref
c     global moments
            dCmacda=(am*cam*dCmacdam0+ac*cac*dCmacdac0)/(aref*lref)
            Cmac0=(am*cam*Cmacm00+ac*cac*Cmacc00)/(aref*lref)
            dCmda=(am*cam*dCmdam0+ac*cac*dCmdac0)/(aref*lref)
            Cm0=(am*cam*Cmm00+ac*cac*Cmc00)/(aref*lref)
         endif
         if(k.eq.1.and.inc(1).ge.0.0)then
            k0=1
c     lift coefficients of main wing
            dCldam0=(cz(2)-cz(1))/(inc(2)-inc(1))
            Clm00=cz(1)+dCldam0*(-inc(1))
c     add flap influence on Clm
            Clm00=Clm00+dClmdtf*tf
c     moment coefficients of the main wing Cmacm and Cmmo
            dCmacdam0=(cq(2)-cq(1))/(inc(2)-inc(1))
            Cmacm00=cq(1)+dCmacdam0*(-inc(1))
            dCmdam0=dCmacdam0-xacm*dCldam0/cam
            Cmm00=Cmacm00-xacm*Clm00/cam
c     add flap influence on Cmc
            Cmac00=Cmac00+dCmmdtf*tf
            Cmm00=Cmm00+dCmmdtf*tf
c     lift coefficient of single canard
            dCldac0=2.0*pi/(1.0+2.0/arc)
            Clc00=dCldac0*tc
c     moment coefficients of single canard Cmacc and Cmco
            dCmacdac0=0.
            Cmacc00=0.
            dCmdac0=dCmacdac0-xacc*dCldac0/cac
            Cmc00=Cmacc00-xacc*Clc00/cac
c     global coefficients wing+canards (ac is area of 2 canards)
c     global lift
            dClda=(am*dCldam0+ac*dCldac0)/aref
            Cl0=(am*Clm00+ac*Clc00)/aref
c     global moments
            dCmacda=(am*cam*dCmacdam0+ac*cac*dCmacdac0)/(aref*lref)
            Cmac0=(am*cam*Cmacm00+ac*cac*Cmacc00)/(aref*lref)
            dCmda=(am*cam*dCmdam0+ac*cac*dCmdac0)/(aref*lref)
            Cm0=(am*cam*Cmm00+ac*cac*Cmc00)/(aref*lref)
         endif
c     aerodynamic center of airplane
         xac=(am*dClmda0*xacm+ac*dClcda0*xacc)
     &        /(am*dClmda0+ac*dClcda0)
         write(14,*)cx(k),cz(k),cq(k),incd
         write(6,1001)k,incd,cz(k),cx(k),cq(k)
 4    continue
      write(6,*)
      write(6,*)'******extrema pointer:'
      write(6,*)'kxtrm(k)=',(kxtrm(k),k=1,kx)
c*****global coefficients
      write(6,*)'**********************************************'
      write(6,*)'zero incidence aerodynamic coefficients at k0=',k0
      write(6,*)'**************lift:'
      write(6,*)' wing:     dCldam0=',dCldam0
     &     ,'   Clm00=',Clm00
      write(6,*)'canar:     dCldac0=',dCldac0
     &     ,'   Clc00=',Clc00
      write(6,*)'total:       dClda=',dClda
     &     ,'     Cl0=',Cl0
      write(6,*)'*****moment at xac:'
      write(6,*)' wing:   dCmacdam0=',dCmacdam0
     &     ,' Cmacm00=',Cmacm00
      write(6,*)'canar:   dCmacdac0=',dCmacdac0
     &     ,' Cmacc00=',Cmacc00
      write(6,*)'total:     dCmacda=',dCmacda
     &     ,'   Cmac0=',Cmac0
      write(6,*)'*****moment at xcg:'
      write(6,*)' wing:     dCmdam0=',dCmdam0
     &     ,'   Cmm00=',Cmm00
      write(6,*)'canar:     dCmdac0=',dCmdac0
     &     ,'   Cmc00=',Cmc00
      write(6,*)'total:       dCmda=',dCmda
     &     ,'     Cm0=',Cm0
      write(6,*)'**estimate aerodynamic center location:'
      write(6,*)'               xac=',xac
     &     ,'(m)  xcg=',xcg,'(m)'
      write(24,*)'**********************************************'
      write(24,*)'zero incidence aerodynamic coefficients at k0=',k0
      write(24,*)'**************lift:'
      write(24,*)' wing:     dCldam0=',dCldam0
     &     ,'   Clm00=',Clm00
      write(24,*)'canar:     dCldac0=',dCldac0
     &     ,'   Clc00=',Clc00
      write(24,*)'total:       dClda=',dClda
     &     ,'     Cl0=',Cl0
      write(24,*)'*****moment at xac:'
      write(24,*)' wing:   dCmacdam0=',dCmacdam0
     &     ,' Cmacm00=',Cmacm00
      write(24,*)'canar:   dCmacdac0=',dCmacdac0
     &     ,' Cmacc00=',Cmacc00
      write(24,*)'total:     dCmacda=',dCmacda
     &     ,'   Cmac0=',Cmac0
      write(24,*)'*****moment at xcg:'
      write(24,*)' wing:     dCmdam0=',dCmdam0
     &     ,'   Cmm00=',Cmm00
      write(24,*)'canar:     dCmdac0=',dCmdac0
     &     ,'   Cmc00=',Cmc00
      write(24,*)'total:       dCmda=',dCmda
     &     ,'     Cm0=',Cm0
      write(24,*)'**estimate aerodynamic center location:'
      write(24,*)'               xac=',xac
     &     ,'(m)  xcg=',xcg,'(m)'
c     center of gravity position if static margin >0
      if(statmarg.gt.0.)then
         xcg=xac-statmarg*lref
      endif
c*****linear model results
      aleq=-(Cmac0+(xcg-xac)*Cl0/lref)
     &     /(dCmacda+(xcg-xac)*dClda/lref)
      Cleq=Cl0+dClda*aleq
      Cweq=Cleq
      if(Cleq.le.0.)then
         write(6,*)'#################no linear solution with Cleq<0'
      endif
      Ueq=sqrt(2.*mg/(rho*aref*Cleq))
      reyeq=Ref*Ueq/100.
      dynaref=0.5*rho*Ueq**2*aref
      call thrust(ndatx,ndat,Ueq,vr,tr,T,dTdv)
      Cteq=rcor*(T+dTdv*Ueq)/dynaref
c     search for point on wing polar and interpolate Cd
         m=1
         prod=aleq+0.5*pi
         do 5 k=2,kx-1
            prod=prod*(aleq-inc(k))
            if(prod.le.-eps)goto 6
            prod=1.
            m=k
 5       continue
 6       continue
         Cdeq=cx(m)+(aleq-inc(m))*(cx(m+1)-cx(m))
     &        /(inc(m+1)-inc(m))
c     wing induced drag estimate
         Cdmeq=Cdeq
         Clmeq=Clm0+dCldam*aleq
         Cdim=Clmeq**2/(pi*em*arm)
         Cdm0=Cdmeq-Cdim
c     calculate fuselage friction drag
         reyf=reyeq*lf/cam
         Cdf0=1.328/sqrt(reyf)
         if(reyf.gt.1.0E5)then
            Cdf0=.072/reyf**.2
         endif
         Cdfeq=Cdf0
c     calculate canard induced drag and friction drag (2 elements)
         Clceq=Clc0+dCldac*aleq
         Cdic=2.0*Clceq**2/(pi*ec*arc)
         reyc=reyeq*cac/cam
         Cdc0=1.328/sqrt(reyc)
         if(reyc.gt.1.0E5)then
            Cdc0=.072/reyc**.2
         endif
         Cdceq=Cdic+2.0*Cdc0
c     total drag
         Cdeq=(am*Cdmeq+af*Cdfeq+ac*Cdceq)/aref
         Cdeq=Cdeq+Cdbrake
         beteq=(Cteq-Cdeq)/Cweq
         Cmac=Cmac0+dCmacda*aleq
         Cmeq=Cm0+dCmda*aleq
         bb(1)=-(Cmeq+xcg*Cweq*cos(aleq+beteq)/lref)
         bb(2)=-(Cleq-Cweq*cos(beteq)+Cteq*sin(aleq))
         bb(3)=-(Cdeq+Cweq*sin(beteq)-Cteq*cos(aleq))
         write(6,1002)
         write(6,*)'****linear model results: m=',m
         write(6,*)'aleq=',aleq,' Ueq=',Ueq,'beteq=',beteq
         write(6,*),'Cleq=',Cleq,'Cdeq=',Cdeq,' Cteq=',Cteq
     &        ,'CM,ac=',Cmac
         write(24,1002)
         write(24,*)'aleq=',aleq,' Ueq=',Ueq,'beteq=',beteq
         write(24,*),'Cleq=',Cleq,'Cdeq=',Cdeq,' Cteq=',Cteq
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
      write(6,1002)
      write(6,*)'begin non-linear equilibrium algorithm' 
      write(24,1002)
      write(24,*)'begin non-linear equilibrium algorithm' 
c*****read canar setting angle tcd (deg)
      write(6,*)'tcd=? (<-99 exit)'
      write(24,*)'tcd=? (<-99 exit)'
      read(5,*)tcd
      write(6,*)'canar tcd=',tcd,'flap tfd=',tfd
      write(24,*)'canar tcd=',tcd,'flap tfd=',tfd
      if(tcd.lt.-99.)goto 200
      tc=degrad*tcd
c*****update with best design input data
      aleq=aleq0
      Ueq=Ueq0
      beteq=beteq0
      Cleq=Cleq0
      Cdeq=Cdeq0
      Cmac=Cmac0
      Cteq=Cdeq
      Cweq=Cleq
      Cmeq=Cmac-xac*Cleq*cos(aleq)/lref
      Cmcg=Cmeq+xcg*Cweq*cos(aleq+beteq)/lref
c     aerodynamic center moment coefficient from design data:
      write(6,*)'********************************************'
      write(6,*)'xcg, xac (xcg calculated or given):'
      write(6,*)'center of gravity:       xcg=',xcg,' (m)'
      write(6,*)'aerodynamic center:      xac=',xac,' (m)'
      write(6,1003)amarg,lref
      write(24,*)'********************************************'
      write(24,*)' xcg, xac (xcg calculated or given):'
      write(24,1002)
      write(24,*)'center of gravity:       xcg=',xcg,' (m)'
      write(24,*)'aerodynamic center:      xac=',xac,' (m)'
      write(24,1003)amarg,lref
c*****equilibrium
c     alpha at equilibrium, lift and moment
      iter=0
c*****iterations
      do 9 it=1,itx
c     engine operating point
         call thrust(ndatx,ndat,Ueq,vr,tr,T,dTdv)
         dynaref=0.5*rho*aref*Ueq**2
         Cteq=rcor*(T+dTdv*Ueq)/dynaref
         Cweq=mg/dynaref
         Reyeq=Ref*Ueq/100.
c     lift coefficient of canar
         dCldac=2.0*pi/(1.0+2.0/arc)
         Clc0=dCldac*tc
         Clceq=Clc0+dCldac*aleq
c     moment coefficients of single canard Cmacc and Cmco
         dCmacdac=0.
         Cmacc0=0.
         dCmdac=dCmacdac-xacc*dCldac/cac
         Cmc0=Cmacc0-xacc*Clc0/cac
c     search for point on polar of wing
         m=1
         prod=aleq+.5*pi
         do 7 k=2,kx-1
            prod=prod*(aleq-inc(k))
            if(prod.le.eps)goto 8
            prod=1.
            m=k
 7       continue
 8       continue
c     lift coefficients of main wing
         dCldam=(cz(m+1)-cz(m))/(inc(m+1)-inc(m))
         Clm0=cz(m)+dCldam*(-inc(m))
c     add canar influence on Clm
         Clm0=Clm0+acw*dCldalind*Clceq/(pi*arc)
c     add flap influence on Clm
         Clm0=Clm0+dClmdtf*tf
c     moment coefficients of the main wing Cmacm and Cmmo
         dCmacdam=(cq(m+1)-cq(m))/(inc(m+1)-inc(m))
         Cmacm0=cq(m)+dCmacdam*(-inc(m))
         dCmdam=dCmacdam-xacm*dCldam/cam
         Cmm0=Cmac0-xacm*Clm0/cam
c     add flap influence on Cmm
         Cmac0=Cmac0+dCmmdtf*tf
         Cmm0=Cmm0+dCmmdtf*tf
c     add influence of wing on canard (downwash)
         Clc0=Clc0+dCldac*awc*Clm0/(pi*arm)
c     global coefficients wing+canard (ac is area of 2 canards)
c     global lift
         dClda=(am*dCldam+ac*dCldac)/aref
         Cl0=(am*Clm0+ac*Clc0)/aref
c     global moments (ac is area of 2 canards)
         dCmacda=(am*cam*dCmacdam+ac*cac*dCmacdac)/(aref*lref)
         Cmac0=(am*cam*Cmacm0+ac*cac*Cmacc0)/(aref*lref)
         dCmda=(am*cam*dCmdam+ac*cac*dCmdac)/(aref*lref)
         Cm0=(am*cam*Cmm0+ac*cac*Cmc0)/(aref*lref)
c     main wing drag
         Cdmeq=cx(m)+(aleq-inc(m))*(cx(m+1)-cx(m))/(inc(m+1)-inc(m))
c     add drag of flap
         Cdmeq=Cdmeq+(dCdtf0+dCdtf1*aleq+dCdtf2*aleq**2)*tf/0.1744
c     wing induced drag estimation
         Clmeq=Clm0+dCldam*aleq
         Cdim=Clmeq**2/(pi*em*arm)
         Cdm0=Cdmeq-Cdim
c     calculate fuselage friction drag
         reyf=reyeq*lf/cam
         Cdf0=1.328/sqrt(reyf)
         if(reyf.gt.1.0E5)then
            Cdf0=.072/reyf**.2
         endif
         Cdfeq=Cdf0
c     calculate canard induced drag and friction drag
         Clceq=Clc0+dCldac*aleq
         Cdic=2.0*Clceq**2/(pi*ec*arc)
         reyc=reyeq*cac/cam
         Cdc0=1.328/sqrt(reyc)
         if(reyc.gt.1.0E5)then
            Cdc0=.072/reyc**.2
         endif
         Cdceq=Cdic+2.0*Cdc0
c     global lift and moment at x=0
         Cleq=Cl0+dClda*aleq
         Cmeq=Cm0+dCmda*aleq
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
         Cdeq=(am*Cdmeq+af*Cdfeq+ac*Cdceq+ar*Cdreq+an*Cdneq)/aref
         Cdeq=Cdeq+Cdbrake
c     relaxation method equation by equation
c     equation 1 for daleq with correction for influence of canard on wing
         daleq=-(Cmeq+xcg*Cweq*cos(aleq+beteq)/lref-zeng*Cteq/lref)
     &        /(dCmda-xcg*Cweq*sin(aleq+beteq)/lref)
         aleq=aleq+omega*daleq
c         daleq=0.
c     update global coefficients
         aleqd=radeg*aleq
         Cleq=Cl0+dClda*aleq
         Cmeq=Cm0+dCmda*aleq
         Cmac=Cmac0+dCmacda*aleq
         Cmcg=Cmeq+xcg*Cweq*cos(aleq+beteq)/lref
c     equation 2 for dUeq
         dUeq=-(Cleq-Cweq*cos(beteq)+Cteq*sin(aleq))
     &        /(-2.*(0.*Cleq-Cweq*cos(beteq)+0.*Cteq*sin(aleq))/Ueq
     &        +0.*dTdv*sin(aleq)/dynaref)
c         dUeq=0.
c     velocity and Reynolds number at equilibrium
         Ueq=Ueq+omega*dUeq
         reyeq=rho*Ueq*cam/amu
c     wing Reynolds number, wing induced drag estimation
         reym=reyeq
         Clmeq=Clm0+dCldam*aleq
         Cdim=Clmeq**2/(pi*em*arm)
         Cdm0=Cdmeq-Cdim
c     fuselage viscous drag
         reyf=reyeq*lf/cam
         Cdf0=1.328/sqrt(reyf)
         if(reyf.gt.1.0E5)then
            Cdf0=.072/reyf**.2
         endif
         Cdfeq=Cdf0
c     canard viscous drag
         reyc=reyeq*cac/cam
         Cdc0=1.328/sqrt(reyc)
         if(reyc.gt.1.0E5)then
            Cdc0=.072/reyc**.2
         endif
         Cdceq=Cdic+2.0*Cdc0
c     rudder viscous drag for 2 sides
         reyr=reyeq*rav/cam
         Cdr0=1.328/sqrt(reyr)
         if(reyr.gt.1.0E5)then
            Cdr0=.072/reyr**.2
         endif
         Cdreq=2.0*Cdr0
c     if winglets are used as rudder, double contribution
c         Cdreq=4.0*Cdr0
c     nacelles viscous drag for 2 engines
         reyn=reyeq*lna/cam
         Cdn0=1.328/sqrt(reyn)
         if(reyn.gt.1.0E5)then
            Cdn0=.072/reyn**.2
         endif
         Cdneq=2.0*Cdn0
c     -total drag
         Cdeq=(am*Cdmeq+af*Cdfeq+ac*Cdceq+ar*Cdreq+an*Cdneq)/aref
         Cdeq=Cdeq+Cdbrake
c     equation 3 for dbeteq
         dbeteq=-(Cdeq+Cweq*sin(beteq)-Cteq*cos(aleq))
     &        /(Cweq*cos(beteq))
c         dbeteq=0.
c     angle of descent
         beteq=beteq+omega*dbeteq
         beteqd=radeg*beteq
c     evaluation of moments
         Cmeq=Cm0+dCmda*aleq
         Cmac=Cmac0+dCmacda*aleq
         inewton=0
         det=1.0
         usdet=1./det
c     residuals
         bb(1)=-(Cmeq+xcg*Cweq*cos(aleq+beteq)/lref-zeng*Cteq/lref)
         bb(2)=-(Cleq-Cweq*cos(beteq)+Cteq*sin(aleq))
         bb(3)=-(Cdeq+Cweq*sin(beteq)-Cteq*cos(aleq))
         reseq=bb(1)**2+bb(2)**2+bb(3)**2
         reseq=sqrt(reseq)
         iter=iter+1
         reseq0=reseq
         if(reseq0.lt.10.*eps)then
c     Newton's Method
            inewton=1
            b1=bb(1)
            b2=bb(2)
            b3=bb(3)
            aa(1,1)=dCmda-xcg*Cweq*sin(aleq+beteq)/lref-zeng*Cteq/lref
            aa(1,2)=2.*bb(1)/Ueq+zeng*dTdv/(dynaref*lref)
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
            det=1.0/usdet
         endif
c     engine operating point
         call thrust(ndatx,ndat,Ueq,vr,tr,T,dTdv)
         write(17,*)iter,reseq,inewton
         if(reseq.lt.epser)goto 10
 9       continue
 10   continue
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
     &        ,' dUeq=',bb(2),'dbeteq=',bb(3)
         write(6,*)' equ1=',b1
     &        ,' equ2=',b2,'  equ3=',b3
      endif
c     results
      theq=aleq+beteq
      theqd=aleqd+beteqd
      weq=Ueq*sin(beteq)
      dyn=0.5*rho*Ueq**2
c     finess coefficient (Cl/Cd)
      ClCdeq=Cleq/Cdeq
c     wing coefficients
      winglift=am*dynaref*Cleq/aref
      Cmacm=dCmacdam*aleq+Cmacm0
      xcpm=xacm-cam*Cmacm/Clmeq
c     canard coefficients
      Clceq=Clc0+dCldac*aleq
      Cdic=2.*Clceq**2/(pi*ec*arc)
      Cdceq=Cdic+2.0*Cdc0
      canarlift=ac*Clceq*dynaref/am
      Cmacc=dCmacdac*aleq+Cmacc0
      xcpc=xacc-cac*Cmacc/Clceq
      reyc=reyeq*cac/cam
      tr0=rcor*(T+dTdv*Ueq)
      write(6,*)'*****************************'
      write(6,*)'        canard setting angle=',tc,' (rd) ='
     &     ,tcd,' (deg)'
      if(it.ge.itx)then
         write(6,*)'equilibrium solution:   iter=',iter
     &     ,'      error=',reseq
         write(6,*)' ',nc
      else
         if(aleq.gt.inc(1))then
            write(6,*)'equilibrium solution:   iter=',iter
     &     ,'      error=',reseq
         else
            write(6,*)'equilibrium solution:   iter=',iter
     &     ,'      error=',reseq
            write(6,*)' ',yw
         endif
      endif
      write(6,*)'global results:              '
     &     ,'         (reference area=',aref,' (m**2))'
      write(6,*)'                   incidence=',aleq,' (rd) ='
     &     ,aleqd,' (deg)'
      write(6,*)'                 climb slope=',beteq,' (rd) ='
     &     ,beteqd,' (deg)'
      write(6,*)'            airplane setting=',theq,' (rd) ='
     &     ,theqd,' (deg)'
      write(6,*)'                    velocity=',Ueq,' (m/s)'
      write(6,*)'                  climb rate=',weq,' (m/s)'
      write(6,*)'            dynamic pressure=',dyn,' (Pa)'
      write(6,*)'       dynamic pressure*aref=',dynaref,' (N)'
      write(6,*)'             Reynolds number=',reyeq
      write(6,*)'                     lift CL=',Cleq
      write(6,*)'       weight coefficient CW=',Cweq
      write(6,*)'                moment CM,ac=',Cmac
      write(6,*)'                     drag CD=',Cdeq
      write(6,*)'       thrust coefficient CT=',Cteq
      write(6,*)'                     Cdbrake=',Cdbrake
      write(6,*)'                       CL/CD=',ClCdeq
      write(6,*)'main wing:                   '
     &     ,'          (ref/wing area=',am,' (m**2))'
      write(6,*)'                    lift CLm=',Clmeq
     &     ,'    Lm=',winglift,' (N)'
      write(6,*)'  wing pitching moment CMacm=',Cmacm
      write(6,*)'               wing Reynolds=',reym
      write(6,*)' estimated induced drag CDim=',Cdim
      write(6,*)'                    drag CDm=',Cdmeq
      write(6,*)'                        xcpm=',xcpm,' (m)'
      write(6,*)'canard:                      '
     &     ,'          (ref/can. area=',ac,' (m**2))'
      write(6,*)'                    lift CLc=',Clceq
     &     ,'    Lc=',canarlift,' (N)'
      write(6,*)'canard pitching moment CMacc=',Cmacc
      write(6,*)'             canard Reynolds=',reyc
      write(6,*)' estimated induced drag CDic=',Cdic
      write(6,*)'                    drag CDc=',Cdceq
      write(6,*)'                        xcpc=',xcpc,' (m)'
      write(6,*)'fuselage:                    '
     &     ,'          (ref/fus. area=',af,' (m**2))'
      write(6,*)'           fuselage Reynolds=',reyf
      write(6,*)'                    drag CDf=',Cdfeq
      write(6,*)'rudder:                      '
     &     ,'          (ref/rud. area=',ac,' (m**2))'
      write(6,*)'             rudder Reynolds=',reyr
      write(6,*)'                    drag Cdr=',Cdreq
      write(6,*)'nacelle:                     '
     &     ,'          (ref/wet. area=',an,' (m**2))'
      write(6,*)'            nacelle Reynolds=',reyn
      write(6,*)'                    drag Cdn=',Cdneq
      write(6,1005)Pcent,tr0
      if(it.lt.itx)then
         imarg=amarg+0.1
         itfd=tfd+0.1*sign(1.0,tfd)
         write(25,*)imarg,itfd,aleq,beteq,Ueq,tr0
      endif
      write(24,*)'*****************************'
      write(24,*)'          flap setting angle=',tf,' (rd) ='
     &     ,tfd,' (deg)'
      if(it.ge.itx)then
         write(24,*)'equilibrium solution:   iter=',iter
     &     ,' error=',reseq,' ',nc
      else
         write(24,*)'equilibrium solution:   iter=',iter
     &     ,' error=',reseq
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
      write(24,*)'       dynamic pressure*aref=',dynaref,' (N)'
      write(24,*)'             Reynolds number=',reyeq
      write(24,*)'                     lift CL=',Cleq
      write(24,*)'       weight coefficient CW=',Cweq
      write(24,*)'                moment CM,ac=',Cmac
      write(24,*)'                     drag CD=',Cdeq
      write(24,*)'       thrust coefficient CT=',Cteq
      write(24,*)'                     Cdbrake=',Cdbrake
      write(24,*)'                       CL/CD=',ClCdeq
      write(24,*)'main wing:                   '
     &     ,'         (ref/wing area=',am,' (m**2))'
      write(24,*)'                    lift CLm=',Clmeq
      write(24,*)'  wing pitching moment CMacm=',Cmacm
      write(24,*)'               wing Reynolds=',reym
      write(24,*)' estimated induced drag CDim=',Cdim
      write(24,*)'                    drag CDm=',Cdmeq
      write(24,*)'                        xcpm=',xcpm,' (m)'
      write(24,*)'canard:                       '
     &     ,'        (ref/canard area=',ac,' (m**2))'
      write(24,*)'                    lift CLc=',Clceq
     &     ,'    Lc=',canarlift,' (N)'
      write(24,*)'canard pitching moment CMacc=',Cmacc
      write(24,*)'             canard Reynolds=',reyc
      write(24,*)' estimated induced drag CDic=',Cdic
      write(24,*)'                    drag CDc=',Cdceq
      write(24,*)'                        xcpc=',xcpc,'(m)'
      write(24,*)'fuselage:                    '
     &     ,'         (ref/fus. area=',af,' (m**2))'
      write(24,*)'           fuselage Reynolds=',reyf
      write(24,*)'                    drag CDf=',Cdfeq
      write(24,*)'rudder:'
      write(24,*)'      (ref/rud. area=',ac,' (m**2)'
      write(24,*)'             rudder Reynolds=',reyr
      write(24,*)'                    drag Cdr=',Cdreq
      write(24,*)'nacelle:                     '
     &     ,'          (ref/wet. area=',an,' (m**2))'
      write(24,*)'            nacelle Reynolds=',reyn
      write(24,*)'                    drag Cdn=',Cdneq
      write(24,1005)Pcent,tr0
      goto 100
 200  continue
      tfd=tf*radeg
      write(21,*)xcg,tfd,aleqd,beteqd,Ueq,Clmeq
c*****loop on incidence
      nalmin=radeg*inc(1)
      nalmax=nalmin+nal
      do 13 n=nalmin,nalmax
         alphd=n
         alph=degrad*alphd
c*****lift drag and moment of wing
c     search for point on polar of wing
         m=1
         prod=alph+.5*pi
         do 11 k=2,kx-1
            prod=prod*(alph-inc(k))
            if(prod.le.eps)goto 12
            prod=1.
            m=k
 11      continue
 12      continue
c     main wing lift
         dCldam=(cz(m+1)-cz(m))/(inc(m+1)-inc(m))
         Clm0=cz(m)+dCldam*(-inc(m))
c     add flap influence on Cl
         Clm0=Clm0+dClmdtf*tf
c     main wing moment coefficients
         Clmeq=Clm0+dCldam*aleq
         dCmacda=(cq(m+1)-cq(m))/(inc(m+1)-inc(m))
         dCmdam=dCmacda-xac*dCldam/cam
         Cmac0=cq(m)+dCmacda*(-inc(m))
         Cmm0=Cmac0-xac*Clm0/cam
c     add flap influence on Cm
         Cmac0=Cmac0+dCmmdtf*tf
         Cmm0=Cmm0+dCmmdtf*tf
c     moments at ac and x=0
         Cmac=Cmac0+dCmacda*aleq
         Cmeq=Cm0+dCmda*aleq
         if(kxtrm(m).ne.0)m=m+1
         if(m.lt.2)m=2
         if(m.gt.kx-1)m=kx-1
         dcxdcz=(cx(m)-cx(m-1))/(cz(m)-cz(m-1))
         d2cxdcz2=((cx(m+1)-cx(m))*(cz(m)-cz(m-1))
     &        -(cx(m)-cx(m-1))*(cz(m+1)-cz(m)))
     &     /((cz(m+1)-cz(m))*(cz(m+1)-cz(m-1))*(cz(m)-cz(m-1)))
         Cdm=cx(m-1)+dcxdcz*(Clmeq-cz(m-1))
     &        +d2cxdcz2*(Clmeq-cz(m-1))*(Clmeq-cz(m))
         dClda=dCldam
         Cl0=Clm0
         dCmda=dCmdam
         Cm0=Cmm0
c*****global lift drag and moment
         Cl=Cl0+dClda*alph
         Cm=Cm0+dCmda*alph
         Cd=(am*Cdm+af*Cdf0)/aref
         Cl3=Cl**3
         Cd2=Cd**2
         ratio1=Cl/Cd
         ratio2=Cl3/Cd2
         incidence=n
         write(16,*)Cdm,Clmeq,Cmmeq,incidence
         write(18,*)Cd,Cl,Cm,incidence
         write(19,*)Cl,Cd,ratio1,Cd2,Cl3,ratio2,incidence
         Cmcg=Cmeq+xcg*Cweq*cos(aleq+beteq)/lref
         write(20,*)incidence,Cmcg
 13   continue
      Cmmeq=Cmeq
      write(22,*)Cdmeq,Clmeq,Cmmeq,aleqd
      write(23,*)Cdeq,Cleq,Cmeq,aleqd,Cmcg
c*****files
      write(6,*)
      write(6,*)'******input file'
      write(6,*)'canarpolar.dat       :polar from prandtline & XFOIL'
     &     ,' (optional)'
      write(6,*)'canareq.data         :parameters file'
      write(6,*)'******output files'
      write(6,*)'canarpolar.prandtline:CDm,CLm,CMm,alpha prandtline'
      write(6,*)'canarpolar.cdclcq    :CDm,CLm,CMm,alpha'
      write(6,*)'canareq.cdclcqm      :CDm,CLm,CMm,alpha'
      write(6,*)'canareq.cdclcq       :CD,CL,CM,alpha'
      write(6,*)'canareq.cdmcln       :CD,CL,CL/CD,CD**2,CL**3,CL3/CD2'
     &     ,',alpha'
      write(6,*)'canareq.itres        :iter,alog10(res)'
      write(6,*)'canareq.cmcg         :alpha,CMcg'
      write(6,*)'canareq.xtabvcl      :xcg,tf,aleq,beteq,Ueq,Clmeq'
      write(6,*)'canareq.cdclcqmeq    :CDmeq,CLmeq,CMmeq,aleq'
      write(6,*)'canareq.cdclcqcmcgeq :CDeq,CLeq,CMeq,aleq,CMcgeq'
      write(6,*)'canareq.list         :listing of results'
      write(6,*)'canareq.stab         :imarg,itfd,aleq,beteq,Ueq,tr0'
 1000 format(18a4)
 1001 format(1x,i4,4x,f8.4,4x,f8.4,4x,f8.4,4x,f8.4)
 1002 format(' #####################################################',
     &     '################')
 1003 format('static margin=',f5.1,' %'
     &     ,'                  (reference lref=',f12.7,'       (m))')
 1004 format(3x,i3,4x,f8.4,4x,e12.4)
 1005 format(' Pcent=',f8.0,'        thrust=',f12.4,'       (N)')
      end

      subroutine thrust(ndatx,ndat,v,vr,tr,T,dTdv)
      implicit none
      integer ndatx,ndat,i,im
      real v,T,dTdv,prod
      real vr(ndatx),tr(ndatx)
      if(abs(v).lt.1.0)then
         dTdv=(tr(2)-tr(1))/(vr(2)-vr(1))
         T=tr(1)+(v-vr(1))*dTdv
         goto 2
      endif
      prod=v-vr(1)
      do 1 i=2,ndat
         prod=prod*(v-vr(i))
         if(prod.le.0.)then
            dTdv=(tr(i)-tr(i-1))/(vr(i)-vr(i-1))
            T=tr(i-1)+dTdv*(-vr(i-1))
            goto 2
         endif
         prod=1.
         im=i
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


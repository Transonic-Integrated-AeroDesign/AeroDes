      program ould
c*****Optimization Under Lift Deformation
      implicit none
      integer jxx,lxx,nxx,ipolar,nx,n,km,kfirst,k,kdum,kp
      integer jx,index,jl,jr,jx2,is,j0,j,jm,iread,ipr,ivis
      integer itx,mj,iter,it,jdx,jdz,jp,ial,ice,iop,iper
      integer jcxw,int
      parameter(jxx=501,lxx=101,nxx=3)
      real eps,pi,radeg,degrad,prod,finex,cloptn,duminc,dihedr
      real czdum,cxdum,dum,cqdum,fines,dcz,dcxm,dcxp,at,bt,al
      real omegades,omegana,avis,cltar,cloptw,eix,cxm,am,dm,dz,g0
      real alphad,cxw,aw,dw,yawd,deltald,deltard,Lref,Vref,Rho,Amu
      real aream,arm,areaw,arw,alpha,yaw,deltal,deltar,reynolds
      real si,dt,tj,yj,zj,etaj,zetaj,ztip,smid,reym,cdm0,vis,alfj
      real czj,dcxdcz,d2cxdcz2,cxj,cost,alagr,cl,A,B,dgx,deta,dzeta
      real gj,ck,dk,piv,dgj,dzx,zm,dy,dzdy,dztip,a1mean,sumj,ratio
      real cdilw,cdirw,cdim,cdvlw,cdvrw,cdvm,d0,cdi,cdv,cd,gmax
      real cle,cdie,cdve,cde,toein,toeind,alfgeo,alfgeod,cmom,cav
      real dumdum,cnom,cnow,clo,sum1,sum2,aj,visc,eme,clm,em,cdlw
      real cdm,cmolw,cdrw,cmorw,e,xlew,tmj,crw,cdil,cllw,clrw,cno
      real gkink,wi
      real y(jxx),z(jxx),eta(jxx),zeta(jxx)
      real s(jxx),ds(jxx)
      real g(jxx),un(jxx)
      real ge(jxx),we(jxx)
      real dcdds(jxx),dcddse(jxx)
      real cmf(jxx),dih(jxx)
      real fy(jxx),fz(jxx)
      real cm(jxx),cm0(jxx),dem(jxx),tm(jxx)
      real alf(jxx),d(jxx),l(jxx)
      real a0(jxx),a1(jxx),a2(jxx)
      real b0(jxx),b1(jxx)
      real c0(jxx),c1(jxx)
      real sbreak(nxx),clopt(nxx)
      real cx(lxx,nxx),cz(lxx,nxx),cq(lxx,nxx),inc(lxx,nxx)
      real incd
      integer m(jxx),polar(jxx)
      integer kx(nxx),kxtrm(lxx,nxx),mxtrm(nxx)
      character*4 bry(18)
      character*4 title(18)
      character*4 typecode(18)      
      data m/jxx*2/
      data sbreak/nxx*1./
      character*9 igeo
      common/para/eps,bt,index
      data igeo/' iterate '/
      open(unit=12,file='polarbl.dat',form='formatted')
      open(unit=13,file='polarbl.cdcl',form='formatted')
      open(unit=14,file='ould.in',form='formatted')
      open(unit=15,file='ould.data',form='formatted')
      open(unit=16,file='ould.yz0',form='formatted')
      open(unit=17,file='ould.yz',form='formatted')
      open(unit=18,file='ould.gun',form='formatted')
      open(unit=19,file='ould.gwe',form='formatted')
      open(unit=20,file='ould.cdy',form='formatted')
      open(unit=21,file='ould.ymf',form='formatted')
      open(unit=22,file='ould.yfz',form='formatted')
      open(unit=23,file='ould.cdt',form='formatted')
      open(unit=24,file='ould.guna',form='formatted')
      open(unit=25,file='ould.gwea',form='formatted')
      open(unit=26,file='ould.zxdtw',form='formatted')
      open(unit=27,file='ould.tcw',form='formatted')
      open(unit=28,file='ould.out',form='formatted')
      open(unit=29,file='ould.clcdcq',form='formatted')
      open(unit=30,file='ould.cdcl',form='formatted')
      open(unit=31,file='ould.list',form='formatted')
c*****constants
      eps=0.5e-7
      pi=2.*asin(1.)
      radeg=180./pi
      degrad=1./radeg
c*****polar data
      sbreak(1)=1000.
      write(6,*)
      write(6,*)'******do you want to use polar data? Y/N=1/0'
      write(31,*)
      write(31,*)'******do you want to use polar data? Y/N=1/0'
      read(5,*)ipolar
      if(ipolar.ne.1)goto 4
      read(12,*)nx
      write(6,*)'******profile polars:'
      write(6,*)' nx= ',nx,' number of polars to be read'
      write(31,*)'******profile polars:'
      write(31,*)' nx= ',nx,' number of polars to be read'
      if(nx.gt.nxx)then
         write(6,*)'!! nx > nxx !!'
         write(6,*)'TOO MANY POLARS: EXITING! CHANGE HARDCODE!'
         write(31,*)'!! nx > nxx !!'
         write(31,*)'TOO MANY POLARS: EXITING! CHANGE HARDCODE!'
         stop
      endif
      do 100 n=1,nx     
      write(6,*)'******n= ',n
      write(31,*)'******n= ',n
      read(12,1001)title
      write(6,1001)title
      write(31,1001)title
      read(12,1001)title
      write(6,1001)title
      write(31,1001)title
      read(12,1001)title
      write(6,1001)title
      write(31,1001)title
      read(12,1001)title
      write(6,1001)title
      write(31,1001)title
      read(12,1001)title
      write(6,1001)title
      write(31,1001)title
      read(12,1001)title
      write(6,1001)title
      write(31,1001)title
      read(12,1001)title
      write(6,1001)title
      write(31,1001)title
      read(12,1001)title
      write(6,1001)title
      write(31,1001)title
      read(12,1001)title
      write(6,1001)title
      write(31,1001)title
      read(12,1001)title
      write(6,1001)title
      write(31,1001)title
      read(12,1001)title
      write(6,1001)title
      write(31,1001)title
      read(12,1001)title
      write(6,1001)title
      write(31,1001)title
      write(6,*)
      write(31,*)
      read(5,1001)bry
      write(6,*)
      write(6,*)'*****extrema of the cl(alpha) function:'
      write(31,*)
      write(31,*)'*****extrema of the cl(alpha) function:'
      prod=1.
      km=1
      kfirst=0
      finex=0.
      cloptn=0.
      do 1 k=1,lxx
         kdum=k
         read(12,*,end=2)duminc,czdum,cxdum,dum,cqdum
         inc(k,n)=duminc
         cx(k,n)=cxdum
         cz(k,n)=czdum
         cq(k,n)=cqdum
         fines=czdum/cxdum
         if(fines.gt.finex)then
            finex=fines
            cloptn=czdum
         endif
         if(duminc.gt.89.)then
            read(12,*)sbreak(n)
            if(n.ge.nx)then
               sbreak(n)=1000.
            endif
            kdum=k+1
            goto 2
         endif
         kxtrm(k,n)=0
         if(k.eq.1)goto 1
         dcz=cz(k,n)-cz(km,n)
         prod=prod*dcz
         if(prod.lt.eps)then
            write(6,*)'kxtrm(',n,')=',km,' cz(kxtrm,n)=',cz(km,n)
            write(31,*)'kxtrm(',n,')=',km,' cz(kxtrm,n)=',cz(km,n)
            kxtrm(km,n)=km
            if(kfirst.le.0)then
               mxtrm(n)=km
               kfirst=1
            endif
         endif
         prod=sign(1.,dcz)
         km=k
 1    continue
 2    continue
      kx(n)=kdum-1
      clopt(n)=cloptn
      if(kx(n).eq.lxx-1)then
         write(6,*)' attention: check if all data has been read;'
     &        ,' continuing/exiting=1/0?'
         write(31,*)' attention: check if all data has been read;'
     &        ,' continuing/exiting=1/0?'
         read(5,*)ice
         if(ice.eq.0)stop
      endif
      write(6,*)
      write(6,*)'*************profile data from Xfoil:'
      write(6,1002)
      write(31,*)
      write(31,*)'*************profile data from Xfoil:'
      write(31,1002)
      do 3 k=1,kx(n)
         kp=k+1
         if(kp.gt.kx(n))kp=kx(n)
         km=k-1
         if(km.lt.1)km=1
         prod=1.
         if(k.gt.1.and.k.lt.kx(n))then
         dcxm=((cx(k,n)-cx(km,n))*(cz(kp,n)-cz(k,n))*(cz(k,n)+cz(kp,n)
     &        -2.*cz(km,n))-(cx(kp,n)-cx(k,n))*(cz(k,n)-cz(km,n))**2)/
     &     ((cz(kp,n)-cz(k,n))*(cz(kp,n)-cz(km,n))*(cz(k,n)-cz(km,n)))
         dcxp=((cx(k,n)-cx(kp,n))*(cz(km,n)-cz(k,n))*(cz(k,n)+cz(km,n)
     &        -2.*cz(kp,n))-(cx(km,n)-cx(k,n))*(cz(k,n)-cz(kp,n))**2)/
     &     ((cz(km,n)-cz(k,n))*(cz(km,n)-cz(kp,n))*(cz(k,n)-cz(kp,n)))
         prod=dcxm*dcxp
         endif
         if(prod.lt.-eps.and.(kxtrm(km,n).ne.0.
     &        or.kxtrm(kp,n).ne.0))then
            write(6,*)'bad data distribution:',
     &           ' interpolate a new data n=',n
            write(31,*)'bad data distribution:',
     &           ' interpolate a new data n=',n
         endif
         incd=inc(k,n)
         inc(k,n)=degrad*inc(k,n)
         write(6,1003)k,incd,cz(k,n),cx(k,n),cq(k,n)
         write(31,1003)k,incd,cz(k,n),cx(k,n),cq(k,n)
         write(13,*)cx(k,n),cz(k,n),incd,cq(k,n)
 3    continue
      write(6,*)
      write(6,*)'******extrema pointer:'
      write(6,*)'kxtrm(k,',n,')=',(kxtrm(k,n),k=1,kx(n))
      write(6,*)
      write(31,*)
      write(31,*)'******extrema pointer:'
      write(31,*)'kxtrm(k,',n,')=',(kxtrm(k,n),k=1,kx(n))
      write(31,*)
100   continue
      write(6,*)'clopt=',(clopt(n),n=1,nx)
      write(6,*)'sbreak=',(sbreak(n),n=1,nx)
      write(31,*)'clopt=',(clopt(n),n=1,nx)
      write(31,*)'sbreak=',(sbreak(n),n=1,nx)
 4    continue
c*****data
      read(15,*)jx
      read(15,*)at
      read(15,*)bt
      read(15,*)index
      read(15,*)al
      read(15,*)omegades
      read(15,*)omegana
      read(15,*)avis
      read(15,*)cltar
      read(15,*)cloptw
      read(15,*)eix
      read(15,*)cxm
      read(15,*)am
      read(15,*)dm
      read(15,*)alphad
      read(15,*)cxw
      read(15,*)aw
      read(15,*)dw
      read(15,*)yawd
      read(15,*)deltald
      read(15,*)deltard
      read(15,*)Lref
      read(15,*)Vref
      read(15,*)Rho
      read(15,*)Amu
      aream=2.*am*cxm
      arm=4./aream
      areaw=(abs(at)+eps)*cxw*(aw+eps)
      arw=at**2/areaw
      alpha=alphad/radeg
      yaw=yawd/radeg
      deltal=deltald/radeg
      deltar=deltard/radeg
      reynolds=Rho*Vref*Lref/Amu
      if(ipolar.eq.1.and.cloptw.lt.0.)then
         cloptw=clopt(1)
      endif
      if(jx.gt.jxx) jx=jxx
      if(abs(at).lt.eps)then
         jl=1
         jr=jx
      else
         int=sqrt(2./abs(at))
         jl=(jx+2)/(3+int)
         jr=jx+1-jl
      endif
      eix=cltar*eix
      write(6,*)
      write(6,*)'***********numerical data:'
      write(6,*)' number of mesh points jx=',jx
      write(6,*)'            left plate jl=',jl,
     &     ' right plate jr=',jr
      write(6,*)'     curve exponent index=',index
      write(6,*)'   Lagrange multiplier al=',al
      write(6,*)' relaxation factor omegad=',omegades
      write(6,*)' relaxation factor omegan=',omegana
      write(6,*)'artificial viscosity avis=',avis
      write(6,*)'          target cl cltar=',cltar
      write(6,*)'           winglet cloptw=',cloptw
      write(6,*)
      write(6,*)'*******main wing geometry:'
      write(6,*)'    E*area moment inertia=',eix
      write(6,*)'           span height bt=',bt
      write(6,*)'           root chord cxm=',cxm
      write(6,*)' wing area coefficient am=',am
      write(6,*)'                wing area=',aream,' (ref. B**2/4)'
      write(6,*)'        wing aspect ratio=',arm
      write(6,*)'  main wing prof rel camb=',dm
      write(6,*)'          incidence alpha=',alpha,' (rad)  ='
     &     ,alphad,' (deg)'
      write(6,*)
      write(6,*)'*********winglet geometry:'
      write(6,*)'        winglet height at=',at
      write(6,*)'   winglet root chord cxw=',cxw
      write(6,*)'    winglet area coef. aw=',aw
      write(6,*)'             winglet area=',areaw,' (ref. B**2/4)'
      write(6,*)'     winglet aspect ratio=',arw
      write(6,*)'    winglet prof rel camb=',dw
      write(6,*)'                      yaw=',yaw,' (rad)  ='
     &     ,yawd,' (deg)'
      write(6,*)'                   deltal=',deltal,' (rad)  ='
     &     ,deltald,' (deg)'
      write(6,*)'                   deltar=',deltar,' (rad)  ='
     &     ,deltard,' (deg)'
      write(6,*)
      write(6,*)'***dimensional parameters:'
      write(6,*)'           half span Lref=',Lref,' (m)'
      write(6,*)'            velocity Vref=',Vref,' (m/s)'
      write(6,*)'          air density Rho=',Rho,' (kg/m**3)'
      write(6,*)'        air viscosity Amu=',Amu,' (kg/(m*s))'
      write(6,*)'                 Reynolds=',reynolds
      write(6,*)
      write(6,*)'###############################normalized variables:'
      write(6,*)'                        Y=0.5*B*y'
      write(6,*)'                        Z=0.5*B*z'
      write(6,*)'                        S=0.5*B*s'
      write(6,*)'                        C=0.5*B*cxm'
      write(6,*)'                      Cav=0.5*B*cav'
      write(6,*)'                       AT=0.5*B*at'
      write(6,*)'                       BT=0.5*B*bt'
      write(6,*)'                        G=0.5*U*B*g'
      write(6,*)'                     Aref=0.5*B**2*cxm*am'
     &     ,'=0.25*B**2*aream'
      write(6,*)'                       Un=U*un (at lifting line)'
      write(6,*)'                       CL=cl'
      write(6,*)'                      CDi=cdi'
      write(6,*)'                        L=0.5*RHO*U**2*Aref*cl'
      write(6,*)'                        D=0.5*RHO*U**2*C**2*cdi'
      write(6,*)'                      M,o=0.5*RHO*U**2*Aref*Cav*cmo'
      write(6,*)'                       Fy=0.5*RHO*U**2*B*C*fy'
      write(6,*)'                       Fz=0.5*RHO*U**2*B*C*fz'
      write(6,*)'                       Mf=0.5*RHO*U**2*Aref*B*cmf'
      write(6,*)'                      EIx=0.25*RHO*U**2*B**3*C*eix'
      write(6,*)'                      EIx=0.5*Weight*B**2*eix/cl'
      write(6,*)'                      N,o=0.25*RHO*U**2*Aref*B*cno'
      write(6,*)'                      L,o=0.25*RHO*U**2*Aref*B*clo'
      write(6,*)
      write(31,*)
      write(31,*)'***********numerical data:'
      write(31,*)' number of mesh points jx=',jx
      write(31,*)'            left plate jl=',jl,
     &     ' right plate jr=',jr
      write(31,*)'     curve exponent index=',index
      write(31,*)'   Lagrange multiplier al=',al
      write(31,*)' relaxation factor omegad=',omegades
      write(31,*)' relaxation factor omegan=',omegana
      write(31,*)'artificial viscosity avis=',avis
      write(31,*)'          target cl cltar=',cltar
      write(31,*)'           winglet cloptw=',cloptw
      write(31,*)
      write(31,*)'*******main wing geometry:'
      write(31,*)'    E*area moment inertia=',eix
      write(31,*)'           span height bt=',bt
      write(31,*)'           root chord cxm=',cxm
      write(31,*)' wing area coefficient am=',am
      write(31,*)'                wing area=',aream,' (ref. B**2/4)'
      write(31,*)'        wing aspect ratio=',arm
      write(31,*)'  main wing prof rel camb=',dm
      write(31,*)'          incidence alpha=',alpha,' (rad)  ='
     &     ,alphad,' (deg)'
      write(31,*)
      write(31,*)'*********winglet geometry:'
      write(31,*)'        winglet height at=',at
      write(31,*)'   winglet root chord cxw=',cxw
      write(31,*)'    winglet area coef. aw=',aw
      write(31,*)'             winglet area=',areaw,' (ref. B**2/4)'
      write(31,*)'     winglet aspect ratio=',arw
      write(31,*)'    winglet prof rel camb=',dw
      write(31,*)'                      yaw=',yaw,' (rad)  ='
     &     ,yawd,' (deg)'
      write(31,*)'                   deltal=',deltal,' (rad)  ='
     &     ,deltald,' (deg)'
      write(31,*)'                   deltar=',deltar,' (rad)  ='
     &     ,deltard,' (deg)'
      write(31,*)
      write(31,*)'***dimensional parameters:'
      write(31,*)'           half span Lref=',Lref,' (m)'
      write(31,*)'            velocity Vref=',Vref,' (m/s)'
      write(31,*)'          air density Rho=',Rho,' (kg/m**3)'
      write(31,*)'        air viscosity Amu=',Amu,' (kg/(m*s))'
      write(31,*)'                 Reynolds=',reynolds
      write(31,*)
      write(31,*)'###############################normalized variables:'
      write(31,*)'                        Y=0.5*B*y'
      write(31,*)'                        Z=0.5*B*z'
      write(31,*)'                        S=0.5*B*s'
      write(31,*)'                        C=0.5*B*cxm'
      write(31,*)'                      Cav=0.5*B*cav'
      write(31,*)'                       AT=0.5*B*at'
      write(31,*)'                       BT=0.5*B*bt'
      write(31,*)'                        G=0.5*U*B*g'
      write(31,*)'                     Aref=0.5*B**2*cxm*am'
     &     ,'=0.25*B**2*aream'
      write(31,*)'                       Un=U*un (at lifting line)'
      write(31,*)'                       CL=cl'
      write(31,*)'                      CDi=cdi'
      write(31,*)'                        L=0.5*RHO*U**2*Aref*cl'
      write(31,*)'                        D=0.5*RHO*U**2*C**2*cdi'
      write(31,*)'                      M,o=0.5*RHO*U**2*Aref*Cav*cmo'
      write(31,*)'                       Fy=0.5*RHO*U**2*B*C*fy'
      write(31,*)'                       Fz=0.5*RHO*U**2*B*C*fz'
      write(31,*)'                       Mf=0.5*RHO*U**2*Aref*B*cmf'
      write(31,*)'                      EIx=0.25*RHO*U**2*B**3*C*eix'
      write(31,*)'                      EIx=0.5*Weight*B**2*eix/cl'
      write(31,*)'                      N,o=0.25*RHO*U**2*Aref*B*cno'
      write(31,*)'                      L,o=0.25*RHO*U**2*Aref*B*clo'
      write(31,*)
      jx2=jx/2
      is=mod(jx,2)
      si=is
      j0=jx2+is
c*****cambered span
      dt=pi/(jr-jl)
      tj=-.5*(1-is)*dt
      y(j0)=0.
      z(j0)=0.
      eta(j0)=0.
      zeta(j0)=0.
      do 5 j=j0+1,jr
         tj=tj+dt
         yj=sin(tj)
         y(j)=yj
         y(jx+1-j)=-yj
         zj=dihedr(yj)
         dih(j)=zj
         z(j)=zj
         dih(jx+1-j)=zj
         z(jx+1-j)=zj
         if(is.ne.0.or.j.ne.j0+1)then
            etaj=sin(tj-.5*dt)
            eta(j-1)=etaj
            eta(jx+1-j)=-etaj
            zetaj=dihedr(etaj)
            zeta(j-1)=zetaj
            zeta(jx+1-j)=zetaj
         endif
 5    continue
      ztip=z(jr)
c*****horizontal winglets
      dt=pi/(jx-jr+eps)
      tj=0.
      do 6 j=jr+1,jx
         tj=tj+dt
         yj=1.+.5*abs(at)*(1.-cos(tj))
         y(j)=yj
         y(jx+1-j)=-yj
         zj=ztip
         z(j)=zj
         z(jx+1-j)=zj
         etaj=1.+.5*abs(at)*(1.-cos(tj-.5*dt))
         eta(j-1)=etaj
         eta(jx+1-j)=-etaj
         zetaj=ztip
         zeta(j-1)=zetaj
         zeta(jx+1-j)=zetaj
 6    continue
      jm=1
      s(1)=0.
      do 7 j=1,jx
         alf(j)=0.
         d(j)=0.
         l(j)=0.
         cm(j)=cxm
         dem(j)=dm
         if(j.lt.jl.or.j.gt.jr)then
            dem(j)=dw
         endif
         tm(j)=0.
         g(j)=0.
         un(j)=0.
         ge(j)=0.
         we(j)=0.
         s(j)=s(jm)+sqrt((y(j)-y(jm))**2+(z(j)-z(jm))**2)
         jm=j
         ds(j)=0.
         if(j.gt.1.and.j.lt.jx)then
            ds(j)=sqrt((eta(j)-eta(j-1))**2+(zeta(j)-zeta(j-1))**2)
         endif
 7    continue
      smid=.5*s(jx)
      do 8 j=1,jx
         s(j)=s(j)-smid
 8    continue
c*****vertical winglets
      do 9 j=jr+1,jx
         zj=sign(1.,at)*(y(j)-1.)+ztip
         yj=1.
         y(j)=yj
         y(jx+1-j)=-yj
         z(j)=zj
         dih(j)=zj
         z(jx+1-j)=zj
         dih(jx+1-j)=zj
         zetaj=sign(1.,at)*(eta(j-1)-1.)+ztip
         etaj=1.
         eta(j-1)=etaj
         eta(jx+1-j)=-etaj
         zeta(j-1)=zetaj
         zeta(jx+1-j)=zetaj
 9    continue
      eta(jx)=eta(jx-1)
      zeta(jx)=zeta(jx-1)
      write(6,*)'read previous result for analysis? y/n=1/0'
      write(31,*)'read previous result for analysis? y/n=1/0'
      read(5,*)iread
      write(6,*)'print coordinates? y=1/n=0'
      write(31,*)'print coordinates? y=1/n=0'
      read(5,*)ipr
      write(6,1004)
      write(31,1004)
      n=1
      do 10 j=1,jx
         if(iread.eq.1)then
            read(14,*)y(j),eta(j),z(j),zeta(j),s(j),ds(j),g(j)
     &        ,cm(j),dem(j),tm(j),alf(j),un(j)
         endif
         if(s(j).gt.sbreak(n)-eps)then
            write(6,1024)sbreak(n)
            write(31,1024)sbreak(n)
            n=n+1
         endif
         if(s(j).gt.sbreak(n)-eps)then
            write(6,1024)sbreak(n)
            write(31,1024)sbreak(n)
            n=n+1
         endif
         polar(j)=n
         if(ipr.eq.1)then
            write(6,1005)y(j),eta(j),z(j),zeta(j),s(j),ds(j),g(j)
     &           ,cm(j),dem(j),tm(j),alf(j),polar(j)
            write(31,1005)y(j),eta(j),z(j),zeta(j),s(j),ds(j),g(j)
     &           ,cm(j),dem(j),tm(j),alf(j),polar(j)
         endif
         reym=reynolds*cxm
         cdm0=1.328/sqrt(reym)
         if(reym.gt.1.0E5)then
            cdm0=.072/reym**.2
         endif
         cdm0=2.*cdm0
         a0(j)=cdm0
         a1(j)=0.
         a2(j)=2.*cdm0
         b0(j)=2.*pi*2.*dem(j)
         b1(j)=2.*pi
         c0(j)=-2.*pi*dem(j)
         c1(j)=-.5*pi
 10   continue
      write(6,*)
      write(6,*)'do you want to skip the optimization and run the'
     &     ,' analysis? y/n=1/0'
      write(31,*)
      write(31,*)'do you want to skip the optimization and run the'
     &     ,' analysis? y/n=1/0'
      read(5,*)iop
      iter=0
      if(iop.eq.1)goto 800
      ivis=0
      vis=0.
 200  continue
      if(ivis.eq.0)then
         write(6,*)'do you want to introduce viscous effects?'
         write(6,*)'viscous/inviscid=1/0   choose ivis='
         write(31,*)'do you want to introduce viscous effects?'
         write(31,*)'viscous/inviscid=1/0   choose ivis='
         read(5,*)ivis
         if(ivis.ne.0)then
            ivis=1
            vis=1.
         endif
      endif
      write(6,*)
      write(6,*)'#############################################'
      write(6,*)'optimization: do you want to iterate? y=1/n=0'
      write(31,*)
      write(31,*)'#############################################'
      write(31,*)'optimization: do you want to iterate? y=1/n=0'
      read(5,*)itx
      if(itx.le.0)goto 700
      if(ipolar.ne.1)then
         do 11 j=1,jx
            reym=reynolds*cxm
            cdm0=1.328/sqrt(reym)
            if(reym.gt.1.0E5)then
               cdm0=.072/reym**.2
            endif
            cdm0=2.*cdm0
            a0(j)=cdm0
            a1(j)=0.
            a2(j)=2.*cdm0
            alfj=(cltar-b0(j))/b1(j)
            cm0(j)=c0(j)+c1(j)*alfj
            l(j)=cltar
            if(j.lt.jl.or.j.gt.jr)then
               alfj=(cloptw-b0(j))/b1(j)
               l(j)=cloptw
               cm0(j)=0.
            endif
            alf(j)=alfj
            d(j)=vis*(a0(j)+a1(j)*l(j)+a2(j)*l(j)**2)
 11      continue
      else
c     search for point on polar and polar coefficients
      do 14 j=1,jx
         n=polar(j)
         czj=cltar
         if(j.lt.jl.or.j.gt.jr)then
            czj=cloptw
         endif
         mj=1
         prod=1.
         do 12 k=2,kx(n)-1
            prod=prod*(czj-cz(k,n))
            if(prod.le.0.)goto 13
            prod=1.
            mj=k
 12      continue
 13   continue
         b1(j)=(cz(mj+1,n)-cz(mj,n))/(inc(mj+1,n)-inc(mj,n))
         b0(j)=cz(mj,n)-b1(j)*inc(mj,n)
         alfj=inc(mj,n)+(inc(mj+1,n)-inc(mj,n))
     &        *(czj-cz(mj,n))/(cz(mj+1,n)-cz(mj,n))
         if(kxtrm(mj,n).ne.0)mj=mj+1
         if(mj.lt.2)mj=2
         if(mj.gt.kx(n)-1)mj=kx(n)-1
         m(j)=mj
         dcxdcz=(cx(mj,n)-cx(mj-1,n))/(cz(mj,n)-cz(mj-1,n))
         d2cxdcz2=((cx(mj+1,n)-cx(mj,n))*(cz(mj,n)-cz(mj-1,n))
     &        -(cx(mj,n)-cx(mj-1,n))*(cz(mj+1,n)-cz(mj,n)))
     &       /((cz(mj+1,n)-cz(mj,n))*(cz(mj+1,n)-cz(mj-1,n))
     &       *(cz(mj,n)-cz(mj-1,n)))
         cxj=cx(mj-1,n)+dcxdcz*(czj-cz(mj-1,n))
     &        +d2cxdcz2*(czj-cz(mj-1,n))*(czj-cz(mj,n))
         a0(j)=cx(mj-1,n)-dcxdcz*cz(mj-1,n)
     &        +d2cxdcz2*cz(mj-1,n)*cz(mj,n)
         a1(j)=dcxdcz-d2cxdcz2*(cz(mj-1,n)+cz(mj,n))
         a2(j)=d2cxdcz2
         alf(j)=alfj
         d(j)=vis*cxj
         l(j)=czj
         cm0(j)=cq(mj,n)+(cq(mj+1,n)-cq(mj,n))*(alfj-inc(mj,n))
     &        /(inc(mj+1,n)-inc(mj,n))-0.25*czj
         if(j.lt.jl.or.j.gt.jr)then
            cm0(j)=0.
         endif
 14   continue
      endif
c      write(6,*)'a0(j)=',(a0(j),j=1,jx)
c      write(6,*)'a1(j)=',(a1(j),j=1,jx)
c      write(6,*)'a2(j)=',(a2(j),j=1,jx)
 300  continue
      iter=0
      do 500 it=1,1000
      iter=iter+1
      cost=0.
      alagr=0.
      cl=0.
      do 16 j=2,jx-1
         A=0.
         B=0.
         jdx=0
         dgx=0.
         deta=eta(j)-eta(j-1)
         dzeta=zeta(j)-zeta(j-1)
         gj=g(j)
         km=1
         do 15 k=1,jx-1
            ck=(deta*(y(j)-eta(k))+dzeta*(z(j)-zeta(k)))/
     &           ((y(j)-eta(k))**2+(z(j)-zeta(k))**2)
            dk=((eta(k)-eta(km))*(y(k)-eta(j-1))+
     &           (zeta(k)-zeta(km))*(z(k)-zeta(j-1)))/
     &           ((y(k)-eta(j-1))**2+(z(k)-zeta(j-1))**2)-
     &           ((eta(k)-eta(km))*(y(k)-eta(j))+
     &           (zeta(k)-zeta(km))*(z(k)-zeta(j)))/
     &           ((y(k)-eta(j))**2+(z(k)-zeta(j))**2)
            if(k.eq.j-1)piv=ck
            if(k.eq.j)then
               piv=piv-ck
               piv=piv+dk
            endif
            A=A+(g(k+1)-g(k))*ck
            B=B+g(k)*dk
            cost=cost+gj*ck*(g(k+1)-g(k))
            km=k
 15      continue
         A=-A/(2.*pi)
         B=-B/(2.*pi)
         alagr=alagr+2.*g(j)*deta
         cost=cost+2.*pi*cm(j)*d(j)*ds(j)
         piv=piv/(2.*pi)+8.*vis*a2(j)*ds(j)/cm(j)
         dgj=(A+B-2.*al*deta
     &        -vis*(2.*a1(j)+8.*a2(j)*g(j)/cm(j))*ds(j)
     &        )/piv
         g(j)=g(j)+omegades*dgj
         if(abs(dgj).gt.abs(dgx))then
            dgx=dgj
            jdx=j
         endif
         cl=cl+g(j)*deta
 16   continue
      cost=cost/(2.*pi)+al*alagr
      cost=0.25*arm*cost
      fy(1)=0.
      fz(1)=0.
      fy(jx)=0.
      fz(jx)=0.
      cmf(1)=0.
      cmf(2)=0.
      cmf(jx)=0.
      do 17 j=2,jx-1
         if(j.eq.jx2+1)then
            fy(j)=fy(j-1)-(1.-si)*g(j)*(zeta(j)-zeta(j-1))
            fz(j)=-fz(j-1)+(1.-si)*g(j)*(eta(j)-eta(j-1))
         else
            fy(j)=fy(j-1)-g(j)*(zeta(j)-zeta(j-1))
            fz(j)=fz(j-1)+g(j)*(eta(j)-eta(j-1))
         endif
         cmf(j+1)=cmf(j)+0.5*fy(j)*(zeta(j+1)-zeta(j))
     &        -0.5*fz(j)*(eta(j+1)-eta(j))
 17   continue
      if(igeo.eq.'converged')goto 400
      dzx=0.
      jdz=0
      zm=0.
      dy=y(jx2+1+is)-y(j0)
      dzdy=0.
      deta=.5*y(jx2+1+is)
      do 18 j=j0+1,jx
         dz=zm+dy*(dzdy-deta*cmf(j-1)/eix)-z(j)+dih(j)
         if(abs(dz).gt.abs(dzx))then
            dzx=dz
            jdz=j
         endif
         z(j)=z(j)+dz
         z(jx+1-j)=z(j)
         zm=z(j)-dih(j)
         jp=j+1
         if(j.eq.jx)jp=jx
         dy=y(jp)-y(j)
         dzdy=(z(j)-dih(j)-z(j-1)+dih(j-1))/(s(j)-s(j-1))
         deta=eta(j)-eta(j-1)
         if(j.eq.jr)dztip=dz
 18   continue
      do 19 j=jl,jr-1
         zeta(j)=z(j)+(z(j+1)-z(j))*(eta(j)-y(j))/(y(j+1)-y(j))
 19   continue
      do 20 j=jr,jx-1
         zeta(j)=zeta(j)+dztip
         zeta(jx-j)=zeta(j)
 20   continue
      zeta(jx)=zeta(jx-1)
      if(abs(dzx).le.eps)igeo='converged'
 400  continue
      if(abs(dgx).lt.0.1*eps)goto 600
 500  continue
 600  continue
c     design of wing and winglets
      g0=g(j0)+eps
      aream=0.
      jm=1
      do 21 j=1,jx
         cm(j)=cxm*g(j)/g0+eps
         if(j.lt.jl.or.j.gt.jr)then
            cm(j)=2.*g(j)/cloptw+eps
c            cm(j)=cxm*g(j)/g0+eps
         endif
         if(j.lt.jx)aream=aream+cm(j)*(eta(j)-eta(jm))
         jm=j
 21   continue
      arm=4./aream
      cl=0.5*arm*cl
      write(6,*)' iter=',iter,' dgx=',dgx,' jdx=',jdx,' cost=',cost
      write(6,*)'  dzx=',dzx,' jdz=',jdz,' igeo=',igeo,' cl=',cl
      write(6,*)
      write(6,*)'************************************'
      write(6,*)'iterate Lagrange multiplier? y=1/n=0'
      write(31,*)' iter=',iter,' dgx=',dgx,' jdx=',jdx,' cost=',cost
      write(31,*)'  dzx=',dzx,' jdz=',jdz,' igeo=',igeo,' cl=',cl
      write(31,*)
      write(31,*)'************************************'
      write(31,*)'iterate Lagrange multiplier? y=1/n=0'
      read(5,*)ial
      if(ial.eq.1)then
         a1mean=a1(j0)
         ratio=0.
         if(abs(cl).gt.eps)ratio=cltar/cl
         al=ratio*al+vis*(ratio-1.)*a1mean
         write(6,*)'al=',al
         write(31,*)'al=',al
         igeo=' iterate '
         do 22 j=1,jx
            g(j)=ratio*g(j)
            un(j)=ratio*un(j)
            fy(j)=ratio*fy(j)
            fz(j)=ratio*fz(j)
            cmf(j)=ratio*cmf(j)
 22      continue
         goto 300
      else
         goto 200
      endif
 700  continue
      cl=0.
      do 24 j=2,jx-1
         deta=(eta(j)-eta(j-1))/ds(j)
         dzeta=(zeta(j)-zeta(j-1))/ds(j)
         sumj=0.
         do 23 k=1,jx-1
            ck=(deta*(y(j)-eta(k))+dzeta*(z(j)-zeta(k)))/
     &           ((y(j)-eta(k))**2+(z(j)-zeta(k))**2)
            sumj=sumj+(g(k+1)-g(k))*ck
 23      continue
         un(j)=-sumj/(4.*pi)
         cl=cl+g(j)*(eta(j)-eta(j-1))
 24   continue
      un(1)=un(2)+(un(3)-un(2))*(s(1)-s(2))/(s(3)-s(2))
      un(jx)=un(jx-1)+
     &     (un(jx-2)-un(jx-1))*(s(jx)-s(jx-1))/(s(jx-2)-s(jx-1))
      cl=0.5*arm*cl
c*****solution for target lift
      ratio=0.
      if(abs(cl).gt.eps)ratio=cltar/cl
      cdilw=0.
      cdirw=0.
      cdim=0.
      cdvlw=0.
      cdvrw=0.
      cdvm=0.
      aream=0.
      areaw=0.
      dztip=ratio*(z(jr)-dih(jr))
      jm=1
      d0=ratio**2*g(2)*un(2)/cm(2)
      if(jl.eq.1)then
         d0=ratio**2*g(j0)*un(j0)/cm(j0)
      endif
      do 25 j=1,jx
         g(j)=ratio*g(j)
         un(j)=ratio*un(j)
         fy(j)=ratio*fy(j)
         fz(j)=ratio*fz(j)
         cmf(j)=ratio*cmf(j)
         if(j.ge.jl.and.j.le.jr)then
            z(j)=dih(j)+ratio*(z(j)-dih(j))
            cdim=cdim-2.*g(j)*un(j)*ds(j)
            cdvm=cdvm+vis*d(j)*cm(j)*ds(j)
         else
            z(j)=dih(j)+dztip
            if(j.lt.jl)then
               cdilw=cdilw-2.*g(j)*un(j)*ds(j)
               cdvlw=cdvlw+vis*d(j)*cm(j)*ds(j)
               areaw=areaw+cm(j)*ds(j)
            else
               cdirw=cdirw-2.*g(j)*un(j)*ds(j)
               cdvrw=cdvrw+vis*d(j)*cm(j)*ds(j)
            endif
         endif
         dcdds(j)=-2.*g(j)*un(j)
         dcddse(j)=0.
         if(j.lt.jx)aream=aream+cm(j)*(eta(j)-eta(jm))
         d(j)=d(j)-(2.*g(j)*un(j)+2.0*eps*d0)/cm(j)
         jm=j
 25   continue
      if(am.ne.1.0)am=aream/(2.*cxm)
      arm=4./aream
      aw=1.
      if(cxw.gt.2.*eps.and.abs(at).gt.eps)then
         aw=areaw/(at*cxw)
      endif
      cdilw=0.25*arm*cdilw
      cdirw=0.25*arm*cdirw
      cdim=0.25*arm*cdim
      cdi=cdilw+cdirw+cdim
      dcdds(1)=0.
      dcdds(jx)=0.
      cdvlw=0.25*arm*cdvlw
      cdvrw=0.25*arm*cdvrw
      cdvm=0.25*arm*cdvm
      cdv=cdvlw+cdvrw+cdvm
      cd=cdi+cdv
      write(6,*)'print results: g,un,alf? y=1/n=0'
      write(31,*)'print results: g,un,alf? y=1/n=0'
      read(5,*)ipr
      if(ipr.eq.1)then
         write(6,*)'******scaled solution:'
         write(6,*)' g(j)=',(g(j),j=1,jx)
c         write(6,*)'un(j)=',(un(j),j=1,jx)
c         write(6,*)'alf(j)=',(alf(j),j=1,jx)
         write(6,*)'m(j)=',(m(j),j=1,jx)
         write(31,*)'******scaled solution:'
         write(31,*)' g(j)=',(g(j),j=1,jx)
c         write(31,*)'un(j)=',(un(j),j=1,jx)
c         write(31,*)'alf(j)=',(alf(j),j=1,jx)
         write(31,*)'m(j)=',(m(j),j=1,jx)
      endif
      do 26 j=1,jx
         write(16,*)y(j),dih(j)
         write(17,*)y(j),z(j)
         write(18,*)s(j),g(j),un(j)
 26   continue
      write(6,*)
      write(6,*)'******scaled solution:'
      write(6,1006)g(jl),g(jx2+1),g(jr)
      write(31,*)
      write(31,*)'******scaled solution:'
      write(31,1006)g(jl),g(jx2+1),g(jr)
      cl=ratio*cl
      cost=cd+al*cl
      write(6,1007)cl,cdi,cmf(jx2+1)
      write(6,1008)cost,al,cdvlw,cdvm,cdvrw,cd
      write(31,1007)cl,cdi,cmf(jx2+1)
      write(31,1008)cost,al,cdvlw,cdvm,cdvrw,cd
c*****exact solution for elliptic wing (no winglet)
      dt=pi/(jr-jl)
      gmax=4.*cl/(pi*arm)
      wi=-.25*gmax
      do 27 j=jl,jr
         tj=(j-jl)*dt
         ge(j)=gmax*sin(tj)
         we(j)=wi
         dcddse(j)=-2.*ge(j)*we(j)
         write(19,*)y(j),ge(j),we(j)
 27   continue
      write(6,*)
      write(6,*)'******exact solution for elliptic wing:'
      write(31,*)
      write(31,*)'******exact solution for elliptic wing:'
      cle=cl
      cdie=cl**2/(pi*arm)
      cdve=cdvm
      cde=cdie+cdve
      write(6,1009)cle,cdie,cdve,cde
      write(6,*)
      write(6,*)'******percent gain in drag:'
      write(31,1009)cle,cdie,cdve,cde
      write(31,*)
      write(31,*)'******percent gain in drag:'
      iper=0
      if(abs(cl).gt.eps)then
         iper=100-100.*cd/cde
      endif
      write(6,1010)iper
      write(31,1010)iper
      do 28 j=1,jx
         write(20,*)s(j),dcdds(j),dcddse(j)
         write(21,*)s(j),cmf(j)
         write(22,*)eta(j),fz(j)
 28   continue
c*****design and analysis
      write(6,*)
      write(6,*)'******design and analysis:'
      write(31,*)
      write(31,*)'******design and analysis:'
      g0=g(j0)+eps
      toein=0.
      if(jl.gt.1)then
         toein=alf(2)
      endif
      toeind=radeg*toein
      alfgeo=alf(j0)-atan2(un(j0),1.)
      alfgeod=radeg*alfgeo
      write(6,1011)cxm,arm,toein,toeind,alfgeo,alfgeod
      write(6,*)
      write(31,1011)cxm,arm,toein,toeind,alfgeo,alfgeod
      write(31,*)
      gkink=g(jl)
      cmom=0.
      cav=eps
      jm=1
      do 29 j=1,jx
         cm(j)=cxm*g(j)/g0
         dem(j)=dm
         tm(j)=0.
         if(j.lt.jl.or.j.gt.jr)then
            cm(j)=2.*g(j)/cloptw+eps
c            cm(j)=cxm*g(j)/g0+eps
            dem(j)=dw
            tm(j)=toein
         endif
         dum=-0.25*cm(j)-0.
         dumdum=0.75*cm(j)-0.
         write(23,*)s(j),dum,dumdum
c     ,dem(j),tm(j)
         if(j.lt.jx)then
            cmom=cmom+(cm0(j)*cm(j)**2
     &           -0.25*(cxm-cm(j))*l(j)*cm(j))*(eta(j)-eta(jm))
            cav=cav+cm(j)**2*(eta(j)-eta(jm))
         endif
         jm=j
 29   continue
      cav=0.25*arm*cav
      cmom=0.25*arm*cmom/cav
      iter=0
 800  continue
      if(iter.eq.0)then
         ivis=0
         vis=0.
      endif
      if(ivis.eq.0)then
         write(6,*)'do you want to introduce viscous effects?'
         write(6,*)'viscous/inviscid=1/0   choose ivis='
         write(31,*)'do you want to introduce viscous effects?'
         write(31,*)'viscous/inviscid=1/0   choose ivis='
         read(5,*)ivis
         if(ivis.ne.0)then
            ivis=1
            vis=1.
         endif
      endif
      jcxw=jl-1
      if(jcxw.lt.1)jcxw=1
      cxw=cm(jcxw)
      write(6,*)
      write(6,*)'analysis: itx=?'
      write(31,*)
      write(31,*)'analysis: itx=?'
      read(5,*)itx
      if(itx.le.0)goto 1000
      do 900 it=1,itx
      iter=iter+1
      if(ipolar.ne.1)then
         do 30 j=1,jx
            reym=reynolds*cxm
            cdm0=1.328/sqrt(reym)
            if(reym.gt.1.0E5)then
               cdm0=.072/reym**.2
            endif
            cdm0=2.*cdm0
            a0(j)=cdm0
            a1(j)=0.
            a2(j)=2.*cdm0
            b0(j)=2.*pi*2.*dem(j)
            b1(j)=2.*pi
            c0(j)=-2.*pi*dem(j)
            c1(j)=-.5*pi
            alfj=alf(j)
            czj=b0(j)+b1(j)*alfj
            d(j)=vis*(a0(j)+a1(j)*czj+a2(j)*czj**2)
            l(j)=czj
            cm0(j)=c0(j)+c1(j)*alfj
            if(j.le.jl.or.j.ge.jr)then
               cm0(j)=0.
            endif
 30      continue
      else
c     search for point on polar and polar coefficients
         do 33 j=1,jx
            n=polar(j)
            alfj=alf(j)
            mj=1
            prod=alfj+.5*pi
            do 31 k=2,kx(n)-1
               prod=prod*(alfj-inc(k,n))
               if(prod.le.eps)goto 32
               prod=1.
               mj=k
 31         continue
 32         continue
            b1(j)=(cz(mj+1,n)-cz(mj,n))/(inc(mj+1,n)-inc(mj,n))
            b0(j)=cz(mj,n)-b1(j)*inc(mj,n)
            czj=b0(j)+b1(j)*alfj
            if(kxtrm(mj,n).ne.0)mj=mj+1
            if(mj.lt.2)mj=2
            if(mj.gt.kx(n)-1)mj=kx(n)-1
            m(j)=mj
            dcxdcz=(cx(mj,n)-cx(mj-1,n))/(cz(mj,n)-cz(mj-1,n))
            d2cxdcz2=((cx(mj+1,n)-cx(mj,n))*(cz(mj,n)-cz(mj-1,n))
     &           -(cx(mj,n)-cx(mj-1,n))*(cz(mj+1,n)-cz(mj,n)))
     &           /((cz(mj+1,n)-cz(mj,n))*(cz(mj+1,n)-cz(mj-1,n))
     &           *(cz(mj,n)-cz(mj-1,n)))
            cxj=cx(mj-1,n)+dcxdcz*(czj-cz(mj-1,n))
     &           +d2cxdcz2*(czj-cz(mj-1,n))*(czj-cz(mj,n))
            a0(j)=cx(mj-1,n)-dcxdcz*cz(mj-1,n)+d2cxdcz2*cz(mj-1,n)
     &           *cz(mj,n)
            a1(j)=dcxdcz-d2cxdcz2*(cz(mj-1,n)+cz(mj,n))
            a2(j)=d2cxdcz2
            l(j)=czj
            d(j)=cxj
            cm0(j)=cq(mj,n)+(cq(mj+1,n)-cq(mj,n))*(czj-cz(mj,n))
     &           /(cz(mj+1,n)-cz(mj,n))
     &           -0.25*czj
            if(j.lt.jl.or.j.gt.jr)then
               cm0(j)=0.
            endif
 33      continue
      endif
      jdx=0
      dgx=0.
      cllw=0.
      clrw=0.
      cl=0.
      cdilw=0.
      cdirw=0.
      cdim=0.
      cdvlw=0.
      cdvrw=0.
      cdvm=0.
      cmom=0.
      cnom=0.
      cnow=0.
      clo=0.
      sum1=0.
      sum2=0.
      aream=0.
      areaw=0.
      cav=eps
      do 35 j=2,jx-1
         deta=(eta(j)-eta(j-1))/ds(j)
         dzeta=(zeta(j)-zeta(j-1))/ds(j)
         sumj=0.
         do 34 k=1,jx-1
            ck=(deta*(y(j)-eta(k))+dzeta*(z(j)-zeta(k)))/
     &           ((y(j)-eta(k))**2+(z(j)-zeta(k))**2)
            sumj=sumj+(g(k+1)-g(k))*ck
            if(k.eq.j-1)piv=ck
            if(k.eq.j)piv=piv-ck
 34      continue
         un(j)=-sumj/(4.*pi)
         aj=alpha
         if(j.lt.jl.or.j.gt.jr)then
            if(j.lt.jl)then
               aj=-yaw+deltal
            else
               aj=yaw-deltar
            endif
         endif
         aj=aj+atan2(20.*un(j),1.)/20.
         aj=aj+atan2(un(j),1.)
         alfj=aj+tm(j)
         czj=b0(j)+b1(j)*alfj
         visc=0.
         if(ipolar.eq.1)then
            if(m(j).ge.mxtrm(n))then
               visc=0.25*cm(j)*b1(j)*piv/(4.*pi*(1.+un(j)**2))
               if(visc.lt.eps)visc=0.
            endif
         endif
         dgj=(0.5*cm(j)*czj-g(j)
     &        +(avis+visc)*(g(j+1)-2.*g(j)+g(j-1)))
     &        /(1.+0.5*cm(j)*b1(j)*piv/(4.*pi*(1.+un(j)**2))
     &        +2.*(avis+visc))
         alf(j)=alfj
         g(j)=g(j)+omegana*dgj
         if(abs(dgx).lt.abs(dgj))then
            dgx=dgj
            jdx=j
         endif
         cl=cl+g(j)*(eta(j)-eta(j-1))
         if(j.lt.jl)then
            cllw=cllw-g(j)*(zeta(j)-zeta(j-1))
            cdilw=cdilw-2.*g(j)*un(j)*ds(j)
            cnow=cnow+2.*g(j)*un(j)*y(j)*ds(j)
            sum2=sum2+toein*g(j)*y(j)*ds(j)
            cdvlw=cdvlw+vis*d(j)*cm(j)*ds(j)
         else
            if(j.gt.jr)then
               clrw=clrw+g(j)*(zeta(j)-zeta(j-1))
               cdirw=cdirw-2.*g(j)*un(j)*ds(j)
               cnow=cnow+2.*g(j)*un(j)*y(j)*ds(j)
               sum2=sum2+toein*g(j)*y(j)*ds(j)
               cdvrw=cdvrw+vis*d(j)*cm(j)*ds(j)
            else
               cdim=cdim-2.*g(j)*un(j)*ds(j)
               cnom=cnom+2.*g(j)*un(j)*y(j)*ds(j)
               clo=clo+g(j)*y(j)*ds(j)
               sum2=sum2+alpha*g(j)*y(j)*ds(j)
               cdvm=cdvm+vis*d(j)*cm(j)*ds(j)
            endif
         endif
         cmom=cmom+(cm0(j)*cm(j)**2
     &        -0.25*(cxm-cm(j))*l(j)*cm(j))*(eta(j)-eta(j-1))
         if(j.lt.jl)then
            sum1=sum1-yaw*g(j)*y(j)*ds(j)
            areaw=areaw+cm(j)*ds(j)
         endif
         if(j.gt.jr)then
            sum1=sum1+yaw*g(j)*y(j)*ds(j)
         endif
         aream=aream+cm(j)*(eta(j)-eta(j-1))
         cav=cav+cm(j)**2*(eta(j)-eta(j-1))
 35   continue
      arm=4./aream
      if(am.ne.1.0)am=0.5*aream/cxm
      cav=0.25*arm*cav
      arw=at*at/(areaw+eps)
      aw=1.
      if(cxw.gt.2.*eps)then
         aw=areaw/((at+eps)*cxw)
      endif
      cl=0.5*arm*cl
      cllw=0.5*arm*cllw
      clrw=0.5*arm*clrw
      cdilw=0.25*arm*cdilw
      cdim=0.25*arm*cdim
      cdirw=0.25*arm*cdirw
      cdi=cdil+cdim+cdirw
      cmom=0.25*arm*cmom/cav
      cnom=0.5*arm*cnom
      cnow=0.5*arm*cnow
      cno=cnom+cnow
      clo=0.5*arm*clo
      cdvm=cdvm/aream
      cdvlw=cdvlw/aream
      cdvrw=cdvrw/aream
      un(1)=un(2)+(un(3)-un(2))*(s(1)-s(2))/(s(3)-s(2))
      un(jx)=un(jx-1)+
     &     (un(jx-2)-un(jx-1))*(s(jx)-s(jx-1))/(s(jx-2)-s(jx-1))
      l(1)=l(2)+(l(3)-l(2))*(s(1)-s(2))/(s(3)-s(2))
      l(jx)=l(jx-1)+
     &     (l(jx-2)-l(jx-1))*(s(jx)-s(jx-1))/(s(jx-2)-s(jx-1))
 900  continue
      write(6,*)' iter=',iter,' dgx=',dgx,' jdx=',jdx
      write(6,*)'m(j)=',(m(j),j=1,jx)
      write(31,*)' iter=',iter,' dgx=',dgx,' jdx=',jdx
      write(31,*)'m(j)=',(m(j),j=1,jx)
      goto 800
1000  continue
      dt=pi/(jr-jl)
      gmax=4.*cl/(pi*arm)
      wi=-.25*gmax
      write(6,*)
      write(6,1012)
      write(31,*)
      write(31,1012)
      do 36 j=1,jx
         write(30,*)d(j),l(j),alf(j)
         write(6,1013)y(j),z(j),cm(j),dem(j),tm(j),g(j),un(j)
     &        ,l(j),d(j),cm0(j),alf(j),polar(j)
         write(31,1013)y(j),z(j),cm(j),dem(j),tm(j),g(j),un(j)
     &        ,l(j),d(j),cm0(j),alf(j),polar(j)
         write(24,*)s(j),g(j),un(j)
         if(j.gt.jl.and.j.lt.jr)then
            tj=(j-jl)*dt
            ge(j)=gmax*sin(tj)
            we(j)=wi
            write(25,*)s(j),ge(j),we(j)
         endif
 36   continue
      cdie=cl*cl/(pi*arm)
      eme=1.
      cdve=cdvm
      clm=cl
      em=1.
      if(cdim.gt.eps)then
         em=clm*clm/(pi*arm*cdim)
      endif
      cde=cdie+cdve
      cdm=cdim+cdvm
      cdlw=cdilw+cdvlw
      cmolw=0.
      cdrw=cdirw+cdvrw
      cmorw=0.
      cdi=cdilw+cdim+cdirw
      cdv=cdvm+cdvlw+cdvrw
      e=1.
      if(cdi.gt.eps)then
         e=cl*cl/(pi*arm*cdi)
      endif
      cd=cdi+cdv
      write(6,*)
      write(6,*)'*******main wing geometry:'
      write(6,*)'    E*area moment inertia=',eix
      write(6,*)'           span height bt=',bt
      write(6,*)'           root chord cxm=',cxm
      write(6,*)' wing area coefficient am=',am
      write(6,*)'                wing area=',aream,' (ref. B**2/4)'
      write(6,*)'        wing aspect ratio=',arm
      write(6,*)'  main wing prof rel camb=',dm
      write(6,*)'          incidence alpha=',alpha,' (rad)  ='
     &     ,alphad,' (deg)'
      write(6,*)
      write(6,*)'*********winglet geometry:'
      write(6,*)'        winglet height at=',at
      write(6,*)'   winglet root chord cxw=',cxw
      write(6,*)'    winglet area coef. aw=',aw
      write(6,*)'             winglet area=',areaw,' (ref. B**2/4)'
      write(6,*)'     winglet aspect ratio=',arw
      write(6,*)'    winglet prof rel camb=',dw
      write(6,*)'                      yaw=',yaw,' (rad)  ='
     &     ,yawd,' (deg)'
      write(6,*)
      write(6,*)'***********************'
      write(6,*)' results elliptic wing:'
      write(6,1014)cdie,eme
      write(6,1015)cdve
      write(6,1016)cde,clm
      write(6,*)
      write(6,*)'***********************'
      write(6,*)'     results main wing:'
      write(6,1014)cdim,em
      write(6,1015)cdvm
      write(6,1016)cdm,clm
      write(6,1017)cmom
      write(6,*)'***********************'
      write(6,*)'  results left winglet:'
      write(6,1014)cdilw
      write(6,1015)cdvlw
      write(6,1016)cdlw,cllw
      write(6,1017)cmolw
      write(6,*)'***********************'
      write(6,*)' results right winglet:'
      write(6,1014)cdirw
      write(6,1015)cdvrw
      write(6,1016)cdrw,clrw
      write(6,1017)cmorw
      write(6,*)'***********************'
      write(6,*)'        global results:'
      write(6,1014)cdi,e
      write(6,1015)cdv
      write(6,1016)cd,cl
      write(6,1017)cmom
      write(6,*)'***********************'
      write(6,*)'   winglet deflections:'
      write(6,1018)deltal,deltald
      write(6,1019)deltar,deltard
      write(6,*)'         effect of yaw:'
      write(6,1020)yaw,yawd,cnom,cnow,cno,clo
      write(6,1021)cdilw,cdim,cdirw
      write(31,*)
      write(31,*)'*******main wing geometry:'
      write(31,*)'    E*area moment inertia=',eix
      write(31,*)'           span height bt=',bt
      write(31,*)'           root chord cxm=',cxm
      write(31,*)' wing area coefficient am=',am
      write(31,*)'                wing area=',aream,' (ref. B**2/4)'
      write(31,*)'        wing aspect ratio=',arm
      write(31,*)'  main wing prof rel camb=',dm
      write(31,*)'          incidence alpha=',alpha,' (rad)  ='
     &     ,alphad,' (deg)'
      write(31,*)
      write(31,*)'*********winglet geometry:'
      write(31,*)'        winglet height at=',at
      write(31,*)'   winglet root chord cxw=',cxw
      write(31,*)'    winglet area coef. aw=',aw
      write(31,*)'             winglet area=',areaw,' (ref. B**2/4)'
      write(31,*)'     winglet aspect ratio=',arw
      write(31,*)'    winglet prof rel camb=',dw
      write(31,*)'                      yaw=',yaw,' (rad)  ='
     &     ,yawd,' (deg)'
      write(31,*)
      write(31,*)'***********************'
      write(31,*)' results elliptic wing:'
      write(31,1014)cdie,eme
      write(31,1015)cdve
      write(31,1016)cde,clm
      write(31,*)
      write(31,*)'***********************'
      write(31,*)'     results main wing:'
      write(31,1014)cdim,em
      write(31,1015)cdvm
      write(31,1016)cdm,clm
      write(31,1017)cmom
      write(31,*)'***********************'
      write(31,*)'  results left winglet:'
      write(31,1014)cdilw
      write(31,1015)cdvlw
      write(31,1016)cdlw,cllw
      write(31,1017)cmolw
      write(31,*)'***********************'
      write(31,*)' results right winglet:'
      write(31,1014)cdirw
      write(31,1015)cdvrw
      write(31,1016)cdrw,clrw
      write(31,1017)cmorw
      write(31,*)'***********************'
      write(31,*)'        global results:'
      write(31,1014)cdi,e
      write(31,1015)cdv
      write(31,1016)cd,cl
      write(31,1017)cmom
      write(31,*)'***********************'
      write(31,*)'   winglet deflections:'
      write(31,1018)deltal,deltald
      write(31,1019)deltar,deltard
      write(31,*)'         effect of yaw:'
      write(31,1020)yaw,yawd,cnom,cnow,cno,clo
      write(31,1021)cdilw,cdim,cdirw
c*****save results
      write(6,1022)
      xlew=0.
      tmj=tm(jr+1)
      write(6,1023)z(jr),xlew,dem(jr),tmj
      write(26,1023)xlew,z(jr),dem(jr),tmj
      dt=pi/(jx-jr+eps)
      tj=0.
      crw=cm(jr+1)
      write(27,*)tj,cxm
      do 37 j=1,jx
         if(j.ge.jr+1)then
            tj=tj+dt
            xlew=1.-cm(j)/crw
            tmj=tm(j)
            zj=z(j)/at
            write(6,1023)zj,xlew,dem(j),tmj
            write(26,1023)xlew,zj,dem(j),tmj
            write(27,*)tj,cm(j)
         endif
         write(28,*)y(j),eta(j),z(j),zeta(j),s(j),ds(j),g(j)
     &        ,cm(j),dem(j),tm(j),alf(j),un(j)
 37   continue
      write(29,*)alphad,cl,cd,cmom
c*****files
      write(6,*)
      write(6,*)'******input files:'
      write(6,*)'polarbl.dat       =profile polars'
      write(6,*)'ould.in           =mesh and circulation for restart'
      write(6,*)'ould.data         =input data'
      write(6,*)'*****output files:'
      write(6,*)'polarbl.cdcl      =cd,cl,inc'
      write(6,*)'ould.yz0          =unloaded given dihedral shape'
      write(6,*)'ould.yz           =loaded dihedral shape'
      write(6,*)'ould.gun          =g(j) & un(j) from optimization'
      write(6,*)'ould.gwe          =ge(j) & we(j) for elliptic loading'
      write(6,*)'ould.cdy          =induced drag distributions'
      write(6,*)'ould.ymf          =flexion moment distribution'
      write(6,*)'ould.yfz          =normal force distribution'
      write(6,*)'ould.cdt          =y(j),cm(j),dem(j),tm(j)'
      write(6,*)'ould.guna         =g(j),un(j) from analysis'
      write(6,*)'ould.gwea         =ge(j),we(j) analysis'
      write(6,*)'ould.zxdtw        =z(j),xle(j),dem(j),tm(j) on winglet'
      write(6,*)'ould.tcw          =theta(j),cm(j) on winglet'
      write(6,*)'ould.out          =y(j),eta(j),z(j),zeta(j),g(j),etc'
      write(6,*)'ould.clcdcq       =alphad,cl,cd,cmom'
      write(6,*)'ould.cdcl         =cd(j),cl(j),alf(j)'
      write(6,*)'ould.list         =listing'
 1001 format(18a4)
 1002 format(3x,'k=',5x,'inc(k)=',6x,'cz(k)=',6x,'cx(k)='
     &     ,6x,'cq(k)=')
 1003 format(1x,i4,4x,f8.4,4x,f8.4,4x,f8.4,4x,f8.4)
 1004 format(4x,'y(j)',4x,'eta(j)',4x,'z(j)',4x,'zeta(j)',3x,'s(j)',
     &     5x,'ds(j)',4x,'g(j)',5x,'cm(j)',3x,'dem(j)',4x,'tm(j)',
     &     3x,'alf(j)',4x,'polar(j)')
 1005 format(1x,f8.4,1x,f8.4,1x,f8.4,1x,f8.4,1x,f8.4,
     &     1x,f8.4,1x,f8.4,1x,f8.4,1x,f8.4,1x,f8.4,1x,f8.4,1x,i5)
 1006 format(' g(jl)=',f10.6,'  g(jx/2)=',f10.6,'      g(jr)=',f10.6)
 1007 format('    CL=',f10.6,'      CDi=',f10.6,'  CMf(jx/2)=',f10.6)
 1008 format('  Cost=',f10.6,'       al=',f10.6,
     &     /,' CDvlw=',f10.6,'     CDvm=',f10.6,'      CDvrw=',f10.6
     &     /,'    CD=',f10.6)
 1009 format('   CLe=',f10.6,'     CDie=',f10.6,'       CDve=',f10.6
     &     /,'   CDe=',f10.6)
 1010 format('  gain=',i4,'%')
 1011 format('     cxm=',f10.6,'      ARm=',f10.6,/
     &      ,'toeinout=',f10.6,' (rd)    =',f10.6,'  (deg)',/
     &      ,'  alfgeo=',f10.6,' (rd)    =',f10.6,'  (deg)')
 1012 format(4x,'y(j)',5x,'z(j)',5x,'cm(j)',3x,'dem(j)',4x,'tm(j)'
     &     ,4x,'g(j)',5x,'un(j)',4x,'cl(j)',4x,'cd(j)',3x,'cm0(j)'
     &     ,3x,'alf(j)',4x,'polar(j)')
 1013 format(1x,f8.4,1x,f8.4,1x,f8.4,1x,f8.4,1x,f8.4,1x,f8.4
     &     ,1x,f8.4,1x,f8.4,1x,f8.4,1x,f8.4,1x,f8.4,1x,i5)
 1014 format(' inviscid contributions: CDi=',f10.6
     &     ,'      wing efficiency e=',f10.4)
 1015 format('  viscous contributions: CDv=',f10.6)
 1016 format('         global results:  CD=',f10.6
     &     ,'    lift coefficient CL=',f10.4)
 1017 format('                                      '
     &     ,' moment coefficient CM,o=',f10.4)
 1018 format('deltal=',f10.6,' (rd) =',f10.6)
 1019 format('deltar=',f10.6,' (rd) =',f10.6)
 1020 format('   yaw=',f10.6,' (rd) =',f10.6,' (deg)   CN,om=',f10.6
     &     ,'  CN,ow=',f10.6,'  CN,o=',f10.6,'  CL,o=',f10.6)
 1021 format(' CDilw=',f10.6,'  CDim=',f10.6,'         CDirw=',f10.6)
 1022 format(2x,'z(j)/at',8x,'xle(j)',9x,'dem(j)',6x,'toe-in(j)')
 1023 format(f10.6,4x,f10.6,4x,f10.6,4x,f10.6)
 1024 format(4x,f10.4,'-----------> break - new polar <------------')
      end

      real function dihedr(yj)
      implicit none
      integer n
      real eps,bt,teta,yj
      common/para/eps,bt,n
      teta=asin((yj)**(n/2.))
      dihedr=bt*(1.+4.*eps/n-(cos(teta)+eps)**(2./n))
      return
      end
      


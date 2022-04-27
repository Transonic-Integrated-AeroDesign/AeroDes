      program prandtline
      implicit none
      integer jxx,lxx,nxx,ipolar,nx,n,km,kfirst,k,kdum,ice,kp,jx
      integer jx2,is,iwing,j,nsteps,ivis,nstep,iter,it,mj,jdx,jm
      integer itx
      parameter(jxx=201,lxx=101,nxx=10)
      real eps,pi,degrad,prod,dcz,dcxm,dcxp,incd,si,omega,avis
      real B,cxm,dm,tmd,Rho,Vinf,Amu,alphad,tm,alpha,Re,amdum
      real camdum,dtet,tetj,yj,etaj,etajm,am,cam,arm,alphain
      real alphafi,alstep,vis,cxj,czj,qj,dgx,sum,wj,atj,attj
      real reg,res0,alogres,cl,cm0,xac,cmac,cd0,sum0
      real rey,cdi,cdv,em,cd,dum,acwash,xcp,Cx0,Rstr0,rstr
      real Rf0,rf,phij,phi0,dwkj,Lambd,lamb,clf
      real c(jxx),g(jxx),dg(jxx),y(jxx),eta(jxx)
      real w(jxx),t(jxx),dem(jxx)
      real a0(jxx),a1(jxx),b0(jxx),b1(jxx),c0(jxx),c1(jxx)
      real l(jxx),d(jxx),q(jxx),at(jxx)
      real cx(lxx,nxx),cz(lxx,nxx),cq(lxx,nxx),inc(lxx,nxx)
      real cmf(jxx),cmt(jxx),fz(jxx)
      real xle(jxx),xte(jxx),wcanar(jxx),xacm(jxx),xiac(jxx)
      real rbreak(nxx)
      integer m(jxx),polar(jxx)
      integer kx(nxx),kxtrm(lxx,nxx),mxtrm(nxx)
      character*4 bry(18)
      character*4 title(18)
      character*4 typcode(18)
      data a0/jxx*0./a1/jxx*0./b0/jxx*0./b1/jxx*0./c0/jxx*0./c1/jxx*0./
      data wcanar/jxx*0./
      data m/jxx*0/rbreak/nxx*2./
      open(unit=13,file='polarbl.dat',form='formatted')
      open(unit=14,file='polarbl.cdcl',form='formatted')
      open(unit=15,file='prandtline.data',form='formatted')
      open(unit=16,file='prandtline.in',form='formatted')
      open(unit=17,file='prandtline.ycdt',form='formatted')
      open(unit=18,file='prandtline.ygw',form='formatted')
      open(unit=19,file='prandtline.yg',form='formatted')
      open(unit=20,file='prandtline.yw',form='formatted')
      open(unit=21,file='prandtline.ycl',form='formatted')
      open(unit=22,file='prandtline.yat',form='formatted')
      open(unit=23,file='prandtline.ycd',form='formatted')
      open(unit=24,file='prandtline.out',form='formatted')
      open(unit=25,file='prandtline.clcdcq',form='formatted')
      open(unit=26,file='prandtline.ycmf',form='formatted')
      open(unit=27,file='prandtline.yfz',form='formatted')
      open(unit=28,file='prandtline.ycmt',form='formatted')
      open(unit=29,file='prandtline.yxlexte',form='formatted')
      open(unit=30,file='prandtline.incld',form='formatted')
      open(unit=31,file='prandtline.itres',form='formatted')
      open(unit=32,file='canarwash.ylwl',form='formatted')
      open(unit=33,file='prandtline.etaxiac',form='formatted')
      open(unit=34,file='prandtline.yxfuse',form='formatted')
      open(unit=35,file='prandtline.cdcl',form='formatted')
c*****constants
      eps=1.e-7
      pi=2.*asin(1.)
      degrad=pi/180.
c*****polar data
      write(6,*)
      write(6,*)'******do you want to use polar data? Y/N=1/0'
      read(5,*)ipolar
      if(ipolar.eq.0)goto 5
      read(13,*)nx
      write(6,*)'******profile polars:'
      write(6,*)' nx= ',nx,' number of polars to be read'
      if(nx.ge.nxx)then
         write(6,*)'!! nx > nxx !!'
         write(6,*)'TOO MANY POLARS: EXITING!'
         stop
      endif
      do 4 n=1,nx     
      write(6,*)'******n= ',n
      read(13,1000)title
      write(6,1000)title
      read(13,1000)title
      write(6,1000)title
      read(13,1000)title
      write(6,1000)title
      read(13,1000)title
      write(6,1000)title
      read(13,1000)title
      write(6,1000)title
      read(13,1000)title
      write(6,1000)title
      read(13,1000)title
      write(6,1000)title
      read(13,1000)title
      write(6,1000)title
      read(13,1000)title
      write(6,1000)title
      read(13,1000)title
      write(6,1000)title
      read(13,1000)title
      write(6,1000)title
      read(13,1000)title
      write(6,1000)title
      write(6,*)
      read(5,1000)bry
      write(6,*)
      write(6,*)'*****extrema of the cl(alpha) function:'
      prod=1.
      km=1
      kfirst=0
      do 1 k=1,lxx
         kdum=k
         read(13,*,end=2)inc(k,n),cz(k,n),cx(k,n),dum,cq(k,n)
         if(inc(k,n).gt.89.)then
            read(13,*)rbreak(n)
            write(6,*)' n=',n,' rbreak(n)=',rbreak(n)
            if(nx.eq.1)then
               read(13,*)rbreak(2)
            endif
            kdum=k+1
            goto 2
         endif
         kxtrm(k,n)=0
         if(k.eq.1)goto 1
         dcz=cz(k,n)-cz(km,n)
         prod=prod*dcz
         if(prod.lt.-eps)then
            write(6,*)'kxtrm(',n,')=',km,' cz(kxtrm,n)=',cz(km,n)
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
      if(kx(n).eq.lxx-1)then
         write(6,*)' attention: check if all data has been read;'
     &        ,' continuing/exiting=1/0?'
         read(5,*)ice
         if(ice.eq.0)stop
      endif
      write(6,*)
      write(6,*)'*************profile data from Xfoil:'
      write(6,1001)
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
         endif
         incd=inc(k,n)
         inc(k,n)=degrad*inc(k,n)
         write(6,*)' k=',k,' inc(k,n)=',inc(k,n),' cz(k,n)=',cz(k,n)
     &        ,' cx(k,n)=',cx(k,n),' cq(k,n)=',cq(k,n)
         write(14,*)cx(k,n),cz(k,n),cq(k,n),incd
 3    continue
      write(6,*)
      write(6,*)'******extrema pointer:'
      write(6,*)'kxtrm(k,',n,')=',(kxtrm(k,n),k=1,kx(n))
      write(6,*)
 4    continue
 5    continue
c*****dimensionless variables
      write(6,*)'************************'
      write(6,*)'dimensionless variables:'
      write(6,*)'                      Y=0.5*B*y'
      write(6,*)'                      C=0.5*B*c'
      write(6,*)'                      A=0.25*B**2*am'
      write(6,*)'                      D=C*d'
      write(6,*)'                      W=U*w'
      write(6,*)'                   GAMA=0.5*U*B*g'
      write(6,*)'                   LIFT=0.5*RHO*U**2*A*Cl'
      write(6,*)'                   DRAG=0.5*RHO*U**2*A*Cd'
      write(6,*)'               MOMENT,0=0.5*RHO*U**2*A*Cam*Cm0'
      write(6,*)'               REYNOLDS=RHO*U*C/AMU'
      write(6,*)'                     Fz=0.5*RHO*U**2*A*fz'
      write(6,*)'                     Mf=0.25*RHO*U**2*A*B*cmf'
      write(6,*)'                     Mt=0.5*RHO*U**2*A*Cam*cmt'
c*****read in data
      read(15,*)jx
      if(jx.gt.jxx)then
         write(6,*)'jx=',jx,'>jxx=',jxx,', change dimension:exiting!'
         stop
      endif
      read(15,*)itx
      read(15,*)omega
      read(15,*)avis
      read(15,*)B
      read(15,*)Cx0
      read(15,*)Lambd
      read(15,*)Rstr0
      read(15,*)Rf0
      read(15,*)dm
      read(15,*)tmd
      read(15,*)iwing
      read(15,*)alphad
      read(15,*)acwash
      read(15,*)Rho
      read(15,*)Vinf
      read(15,*)Amu
      cxm=2.0*Cx0/B
      lamb=degrad*Lambd
      rstr=2.0*Rstr0/B
      rf=2.0*Rf0/B
      tm=degrad*tmd
      alpha=degrad*alphad
c*****initialization
c     downwash due to canard on wing for Clc/(pi*arc)=0.1
      if(acwash.ne.0.0)then
	read(32,*)(dum,wcanar(j),j=1,jx)
      endif
c*****mesh, geometry and flow
      write(6,*)' '
      write(6,*)'******point distribution, geometry and flow:'
      write(6,1002)
      amdum=0.
      camdum=0.
      dtet=pi/(jx-1)
      n=1
      do 6 j=1,jx
         tetj=(j-1)*dtet
         yj=-cos(tetj)
         y(j)=yj
         etaj=-cos(tetj+.5*dtet)
         eta(j)=etaj
         if(j.eq.1)then
            etajm=-1.
         else
            etajm=eta(j-1)
         endif
         if(j.eq.jx-1)then
            etaj=1.0
         endif
         dem(j)=dm
         t(j)=tm
         g(j)=0.
         w(j)=0.
         at(j)=0.
         a0(j)=-pi*dm
         a1(j)=0.
         b0(j)=2.*pi*(2.*dem(j))
         b1(j)=2.*pi
c*****elliptic wing
         if(iwing.eq.0)then
            c(j)=cxm*sin(tetj)
            xacm(j)=0.25*cxm
            xacm(j)=xacm(j)+tan(lamb)*abs(y(j))
            xle(j)=xacm(j)-0.25*c(j)
            xte(j)=xle(j)+c(j)
            if(j.ge.2)xiac(j-1)=0.5*(xacm(j)+xacm(j-1))
            amdum=amdum+c(j)*(etaj-etajm)
            camdum=camdum+c(j)**2*(etaj-etajm)
         endif
c*****rectangular wing
         if(iwing.eq.1)then
            c(j)=cxm
            xle(j)=0.0
            xte(j)=cxm
            xacm(j)=0.25*cxm
            xiac(j)=0.25*cxm
            amdum=amdum+c(j)*(etaj-etajm)
            camdum=camdum+c(j)**2*(etaj-etajm)
         endif
c******tailess configuration
         if(iwing.eq.2)then
            if(abs(yj).ge.rf)then
               if(abs(yj).ge.rstr)then
                  xle(j)=cxm+0.01589-0.468*(1.0-abs(yj))
                  xte(j)=cxm+0.1269-0.181*(1.0-abs(yj))
               else
                  xle(j)=cxm-0.31111-1.0*(0.3-abs(yj))
                  xte(j)=cxm
               endif
               c(j)=xte(j)-xle(j)
               xacm(j)=xle(j)+0.25*c(j)
               if(j.ge.2)then
                  xiac(j-1)=0.5*(xacm(j)+xacm(j-1))
               endif
               amdum=amdum+c(j)*(etaj-etajm)
               camdum=camdum+c(j)**2*(etaj-etajm)
	       a0(j)=0.
            else
c     slender body treatment of fuselage
               phij=acos(yj/rf)
               xle(j)=rf*(1.0-sin(phij))
               xte(j)=cxm
               c(j)=xte(j)-xle(j)
               xacm(j)=xacm(j-1)+0.033333*((yj/0.11111)**2
     &              -(y(j-1)/0.11111)**2)
               if(j.ge.2)then
                  xiac(j-1)=0.5*(xacm(j)+xacm(j-1))
               endif
               amdum=amdum+c(j)*(etaj-etajm)
               camdum=camdum+c(j)**2*(etaj-etajm)
               a0(j)=0.
            endif
         endif
         if(ipolar.ne.0)then
            if(nx.gt.1)then
               polar(j)=n
               if(y(j).gt.rbreak(n)-eps)then
                  write(6,*)'rbreak(n)=',rbreak(n)
                  n=n+1
                  polar(j)=n
               endif
            else
               polar(j)=n
               if(y(j).gt.rbreak(n)-eps)then
                  polar(j)=0
               endif
               if(y(j).gt.rbreak(2)-eps)then
                  polar(j)=1
               endif
            endif
         else
            polar(j)=0
         endif
c         write(6,1003)y(j),eta(j),c(j),t(j),dem(j),g(j),w(j),at(j)
c     &        ,polar(j)
         write(17,*)y(j),c(j),dem(j),t(j)
         write(29,*)y(j),xle(j),xte(j),xacm(j)
         if(j.ge.2)write(33,*)eta(j-1),xiac(j-1)
 6    continue
      write(6,*)
      eta(jx)=eta(jx-1)
      xiac(jx)=xiac(jx-1)
      am=amdum
      cam=camdum/am
      arm=4./am
      Re=0.5*Rho*Vinf*B/Amu
      write(6,*)'***************'
      write(6,*)'numerical data:'
      write(6,*)'          number of points jx=',jx
      write(6,*)' max number of iterations itx=',itx
      write(6,*)'                        omega=',omega
      write(6,*)'   viscosity coefficient avis=',avis
      write(6,*)'main wing data:'
      write(6,*)'                  wing span B=',B,' (m)'
      write(6,*)'   maximum chord/fuselage Cx0=',Cx0,' (m)'
      write(6,*)'      a. c. sweep angle Lambd=',Lambd,'(deg)'
      write(6,*)'    half span of strake Rstr0=',Rstr0,' (m)'
      write(6,*)'          fuselage radius Rf0=',Rf0,' (m)'
      write(6,*)'    relative camber height dm=',dm,' (ref. C)'
      write(6,*)'        wing setting angle tm=',tm,' (rd) ='
     &     ,tmd,' (deg)'
      write(6,*)'             wing shape 0/1/2=',iwing
      write(6,*)'    downwash of canard acwash=',acwash
      write(6,*)'air data:'
      write(6,*)'              air density Rho=',Rho,' (kg/m**3)'
      write(6,*)'           wind velocity Vinf=',Vinf,' (m/s)'
      write(6,*)'        dynamic viscosity Amu=',Amu,' (kg/(m*s))'
      write(6,*)'        reference Reynolds Re=',Re,' (ref. Cam)'
      write(6,*)'****************'
      write(6,*)'calculated data:'
      write(6,*)'   maximum chord/fuselage cxm=',cxm,' ref. B/2)'
      write(6,*)'   wing+fuse planform area am=',am,' (ref. B**2/4)'
      write(6,*)'   wing+fuse aspect ratio arm=',arm
      write(6,*)'average aerodynamic chord cam=',cam,' (ref. B/2)'
      write(6,*)'           fuselage radius rf=',rf,'(ref. B/2)'
      write(6,*)' '
c*****iterations
      write(6,*)' alphain alphafi alstep ?'
      read(5,*)alphain,alphafi,alstep
      if(alstep.eq.0)then
         nsteps=1
      else
         nsteps=1+(alphafi-alphain)/alstep
      endif
      alphad=alphain-alstep
      write(6,*)'alphain=',alphain,'alphafi=',alphafi,'alstep=',alstep
      write(6,*)'nsteps=',nsteps
      if(nsteps.le.0.or.nsteps.gt.101)write(6,*)'bad sequence, exit'
      if(nsteps.le.0.or.nsteps.gt.101)stop
      write(6,*)'name of saved polar?'
      read(5,1000)title
      write(25,1000)title
      write(30,1000)title
      ivis=0
      vis=0.
      do 400 nstep=1,nsteps
      alphad=alphad+alstep
      alpha=degrad*alphad
      if(abs(alphad).ge.91.0)then
         write(6,*)'alphad>91 deg, stop'
         stop
      endif
      iter=0
      write(6,*)
      write(6,*)'****************solution:'
      if(ivis.eq.0)then
         write(6,*)'do you want to introduce viscous effects?'
         write(6,*)'viscous/inviscid=1/0   choose ivis='
         read(5,*)ivis
         if(ivis.ne.0)then
            ivis=1
            vis=1.
         endif
      else
         ivis=1
         vis=1.
      endif
      do 200 it=1,itx
      iter=iter+1
      if(ipolar.ne.1)goto 10
c     search for point on polar and polar coefficients
      do 9 j=2,jx-1
         n=polar(j)
         if(n.ne.0)then
            atj=at(j)
            atj=atj+acwash*wcanar(j)
            mj=2
            prod=1.57-atj
            do 7 k=2,kx(n)-1
               prod=prod*(inc(k,n)-atj)
               if(prod.ge.-eps)goto 8
               prod=1.
               mj=k
 7          continue
 8          continue
            a1(j)=(cx(mj+1,n)-cx(mj,n))/(inc(mj+1,n)-inc(mj,n))
            a0(j)=cx(mj,n)-a1(j)*inc(mj,n)
            b1(j)=(cz(mj+1,n)-cz(mj,n))/(inc(mj+1,n)-inc(mj,n))
            b0(j)=cz(mj,n)-b1(j)*inc(mj,n)
            c1(j)=(cq(mj+1,n)-cq(mj,n))/(inc(mj+1,n)-inc(mj,n))
            c0(j)=cq(mj,n)-c1(j)*inc(mj,n)
            cxj=a0(j)+a1(j)*atj
            czj=b0(j)+b1(j)*atj
            qj=c0(j)+c1(j)*atj
            m(j)=mj
            l(j)=czj
            d(j)=vis*cxj
            q(j)=qj
         else
            l(j)=2.0*g(j)/c(j)
            d(j)=0.0
            q(j)=0.0
            c1(j)=0.0
            c0(j)=-pi*dm
         endif
 9    continue
      m(1)=m(2)
      l(1)=l(2)+(l(3)-l(2))*(y(1)-y(2))/(y(3)-y(2))
      d(1)=d(2)+(d(3)-d(2))*(y(1)-y(2))/(y(3)-y(2))
      q(1)=q(2)+(q(3)-q(2))*(y(1)-y(2))/(y(3)-y(2))
      m(jx)=m(jx-1)
      l(jx)=l(jx-1)+(l(jx-1)-l(jx-2))*(y(jx)-y(jx-1))
     &     /(y(jx-1)-y(jx-2))
      d(jx)=d(jx-1)+(d(jx-1)-d(jx-2))*(y(jx)-y(jx-1))
     &     /(y(jx-1)-y(jx-2))
      q(jx)=q(jx-1)+(q(jx-1)-q(jx-2))*(y(jx)-y(jx-1))
     &     /(y(jx-1)-y(jx-2))
 10   continue
c     fixed point iteration
      g(1)=0.
      g(jx)=0.
      dgx=0.
      jdx=0
      do 12 j=2,jx-1
         sum=0.
         do 11 k=1,jx-1
c     downwash for non-straight lifting line
c     trailed vorticity
            phi0=sign(1.0,y(j)-eps)
     &           *atan((xacm(j)-xiac(k))/(y(j)-eta(k)))
            sum=sum+(g(k+1)-g(k))*(1.0-sin(phi0))/(y(j)-eta(k))
c     downwash due to lifting line
            if(k+1.ne.j.and.k.lt.jx-1)then
               dwkj=-g(k+1)*((xiac(k+1)-xiac(k))*(y(j)-y(k+1))
     &              -(xacm(j)-xacm(k+1))*(eta(k+1)-eta(k)))
     &              /((xacm(j)-xacm(k+1))**2+(y(j)-y(k+1))**2
     &              +10.0*(eta(k+1)-eta(k)))**1.5
               sum=sum+dwkj
            endif
 11      continue
         wj=-sum/(4.*pi)
         atj=alpha+atan2(wj,1.)
         atj=atj+acwash*wcanar(j)
         attj=atj+t(j)
         czj=b0(j)+b1(j)*attj
         if(ipolar.eq.1)then
            reg=0.
            if(m(j).ge.mxtrm(n))then
               reg=-c(j)*b1(j)
     &              *(1./(y(j)-eta(j-1))-1./(y(j)-eta(j)))
     &              /(16.*pi*(1.0+wj*wj))
               if(reg.lt.eps)reg=0.
            endif
         endif
         dg(j)=(0.5*c(j)*czj-g(j)
     &        +(avis+reg)*(g(j+1)-2.*g(j)+g(j-1)))
     &        /(1.+c(j)*b1(j)*(1./(y(j)-eta(j-1))-1./(y(j)-eta(j)))
     &        /(8.*pi*(1.0+wj*wj))+2.*(avis+reg))
         w(j)=wj
         at(j)=attj
         if(iwing.eq.2.and.abs(y(j)).le.rf)then
           dg(j)=g(j-1)+2.0*rf*sin(3.0*alpha)/3.0*((1.0-(y(j)/rf)**2)
     &           /(1.0+(y(j)/rf)**2)-(1.0-(y(j-1)/rf)**2)
     &           /(1.0+(y(j-1)/rf)**2))-g(j)
         endif
         g(j)=g(j)+omega*dg(j)
         if(abs(dgx).lt.abs(dg(j)))then
            dgx=dg(j)
            jdx=j
         endif
 12   continue
      if(iter.eq.1)res0=dgx
      alogres=alog10(abs(dgx/res0)+eps)
      write(31,*)iter,alogres
      w(1)=w(2)
     &     +(w(3)-w(2))*(y(1)-y(2))/(y(3)-y(2))
      w(jx)=w(jx-1)
     &     +(w(jx-2)-w(jx-1))*(y(jx)-y(jx-1))/(y(jx-2)-y(jx-1))
      at(1)=at(2)
     &     +(at(3)-at(2))*(y(1)-y(2))/(y(3)-y(2))
      at(jx)=at(jx-1)
     &     +(at(jx-2)-at(jx-1))*(y(jx)-y(jx-1))/(y(jx-2)-y(jx-1))
      if(abs(dgx).lt.eps)goto 300
 200  continue
      write(6,*)'************NOT CONVERGED!!!!!'
 300  continue
c*****results
      write(6,*)'alphad=',alphad,' deg'
      write(6,*)'iter=',iter,' dgx=',dgx,' jdx=',jdx
      write(6,*)'m(j)=',(m(j),j=1,jx)
      jx2=jx/2
      is=mod(jx,2)
      si=is
      cl=0.
      clf=0.
      xac=0.
      cmac=0.
      fz(1)=0.
      fz(jx)=0.
      cmf(1)=0.
      cmf(2)=0.
      cmf(jx)=0.
      cmt(1)=0.
      cmt(jx)=0.
      do 13 j=2,jx-1
         if(j.eq.2)eta(j-1)=-1.0
         if(j.eq.jx-1)eta(j)=1.0
         cl=cl+g(j)*(eta(j)-eta(j-1))
         if(abs(y(j)).lt.rf)then
            clf=clf+g(j)*(eta(j)-eta(j-1))
         endif
         if(ipolar.ne.1)then
            at(j)=alpha+atan2(w(j),1.)+t(j)
            q(j)=c0(j)+c1(j)*at(j)
            xac=xac+xacm(j)*c(j)*(eta(j)-eta(j-1))
            cmac=cmac+c(j)**2*q(j)*(eta(j)-eta(j-1))
            if(j.eq.jx2+1)then
               cmt(j)=-cmt(j-1)+(1.-si)*c(j)**2*q(j)*(eta(j)-eta(j-1))
            else
               cmt(j)=cmt(j-1)+c(j)**2*q(j)*(eta(j)-eta(j-1))
            endif
         else
            xac=xac+xacm(j)*c(j)*(eta(j)-eta(j-1))
            cmac=cmac+c(j)**2*q(j)*(eta(j)-eta(j-1))
         endif
         if(j.eq.jx2+1)then
            fz(j)=-fz(j-1)+(1.-si)*g(j)*(eta(j)-eta(j-1))
            cmt(j)=-cmt(j-1)+(1.-si)*c(j)**2
     &           *(q(j)+(xacm(j+1)-xacm(j))*fz(j)/cam)
     &           *(eta(j)-eta(j-1))
         else
            fz(j)=fz(j-1)+g(j)*(eta(j)-eta(j-1))
            cmt(j)=cmt(j-1)+c(j)**2
     &           *(q(j)+(xacm(j+1)-xacm(j))*fz(j)/cam)
     &           *(eta(j)-eta(j-1))
         endif
         cmf(j+1)=cmf(j)-fz(j)*(eta(j+1)-eta(j))
 13   continue
      at(1)=at(2)
     &     +(at(3)-at(2))*(y(1)-y(2))/(y(3)-y(2))
      at(jx)=at(jx-1)
     &     +(at(jx-2)-at(jx-1))*(y(jx)-y(jx-1))/(y(jx-2)-y(jx-1))
      cl=0.5*arm*cl
      clf=2.0*clf/am
      xac=xac/am
      cmac=cmac/(am*cam)
      cm0=cmac-xac*cl/cam
      xcp=xac-cam*cmac/cl
      cd0=0.
      sum=0.
      sum0=0.
      do 14 j=1,jx
         jm=j-1
         if(j.eq.1)jm=1
         czj=b0(j)+b1(j)*at(j)
         sum=sum+g(j)*w(j)*(eta(j)-eta(jm))
         if(ipolar.ne.1)then
            rey=c(j)*Re
            if(rey.lt.1.0e5)rey=1.0e5
            cd0=1.328/sqrt(rey)
            if(rey.gt.1.0e5)cd0=.072/rey**.2
            cd0=2.*cd0
            sum0=sum0+cd0*(eta(j)-eta(jm))
            l(j)=b0(j)+b1(j)*at(j)
            d(j)=vis*cd0
         else
            rey=c(j)*Re
            if(rey.lt.1.0e5)rey=1.0e5
            cd0=1.328/sqrt(rey)
            if(rey.gt.1.0e5)cd0=.072/rey**.2
            sum0=sum0+c(j)*d(j)*(eta(j)-eta(jm))
         endif
 14   continue
      cdi=-0.5*arm*sum
      cdv=vis*(0.25*arm*sum0+2.0*pi*rf*cxm*cd0/am)
      if(abs(cdi).lt.eps)then
         em=1.0
      else
         em=cl*cl/(pi*arm*cdi)
      endif
      cd=cdi+cdv
c*****distributions
      write(6,1004)
      do 15 j=1,jx
         write(6,1003)y(j),c(j),t(j),dem(j),g(j),w(j),l(j),d(j),polar(j)
         write(18,*)y(j),g(j),w(j)
         write(19,*)y(j),g(j)
         write(20,*)y(j),w(j)
         write(21,*)y(j),l(j)
         write(22,*)y(j),at(j)
         write(23,*)y(j),d(j)
         write(24,*)y(j),eta(j),c(j),t(j),dem(j),g(j),w(j)
 15   continue
c*****force and moment
      do 16 j=1,jx
         fz(j)=0.5*arm*fz(j)
         cmf(j)=0.5*arm*cmf(j)
         cmt(j)=0.25*arm*cmt(j)/cam
         write(26,*)y(j),cmf(j)
         write(27,*)eta(j),fz(j)
         write(28,*)eta(j),cmt(j)
 16   continue
      write(6,*)
      write(6,*)'********'
      write(6,*)'results:'
      write(6,1005)cdi,em
      write(6,1006)cdv
      write(6,1007)cd,cl,clf
      write(6,1008)cm0,cmac,xac
      write(6,1011)xcp
      write(6,1009)-cmf(jx2+1)
      write(6,1010)-cmt(jx2+1)
      dum=0.
      write(25,*)alphad,cl,cd,dum,cmac
      write(30,*)alphad,cl/cd
      write(6,*)
      write(6,*)'at(j)=',(at(j),j=1,jx)
      write(6,*)
 400  continue
      y(47)=-0.11111
      xle(47)=1.3
      y(48)=-0.11110
      xle(48)=cxm
      y(54)=0.11110
      xle(54)=cxm
      y(55)=0.11111
      xle(55)=1.3
      write(34,*)y(47),xle(47)
      write(34,*)y(48),xle(48)
      write(34,*)y(54),xle(54)
      write(34,*)y(55),xle(55)
      write(35,*)cd,cl
      write(6,*)
c*****files
      write(6,*)'******data file:'
      write(6,*)'polarbl.dat           :profile polar from Xfoil'
      write(6,*)'prandtline.data       :geometric data for wing'
      write(6,*)'prandtline.in         :restart file'
      write(6,*)'******output files:'
      write(6,*)'polarbl.cdcl          :profile polar cd,cl,cq,inc'
      write(6,*)'prandtline.ycdt       :y,c(y),d(y) and t(y)'
      write(6,*)'prandtline.ygw        :y,g(y) and w(y)'
      write(6,*)'prandtline.yg         :y,g(y) circulation gama'
      write(6,*)'prandtline.yw         :y,w(y) downwash'
      write(6,*)'prandtline.ycl        :y,cl(y)local lift coefficient' 
      write(6,*)'prandtline.yat        :y,at(y) angle of attack'
      write(6,*)'prandtline.ycd        :y(j),cd(j)'
      write(6,*)'prandtline.clcdcq     :alpha,CL,CD,CM,0'
      write(6,*)'prandtline.out        :save results for restart'
      write(6,*)'prandtline.ycmf       :y,cmf(y)'
      write(6,*)'prandtline.yfz        :y,fz(y)'
      write(6,*)'prandtline.ycmt       :y,cmt(y)'
      write(6,*)'prandtline.yxlexte    :y,xle(y),xte(y)'
      write(6,*)'prandtline.incld      :inc,L/D'
      write(6,*)'prandtline.itres      :iter,alog(residue)' 
 1000 format(18a4)
 1001 format(1x,i4,4x,f8.4,4x,f8.4,4x,f8.4,4x,f8.4)
 1002 format(7x,'y(j)=',5x,'eta(j)=',7x,'c(j)=',7x,'t(j)='
     &     ,7x,'d(j)=',7x,'g(j)=',7x,'w(j)=',6x,'at(j)=',3x,'polar(j)=')
 1003 format(8f12.4,i12)
 1004 format(7x,'y(j)=',7x,'c(j)=',7x,'t(j)=',7x,'d(j)=',7x,'g(j)=',
     &     7x,'w(j)=',7x,'  Cl=',7x,'  Cd=',3x,'polar(j)=')
 1005 format(' inviscid contribution: CDi=',f10.6
     &     ,'    Oswald efficiency e=',f10.4)
 1006 format('  viscous contribution: CDv=',f10.6)
 1007 format('        global results:  CD=',f10.6
     &     ,'    lift coefficient CL=',f10.4
     &     ,'   fuselage lift CLf=',f10.4)
 1008 format('                             '
     &     ,'pitching moment coefficient CM,0=',f10.4
     &     ,'               CM,ac=',f10.4
     &     ,/,'                                      '
     &     ,'aerodynamic center x,ac=',f10.4)
 1009 format('                             '    
     &     ,'  root bending moment coef. CM,x=',f10.4)
 1010 format('                             '
     &     ,'  root torsion moment coef. CM,y=',f10.4)
 1011 format('                                      '
     &     ,'center of pressure x,cp=',f10.4)
      end

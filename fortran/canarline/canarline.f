      program canarline
      implicit none
      integer jxx,lxx,nxx,ipolar,nx,n,km,kfirst,k,kdum,ice,kp,jx
      integer jx2,is,icanar,j,nsteps,ivis,nstep,iter,it,mj,jdx,jm
      integer itx,jxs2,jc,ixx,i,ix,ixw
      parameter(ixx=201,jxx=102,lxx=102,nxx=10)
      real eps,pi,degrad,prod,dcz,dcxm,dcxp,incd,si,omega,avis
      real B,cxc,dc,tcd,Rho,Vinf,Amu,alphad,tc,alpha,Re,acdum
      real cacdum,dtet,tetj,yj,etaj,etajm,ac,cac,arc,alphain
      real alphafi,alstep,vis,cxj,czj,qj,dgx,sum,wj,atj,attj
      real reg,res0,alogres,cl,cm0,xac,cmac,cd0,sum0,sum1,sum2
      real rey,cdi,cdv,ec,cd,dum,xcp,Bc0,bc,Cc0,cl0,cl1
      real Rf0,rf,phij,phi0,dwkj,Lambd,lamb,dClcda0,arceff
      real Dx0,xi,str,dxm,Lf0,lf,Zc0,zcanar,xcim,zcim,zwake,eceff
      real Zw0,zwing,xacmstr,yn,wn,aleqd,aleq,theqd,theq
      real c(jxx),g(jxx),dg(jxx),y(jxx),eta(jxx)
      real w(jxx),t(jxx),dec(jxx)
      real a0(jxx),a1(jxx),b0(jxx),b1(jxx),c0(jxx),c1(jxx)
      real l(jxx),d(jxx),q(jxx),at(jxx)
      real cx(lxx,nxx),cz(lxx,nxx),cq(lxx,nxx),inc(lxx,nxx)
      real cmf(jxx),cmt(jxx),fz(jxx)
      real xle(jxx),xte(jxx),wcanar(jxx),xacc(jxx),xiac(jxx),xacw(jxx)
      real xc(ixx),zc(ixx)
      real rbreak(nxx)
      integer m(jxx),polar(jxx)
      integer kx(nxx),kxtrm(lxx,nxx),mxtrm(nxx)
      character*4 bry(18)
      character*4 title(18)
      character*4 typcode(18)
      data a0/jxx*0./a1/jxx*0./b0/jxx*0./b1/jxx*0./c0/jxx*0./c1/jxx*0./
      data m/jxx*0/
      data rbreak/nxx*2./
      open(unit=11,file='wing.yxlexte',form='formatted')
      open(unit=12,file='geocanard.xzmses',form='formatted')
      open(unit=13,file='polarbl.dat',form='formatted')
      open(unit=14,file='polarbl.cdcl',form='formatted')
      open(unit=15,file='canarline.data',form='formatted')
      open(unit=16,file='canarline.in',form='formatted')
      open(unit=17,file='canarline.ycdt',form='formatted')
      open(unit=18,file='canarline.ygw',form='formatted')
      open(unit=19,file='canarline.yg',form='formatted')
      open(unit=20,file='canarline.yw',form='formatted')
      open(unit=21,file='canarline.ycl',form='formatted')
      open(unit=22,file='canarline.yat',form='formatted')
      open(unit=23,file='canarline.ycd',form='formatted')
      open(unit=24,file='canarline.out',form='formatted')
      open(unit=25,file='canarline.clcdcq',form='formatted')
      open(unit=26,file='canarline.ycmf',form='formatted')
      open(unit=27,file='canarline.yfz',form='formatted')
      open(unit=28,file='canarline.ycmt',form='formatted')
      open(unit=29,file='canarline.yxlexte',form='formatted')
      open(unit=30,file='canarline.incld',form='formatted')
      open(unit=31,file='canarline.itres',form='formatted')
      open(unit=32,file='canarwash.ylwl',form='formatted')
      open(unit=33,file='canarline.etaxiac',form='formatted')
      open(unit=34,file='canarline.yxfuse',form='formatted')
      open(unit=35,file='prandtline.cdcl',form='formatted')
      open(unit=36,file='canarwake.xz',form='formatted')
      open(unit=37,file='canaredge.xy',form='formatted')
c*****constants
      eps=1.e-7
      pi=2.*asin(1.)
      degrad=pi/180.
c*****polar data
      write(6,*)
      write(6,*)'******do you want to use polar data? Y/N=1/0'
      read(5,*)ipolar
      if(ipolar.ne.1)goto 5
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
            if(n.ge.nx)then
               rbreak(n)=1.+eps
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
      write(6,*)'                      Y=0.5*Bc0*y'
      write(6,*)'                      C=0.5*Bc0*c'
      write(6,*)'                      A=0.25*Bc0**2*ac'
      write(6,*)'                      D=C*d'
      write(6,*)'                      W=U*w'
      write(6,*)'                   GAMA=0.5*U*Bc0*g'
      write(6,*)'                   LIFT=0.5*RHO*U**2*A*Cl'
      write(6,*)'                   DRAG=0.5*RHO*U**2*A*Cd'
      write(6,*)'               MOMENT,0=0.5*RHO*U**2*A*Cac*Cm0'
      write(6,*)'               REYNOLDS=RHO*U*C/AMU'
      write(6,*)'                     Fz=0.5*RHO*U**2*A*fz'
      write(6,*)'                     Mf=0.25*RHO*U**2*A*Bc0*cmf'
      write(6,*)'                     Mt=0.5*RHO*U**2*A*Cac*cmt'
c*****read in data
      read(15,*)jxs2
      jx=2*jxs2
      if(jx.gt.jxx)then
         write(6,*)'jx=',jx,'>jxx=',jxx,', change dimension:exiting!'
         stop
      endif
      read(15,*)itx
      read(15,*)omega
      read(15,*)avis
      read(15,*)tcd
      read(15,*)theqd
      read(15,*)Bc0
      read(15,*)Cc0
      read(15,*)Zc0
      read(15,*)Lambd
      read(15,*)Rf0
      read(15,*)dc
      read(15,*)arceff
      read(15,*)icanar
      read(15,*)B
      read(15,*)xacmstr
      read(15,*)Zw0
      read(15,*)Lf0
      read(15,*)Dx0
      read(15,*)str
      read(15,*)Rho
      read(15,*)Vinf
      read(15,*)Amu
c     reference Bc0
      tc=degrad*tcd
      theq=degrad*theqd
      bc=2.0
      cxc=2.0*Cc0/Bc0
      zcanar=2.0*Zc0/B
      lamb=degrad*Lambd
      rf=2.0*Rf0/Bc0
      zwing=2.0*Zw0/B
      lf=2.0*Lf0/B
      dxm=2.0*Dx0/B
c*****initialization
c*****mesh, geometry and flow
      write(6,*)' '
      write(6,*)'******point distribution, geometry and flow:'
      write(6,1002)
      acdum=0.
      cacdum=0.
      dtet=pi/(jxs2-1)
      n=1
      do 6 j=1,jxs2
         tetj=(j-1)*dtet
         yj=-1.0+(1.0-cos(tetj))*(1.0-rf)/2.0
         y(j)=yj
         etaj=-1.0+(1.0-cos(tetj+.5*dtet))*(1.0-rf)/2.0
         if(j.eq.1)then
            etajm=-1.0
         else
            etajm=eta(j-1)
         endif
         if(j.eq.jxs2)then
            etaj=y(j)
         endif
         eta(j)=etaj
         dec(j)=dc
         t(j)=tc
         g(j)=0.
         w(j)=0.
         at(j)=0.
         a0(j)=-pi*dc
         a1(j)=0.
         b0(j)=2.*pi*(2.*dec(j))
         b1(j)=2.*pi
c*****elliptic canard
         if(icanar.eq.0)then
            c(j)=cxc*sin(tetj)
            xacc(j)=0.25*cxc
            xacc(j)=xacc(j)+tan(lamb)*abs(y(j))
            xle(j)=xacc(j)-0.25*c(j)
            xte(j)=xle(j)+c(j)
            if(j.ge.2)xiac(j-1)=0.5*(xacc(j)+xacc(j-1))
            acdum=acdum+c(j)*(etaj-etajm)
            cacdum=cacdum+c(j)**2*(etaj-etajm)
         endif
c*****rectangular canard
         if(icanar.eq.1)then
            c(j)=cxc
            xle(j)=0.0
            xte(j)=cxc
            xacc(j)=0.5*cxc
            xiac(j)=0.5*cxc
            acdum=acdum+c(j)*(etaj-etajm)
            cacdum=cacdum+c(j)**2*(etaj-etajm)
         endif
c******canard geometry
         if(icanar.eq.2)then
            if(abs(yj).ge.rf-eps)then
               xle(j)=0.8497-0.466*(1.0-abs(y(j)))
               xte(j)=0.9747-0.266*(1.0-abs(y(j)))
               c(j)=xte(j)-xle(j)
               xacc(j)=xle(j)+0.25*c(j)
               if(j.ge.2)then
                  xiac(j-1)=0.5*(xacc(j)+xacc(j-1))
               endif
               acdum=acdum+c(j)*(etaj-etajm)
               cacdum=cacdum+c(j)**2*(etaj-etajm)
	       a0(j)=0.
            endif
         endif
         if(ipolar.eq.1)then
            if(y(j).gt.rbreak(n)-eps)then
               write(6,*)'rbreak(n)=',rbreak(n)
               n=n+1
            endif
         else
            n=0
         endif
         polar(j)=n
         write(6,1003)y(j),eta(j),c(j),t(j),dec(j),g(j),w(j),at(j)
     &        ,polar(j)
         write(17,*)y(j),c(j),dec(j),t(j)
         write(29,*)y(j),xle(j),xte(j),xacc(j)
         if(j.ge.2)write(33,*)eta(j-1),xiac(j-1)
         y(jx+1-j)=-y(j)
         if(j.lt.jxs2)eta(jx-j)=-eta(j)
 6    continue
      eta(jxs2)=0.0
      eta(jx)=eta(jx-1)
      do 61 jc=1,jxs2
         j=jxs2+jc
         xle(j)=xle(jxs2+1-jc)
         xte(j)=xte(jxs2+1-jc)
         c(j)=c(jxs2+1-jc)
         xacc(j)=xle(j)+0.25*c(j)
         if(j.ge.jxs2+2)xiac(j-1)=0.5*(xacc(j)+xacc(j-1))
         write(29,*)y(j),xle(j),xte(j),xacc(j)
         if(j.ge.jxs2+2)write(33,*)eta(j-1),xiac(j-1)
 61   continue
      xiac(jxs2)=xiac(jxs2-1)
      xiac(jx)=xiac(jx-1)
      write(6,*)' '
      ac=acdum
      cac=cacdum/ac
      arc=0.5*bc**2/ac
      Re=Rho*Vinf*Cc0*cac/Amu
      write(6,*)'***************'
      write(6,*)'numerical data:'
      write(6,*)'          number of points jx=',jx
      write(6,*)' max number of iterations itx=',itx
      write(6,*)'                        omega=',omega
      write(6,*)'   viscosity coefficient avis=',avis
      write(6,*)'     canard setting angle tcd=',tcd,' (deg) ='
     &     ,tc,' (rd)'
      write(6,*)' airplane setting angle theqd=',theqd,' (deg) ='
     &     ,theq,' (rd)'
      write(6,*)'canard data:'
      write(6,*)'              canard span Bc0=',Bc0,' (m)'
      write(6,*)'        canard root chord Cc0=',Cc0,' (m)'
      write(6,*)'        canard z location Zc0=',Zc0,' (m)'
      write(6,*)'   0/1 a.c. sweep angle Lambd=',Lambd,'(deg)'
      write(6,*)'          fuselage radius Rf0=',Rf0,' (m)'
      write(6,*)'    relative camber height dc=',dc,' (ref. C)'
      write(6,*)'     canard efficiency arceff=',arceff
     &     ,' initial value 1 to be recalculated'
      write(6,*)'           canard shape 0/1/2=',icanar
      write(6,*)'                  wing span B=',B,' (m)'
      write(6,*)'wing strake end x-loc xacmstr=',xacmstr,' (ref. B/2)'
      write(6,*)'wing  z-location wrt fuse Zw0=',Zw0,' (m)'
      write(6,*)'          fuselage length Lf0=',Lf0,' (m)'
      write(6,*)'    inital wake mesh step Dx0=',Dx0,' (m)'
      write(6,*)'  wake stretch paprameter str=',str
      write(6,*)'air data:'
      write(6,*)'              air density Rho=',Rho,' (kg/m**3)'
      write(6,*)'           wind velocity Vinf=',Vinf,' (m/s)'
      write(6,*)'        dynamic viscosity Amu=',Amu,' (kg/(m*s))'
      write(6,*)'        reference Reynolds Re=',Re
      write(6,*)'****************'
      write(6,*)'calculated data:'
      write(6,*)'        single canard area ac=',ac,' (ref. Bc0**2/4)'
      write(6,*)'theo. canard aspect ratio arc=',arc
      write(6,*)'average aerodynamic chord cac=',cac,' (ref. Bc0/2)'
      write(6,*)' '
c*****iterations
      write(6,*)' alphain alphafi alstep'
     &     ,'or aleqd,aleqd,0 from canareq.f for wake calculation?'
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
      alpha=degrad*alphad+eps
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
      do 9 j=2,jxs2-1
         n=polar(j)
         atj=at(j)
         mj=2
         prod=1.57-atj
         do 7 k=2,kx(n)-1
            prod=prod*(inc(k,n)-atj)
            if(prod.ge.-eps)goto 8
            prod=1.
            mj=k
 7       continue
 8       continue
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
         l(jx+1-j)=l(j)
         d(jx+1-j)=d(j)
         q(jx+1-j)=q(j)
 9    continue
      m(1)=m(2)
      l(1)=l(2)+(l(3)-l(2))*(y(1)-y(2))/(y(3)-y(2))
      d(1)=d(2)+(d(3)-d(2))*(y(1)-y(2))/(y(3)-y(2))
      q(1)=q(2)+(q(3)-q(2))*(y(1)-y(2))/(y(3)-y(2))
      m(jxs2)=m(jxs2-1)
      l(jxs2)=l(jxs2-1)+(l(jxs2-1)-l(jxs2-2))*(y(jxs2)-y(jxs2-1))
     &     /(y(jx-1)-y(jx-2))
      d(jxs2)=d(jxs2-1)+(d(jxs2-1)-d(jxs2-2))*(y(jx)-y(jx-1))
     &     /(y(jx-1)-y(jx-2))
      q(jxs2)=q(jxs2-1)+(q(jxs2-1)-q(jxs2-2))*(y(jx)-y(jx-1))
     &     /(y(jx-1)-y(jx-2))
      l(jxs2+1)=l(jxs2)
      d(jxs2+1)=d(jxs2)
      q(jxs2+1)=q(jxs2)
      l(jx)=l(1)
      d(jx)=d(1)
      q(jx)=q(1)
 10   continue
c     fixed point iteration
      g(1)=0.
      g(jxs2)=0.
      dgx=0.
      jdx=0
      do 12 j=2,jxs2-1
         sum=0.
         do 11 k=1,jx-1
            if(k.eq.jxs2)goto 11
c     downwash for non-straight lifting line
c     trailed vorticity
            phi0=sign(1.0,y(j)-eps)
     &           *atan((xacc(j)-xiac(k))/(y(j)-eta(k)))
            sum=sum+(g(k+1)-g(k))*(1.0-sin(phi0))/(y(j)-eta(k))
c     downwash due to lifting line
            if(k+1.ne.j.and.k.lt.jx-1)then
               dwkj=-g(k+1)*((xiac(k+1)-xiac(k))*(y(j)-y(k+1))
     &              -(xacc(j)-xacc(k+1))*(eta(k+1)-eta(k)))
     &              /((xacc(j)-xacc(k+1))**2+(y(j)-y(k+1))**2
     &              +10.0*(eta(k+1)-eta(k)))**1.5
               sum=sum+dwkj
            endif
 11      continue
         wj=-sum/(4.*pi)
         atj=alpha+atan2(wj,1.)
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
         if(icanar.eq.2.and.abs(y(j)).le.rf)then
           dg(j)=g(j-1)+2.0*rf*sin(3.0*alpha)/3.0*((1.0-(y(j)/rf)**2)
     &           /(1.0+(y(j)/rf)**2)-(1.0-(y(j-1)/rf)**2)
     &           /(1.0+(y(j-1)/rf)**2))-g(j)
         endif
         g(j)=g(j)+omega*dg(j)
         if(abs(dgx).lt.abs(dg(j)))then
            dgx=dg(j)
            jdx=j
         endif
         g(jx+1-j)=g(j)
         w(jx+1-j)=w(j)
         at(jx+1-j)=at(j)
 12   continue
      g(jxs2+1)=0.
      g(jx)=0.
      if(iter.eq.1)res0=dgx
      alogres=alog10(abs(dgx/res0)+eps)
      write(31,*)iter,alogres
      w(1)=w(2)
     &     +(w(3)-w(2))*(y(1)-y(2))/(y(3)-y(2))
      w(jxs2)=w(jxs2-1)
     &     +(w(jxs2-2)-w(jxs2-1))*(y(jxs2)-y(jxs2-1))
     &     /(y(jxs2-2)-y(jxs2-1))
      w(jxs2+1)=w(jxs2)
      w(jx)=w(1)
      at(1)=at(2)
     &     +(at(3)-at(2))*(y(1)-y(2))/(y(3)-y(2))
      at(jxs2)=at(jxs2-1)
     &     +(at(jxs2-2)-at(jxs2-1))*(y(jxs2)-y(jxs2-1))
     &     /(y(jxs2-2)-y(jxs2-1))
      at(jxs2+1)=at(jxs2)
      at(jx)=at(1)
      if(abs(dgx).lt.eps)goto 300
 200  continue
      write(6,*)'************NOT CONVERGED!!!!!'
 300  continue
c*****results
      write(6,*)
      write(6,*)'alphad=',alphad,' (deg)','   alpha=',alpha,' (rd)'
      write(6,*)'  iter=',iter,'                dgx=',dgx,' jdx=',jdx
      write(6,*)
      write(6,*)'m(j)=',(m(j),j=1,jxs2)
      cl=0.
      cm0=0.
      xac=0.
      cmac=0.
      fz(1)=0.
      cmf(1)=0.
      cmf(2)=0.
      cmt(1)=0.
      do 13 j=2,jxs2-1
         cl=cl+g(j)*(eta(j)-eta(j-1))
         if(ipolar.ne.1)then
            at(j)=alpha+atan2(w(j),1.)+t(j)
            q(j)=c0(j)+c1(j)*at(j)
            xac=xac+xacc(j)*c(j)*(eta(j)-eta(j-1))
            cmac=cmac+c(j)**2*q(j)*(eta(j)-eta(j-1))
            cmt(j)=cmt(j-1)+c(j)**2*q(j)*(eta(j)-eta(j-1))
         else
            xac=xac+xacc(j)*c(j)*(eta(j)-eta(j-1))
            cmac=cmac+c(j)**2*q(j)*(eta(j)-eta(j-1))
         endif
            fz(j)=fz(j-1)+g(j)*(eta(j)-eta(j-1))
            cmt(j)=cmt(j-1)+c(j)**2
     &           *(q(j)+(xacc(j+1)-xacc(j))*fz(j)/cac)
     &           *(eta(j)-eta(j-1))
            cmf(j+1)=cmf(j)-fz(j)*(eta(j+1)-eta(j))
 13   continue
      fz(jxs2)=fz(jxs2-1)
     &     +(fz(jxs2-2)-fz(jxs2-1))*(y(jxs2)-y(jxs2-1))
     &     /(y(jxs2-2)-y(jxs2-1))
      cmt(jxs2)=cmt(jxs2-1)
     &     +(cmt(jxs2-2)-cmt(jxs2-1))*(y(jxs2)-y(jxs2-1))
     &     /(y(jxs2-2)-y(jxs2-1))
      cl=arc*cl
      if(alphad.eq.0.)then
          cl0=cl+eps
      endif
      if(alphad.eq.1.)then
          cl1=cl
      endif
      xac=xac/ac
      cmac=cmac/(ac*cac)
      cm0=cmac-xac*cl/cac
      xcp=xac-cac*cmac/cl
      cd0=0.
      sum=0.
      sum0=0.
      sum1=0.
      sum2=0.
      do 14 j=1,jxs2
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
            sum1=0.
            sum2=0.
            l(j)=b0(j)+b1(j)*at(j)
            d(j)=vis*cd0
         else
            sum0=sum0+c(j)*d(j)*(eta(j)-eta(jm))
            sum1=0.
            sum2=0.
         endif
 14   continue
      cdi=-arc*sum
      cdv=vis*0.25*arc*(sum0+sum1+sum2)
      if(abs(cdi).lt.eps)then
         ec=1.0
      else
         ec=cl*cl/(pi*arc*cdi)
      endif
      cd=cdi+cdv
c*****distributions
      write(6,1004)
      do 15 j=1,jx
         write(6,1003)y(j),c(j),t(j),dec(j),g(j),w(j),l(j),d(j),polar(j)
         write(18,*)y(j),g(j),w(j)
         write(19,*)y(j),g(j)
         write(20,*)y(j),w(j)
         write(21,*)y(j),l(j)
         write(22,*)y(j),at(j)
         write(23,*)y(j),d(j)
         write(24,*)y(j),eta(j),c(j),t(j),dec(j),g(j),w(j)
 15   continue
c*****force and moment
      do 16 j=1,jxs2
         fz(j)=0.5*arc*fz(j)
         cmf(j)=0.5*arc*cmf(j)
         cmt(j)=0.25*arc*cmt(j)/cac
         write(26,*)y(j),cmf(j)
         write(27,*)eta(j),fz(j)
         write(28,*)eta(j),cmt(j)
 16   continue
      write(6,*)
      write(6,*)'********'
      write(6,*)'results:'
      write(6,1005)cdi,ec
      write(6,1006)cdv
      write(6,1007)cd,cl
      write(6,1008)cm0,cmac,xac
      write(6,1011)xcp
      write(6,1009)-cmf(jxs2)
      write(6,1010)-cmt(jxs2)
      dum=0.
      write(25,*)alphad,cl,cd,dum,cmac
      write(30,*)alphad,cl/cd
      write(6,*)
      write(6,*)'at(j)=',(at(j),j=1,jx)
      write(6,*)
      y(47)=-0.11111
      xle(47)=1.3
      y(48)=-0.11110
      xle(48)=cxc
      y(54)=0.11110
      xle(54)=cxc
      y(55)=0.11111
      xle(55)=1.3
      write(34,*)y(47),xle(47)
      write(34,*)y(48),xle(48)
      write(34,*)y(54),xle(54)
      write(34,*)y(55),xle(55)
      write(35,*)cd,cl
 400  continue
      write(6,*)
      write(6,*)'run alpha=0 to 1 deg to get effective'
     &     ,' parameters of canard dClcda0, arceff and eceff'
      if(alstep.lt.eps)then
         write(6,*)'********************************'
     &        ,'not enough info to calculate dClcda0, arceff and eceff'
      else
         dClcda0=(cl1-cl0)/degrad
         arceff=2.*dClcda0/(2.*pi-dClcda0)
         eceff=cl**2/(pi*arceff*cdi)
         write(6,*)'dClcda0=',dClcda0,'arceff=',arceff,'eceff=',eceff
      endif
c     canard profile
      do 17 i=1,ixx
         read(12,*,end=18)xc(i),zc(i)
         ix=i
 17   continue
 18   continue
      write(6,*)' ix=',ix
c     integration of wake trajectory in airplane coordinates (ref Bc0/2)
      ixw=1+alog(1.0+(str-1.0)*((lf-cxc)/dxm))/alog(str)
      write(6,*)'ixw=',ixw
      do 19 i=1,ix
         xc(i)=Bc0*(xle(jxs2)+cxc*xc(i))/B
         zc(i)=Bc0*cxc*zc(i)/B
 19   continue
      xcim=0.0
      zcim=0.0
      dtet=pi/(jxs2-1)
      do 21 i=ix+1,ix+ixw
         xi=dxm*(1.0-str**(i-ix))/(1.0-str)
         sum=0.
         do 20 j=1,jxs2-1
            tetj=(j-1)*dtet
            sum=sum+dtet/sqrt(1.+(cos(tetj)/xi)**2)
 20      continue
         xc(i)=xi
         sum=theq-cl*(1.0+sum/pi)/(pi*arc)
         zc(i)=zcim+sum*(xc(i)-xcim)
         xcim=xc(i)
         zcim=zc(i)
c     canard coordinates
 21   continue
      do 22 i=1,ix
         zc(i)=zcanar+zc(i)
         write(36,*)xc(i),zc(i)
         zc(i)=zc(i)-zcanar
 22   continue
c     attach wake to the canard trailing edge
      zwake=0.
      do 23 i=ix+1,ix+ixw
         xc(i)=xc(ix)+xc(i)
         zc(i)=zcanar+zc(i)
         if(xc(i-1).lt.xacmstr.and.xc(i).gt.xacmstr)then
            zwake=zc(i)-zwing
            write(6,*)'i=',i,'xacmstr=',xacmstr,'zc(i)=',zc(i)
     &           ,'zcanar=',zcanar,'zwake=',zwake
         endif
         write(36,*)xc(i),zc(i)
 23   continue
      write(6,*)'zwake=',zwake,' (ref. B/2)'
c     downwash on wing in airplane coordinate
      nx=101
      dtet=pi/(nx-1)
      do 25 n=1,nx
         read(11,*)dum,dum,dum,xacw(n)
         tetj=(n-1)*dtet
         yn=-cos(tetj)
         sum=0.
         do 24 k=1,jx-1
            if(k.lt.jxs2)then
               phi0=(xacw(n)-Bc0*xiac(k)/B)
     &           /sqrt(zwake**2+(yn-Bc0*eta(k)/B)**2)
            else
               phi0=(xacw(n)-Bc0*xiac(jx+1-k)/B)
     &           /sqrt(zwake**2+(yn+Bc0*eta(jx+1-k)/B)**2)
            endif
            phi0=atan(phi0)
            sum=sum+(g(k+1)-g(k))*(1.0+sin(phi0))*(yn-Bc0*eta(k)/B)
     &           /(zwake**2+(yn-Bc0*eta(k)/B)**2)
 24      continue
         wn=-sum/(2.*pi)
         write(32,*)yn,wn
 25   continue
      write(6,*)
c     point A
      dum=-1.
      acdum=0.8497
c     point B
      cacdum=0.9747
      write(37,*)dum,acdum
      write(37,*)dum,cacdum
c     point C
      dum=-0.25
      cacdum=0.7752
c     point D
      acdum=0.5
      write(37,*)dum,cacdum
      write(37,*)dum,acdum
      write(37,*)-dum,acdum
      write(37,*)-dum,cacdum
c     point E
      dum=0.25
      acdum=0.5
c     point F
      cacdum=0.7752
      write(37,*)dum,acdum
      write(37,*)dum,cacdum
c     point G
      dum=1.
      acdum=0.9747
c     point H
      cacdum=0.8497
      write(37,*)dum,acdum
      write(37,*)dum,cacdum
c*****files
      write(6,*)'******data file:'
      write(6,*)'polarbl.dat           :profile polar from Xfoil'
      write(6,*)'prandtline.data       :geometric data for wing'
      write(6,*)'prandtline.in         :restart file'
      write(6,*)'******output files:'
      write(6,*)'polarbl.cdcl          :profile polar cd,cl,cq,inc'
      write(6,*)'canarline.ycdt        :y,c(y),d(y) and t(y)'
      write(6,*)'canarline.ygw         :y,g(y) and w(y)'
      write(6,*)'canarline.yg          :y,g(y) circulation gama'
      write(6,*)'canarline.yw          :y,w(y) downwash'
      write(6,*)'canarline.ycl         :y,cl(y)local lift coefficient' 
      write(6,*)'canarline.yat         :y,at(y) angle of attack'
      write(6,*)'canarline.ycd         :y(j),cd(j)'
      write(6,*)'canarline.clcdcq      :alpha,CL,CD,CM,0'
      write(6,*)'canarline.out         :save results for restart'
      write(6,*)'canarline.ycmf        :y,cmf(y)'
      write(6,*)'canarline.yfz         :y,fz(y)'
      write(6,*)'canarline.ycmt        :y,cmt(y)'
      write(6,*)'canarline.yxlexte     :y,xle(y),xte(y)'
      write(6,*)'canarline.incld       :inc,L/D'
      write(6,*)'canarline.itres       :iter,alog(residue)' 
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
     &     ,'    lift coefficient CL=',f10.4)
 1008 format('                             '
     &     ,'pitching moment coefficient CM,0=',f10.4,' CM,ac=',f10.4
     &     ,/,'                                      '
     &     ,'aerodynamic center x,ac=',f10.4)
 1009 format('                             '    
     &     ,'  root bending moment coef. CM,x=',f10.4)
 1010 format('                             '
     &     ,'  root torsion moment coef. CM,y=',f10.4)
 1011 format('                                      '
     &     ,'center of pressure x,cp=',f10.4)
      end

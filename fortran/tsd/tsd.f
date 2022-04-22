      program tsd
      implicit none
      integer ixx,kxx,ikxx,ithick,ile,klo,kup,kpt,kx,iprof,ite
      integer ipt,ix,k,kc,i,ic,iter,inflow,itx,idx,kdx,iwrite
      integer it,jxx,ijkxx,j,jx,ijxx,ixdum,jxdum,kxdum,idum,inprof
      integer jtip,jpt,jc,jstr,jdx,jj,jtipp,n
      parameter(ixx=201,jxx=41,kxx=101,ijxx=ixx*jxx
     &     ,ijkxx=ixx*jxx*kxx,ikxx=ixx*kxx)
      real pi,eps,gamp,mach0,ucr,bet0,gamach
      real usdpi,degrad,dx0,dz0,xmin,xmax,zmin,zmax,str
      real omega,dm,em,alphad,alpha,zk,xii,dtet
      real um,ui,rex,rex1,tenlog,cd,cm0,cl,cpcr
      real dy0,ymin,ymax,dmplus,dum,lamb,swp,ratio,stop
      real yj,bs2,ystr,piv,dga,yjm,am,cdum,cmum,cav,cdw
      real x(ixx),y(jxx),z(kxx),xi(ixx,jxx)
      real ph(ixx,jxx,kxx)
      real aa(kxx),bb(kxx),cc(kxx),dd(kxx)
      real d(ixx),e(ixx),dp(ixx,jxx),ep(ixx,jxx)
      real cpo(ixx,jxx),cpu(ixx,jxx),gp(ixx,jxx)
      real pho(ixx,jxx),phu(ixx,jxx)
      real cpwo(ixx,jxx),cpwu(ixx,jxx)
      real zu(ixx,jxx),zo(ixx,jxx)
      real cp(ixx,jxx,kxx),u(ixx,jxx,kxx)
      real ax(ixx),ay(ixx),xle(jxx),xte(jxx),c(jxx),ga(jxx)
      real cz(jxx),cx(jxx),cmo(jxx),xcp(jxx)
      data ph/ijkxx*0.0/cp/ijkxx*0.0/u/ijkxx*0.0/
      data cpo/ijxx*0.0/cpu/ijxx*0.0/cpwo/ijxx*0.0/cpwu/ijxx*0.0/
      data gp/ijxx*0.0/zu/ijxx*0.0/zo/ijxx*0.0/dp/ijxx*0.0/ep/ijxx*0.0/
      data ga/jxx*0.0/cz/jxx*0.0/cx/jxx*0.0/cmo/jxx*0.0/xcp/jxx*0.0/
      common/constants/pi,eps,gamp,mach0,ucr,bet0,gamach,piv,cdw
c*****files
      open(unit=14,form='formatted',file='geoprofortsd.xde')
      open(unit=15,form='formatted',file='tsd.data')
      open(unit=16,form='formatted',file='tsd.itr')
      open(unit=17,form='formatted',file='tsd.ix')
      open(unit=18,form='formatted',file='tsd.kz')
      open(unit=19,form='formatted',file='tsd.cpo')
      open(unit=20,form='formatted',file='tsd.cpu')
      open(unit=21,form='formatted',file='tsd.gp')
      open(unit=22,form='formatted',file='tsd.cpwo')
      open(unit=23,form='formatted',file='tsd.cpwu')
      open(unit=24,form='formatted',file='tsd.in')
      open(unit=25,form='formatted',file='tsd.out')
      open(unit=26,form='formatted',file='tsd.cpcr')
      open(unit=27,form='formatted',file='tsd.itcl')
      open(unit=28,form='formatted',file='tsd.cpcon')
      open(unit=29,form='formatted',file='tsd.vcon')
      open(unit=30,form='formatted',file='tsd.xzmses')
      open(unit=31,form='formatted',file='tsd.dom')
      open(unit=32,form='formatted',file='tsd.jy')
      open(unit=33,form='formatted',file='tsd.yxlexte')
      open(unit=34,form='formatted',file='tsd.xymesh1')
      open(unit=35,form='formatted',file='tsd.xymesh2')
c*****constants
      pi=2.0*asin(1.0)
      usdpi=1.0/(2.0*pi)
      degrad=pi/180.0
      eps=1.e-6
c*****definition of mesh system
      read(15,*)dx0
      read(15,*)dy0
      read(15,*)dz0
      read(15,*)xmin
      read(15,*)xmax
      read(15,*)bs2
      read(15,*)ystr
      read(15,*)ymin
      read(15,*)ymax
      read(15,*)zmin
      read(15,*)zmax
      read(15,*)str
      read(15,*)omega
      read(15,*)piv
      read(15,*)dm
      read(15,*)dmplus
      read(15,*)em
      read(15,*)ithick
      read(15,*)inprof
      read(15,*)ratio
      read(15,*)alphad
      read(15,*)gamp
      read(15,*)mach0
      read(15,*)lamb
      if(str.lt.1.0-eps)then
         write(6,*)
         write(6,*)'######parameter str<1.0 exiting'
         write(6,*)
         stop
      endif
      alpha=degrad*alphad
      gamach=gamp*mach0**2
      swp=0.0
c     for infinite swept wing with lamb non-zero
      if(abs(lamb).gt.eps)then
c        2d case
         gamach=gamp*(mach0*cos(pi*lamb/180.0))**2
         swp=mach0*sin(pi*lamb/180.0)
      endif
      bet0=1.0-mach0**2+swp**2
      if(mach0.lt.eps)mach0=eps
      ucr=bet0/gamach
      write(6,*)
      write(6,*)'gamach=',gamach,'  bet0=',bet0,'   ucr=',ucr
      if(str-1.0-eps.lt.0.0)then
         ile=1+eps-xmin/dx0
         klo=1+eps-zmin/dz0
         kup=klo+1
         kpt=1+eps+zmax/dz0
         kx=klo+kpt
      else
         ile=1+0.5*pi*sqrt(-xmin/(2.0*dx0))
         klo=1+alog(1.0+(str-1.0)*(0.5-zmin/dz0))/alog(str)
         kup=klo+1
         kpt=1+alog(1.0+(str-1.0)*(0.5+zmax/dz0))/alog(str)
         kx=klo+kpt
      endif
      write(6,*)
      if(str-1.0-eps.lt.0.0)then
         iprof=eps+1.0/dx0
         ite=ile+iprof
         ipt=(xmax-1.0)/dx0+eps
         ix=ite+ipt
      else
         iprof=1+0.5*pi/sqrt(dx0)
         ite=ile+iprof
         ipt=1+0.5*pi*sqrt((xmax-1.0)/(2.0*dx0))
         ix=ite+ipt
      endif
      if(ix.gt.ixx.or.kx.gt.kxx)then
         write(6,*)'ix=',ix,'kx=',kx,' ix or kx too large: exiting'
         stop
      endif
      if(str-1.0-eps.lt.0.0)then
         jtip=eps+1+(bs2+0.5*dy0)/dy0
      else
         jtip=1+0.25*pi*sqrt(0.5*bs2/dy0)
         jtip=2*jtip
         dtet=pi/(2*jtip-1)
         jtip=jtip+1
      endif
c     set mesh system
      jstr=0
      do 1 j=1,jtip
         if(str-1.0-eps.lt.0.)then
            y(j)=(j-1.5)*dy0
         else
            y(j)=bs2*cos((j-jtip)*dtet)
         endif
         if(y(j).gt.ystr)then
            if(abs(lamb).gt.eps)then
               xle(j)=tan(pi*lamb/180.0)*y(j)
               xte(j)=xle(j)+1.0
            else
               xle(j)=1.0302+0.466*(abs(y(j))-1.7778)
               xte(j)=1.218+0.156*(abs(y(j))-1.7778)
            endif
         else
            if(abs(lamb).gt.eps)then
               xle(j)=tan(pi*lamb/180.0)*y(j)
               xte(j)=xle(j)+1.0
            else
               xle(j)=abs(y(j))
               xte(j)=1.0
            endif
            jstr=j
         endif
         c(j)=xte(j)-xle(j)
         write(32,*)j,y(j)
         write(33,*)y(j),xle(j),xte(j)
 1    continue
      write(6,*)'  jstr=',jstr
      if(str-1.0-eps.lt.0.0)then
         jpt=(ymax-y(jstr))/dy0+eps
      else
         jpt=1+0.5*pi*sqrt((ymax-bs2)/(2.0*dy0))
      endif
      jx=jtip+jpt
      if(jx.gt.jxx)write(6,*)'jx too large, increase jxx, exiting'
     &     ,stop
      dtet=0.5*pi/(jx-jtip)
      do 2 j=jtip+1,jx
         jc=j-jtip
         if(str-1.0-eps.lt.0.0)then
            yj=y(jtip)+jc*dy0
         else
            yj=y(jtip)+(ymax-bs2)*(1.0-cos((j-jtip)*dtet))
         endif
         y(j)=yj
         if(abs(lamb).gt.eps)then
            xle(j)=tan(pi*lamb/180.0)*y(j)
            xte(j)=xle(j)+1.0
         else
            xle(j)=xle(jtip)
            xte(j)=xte(jtip)
         endif
         c(j)=xte(j)-xle(j)
         write(32,*)j,y(j)
         write(33,*)y(j),xle(j),xte(j)
 2    continue
      if(abs(lamb).gt.eps)then
         jx=3
         jtip=jx
      endif
      am=0.0
      yjm=0.5*(y(2)+y(1))
      do 3 j=2,jtip
         if(j.eq.jtip)then
            am=am+c(j)*(y(jtip)-yjm)
         else
            am=am+c(j)*(0.5*(y(j+1)+y(j))-yjm)
         endif
         yjm=0.5*(y(j+1)+y(j))
 3    continue
      am=am
c*****output on listing
      write(6,*)'*************************'
      write(6,*)'dimensionless parameters:'
      write(6,*)'                           C=root chord'
      write(6,*)'                           X=C*x'
      write(6,*)'                           Z=C*z'
      write(6,*)'                           G=U*C*ga'
      write(6,*)'                      P-Pinf=0.5*Rho*V**2*Cp'
      write(6,*)'                           L=0.5*Rho*V**2*Am*Cl'
      write(6,*)'                           D=0.5*Rho*V**2*Am*Cd'
      write(6,*)'                         M,o=0.5*Rho*V**2*Am*Cam*Cm,o'
      write(6,*)
      write(6,*)'****************'
      write(6,*)'mesh parameters:'
      write(6,*)'                         dx0=',dx0,'  dy0=',dy0
     &     ,'   dz0=',dz0
      write(6,*)'                        xmin=',xmin,' xmax=',xmax
      write(6,*)'                        ymin=',ymin,' ymax=',ymax
      write(6,*)'                        zmin=',zmin,' zmax=',zmax
      write(6,*)'                         str=',str
      write(6,*)'                         ile=',ile,'       ite=',ite
     &     ,'    ix=',ix
      write(6,*)'                        jstr=',jstr,'      jtip=',jtip
     &     ,'    jx=',jx
      write(6,*)'                         klo=',klo,'       kup=',kup
     &     ,'    kx=',kx
      write(6,*)'******************'
      write(6,*)'relaxation method:'
      write(6,*)'                       omega=',omega,'piv=',piv
      write(6,*)'******************'
      write(6,*)'main profile data:'
      write(6,*)'               maximum chord=   1.0'
      write(6,*)'               half-span bs2=',y(jtip)
     &     ,' (ref. root chord)'
      write(6,*)'          end of strake ystr=',ystr
     &     ,' (ref. root chord)'
      write(6,*)'                wing area am=',am,' (ref. C**2)'
      write(6,*)'             relative camber=',dm,' <0 for flying wing'
      write(6,*)' added parabolic camb.dmplus=',dmplus
      write(6,*)'          relative thickness=',em
      write(6,*)'      thickness distribution=',ithick,
     &     '     (ellip/semicubic/q-j/naca00xx/selig/s-selig/biconvex)'
      write(6,*)'read-in discrete data inprof=',inprof
      write(6,*)'change thickness ratio ratio=',ratio
      write(6,*)'***********************'
      write(6,*)'aerodynamic parameters:'
      write(6,*)'             angle of attack=',alpha,' (rd) ='
     &     ,alphad,' (deg)'
      write(6,*)'   coefficient (gama+1) gamp=',gamp
      write(6,*)'     incoming Mach number M0=',mach0
      write(6,*)'             wing sweep lamb=',lamb,' (deg)'
      write(6,*)'           tangent(lamb) swp=',swp
      write(6,*)'         1-M0**2+swp**2 bet0=',bet0
      write(6,*)'                         ucr=',ucr
      write(6,*)
c*****initialization, mesh
      do 4 k=1,klo
         kc=klo+1-k
         if(str-1.0-eps.lt.0.0)then
            zk=-0.5*dz0-(kc-1)*dz0
         else
            zk=-0.5*dz0-dz0*(1.0-str**(kc-1))/(1.0-str)
         endif
         z(k)=zk
 4    continue
      do 5 k=kup,kx
         kc=k-klo
         if(str-1.0-eps.lt.0.0)then
            zk=0.5*dz0+(kc-1)*dz0
         else
            zk=0.5*dz0+dz0*(1.0-str**(kc-1))/(1.0-str)
         endif
         z(k)=zk
 5    continue
      dtet=pi/(2.0*(ile-1))
      do 8 i=1,ile
         if(str-1.0-eps.lt.0.0)then
            xii=-(ile-i)*dx0
            x(i)=xii
            do 6 j=1,jx
               xi(i,j)=x(i)+xle(j)
 6       continue
         else
            xii=xmin*(1.0-cos((ile-i)*dtet))
            x(i)=xii
            do 7 j=1,jx
               xi(i,j)=x(i)+xle(j)
 7       continue
         endif
         d(i)=0.0
         e(i)=0.0
         dp(i,2)=0.0
         ep(i,2)=0.0
 8    continue
      dtet=pi/(ite-ile)
      do 11 i=ile+1,ite
         if(str-1.0-eps.lt.0.0)then
            xii=(i-ile)*dx0
            x(i)=xii
            do 9 j=1,jx
               xi(i,j)=xle(j)+x(i)*(xte(j)-xle(j))
 9       continue
         else
            ic=i-ile
            xii=0.5*(1.0-cos(ic*dtet))
            x(i)=xii
            do 10 j=1,jx
               xi(i,j)=xle(j)+x(i)*(xte(j)-xle(j))
 10      continue
         endif
         d(i)=0.0
         e(i)=0.0
         dp(i,2)=0.0
         ep(i,2)=0.0
 11   continue
      dtet=pi/(2.0*(ix-ite))
      do 14 i=ite+1,ix
         if(str-1.0-eps.lt.0.0)then
            xii=1.0+(i-ite)*dx0
            x(i)=xii
            do 12 j=1,jx
               xi(i,j)=x(i)+xte(j)-1.0
 12         continue
         else
            xii=1.0+xmax*(1.0-cos((i-ite)*dtet))
            x(i)=xii
            do 13 j=1,jx
               xi(i,j)=x(i)+xte(j)-1.0
 13         continue
         endif
         d(i)=0.0
         e(i)=0.0
         dp(i,2)=0.0
         ep(i,2)=0.0
 14   continue
      do 16 j=1,jx
         do 15 i=1,ix
            ic=ix+1-i
            if(mod(j,2).eq.1)then
               write(34,*)xi(i,j),y(j)
            else
               write(34,*)xi(ic,j),y(j)
            endif
 15      continue
 16   continue
      do 18 i=1,ix
         do 17 j=1,jx
            jc=jx+1-j
            if(mod(i,2).eq.1)then
               write(35,*)xi(i,j),y(j)
            else
               write(35,*)xi(i,jc),y(jc)
            endif
 17      continue
 18   continue
c     profile geometry
      do 19 i=ile,ite-1
         if(inprof.ne.0)then
            read(14,*)dum,d(i),e(i)
         endif
 19   continue
      dtet=pi/(ite-ile)
      xii=1.0
      zu(ite,2)=0.
      zo(ite,2)=0.
      write(30,*)xii,zo(ite,2)
      do 20 ic=ile,ite-1
         i=ite-1+ile-ic
         xii=0.5*(1.0-cos((i-ile+0.5)*dtet))
         if(inprof.eq.0)then
            call camberdis(1.0,dm,dmplus,xii,d(i))
            call thickdis(ithick,1.0,em,xii,e(i))
         endif
         zu(i,2)=d(i)+e(i)
         zo(i,2)=d(i)-e(i)
         write(30,*)xii,zo(i,2)
 20   continue
      xii=0.0
      zu(ile,2)=0.0
      zo(ile,2)=0.0
      write(30,*)xii,zo(ile,2)
      do 21 i=ile,ite-1
         zu(i,2)=d(i)+e(i)
         zo(i,2)=d(i)-e(i)
         xii=0.5*(1.0-cos((i-ile+0.5)*dtet))
         write(30,*)xii,zu(i,2)
 21   continue
      xii=1.0
      zu(ite,2)=0.
      zo(ite,2)=0.
      write(30,*)xii,zu(ite,2)
      write(6,*)
      write(6,*)'******do you want to write the geometry? Y/N=1/0'
      read(5,*)iwrite
      write(6,*)
      write(6,1003)
      dp(1,2)=0.
      ep(1,2)=0.
      if(iwrite.eq.1)write(6,1004)1,x(1),d(1),dp(1,2),e(1),ep(1,2)
      do 23 i=2,ix-1
         do 22 j=1,jtip
            dp(i,j)=2.0*ratio*(d(i)-d(i-1))/(x(i+1)-x(i-1))
            ep(i,j)=2.0*ratio*(e(i)-e(i-1))/(x(i+1)-x(i-1))
            if(i.eq.ite)then
               dp(i,j)=dp(i-1,j)
               ep(i,j)=ep(i-1,j)
            endif
 22   continue
         if(iwrite.eq.1)write(6,1004)i,x(i),d(i),dp(i,2)
     &        ,e(i),ep(i,2)
 23   continue
      xii=1.0
      dp(ix,2)=0.
      ep(ix,2)=0.
      if(iwrite.eq.1)write(6,1004)ix,x(ix),d(ix),dp(ix,2),e(ix),ep(ix,2)
      write(6,*)
      iter=0
      write(6,*)'******do you want to read-in the flow? Y/N=1/0'
      read(5,*)inflow
      if(inflow.eq.1)then
         write(6,*)
         read(24,*)ixdum,jxdum,kxdum
         write(6,*)'ix=',ixdum,'jx=',jxdum,'kx=',kxdum
         do i=1,ix
            do j=1,jx
                do k=1,kx
                    read(24,*)ph(i,j,k)
                enddo
            enddo
         enddo
c         read(24,*)(((ph(i,j,k),i=1,ix),j=1,jx),k=1,kx)
         read(24,*)iter,rex,(ga(j),j=1,jx)
         write(6,*)
         do i=1,ix
            do j=1,jx
                do k=1,kx
                    write(6,*)'ph = ',ph(i,j,k)
                enddo
            enddo
         enddo
         write(6,*)'(ga(j),j=1,jx)',(ga(j),j=1,jx)
         write(6,*)
         write(6,*)'iter=',iter,' rex=',rex
      endif
      write(6,*)
 100  continue
      write(6,*)'itx=?'
      read(5,*)itx
      if(itx.le.0)goto 400
      do 300 it=1,itx
      iter=iter+1
      cdw=0.0
c*****scheme
c     x-sweep
      rex=0.
      idx=0
      jdx=0
      kdx=0
      write(6,*)'kx = ',kx
c     y-sweep
      do 200 jj=2,jx-1
         if(mod(iter,2).eq.1)then
            j=jj
         else
            j=jx+1-jj
            j=jj
         endif
c*****bc at i=1
      do 24 k=1,kx
         if(mach0.lt.1.0-eps)then
            ph(1,j,k)=usdpi*ga(j)*atan2(z(k),-(xi(1,j)-xle(j)))
     &           /sqrt(bet0)
         else
            ph(1,j,k)=0.0
         endif
 24   continue
      do 27 n=1,1
      do 27 i=2,ix-1
c     bc at k=1 and k=kx
         if(mach0.lt.1.0-eps)then
            aa(1)=0.0
            bb(1)=1.0
            cc(1)=0.0
            dd(1)=0.0
            ph(i,j,1)=usdpi*ga(j)*atan2(z(1),-(xi(i,j)-xle(j)))
     &           /sqrt(bet0)
            aa(kx)=0.0
            bb(kx)=1.0
            cc(kx)=0.0
            dd(kx)=0.0
            ph(i,j,kx)=usdpi*ga(j)*atan2(z(kx),-(xi(i,j)-xle(j)))
     &           /sqrt(bet0)
         else
            aa(1)=0.0
            bb(1)=1.0
            cc(1)=0.0
            dd(1)=0.0
            aa(kx)=0.0
            bb(kx)=1.0
            cc(kx)=0.0
            dd(kx)=0.0
         endif
c     z-sweep interior points
         do 25 k=2,kx-1
            um=(ph(i,j,k)-ph(i-1,j,k))/(xi(i,j)-xi(i-1,j))
            if(i.gt.2)then
               um=0.5*(um
     &           +(ph(i-1,j,k)-ph(i-2,j,k))/(xi(i-1,j)-xi(i-2,j)))
            endif
            u(i-1,j,k)=um
            ui=0.5*((ph(i+1,j,k)-ph(i,j,k))/(xi(i+1,j)-xi(i,j))
     &        +(ph(i,j,k)-ph(i-1,j,k))/(xi(i,j)-xi(i-1,j)))
            u(i,j,k)=ui
c            write(6,*)'i = ',i,' j = ',j,'k = ',k,' ui = ',ui
c            write(6,*)'ph+1',ph(i+1,j,k)
c            write(6,*)'ph',ph(i,j,k)
c            write(6,*)'ph-1',ph(i-1,j,k)
c            write(6,*)'xi+1',xi(i+1,j)
c            write(6,*)'xi',xi(i,j)
            if(k.lt.klo.or.k.gt.kup.or.i.lt.ile.or.j.gt.jtip)then
               call jjscheme(ixx,jxx,kxx,i,j,k,um,ui,xi,y,z,ph
     &              ,aa,bb,cc,dd)
c               write(6,*)'hi 0'
            else
               if(k.eq.klo.and.i.le.ite.and.j.le.jtip)then
               call jjscheme(ixx,jxx,kxx,i,j,k,um,ui,xi,y,z,ph
     &                 ,aa,bb,cc,dd)
                  bb(k)=bb(k)
     &                 -(1.0/(z(k+1)-z(k)))
     &                 *0.5*(xi(i+1,j)-xi(i-1,j))
                  cc(k)=0.0
                  dd(k)=dd(k)
     &                 +(dp(i,j)-ep(i,j)-alpha
     &                 -(ph(i,j,k+1)-ph(i,j,k))/(z(k+1)-z(k)))
     &                 *0.5*(xi(i+1,j)-xi(i-1,j))
c                  write(6,*)'hi 1'
               endif
               if(k.eq.klo.and.i.gt.ite)then
                  call jjscheme(ixx,jxx,kxx,i,j,k,um,ui,xi,y,z,ph
     &                 ,aa,bb,cc,dd)
                  dd(k)=dd(k)
     &                 +(-ga(j)/(z(k+1)-z(k)))
     &                 *0.5*(xi(i+1,j)-xi(i-1,j))
c                  write(6,*)'hi 2'
               endif
               if(k.eq.kup.and.i.le.ite.and.j.le.jtip)then
                  call jjscheme(ixx,jxx,kxx,i,j,k,um,ui,xi,y,z,ph
     &                 ,aa,bb,cc,dd)
                  aa(k)=0.0
                  bb(k)=bb(k)
     &                 -(1.0/(z(k)-z(k-1)))
     &                 *0.5*(xi(i+1,j)-xi(i-1,j))
                  dd(k)=dd(k)
     &                 +(-(dp(i,j)+ep(i,j)-alpha)
     &                 +(ph(i,j,k)-ph(i,j,k-1))/(z(k)-z(k-1)))
     &                 *0.5*(xi(i+1,j)-xi(i-1,j))
c                  write(6,*)'hi 3'
               endif
               if(k.eq.kup.and.i.gt.ite)then
                  call jjscheme(ixx,jxx,kxx,i,j,k,um,ui,xi,y,z,ph
     &                 ,aa,bb,cc,dd)
                  dd(k)=dd(k)
     &                 +(ga(j)/(z(k)-z(k-1)))
     &                 *0.5*(xi(i+1,j)-xi(i-1,j))
c                  write(6,*)'hi 4'
               endif
            endif
            if(ui.lt.ucr-eps)then
               dd(k)=omega*dd(k)
            else
               if(i.gt.2.and.i.lt.ix-1)then
                  dd(k)=dd(k)-0.0*(u(i-1,j,k)-2.0*u(i,j,k)
     &                 +u(i+1,j,k))
               endif
            endif
            if(abs(dd(k)).gt.abs(rex))then
               rex=dd(k)
               idx=i
               jdx=j
               kdx=k
            endif
c            write(6,*)'aa(k)',aa(k)
c            write(6,*)'bb(k)',bb(k)
c            write(6,*)'cc(k)',cc(k)
c            write(6,*)'dd(k)',dd(k)
c            write(6,*)''
 25      continue
         call tridiag(aa,bb,cc,dd,1,kx)
         do 26 k=1,kx
c            write(6,*)''
c            write(6,*)'i = ',i,' j = ',j,'k = ',k,' d* = ',dd(k)
c            write(6,*)'i = ',i,' j = ',j,'k = ',k,' ph* = ',ph(i,j,k)
            ph(i,j,k)=ph(i,j,k)+dd(k)
c            write(6,*)'i = ',i,' j = ',j,'k = ',k,' ph = ',ph(i,j,k)
c            write(6,*)''
 26      continue
         if(i.eq.ite.and.j.le.jtip)then
            dga=ph(ite,j,kup)-ph(ite,j,klo)
     &          +(ph(ite,j,kup+1)-ph(ite,j,kup))*(0.0-z(kup))
     &          /(z(kup+1)-z(kup))
     &          -(ph(ite,j,klo)-ph(ite,j,klo-1))*(0.0-z(klo))
     &          /(z(klo)-z(klo-1))-ga(j)
            ga(j)=ga(j)+omega*dga
         endif
 27   continue
c     bc at i=ix
      do 28 k=1,kx
         ph(ix,j,k)=ph(ix-1,j,k)
 28   continue
 200  continue
c     first j-plane  and last j-plane
      do 30 i=1,ix
         do 29 k=1,kx
            ph(i,1,k)=ph(i,2,k)
            ph(i,jx,k)=ph(i,jx-1,k)
 29      continue
 30   continue
      ga(1)=ga(2)
      ga(jx)=ga(jx-1)
      do 31 j=jtip+1,jx
         ga(j)=0.0
         cx(j)=0.0
         cmo(j)=0.0
 31   continue
c     calculate lift, drag and moment
      if(iter.eq.1)then
         rex1=rex
      endif
      tenlog=log10(abs(rex/rex1)+eps*eps)
      write(16,*)iter,tenlog
      cl=0.0
      yjm=0.0
      do 32 j=2,jtip
         if(j.eq.jtip)then
            cl=cl+ga(j)*(y(jtip)-yjm)
         else
            cl=cl+ga(j)*(0.5*(y(j+1)+y(j))-yjm)
         endif
         yjm=0.5*(y(j+1)+y(j))
 32   continue
      cl=2.0*cl/am
      cdw=-gamach*cdw/(6.0*am)
      write(27,*)iter,cl
 300  continue
      write(6,1000)iter,rex,idx,jdx,kdx,cl,cdw
      goto 100
 400  continue
      do 34 j=1,jx
         do 33 i=1,ix
            pho(i,j)=ph(i,j,klo)+(ph(i,j,klo)-ph(i,j,klo-1))
     &           *(0.0-z(klo))/(z(klo)-z(klo-1))
            phu(i,j)=ph(i,j,kup)+(ph(i,j,kup+1)-ph(i,j,kup))
     &           *(0.0-z(kup))/(z(kup+1)-z(kup))
         if(i.ge.ite.or.j.gt.jtip)then
            pho(i,j)=0.5*(pho(i,j)+phu(i,j))
            phu(i,j)=pho(i,j)+0.5*ga(j)
            pho(i,j)=pho(i,j)-0.5*ga(j)
         endif
         if(j.eq.2)write(17,*)i,x(i)
 33   continue
 34   continue
      do 35 k=1,kx
         write(18,*)k,z(k)
 35   continue
      do 37 j=1,jx
         do 36 i=2,ix-1
            cpo(i,j)=2.0*(pho(i+1,j)-pho(i,j))/(xi(i+1,j)-xi(i,j))
            cpu(i,j)=2.0*(phu(i+1,j)-phu(i,j))/(xi(i+1,j)-xi(i,j))
            gp(i,j)=0.5*(cpu(i,j)-cpo(i,j))
            if(i.lt.ile.or.i.ge.ite.or.j.gt.jtip)then
               cpo(i,j)=0.5*(cpo(i,j)+cpu(i,j))
               cpu(i,j)=cpo(i,j)
               gp(i,j)=0.0
            endif
            cpwo(i,j)=2.0*(ph(i+1,j,1)-ph(i-1,j,1))
     &           /(xi(i+1,j)-xi(i-1,j))
            cpwu(i,j)=2.0*(ph(i+1,j,kx)-ph(i-1,j,kx))
     &           /(xi(i+1,j)-xi(i-1,j))
 36      continue
 37   continue
      do 38 j=1,jx
         cpo(1,j)=cpo(2,j)
         cpu(1,j)=cpu(2,j)
         cpwo(1,j)=cpwo(2,j)+(cpwo(3,j)-cpwo(2,j))*(xi(1,j)-xi(2,j))
     &        /(xi(3,j)-xi(2,j))
         cpwu(1,j)=cpwu(2,j)+(cpwu(3,j)-cpwu(2,j))*(xi(1,j)-xi(2,j))
     &        /(xi(3,j)-xi(2,j))
         cpwo(ix,j)=cpwo(ix-1,j)+(cpwo(ix-1,j)
     &        -cpwo(ix-2,j))*(xi(ix,j)-xi(ix-1,j))
     &        /(xi(ix-1,j)-xi(ix-2,j))
         cpwu(ix,j)=cpwu(ix-1,j)+(cpwu(ix-1,j)
     &        -cpwu(ix-2,j))*(x(ix)-x(ix-1))
     &        /(xi(ix-1,j)-xi(ix-2,j))
         cpo(ix,j)=0.0
         cpu(ix,j)=0.0
         gp(1,j)=0.0
         gp(ix,j)=0.0
 38   continue
      write(25,*)ix,jx,kx
      if(abs(lamb).gt.eps)then
         jtipp=jtip-1
      else
         jtipp=jtip
      endif
      do 40 j=2,jtipp
         do 39 i=1,ix
            ic=ix+1-i
            if(mod(j,2).eq.1)then
               write(19,*)xi(i,j),cpo(i,j)
               write(20,*)xi(i,j),cpu(i,j)
               write(21,*)xi(i,j),gp(i,j)
               write(22,*)xi(i,j),cpwo(i,j)
               write(23,*)xi(i,j),cpwu(i,j)
            else
               write(19,*)xi(ic,j),cpo(ic,j)
               write(20,*)xi(ic,j),cpu(ic,j)
               write(21,*)xi(ic,j),gp(ic,j)
               write(22,*)xi(ic,j),cpwo(ic,j)
               write(23,*)xi(ic,j),cpwu(ic,j)
            endif
            write(40+j,*)xi(i,j),cpo(i,j),cpu(i,j)
 39      continue
 40   continue
      do i=1,ix
        do j=1,jx
            do k=1,kx
                write(25,*)ph(i,j,k)
            enddo
        enddo
      enddo
c      write(25,*)(((ph(i,j,k),i=1,ix),j=1,jx),k=1,kx)
      write(25,*)iter,rex,(ga(j),j=1,jx)
      do 42 k=2,kx-1
         do 41 i=2,ix-1
            j=2
            cp(i,j,k)=-2.0*(ph(i+1,j,k)-ph(i-1,j,k))
     &           /(xi(i+1,j)-xi(i-1,j))
 41      continue
 42   continue
      do 43 k=1,kx
         write(28,*)((cp(i,j,k),i=1,ix),j=2,2)
 43   continue
      cl=0.0
      cav=0.0
      cdum=0.0
      cmum=0.0
      yjm=0.0
      do 45 j=1,jtip
         cd=0.0
         cm0=0.0
         do 44 i=ile,ite
            cd=cd+0.5*((cpo(i-1,j)+cpo(i,j))*(dp(i,j)-ep(i,j))
     &           -(cpu(i-1,j)+cpu(i,j))*(dp(i,j)+ep(i,j)))
     &           *0.5*(xi(i+1,j)-xi(i-1,j))
            cm0=cm0-(gp(i-1,j)+gp(i,j))*xi(i,j)
     &           *0.5*(xi(i+1,j)-xi(i-1,j))
 44      continue
         cz(j)=2.0*ga(j)/c(j)
         cx(j)=cd
         cmo(j)=cm0
         if(abs(cz(j)).gt.eps)then
            xcp(j)=-cm0/cz(j)
         else
            xcp(j)=1.0/eps
         endif
         if(j.eq.jtip)then
            cl=cl+ga(j)*(y(jtip)-yjm)
            cdum=cdum+cx(j)*(y(jtip)-yjm)
            cmum=cmum+cmo(j)*(y(jtip)-yjm)
            cav=cav+c(j)**2*(y(jtip)-yjm)
         else
            cl=cl+ga(j)*(0.5*(y(j+1)+y(j))-yjm)
            cdum=cdum+cx(j)*(0.5*(y(j+1)+y(j))-yjm)
            cmum=cmum+cmo(j)*(0.5*(y(j+1)+y(j))-yjm)
            cav=cav+c(j)**2*(0.5*(y(j+1)+y(j))-yjm)
         endif
         yjm=0.5*(y(j+1)+y(j))
 45   continue
      cl=2.0*cl/am
      cav=cav/am
      cdum=cdum/am
      cmum=cmum/(cav*am)
      write(6,*)
      write(6,*)'******global results:'
      write(6,*)'   M0=',mach0,'alpha=',alphad,' lamb=',lamb
     &     ,'ratio=',ratio,'cav=',cav
      write(6,*)'   Cl=',cl,'   Cd=',cdw
      write(6,*)' Cm,o=',cmum
      write(6,*)
      cpcr=2.0*ucr
      write(26,*)xmin,cpcr
      write(26,*)xmax,cpcr
      write(6,*)
      write(6,*)'(ga(j),j=1,jx)',(ga(j),j=1,jx)
      write(6,*)
      write(6,*)'(cx(j),j=1,jx)',(cx(j),j=1,jx)
      write(6,*)
      write(6,*)'(cmo(j),j=1,jx)',(cmo(j),j=1,jx)
      do 46 k=1,kx
         write(29,*)((u(i,j,k),i=1,ix),j=2,2)
 46   continue
      write(6,*)
      ax(1)=0.0
      ay(1)=0.0
      ax(2)=xi(ile,jstr)
      ay(2)=y(jstr)
      ax(3)=xi(ile,jtip)
      ay(3)=y(jtip)
      ax(4)=xi(ite,jtip)
      ay(4)=y(jtip)
      ax(5)=xi(ite,jstr)
      ay(5)=y(jstr)
      ax(6)=1.0
      ay(6)=0.0
      ax(7)=x(1)
      ay(7)=0.0
      do 47 i=1,7
         write(31,*)ay(i),ax(i)
 47   continue
c*****input files
      write(6,*)'********input files:'
      write(6,*)'tsd.data    : parameters'
c*****output files
      write(6,*)'*******output files:'
      write(6,*)'oneram6.xz  : inprof.ne.0 profile points'
      write(6,*)'tsd.itr     : residual vs. iteration'
      write(6,*)'tsd.ix      : mesh distribution in x'
      write(6,*)'tsd.kz      : mesh distribution in z'
      write(6,*)'tsd.cpo     : lower pressure coefficient'
      write(6,*)'tsd.cpu     : upper pressure coefficient'
      write(6,*)'tsd.gp      : vorticity distribution'
      write(6,*)'tsd.cpwo    : lower boundary Cp coefficient'
      write(6,*)'tsd.cpwu    : upper boundary Cp coefficient'
      write(6,*)'tsd.in      : potential from calculation'
      write(6,*)'tsd.out     : potential for restart'
      write(6,*)'tsd.cpcr    : critical Cp value'
      write(6,*)'tsd.itcl    : cl vs. iteration'
      write(6,*)'tsd.cpcon   : ix,jx,x(i),z(k),cp(i,k) for contour plot'
      write(6,*)'tsd.xzmses  : profile geopmetry'
      write(6,*)'tsd.dom     : projected planform of domain and wing'
 1000 format(1x,'iter=',i10,'   resx=',e17.8,'   idx=',i3,'   jdx=',i3
     &     ,'   kdx=',i3,'   cl=',f8.4,'   cdw=',f8.4)
 1001 format((1x,10f10.4))
 1002 format(' z(k)')
 1003 format(1x,' i=    ','    x(i)=    ','    d(i)=    '
     &     ,'    dp(i)=   ','    e(i)=    ','    ep(i)=')
 1004    format(1x,i4,5f13.4)
 1005 format(' ph(i,k)')
      end

      subroutine jjscheme(ixx,jxx,kxx,i,j,k,um,ui,xi,y,z,ph,aa,bb,cc,dd)
      implicit none
      integer ixx,jxx,kxx,i,j,k,jtip,ile,ite
      real pi,eps,gamp,mach0,ucr,bet0,gamach,um,ui,swp,piv,cdw,ddkm
      real xi(ixx,jxx),y(jxx),z(kxx)
      real ph(ixx,jxx,kxx)
      real aa(kxx),bb(kxx),cc(kxx),dd(kxx)
      common/constants/pi,eps,gamp,mach0,ucr,bet0,gamach,piv,cdw
      if(um.gt.ucr+eps)then
         if(ui.gt.ucr+eps)then
c     supersonic point
            aa(k)=-1.0/(z(k)-z(k-1))
     &           *0.5*(xi(i+1,j)-xi(i-1,j))
            bb(k)=-(bet0-gamach*um)/(xi(i,j)-xi(i-1,j))
     &           /(xi(i,j)-xi(i-1,j))
     &           *0.25*(xi(i+1,j)-xi(i-1,j))*(z(k+1)-z(k-1))
     &           +(1.0/(y(j+1)-y(j))+1.0/(y(j)-y(j-1)))
     &           /(y(j+1)-y(j-1))
     &           *0.5*(xi(i+1,j)-xi(i-1,j))*(z(k+1)-z(k-1))
     &           +(1.0/(z(k+1)-z(k))+1.0/(z(k)-z(k-1)))
     &           *0.5*(xi(i+1,j)-xi(i-1,j))
     &           +piv/(xi(i,j)-xi(i-1,j))
     &           *0.25*(xi(i+1,j)-xi(i-1,j))*(z(k+1)-z(k-1))
            cc(k)=-1.0/(z(k+1)-z(k))
     &           *0.5*(xi(i+1,j)-xi(i-1,j))
            ddkm=dd(k)
            dd(k)=0.0
            if(i.gt.2)then
               bb(k)=bb(k)-(bet0-gamach*um)/(xi(i,j)-xi(i-1,j))
     &              /(xi(i,j)-xi(i-2,j))
     &              *0.5*(xi(i+1,j)-xi(i-1,j))*(z(k+1)-z(k-1))
     (              +(bet0-gamach*um)/(xi(i,j)-xi(i-1,j))
     &              /(xi(i,j)-xi(i-1,j))
     &              *0.25*(xi(i+1,j)-xi(i-1,j))*(z(k+1)-z(k-1))
               dd(k)=(bet0-gamach*um)*((ph(i,j,k)-ph(i-1,j,k))
     &              /(xi(i,j)-xi(i-1,j))
     &              -(ph(i-1,j,k)-ph(i-2,j,k))/(xi(i-1,j)-xi(i-2,j)))
     &              /(xi(i,j)-xi(i-2,j))
     &              *0.5*(xi(i+1,j)-xi(i-1,j))*(z(k+1)-z(k-1))
            endif
            dd(k)=dd(k)
     &           +((ph(i,j+1,k)-ph(i,j,k))/(y(j+1)-y(j))
     &           -(ph(i,j,k)-ph(i,j-1,k))/(y(j)-y(j-1)))
     &           /(y(j+1)-y(j-1))
     &           *0.5*(xi(i,j)-xi(i-1,j))*(z(k+1)-z(k-1))
     &           +((ph(i,j,k+1)-ph(i,j,k))/(z(k+1)-z(k))
     &           -(ph(i,j,k)-ph(i,j,k-1))/(z(k)-z(k-1)))
     &           *0.5*(xi(i+1,j)-xi(i-1,j))
     &           -((xi(i,j+1)-xi(i,j))/(y(j+1)-y(j))
     &           +(xi(i,j)-xi(i,j-1))/(y(j)-y(j-1)))
     &           *(ph(i+1,j+1,k)-ph(i+1,j-1,k)
     &           -ph(i-1,j+1,k)+ph(i-1,j-1,k))
     &           /(y(j+1)-y(j-1))
     &           /(xi(i+1,j)-xi(i-1,j)-(xi(i+1,j+1)-xi(i+1,j-1)
     &           -xi(i-1,j+1)+xi(i-1,j-1))*y(j)/(y(j+1)-y(j-1)))
     &           *0.25*(xi(i+1,j)-xi(i-1,j))*(z(k+1)-z(k-1))
     &           -2.0*((xi(i,j+1)-xi(i,j))/(y(j+1)-y(j))
     &           -(xi(i,j)-xi(i,j-1))/(y(j)-y(j-1)))
     &           /(y(j+1)-y(j-1))*(ph(i+1,j,k)-ph(i-1,j,k))
     &           /(xi(i+1,j)-xi(i-1,j)-(xi(i+1,j+1)-xi(i+1,j-1)
     &           -xi(i-1,j+1)+xi(i-1,j-1))*y(j)/(y(j+1)-y(j-1)))
     &           *0.25*(xi(i+1,j)-xi(i-1,j))*(z(k+1)-z(k-1))
     &           +piv*ddkm/(xi(i,j)-xi(i-1,j))
     &           *0.25*(xi(i+1,j)-xi(i-1,j))*(z(k+1)-z(k-1))
c            write(6,*)'case 1'
         else
c     shock point
            aa(k)=-1.0/(z(k)-z(k-1))
     &           *0.5*(xi(i+1,j)-xi(i-1,j))
            bb(k)=-(bet0-gamach*um)/(xi(i,j)-xi(i-1,j))
     &           /(xi(i,j)-xi(i-1,j))
     &           *0.25*(xi(i+1,j)-xi(i-1,j))*(z(k+1)-z(k-1))
     &           +(1.0/(y(j+1)-y(j))+1.0/(y(j)-y(j-1)))
     &           /(y(j+1)-y(j-1))
     &           *0.5*(xi(i+1,j)-xi(i-1,j))*(z(k+1)-z(k-1))
     &           +(1.0/(z(k+1)-z(k))+1.0/(z(k)-z(k-1)))
     &           *0.5*(xi(i+1,j)-xi(i-1,j))
     &           +piv/(xi(i,j)-xi(i-1,j))
     &           *0.25*(xi(i+1,j)-xi(i-1,j))*(z(k+1)-z(k-1))
            cc(k)=-1.0/(z(k+1)-z(k))
     &           *0.5*(xi(i+1,j)-xi(i-1,j))
            ddkm=dd(k)
            dd(k)=0.0
            if(i.gt.2)then
               bb(k)=bb(k)-(bet0-gamach*um)/(xi(i,j)-xi(i-1,j))
     &              /(xi(i,j)-xi(i-2,j))
     &              *0.5*(xi(i+1,j)-xi(i-1,j))*(z(k+1)-z(k-1))
     (              +(bet0-gamach*um)/(xi(i,j)-xi(i-1,j))
     &              /(xi(i,j)-xi(i-1,j))
     &              *0.25*(xi(i+1,j)-xi(i-1,j))*(z(k+1)-z(k-1))
               dd(k)=(bet0-gamach*um)*((ph(i,j,k)-ph(i-1,j,k))
     &              /(xi(i,j)-xi(i-1,j))
     &              -(ph(i-1,j,k)-ph(i-2,j,k))/(xi(i-1,j)-xi(i-2,j)))
     &              /(xi(i,j)-xi(i-2,j))
     &              *0.5*(xi(i+1,j)-xi(i-1,j))*(z(k+1)-z(k-1))
            endif
            dd(k)=dd(k)
     &           +(bet0-gamach*ui)*((ph(i+1,j,k)-ph(i,j,k))
     &           /(xi(i+1,j)-xi(i,j))-(ph(i,j,k)-ph(i-1,j,k))
     &           /(xi(i,j)-xi(i-1,j)))*0.5*(z(k+1)-z(k-1))
     &           +((ph(i,j+1,k)-ph(i,j,k))/(y(j+1)-y(j))
     &           -(ph(i,j,k)-ph(i,j-1,k))/(y(j)-y(j-1)))
     &           /(y(j+1)-y(j-1))
     &           *0.5*(xi(i+1,j)-xi(i-1,j))*(z(k+1)-z(k-1))
     &           +((ph(i,j,k+1)-ph(i,j,k))/(z(k+1)-z(k))
     &           -(ph(i,j,k)-ph(i,j,k-1))/(z(k)-z(k-1)))
     &           *0.5*(xi(i+1,j)-xi(i-1,j))
     &           -((xi(i,j+1)-xi(i,j))/(y(j+1)-y(j))
     &           +(xi(i,j)-xi(i,j-1))/(y(j)-y(j-1)))
     &           *(ph(i+1,j+1,k)-ph(i+1,j-1,k)
     &           -ph(i-1,j+1,k)+ph(i-1,j-1,k))
     &           /(y(j+1)-y(j-1))
     &           /(xi(i+1,j)-xi(i-1,j)-(xi(i+1,j+1)-xi(i+1,j-1)
     &           -xi(i-1,j+1)+xi(i-1,j-1))*y(j)/(y(j+1)-y(j-1)))
     &           *0.25*(xi(i+1,j)-xi(i-1,j))*(z(k+1)-z(k-1))
     &           -2.0*((xi(i,j+1)-xi(i,j))/(y(j+1)-y(j))
     &           -(xi(i,j)-xi(i,j-1))/(y(j)-y(j-1)))
     &           /(y(j+1)-y(j-1))*(ph(i+1,j,k)-ph(i-1,j,k))
     &           /(xi(i+1,j)-xi(i-1,j)-(xi(i+1,j+1)-xi(i+1,j-1)
     &           -xi(i-1,j+1)+xi(i-1,j-1))*y(j)/(y(j+1)-y(j-1)))
     &           *0.25*(xi(i+1,j)-xi(i-1,j))*(z(k+1)-z(k-1))
     &           +piv*ddkm/(xi(i,j)-xi(i-1,j))
     &           *0.25*(xi(i+1,j)-xi(i-1,j))*(z(k+1)-z(k-1))
            cdw=cdw+2.0*(ucr-um)**3*(z(k+1)-z(k-1))*(y(j+1)-y(j-1))
c            write(6,*)'case 2'
         endif
      else
         if(ui.gt.ucr+eps)then
c     sonic point
            aa(k)=-1.0/(z(k)-z(k-1))
     &           *0.5*(xi(i+1,j)-xi(i-1,j))
            bb(k)=gamach*(ui-um)/(xi(i,j)-xi(i-1,j))
     &           /(xi(i,j)-xi(i-1,j))
     &           *0.5*(xi(i+1,j)-xi(i-1,j))*(z(k+1)-z(k-1))
     &           +(1.0/(y(j+1)-y(j))+1.0/(y(j)-y(j-1)))
     &           /(y(j+1)-y(j-1))
     &           *0.5*(xi(i+1,j)-xi(i-1,j))*(z(k+1)-z(k-1))
     &           +(1.0/(z(k+1)-z(k))+1.0/(z(k)-z(k-1)))
     &           *0.5*(xi(i+1,j)-xi(i-1,j))
     &           +piv/(xi(i,j)-xi(i-1,j))
     &           *0.25*(xi(i+1,j)-xi(i-1,j))*(z(k+1)-z(k-1))
            cc(k)=-1.0/(z(k+1)-z(k))
     &           *0.5*(xi(i+1,j)-xi(i-1,j))
            ddkm=dd(k)
            dd(k)=(bet0-gamach*(ph(i,j,k)-ph(i-1,j,k))
     &           /(xi(i,j)-xi(i-1,j)))
     &           *(ui-um)/(xi(i,j)-xi(i-1,j))
     &           *0.25*(xi(i+1,j)-xi(i-1,j))*(z(k+1)-z(k-1))
     &           +((ph(i,j+1,k)-ph(i,j,k))/(y(j+1)-y(j))
     &           -(ph(i,j,k)-ph(i,j-1,k))/(y(j)-y(j-1)))
     &           /(y(j+1)-y(j-1))
     &           *0.5*(xi(i+1,j)-xi(i-1,j))*(z(k+1)-z(k-1))
     &           +((ph(i,j,k+1)-ph(i,j,k))/(z(k+1)-z(k))
     &           -(ph(i,j,k)-ph(i,j,k-1))/(z(k)-z(k-1)))
     &           *0.5*(xi(i+1,j)-xi(i-1,j))
     &           -((xi(i,j+1)-xi(i,j))/(y(j+1)-y(j))
     &           +(xi(i,j)-xi(i,j-1))/(y(j)-y(j-1)))
     &           *(ph(i+1,j+1,k)-ph(i+1,j-1,k)
     &           -ph(i-1,j+1,k)+ph(i-1,j-1,k))
     &           /(y(j+1)-y(j-1))
     &           /(xi(i+1,j)-xi(i-1,j)-(xi(i+1,j+1)-xi(i+1,j-1)
     &           -xi(i-1,j+1)+xi(i-1,j-1))*y(j)/(y(j+1)-y(j-1)))
     &           *0.25*(xi(i+1,j)-xi(i-1,j))*(z(k+1)-z(k-1))
     &           -2.0*((xi(i,j+1)-xi(i,j))/(y(j+1)-y(j))
     &           -(xi(i,j)-xi(i,j-1))/(y(j)-y(j-1)))
     &           /(y(j+1)-y(j-1))*(ph(i+1,j,k)-ph(i-1,j,k))
     &           /(xi(i+1,j)-xi(i-1,j)-(xi(i+1,j+1)-xi(i+1,j-1)
     &           -xi(i-1,j+1)+xi(i-1,j-1))*y(j)/(y(j+1)-y(j-1)))
     &           *0.25*(xi(i+1,j)-xi(i-1,j))*(z(k+1)-z(k-1))
     &           +piv*ddkm/(xi(i,j)-xi(i-1,j))
     &           *0.25*(xi(i+1,j)-xi(i-1,j))*(z(k+1)-z(k-1))
c            write(6,*)'case 3'
         else
c     subsonic point
c            write(6,*)'gamach = ',gamach,' bet0 = ',bet0,' ui = ',ui
            aa(k)=-1.0/(z(k)-z(k-1))
     &           *0.5*(xi(i+1,j)-xi(i-1,j))
            bb(k)=(bet0-gamach*ui)
     &           *(1.0/(xi(i+1,j)-xi(i,j))+1.0/(xi(i,j)-xi(i-1,j)))
     &           *0.5*(z(k+1)-z(k-1))
     &           +(1.0/(y(j+1)-y(j))+1.0/(y(j)-y(j-1)))
     &           /(y(j+1)-y(j-1))
     &           *0.25*(xi(i+1,j)-xi(i-1,j))*(z(k+1)-z(k-1))
     &           +(1.0/(z(k+1)-z(k))+1.0/(z(k)-z(k-1)))
     &           *0.5*(xi(i+1,j)-xi(i-1,j))
     &           +piv/(xi(i,j)-xi(i-1,j))
     &           *0.25*(xi(i+1,j)-xi(i-1,j))*(z(k+1)-z(k-1))
            cc(k)=-1.0/(z(k+1)-z(k))
     &           *0.5*(xi(i+1,j)-xi(i-1,j))
            ddkm=dd(k)
            dd(k)=(bet0-gamach*ui)*((ph(i+1,j,k)-ph(i,j,k))
     &           /(xi(i+1,j)-xi(i,j))
     &           -(ph(i,j,k)-ph(i-1,j,k))/(xi(i,j)-xi(i-1,j)))
     &           *0.5*(z(k+1)-z(k-1))
     &           +((ph(i,j+1,k)-ph(i,j,k))/(y(j+1)-y(j))
     &           -(ph(i,j,k)-ph(i,j-1,k))/(y(j)-y(j-1)))
     &           /(y(j+1)-y(j-1))
     &           *0.5*(xi(i+1,j)-xi(i-1,j))*(z(k+1)-z(k-1))
     &           +((ph(i,j,k+1)-ph(i,j,k))/(z(k+1)-z(k))
     &           -(ph(i,j,k)-ph(i,j,k-1))/(z(k)-z(k-1)))
     &           *0.5*(xi(i+1,j)-xi(i-1,j))
     &           -((xi(i,j+1)-xi(i,j))/(y(j+1)-y(j))
     &           +(xi(i,j)-xi(i,j-1))/(y(j)-y(j-1)))
     &           *(ph(i+1,j+1,k)-ph(i+1,j-1,k)
     &           -ph(i-1,j+1,k)+ph(i-1,j-1,k))
     &           /(y(j+1)-y(j-1))
     &           /(xi(i+1,j)-xi(i-1,j)-(xi(i+1,j+1)-xi(i+1,j-1)
     &           -xi(i-1,j+1)+xi(i-1,j-1))*y(j)/(y(j+1)-y(j-1)))
     &           *0.25*(xi(i+1,j)-xi(i-1,j))*(z(k+1)-z(k-1))
     &           -2.0*((xi(i,j+1)-xi(i,j))/(y(j+1)-y(j))
     &           -(xi(i,j)-xi(i,j-1))/(y(j)-y(j-1)))
     &           /(y(j+1)-y(j-1))*(ph(i+1,j,k)-ph(i-1,j,k))
     &           /(xi(i+1,j)-xi(i-1,j)-(xi(i+1,j+1)-xi(i+1,j-1)
     &           -xi(i-1,j+1)+xi(i-1,j-1))*y(j)/(y(j+1)-y(j-1)))
     &           *0.25*(xi(i+1,j)-xi(i-1,j))*(z(k+1)-z(k-1))

c     &           +piv*ddkm/(xi(i,j)-xi(i-1,j))
c     &           *0.25*(xi(i+1,j)-xi(i-1,j))*(z(k+1)-z(k-1))
c            write(6,*)'case 4 bb_k', bb(k),'dd_k',dd(k)
         endif
      endif
c      write(6,*)'aa(k) = ',aa(k),' bb(k) = ',bb(k),' cc(k) = ',cc(k)
      return
      end
       

      subroutine tridiag(aa,bb,cc,ff,n1,n)
      implicit none
      integer n1,n,n2,n1n,k,k1
      real aa(n),bb(n),cc(n),ff(n)
      bb(n1)=1./bb(n1)
      aa(n1)=ff(n1)*bb(n1)
      n2=n1+1
      n1n=n1+n
      do k=n2,n
         k1=k-1
         cc(k1)=cc(k1)*bb(k1)
         bb(k)=bb(k)-aa(k)*cc(k1)
         bb(k)=1./bb(k)
         aa(k)=(ff(k)-aa(k)*aa(k1))*bb(k)
c         write(6,*)'aaT(k)',aa(k)
c         write(6,*)'bbT(k)',bb(k)
c         write(6,*)'ccT(k)',cc(k)
c         write(6,*)'ddT(k)',ff(k)
c         write(6,*)''
      enddo
c     back substitution
      ff(n)=aa(n)
      do k1=n2,n
         k=n1n-k1
         ff(k)=aa(k)-cc(k)*ff(k+1)
c         write(6,*)'k',k,'aaT(k)',aa(k)
c         write(6,*)'k',k,'bbT(k)',bb(k)
c         write(6,*)'k',k,'ccT(k)',cc(k)
c         write(6,*)'k',k,'ddT(k)',ff(k)
c         write(6,*)''
      enddo
      end

      subroutine camberdis(cxm,dm,dmplus,xi,fim)
      implicit none
      real cxm,dm,xi,fim,dmplus
c     parabolic relative camber dm>0 (alfades=0, Cldes=4*pi*dm, Cmac=-pi*dm)
      fim=4.*dm*xi*(1.-xi/cxm)
      if(dm.ge.0)return
c     Profile for tailless airplane cubic+parabolic dm<0 (alfades=-dm/3.0
c     Cldes=pi*(-dm+4.0*dmplus), Cmac=pi*dmplus)
      fim=-dm/3.*xi*(7.0-8.0*xi/cxm)*(1.-xi/cxm)
     &     +4.*dmplus*xi*(1.-xi/cxm)
      return
      end

      subroutine thickdis(ithick,cxm,em,xi,eim)
      implicit none
      integer ithick
      real cxm,em,xi,eim,fasc,faqj,fasl,fass,teti
      data fasc/.4592793/faqj/.3849002/fasl/.3088162/fass/.2414953/
c
c     1=elliptic distribution
c     2=semi-cubic distribution
c     3=quasi-joukowski distribution
c     4=naca00em distribution
c     5=selig distribution
c     6=super-selig distribution
c     7=biconvex distribution
c
      teti=acos(1.0-2.0*xi)
      goto(1,2,3,4,5,6,7)ithick
 1    eim=.5*em*cxm*sin(teti)
      return
 2    eim=fasc*em*cxm*sqrt(1.+cos(teti))*sin(teti)
      return
 3    eim=faqj*em*cxm*(1.+cos(teti))*sin(teti)
      return
 4    eim=5.*em*cxm*(.2969*sqrt(xi/cxm)-.126*xi/cxm-
     &.3537*(xi/cxm)**2+.2843*(xi/cxm)**3-.1015*(xi/cxm)**4)
      return
 5    eim=fasl*em*cxm*sin(teti)*(1.+cos(teti))**1.5
      return
 6    eim=fass*em*cxm*sin(teti)*(1.+cos(teti))**2
      return
 7    eim=2.0*em*cxm*xi*(1.0-xi)
      return
      end



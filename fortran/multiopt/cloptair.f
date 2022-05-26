      program cloptair
      implicit none
      integer ix,ixx,ic,n
      parameter(ix=51,ixx=ix+ix,n=ixx*ixx)
      doubleprecision fact(4),ex(4)
      doubleprecision x1(ix),axi1(ix),azi1(ix)
      doubleprecision x2(ix),axi2(ix),azi2(ix)
      doubleprecision f1(ix),fp1(ix),s1(ix)
      doubleprecision f2(ix),fp2(ix),s2(ix)
      doubleprecision e1(ix),ep1(ix)
      doubleprecision e2(ix),ep2(ix)
      doubleprecision g1(ix),gp1(ix)
      doubleprecision g2(ix),gp2(ix)
      doubleprecision u1(ix),w1(ix)
      doubleprecision u2(ix),w2(ix)
      doubleprecision uru1(ix),wru1(ix),cpu1(ix),fu1(ix)
      doubleprecision uro1(ix),wro1(ix),cpo1(ix),fo1(ix)
      doubleprecision uru2(ix),wru2(ix),cpu2(ix),fu2(ix)
      doubleprecision uro2(ix),wro2(ix),cpo2(ix),fo2(ix)
      doubleprecision veu(ix),xxu(ix),vgradu(ix),thetau(ix),cfru(ix)
      doubleprecision veo(ix),xxo(ix),vgrado(ix),thetao(ix),cfro(ix)
      doubleprecision ama(ixx,ixx),rhs(ixx)
      doubleprecision XINL,XOUT,YBOT,YTOP
      integer ipvt(ixx)
      character*17 thick(4)
c*****constants
      doubleprecision eps, pi, degrad, dtet
c     semi-cubic thickness parameter
      doubleprecision exsc, fasc
c     quasi-joukowski thickness distribution
      doubleprecision exqj, faqj
c     selig thickness distribution
      doubleprecision exsg, fasg
c     superselig thickness distribution
      doubleprecision exss, fass
c*****parameters
c*****iteration parameters
      integer itx
      doubleprecision omega, omgeo, omak
c*****design cl
      doubleprecision cld, ak
c*****location of te1
      doubleprecision m1, xte1
c*****thickness of profile 1
      integer k1
      doubleprecision ex1, fact1, n1, em1
c*****location of le2
      doubleprecision m2, xle2
c*****thickness of profile 2
      integer k2
      doubleprecision ex2, fact2, n2, em2, tethick
c*****slot height
      doubleprecision  gap
c*****initialization
      doubleprecision axii1, axii2, eim1, eim2, axiim1, axiim2
      doubleprecision teti, xi1, xi2
      integer i, iter, nopt
      doubleprecision dm, dmplus, fim1, alphad, alpha, dum, alpheffd
     &	,alpheff
c*****iteration loop
      integer itd, it, idx
      doubleprecision dgx, sum, sum1, sum2
      integer j
      doubleprecision aij, bij, cij, dij
c*****boundary conditions
      integer info, k, jc, im
      doubleprecision dgx0, coefj, f1m, f2m, reduc, cl1, cl2, clt, dak
c*****compute velocity due to thickness
c*****end of iteration-begin post processing
c*****compute fp
      integer ip
c*****regularize velocity field due to thickness
      doubleprecision gpm1, gpm2, dx, dz, axi, gpe
c*****output on listing
c*****output for nxyplot
      doubleprecision f1e
cccccccccccccccccccccccccccccccccc
      data ama/n*0.0e0/rhs/ixx*0.0e0/w1/ix*0.0e0/w2/ix*0.0e0/
      data g1/ix*0.0e0/g2/ix*0.0e0/
      open(unit=15,file='cloptair.data',form='formatted')
      open(unit=16,file='cloptcp1.dat',form='formatted')
      open(unit=17,file='cloptcp2.dat',form='formatted')
      open(unit=18,file='cloptxf1.dat',form='formatted')
      open(unit=19,file='cloptxf2.dat',form='formatted')
      open(unit=20,file='cloptcf1.dat',form='formatted')
      open(unit=21,file='cloptcf2.dat',form='formatted')
      open(unit=22,file='cloptxs1.dat',form='formatted')
      open(unit=23,file='cloptxs2.dat',form='formatted')
      open(unit=24,file='cloptxfe.dat',form='formatted')
      open(unit=25,file='cloptcpe.dat',form='formatted')
      open(unit=26,file='cloptair.out',form='formatted')
      open(unit=27,file='cloptxz1.dat',form='formatted')
      open(unit=28,file='cloptxz2.dat',form='formatted')
      open(unit=29,file='cloptxg1.dat',form='formatted')
      open(unit=30,file='cloptxg2.dat',form='formatted')
      open(unit=31,file='blade.profile',form='formatted')
c*****variables
      write(6,*)
      write(6,*)'*dimensionless variables:'
      write(6,*)'                       X=C*x'
      write(6,*)'                       Z=C*z'
      write(6,*)'                       G=Vinf*C*g'
      write(6,*)'                       U=Vinf*(1+u)'
      write(6,*)'                       W=Vinf*w'
      write(6,*)'                       D=0.5*Rho*Vinf**2*C*cd'
      write(6,*)'                       L=0.5*Rho*Vinf**2*C*cl'
      write(6,*)'                      Re=Rho*Vinf*C/Amu'
      write(26,*)
      write(26,*)'*dimensionless variables:'
      write(26,*)'                       X=C*x'
      write(26,*)'                       Z=C*z'
      write(26,*)'                       G=Vinf*C*g'
      write(26,*)'                       U=Vinf*(1+u)'
      write(26,*)'                       W=Vinf*w'
      write(26,*)'                       D=0.5*Rho*Vinf**2*C*cd'
      write(26,*)'                       L=0.5*Rho*Vinf**2*C*cl'
      write(26,*)'                      Re=Rho*Vinf*C/Amu'
c*****constants
      eps=1.e-10
      pi=2.*dasin(1.0d0)
      degrad=pi/180.
      dtet=pi/(ix-1)
c     semi-cubic thickness parameter
      exsc=0.5
      fasc=.125*3.*dsqrt(3.0d0/2.0d0)
c     quasi-joukowski thickness distribution
      exqj=1.0
      faqj=2./(3.*dsqrt(3.0d0))
c     selig thickness distribution
      exsg=1.5
      fasg=25.*dsqrt(2.5d0)/128.
c     superselig thickness distribution
      exss=2.
      fass=27./(50.*dsqrt(5.0d0))
c*****parameters
      ex(1)=exsc
      ex(2)=exqj
      ex(3)=exsg
      ex(4)=exss
      fact(1)=fasc
      fact(2)=faqj
      fact(3)=fasg
      fact(4)=fass
      thick(1)='  semi-cubic'
      thick(2)='  quasi-joukowski'
      thick(3)='  selig'
      thick(4)='  super-selig'
c*****iteration parameters
      read(15,*)itx
      read(15,*)omega
      read(15,*)omgeo
      read(15,*)omak
c*****reading geometric coefficients
c*****design cl
      read(15,*)cld
      ak=2.*cld/pi
c*****location of te1
      read(15,*)m1
      xte1=.01*m1
c*****thickness of profile 1
      read(15,*)k1
      ex1=ex(k1)
      fact1=fact(k1)
      read(15,*)n1
      em1=.01*n1*xte1
c*****location of le2
      read(15,*)m2
      xle2=.01*m2
c*****thickness of profile 2
      read(15,*)k2
      ex2=ex(k2)
      fact2=fact(k2)
      read(15,*)n2
      em2=.01*n2*(1.-xle2)
c*****trailing edge thickness
      read(15,*)tethick
c*****slot height
      read(15,*)gap
c*****flywingplus profile
      read(15,*)dm
      read(15,*)dmplus
      read(15,*)nopt
      read(15,*)alphad
      alpha=pi*alphad/180.
      alpheff=gap+xte1*alpha
      alpheffd=180.*alpheff/pi
c*****boundary domain
      read(15,*)XINL
      read(15,*)XOUT
      read(15,*)YBOT
      read(15,*)YTOP
c*****output on listing
      write(6,*)
      write(6,*)'***************'
      write(6,*)'iteration parameters:'
      write(6,*)'      maximum iterations itx=',itx
      write(6,*)'                       omega=',omega
      write(6,*)'     relaxation factor omgeo=',omgeo
      write(6,*)'      relaxation factor omak=',omak
      write(6,*)'***************'
      write(6,*)'profile data:'
      write(6,*)'     profile 1 trailing edge=',xte1
      write(6,*)'      thickness distribution=',thick(k1)
      write(6,*)'             thickness ratio=',n1,' (%)'
      write(6,*)'     profile 2  leading edge=',xle2
      write(6,*)'      thickness distribution=',thick(k2)
      write(6,*)'             thickness ratio=',n2,' (%)'
      write(6,*)' 1/2 trailing edge thickness=',tethick
      write(6,*)'                 slot height=',gap
      write(6,*)'                   design cl=',cld
      write(6,*)'           flywing camber dm=',dm
      write(6,*)'reversed parab camber dmplus=',dmplus
      write(6,*)'opt1 element/2 elements nopt=',nopt
      write(6,*)'  incidence of elem 1 alphad=',alphad
     &     ,' (deg)   alpha=',alpha,' (rd)'
      write(6,*)'effective incidence alpheffd=',alpheffd
     &     ,' (deg) alpheff=',alpheff,' (rd)'
      write(6,*)'inlet boundary location XINL=',XINL
      write(6,*)'        outlet location XOUT=',XOUT
      write(6,*)'        bottom location YBOT=',YBOT
      write(6,*)'  top boundary location YTOP=',YTOP
      write(26,*)
      write(26,*)'***************'
      write(26,*)'iteration parameters:'
      write(26,*)'      maximum iterations itx=',itx
      write(26,*)'                       omega=',omega
      write(26,*)'     relaxation factor omgeo=',omgeo
      write(26,*)'      relaxation factor omak=',omak
      write(26,*)'***************'
      write(26,*)'profile data:'
      write(26,*)'     profile 1 trailing edge=',xte1
      write(26,*)'      thickness distribution=',thick(k1)
      write(26,*)'             thickness ratio=',n1,' (%)'
      write(26,*)'     profile 2  leading edge=',xle2
      write(26,*)'      thickness distribution=',thick(k2)
      write(26,*)'             thickness ratio=',n2,' (%)'
      write(26,*)' 1/2 trailing edge thickness=',tethick
      write(26,*)'                 slot height=',gap
      write(26,*)'                   design cl=',cld
      write(26,*)'           flywing camber dm=',dm
      write(26,*)'reversed parab camber dmplus=',dmplus
      write(26,*)'opt1 element/2 elements nopt=',nopt
      write(26,*)'  incidence of elem 1 alphad=',alphad
     &     ,' (deg)   alpha=',alpha,' (rd)'
      write(26,*)'effective incidence alpheffd=',alpheffd
     &     ,' (deg) alpheff=',alpheff,' (rd)'
      write(26,*)'inlet boundary location XINL=',XINL
      write(26,*)'        outlet location XOUT=',XOUT
      write(26,*)'        bottom location YBOT=',YBOT
      write(26,*)'  top boundary location YTOP=',YTOP
c*****initialization
      write(6,*)' '
      write(6,*)'  i ','    x1(i) ','    f1(i) ','   fp1(i) '
     &     ,'    e1(i) ','   ep1(i) ','    x2(i) ','    f2(i) '
     &     ,'   fp2(i) ','    e2(i) ','   ep2(i)'
      write(26,*)'  i ','    x1(i) ','    f1(i) ','   fp1(i) '
     &     ,'    e1(i) ','   ep1(i) ','    x2(i) ','    f2(i) '
     &     ,'   fp2(i) ','    e2(i) ','   ep2(i)'
      write(26,*)' '
      write(26,*)'  i ','    x1(i) ','    f1(i) ','    e1(i) '
     &     ,'    x2(i) ','    f2(i) ','    e2(i) '
      axii1=-.5*xte1*(1.-dcos(.5*dtet))
      axii2=xle2-.5*(1.-xle2)*(1.-dcos(.5*dtet))
      fim1=gap
c     if element 1 fixed
      if(nopt.eq.1)then
         fim1=gap+dm*axii1*(7.-8.*axii1/xte1)*(1.-axii1/xte1)
     &        +dmplus*axii1*(1.-axii1/xte1)
         fim1=fim1+xte1*(1.-axii1/xte1)*alpha
      endif
      eim1=-fact1*em1*(1.+dcos(.5*dtet))**ex1*dsin(.5*dtet)
      eim2=-fact2*em2*(1.+dcos(.5*dtet))**ex2*dsin(.5*dtet)
      axiim1=axii1
      axiim2=axii2
      do 1 i=1,ix
         teti=(i-1)*dtet
         xi1=.5*xte1*(1.-dcos(teti))
         axii1=.5*xte1*(1.-dcos(teti+.5*dtet))
         x1(i)=xi1
         axi1(i)=axii1
         xi2=xle2+.5*(1.-xle2)*(1.-dcos(teti))
         axii2=xle2+.5*(1.-xle2)*(1.-dcos(teti+.5*dtet))
         x2(i)=xi2
         axi2(i)=axii2
         f1(i)=gap
         f2(i)=0.
         fp1(i)=0.
         fp2(i)=0.
c     if element 1 fixed
         if(nopt.eq.1)then
            f1(i)=gap+dm*axii1*(7.-8.*axii1/xte1)*(1.-axii1/xte1)/3.
     &           +dmplus*axii1*(1.-axii1/xte1)
            f1(i)=f1(i)+xte1*(1.-axii1/xte1)*alpha
            azi1(i)=f1(i)
            fp1(i)=(f1(i)-fim1)/(axii1-axiim1)
            fim1=f1(i)
            f1(i)=gap+dm*xi1*(7.-8.*xi1/xte1)*(1.-xi1/xte1)/3.
     &           +dmplus*xi1*(1.-xi1/xte1)
            f1(i)=f1(i)+xte1*(1.-xi1/xte1)*alpha
         endif
         e1(i)=fact1*em1*(1.+dcos(teti+.5*dtet))**ex1*dsin(teti+.5*dtet)
         e2(i)=fact2*em2*(1.+dcos(teti+.5*dtet))**ex2*dsin(teti+.5*dtet)
         ep1(i)=(e1(i)-eim1)/(axii1-axiim1)
         ep2(i)=(e2(i)-eim2)/(axii2-axiim2)
         eim1=e1(i)
         eim2=e2(i)
         e1(i)=fact1*em1*(1.+dcos(teti))**ex1*dsin(teti)
         e2(i)=fact2*em2*(1.+dcos(teti))**ex2*dsin(teti)
c     if element 1 fixed
         if(nopt.eq.2)then
            azi1(i)=f1(i)
         endif
         azi2(i)=f2(i)
         axiim1=axii1
         axiim2=axii2
         if(i.eq.ix-1)then
            eim1=-3.*eim1
            axiim1=2.*xte1-axii1
            eim2=-3.*eim2
            axiim2=2.-axii2
         endif
         gp1(i)=0.
         g1(i)=0.
         gp2(i)=0.
         g2(i)=0.
         s1(i)=f1(i)
         s2(i)=f2(i)
         if(i.lt.ix)write(6,1000)i,x1(i),f1(i),fp1(i),e1(i),ep1(i)
     &        ,x2(i),f2(i),fp2(i),e2(i),ep2(i)
         if(i.lt.ix)write(26,1000)i,x1(i),f1(i),fp1(i),e1(i),ep1(i)
     &        ,x2(i),f2(i),fp2(i),e2(i),ep2(i)
 1    continue
      axi1(ix)=axi1(ix-1)
      axi2(ix)=axi2(ix-1)
      azi1(ix)=azi1(ix-1)
      azi2(ix)=azi2(ix-1)
      fp1(ix)=fp1(ix-1)
      fp2(ix)=fp2(ix-1)
      write(6,1000)ix,x1(ix),f1(ix),fp1(ix),e1(ix),ep1(ix)
     &     ,x2(ix),f2(ix),fp2(ix),e2(ix),ep2(ix)
      write(26,1000)ix,x1(ix),f1(ix),fp1(ix),e1(ix),ep1(ix)
     &     ,x2(ix),f2(ix),fp2(ix),e2(ix),ep2(ix)
      iter=0
c*****iteration loop
 100  continue
      write(6,*)
      write(6,*)' iterations?(Y=1/N=0)'
      write(26,*)
      write(26,*)' iterations?(Y=1/N=0)'
      read(5,*)itd 
      if(itd.le.0)goto 400
      do 200 it=1,itx
      iter=iter+1
      idx=0
      dgx=0.
      sum=0.
      do 3 i=2,ix
         sum1=0.
         sum2=0.
         do 2 j=2,ix
            aij=-((x1(i-1)-axi1(j-1))/
     &  ((x1(i-1)-axi1(j-1))**2+(f1(i-1)-azi1(j-1))**2)
     &            -(x1(i)-axi1(j-1))/
     &  ((x1(i)-axi1(j-1))**2+(f1(i)-azi1(j-1))**2))/(x1(i)-x1(i-1))
     &          -((x1(j-1)-axi1(i-1))/
     &  ((x1(j-1)-axi1(i-1))**2+(f1(j-1)-azi1(i-1))**2)
     &            -(x1(j)-axi1(i-1))/
     &  ((x1(j)-axi1(i-1))**2+(f1(j)-azi1(i-1))**2))/(x1(j)-x1(j-1))
            bij=-((x1(i-1)-axi2(j-1))/
     &  ((x1(i-1)-axi2(j-1))**2+(f1(i-1)-azi2(j-1))**2)
     &            -(x1(i)-axi2(j-1))/
     &  ((x1(i)-axi2(j-1))**2+(f1(i)-azi2(j-1))**2))/(x1(i)-x1(i-1))
     &          -((x2(j-1)-axi1(i-1))/
     &  ((x2(j-1)-axi1(i-1))**2+(f2(j-1)-azi1(i-1))**2)
     &            -(x2(j)-axi1(i-1))/
     &  ((x2(j)-axi1(i-1))**2+(f2(j)-azi1(i-1))**2))/(x2(j)-x2(j-1))
            cij=-((x2(i-1)-axi1(j-1))/
     &  ((x2(i-1)-axi1(j-1))**2+(f2(i-1)-azi1(j-1))**2)
     &            -(x2(i)-axi1(j-1))/
     &  ((x2(i)-axi1(j-1))**2+(f2(i)-azi1(j-1))**2))/(x2(i)-x2(i-1))
     &          -((x1(j-1)-axi2(i-1))/
     &  ((x1(j-1)-axi2(i-1))**2+(f1(j-1)-azi2(i-1))**2)
     &            -(x1(j)-axi2(i-1))/
     &  ((x1(j)-axi2(i-1))**2+(f1(j)-azi2(i-1))**2))/(x1(j)-x1(j-1))
            dij=-((x2(i-1)-axi2(j-1))/
     &  ((x2(i-1)-axi2(j-1))**2+(f2(i-1)-azi2(j-1))**2)
     &            -(x2(i)-axi2(j-1))/
     &  ((x2(i)-axi2(j-1))**2+(f2(i)-azi2(j-1))**2))/(x2(i)-x2(i-1))
     &          -((x2(j-1)-axi2(i-1))/
     &  ((x2(j-1)-axi2(i-1))**2+(f2(j-1)-azi2(i-1))**2)
     &            -(x2(j)-axi2(i-1))/
     &  ((x2(j)-axi2(i-1))**2+(f2(j)-azi2(i-1))**2))/(x2(j)-x2(j-1))
            if(i.lt.ix)then
               aij=aij-((x1(j)-axi1(i))/
     &  ((x1(j)-axi1(i))**2+(f1(j)-azi1(i))**2)
     &                  -(x1(j-1)-axi1(i))/
     &  ((x1(j-1)-axi1(i))**2+(f1(j-1)-azi1(i))**2))/(x1(j)-x1(j-1))
     &                -((x1(i+1)-axi1(j-1))/
     &  ((x1(i+1)-axi1(j-1))**2+(f1(i+1)-azi1(j-1))**2)
     &                  -(x1(i)-axi1(j-1))/
     &  ((x1(i)-axi1(j-1))**2+(f1(i)-azi1(j-1))**2))/(x1(i+1)-x1(i))
               bij=bij-((x2(j)-axi1(i))/
     &  ((x2(j)-axi1(i))**2+(f2(j)-azi1(i))**2)
     &                  -(x2(j-1)-axi1(i))/
     &  ((x2(j-1)-axi1(i))**2+(f2(j-1)-azi1(i))**2))/(x2(j)-x2(j-1))
     &                -((x1(i+1)-axi2(j-1))/
     &  ((x1(i+1)-axi2(j-1))**2+(f1(i+1)-azi2(j-1))**2)
     &                  -(x1(i)-axi2(j-1))/
     &  ((x1(i)-axi2(j-1))**2+(f1(i)-azi2(j-1))**2))/(x1(i+1)-x1(i))
               cij=cij-((x1(j)-axi2(i))/
     &  ((x1(j)-axi2(i))**2+(f1(j)-azi2(i))**2)
     &                  -(x1(j-1)-axi2(i))/
     &  ((x1(j-1)-axi2(i))**2+(f1(j-1)-azi2(i))**2))/(x1(j)-x1(j-1))
     &                -((x2(i+1)-axi1(j-1))/
     &  ((x2(i+1)-axi1(j-1))**2+(f2(i+1)-azi1(j-1))**2)
     &                  -(x2(i)-axi1(j-1))/
     &  ((x2(i)-axi1(j-1))**2+(f2(i)-azi1(j-1))**2))/(x2(i+1)-x2(i))
               dij=dij-((x2(j)-axi2(i))/
     &  ((x2(j)-axi2(i))**2+(f2(j)-azi2(i))**2)
     &                  -(x2(j-1)-axi2(i))/
     &  ((x2(j-1)-axi2(i))**2+(f2(j-1)-azi2(i))**2))/(x2(j)-x2(j-1))
     &                -((x2(i+1)-axi2(j-1))/
     &  ((x2(i+1)-axi2(j-1))**2+(f2(i+1)-azi2(j-1))**2)
     &                  -(x2(i)-axi2(j-1))/
     &  ((x2(i)-axi2(j-1))**2+(f2(i)-azi2(j-1))**2))/(x2(i+1)-x2(i))
            endif
            if(j.lt.ix)then
               aij=aij-((x1(i)-axi1(j))/
     &  ((x1(i)-axi1(j))**2+(f1(i)-azi1(j))**2)
     &                  -(x1(i-1)-axi1(j))/
     &  ((x1(i-1)-axi1(j))**2+(f1(i-1)-azi1(j))**2))/(x1(i)-x1(i-1))
     &                -((x1(j+1)-axi1(i-1))/
     &  ((x1(j+1)-axi1(i-1))**2+(f1(j+1)-azi1(i-1))**2)
     &                  -(x1(j)-axi1(i-1))/
     &  ((x1(j)-axi1(i-1))**2+(f1(j)-azi1(i-1))**2))/(x1(j+1)-x1(j))
               bij=bij-((x1(i)-axi2(j))/
     &  ((x1(i)-axi2(j))**2+(f1(i)-azi2(j))**2)
     &                  -(x1(i-1)-axi2(j))/
     &  ((x1(i-1)-axi2(j))**2+(f1(i-1)-azi2(j))**2))/(x1(i)-x1(i-1))
     &                -((x2(j+1)-axi1(i-1))/
     &  ((x2(j+1)-axi1(i-1))**2+(f2(j+1)-azi1(i-1))**2)
     &                  -(x2(j)-axi1(i-1))/
     &  ((x2(j)-axi1(i-1))**2+(f2(j)-azi1(i-1))**2))/(x2(j+1)-x2(j))
               cij=cij-((x2(i)-axi1(j))/
     &  ((x2(i)-axi1(j))**2+(f2(i)-azi1(j))**2)
     &                  -(x2(i-1)-axi1(j))/
     &  ((x2(i-1)-axi1(j))**2+(f2(i-1)-azi1(j))**2))/(x2(i)-x2(i-1))
     &                -((x1(j+1)-axi2(i-1))/
     &  ((x1(j+1)-axi2(i-1))**2+(f1(j+1)-azi2(i-1))**2)
     &                  -(x1(j)-axi2(i-1))/
     &  ((x1(j)-axi2(i-1))**2+(f1(j)-azi2(i-1))**2))/(x1(j+1)-x1(j))
               dij=dij-((x2(i)-axi2(j))/
     &  ((x2(i)-axi2(j))**2+(f2(i)-azi2(j))**2)
     &                  -(x2(i-1)-axi2(j))/
     &  ((x2(i-1)-axi2(j))**2+(f2(i-1)-azi2(j))**2))/(x2(i)-x2(i-1))
     &                -((x2(j+1)-axi2(i-1))/
     &  ((x2(j+1)-axi2(i-1))**2+(f2(j+1)-azi2(i-1))**2)
     &                  -(x2(j)-axi2(i-1))/
     &  ((x2(j)-axi2(i-1))**2+(f2(j)-azi2(i-1))**2))/(x2(j+1)-x2(j))
            endif
            if(i.lt.ix.and.j.lt.ix)then
               aij=aij-((x1(i)-axi1(j))/
     &  ((x1(i)-axi1(j))**2+(f1(i)-azi1(j))**2)
     &                  -(x1(i+1)-axi1(j))/
     &  ((x1(i+1)-axi1(j))**2+(f1(i+1)-azi1(j))**2))/(x1(i+1)-x1(i))
     &                -((x1(j)-axi1(i))/
     &  ((x1(j)-axi1(i))**2+(f1(j)-azi1(i))**2)
     &                  -(x1(j+1)-axi1(i))/
     &  ((x1(j+1)-axi1(i))**2+(f1(j+1)-azi1(i))**2))/(x1(j+1)-x1(j))
               bij=bij-((x1(i)-axi2(j))/
     &  ((x1(i)-axi2(j))**2+(f1(i)-azi2(j))**2)
     &                  -(x1(i+1)-axi2(j))/
     &  ((x1(i+1)-axi2(j))**2+(f1(i+1)-azi2(j))**2))/(x1(i+1)-x1(i))
     &                -((x2(j)-axi1(i))/
     &  ((x2(j)-axi1(i))**2+(f2(j)-azi1(i))**2)
     &                  -(x2(j+1)-axi1(i))/
     &  ((x2(j+1)-axi1(i))**2+(f2(j+1)-azi1(i))**2))/(x2(j+1)-x2(j))
               cij=cij-((x2(i)-axi1(j))/
     &  ((x2(i)-axi1(j))**2+(f2(i)-azi1(j))**2)
     &                  -(x2(i+1)-axi1(j))/
     &  ((x2(i+1)-axi1(j))**2+(f2(i+1)-azi1(j))**2))/(x2(i+1)-x2(i))
     &                -((x1(j)-axi2(i))/
     &  ((x1(j)-axi2(i))**2+(f1(j)-azi2(i))**2)
     &                  -(x1(j+1)-axi2(i))/
     &  ((x1(j+1)-axi2(i))**2+(f1(j+1)-azi2(i))**2))/(x1(j+1)-x1(j))
               dij=dij-((x2(i)-axi2(j))/
     &  ((x2(i)-axi2(j))**2+(f2(i)-azi2(j))**2)
     &                  -(x2(i+1)-axi2(j))/
     &  ((x2(i+1)-axi2(j))**2+(f2(i+1)-azi2(j))**2))/(x2(i+1)-x2(i))
     &                -((x2(j)-axi2(i))/
     &  ((x2(j)-axi2(i))**2+(f2(j)-azi2(i))**2)
     &                  -(x2(j+1)-axi2(i))/
     &  ((x2(j+1)-axi2(i))**2+(f2(j+1)-azi2(i))**2))/(x2(j+1)-x2(j))
            endif
            if(m2.ge.99)then
               bij=0.
               cij=0.
               dij=0.
               if(j.eq.i)dij=1.
            endif
            ama(i,j)=aij
c     if element 1 fixed
            if(nopt.eq.1)then
               if(i.eq.ix)bij=0.
            endif
            ama(i,ix+j)=bij
            ama(ix+i,j)=cij
            ama(ix+i,ix+j)=dij
            sum1=sum1+aij*g1(j)+bij*g2(j)
            sum2=sum2+cij*g1(j)+dij*g2(j)
 2       continue
c     if element 1 fixed
         if(nopt.eq.1)then
            sum1=0.
         endif
         sum=sum+sum1*g1(i)+sum2*g2(i)
         if(i.eq.ix)then
            sum1=sum1-4.*pi*ak
c     if element 1 fixed
            if(m2.lt.99)then
               sum2=sum2-4.*pi*ak
            endif
         endif
         if(nopt.eq.1)sum1=0.
         rhs(i)=-sum1
         rhs(ix+i)=-sum2
 3    continue
      sum=sum/(2.*pi)
c*****boundary conditions
      ama(1,1)=1.
      rhs(1)=0.
      ama(ix+1,ix+1)=1.
      rhs(ix+1)=0.
c     if element 1 fixed
      if(nopt.eq.1)then
         do 600 i=1,ix
            dum=0.
            do 500 j=1,ix-1
               aij=0.
               bij=0.
               if(i.gt.1.and.i.lt.ix)then
                  aij=(x1(i)-axi1(j))/
     &                 ((x1(i)-axi1(j))**2+(f1(i)-azi1(j))**2)
                  bij=(x1(i)-axi2(j))/
     &                 ((x1(i)-axi2(j))**2+(f1(i)-azi2(j))**2)
                  if(j.gt.1)then
                     aij=aij-(x1(i)-axi1(j-1))/
     &                    ((x1(i)-axi1(j-1))**2+(f1(i)-azi1(j-1))**2)
                     bij=bij-(x1(i)-axi2(j-1))/
     &                    ((x1(i)-axi2(j-1))**2+(f1(i)-azi2(j-1))**2)
                  endif
                  dum=dum+bij*g2(j)
               endif
               ama(i,j)=aij
               ama(i,ix+j)=0.*bij
 500        continue
            bij=-(x1(i)-axi2(ix-1))
     &           /((x1(i)-axi2(ix-1))**2+(f1(i)-azi2(ix-1))**2)
            dum=dum+bij*g2(i)
            ama(1,i)=0.
            ama(ix,i)=0.
            if(i.gt.1)ama(i,ix)=-(x1(i)-axi1(ix-1))
     &           /((x1(i)-axi1(ix-1))**2+(f1(i)-azi1(ix-1))**2)
            ama(ix,ix+i)=0.
            if(i.gt.1.and.i.lt.ix)then
                rhs(i)=2.*pi*fp1(i)-dum
                ama(i,ix+ix)=0.
            endif
 600     continue
         ama(1,1)=1.
         ama(ix,1)=0.
      endif
c     if element 1 fixed
      if(nopt.eq.1)then
         ama(ix,ix-2)=(x1(ix)-x1(ix-1))**2
     &        /((x1(ix)-x1(ix-2))**2-(x1(ix)-x1(ix-1))**2)
         ama(ix,ix-1)=-(x1(ix)-x1(ix-2))**2
     &        /((x1(ix)-x1(ix-2))**2-(x1(ix)-x1(ix-1))**2)
         ama(ix,ix)=1.
         rhs(ix)=0.
      endif
      call dgefa(ama,ixx,ixx,ipvt,info)
      if(info.ne.0)then
         write(6,*)' zero pivot, exiting'
         stop
      endif
      call dgesl(ama,ixx,ixx,ipvt,rhs,info)
      do 4 j=1,ix
         g1(j)=g1(j)+omega*rhs(j)
         g2(j)=g2(j)+omega*rhs(ix+j)
c     if element 1 fixed
         if(nopt.eq.1)g1(j)=rhs(j)
         if(m2.ge.99)then
            rhs(ix+j)=0.
         endif
         if(dabs(rhs(j)).gt.dabs(dgx))then
            dgx=rhs(j)
            if(it.eq.1)dgx0=dabs(dgx)
            idx=j
         endif
         if(dabs(rhs(ix+j)).gt.dabs(dgx))then
            dgx=rhs(ix+j)
            if(it.eq.1)dgx0=dabs(dgx)
            idx=j
         endif
         do 4 k=1,ix
            ama(j,k)=0.
            ama(j,ix+k)=0.
            ama(ix+j,k)=0.
            ama(ix+j,ix+k)=0.
 4    continue
      do 7 j=2,ix-1
         sum1=0.
         sum2=0.
         do 5 k=1,ix-1
            coefj=(x1(j)-axi1(k))/
     &           ((x1(j)-axi1(k))**2+(f1(j)-azi1(k))**2)
            sum1=sum1+(g1(k+1)-g1(k))*coefj
            coefj=(x2(j)-axi2(k))/
     &           ((x2(j)-axi2(k))**2+(f2(j)-azi2(k))**2)
            sum2=sum2+(g2(k+1)-g2(k))*coefj
 5       continue
         do 6 k=1,ix-1
            coefj=(x1(j)-axi2(k))/
     &           ((x1(j)-axi2(k))**2+(f1(j)-azi2(k))**2)
            sum1=sum1+(g2(k+1)-g2(k))*coefj
            coefj=(x2(j)-axi1(k))/
     &           ((x2(j)-axi1(k))**2+(f2(j)-azi1(k))**2)
            sum2=sum2+(g1(k+1)-g1(k))*coefj
 6       continue
         w1(j)=-sum1/(2.*pi)
         w2(j)=-sum2/(2.*pi)
 7    continue
      w1(1)=w1(2)
     &     +(w1(3)-w1(2))*(x1(1)-x1(2))/(x1(3)-x1(2))
      w1(ix)=w1(ix-1)
     &     +(w1(ix-2)-w1(ix-1))*(x1(ix)-x1(ix-1))/(x1(ix-2)-x1(ix-1))
      w2(1)=w2(2)
     &     +(w2(3)-w2(2))*(x2(1)-x2(2))/(x2(3)-x2(2))
      w2(ix)=w2(ix-1)
     &     +(w2(ix-2)-w2(ix-1))*(x2(ix)-x2(ix-1))/(x2(ix-2)-x2(ix-1))
      sum1=0.
      s1(ix)=f1(ix)
      sum2=0.
      s2(1)=f2(1)
      do 8 j=1,ix-1
         jc=ix-j
c     if element 1 fixed
         if(nopt.eq.2)then
            sum1=sum1-.5*(w1(jc)+w1(jc+1))*(x1(jc+1)-x1(jc))
            s1(jc)=s1(jc)+omgeo*(s1(ix)+sum1-s1(jc))
         endif
         sum2=sum2+.5*(w2(j)+w2(j+1))*(x2(j+1)-x2(j))
         s2(j+1)=s2(j+1)+omgeo*(s2(1)+sum2-s2(j+1))
 8    continue
      im=1
      if(m2.ge.99)s2(ix)=s1(ix)
      f1m=s1(1)-s2(ix)
      f2m=s2(1)-s2(ix)
      do 9 i=1,ix
c     if element 1 fixed
         if(nopt.eq.2)then
            s1(i)=s1(i)-s2(ix)
            f1(i)=s1(i)
         endif
         s2(i)=s2(i)-s2(ix)
         f2(i)=s2(i)
c     if element 1 fixed
         if(nopt.eq.2)then
            azi1(im)=.5*(f1m+f1(i))
         endif
         azi2(im)=.5*(f2m+f2(i))
         f1m=f1(i)
         f2m=f2(i)
         im=i
 9    continue
      cl1=2.*g1(ix)
      cl2=2.*g2(ix)
      clt=cl1+cl2
      reduc=dabs(dgx)/(dgx0+eps)
      if(abs(dgx).lt.eps)then
         write(6,*)'iter=',iter,' dgx=',dgx,' idx=',idx
     &        ,' sum=',sum,' reduc=',reduc
         write(26,*)'iter=',iter,' dgx=',dgx,' idx=',idx
     &        ,' sum=',sum,' reduc=',reduc
         goto 300
      endif
 200  continue
 300  continue
      dak=ak*cld/clt-ak
      ak=ak+omak*dak
      write(6,*)'it=',it,'ak=',ak,' clt=',clt
      write(26,*)'it=',it,'ak=',ak,' clt=',clt
      goto 100
 400  continue
c*****compute velocity due to thickness
      do 11 i=2,ix-1
         sum1=0.
         sum2=0.
         do 10 j=1,ix-1
            coefj=1./(x1(i)-axi1(j))
            sum1=sum1+(e1(j+1)-e1(j))*coefj
            coefj=(x1(i)-axi2(j))/
     &           ((x1(i)-axi2(j))**2+(f1(i)-azi2(j))**2)
            sum1=sum1+(e2(j+1)-e2(j))*coefj
            coefj=1./(x2(i)-axi2(j))
            sum2=sum2+(e2(j+1)-e2(j))*coefj
            coefj=(x2(i)-axi1(j))/
     &           ((x2(i)-axi1(j))**2+(f2(i)-azi1(j))**2)
            sum2=sum2+(e1(j+1)-e1(j))*coefj
 10      continue
         u1(i)=sum1/pi
         u2(i)=sum2/pi
         if(m2.ge.99)u2(i)=0.
 11   continue
c*****end of iteration-begin post processing
      u1(1)=u1(2)
     &     +(u1(3)-u1(2))*(x1(1)-x1(2))/(x1(3)-x1(2))
      u1(ix)=u1(ix-1)
     &     +(u1(ix-2)-u1(ix-1))*(x1(ix)-x1(ix-1))/(x1(ix-2)-x1(ix-1))
      u2(1)=u2(2)
     &     +(u2(3)-u2(2))*(x2(1)-x2(2))/(x2(3)-x2(2))
      u2(ix)=u2(ix-1)
     &     +(u2(ix-2)-u2(ix-1))*(x2(ix)-x2(ix-1))/(x2(ix-2)-x2(ix-1))
      do 12 i=1,ix-1
         gp1(i)=(g1(i+1)-g1(i))/(x1(i+1)-x1(i))
         gp2(i)=(g2(i+1)-g2(i))/(x2(i+1)-x2(i))
 12   continue
      gp1(ix)=gp1(ix-1)
      gp2(ix)=gp2(ix-1)
c*****compute fp
      do 13 i=1,ix
         im=i-1
         if(i.eq.1)im=1
         ip=i+1
         if(i.eq.ix)ip=ix
         fp1(i)=(f1(ip)-f1(im))/(x1(ip)-x1(im))
         fp2(i)=(f2(ip)-f2(im))/(x2(ip)-x2(im))
 13   continue
c*****regularize velocity field due to thickness
      axiim1=0.
      axiim2=xle2
      gpm1=0.
      gpm2=0.
      do 14 i=1,ix
         dx=axi1(i)-axiim1
         dz=(fp1(i)+ep1(i))*dx
         uru1(i)=(1.+u1(i)+.25*(gpm1+gp1(i)))*dx*dx/(dx*dx+dz*dz)
         wru1(i)=uru1(i)*dz/dx
         uru1(i)=uru1(i)-1.
         dz=(fp1(i)-ep1(i))*dx
         uro1(i)=(1.+u1(i)-.25*(gpm1+gp1(i)))*dx*dx/(dx*dx+dz*dz)
         wro1(i)=uro1(i)*dz/dx
         uro1(i)=uro1(i)-1.
         fu1(i)=f1(i)+e1(i)
         fo1(i)=f1(i)-e1(i)
         axiim1=axi1(i)
         gpm1=gp1(i)
         if(i.eq.ix-1)then
            axiim1=-xte1+2.*axiim1
         endif
         dx=axi2(i)-axiim2
         dz=(fp2(i)+ep2(i))*dx
         uru2(i)=(1.+u2(i)+.25*(gpm2+gp2(i)))*dx*dx/(dx*dx+dz*dz)
         wru2(i)=uru2(i)*dz/dx
         uru2(i)=uru2(i)-1.
         dz=(fp2(i)-ep2(i))*dx
         uro2(i)=(1.+u2(i)-.25*(gpm2+gp2(i)))*dx*dx/(dx*dx+dz*dz)
         wro2(i)=uro2(i)*dz/dx
         uro2(i)=uro2(i)-1.
         fu2(i)=f2(i)+e2(i)
         fo2(i)=f2(i)-e2(i)
         axiim2=axi2(i)
         gpm2=gp2(i)
         if(i.eq.ix-1)then
            axiim2=-1.+2.*axiim2
         endif
 14   continue
      do 15 i=1,ix
         axi=axi1(i)
         gpe=4.*cld*dsqrt(axi*(1.-axi))/pi
         cpu1(i)=2.*uru1(i)
         cpo1(i)=2.*uro1(i)
         cpu2(i)=2.*uru2(i)
         cpo2(i)=2.*uro2(i)
         write(16,*)axi1(i),cpu1(i),cpo1(i)
         write(17,*)axi2(i),cpu2(i),cpo2(i)
         write(25,*)axi1(i),gpe,-gpe
 15   continue
      cl1=2.*g1(ix)
      cl2=2.*g2(ix)
      clt=cl1+cl2
c*****output on listing
      write(6,*)
      write(6,*)'cl1=',cl1
      write(6,*)'cl2=',cl2
      write(6,*)' '
      write(6,*)'  i ','    x1(i) ','    s1(i) ','    w1(i) '
     &     ,'    g1(i)','     x2(i) ','    s2(i) ','    w2(i) '
     &     ,'    g2(i)'
      write(26,*)
      write(26,*)'cl1=',cl1
      write(26,*)'cl2=',cl2
      write(26,*)' '
      write(26,*)'  i ','    x1(i) ','    s1(i) ','    w1(i) '
     &     ,'    g1(i) ','     x2(i) ','    s2(i) ','    w2(i) '
     &     ,'    g2(i)'
c*****output for nxyplot
      do 16 i=1,ix
         axi=x1(i)
         f1e=cld*axi*(1.-axi)/pi
         write(6,1000)i,x1(i),s1(i),w1(i),g1(i),x2(i),s2(i),w2(i),g2(i)
         write(26,1000)i,x1(i),s1(i),w1(i),g1(i),x2(i),s2(i),w2(i),g2(i)
         write(18,*)x1(i),fu1(i),fo1(i)
         write(19,*)x2(i),fu2(i),fo2(i)
         write(22,*)x1(i),s1(i)
         write(23,*)x2(i),s2(i)
         write(24,*)x1(i),f1e
 16   continue
c     add thickness to camber line and small te thickness
!     write(27,*)'x ','z'
!     write(28,*)'x','z'
      do 17 i=1,ix
         ic=ix+1-i
         fu1(ic)=f1(ic)+e1(ic)+tethick*x1(ic)/xte1
         fu2(ic)=f2(ic)+e2(ic)+tethick*(x2(ic)-xle2)/(1.0-xle2)
         write(27,*)x1(ic),fu1(ic)
         write(28,*)x2(ic),fu2(ic)
 17   continue
      write(29,*)x1(1),g1(1)
      write(30,*)x2(1),g2(1)
      do 18 i=2,ix
         fo1(i)=f1(i)-e1(i)-tethick*x1(i)/xte1
         fo2(i)=f2(i)-e2(i)-tethick*(x2(i)-xle2)/(1.0-xle2)
         write(27,*)x1(i),fo1(i)
         write(28,*)x2(i),fo2(i)
         write(29,*)x1(i),g1(i)
         write(30,*)x2(i),g2(i)
 18   continue
c
c      write mses profile
c
      write(31,*)'profile.name'
      write(31,*)XINL,' ',XOUT,' ',YBOT,' ',YTOP
      do 19 i=1,ix
         ic=ix+1-i
         fu1(ic)=f1(ic)+e1(ic)+tethick*x1(ic)/xte1
         write(31,*)x1(ic),fu1(ic)
  19   continue
      do 20 i=2,ix
         fo1(i)=f1(i)-e1(i)-tethick*x1(i)/xte1
         write(31,*)x1(i),fo1(i)
  20   continue
      write(31,*)'999.0 999.0'
      do 21 i=1,ix
         ic=ix+1-i
         fu2(ic)=f2(ic)+e2(ic)+tethick*(x2(ic)-xle2)/(1.0-xle2)
         write(31,*)x2(ic),fu2(ic)
  21   continue
      do 22 i=2,ix
         fo2(i)=f2(i)-e2(i)-tethick*(x2(i)-xle2)/(1.0-xle2)
         write(31,*)x2(i),fo2(i)
  22   continue
c
c
c
      write(6,*)' '
      write(6,*)'ak=',ak,' cl1=',cl1,' cl2=',cl2,' clt=',clt
      write(26,*)' '
      write(26,*)'ak=',ak,' cl1=',cl1,' cl2=',cl2,' clt=',clt
 1000 format(1x,i3,2x,f8.4,2x,f8.4,2x,f8.4,2x,f8.4,2x,f8.4,2x,f8.4
     &     ,2x,f8.4,2x,f8.4,2x,f8.4,2x,f8.4)
 1001 format(f10.2,f10.2,f10.2,f10.2,f10.2,f10.2,f10.2
     &     ,f10.2,f10.2,f10.2,f10.2,f10.2)
c*****input files
      write(6,*)'******input files:'
      write(6,*)'cloptair.data    :data file'
c*****output files
      write(6,*)'******output files:'
      write(6,*)'cloptcp1.dat      :x,Cp1(x) reg. pressure coef'
      write(6,*)'cloptcp2.dat      :x,Cp2(x) reg. pressure coef'
      write(6,*)'cloptxf1.dat      :x,fu1(x),fo1(x)'
      write(6,*)'cloptxf2.dat      :x,fu2(x),fo2(x)'
      write(6,*)'cloptxs1.dat      :x,s1(x)'
      write(6,*)'cloptxs2.dat      :x,s2(x)'
      write(6,*)'cloptxfe.dat      :x,fe(x)'
      write(6,*)'cloptcpe.dat      :x,Cpe(x)'
      write(6,*)'cloptair.out      :listing'
      write(6,*)'cloptxz1.dat      :coordinates of elem 1 for mses'
      write(6,*)'cloptxz2.dat      :coordinates of elem 2 for mses'
      write(6,*)'blade.profile     :coordinates for mses'
      write(6,*)'cloptxg1.dat      :circulation of elem 1'
      write(6,*)'cloptxg2.dat      :circulation of elem 2'
      write(6,*)'blade.profile     :profile geometry for Mses'
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      stop
	end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine dgefa(a,lda,n,ipvt,info)
      integer lda,n,ipvt(1),info
      double precision a(lda,1)
c
c     dgefa factors a double precision matrix by gaussian elimination.
c
c     dgefa is usually called by dgeco, but it can be called
c     directly with a saving in time if  rcond  is not needed.
c     (time for dgeco) = (1 + 9/n)*(time for dgefa) .
c
c     on entry
c
c        a       double precision(lda, n)
c                the matrix to be factored.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c     on return
c
c        a       an upper triangular matrix and the multipliers
c                which were used to obtain it.
c                the factorization can be written  a = l*u  where
c                l  is a product of permutation and unit lower
c                triangular matrices and  u  is upper triangular.
c
c        ipvt    integer(n)
c                an integer vector of pivot indices.
c
c        info    integer
c                = 0  normal value.
c                = k  if  u(k,k) .eq. 0.0 .  this is not an error
c                     condition for this subroutine, but it does
c                     indicate that dgesl or dgedi will divide by zero
c                     if called.  use  rcond  in dgeco for a reliable
c                     indication of singularity.
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas daxpy,dscal,idamax
c
c     internal variables
c
      double precision t
      integer idamax,j,k,kp1,l,nm1
c
c
c     gaussian elimination with partial pivoting
c
      info = 0
      nm1 = n - 1
      if (nm1 .lt. 1) go to 70
      do 60 k = 1, nm1
         kp1 = k + 1
c
c        find l = pivot index
c
         l = idamax(n-k+1,a(k,k),1) + k - 1
         ipvt(k) = l
c
c        zero pivot implies this column already triangularized
c
         if (a(l,k) .eq. 0.0d0) go to 40
c
c           interchange if necessary
c
            if (l .eq. k) go to 10
               t = a(l,k)
               a(l,k) = a(k,k)
               a(k,k) = t
   10       continue
c
c           compute multipliers
c
            t = -1.0d0/a(k,k)
            call dscal(n-k,t,a(k+1,k),1)
c
c           row elimination with column indexing
c
            do 30 j = kp1, n
               t = a(l,j)
               if (l .eq. k) go to 20
                  a(l,j) = a(k,j)
                  a(k,j) = t
   20          continue
               call daxpy(n-k,t,a(k+1,k),1,a(k+1,j),1)
   30       continue
         go to 50
   40    continue
            info = k
   50    continue
   60 continue
   70 continue
      ipvt(n) = n
      if (a(n,n) .eq. 0.0d0) info = n
      return
      end

	subroutine dgesl(a,lda,n,ipvt,b,job)
      integer lda,n,ipvt(1),job
      double precision a(lda,1),b(1)
c
c     dgesl solves the double precision system
c     a * x = b  or  trans(a) * x = b
c     using the factors computed by dgeco or dgefa.
c
c     on entry
c
c        a       double precision(lda, n)
c                the output from dgeco or dgefa.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c        ipvt    integer(n)
c                the pivot vector from dgeco or dgefa.
c
c        b       double precision(n)
c                the right hand side vector.
c
c        job     integer
c                = 0         to solve  a*x = b ,
c                = nonzero   to solve  trans(a)*x = b  where
c                            trans(a)  is the transpose.
c
c     on return
c
c        b       the solution vector  x .
c
c     error condition
c
c        a division by zero will occur if the input factor contains a
c        zero on the diagonal.  technically this indicates singularity
c        but it is often caused by improper arguments or improper
c        setting of lda .  it will not occur if the subroutines are
c        called correctly and if dgeco has set rcond .gt. 0.0
c        or dgefa has set info .eq. 0 .
c
c     to compute  inverse(a) * c  where  c  is a matrix
c     with  p  columns
c           call dgeco(a,lda,n,ipvt,rcond,z)
c           if (rcond is too small) go to ...
c           do 10 j = 1, p
c              call dgesl(a,lda,n,ipvt,c(1,j),0)
c        10 continue
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas daxpy,ddot
c
c     internal variables
c
      double precision ddot,t
      integer k,kb,l,nm1
c
      nm1 = n - 1
      if (job .ne. 0) go to 50
c
c        job = 0 , solve  a * x = b
c        first solve  l*y = b
c
         if (nm1 .lt. 1) go to 30
         do 20 k = 1, nm1
            l = ipvt(k)
            t = b(l)
            if (l .eq. k) go to 10
               b(l) = b(k)
               b(k) = t
   10       continue
            call daxpy(n-k,t,a(k+1,k),1,b(k+1),1)
   20    continue
   30    continue
c
c        now solve  u*x = y
c
         do 40 kb = 1, n
            k = n + 1 - kb
            b(k) = b(k)/a(k,k)
            t = -b(k)
            call daxpy(k-1,t,a(1,k),1,b(1),1)
   40    continue
      go to 100
   50 continue
c
c        job = nonzero, solve  trans(a) * x = b
c        first solve  trans(u)*y = b
c
         do 60 k = 1, n
            t = ddot(k-1,a(1,k),1,b(1),1)
            b(k) = (b(k) - t)/a(k,k)
   60    continue
c
c        now solve trans(l)*x = y
c
         if (nm1 .lt. 1) go to 90
         do 80 kb = 1, nm1
            k = n - kb
            b(k) = b(k) + ddot(n-k,a(k+1,k),1,b(k+1),1)
            l = ipvt(k)
            if (l .eq. k) go to 70
               t = b(l)
               b(l) = b(k)
               b(k) = t
   70       continue
   80    continue
   90    continue
  100 continue
      return
      end

ccccccccccccccccccccccccccccccccccccc
      DOUBLE PRECISION FUNCTION DDOT(N,DX,INCX,DY,INCY)
*     .. Scalar Arguments ..
      INTEGER INCX,INCY,N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION DX(*),DY(*)
*     ..
*
*  Purpose
*  =======
*
*     forms the dot product of two vectors.
*     uses unrolled loops for increments equal to one.
*     jack dongarra, linpack, 3/11/78.
*     modified 12/3/93, array(1) declarations changed to array(*)
*
*
*     .. Local Scalars ..
      DOUBLE PRECISION DTEMP
      INTEGER I,IX,IY,M,MP1
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC MOD
*     ..
      DDOT = 0.0d0
      DTEMP = 0.0d0
      IF (N.LE.0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) GO TO 20
*
*        code for unequal increments or equal increments
*          not equal to 1
*
      IX = 1
      IY = 1
      IF (INCX.LT.0) IX = (-N+1)*INCX + 1
      IF (INCY.LT.0) IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
          DTEMP = DTEMP + DX(IX)*DY(IY)
          IX = IX + INCX
          IY = IY + INCY
   10 CONTINUE
      DDOT = DTEMP
      RETURN
*
*        code for both increments equal to 1
*
*
*        clean-up loop
*
   20 M = MOD(N,5)
      IF (M.EQ.0) GO TO 40
      DO 30 I = 1,M
          DTEMP = DTEMP + DX(I)*DY(I)
   30 CONTINUE
      IF (N.LT.5) GO TO 60
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
          DTEMP = DTEMP + DX(I)*DY(I) + DX(I+1)*DY(I+1) +
     +            DX(I+2)*DY(I+2) + DX(I+3)*DY(I+3) + DX(I+4)*DY(I+4)
   50 CONTINUE
   60 DDOT = DTEMP
      RETURN
      END
cccccccccccccccccccccccccccccccccccccccccccccccccc
      INTEGER FUNCTION IDAMAX(N,DX,INCX)
*     .. Scalar Arguments ..
      INTEGER INCX,N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION DX(*)
*     ..
*
*  Purpose
*  =======
*
*     finds the index of element having max. absolute value.
*     jack dongarra, linpack, 3/11/78.
*     modified 3/93 to return if incx .le. 0.
*     modified 12/3/93, array(1) declarations changed to array(*)
*
*
*     .. Local Scalars ..
      DOUBLE PRECISION DMAX
      INTEGER I,IX
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC DABS
*     ..
      IDAMAX = 0
      IF (N.LT.1 .OR. INCX.LE.0) RETURN
      IDAMAX = 1
      IF (N.EQ.1) RETURN
      IF (INCX.EQ.1) GO TO 20
*
*        code for increment not equal to 1
*
      IX = 1
      DMAX = DABS(DX(1))
      IX = IX + INCX
      DO 10 I = 2,N
          IF (DABS(DX(IX)).LE.DMAX) GO TO 5
          IDAMAX = I
          DMAX = DABS(DX(IX))
    5     IX = IX + INCX
   10 CONTINUE
      RETURN
*
*        code for increment equal to 1
*
   20 DMAX = DABS(DX(1))
      DO 30 I = 2,N
          IF (DABS(DX(I)).LE.DMAX) GO TO 30
          IDAMAX = I
          DMAX = DABS(DX(I))
   30 CONTINUE
      RETURN
      END

cccccccccccccccccccccccccccccccc
      SUBROUTINE DAXPY(N,DA,DX,INCX,DY,INCY)
*     .. Scalar Arguments ..
      DOUBLE PRECISION DA
      INTEGER INCX,INCY,N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION DX(*),DY(*)
*     ..
*
*  Purpose
*  =======
*
*     constant times a vector plus a vector.
*     uses unrolled loops for increments equal to one.
*     jack dongarra, linpack, 3/11/78.
*     modified 12/3/93, array(1) declarations changed to array(*)
*
*
*     .. Local Scalars ..
      INTEGER I,IX,IY,M,MP1
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC MOD
*     ..
      IF (N.LE.0) RETURN
      IF (DA.EQ.0.0d0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) GO TO 20
*
*        code for unequal increments or equal increments
*          not equal to 1
*
      IX = 1
      IY = 1
      IF (INCX.LT.0) IX = (-N+1)*INCX + 1
      IF (INCY.LT.0) IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
          DY(IY) = DY(IY) + DA*DX(IX)
          IX = IX + INCX
          IY = IY + INCY
   10 CONTINUE
      RETURN
*
*        code for both increments equal to 1
*
*
*        clean-up loop
*
   20 M = MOD(N,4)
      IF (M.EQ.0) GO TO 40
      DO 30 I = 1,M
          DY(I) = DY(I) + DA*DX(I)
   30 CONTINUE
      IF (N.LT.4) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,4
          DY(I) = DY(I) + DA*DX(I)
          DY(I+1) = DY(I+1) + DA*DX(I+1)
          DY(I+2) = DY(I+2) + DA*DX(I+2)
          DY(I+3) = DY(I+3) + DA*DX(I+3)
   50 CONTINUE
      RETURN
      END

cccccccccccccccccccccccccccccccc
      SUBROUTINE DSCAL(N,DA,DX,INCX)
*     .. Scalar Arguments ..
      DOUBLE PRECISION DA
      INTEGER INCX,N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION DX(*)
*     ..
*
*  Purpose
*  =======
**
*     scales a vector by a constant.
*     uses unrolled loops for increment equal to one.
*     jack dongarra, linpack, 3/11/78.
*     modified 3/93 to return if incx .le. 0.
*     modified 12/3/93, array(1) declarations changed to array(*)
*
*
*     .. Local Scalars ..
      INTEGER I,M,MP1,NINCX
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC MOD
*     ..
      IF (N.LE.0 .OR. INCX.LE.0) RETURN
      IF (INCX.EQ.1) GO TO 20
*
*        code for increment not equal to 1
*
      NINCX = N*INCX
      DO 10 I = 1,NINCX,INCX
          DX(I) = DA*DX(I)
   10 CONTINUE
      RETURN
*
*        code for increment equal to 1
*
*
*        clean-up loop
*
   20 M = MOD(N,5)
      IF (M.EQ.0) GO TO 40
      DO 30 I = 1,M
          DX(I) = DA*DX(I)
   30 CONTINUE
      IF (N.LT.5) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
          DX(I) = DA*DX(I)
          DX(I+1) = DA*DX(I+1)
          DX(I+2) = DA*DX(I+2)
          DX(I+3) = DA*DX(I+3)
          DX(I+4) = DA*DX(I+4)
   50 CONTINUE
      RETURN
      END

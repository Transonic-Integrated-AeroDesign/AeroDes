      program smoothpolar
      implicit none
      integer lxx,kx,inpolar,km,kfirst,kdum,kmfirst,ice,kxold,kxx
      integer kskip,kp,lc,k,l
      parameter(lxx=101)
      integer kxtrm(lxx)
      real cx(lxx),cz(lxx),cdp(lxx),cq(lxx),inc(lxx)
      real eps,pi,degrad,prod,dcz,czfirst,dczm,dczp,incdum,incd
      real czdum,cxdum,cqdum
      character*4 title(18)
      data cdp/lxx*0./
      data kxtrm/lxx*0/
      open(unit=13,file='smoothpolar.in',form='formatted')
      open(unit=14,file='smoothpolar.out',form='formatted')
      open(unit=15,file='smoothpolar.alcl',form='formatted')
      open(unit=16,file='smoothpolar.alcq',form='formatted')
      open(unit=17,file='smoothpolar.alcd',form='formatted')
c*****constants
      eps=1.e-7
      pi=2.*asin(1.)
      degrad=pi/180.
      write(6,*)' is the polar from Xfoil? Y/N=1/0'
      read(5,*)inpolar
      write(6,*)' inpolar=',inpolar
      write(6,*)
      write(6,*)'******rough polar data from input file'
      if(inpolar.ne.1)goto 1
c*****polar data
      write(6,*)'******profile polar:'
      write(6,*)
      read(13,1000)title
      write(6,1000)title
      write(14,1000)title
      read(13,1000)title
      write(6,1000)title
      write(14,1000)title
      read(13,1000)title
      write(6,1000)title
      write(14,1000)title
      read(13,1000)title
      write(6,1000)title
      write(14,1000)title
      read(13,1000)title
      write(6,1000)title
      write(14,1000)title
      read(13,1000)title
      write(6,1000)title
      write(14,1000)title
      read(13,1000)title
      write(6,1000)title
      write(14,1000)title
      read(13,1000)title
      write(6,1000)title
      write(14,1000)title
      read(13,1000)title
      write(6,1000)title
      write(14,1000)title
      read(13,1000)title
      write(6,1000)title
      write(14,1000)title
      read(13,1000)title
      write(6,1000)title
      write(14,1000)title
 1    continue
      read(13,1000)title
      write(6,1000)title
      write(14,1000)title
      prod=1.
      km=1
      kfirst=0
      do 2 k=1,lxx
         kdum=k
         read(13,*,end=3)inc(k),cz(k),cx(k),cdp(k),cq(k)
         write(6,1001)inc(k),cz(k),cx(k),cdp(k),cq(k)                  
         kxtrm(k)=0
         if(k.eq.1)goto 2
         dcz=cz(k)-cz(km)
         prod=prod*dcz
         if(prod.lt.-eps)then
            kxtrm(km)=km
            if(kfirst.eq.0)then
               kfirst=1
               kmfirst=km
               czfirst=cz(km)
            endif
         endif
         prod=sign(1.,dcz)
         km=k
 2    continue
 3    continue
      write(6,*)
      write(6,*)'******first extremum of polar'
      write(6,*)' kxtrm=',kmfirst,' cz(kxtrm)=',czfirst
      write(6,*)
      kx=kdum-1
      if(kx.eq.lxx-1)then
         write(6,*)' attention: check if all data has been read;'
     &        ,' continuing/exiting=1/0?'
         read(5,*)ice
         if(ice.eq.0)stop
      endif
      write(6,*)
      write(6,*)'******smooth polar with interpolated data at extrema'
      kxold=kx
      kxx=kx
      kskip=0
      kxtrm(kx+1)=0
      write(6,1000)title
      do 5 k=1,kx
         if(kskip.eq.1)then
            kskip=0
            kxtrm(k)=k
            kxtrm(k-1)=0
            goto 5
         endif
         kp=k+1
         if(kp.gt.kx)kp=kx
         km=k-1
         if(km.lt.1)km=1
         prod=1.
         if(k.gt.1.and.k.lt.kx)then
            dczm=cz(k)-cz(km)
            dczp=cz(kp)-cz(k)
            prod=dczm*dczp
         endif
         if(prod.lt.-eps.and.kxtrm(k).ne.0)then
            incdum=.5*(inc(km)+inc(k))
            czdum=.5*(cz(km)+cz(k))
            cxdum=.5*(cx(km)+cx(k))
            cqdum=.5*(cq(km)+cq(k))
            do 4 lc=k+1,kx
               l=kx+k+1-lc
               inc(l+2)=inc(l)
               cz(l+2)=cz(l)
               cx(l+2)=cx(l)
               cq(l+2)=cq(l)
 4          continue
            inc(k+2)=.5*(inc(k)+inc(k+3))
            cz(k+2)=.5*(cz(k)+cz(k+3))
            cx(k+2)=.5*(cx(k)+cx(k+3))
            cq(k+2)=.5*(cq(k)+cq(k+3))
            inc(k+1)=inc(k)
            cz(k+1)=cz(k)
            cx(k+1)=cx(k)
            cq(k+1)=cq(k)
            inc(k)=incdum
            cz(k)=czdum
            cx(k)=cxdum
            cq(k)=cqdum
            write(6,*)' k=',k,' inc(k)=',inc(k),' cz(k)=',cz(k)
     &           ,' cx(k)=',cx(k),' cq(k)=',cq(k)
            write(14,1001)inc(k),cz(k),cx(k),cdp(k),cq(k)
            write(6,*)' k=',k+1,' inc(k)=',inc(k+1),' cz(k)=',cz(k+1)
     &           ,' cx(k)=',cx(k+1),' cq(k)=',cq(k+1)
            write(14,1001)inc(k+1),cz(k+1),cx(k+1),cdp(k+1),cq(k+1)
            kxx=kx+2
            kskip=1
            goto 5
         endif
         incd=inc(k)
         write(6,*)' k=',k,' inc(k)=',inc(k),' cz(k)=',cz(k)
     &        ,' cx(k)=',cx(k),' cq(k)=',cq(k)
         write(14,1001)inc(k),cz(k),cx(k),cdp(k),cq(k)
 5    continue
      if(inc(kxx).ge.89.)then
         kx=kxx
         goto 6
      endif
      kx=kxx+2
      inc(kx)=90.
      cz(kx)=0.
      cx(kx)=1.
      cq(kx)=0.
      inc(kx-1)=2.*inc(kx-2)-inc(kx-3)
      cz(kx-1)=cz(kx-2)+(cz(kx)-cz(kx-2))
     &     *(inc(kx-1)-inc(kx-2))/(inc(kx)-inc(kx-2))
      cx(kx-1)=cx(kx-2)+(cx(kx)-cx(kx-2))
     &     *(inc(kx-1)-inc(kx-2))/(inc(kx)-inc(kx-2))
      cq(kx-1)=cq(kx-2)+(cq(kx)-cq(kx-2))
     &     *(inc(kx-1)-inc(kx-2))/(inc(kx)-inc(kx-2))
 6    continue
      do 7 k=kxold+1,kx
         write(6,*)' k=',k,' inc(k)=',inc(k),' cz(k)=',cz(k)
     &        ,' cx(k)=',cx(k),' cq(k)=',cq(k)
         write(14,1001)inc(k),cz(k),cx(k),cdp(k),cq(k)
 7    continue
      kxtrm(kx)=0
      write(6,*)'******extrema pointer:'
      write(6,*)'kxtrm(k)=',(kxtrm(k),k=1,kx)
      write(6,*)'******kx=',kx
      do 8 k=1,kx
         write(15,*)inc(k),cz(k)
         write(16,*)inc(k),cq(k)
         write(17,*)inc(k),cx(k)
 8    continue
 1000 format(18a4)
 1001 format(1x,f7.3,3x,f6.4,3x,f6.4,3x,f7.5,3x,f7.4)
      end

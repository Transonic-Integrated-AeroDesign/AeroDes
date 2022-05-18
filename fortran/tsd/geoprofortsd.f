      program geoprofortsd
c     geoprofortsd.f reads a profile thickness distribution and 
c     adds a camber distribution that corresponds to a flywing
c     profile plus a negative parabolic camber and writes the
c     camber and thickness at intermediate points to be used 
c     in the code tsd.f between the leading edge xle(j) and 
c     the trailing edge xte(j) of a 3-D wing
      implicit none
      integer ixx,lxx,ile,ite,i,ix,ic,lx,l
      real pi,eps,dtet,xl,el,xlm,dm,dmplus,dl
      parameter(ixx=101,lxx=201)
      real xi(ixx),zu(ixx),zo(ixx)
      real x(lxx),d(lxx),e(lxx)
      open(unit=11,file='geoprofortsd.data',form='formatted')
      open(unit=12,file='geoprofortsd.xzu',form='formatted')
      open(unit=13,file='geoprofortsd.xzmses',form='formatted')
      open(unit=14,file='geoprofortsd.xde',form='formatted')
      pi=2.0*asin(1.0)
      eps=1.0e-6
c     read indices of leading edge and trailing edge of tsd wing ile, ite
c     read flywing camber dm (>0) and extra parabolic camber dmplus (<0)
      read(11,*)ile,ite
      read(11,*)dm
      read(11,*)dmplus
      write(6,*)
      write(6,*)' ile=',ile,'ite=',ite
      write(6,*)
      write(6,*)' added camber distribution of flywing dm=',dm
      write(6,*)'  added negative parabolic camber dmplus=',dmplus
      write(6,*)
c     read thickness distribution of wing profile (symmetric profile)
      do 1 i=1,ixx
         read(12,*,end=2)xi(i),zu(i)
c         write(6,*)' i=',i,'xi(i)=',xi(i),'zu(i)=',zu(i)
         ix=i
 1    continue
 2    continue
      do 3 ic=1,ix
         i=ix+1-ic
         zo(i)=-zu(i)
         write(13,*)xi(i),zo(i)
 3    continue
      do 4 i=2,ix
         write(13,*)xi(i),zu(i)
 4    continue
      lx=ite-ile+1
      write(6,*)' ix=',ix,'  lx=',lx
      dtet=pi/(ite-ile)
      xlm=0.0
      do 6 l=1,lx-1
         xl=0.5*(1.0-cos((l-1+0.5)*dtet))
         x(l)=xl
         xlm=xl
         do 5 i=2,ix
            if(x(l).ge.xi(i-1).and.x(l).le.xi(i))then
               dl=dm/3.*x(l)*(7.-8.*x(l))*(1.-x(l))
     &              +4.*dmplus*x(l)*(1.-x(l))
               d(l)=dl
               el=zu(i-1)+(x(l)-xi(i-1))*(zu(i)-zu(i-1))/(xi(i)-xi(i-1))
               e(l)=el
            endif
 5       continue
         write(14,*)x(l),d(l),e(l)
 6    continue
      write(6,*)
      write(6,*)'geoprofortsd.data          :data file'
      write(6,*)'geoprofortsd.xzu           :profile thickness input'
      write(6,*)'geoprofortsd.xzmses        :profile geometry for Xfoil'
      write(6,*)'geoprofortsd.xde           :profile camber and'
     &,' thickness for tsd at intermediate points'
      end




      


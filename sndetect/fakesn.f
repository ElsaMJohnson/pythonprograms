c       Fakesn creates a 1000x1000 grid with n number of
c       stars (SNe) strewn randomly across the image.
c       The peak brightness of the psf is chosen randomly from a 
c       0 to 500.
c       Outputs 2 files. Use fort.91 and rtextimage to make a fake
c       image in iraf. make sure you pick no for header option
c       

      PROGRAM FAKESN

      real*8 grid(1000,1000),sn(1000,1000)
      real*8 dcsky,noise
      real*8 gasdev,sz,const,seed,seed2,nsn
      real*8 lm,SNR
      integer px(1000),py(1000),tsteps,gx(1000)
      integer gy(1000),sty,stx,xx,yy
    
      write(*,*)"input seed to start"
      read(*,*) seed
      seed2=5*seed+192
      seed3=sqrt(seed2)
      seed4=seed**2-1
c     parameters:
c     average value of background sky. Noise will be added later.
      dcsky=100
c     number of events to generate
      nsn=100
c     size of search box for events
      sz=10
c     size of gaussian 
      sig=1
      pi=3.141592653589793
c     number of times area is randomly searched.
      tsteps=10

c     zero out arrays      
      do i=1,1000
      do j=1,1000
        sn(i,j)=0.
        grid(i,j)=0.
      end do
      end do
      nd=0.
      sumsq=0.
c     Loop for number of observations (1 per night)
      do t=1,tsteps
c     randomly create SN and get rid of first pick b/c it sucks
      tmp1=1+1000*ran3(int(seed))
      tmp2=1+1000*ran3(int(seed2))
c     SN loop for # of events px and py are the randomly selected
c     coordinates. I have avoided generating events on the edges
c     of the image hence the 'goto' statements
      do k=1,nsn
 50     px(k)=1+1000*ran3(int(seed))
 51     py(k)=1+1000*ran3(int(seed2))
      if (px(k).ge.(999-sz).or.px(k).le.(sz+1)) then
         goto 50
      else if (py(k).ge.(999-sz).or.py(k).le.(sz+1)) then 
         goto 51 
      end if  
c     select a peak brightness for SN
      const=500*ran3(int(seed2))
      lm=2*sz
c     psf loop: generating pixel values for psf
c     psf is a double gaussian (1st gaussian has a sigma = 1, 2nd
c     has sigma=2)
      do x=1,lm
      do y=1,lm
        c1=px(k)-sz+x
        c2=py(k)-sz+y
        
        sn(c1,c2)=sn(c1,c2)+const*(1/sig)*(1/(2*pi))*
     *   (exp(-((x-sz)**2+(y-sz)**2)/(2*sig**2))
     *   + (exp(-((x-sz)**2+(y-sz)**2)/(8*sig**2))/2))
c       write(*,*) c1,c2,sn(c1,c2)
       end do
       end do
       end do
c      pick a random 10x10 grid to observe
 52     gx(t)=1+1000*ran3(int(seed3))
 53     gy(t)=1+1000*ran3(int(seed4))
      if (gx(t).ge.(994).or.gx(t).le.(6)) then
         goto 52
      else if (gy(t).ge.(994).or.gy(t).le.(6)) then
         goto 53
      end if 
c      Sum the noise, sky and possible SN counts in each grid
        write(*,*) made it here
       do j=1,1000
       do i=1,1000
c       adding noise, sky and sn counts to each pixel
        call gausran3(int(seed),gasdev)
        noise=gasdev*sqrt(dcsky)
        grid(i,j)=dcsky+noise+sn(i,j)
c       use this for gnuplot
        write(90,*) i,j,grid(i,j)
c       use this below for rtextimage in iraf.
        write(91,*) grid(i,j)
        end do
        end do
c       Sum counts in detector grid. A detection means
c       an SNR of 3 or greater.       
        sumsq=0.
        do j=1,10
        do i=1,10
        stx = gx(t)-6
        sty = gy(t)-6
        xx=stx+i
        yy=sty+j
        sumsq=grid(xx,yy)+sumsq 
        end do
        end do
        write(*,*)'flux',sumsq-10000
c       The background to subtract is 100
        SNR= (sumsq-100*10*10)/sqrt(sumsq)
        write(*,*)"SNR",SNR
        if (SNR.gt.5) then
        write(*,*) "detection",sumsq,gx(t),gy(t)
        nd=nd+1
        end if
        end do 
        write(*,*) "number of detections",nd
       END

c**************************************************

        SUBROUTINE GAUSRAN3(DUMMY,GASDEV)
c       This program of gaussian generated numbers from
c       Numerical recipies called gasdev

        REAL*8 fac,gset,rsq,v1,v2,gasdev
c        REAL*8 ran3
        INTEGER iset,dummy
        SAVE iset,gset
        DATA iset/0/
        if (dummy.lt.0) iset=0
 1      if (iset.eq.0) then
          v1=(2.*ran3(dummy)-1)
          v2=(2.*ran3(dummy)-1)
          rsq=v1**2+v2**2
          if (rsq.gt.1..or.rs1.eq.0.) goto 1
          fac=sqrt(-2.*log(rsq)/rsq)
          gset=v1*fac
          gasdev=v2*fac
          iset=1
        else
          gasdev=gset
          iset=0
        end if
        RETURN
        END

      Real*8 FUNCTION RAN3(IDUM)
      Implicit None
      Integer iff,idum,i,inext,inextp,k,ii
      Real*8 mbig, mseed, mz, fac, ma(55), mj, mk
      save inext,INEXTP, ma

      PARAMETER (MBIG=4000000.D0,MSEED=1618033.D0,MZ=0.D0,FAC=2.5D-7)
CC      PARAMETER (MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1.E-9)
      DATA IFF /0/
        IF(IDUM.LT.0.OR.IFF.EQ.0)THEN
         IFF=1
         MJ=MSEED-dble(IABS(IDUM))
         MJ=MOD(MJ,MBIG)
         MA(55)=MJ
         MK=1
         DO I=1,54,1
            II=MOD(21*I,55)
            MA(II)=MK
            MK=MJ-MK
            IF(MK.LT.MZ)MK=MK+MBIG
            MJ=MA(II)
         END DO
         DO K=1,4,1
            DO I=1,55,1
               MA(I)=MA(I)-MA(1+MOD(I+30,55))
               IF(MA(I).LT.MZ)MA(I)=MA(I)+MBIG
            END DO
         END DO
         INEXT=0
         INEXTP=31
         IDUM=1
      ENDIF
      INEXT=INEXT+1
      IF(INEXT.ge.56)INEXT=1
      INEXTP=INEXTP+1
      IF(INEXTP.ge.56)INEXTP=1
      MJ=MA(INEXT)-MA(INEXTP)
      IF(MJ.LT.MZ)MJ=MJ+MBIG
      MA(INEXT)=MJ
      RAN3=MJ*FAC
      RETURN
      END

    
    


     

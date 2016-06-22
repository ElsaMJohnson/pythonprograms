        PROGRAM READPOLYPHOT
        
c       This is based off of an old program I wrote to analyze dr7
c       magnitudes for SDSS images (asinhmag).
c       This uses the flux values in the polyphot output file rather
c       than the calculated magnitudes because some can be undefined
c       due to the low flux when doing grid analysis.

c       bmag is the bluer of the two colors and rmag is the redder
c       Input is a parameter file of the following (1 item per line)
c       bluer file; 
c       redder file;
c       bluer zeropoint; 
c       redder zeropoint; 
c       bluer extinction (Ab) 
c       redder extinction (Ar)
c       name of output file
c       All of my zeropoints contain the exposure time information. 
c       If yours do not, then you will have to add it to the program
c       unless the images have been calibrated to counts/s

        real*8 Ab,Ar,numb,bx(400000),by(400000)
        real*8 rx(400000),ry(400000),bsky(400000),bsig(400000)
        real*8 rsky(400000),rsig(400000),rtmp,btmp
        real*8 bsum(400000),bar(400000),bflx(400000),rsum(400000)
        real*8 rar(400000),rflx(400000),bexp,rexp,bzp,rzp
        real*8 berr(400000),rerr(400000),bmag(400000),rmag(400000)
        character*28 bfile,rfile,outfile,a(800000),galfile

        write(*,*) 'enter in galaxy info file'
        read(*,150) galfile       
 150    format(A28)
        open(unit=19,file=galfile,form='formatted',status='old')
        
        read(19,*) bfile
        read(19,*) rfile
C       Note the zeropoints should contain the time encoded in the data
C       too so you need to leave it out of the equation.
c       read(19,*) bexp
c       read(19,*) rexp
        read(19,*) bzp
        read(19,*) rzp
        read(19,*) Ab
        read(19,*) Ar
        read(19,*) outfile

        open(unit=9,file=bfile,status='old')
        open(unit=11,file=rfile,status='old')
        open(unit = 16,file=outfile,status='unknown')
        write(16,*)'X    Y    B    R    Berr   Rerr'

        nn=0
 4      continue
        read(9,*,end=5000) a(nn)
        nn=nn+1
        goto 4
 5000   continue
        rewind(9)

        write(*,*) nn
c       Number of data (10 lines per data chunk)
        numb=(nn-82)/10. 
        write(*,*) numb
        read(9,101)
 101    format(/////////////////////////////////////////)
        read(9,107)
 107    format(////////////////////////////////////////)
        read(11,101)
        read(11,107)
        do i=1,numb
           write(*,*) 'i',i
 103       format(f14.3,f11.3)
           read(9,103) bx(i),by(i)
           write(*,*) bx(i),by(i)
           read(9,104) bsky(i),bsig(i)
           read(11,103) rx(i),ry(i)
           read(11,104) rsky(i),rsig(i)
 104       format(f10.2,f10.1/)
           read(9,*) bsum(i),bar(i),bflx(i)
           read(9,106)
           read(11,*) rsum(i),rar(i),rflx(i)
           read(11,106)
c           write(*,*) bsky(i),rsky(i),bsig(i),rsig(i)
 106       format(/////)

        gain=10.
        if (bflx(i).gt.0) then
         bmag(i) = -2.5*log10(bflx(i))+bzp-Ab
         btmp = ((bflx(i)/gain)+bar(i)*bsig(i)**2+(bar(i)*bsig(i))**2)
         berr(i)=1.0857*(btmp**.5)/bflx(i)
        else if (bflx(i).le.0) then
         bmag(i) = 50.0
         berr(i) = 10.0
        end if

        if (rflx(i).gt.0) then
         rmag(i) = -2.5*log10(rflx(i))+rzp-Ar
         rtmp = (rflx(i)/gain) +rar(i)*rsig(i)**2+(rar(i)*rsig(i))**2.
         rerr(i) = 1.0857*(rtmp**.5)/rflx(i)
        else if (rflx(i).le.0) then
         rmag(i) = 40.0
         rerr(i) = 10.0
        end if     

 110    format(2f9.3,1x,4(f9.5,1x))
        write(16,110) rx(i),ry(i),bmag(i),rmag(i),berr(i),rerr(i)
        
        end do
      
        END


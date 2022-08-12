      program dopcpisofcst
c
c Purpose: combine SVR iso components and original iso component into 
c          a new iso components which can be used for retrieving 
c          PCP forecast
c
      parameter (nx=144,my=72,nt=504)
      parameter (ntall=516,nv=20,md=nx*my)
      character*200 pcpisosvr,eofano_ana2_file,
     &              eofpca_ana2_file,eofev_ana2_file,
     &              ydatafriso_ana2_file,
     &              avemonfile
      integer nbegin
      real*4 gh(nx,my)
      real pp(nt)
      real pca(ntall,nv),weigen(nt),treof(md,nt),origpca(nt,nt)
c      real zeall2(nv,ntall)
      real zeall2(nt,ntall)
      real y(md,ntall)
      real rlon1,rlon2,rlat1,rlat2
c
      namelist /dopcpisofcstlst/ pcpisosvr,eofano_ana2_file,
     >                     eofpca_ana2_file,eofev_ana2_file,
     >                     ydatafriso_ana2_file,avemonfile,
     >                     rlon1,rlon2,rlat1,rlat2,
     >                     nbegin
c
c======================================================================
c
      open(66,file='dopcpisofcst.nmlst',form='formatted')
      read(66,dopcpisofcstlst,end=110)
  110 continue
      print dopcpisofcstlst
c
c Pick up the selected domain
c
      rlondiff=360./float(nx)
      rlatdiff=177.5/float(my-1)
      npx1=nint((rlon1-1.25)/rlondiff)+1
      npx2=nint((rlon2-1.25)/rlondiff)+1
      mpy1=nint((rlat1+88.75)/rlatdiff)+1
      mpy2=nint((rlat2+88.75)/rlatdiff)+1
      print *,'npx1,npx2= ',npx1,npx2
      print *,'mpy1,mpy2= ',mpy1,mpy2
      newpx=(npx2-npx1+1)
      newpy=(mpy2-mpy1+1)
      print *,'newpx,newpy= ',newpx,newpy
      newmd=newpx*newpy
      print *,'new spatial dimension, newmd = ',newmd
c
c      pcpisosvr='/geps2/chtseng/prog/PCP_EAmonsoon_svr_fcst_516.dat'
      newpxpy4=newpx*newpy*4
      nt4=nt*4
      open(11,file=trim(pcpisosvr),form='formatted',
     &     status='old')
      open(12,file=trim(eofano_ana2_file),access='direct',
     &     form='unformatted',recl=newpxpy4,status='old')
      open(14,file=trim(eofpca_ana2_file),access='direct',
     &     form='unformatted',recl=nt4,status='old')
      open(15,file=trim(eofev_ana2_file),
     &     form='formatted',status='old')
      open(16,file=trim(ydatafriso_ana2_file),access='direct',
     &     form='unformatted',recl=newpxpy4,status='unknown')
c
c      open(36,file='testread.dat',form='formatted',status='unknown')
c
      do n=1,ntall
      read(11,*) (pca(n,iv),iv=1,nv)
c      write(36,*) (pca(n,iv),iv=1,nv)
      enddo
c
      call getflds2(12,treof,nt,md,newmd,gh,newpx,newpy)
c
      do nc=1,nt
      call getpca(14,pp,nt,nc)
      do n=1,nt
      origpca(n,nc)=pp(n)
      enddo
      enddo
c
      do n=1,nt
      read (15,155) weigen(n)
  155 format(1x,f20.8)
      enddo
c
      print *,'before pca(505,1)= ',pca(505,1)
      id=1
      do n=nbegin,ntall
c      print *,'n-nbegin+1=',n-nbegin+1
      do nc=1,nv
      pca(n,nc)=pca(n,nc)/sqrt(weigen(nc))
      zeall2(nc,n-nbegin+1)=pca(n,nc)
      enddo
      do nc=nv+1,nt
      if (weigen(nc).le.0.0) then
      zeall2(nc,n-nbegin+1)=0.0
      else
      zeall2(nc,n-nbegin+1)=0.1*origpca(n,nc)/sqrt(weigen(nc))
      endif
      enddo
      id=id+1
      enddo
c
      print *,'weigen(1)= ',weigen(1)
      print *,'pca(505,1)= ',pca(505,1)
      print *,'zeall2(1,10)= ',zeall2(1,10)
      print *,'treof(10,1)= ',treof(10,1)
c
c recompose the new data matrix y
c
      newnt=ntall-nbegin+1
      print *,'newnt= ',newnt
      call mtxmlp2 (treof,md,newmd,nt,nt,zeall2,nt,nt,ntall,newnt,y)
c 
      call outpcp(16,y,md,newmd,ntall,newnt,gh,newpx,newpy,avemonfile)
c
      close(11)
      close(12)
      close(13)
      close(14)
      close(15)
      close(16)
c
      stop
      end
c
c
      subroutine getflds2(ifid,reof,nt,md,newmd,gh,nx,my)
      integer nt,md,nx,my,ifid
      real reof(md,nt),gh(nx,my)
      real xx(newmd,nt)
      character*120 eofano_ana2_file
c
c      nxmy4=nx*my*4
c      open(12,file=trim(eofano_ana2_file),access='direct',
c     &     form='unformatted',recl=nxmy4,status='unknown')
c
      do n=1,nt
      nrec=n
      read(ifid,rec=nrec) gh
c
      m=0
      do j=1,my
      do i=1,nx
        m=m+1
        xx(m,n)=gh(i,j)
      enddo
      enddo
      enddo ! end loop of n, time series
c
      do n=1,nt
      do m=1,newmd
      reof(m,n)=xx(m,n)
      enddo
      enddo
      print *,'reof(1,2)=',reof(1,2)
      print *,'xx(1,2)=',xx(1,2)
c
      return
      end
c
c
      subroutine getpca(ifid,pca,nt,ncheck)
      integer ncheck,nt,ifid
      real pca(nt)
c
      read(ifid,rec=ncheck) pca
c
      return
      end
c
c
      subroutine mtxmlp2 (a,malarge,ma,nalarge,na,
     &                    b,mblarge,mb,nblarge,nb,
     &                    c)
c
c  matrix multiply routine
c
c  ***input***
c
c  a:  first array in matrix product
c  ma: row order of matrix a
c  na: column order of matrix a
c  b:  second array in matrix product
c  mb: row order of matrix b
c  nb: column order of matrix b
c
c  NOTE: na must be equal to mb
c
c  malarge,nalarge,mblarge,nblarge
c        : total dimension declaration in main program
c          other than working dimension
c          for preventing Fortran get the wrong memory address
c
c
c ***output***
c
c  c: matrix product
c
c **************************************************************
c
      integer ma,na,mb,nb
      integer malarge,nalarge,mblarge,nblarge
c      dimension a(ma,na),b(mb,nb),c(ma,nb)
      dimension a(malarge,nalarge),b(mblarge,nblarge),
     &          c(malarge,nblarge)
c
c      call zilch(c,ma*nb)
      call zilch(c,malarge*nblarge)
c
c check!
c
      if (na.ne.mb) then
        print *,'Cannot execute matrix mulitiplication!'
        print *,'Check the dimension of Matrix!'
        stop
      endif
c
      do 3 k=1,na
      do 2 j=1,nb
      do 1 i=1,ma
c      if (a(i,k).eq.-999.) goto 1
c      if (b(k,j).eq.-999.) goto 2
      if (a(i,k).ne.-999.and.b(k,j).ne.-999.) then
        c(i,j)= c(i,j)+a(i,k)*b(k,j)
      endif
    1 continue
    2 continue
    3 continue
      return
      end
c
c
      subroutine zilch (x,m)
c
c  subroutine to zero out arrays
c
c  ***input***
c
c  x: input array to zilch
c  m: number of elements to zilch
c
c  ***output***
c
c  x: zeroed array
c
c ****************************************************************
c
      real x(m)
      do 1 i=1,m
      x(i)= 0.0
    1 continue
      return
      end
c
c
      subroutine outpcp(ifid,y,md,newmd,nt,newnt,gh,nx,my,avemonfile)
      integer nt,newnt,md,newmd,nx,my,ifid
      real y(md,nt),gh(nx,my),avemon(nx,my)
      character*200 avemonfile
c
      nxmy4=nx*my*4
c      open(12,file=trim(eofano_ana2_file),access='direct',
c     &     form='unformatted',recl=nxmy4,status='unknown')
      open(36,file=trim(avemonfile),access='direct',
     &     form='unformatted',recl=nxmy4,status='old')
c
      do n=1,newnt
      nrec=n
      read(36,rec=nrec) avemon
c
      m=0
      do j=1,my
      do i=1,nx
        m=m+1
        gh(i,j)=y(m,n)+avemon(i,j)
      enddo
      enddo
c
c find the minimum value
c
      valmin=0.
      do j=1,my
      do i=1,nx
      if (gh(i,j).lt.valmin) then
      valmin=gh(i,j)
      endif
      enddo
      enddo
c
      print *,'the minumum pcp=',valmin
c
      do j=1,my
      do i=1,nx
      gh(i,j)=gh(i,j)-valmin
      enddo
      enddo
c
      print *,'write out record nrec= ',nrec
      print *,'y(20,2)= ',y(20,2)
      print *,'pcp(10,10)= ',gh(10,10)
      write(ifid,rec=nrec) gh
      enddo ! end loop of n, time series
c
      return
      end
c
c

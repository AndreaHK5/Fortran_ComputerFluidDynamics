c==================================================
           subroutine calccp1(ill,ifa,ind0,i1,i2,ivv,j1,j2,j3,ipf)
c------------------------------------------
c calculate interpolating coefficients using
c 1. third order Lagrangian interpolation scheme
c    for points in the physical cells
c 2. second order Lagrangian interpolation scheme
c    for points in ghost cells
c-----------------------------------------------------
       implicit double precision (a-h,o-z)
       include 'main.param'
c----------------------------------------
c ncx, ncy, ncz numero di celle fisiche del blocco
c ng numero di blocchi
      parameter(ncx = il-2, ncy = jl-2, ncz = kl-2)
c nx, ny, nz numero max di celle totali (fisiche + ghost) del blocco
      parameter(nx=ncx+2,ny=ncy+2,nz=ncz+2)
c nx2, ny2, nz2 numero di vertici del dominio orlato
      parameter(nx2 = 2*nx+1, ny2 = 2*ny+1, nz2 = 2*nz+1)

c  +++ global variables
c x2 coordinate dei vertici del dominio (nx2 dispari) o
c dei punti intermedi dei segmenti (nx2 pari)
c
       double precision x2(ng,nx2,ny2,nz2,3)
c dimensioni dei singoli blocchi
       integer lg,ilg(ng),jlg(ng),klg(ng)
       integer nilg(ng),njlg(ng),nklg(ng)
       integer iv(ng),jv(ng),kv(ng)
       integer ivert(ng),jvert(ng),kvert(ng)
       integer icord(ng)
       integer nby(ng),ngcouple(6,ng,3)
       double precision cp(31,ncm,ncm,3,6,ng)
c
c +++ local variables +++
       double precision xmin(ng),xmax(ng),ymin(ng),ymax(ng)
       double precision zmin(ng),zmax(ng)
       double precision xb(3),xxp(3)
       double precision xbb(27,3)
       double precision delx,dely,delz
       double precision a(3),b(3),c(3)
       double precision p1x,p2x,p3x,p1,p2,p3
       double precision q1x,q2x,q3x,q1,q2,q3
       double precision t1x,t2x,t3x,t1,t2,t3
       double precision xx(3),yy(3),zz(3),xi,yi,zi

      common /met1/x2,ilg,jlg,klg,lg
      common /met2/xix,xiy,xiz,etx,ety,etz,ztx,zty,ztz,
     & jac,g11,g22,g33
      common /met3/nby,ngcouple,cp,icord

      pi=4.*atan(1.)
c inizializzazione valori cp
      do i=1,27
      cp(i,i1,i2,ivv,ifa,ill)=0.
      enddo
      cp(28,i1,i2,ivv,ifa,ill)=j1
      cp(29,i1,i2,ivv,ifa,ill)=j2
      cp(30,i1,i2,ivv,ifa,ill)=j3
      cp(31,i1,i2,ivv,ifa,ill)=ipf
c punti di boundary numerica mista esterni al dominio di accoppiamento
      if(ipf.eq.0)return
c igh blocco da cui interpolare 
      igh=ipf
c====================================================== 
c coordinate del punto da interpolare (il punto cambia al cambiare di ivv
c====================================================== 
      if(ivv.eq.2)ip=0
      if(ivv.eq.1)ip=1
      if(ivv.eq.3)ip=2
c inversione ip se faccia 3 o 2
      if(ifa.eq.3.or.ifa.eq.2.or.ifa.eq.6)ip=-ip
      if(ifa.eq.1.or.ifa.eq.3)then
      xb(1)=x2(ill,2*ind0+ip,2*i1,2*i2,1)
      xb(2)=x2(ill,2*ind0+ip,2*i1,2*i2,2)
      xb(3)=x2(ill,2*ind0+ip,2*i1,2*i2,3)
      endif
      if(ifa.eq.2.or.ifa.eq.4)then
      xb(1)=x2(ill,2*i1,2*ind0+ip,2*i2,1)
      xb(2)=x2(ill,2*i1,2*ind0+ip,2*i2,2)
      xb(3)=x2(ill,2*i1,2*ind0+ip,2*i2,3)
      endif
      if(ifa.eq.5.or.ifa.eq.6)then
      xb(1)=x2(ill,2*i1,2*i2,2*ind0+ip,1)
      xb(2)=x2(ill,2*i1,2*i2,2*ind0+ip,2)
      xb(3)=x2(ill,2*i1,2*i2,2*ind0+ip,3)
      endif
c------------------------------------------------------
c conversione coordinate punto da interpolare se blocco cilindrico
c------------------------------------------------------
      if(icord(igh).eq.2)then
       xi=(xb(1)**2.+xb(2)**2.)**0.5
       yi=angle(xb(1),xb(2))
       else
c coordinate cartesiane
       xi=xb(1)
       yi=xb(2)
      endif
       zi=xb(3)
       xxp(1)=xi
       xxp(2)=yi
       xxp(3)=zi

c------------------------------------------------------
c calcolo dei coefficienti: punto interno
c------------------------------------------------------
c direzione x
c-----------------------
      if(j1.gt.1.and.j1.lt.ilg(igh))then
c cella contenente ghost del blocco ill e' cella interna del blocco igh
      do id=1,3
      xbb(13,id)=x2(igh,2*(j1-1),2*j2,2*j3,id)
      xbb(14,id)=x2(igh,2*j1,2*j2,2*j3,id)
      xbb(15,id)=x2(igh,2*(j1+1),2*j2,2*j3,id)
      enddo
c per blocco cilindrico, interpolazione usando coordinate cilindriche
      if(icord(igh).eq.2)then
       xx(1)=(xbb(13,1)**2.+xbb(13,2)**2.)**0.5
       xx(2)=(xbb(14,1)**2.+xbb(14,2)**2.)**0.5
       xx(3)=(xbb(15,1)**2.+xbb(15,2)**2.)**0.5
       else
c coordinate cartesiane
       xx(1)=xbb(13,1)
       xx(2)=xbb(14,1)
       xx(3)=xbb(15,1)
      endif
      call lagran3o(xx,xxp(1),a)
      else
c punto sul bordo indice min
        if(j1.eq.1)then
            do id=1,3
            xbb(14,id)=x2(igh,2*j1,2*j2,2*j3,id)
            xbb(15,id)=x2(igh,2*(j1+1),2*j2,2*j3,id)
            enddo
      if(icord(igh).eq.2)then
       xx(1)=(xbb(14,1)**2.+xbb(14,2)**2.)**0.5
       xx(2)=(xbb(15,1)**2.+xbb(15,2)**2.)**0.5
       xx(3)=1
       else
c coordinate cartesiane
       xx(1)=xbb(14,1)
       xx(2)=xbb(15,1)
       xx(3)=1
      endif
c punto sul bordo indice max
         else
c j1=ilg(igh)  
            do id=1,3
            xbb(13,id)=x2(igh,2*(j1-1),2*j2,2*j3,id)
            xbb(14,id)=x2(igh,2*j1,2*j2,2*j3,id)
            enddo
      if(icord(igh).eq.2)then
       xx(1)=(xbb(13,1)**2.+xbb(13,2)**2.)**0.5
       xx(2)=(xbb(14,1)**2.+xbb(14,2)**2.)**0.5
       xx(3)=3
       else
c coordinate cartesiane
       xx(1)=xbb(13,1)
       xx(2)=xbb(14,1)
       xx(3)=3
      endif
        endif
      call lagran2o(xx,xxp(1),a)
      endif
c fine direzione x
c------------------------------------------------------
c direzione y
c-----------------------
      if(j2.gt.1.and.j2.lt.jlg(igh))then
c cella contenente ghost del blocco ill e' cella interna del blocco igh
      do id=1,3
      xbb(11,id)=x2(igh,2*j1,2*(j2-1),2*j3,id)
      xbb(14,id)=x2(igh,2*j1,2*j2,2*j3,id)
      xbb(17,id)=x2(igh,2*j1,2*(j2+1),2*j3,id)
      enddo
c per blocco cilindrico, interpolazione usando coordinate cilindriche
      if(icord(igh).eq.2)then
       xx(1)=angle(xbb(11,1),xbb(11,2))
       xx(2)=angle(xbb(14,1),xbb(14,2))
       xx(3)=angle(xbb(17,1),xbb(17,2))
       else
c coordinate cartesiane
       xx(1)=xbb(11,2)
       xx(2)=xbb(14,2)
       xx(3)=xbb(17,2)
      endif
      call lagran3o(xx,xxp(2),b)
      else
c punto sul bordo indice min
        if(j2.eq.1)then
            do id=1,3
            xbb(14,id)=x2(igh,2*j1,2*j2,2*j3,id)
            xbb(17,id)=x2(igh,2*j1,2*(j2+1),2*j3,id)
            enddo
      if(icord(igh).eq.2)then
       xx(1)=angle(xbb(14,1),xbb(14,2))
       xx(2)=angle(xbb(17,1),xbb(17,2))
       xx(3)=1
       else
c coordinate cartesiane
       xx(1)=xbb(14,2)
       xx(2)=xbb(17,2)
       xx(3)=1
      endif
c punto sul bordo indice max
         else
c j2=jlg(igh) 
            do id=1,3
            xbb(11,id)=x2(igh,2*j1,2*(j2-1),2*j3,id)
            xbb(14,id)=x2(igh,2*j1,2*j2,2*j3,id)
            enddo
      if(icord(igh).eq.2)then
       xx(1)=angle(xbb(11,1),xbb(11,2))
       xx(2)=angle(xbb(14,1),xbb(14,2))
       xx(3)=3
       else
c coordinate cartesiane
       xx(1)=xbb(11,2)
       xx(2)=xbb(14,2)
       xx(3)=3
      endif
        endif
      call lagran2o(xx,xxp(2),b)
      endif
c fine direzione y
c------------------------------------------------------
c direzione z
c-----------------------
      if(j3.gt.1.and.j3.lt.klg(igh))then
c cella contenente ghost del blocco ill e' cella interna del blocco igh
      do id=1,3
      xbb(5,id)=x2(igh,2*j1,2*j2,2*(j3-1),id)
      xbb(14,id)=x2(igh,2*j1,2*j2,2*j3,id)
      xbb(23,id)=x2(igh,2*j1,2*j2,2*(j3+1),id)
      enddo
c coordinate non cambiano con tipo di blocco lungo z
       xx(1)=xbb(5,3)
       xx(2)=xbb(14,3)
       xx(3)=xbb(23,3)
      call lagran3o(xx,xxp(3),c)
      else
c punto sul bordo indice min
        if(j3.eq.1)then
            do id=1,3
            xbb(14,id)=x2(igh,2*j1,2*j2,2*j3,id)
            xbb(23,id)=x2(igh,2*j1,2*j2,2*(j3+1),id)
            enddo
       xx(1)=xbb(14,3)
       xx(2)=xbb(23,3)
       xx(3)=1
c punto sul bordo indice max
         else
c j2=jlg(igh)
            do id=1,3
            xbb(5,id)=x2(igh,2*j1,2*j2,2*(j3-1),id)
            xbb(14,id)=x2(igh,2*j1,2*j2,2*j3,id)
            enddo
c coordinate cartesiane
       xx(1)=xbb(5,2)
       xx(2)=xbb(14,2)
       xx(3)=3
        endif
      call lagran2o(xx,xxp(3),c)
      endif
c fine direzione z
c-------------------------------------------------
c assegnazione cp
c----------------------------------------------
      do idd=1,3
         do ipo=1,3
           do ipp=1,3
           cp(ipo+(idd-1)*3+(ipp-1)*9,i1,i2,ivv,ifa,ill)=
     &       a(ipo)*b(idd)*c(ipp)
            write(230,*)cp(ipo+(idd-1)*3+(ipp-1)*9,i1,i2,ivv,ifa,ill)
           enddo
         enddo
      enddo
      write(230,*)cp(28,i1,i2,ivv,ifa,ill)
      write(230,*)cp(29,i1,i2,ivv,ifa,ill)
      write(230,*)cp(30,i1,i2,ivv,ifa,ill)
      write(230,*)cp(31,i1,i2,ivv,ifa,ill)
      return
      end






c==========================================
       subroutine lagran2o(xx,xi,a)
c==========================================
c subroutine che calcola il contributo dei tre punti xx nella
c posizione xxp, valori messi in a
c==========================================
c interpolatore Lagrangiano terzo ordine, griglia non uniforme
c==========================================
      implicit double precision (a-h,o-z)
      double precision a(3),xx(3),xi
      p1x=xi-xx(2)
      p2x=xi-xx(1)
      p1=xx(1)-xx(2)

c calcolo coefficienti lungo x
      if(p1.eq.0)then
        print *,'p1', p1,j1,j2,j3
        print *,xx(1),xx(2)
        stop
      endif
      aa1=p1x/p1
      aa2=-p2x/p1
      if(xx(3).eq.1)then
        a(1)=0.
        a(2)=aa1
        a(3)=aa2
       else
        a(1)=aa1
        a(2)=aa2
        a(3)=0.
      endif
      asum=a(1)+a(2)+a(3)
      if(abs(asum-1.).gt.0.01)then
      print *,(a(i),i=1,3)
      print *,'valori coefficiente interpolazione scorretti 3 points' 
      stop
      endif
      return
      end
c------------------------------------------------




c==========================================
       subroutine lagran3o(xx,xi,a)
c==========================================
c subroutine che calcola il contributo dei tre punti xx nella
c posizione xxp, valori messi in a
c==========================================
c interpolatore Lagrangiano terzo ordine, griglia non uniforme
c==========================================
      implicit double precision (a-h,o-z)
      double precision a(3),xx(3),xi
      p1x=xi-xx(2)
      p2x=xi-xx(3)
      p3x=xi-xx(1)
      p1=xx(1)-xx(2)
      p2=xx(2)-xx(3)
      p3=xx(3)-xx(1)

c calcolo coefficienti lungo x
      if(p1.eq.0.or.p2.eq.0.or.p3.eq.0)then
        print *,'p1,p2,p3', p1,p2,p3,j1,j2,j3
        print *,xx(1),xx(2),xx(3)
        stop
      endif
      a(1)=-p1x*p2x/(p1*p3)
      a(2)=-p2x*p3x/(p2*p1)
      a(3)=-p3x*p1x/(p2*p3)
      asum=a(1)+a(2)+a(3)
      if(abs(asum-1.).gt.0.01)then
      print *,(a(i),i=1,3)
      print *,'valori coefficiente interpolazione scorretti 2 points'
      stop
      endif
      return
      end

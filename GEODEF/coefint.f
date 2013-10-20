      subroutine coefint
c------------------------------------------
c
c   subroutine to get interpolation coefficient for multidomain
c
c------------------------------------------
c
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
       integer icord(ng),ijm(ng)
       integer nby(ng),ngcouple(6,ng,3)
       double precision cp(31,ncm,ncm,3,6,ng)
c
c +++ local variables +++
       double precision xmin(ng),xmax(ng),ymin(ng),ymax(ng)
       double precision zmin(ng),zmax(ng)
       double precision coor(3,ncx+1,ncy+1,ncz+1)
       double precision delx,dely,delz

      common /met1/x2,ilg,jlg,klg,lg
      common /met2/xix,xiy,xiz,etx,ety,etz,ztx,zty,ztz,
     & jac,g11,g22,g33
      common /met3/nby,ngcouple,cp,icord

      pi=4.*atan(1.)
c---------------------------------------------------------------------
c
c initialization nby, ngcouple
c
c---------------------------------------------------------------------
      do ig=1,ng
      nby(ig)=0
        do iface=1,6
          do ipo=1,3
        ngcouple(iface,ig,ipo)=0
          enddo
        enddo
      enddo

c---------------------------------------------------------------------
c
c Input data: da file cmin informazioni su accoppiamenti
c----------------------------------
c
      read(10,*)
      do ig=1,lg
c nby(ip) = numero di boundary numeriche
      read(10,*)
      read(10,*) ip,nby(ip)
      print *,ip,nby(ip)
      if(nby(ip).eq.0)goto 11
      read(10,*)
        do iface=1,nby(ip)
        read(10,*) if,(ngcouple(if,ip,ipo),ipo=1,3)
        enddo
 11   continue
      enddo
c----------------------------------------------
c ciclo sui blocchi
      do 1,ill=1,lg
        if(nby(ill).eq.0)goto 1
c se sono definiti accoppiamenti numerici, cerca la faccia
        do 2 ifa=1,6
          if(ngcouple(ifa,ill,1).eq.0)goto 2
          print *,'doing blocco',ill,'faccia',ifa
c direzione x
          if(ifa.eq.1.or.ifa.eq.3)then
           ind1=jlg(ill)-1
           ind2=klg(ill)-1
c assegnazione indice i cella ghost da interpolare
           if(ifa.eq.1)then
             ind0=1
           else
             ind0=ilg(ill)-1
           endif
          endif
c direzione y
          if(ifa.eq.2.or.ifa.eq.4)then
           ind1=ilg(ill)-1
           ind2=klg(ill)-1
           if(ifa.eq.4)then
             ind0=1
           else
             ind0=jlg(ill)-1
           endif
          endif
c direzione z
          if(ifa.eq.5.or.ifa.eq.6)then
           ind1=ilg(ill)-1
           ind2=jlg(ill)-1
           if(ifa.eq.5)then
             ind0=1
           else
             ind0=klg(ill)-1
           endif
          endif
           do i1=1,ind1
            do i2=1,ind2
             do ivv=1,3
c call ripetuta per punti adiacenti al ghost usati per il calcolo
c delle derivate
            do 10,ipo=1,3
            ipf=ngcouple(ifa,ill,ipo)
            if(ipf.eq.0)goto 10
           call locate(ill,ifa,ind0,i1,i2,ivv,j1,j2,j3,ipf)
            if(j1.ne.0.and.j2.ne.0.and.j3.ne.0)goto 12
 10         continue
            ipf=0
 12         continue
           call calccp1(ill,ifa,ind0,i1,i2,ivv,j1,j2,j3,ipf)
             enddo
            enddo
           enddo
 2       continue
 1       continue
c---------------------------------------------------
c scrittura file coef.in
c ---------------------------------------------------
          do ill=1,lg
          ijm(ill)=max0(ilg(ill)-1,jlg(ill)-1,klg(ill)-1)
          write(6,*)'blocco',ill,'max dim',ijm(ill),
     &     'numero bnum',nby(ill)
          if(nby(ill).ne.0)then
          write(22)(((((cp(i,j,k,ip,nyy,ill),i=1,31),j=1,ijm(ill)),
     &     k=1,ijm(ill)),ip=1,3),nyy=1,nby(ill))      
c           write(229,*)'blocco',ill,'nbnum',nby(ill)
c          write(229,800)(((((cp(i,j,k,ip,nyy,ill),i=1,31),j=1,ijm(ill)),
c     &     k=1,ijm(ill)),ip=1,3),nyy=1,nby(ill))      
c 800      format(f10.5)
          endif
          enddo
         return
         end


         subroutine locate(ill,ifa,ind0,i1,i2,ivv,index1,
     & index2,index3,ipf)
c------------------------------------------
c
c   subroutine to localize point of the boundary in the interior of
c   the coupled domain
c
c   variables: ill, domain number
c              ifa, face number
c              ipo, first, second or third block
c              ind0, value of the grid fixed on the boundary 
c              i1, first index in boundary face
c              i2, second index in boundary face
c              ivv, index for point position
c              index1, first direction index for location in the 2nd domain
c              index2, second direction index for location in the 2nd domain
c              index3, third direction index for location in the 2nd domain
c j1, j2 and j3 are assigned to cp(28,,,,), cp(29,,,,) and cp(30,,,,,)
c cp(31,,,,) contains the number of the block coupled to the point
c------------------------------------------
c
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
       double precision xbound(3)
       double precision delx,dely,delz

      common /met1/x2,ilg,jlg,klg,lg
      common /met2/xix,xiy,xiz,etx,ety,etz,ztx,zty,ztz,
     & jac,g11,g22,g33
      common /met3/nby,ngcouple,cp,icord

      pi=4.*atan(1.)
      index1=0
      index2=0
      index3=0
c assegnazione indici per localizzazione punto in base a ivv:
c ivv=2 centro cella ghost
c ivv=1 punto sulla faccia di boundary
c ivv=3 primo punto interno

      if(ivv.eq.2)ip=0
      if(ivv.eq.1)ip=1
      if(ivv.eq.3)ip=2
c inversione indice se ultima cella ghost: calcolo va fatto per celle interne
      if(ifa.eq.3.or.ifa.eq.2.or.ifa.eq.6)ip=-ip
c assegnazione coordinata
        if(ifa.eq.1.or.ifa.eq.3)then   
      xbound(1)=x2(ill,2*ind0+ip,2*i1,2*i2,1)
      xbound(2)=x2(ill,2*ind0+ip,2*i1,2*i2,2)
      xbound(3)=x2(ill,2*ind0+ip,2*i1,2*i2,3)
        endif
        if(ifa.eq.2.or.ifa.eq.4)then   
      xbound(1)=x2(ill,2*i1,2*ind0+ip,2*i2,1)
      xbound(2)=x2(ill,2*i1,2*ind0+ip,2*i2,2)
      xbound(3)=x2(ill,2*i1,2*ind0+ip,2*i2,3)
        endif
        if(ifa.eq.5.or.ifa.eq.6)then   
      xbound(1)=x2(ill,2*i1,2*i2,2*ind0+ip,1)
      xbound(2)=x2(ill,2*i1,2*i2,2*ind0+ip,2)
      xbound(3)=x2(ill,2*i1,2*i2,2*ind0+ip,3)
        endif

c igh blocco accoppiato alla faccia ifa del blocco ill
c se esistono piu' blocchi accoppiati alla stessa faccia, 
c call a locate e' ripetuta per ogni blocco possibile
      igh=ipf
c per blocco igh definito cilindrico, cambio di variabili
       if(icord(igh).eq.2)then
         varx=(xbound(1)**2.+xbound(2)**2.)**0.5
         vary=angle(xbound(1),xbound(2))
        else
        varx=xbound(1)
        vary=xbound(2)
       endif
        varz=xbound(3)
c confronto coordinate: cartesiane o cilindriche in base al blocco in
c cui si cerca
c cerca solo su celle fisiche, perche' serve una cella in piu' da ogni lato
c per interpolazione
      do i=1,ilg(igh)-1
        i22=2*i-1
        i3=2*i+1
c assumo distribuzione omogenea nelle direzioni all'interno dei blocchi
c primo indice= coordinata radiale
        if(icord(igh).eq.2)then
          hvarx1=(x2(igh,i22,2,2,1)**2.+x2(igh,i22,2,2,2)**2.)**0.5 
          hvarx2=(x2(igh,i3,2,2,1)**2.+x2(igh,i3,2,2,2)**2.)**0.5 
        else
c o coordinata x
          hvarx1=x2(igh,i22,2,2,1) 
          hvarx2=x2(igh,i3,2,2,1) 
        endif
        if(varx.ge.hvarx1.and.varx.le. hvarx2)then
         index1=i
         goto 2
        endif
      enddo
      print *,'index i not found'
      write(888,*)ill,ifa,ind0,i1,i2,ivv
      return
 2    continue
        do j=1,jlg(igh)-1
        j22=2*j-1
        j3=2*j+1
c secondo indice: coordinata theta
        if(icord(igh).eq.2)then
          hvary1=angle(x2(igh,2*index1,j22,2,1),
     &        x2(igh,2*index1,j22,2,2))
          hvary2=angle(x2(igh,2*index1,j3,2,1),x2(igh,2*index1,j3,2,2))
          if(hvary1.gt.hvary2)hvary2=hvary2+2.*pi
        else
c o coordinata y 
          hvary1=x2(igh,2*index1,j22,2,2)
          hvary2=x2(igh,2*index1,j3,2,2)
        endif
        if(vary.gt.hvary1.and.vary.le. hvary2)then
         index2=j
         goto 3
        endif
      enddo
      print *,'index j not found'
      return
 3    continue
c accoppiamento su piano z-normale
        do k=1,klg(igh)-1
        k22=2*k-1
        k3=2*k+1
c secondo indice: coordinata theta
          hvarz1=x2(igh,2*index1,index2,k22,3)
          hvarz2=x2(igh,2*index1,index2,k3,3)
        if(varz.gt.hvarz1.and.varz.le. hvarz2)then
         index3=k
         goto 4
        endif
        enddo
      print *,'index k not found'
      return
 4    continue

      return
      end

c-----------------------------------------
         function angle(xx,yy)
c-----------------------------------------
         implicit double precision (a-h,o-z)
         pi=4.*atan(1.)
         if(xx.eq.0)then
           angle=pi/2.
           if(yy.lt.0) angle=3.*pi/2.
         else
         angle=atan(yy/xx)
         if(yy.gt.0.and.xx.lt.0)angle=angle+pi
         if(yy.lt.0.and.xx.lt.0)angle=angle+pi
         if(yy.lt.0.and.xx.gt.0)angle=angle+2*pi
         endif 
         end




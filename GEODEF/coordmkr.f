c--------------------------------------------------
c                a Gino's production
c
c   
c   subroutine!!!!! calcolo x2
c-----------------------------------------------------

      subroutine coordmkr
      implicit double precision (a-h,o-z)     
      include 'main.param'
      parameter( nmaxx=il-2, nmaxy=jl-2, nmaxz=kl-2)
      parameter( nx2=2*il+1, ny2=2*jl+1, nz2=2*kl+1 )
      integer ilg(ng), jlg(ng), klg(ng),lg
      integer iv(ng), jv (ng), kv (ng)
      integer nilg(ng),njlg(ng),nklg(ng),icord(ng)
      double precision coor(3,nmaxx+1,nmaxy+1,nmaxz+1)
      double precision x2(ng,nx2,ny2,nz2,3)   
      double precision xmin(ng),xmax(ng)
      double precision ymin(ng),ymax(ng)
      double precision zmin(ng),zmax(ng)
      integer nx(2,ng),ny(2,ng),nz(2,ng)

c   main.param contiene i valori massimi per il problema proposto
c  ng numero di blocchi
c  nmaxx,nmaxy,nmaxz max assignable cells in x,y,z direction
c  ilg,jlg,klg numero nodi (ill)
c  iv,jv,kv numero totale di vertici all inclusive (vertici+nodi di celle normali e
c     ghost)
c  nx2 numero massimo di vertici su x      
c  Attenzione, coor e x2 hanno definizioni diverse!!!
      
      common /met1/x2,iv,jv,kv,lg

c in questa common ci sono:
c x2 che è il matricione della griglia
c iv,jv,kv che sono i numeri di vertici sugli assi   
c ilg,jlg,klg che sono il numero di faccie sugli assi      
c
c  la common met 2 nn è stata dichiarata qui in quanto contiene i paramentri per
c  la controvarianza che qui non servono 
      
      open(unit=10, file='inputgr')
      
c  inputgr ha i paramentri richiesti dall'utente       
 
 

c ------------------
c lettura inputgr e controlli di coerenza
c ----------------
  
      read(10,*)
      read(10,*) lg
      if (lg.gt.ng)then 
        print *,'to many blocks, recopile with greater ng in main.param'
	print *,'current max block number ',ng
	stop
      endif
      do 1,ill=1,lg
       read(10,*)
       read(10,*) xmin(ill),xmax(ill)
       read(10,*) ymin(ill),ymax(ill)
       read(10,*) zmin(ill),zmax(ill)
       read(10,*)
       read(10,*) nx(2,ill),nx(1,ill)
       nilg(ill)=nx(1,ill)+nx(2,ill)
       if (nilg(ill).ge.nmaxx)then
         print *,'too many asked cells in x directions'
	 print *,'current max cell number in x direction',nmaxx
	 stop
       endif
       
       read(10,*)
       read(10,*) ny(2,ill),ny(1,ill)
       njlg(ill)=ny(1,ill)+ny(2,ill)       
       if (njlg(ill).ge.nmaxy)then
         print *,'too many asked cells in y directions'
	 print *,'current max cell number in y direction',nmaxy
	 stop
       endif
      
       read(10,*)
       read(10,*) nz(2,ill),nz(1,ill)
        nklg(ill)=nz(1,ill)+nz(2,ill)
      if (nklg(ill).ge.nmaxz)then
         print *,'too many asked cells in z directions'
	 print *,'current max cell number in z direction',nmaxz
	 stop
       endif
      read(10,*)
      read(10,*)
 1    continue


c-----------------------
c fine lettura inputgr
c inizio calcolo x2
c-----------------------

       do 2,ill=1,lg
         
	 icord(ill)=1
	 
	 delx=(xmax(ill)-xmin(ill))/(1.*nilg(ill))
         dely=(ymax(ill)-ymin(ill))/(1.*njlg(ill))
         delz=(zmax(ill)-zmin(ill))/(1.*nklg(ill))

         ilg(ill)=nilg(ill)+1
         jlg(ill)=njlg(ill)+1
         klg(ill)=nklg(ill)+1
c----------------
c calcolo coor, solo vertici
c Questi ilg,jlg,klg in alto sono il numero dei vertici delle celle
c------------------------
         
	 print *,''
	 print *,'Phisical dimension of block #',ill
	 print *,xmin(ill),' ->',xmax(ill),' on x'
	 print *,ymin(ill),' ->',ymax(ill),' on y'
	 print *,zmin(ill),' ->',zmax(ill),' on z'
	 print *,'Phisical cells in block #',ill,' =',
     & nilg(ill),njlg(ill),nklg(ill)

c--------------------------
c coordinate cartesiane
c-----------------------
         if(icord(ill).eq.1)then
         do 3,i=1,ilg(ill)
         do 3,j=1,jlg(ill)
         do 3,k=1,klg(ill)
         coor(1,i,j,k)=xmin(ill)+dfloat(i-1)*delx
         coor(2,i,j,k)=ymin(ill)+dfloat(j-1)*dely
         coor(3,i,j,k)=zmin(ill)+dfloat(k-1)*delz
 3       continue
         else
c--------------------------
c coordinate cilindriche: r, theta, z
c-----------------------
         do 33,i=1,ilg(ill)
         do 33,j=1,jlg(ill)
         do 33,k=1,klg(ill)
         erre=xmin(ill)+(i-1)*delx
         theta=ymin(ill)+(j-1)*dely
         thetarad=theta*pi/180.
         coor(1,i,j,k)=erre*cos(thetarad)
         coor(2,i,j,k)=erre*sin(thetarad)
         coor(3,i,j,k)=zmin(ill)+(k-1)*delz
 33       continue
         endif
c --------------         
c  copio coor in x2
c  ma prima trovo le dimensioni di x2 per il caso proposto
c ---------------
      iv(ill)=2*(nilg(ill)+2)+1
      jv(ill)=2*(njlg(ill)+2)+1
      kv(ill)=2*(nklg(ill)+2)+1
      print *,'Total vertices in block #',ill,' =',
     & iv(ill),jv(ill),kv(ill)
      print *,''
            
      do 4,i=1,ilg(ill)
      do 4,j=1,jlg(ill)
      do 4,k=1,klg(ill)
        i2=2*i+1
        j2=2*j+1
        k2=2*k+1
        x2(ill,i2,j2,k2,1)=coor(1,i,j,k)
        x2(ill,i2,j2,k2,2)=coor(2,i,j,k)
        x2(ill,i2,j2,k2,3)=coor(3,i,j,k)
4     continue      
      
c -----------------
c inziamo ad orlare con le ghost su x
c --------------------

         do 5,j=1,jlg(ill)
         do 5,k=1,klg(ill)
         i2=1
         j2=2*j+1
         k2=2*k+1
	 x2(ill,i2,j2,k2,1)=2.*x2(ill,i2+2,j2,k2,1)-x2(ill,i2+4,j2,k2,1)
         x2(ill,i2,j2,k2,2)=2.*x2(ill,i2+2,j2,k2,2)-x2(ill,i2+4,j2,k2,2)
         x2(ill,i2,j2,k2,3)=2.*x2(ill,i2+2,j2,k2,3)-x2(ill,i2+4,j2,k2,3)
         i2=iv(ill)
         x2(ill,i2,j2,k2,1)=2.*x2(ill,i2-2,j2,k2,1)-x2(ill,i2-4,j2,k2,1)
         x2(ill,i2,j2,k2,2)=2.*x2(ill,i2-2,j2,k2,2)-x2(ill,i2-4,j2,k2,2)
         x2(ill,i2,j2,k2,3)=2.*x2(ill,i2-2,j2,k2,3)-x2(ill,i2-4,j2,k2,3)
5       continue
         
c -----------------
c   ghost su y
c   il contatore i ha gli estremi aggiornati ai vetrici appena aggiunti
c --------------------
 
	 do 6,i=0,ilg(ill)+1
         do 6,k=1,klg(ill)
         j2=1
         i2=2*i+1
         k2=2*k+1
         x2(ill,i2,j2,k2,1)=2.*x2(ill,i2,j2+2,k2,1)-x2(ill,i2,j2+4,k2,1)
         x2(ill,i2,j2,k2,2)=2.*x2(ill,i2,j2+2,k2,2)-x2(ill,i2,j2+4,k2,2)
         x2(ill,i2,j2,k2,3)=2.*x2(ill,i2,j2+2,k2,3)-x2(ill,i2,j2+4,k2,3)
         j2=jv(ill)
         x2(ill,i2,j2,k2,1)=2.*x2(ill,i2,j2-2,k2,1)-x2(ill,i2,j2-4,k2,1)
         x2(ill,i2,j2,k2,2)=2.*x2(ill,i2,j2-2,k2,2)-x2(ill,i2,j2-4,k2,2)
         x2(ill,i2,j2,k2,3)=2.*x2(ill,i2,j2-2,k2,3)-x2(ill,i2,j2-4,k2,3)
 6       continue
c --------------------
c ghost direzione z
c i contatori i e j sono aggiornati e permettono di completare 
c   la griglia delle phantom coi punti mancanti
c --------------------
         do 7,j=0,jlg(ill)+1
         do 7,i=0,ilg(ill)+1
         k2=1
         j2=2*j+1
         i2=2*i+1
         x2(ill,i2,j2,k2,1)=2.*x2(ill,i2,j2,k2+2,1)-x2(ill,i2,j2,k2+4,1)
         x2(ill,i2,j2,k2,2)=2.*x2(ill,i2,j2,k2+2,2)-x2(ill,i2,j2,k2+4,2)
         x2(ill,i2,j2,k2,3)=2.*x2(ill,i2,j2,k2+2,3)-x2(ill,i2,j2,k2+4,3)
         k2=kv(ill)
         x2(ill,i2,j2,k2,1)=2.*x2(ill,i2,j2,k2-2,1)-x2(ill,i2,j2,k2-4,1)
         x2(ill,i2,j2,k2,2)=2.*x2(ill,i2,j2,k2-2,2)-x2(ill,i2,j2,k2-4,2)
         x2(ill,i2,j2,k2,3)=2.*x2(ill,i2,j2,k2-2,3)-x2(ill,i2,j2,k2-4,3)
7       continue

c --------------
c aggiungiamo i centri segmenti
c prima lungo x (punti ad indice i pari)
c ---------------

      do 8 k=1,kv(ill),2
      do 8 j=1,jv(ill),2
      do 8 i=2,iv(ill),2
        x2(ill,i,j,k,1)=0.5*(x2(ill,i-1,j,k,1)+
     &                              x2(ill,i+1,j,k,1))
        x2(ill,i,j,k,2)=0.5*(x2(ill,i-1,j,k,2)+
     &                             x2(ill,i+1,j,k,2))
        x2(ill,i,j,k,3)=0.5*(x2(ill,i-1,j,k,3)+
     &                              x2(ill,i+1,j,k,3))
8     continue

c----------------------       
c  ora su y (indice j pari)
c -----------------       
       
       do 18 k=1,kv(ill),2
       do 18 j=2,jv(ill),2
       do 18 i=1,iv(ill),2
         x2(ill,i,j,k,1)=0.5*(x2(ill,i,j-1,k,1)+
     &                              x2(ill,i,j+1,k,1))
         x2(ill,i,j,k,2)=0.5*(x2(ill,i,j-1,k,2)+
     &                              x2(ill,i,j+1,k,2))
         x2(ill,i,j,k,3)=0.5*(x2(ill,i,j-1,k,3)+
     &                              x2(ill,i,j+1,k,3))
18     continue
       
c ----------------       
c  e infine su z Indice k pari)
c  -----------------------       
      do 28 k=2,kv(ill),2
      do 28 j=1,jv(ill),2
      do 28 i=1,iv(ill),2
        x2(ill,i,j,k,1)=0.5*(x2(ill,i,j,k-1,1)+
     &                              x2(ill,i,j,k+1,1))
         x2(ill,i,j,k,2)=0.5*(x2(ill,i,j,k-1,2)+
     &                              x2(ill,i,j,k+1,2))
         x2(ill,i,j,k,3)=0.5*(x2(ill,i,j,k-1,3)+
     &                              x2(ill,i,j,k+1,3))
28     continue

c --------------------------
c calcolo centri faccie
c facce orizzontali
c ---------------------------

      do 9 k=1,kv(ill),2
      do 9 j=2,jv(ill),2
      do 9 i=2,iv(ill),2
        x2(ill,i,j,k,1)=0.5*(x2(ill,i,j-1,k,1) +
     &                                x2(ill,i,j+1,k,1))
        x2(ill,i,j,k,2)=0.5*(x2(ill,i,j-1,k,2) +
     &                                x2(ill,i,j+1,k,2))
        x2(ill,i,j,k,3)=0.5*(x2(ill,i,j-1,k,3) +
     &                                x2(ill,i,j+1,k,3))
9       continue

c -----------------      
c  facce parallele al laterale
c---------------
      
      do 19 k=2,kv(ill),2
      do 19 j=1,jv(ill),2
      do 19 i=2,iv(ill),2
        x2(ill,i,j,k,1)=0.5*(x2(ill,i+1,j,k,1) +
     &                                x2(ill,i-1,j,k,1))
        x2(ill,i,j,k,2)=0.5*(x2(ill,i+1,j,k,2) +
     &                                x2(ill,i-1,j,k,2))
        x2(ill,i,j,k,3)=0.5*(x2(ill,i+1,j,k,3) +
     &                                x2(ill,i-1,j,k,3))
19       continue

c -------------
c  facce parallele all'ingresso
c -----------------

      do 29 k=2,kv(ill),2
      do 29 j=2,jv(ill),2
      do 29 i=1,iv(ill),2
        x2(ill,i,j,k,1)=0.5*(x2(ill,i,j,k-1,1) +
     &                                x2(ill,i,j,k+1,1))
        x2(ill,i,j,k,2)=0.5*(x2(ill,i,j,k-1,2) +
     &                                x2(ill,i,j,k+1,2))
        x2(ill,i,j,k,3)=0.5*(x2(ill,i,j,k-1,3) +
     &                                x2(ill,i,j,k+1,3))
29     continue

c --------------
c  ed infine i centri cella
c ------------------

      do 10 k=2,kv(ill),2
      do 10 j=2,jv(ill),2
      do 10 i=2,iv(ill),2
        x2(ill,i,j,k,1)=0.5*(x2(ill,i,j,k-1,1) +
     &                                  x2(ill,i,j,k+1,1))
        x2(ill,i,j,k,2)=0.5*(x2(ill,i,j,k-1,2) +
     &                                  x2(ill,i,j,k+1,2))
        x2(ill,i,j,k,3)=0.5*(x2(ill,i,j,k-1,3) +
     &                                  x2(ill,i,j,k+1,3))
10     continue
2       continue  
      return
         end














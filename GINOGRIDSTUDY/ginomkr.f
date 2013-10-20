c--------------------------------------------------
c                a Gino's production
c
c   
c   subroutine!!!!! calcolo x2
c-----------------------------------------------------

      subroutine ginomkr 
      implicit double precision (a-h,o-z)     
      include 'main.param'
      parameter( nmaxx=il-2, nmaxy=jl-2, nmaxz=kl-2)
      parameter( nx2=2*il+1, ny2=2*jl+1, nz2=2*kl+1 )
      double precision xmin(ng),xmax(ng)
      double precision ymin(ng),ymax(ng)
      double precision zmin(ng),zmax(ng)
      integer ill,sg
      integer nx(2,ng),ny(2,ng),nz(2,ng)
      integer nxt,nyt,nzt,nxs,nys,nzs
      double precision Dx,Dy,Dz
      double precision x2(ng,nx2,ny2,nz2,3)
      double precision coor(3,nmaxx+1,nmaxy+1,nmaxz+1)
      double precision gx(nmaxx+1),gy(nmaxy+1),gz(nmaxz+1)
      double precision lxtemp,lytemp,lztemp,ax,ay,az
      integer axr(ng),ayr(ng),azr(ng)
      integer iv(ng), jv (ng), kv (ng)
      integer ilg(ng), jlg(ng), klg(ng)
       
       double precision DDX
      
      
c --------------------------
c main.patram è il nostro glorioso sponson e nn necessita presentazioni
c ill è un altro nostro cliente di vecchia data
c sg è l'idice di sottoblocco (per fare la gino griglia)
c nx, ny, nz sono i numeri di celle MA sono sul semiasse negativo (2,ill)
c    e positivo (1,ill) 
c nxy ed nxs sono delle variabili temporanee per le iterazioni in sottoblocco
c Dx, Dy, Dz serve dire??
c x2 è l'output di questa subroutine, along with iv,jv,kv,lg
c  gx(i) contiene nelle iterazioni di sottoblocco la divisione del semiasse
c  lxtemp-ax sono all'interno dei sottoblocchi per creare gx
c  axr è il ginoparam della situazione
c --------------------------------  

      common /met1/x2,iv,jv,kv,lg
     
c-------------
c in questa common ci sono:
c x2 che è il matricione della griglia
c iv,jv,kv che sono i numeri di vertici sugli assi   
c ilg,jlg,klg che sono il numero di faccie sugli assi      
c
c  la common met 2 nn è stata dichiarata qui in quanto contiene i paramentri per
c  la controvarianza che qui non servono 
c -------------    
     
      
      open(unit=10, file='inputgr')
      open(unit=41, file='usrsegm')
c ----------------------  
c  inputgr contiene i parametri in lettura
c  phisicell ha solo la cella fisica
c ----------------------

c ----------------
c lettura inputgr e controlo di coerenza
c ----------------
     
      read(10,*)
      read(10,*) lg
      if (lg.gt.ng)then 
        print *,'to many blocks, recopile with greater ng in main.param'
	print *,'current max block number ',ng
	stop
      endif
      
      do 1 ill=1,lg
      read(10,*)
      read(10,*) xmin(ill),xmax(ill)
      read(10,*) ymin(ill),ymax(ill)
      read(10,*) zmin(ill),zmax(ill)
      if (xmin(ill).gt.0.0.or. xmax(ill).lt.0.0.or.
     & ymin(ill).gt.0.0.or. ymax(ill).lt.0.0.or.
     & zmin(ill).gt.0.0.or. zmax(ill).lt.0.0) then
        print *,''
	print *,'error in inputgr definition, bad asked dimension(s)'
	print *,'Cant make a "Ginogrid"(pat.pend)!'
	print *,''
	stop
      endif
      read (10,*)
      read (10,*) nx(2,ill),nx(1,ill)
      read (10,*)
      read (10,*) ny(2,ill),ny(1,ill)
      read (10,*)
      read (10,*) nz(2,ill),nz(1,ill)
      read (10,*)
      read (10,*) axr(ill),ayr(ill),azr(ill)
 

        if ((nx(2,ill)+nx(1,ill)).ge.nmaxx)then
         print *,'too many asked cells in x directions'
	 print *,'current max cell number in x direction',nmaxx
	 print *,'Cant make a "Ginogrid"(pat.pend)!'
	 print *,''
	 stop
       endif
        if ((ny(2,ill)+ny(1,ill)).ge.nmaxy)then
         print *,'too many asked cells in y directions'
	 print *,'current max cell number in y direction',nmaxy
	 print *,'Cant make a "Ginogrid"(pat.pend)!'
	 print *,''
	 stop
       endif
      if ((nz(2,ill)+nz(1,ill)).ge.nmaxz)then
         print *,'too many asked cells in z directions'
	 print *,'current max cell number in z direction',nmaxz
	 print *,'Cant make a "Ginogrid"(pat.pend)!'
	 print *,''
	 stop
       endif
 
 1     continue
      
      
      do 2 ill=1,lg     
        
	 print *,''
	 print *,'Phisical dimension of block #',ill
	 print *,xmin(ill),' ->',xmax(ill),' on x'
	 print *,ymin(ill),' ->',ymax(ill),' on y'
	 print *,zmin(ill),' ->',zmax(ill),' on z'
	 print *,'Phisical cells in block #',ill,' =',
     & (nx(1,ill)+nx(2,ill)),(ny(1,ill)+ny(2,ill)),
     & (nz(1,ill)+nz(2,ill))

    
c ------------------------
c controlli su correttezza numeri di griglie su semiassi negativi
c ------------------------      
      
      if (xmin(ill).eq.0.0.and.nx(2,ill).ne.0)then
         print *,'bad-x-neg-cell-number value, changed to zero!'
	 print *,'change to zero also in inputgr or initial wont work'
	 print *,''
	 nx(2,ill)=0
      endif
       
      if (ymin(ill).eq.0.0.and.ny(2,ill).ne.0)then
         print *,'bad-y-neg-cell-number value, changed to zero!'
	 print *,'change to zero also in inputgr or initial wont work'
	 print *,''
	 ny(2,ill)=0
      endif
        
      if (zmin(ill).eq.0.0.and.nz(2,ill).ne.0)then
	print *,'bad z-neg-cell-number value, changed to zero!'
	print *,'change to zero also in inputgr or initial wont work'
	print *,''
	nz(2,ill)=0
      endif

c ------------------
c iterazione per gli otto blocchi
c  acquisizione dati, calcolo delle divisioni delgi assi
c  scrittura a mitraglia su chi di dovere 
c
c ----------------- 
      	
	
	ax = axr(ill)/100.
	ay = ayr(ill)/100.
	az = azr(ill)/100.
      
      do 3 sg=1,8
      
        
      if (sg.eq.1) then	  

	  nxt=nx(1,ill)
	  nyt=ny(1,ill)
	  nzt=nz(1,ill)
	  
	  nxs=nx(2,ill)
	  nys=ny(2,ill)
	  nzs=nz(2,ill)
	  
	  Dx = (xmax(ill))/(1.*nxt)
	  Dy = (ymax(ill))/(1.*nyt)
	  Dz = (zmax(ill))/(1.*nzt)     
          
	  lxtemp = 0.0
          gx(1)=0.0
	  do 10 i=2,nxt+1
            lxtemp=lxtemp+(ax+2.*(i-2)*(1-ax)/(nxt-1))*Dx   
            gx(i)=lxtemp
 10     continue

	  lytemp = 0.0
	  gy(1)=0.0
	  do 11 j=2,nyt+1
            lytemp=lytemp+(ay+2.*(j-2)*(1-ay)/(nyt-1))*Dy   
            gy(j)=lytemp
 11       continue
	  

	  lztemp = 0.0
	  gz(1)=0.0
	  do 12 k=2,nzt+1
	    lztemp=lztemp+(az+2.*(k-2)*(1-az)/(nzt-1))*Dz   
	    gz(k)=lztemp
 12       continue
      
        do 13 i=1,nxt+1
        do 13 j=1,nyt+1
	do 13 k=1,nzt+1
	  coor(1,nxs+i,nys+j,nzs+k)=gx(i)
 	  coor(2,nxs+i,nys+j,nzs+k)=gy(j)     
 	  coor(3,nxs+i,nys+j,nzs+k)=gz(k)     
 13     continue	


c------------------
c  estrapolo i gx, gy, gz e le dimensioni max e min delle celle 
c    per ul file usr segm
c --------------------      

      write (41,*) 
     &'Dimensione della più piccola e più grande cella su x > 0'
      write (41,*)  gx(2)-gx(1),gx(nxt+1)-gx(nxt)    
      write (41,*) 
     &'Dimensione della più piccola e più grande cella su y > 0'
      write (41,*)  gy(2)-gy(1),gy(nyt+1)-gy(nyt) 
      write (41,*) 
     &'Dimensione della più piccola e più grande cella su z > 0'
      write (41,*)  gz(2)-gz(1),gz(nzt+1)-gz(nzt)
      
      endif
            


      if (sg.eq.2) then
	  
	  nxt=nx(2,ill)
	  nyt=ny(1,ill)
	  nzt=nz(1,ill)
	  
	  nxs=0
	  nys=ny(2,ill)
	  nzs=nz(2,ill)
	  
	  Dx = (xmin(ill))/(1.*nxt)
	  Dy = (ymax(ill))/(1.*nyt)
	  Dz = (zmax(ill))/(1.*nzt)

      	  lxtemp = 0.0
          gx(1)=0.0
	  do 20 i=2,nxt+1
            lxtemp=lxtemp+(ax+2.*(i-2)*(1-ax)/(nxt-1))*Dx   
            gx(i)=lxtemp
 20     continue

	  lytemp = 0.0
	  gy(1)=0.0
	  do 21 j=2,nyt+1
            lytemp=lytemp+(ay+2.*(j-2)*(1-ay)/(nyt-1))*Dy   
            gy(j)=lytemp
 21       continue
	  

	  lztemp = 0.0
	  gz(1)=0.0
	  do 22 k=2,nzt+1
	    lztemp=lztemp+(az+2.*(k-2)*(1-az)/(nzt-1))*Dz   
	    gz(k)=lztemp
 22       continue

	do 23 i=1,nxt+1
        do 23 j=1,nyt+1
	do 23 k=1,nzt+1
	  coor(1,nxt-i+2,nys+j,nzs+k)=gx(i)
 	  coor(2,nxt-i+2,nys+j,nzs+k)=gy(j)     
 	  coor(3,nxt-i+2,nys+j,nzs+k)=gz(k)     
 23     continue	
      endif

      
      if (sg.eq.3) then
	  
	  nxt=nx(2,ill)
	  nyt=ny(2,ill)
	  nzt=nz(1,ill)
	  
	  nxs=0
	  nys=0
	  nzs=nz(2,ill)
	    
	  
	  Dx = (xmin(ill))/(1.*nx(2,ill))
	  Dy = (ymin(ill))/(1.*ny(2,ill))
	  Dz = (zmax(ill))/(1.*nz(1,ill))
      
          lxtemp = 0.0
          gx(1)=0.0
	  do 30 i=2,nxt+1
            lxtemp=lxtemp+(ax+2.*(i-2)*(1-ax)/(nxt-1))*Dx   
            gx(i)=lxtemp
 30     continue

	  lytemp = 0.0
	  gy(1)=0.0
	  do 31 j=2,nyt+1
            lytemp=lytemp+(ay+2.*(j-2)*(1-ay)/(nyt-1))*Dy   
            gy(j)=lytemp
 31       continue
	  

	  lztemp = 0.0
	  gz(1)=0.0
	  do 32 k=2,nzt+1
	    lztemp=lztemp+(az+2.*(k-2)*(1-az)/(nzt-1))*Dz   
	    gz(k)=lztemp
 32       continue
	
	
	do 33 i=1,nxt+1
        do 33 j=1,nyt+1
	do 33 k=1,nzt+1
	  coor(1,nxt-i+2,nyt-j+2,nzs+k)=gx(i)
 	  coor(2,nxt-i+2,nyt-j+2,nzs+k)=gy(j)     
 	  coor(3,nxt-i+2,nyt-j+2,nzs+k)=gz(k)     
 33     continue
      endif
      
      
      if (sg.eq.4) then
	  
	  Dx = (xmax(ill))/(1.*nx(1,ill))
	  Dy = (ymin(ill))/(1.*ny(2,ill))
	  Dz = (zmax(ill))/(1.*nz(1,ill))
      	  
	  nxt=nx(1,ill)
	  nyt=ny(2,ill)
	  nzt=nz(1,ill)
	  
	  nxs=nx(2,ill)
	  nys=0
	  nzs=nz(2,ill)

      
          lxtemp = 0.0
          gx(1)=0.0
	  do 40 i=2,nxt+1
            lxtemp=lxtemp+(ax+2.*(i-2)*(1-ax)/(nxt-1))*Dx   
            gx(i)=lxtemp
 40     continue

	  lytemp = 0.0
	  gy(1)=0.0
	  do 41 j=2,nyt+1
            lytemp=lytemp+(ay+2.*(j-2)*(1-ay)/(nyt-1))*Dy   
            gy(j)=lytemp
 41       continue
	  

	  lztemp = 0.0
	  gz(1)=0.0
	  do 42 k=2,nzt+1
	    lztemp=lztemp+(az+2.*(k-2)*(1-az)/(nzt-1))*Dz   
	    gz(k)=lztemp
 42       continue
	
	
	do 43 i=1,nxt+1
        do 43 j=1,nyt+1
	do 43 k=1,nzt+1
	  coor(1,nxs+i,nyt-j+2,nzs+k)=gx(i)
 	  coor(2,nxs+i,nyt-j+2,nzs+k)=gy(j)     
 	  coor(3,nxs+i,nyt-j+2,nzs+k)=gz(k)     
 43     continue
      endif

      
      if (sg.eq.5) then	  

	  nxt=nx(1,ill)
	  nyt=ny(1,ill)
	  nzt=nz(2,ill)
	  
	  nxs=nx(2,ill)
	  nys=ny(2,ill)
	  nzs=0
	  
	  Dx = (xmax(ill))/(1.*nxt)
	  Dy = (ymax(ill))/(1.*nyt)
	  Dz = (zmin(ill))/(1.*nzt)     
          
	  lxtemp = 0.0
          gx(1)=0.0
	  do 50 i=2,nxt+1
            lxtemp=lxtemp+(ax+2.*(i-2)*(1-ax)/(nxt-1))*Dx   
            gx(i)=lxtemp
 50     continue

	  lytemp = 0.0
	  gy(1)=0.0
	  do 51 j=2,nyt+1
            lytemp=lytemp+(ay+2.*(j-2)*(1-ay)/(nyt-1))*Dy   
            gy(j)=lytemp
 51       continue
	  

	  lztemp = 0.0
	  gz(1)=0.0
	  do 52 k=2,nzt+1
	    lztemp=lztemp+(az+2.*(k-2)*(1-az)/(nzt-1))*Dz   
	    gz(k)=lztemp
 52       continue
	

        do 53 i=1,nxt+1
        do 53 j=1,nyt+1
	do 53 k=1,nzt+1
	  coor(1,nxs+i,nys+j,nzt-k+2)=gx(i)
 	  coor(2,nxs+i,nys+j,nzt-k+2)=gy(j)     
 	  coor(3,nxs+i,nys+j,nzt-k+2)=gz(k)     
 53     continue	
      endif
      
      
      if (sg.eq.6) then
	  
	  nxt=nx(2,ill)
	  nyt=ny(1,ill)
	  nzt=nz(2,ill)
	  
	  nxs=0
	  nys=ny(2,ill)
	  nzs=0
	  
	  Dx = (xmin(ill))/(1.*nxt)
	  Dy = (ymax(ill))/(1.*nyt)
	  Dz = (zmin(ill))/(1.*nzt)

      	  lxtemp = 0.0
          gx(1)=0.0
	  do 60 i=2,nxt+1
            lxtemp=lxtemp+(ax+2.*(i-2)*(1-ax)/(nxt-1))*Dx   
            gx(i)=lxtemp
 60     continue

	  lytemp = 0.0
	  gy(1)=0.0
	  do 61 j=2,nyt+1
            lytemp=lytemp+(ay+2.*(j-2)*(1-ay)/(nyt-1))*Dy   
            gy(j)=lytemp
 61       continue
	  

	  lztemp = 0.0
	  gz(1)=0.0
	  do 62 k=2,nzt+1
	    lztemp=lztemp+(az+2.*(k-2)*(1-az)/(nzt-1))*Dz   
	    gz(k)=lztemp
 62       continue

	do 63 i=1,nxt+1
        do 63 j=1,nyt+1
	do 63 k=1,nzt+1
	  coor(1,nxt-i+2,nys+j,nzt-k+2)=gx(i)
 	  coor(2,nxt-i+2,nys+j,nzt-k+2)=gy(j)     
 	  coor(3,nxt-i+2,nys+j,nzt-k+2)=gz(k)     
 63     continue
      endif


      if (sg.eq.7) then
	  
	  nxt=nx(2,ill)
	  nyt=ny(2,ill)
	  nzt=nz(2,ill)
	  
	  nxs=0
	  nys=0
	  nzs=0
	    
	  Dx = (xmin(ill))/(1.*nx(2,ill))
	  Dy = (ymin(ill))/(1.*ny(2,ill))
	  Dz = (zmin(ill))/(1.*nz(1,ill))
      
          lxtemp = 0.0
          gx(1)=0.0
	  do 70 i=2,nxt+1
            lxtemp=lxtemp+(ax+2.*(i-2)*(1-ax)/(nxt-1))*Dx   
            gx(i)=lxtemp
 70     continue

	  lytemp = 0.0
	  gy(1)=0.0
	  do 71 j=2,nyt+1
            lytemp=lytemp+(ay+2.*(j-2)*(1-ay)/(nyt-1))*Dy   
            gy(j)=lytemp
 71       continue
	  

	  lztemp = 0.0
	  gz(1)=0.0
	  do 72 k=2,nzt+1
	    lztemp=lztemp+(az+2.*(k-2)*(1-az)/(nzt-1))*Dz   
	    gz(k)=lztemp
 72       continue
	
	
	do 73 i=1,nxt+1
        do 73 j=1,nyt+1
	do 73 k=1,nzt+1
	  coor(1,nxt-i+2,nyt-j+2,nzt-k+2)=gx(i)
 	  coor(2,nxt-i+2,nyt-j+2,nzt-k+2)=gy(j)     
 	  coor(3,nxt-i+2,nyt-j+2,nzt-k+2)=gz(k)     
 73     continue

      if (xmin(ill).ne.0.0) then
         write (41,*) 
     &'Dimensione della più piccola e più grande cella su x < 0'
         write (41,*)  gx(1)-gx(2),gx(nxt)-gx(nxt+1)
       else 
	  write(41,*)'No cells on x < 0'
	endif
       
       if (ymin(ill).ne.0.0) then
         write (41,*) 
     &'Dimensione della più piccola e più grande cella su y < 0'
         write (41,*)  gy(1)-gy(2),gy(nyt)-gy(nyt+1)
       else
         write (41,*) 'No cells on y < 0'
	endif
	
	if (zmin(ill).ne.0.0) then
         write (41,*) 
     &'Dimensione della più piccola e più grande cella su z < 0'
         write (41,*)  gz(1)-gz(2),gz(nzt)-gz(nzt+1)
	else 
	  write(41,*)'No cells on z < 0'
	  endif
      
      
      
      
      
      
      
      endif
      
      
      if (sg.eq.8) then
	  
	  Dx = (xmax(ill))/(1.*nx(1,ill))
	  Dy = (ymin(ill))/(1.*ny(2,ill))
	  Dz = (zmin(ill))/(1.*nz(1,ill))
      	  
	  nxt=nx(1,ill)
	  nyt=ny(2,ill)
	  nzt=nz(2,ill)
	  
	  nxs=nx(2,ill)
	  nys=0
	  nzs=0

      
          lxtemp = 0.0
          gx(1)=0.0
	  do 80 i=2,nxt+1
            lxtemp=lxtemp+(ax+2.*(i-2)*(1-ax)/(nxt-1))*Dx   
            gx(i)=lxtemp
 80     continue

	  lytemp = 0.0
	  gy(1)=0.0
	  do 81 j=2,nyt+1
            lytemp=lytemp+(ay+2.*(j-2)*(1-ay)/(nyt-1))*Dy   
            gy(j)=lytemp
 81       continue
	  

	  lztemp = 0.0
	  gz(1)=0.0
	  do 82 k=2,nzt+1
	    lztemp=lztemp+(az+2.*(k-2)*(1-az)/(nzt-1))*Dz   
	    gz(k)=lztemp
 82       continue
	
	
	do 83 i=1,nxt+1
        do 83 j=1,nyt+1
	do 83 k=1,nzt+1
	  coor(1,nxs+i,nyt-j+2,nzt-k+2)=gx(i)
 	  coor(2,nxs+i,nyt-j+2,nzt-k+2)=gy(j)     
 	  coor(3,nxs+i,nyt-j+2,nzt-k+2)=gz(k)     
 83     continue
       endif  

c ---------------------
c controlli nel caso uno dei minimi del dominio venga posto a zero
c sovrascrittura eventuale errore nella definizione del numero di celle
c
c ---------------------       


 3      continue
        
     
      
      ilg(ill)=nx(2,ill)+nx(1,ill)+1
      jlg(ill)=ny(2,ill)+ny(1,ill)+1
      klg(ill)=nz(2,ill)+nz(1,ill)+1


c --------------         
c  copio coor in x2
c  ma prima trovo le dimensioni di x2 per il caso proposto
c ---------------

      iv(ill)=2*(nx(2,ill)+nx(1,ill)+2)+1
      jv(ill)=2*(ny(2,ill)+ny(1,ill)+2)+1
      kv(ill)=2*(nz(2,ill)+nz(1,ill)+2)+1
      print *,'Total vertices in block #',ill,' =',
     & iv(ill),jv(ill),kv(ill)
      print *,''

c ----------------
c copio coor in x2
c ---------------
 
             
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

      do 15 k=2,kv(ill),2
      do 15 j=2,jv(ill),2
      do 15 i=2,iv(ill),2
        x2(ill,i,j,k,1)=0.5*(x2(ill,i,j,k-1,1) +
     &                                  x2(ill,i,j,k+1,1))
        x2(ill,i,j,k,2)=0.5*(x2(ill,i,j,k-1,2) +
     &                                  x2(ill,i,j,k+1,2))
        x2(ill,i,j,k,3)=0.5*(x2(ill,i,j,k-1,3) +
     &                                  x2(ill,i,j,k+1,3))
 15     continue

 2    continue      
      return
      end

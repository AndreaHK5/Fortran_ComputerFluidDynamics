      program tracker
c ------------------------
c   programma per il tracciamento dele particielle
c   (un'altra gino produxction, 'o duri!)
c
c  L'output viene fatto su file tipo pvp_?????, dove pvp è per
c   Position-Velocity-Particles e ????? è il time step relativo
c  Il programma cercherà il file pvp_????? con ????? numero di step
c   di input tr, se INIT è 0. 
c
c -------------------------

      implicit double precision (a-h,o-z)
      include 'main.param'
      double precision ROP,ROF,KVF,DD,GAx,GAy,GAz,EPS,DT,maxR,str
      integer nm,SI,INIT,first,lest,step,th
      integer ilg(ng),jlg(ng),klg(ng)
      double precision G(il,jl,kl,3)
      double precision xmi(ng),xma(ng)
      double precision ymi(ng),yma(ng)
      double precision zmi(ng),zma(ng)
      double precision XMIN,XMAX,YMIN,YMAX,ZMIN,ZMAX
      double precision P0(3,maxpart),P1(3,maxpart)
      double precision U0(3,maxpart),U1(3,maxpart)
      double precision A(il,jl,kl,3)
      double precision Ptemp(3),Pnew(3),Utemp(3),Vtemp(3),Unew(3)
      integer nstep
      double precision K01,K02,K04
      integer R1,R2,R3,R4,R5,R6
c --------------------------------
c  P0 posizione a t-1
c  U0 velocità particielle t-1
c  A è il campo dele velocità letto da rstfile
c  V0 è la velocità del fluido a t-1
c  nstep è il contatore per l'iterazione sui numeri di step
c --------------------------------

      character *9 rstfile
      character *5 stepcha
      character *9 pvpfile

      common /general/G,A,ilg,jlg,klg
      common /param/DD,DT,KVF,K02,K04,GAx,GAy,GAz
      common /rebo/XMAX,XMIN,YMAX,YMIN,ZMAX,ZMIN
     &                                ,R1,R2,R3,R4,R5,R6
     
      
      print *,'   '
      print *,'    --------------------------'
      print *,'    >> a Gino''s production <<'
      print *,'    --------------------------' 
      print *,'     (che credevate, o duri!)'
      print *,''





c -------------------------------
c iniziamo col legegre il file con tuuuuutti i parametri dell'ammmore
c cche si trovano in inputtr
c -------------------------------

       open (unit=2,file='inputtr',form='formatted')
      

       read (2,*) ROP
       read (2,*) ROF
       read (2,*) KVF
       read (2,*) DD
       read (2,*) GAx
       read (2,*) GAy
       read (2,*) GAz
       read (2,*) EPS
       read (2,*) DT
       read (2,*) nm
       read (2,*) maxR
       read (2,*) str
       read (2,*) th
       read (2,*) SI
       read (2,*) INIT
       read (2,*) first
       read (2,*) last
       read (2,*) step
       read (2,*) R1,R2,R3,R4,R5,R6

       close (2)
c --------------------------
c Se si vuole il force tracking ci vogliono altri aprametri
c
c  Lettura inputgr per i valori mi max e min del dominio
c -------------------



      open(unit=10, file='inputgr')
      
      read(10,*)
      read(10,*) lg
      do 2 ig=1,lg
      read(10,*)
      read(10,*) xmi(ig),xma(ig)
      read(10,*) ymi(ig),yma(ig)
      read(10,*) zmi(ig),zma(ig)
      read (10,*)
      read (10,*) 
      read (10,*)
      read (10,*) 
      read (10,*)
      read (10,*) 
      read (10,*)
      read (10,*) 

  2   continue
     
c ---------------------------
c Adesso passiamo alla lettura dei centri griglia 
c (importantissimi x le velocità)
c E delle loro velocità (A)
c ---------------------------      
      
      open(9,file = 'xyu',form = 'unformatted')
      read(9)lg
      
      write (stepcha,999) first
      rstfile='rst_'//stepcha
 999  format (i5.5)

      open (11,file=rstfile, form='unformatted')
      read (11) nt,ntime
     
      do 1 ig=1,lg
       read(9)ilg(ig), jlg(ig), klg(ig)
       read(9)((((G(i,j,k,l),i = 1,ilg(ig)),
     &          j = 1,jlg(ig)), k = 1,klg(ig)), l = 1,3)

       
       read(11) ((((A(i,j,k,l),i=1,ilg(ig)),j=1,jlg(ig)),
     &      k=1,klg(ig)),l=1,3)
       read(11)pp
       read(11)pp
       read(11)pp
       read(11)pp
       read(11)pp
       read(11)pp
       read(11)pp
       read(11)pp
       read(11)pp

       close(11)

c ----------------------------
c inizializzazione del campo di moto e delle posizioni, se 
c  è il primo step
c ----------------------------       

      XMIN=xmi(ig)
      XMAX=xma(ig)
      YMIN=ymi(ig)
      YMAX=yma(ig)
      ZMIN=zmi(ig)
      ZMAX=zma(ig)
	
      if(INIT.EQ.1) then
      
      print *,'Randomly positioning the particles'
      call posin(P0,maxR,str,th,nm,DD)


      if (SI.eq.0) then
       print *,'Velocità iniziale nulla per le particielle'
       print *,''
        do  l=1,3
        do  n=1,nm
          U0(1,n)=0.
        enddo
        enddo
       endif

       if (SI.ne.0) then
        print *,'Particielle in moto con il jet,inizializzazione' 
        do  n=1,nm
          Ptemp(1)=P0(1,n)
	  Ptemp(2)=P0(2,n)
          Ptemp(3)=P0(3,n)
        call velxyz(Ptemp,Utemp)
          U0(1,n)=Utemp(1)
          U0(2,n)=Utemp(2)
          U0(3,n)=Utemp(3)

	
	enddo
       endif

      endif


      if (INIT.eq.0)then
        pvpfile='pvp_'//stepcha     
        open (12,file='pvpstart',form='unformatted')
        print *,
     &'Reading initial position and velocity from file "pvpstart"'
        print *,''
        do n=1,nm
         read (12) (P0(l,n),l=1,3),(U0(h,n),h=1,3)
        enddo

        close(12)
      endif


c ----------------------------
c calcolo di qualche paramentro interessante che poi andarà in
c  qualche common
c ----------------------------	

     	K04=ROF/ROP
        K01=KVF*18/(DD*DD)
        K02=K01*K04

c ---------------------------
c Inizio del ciclo sui time steps
c ---------------------------


      do 3 nstep=first,last,step
        
	print *,''
	print *,''      
	print *,'>>>>>>>>>step #',nstep,'<<<<<<<<<<<<'
  
c --------------------------
c Se è un entrata successiva legge P0 e V0 e le velocità del fluido
c  da file precedente
c ---------------------------


       if (nstep.gt.first) then
         
         write (stepcha,999) nstep-step 
         pvpfile='pvp_'//stepcha     
         open (12,file=pvpfile,form='unformatted')
         print *,'Reading position and velocity from ',pvpfile
         do n=1,nm
           read (12) (P0(l,n),l=1,3),(U0(h,n),h=1,3)
         enddo
         close(12)

	 write (stepcha,999) nstep
	 rstfile='rst_'//stepcha
         print *,'Reading fluid velocity from ',rstfile
        

         open (11,file=rstfile, form='unformatted')
         read (11) nt,ntime
       
         read(11) ((((A(i,j,k,l),i=1,ilg(ig)),j=1,jlg(ig)),
     &      k=1,klg(ig)),l=1,3)
         read(11)pp
         read(11)pp
         read(11)pp
         read(11)pp
         read(11)pp
         read(11)pp
         read(11)pp
         read(11)pp
         read(11)pp

        close(11)

      endif

c -------------------------
c Calcolo velocità fluido nelle posizioni delle particielle,
c  delle nuove velocità
c  e successivamente calcolo della nuova posizione della particiella
c -------------------------

       do  n=1,nm
          Ptemp(1)=P0(1,n)
	  Ptemp(2)=P0(2,n)
          Ptemp(3)=P0(3,n)
          Utemp(1)=U0(1,n)
          Utemp(2)=U0(2,n)
          Utemp(3)=U0(3,n)
	  
c ----------------------------
c  Se la particiella è sul bordo vuol dire che è uscirta, skippo
c   le operazioni di calcolo...
c  L'unico dato trattato diversamente è Zmin, xche se nò nn mi prende le cond a
c    contorno e me le elimina...
c ----------------------------	  
	  
	  if (Ptemp(1).eq.XMIN.or.Ptemp(1).eq.XMAX.or.
     &               Ptemp(2).eq.YMIN.or.Ptemp(2).eq.YMAX.or.
     &               Ptemp(3).eq.ZMAX) then  
           P1(1,n)=Ptemp(1)
	   P1(2,n)=Ptemp(2)
	   P1(3,n)=Ptemp(3)
	   U1(1,n)=Utemp(1)
	   U1(2,n)=Utemp(2)
	   U1(3,n)=Utemp(3)
c	   print *,'Particella uscita dal dominio'
	   goto 5
	  endif

c -----------------------
c In quiesta maniera puàò essere fissato come uscita anche 
c il bordo inferiore
c ----------------------
         
	 if (Ptemp(3).eq.ZMIN.and.
     &           sqrt(Ptemp(1)**2+Ptemp(2)**2).gt.maxR)then
	  P1(1,n)=Ptemp(1)
	   P1(2,n)=Ptemp(2)
	   P1(3,n)=Ptemp(3)
	   U1(1,n)=Utemp(1)
	   U1(2,n)=Utemp(2)
	   U1(3,n)=Utemp(3)
c	   print *,'Particella uscita dal dominio'
	   goto 5
	  endif

c ----------------	  
c Per le particielle che hanno punto di partenza interno (strettamente) al
c  dominio inizia il calcolo	  
c ---------------------	  
	  call velxyz(Ptemp,Vtemp)
	  call velpar(Utemp,Vtemp,Unew)
	  call pospar(Ptemp,Pnew,Utemp,Unew)
      	  U1(1,n)=Unew(1)
	  U1(2,n)=Unew(2)
	  U1(3,n)=Unew(3)
      	  P1(1,n)=Pnew(1)
	  P1(2,n)=Pnew(2)
	  P1(3,n)=Pnew(3)	  
 5      continue
	enddo     
        
	
	write (stepcha,999) nstep 
	pvpfile='pvp_'//stepcha     
        open (13,file=pvpfile,form='unformatted')
	print *,'Writing velocity on ',pvpfile
	do 4 n=1,nm
         write (13) (P1(l,n),l=1,3),(U1(h,n),h=1,3)
   4   continue
        close(13)
   3    continue
   
   
   
   1    continue
      
      print *,''
      print *,'End of records,'
      print *,'Have a nice day'
      print *,''
      stop
 
      end

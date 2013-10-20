c---------------------------------------------------------------
c   program to compute the initial and boundary velocity conditions
c   for the input of 3-D Composite Navier-Stokes code
c
c  un'altra Gino's production
c----------------------------------------------------------------

      
      program boundary
      implicit double precision (a-h,o-z)
      include 'main.param'
      parameter(ijm = 160)
      double precision Ucf,Ujet,Rjet,Dcf,Djet,ux
      double precision Rlim
      integer ncx(ng),ncy(ng),ncz(ng),idc,pos
      double precision xp(il,jl,kl,3,ng)
      double precision sini(il,jl,kl,ng)
      double precision uby(ijm,ijm,3,6,ng),
     & ubyconvold(ijm,ijm,3,6,ng,3),
     & ubyconvold2(ijm,ijm,3,6,ng,3),Deta(ijm,ijm,3,ng)
      logical bsoli(il,ng),bfree(il,ng)

c ------------------------
c	real Uinflo,Uinflo1,ux
c rivedere xche in caso sono da ri dichiarare Ujet e amici vari
c ------------------------


c -----------------------------
c ijm dimensione max del blocco (usata per assegnare boundary su faccia)
c ncx,ncy,ncz numero di celle
c xp è la matrice con dentro i centri volume
c Rlim è il raggio a cui parte lo strato limite
c idc= indice per la ricorsione di direzione in calcolodi DETA
c sini è un dummy, saranno dele condizioni di turbolenza (inizializzato a zero e
c         scrito)
c  uby velocita' sulla frontiera: 2 indici di faccia, componente,
c   faccia del blocco, blocco di appartenenza
c -----------------------------




c -----------
c apertura file's
c
c
c ------------------


      open(9, file = 'inputbo')
      open(11,file = 'xyu',form = 'unformatted')
      open(22,file = 'uby.in',form = 'unformatted')
      open(23,file = 'ublog',form = 'formatted')
      open(24,file = 'ubyconv',form = 'unformatted')
      open(25,file = 'Detaconv',form = 'unformatted')
c      open(20,file = 'usrvel')

      open(55,file = 'usrjetcheck')      
      open(56,file = 'usrincheck')
c ----------
c ubin ha i valori iniziali di vel cf e varie
c uby.in ha dentro i valori delle vel sulle facce
c il file 23 aveva delle scritture (commentate) di bsoli
c   e bfree, che sn state saldobrasate
c ubyconv contiene le condizioni di partenza
c Detaconv le info per porre le condizioni di differenziale 
c    zero in uscita
c  pippo è in fase di definizione
c-----------------
     
      print *,'   '
      print *,'    --------------------------'
      print *,'    >> a Gino''s production <<'
      print *,'    --------------------------' 
      print *,''
      
      
c ---------------------
c  Lettura di ubin
c ---------------------
      
      read (9,*)
      read (9,*) Ucf
      read (9,*)
      read (9,*) Ujet
      read (9,*)
      read (9,*) Rjet
      read (9,*)
      read (9,*) Dcf
      read (9,*)
      read (9,*) Djet


      print *,'Riepilogo dei Goal dall"Olimpico in Roma'
      print *,''
      print *,'         Velocità CF',Ucf
      print *,'        Velocità jet',Ujet
      print *,'          Raggio jet',Rjet
      print *,'   Spess str.lim. CF',Dcf
      print *,'  Spess str.lim. jet',Djet
      print *,''
      print *,'Linea allo studio per i secondi tempi' 
      print *,''

      print *,''
      
      
      Rlim=Rjet-Djet


      read(11)lg
      do 16 ig = 1,lg
       read(11)ncx(ig),ncy(ig),ncz(ig)
       read(11)((((xp(i,j,k,l,ig), i = 1,ncx(ig)), j = 1,ncy(ig)),
     &           k = 1,ncz(ig)), l = 1,3)
 16   continue

c ----------------------------------------
c delimitazione delle zone con condizioni al contorno miste
c   credo bsoli e b free siano duna sorta di "flag" x il solutore
c -------------------------------------

       do 100 ig=1,lg
       do 100 i=1,ncx(ig)
           bsoli(i,ig)=.true.
           bfree(i,ig)=.false.
  100    continue




c-------------------------------------
c condizioni convettive su sezione normale z
c inizializzazione a zero della matrice uby
c-------------------------------------
      
      do 10 ig = 1,lg
      do 15 mb = 1,6
      do 21 ni = 1,3
      do 25 k  = 1,ijm
      do 30 j  = 1,ijm
        uby(j,k,ni,mb,ig) = 0.d0
        do pos=1,3
          ubyconvold(j,k,ni,mb,ig,pos) = 0.d0
          ubyconvold2(j,k,ni,mb,ig,pos) = 0.d0
        enddo 
  30  continue
  25  continue
  21  continue
  15  continue
  10  continue

c --------------------------
c assegnazione velocita' a seconda della faccia di appartenenza
c --------------------------
       
       do 1 ig=1,lg
        do mb=1,6

c --------------
c facce 1 e 3
c --------------
         
	 if(mb.eq.1.or.mb.eq.3)then
          if (mb.eq.1) then 
           ux= Ucf
          else
           ux= 0.d0
          endif
          nparab=0
          do 2,k=1,ncz(ig)
           z=xp(1,1,k,3,ig)
           if(z.le.Dcf)then
            parabz=(-z**2.+2*Dcf*z)/Dcf**2.
            nparab=nparab+1
           else
            parabz=1.
           endif
           do 2,j=1,ncy(ig)
             uby(j,k,1,mb,ig)=ux*parabz
             do 2,pos=1,3
               ubyconvold(j,k,1,mb,ig,pos) = ux*parabz
               ubyconvold2(j,k,1,mb,ig,pos) = ux*parabz
 2           continue 
          endif
          
	  if (mb.eq.1) then   
	     print *,'Centri faccia in str.lim.C.F.=',nparab
             print *,''
	  endif

c ---------------
c facce 2 e 4
c ---------------

         if(mb.eq.2.or.mb.eq.4)then
           ux= 0.0
             do 3,k=1,ncz(ig)
             do 3,i=1,ncx(ig)
              uby(j,k,2,mb,ig)=ux
              do 3,pos=1,3
               ubyconvold(i,k,2,mb,ig,pos) = ux
               ubyconvold2(i,k,2,mb,ig,pos) = ux
 3           continue
          endif

c --------------
c facce 5 e 6
c --------------

         if(mb.eq.5.or.mb.eq.6)then
           if(mb.eq.6)then
           ux= 0.d0
             do 4,j=1,ncy(ig)
             do 4,i=1,ncx(ig)
              uby(i,j,1,mb,ig)=ux
             do 4,pos=1,3
               ubyconvold(i,j,3,mb,ig,pos) = ux
               ubyconvold2(i,j,3,mb,ig,pos) = ux
 4           continue
            else
	     nmixte=0
             do 5,j=1,ncy(ig)
             do 5,i=1,ncx(ig)

c ------------------
c delimitazione celle in condizione di inlet
c ------------------

              erre=(xp(i,j,1,1,ig)**2.+xp(i,j,1,2,ig)**2.)**0.5
              if(erre.lt.Rlim)then
                nmixte=nmixte+1
                ux=Ujet*(1.-tanh((Rlim/Djet)/4.*(2*erre-1./(2*erre))))/2.
              elseif(Rlim.le.erre.and.erre.lt.Rjet)then
                nmixte=nmixte+1
	        ux=Ujet*(1.-(erre-Rlim)/(Rjet-Rlim))
              else
                ux=0. 
              endif
              uby(i,j,3,mb,ig)=ux
              if(erre.lt.Rjet)then
c                print *,'Point',2*i,2*j,' (i,j) defined as inlet with
c     & uby= ',uby(i,j,3,mb,ig)
              endif
              do pos=1,3
               ubyconvold(i,j,3,mb,ig,pos) = ux
               ubyconvold2(i,j,3,mb,ig,pos) = ux
              enddo
 5           continue
            endif
           endif
          enddo
 1     continue

       
c ---------------------       
c Rapido recap su tutti sti endif ed enddo: 
c  il primo enddo chiude il do su ubconvold     
c  il continue 5 chiude la ricorsione sul "fondo per i e j
c  endif che chiude l'if su mb= 5 o 6
c  endif che chiude sul controllo di tutti gli mb
c  enddo sulla ricorsione di mb
c  1 che chiude su i per il multiblocco
c ------------------------       
       
       if(nmixte.eq.0)then
         print *,'no cell face defined as mixed boundary'
	 print *,'hence this program won"t work, no inlet...'
         print *,'check cmin and ubin!!'
         print *,''
       else
         print *,''
	 print *,nmixte,' cell faces defined as jet inlet'
         print *,''
       endif

c ---------------------------------
c scrittura valori di output
c 1. uby sulle faccie del dominio (file uby.in)
c 2. sini nei centri cella        (file uby.in)
c 3. bsoli, bfree                 (file ubylog)
c 
c e in più c'è anche l'inizializzazione di sini a zero
c -----------------------------------

      do 3000 ig = 1,lg


       do 425 i=1,ncx(ig)
       do 425 j=1,ncy(ig)
       do 425 k=1,ncz(ig)
           
	 sini(i,j,k,ig)=0.0
       
 425    continue
	
	
	nim = max0(ncx(ig),ncy(ig),ncz(ig))

        write(22)((((uby(i,j,mv,mb,ig), i = 1,nim), j = 1,nim),
     &            mv = 1,3), mb = 1,6)
        write(22)(((sini(i,j,k,ig),i = 1,ncx(ig)), j = 1,ncy(ig)),
     &                              k = 1,ncz(ig))

        write(23,*) (bsoli(i,ig),i= 1,ncx(ig))
        write(23,*) (bfree(i,ig),i= 1,ncx(ig))


c --------------------	
c scrittura su fantomatico file 20 (formattato), che cos'è??
c ----------------------	
	
c	do mb=1,6
c         do i=1,nim
c          do j=1,nim
c           write(20,202)mb,i,j,(uby(i,j,mv,mb,ig),mv=1,3)
c 202         format(3(1x,i4),3(1x,f12.8))
c          enddo
c          write(20,*)
c         enddo
c        write(20,*)
c        enddo
3000   continue


c ----------------------------------
c Calcolo della distanza tra punti griglia in direzione normale
c all'outflow convettivo: la distanza e' pari a 2 passi griglia
c (vedi espressione per dv/deta su Ferziger (3.33 p 51)
c La distanza viene ricostriuta dalla posizione dei centrocella
c prima calcolo il punto medio L tra nj (ultimo) e nj-1 (penultimo)
c ----------------------------------

      do 4000 ig=1,lg 
        nig = ncx(ig)
        njg = ncy(ig)
	nkg = ncz(ig)
        
	do idc=1,3
 
c --------------------
c direzione x normal (faccia 3)
c --------------------         
	  
	if(idc.eq.1)then
          do 4002 j = 1,njg
          do 4002 k = 1,nkg

            xL=0.5*(xp(nig-1,j,k,1,ig)+xp(nig,j,k,1,ig))
            yL=0.5*(xp(nig-1,j,k,2,ig)+xp(nig,j,k,2,ig))
            zL=0.5*(xp(nig-1,j,k,3,ig)+xp(nig,j,k,3,ig))

c -----------------------
c poi calcolo il punto medio H tra nj-2 e nj-3
c ----------------------

            xH=0.5*(xp(nig-3,j,k,1,ig)+xp(nig-2,j,k,1,ig))
            yH=0.5*(xp(nig-3,j,k,2,ig)+xp(nig-2,j,k,2,ig))
            zH=0.5*(xp(nig-3,j,k,3,ig)+xp(nig-2,j,k,3,ig))

c -----------------
c quindi calcolo Deta= 2 delta z
c -----------------

	    Deta(j,k,idc,ig)=sqrt((xL-xH)**2+(yL-yH)**2+(zL-zH)**2)
        if(Deta(j,k,idc,ig).EQ.0) print*,'DETA =0 
     &                       for face x direction ',j,k,'that"s bad!'
 4002   continue
       endif

c ------------------
c direzione 2 (profondità)
c ------------------

       if(idc.eq.2)then
        do 4003 i = 1,nig
        do 4003 k = 1,nkg
          xL=0.5*(xp(i,njg-1,k,1,ig)+xp(i,njg,k,1,ig))
          yL=0.5*(xp(i,njg-1,k,2,ig)+xp(i,njg,k,2,ig))
          zL=0.5*(xp(i,njg-1,k,3,ig)+xp(i,njg,k,3,ig))
          xH=0.5*(xp(i,njg-3,k,1,ig)+xp(i,njg-2,k,1,ig))
          yH=0.5*(xp(i,njg-3,k,2,ig)+xp(i,njg-2,k,2,ig))
          zH=0.5*(xp(i,njg-3,k,3,ig)+xp(i,njg-2,k,3,ig))
          Deta(i,k,idc,ig)=sqrt((xL-xH)**2+(yL-yH)**2+(zL-zH)**2)
          if(Deta(i,k,idc,ig).EQ.0) print*,'DETA =0
     &                     for face y direction ',i,k,'that"s bad!'
 4003  continue
       endif

c --------------------
c direzione 3 (verticale)
c ---------------------
        
	if(idc.eq.3)then
	  do 4001 i = 1,nig
          do 4001 j = 1,njg      
           xL=0.5*(xp(i,j,nkg-1,1,ig)+xp(i,j,nkg,1,ig))
	   yL=0.5*(xp(i,j,nkg-1,2,ig)+xp(i,j,nkg,2,ig))
	   zL=0.5*(xp(i,j,nkg-1,3,ig)+xp(i,j,nkg,3,ig))
           xH=0.5*(xp(i,j,nkg-3,1,ig)+xp(i,j,nkg-2,1,ig))
	   yH=0.5*(xp(i,j,nkg-3,2,ig)+xp(i,j,nkg-2,2,ig))
	   zH=0.5*(xp(i,j,nkg-3,3,ig)+xp(i,j,nkg-2,3,ig))
           Deta(i,j,idc,ig)=sqrt((xL-xH)**2+(yL-yH)**2+(zL-zH)**2)
           if(Deta(i,j,idc,ig).EQ.0) print*,'DETA =0
     &                      for face z direction',i,j,'that"s bad!'
 4001   continue
       endif

       enddo
4000  continue

c ----------------
c enddo su ricorsione fino a 3 per le direzioni
c continue 4000 sul numero di blocchi
c -----------------------


c------------------------------------------
c now write to files 
c ubyconvold, ubyconvold2 (file ubyconv)
c deta                    (file Detaconv)
c------------------------------------------
      do 4100 ig = 1,lg
         nig = ncx(ig)
         njg = ncy(ig)
         nkg = ncz(ig)
	 nim = max0(ncx(ig),ncy(ig),ncz(ig))

         write(24)(((((ubyconvold(i,j,l,mb,ig,pos), i=1,nim), j=1,nim),
     &          l=1,3), mb=1,6), pos=1,3)

         write(24)(((((ubyconvold2(i,j,l,mb,ig,pos),i=1,nim), j=1,nim),
     &          l=1,3), mb=1,6), pos=1,3) 

         write(25)(((Deta(i,j,idc,ig),i=1,nim),j=1,nim),idc=1,3)


c ----------------------
c  scrittura file usrjetcheck
c  contiene le velocità sulla faccia inferiore
c -------------


      do 2710 i=1,nig
      do 2710 j=1,njg
        if (uby(i,j,3,5,ig).ne.0.0)then
         write(55,1001) xp(i,j,1,1,ig),xp(i,j,1,2,ig),uby(i,j,3,5,ig)
        endif
2710    continue

      do 3011 j=1,njg
      do 3011 k=1,nkg
         
	 write(56,1001) xp(1,j,k,2,ig),xp(1,j,k,3,ig),uby(j,k,1,1,ig)

 3011   continue

4100	continue

 1001 format (3(1x,f12.8))

      end







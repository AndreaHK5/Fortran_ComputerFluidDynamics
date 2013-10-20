             program initial
c--------------------------------------------------
c
c   Initialization to an homogeneous velocity field at
c      ambient pressure for multi domain
c
c
c --------------------------------------
      
      implicit double precision (a-h,o-z)
      include 'main.param'
      double precision uinf(ng),vinf(ng),winf(ng),pinf(ng)
      double precision u(il,jl,kl,3),p(il,jl,kl)
      double precision xp(il,jl,kl,3)
      double precision uxi(il,jl,kl),uej(il,jl,kl),uzk(il,jl,kl)
      double precision xix(il,jl,kl), etx(il,jl,kl),
     &                 xiy(il,jl,kl), ety(il,jl,kl),
     &                 xiz(il,jl,kl), etz(il,jl,kl),
     &                 ztx(il,jl,kl), zty(il,jl,kl),
     &                 ztz(il,jl,kl), jac(il,jl,kl)
      double precision dummy(il,jl,kl,3)
      double precision dummy2(il,jl,kl)
      double precision dummy3(nmax,nmax,4,6)

      integer ilg(ng),jlg(ng),klg(ng),icord(ng)
      integer nilg(ng),njlg(ng),nklg(ng)
      integer ncx,ncy,ncz
      integer nx(2,ng),ny(2,ng),nz(2,ng)


c --------------------------
c  xp è dove immagazzino tutti i centri cella
c  dummy è un pò un brodo...
c  xix e compagnia varia nn è il diciannovesimo secolo ma i
c    ns amici vattori controvarianti
c  dummy 2 serve a leggere (e scartare) i g11,g22,g33 che sono scritti in
c         metric
c  In generale le variabili dummy servono a "scrivere" nei file degli zeri
c    dove servono per gli spazi voluto (se vai sotto capisci che cosa significa)  
c --------------------------
      
      open(unit=10, file='inputgr')
      open(unit=1,file='inputin',
     &     status='old',form='formatted')
      open(12,file = 'metric',form = 'unformatted')
      open(22,file = 'rstore',  form = 'unformatted')
      open(55,file = 'dummy')
      
c --------------------
c  
c  cmin is the same input file for griglia (hence data MUST be correct)
c  inprini contains the input of velocity and speed for the domain
c  metric ha i calcoli di griglia.f
c  rstore ancora nn è chiaro sosa faccia, ci inizio scrivendoci uno zero e poi
c    le velocità di centro griglia (uhmmm, griglia!) 
c
c
c ----------------------

    
      print *,'   '
      print *,'    --------------------------'
      print *,'    >> a Gino''s production <<'
      print *,'    --------------------------' 
      print *,''



c-------------------------------------------
c Read the number of grid points
c
c it was done by reading datain, a file written by the gril mkr
c  instead, to get a rid of some files ad to obtain a lighter application,
c  it is done now via cmin readin (as in coordmkr) but wout the check on max
c  available cells (already done) 
c-------------------------------------------
      
      read(10,*)
      read(10,*) lg
      do 1,ill=1,lg
       read(10,*)
       read(10,*) 
       read(10,*) 
       read(10,*)
       read(10,*) 
       read(10,*) nx(2,ill),nx(1,ill)
       read(10,*)
       read(10,*) ny(2,ill),ny(1,ill)
       read(10,*)
       read(10,*) nz(2,ill),nz(1,ill)
       read(10,*)
       read(10,*)

 1    continue


c-------------------------------------------------------
c Read velocity components at infinity and ambient pressure
c da mettere per ogni blocco
c-------------------------------------------------------

       do ig=1,lg 
        
	nilg(ig)=nx(1,ig)+nx(2,ig)
	njlg(ig)=ny(1,ig)+ny(2,ig)
	nklg(ig)=nz(1,ig)+nz(2,ig)
	
	read (1,*)
	read (1,*) uinf(ig), vinf(ig), winf(ig)
        print *,'>>Block #', ig,'<<'
	print *,'Uinf, Vinf, Winf ',
     &              uinf(ig), vinf(ig), winf(ig)
        read (1,*)
	read (1,*) pinf(ig)
        print *,'Initial pressure',pinf(ig)
       enddo
       
c  ----------------------------------
c Writing the file for restart (che vor dì??)
c
c nt = n step Iterazione)
c rtime = instante di calcolo
c
c questo xche il file rstore è il source per rstart che è il file usato dal
c   solver per riniziare il calcolo 
c ----------------------------

      nt=0
      rtime=0.
      write(22)nt, rtime

c ------------------------------------------
c Read metrics;
c data used are contravariant component of velocity to
c calculate fluxes from physical velocities
c
c  ncx obtained from requested number of cels +2 (phantom) 
c
c -----------------------------------------

      do 100 ig=1,lg
        ncx=nilg(ig)+2
        ncy=njlg(ig)+2
        ncz=nklg(ig)+2
      print *,'Total cells',ncx,ncy,ncz

c ---------------------------
c efettua la lettura dei nodi di griglia
c ---------------------------

      read(12)((((xp(i,j,k,l), i = 1,ncx),
     &              j = 1,ncy), k = 1,ncz), l = 1,3)

c --------------------------
c acquisisce le componenti dei vettori controvarianti
c e il volume
c---------------------------

      
      read(12)(((xix(i,j,k), i = 1,ncx), j = 1,ncy),
     &                        k = 1,ncz)
      read(12)(((xiy(i,j,k), i = 1,ncx), j = 1,ncy),
     &                        k = 1,ncz)
      read(12)(((xiz(i,j,k), i = 1,ncx), j = 1,ncy),
     &                        k = 1,ncz)

      read(12)(((etx(i,j,k), i = 1,ncx), j = 1,ncy),
     &                        k = 1,ncz)
      read(12)(((ety(i,j,k), i = 1,ncx), j = 1,ncy),
     &                        k = 1,ncz)
      read(12)(((etz(i,j,k), i = 1,ncx), j = 1,ncy),
     &                        k = 1,ncz)

      read(12)(((ztx(i,j,k), i = 1,ncx), j = 1,ncy),
     &                        k = 1,ncz)
      read(12)(((zty(i,j,k), i = 1,ncx), j = 1,ncy),
     &                        k = 1,ncz)
      read(12)(((ztz(i,j,k), i = 1,ncx), j = 1,ncy),
     &                        k = 1,ncz)

      read(12)(((jac(i,j,k), i = 1,ncx), j = 1,ncy),
     &                        k = 1,ncz)
      read(12)(((dummy2(i,j,k), i = 1,ncx),j = 1,ncy), 
     &                        k = 1,ncz)
      read(12)(((dummy2(i,j,k), i = 1,ncx),j = 1,ncy), 
     &                        k = 1,ncz)
      read(12)(((dummy2(i,j,k), i = 1,ncx),j = 1,ncy), 
     &                        k = 1,ncz)


c ----------------------
c Setting the velocity field constant equal to the farfield velocity
c and pressure equal to p inf
c  and compute bnigger and smaller cell
c --------------------------
      rjacmax=0.
      rjacmin=1000.
      do k=1,ncz
      do j=1,ncy
      do i=1,ncx
        u(i,j,k,1)=uinf(ig)
        u(i,j,k,2)=vinf(ig)
        u(i,j,k,3)=winf(ig)
        p(i,j,k)=pinf(ig)
        
	rjacmax=max(dabs(jac(i,j,k)),rjacmax)
        rjacmin=min(dabs(jac(i,j,k)),rjacmin)
      enddo
      enddo
      enddo

      print *,'Bigger & Smaller cell'
      print *,'JACMAX=', rjacmax, '(Gianni); JACMIN=', 
     &                               rjacmin,' (Pinotto) '
      print *,''
c -------------------------------------
c Inizializzazione dei flussi controvarianti
c su tutte le celle
c -------------------------------------

      do k=1,ncz
      do j=1,ncy
      do i=1,ncx
        uxi(i,j,k)=(xix(i,j,k)*u(i,j,k,1)+xiy(i,j,k)*
     &           u(i,j,k,2)+xiz(i,j,k)*u(i,j,k,3))
        uej(i,j,k)=(etx(i,j,k)*u(i,j,k,1)+ety(i,j,k)*
     &           u(i,j,k,2)+etz(i,j,k)*u(i,j,k,3))
        uzk(i,j,k)=(ztx(i,j,k)*u(i,j,k,1)+zty(i,j,k)*
     &           u(i,j,k,2)+ztz(i,j,k)*u(i,j,k,3))
      enddo
      enddo
      enddo


c ------------------------
c Setting the dummy variables to zero
c ------------------------
      do k=1,ncz
      do j=1,ncy
      do i=1,ncx
       do n=1,3
        dummy(i,j,k,n)=0.
       enddo
       dummy2(i,j,k)=0.
      enddo
      enddo
      enddo
      nijk=max0(ncx,ncy,ncz) 
      
      do i=1,nijk
      do j=1,nijk
      do id=1,4
      do jd=1,6
        dummy3(i,j,id,jd)=0.
      enddo
      enddo
      enddo
      enddo

c---------------------------------
c Writing the file for restart
c Restart file puo' essere un file fittizio, come questo creato 
c all'inizio della simulazione
c---------------------------------
         print *,'B sure that 2nd line in datainS in', ncx,ncy,ncz
         print*,''
	 
         write(22)((((u(i,j,k,l),i=1,ncx),j=1,ncy),k=1,ncz),l=1,3)
         write(22)((((dummy(i,j,k,l),i=1,ncx),j=1,ncy),k=1,ncz),l=1,3)
         write(22)(((uxi(i,j,k),i=1,ncx),j=1,ncy),k=1,ncz)
         write(22)(((uej(i,j,k),i=1,ncx),j=1,ncy),k=1,ncz)
         write(22)(((uzk(i,j,k),i=1,ncx),j=1,ncy),k=1,ncz)
c dp
         write(22)(((p(i,j,k),i=1,ncx),j=1,ncy),k=1,ncz)
c s
         write(22)(((dummy2(i,j,k),i=1,ncx),j=1,ncy),k=1,ncz)
c hbs

         write(22)(((dummy2(i,j,k),i=1,ncx),j=1,ncy),k=1,ncz)

c usb
         nijk=max0(ncx,ncy,ncz) 
	 
	 
	 write(22)((((dummy3(i,j,k,l),i=1,nijk),j=1,nijk)
     &                                             ,k=1,4),l=1,6)
c vst
         write(22)(((dummy2(i,j,k),i=1,ncx),j=1,ncy),k=1,ncz)
c
 100   CONTINUE

       stop
       end   

c--------------------------------------------------
c                a Gino's production
c
c   creiamo la griglia, ma con la subroutine 
c   occhio alle common!
c-----------------------------------------------------

      program griglia
      implicit double precision (a-h,o-z)     
      include 'main.param'
      parameter( nx2=2*il+1, ny2=2*jl+1, nz2=2*kl+1 )
      integer iv(ng), jv(ng), kv(ng),lg
      integer ncx,ncy,ncz
c      integer nby(ng),ngcouple(6,ng,3)
c      integer icord(ng)
c      double precision cp(31,ncm,ncm,3,6,ng)
      double precision x2(ng,nx2,ny2,nz2,3)   
      double precision xix(il,jl,kl), etx(il,jl,kl),
     &                 xiy(il,jl,kl), ety(il,jl,kl),
     &                 xiz(il,jl,kl), etz(il,jl,kl),
     &                 ztx(il,jl,kl), zty(il,jl,kl),
     &                 ztz(il,jl,kl), jac(il,jl,kl),
     &                 g11(il,jl,kl), g22(il,jl,kl),
     &                 g33(il,jl,kl)

      
c   main.param contiene i valori massimi per il problema proposto
c  lg numero di blocchi
c  iv,jv,kv numero totale di vertici all inclusive (vertici+nodi di celle normali e
c     ghost)
c  nx2 numero massimo di vertici su x      
c
c xix e varie sono i vettori controvarianti, jac è il det della jacob (vol
c  cella), g11 e varie sono i coseni direttori 
c nby,ngcouple,icord,cp sono copiati in maniera pecoreccia per la subroutine che 
c   calcola i coeff x il multiblocco

      
      common /met1/x2,iv,jv,kv,lg
      common /met2/xix,xiy,xiz,etx,ety,etz,ztx,zty,ztz,
     & jac,g11,g22,g33
c      common /met3/nby,ngcouple,cp,icord
     
          
      open (11,file='xyu',form='unformatted')
      open (12,file='metric',form='unformatted')
      open (14,file = 'met',  form = 'unformatted')
      open (21,file='usrckmetric')
      open (33,file='usrcell')
      
      
        
c ----------
c 33 file formattato per controllo dell'utente sunedit lla griglia
c 12 è il file che ha dentro i vettori controvarianti (è scritto
c    da metric
c 21 file di controllo valori max coseni direttori metricout
c ----------
     
      print *,'   '
      print *,'    --------------------------'
      print *,'    >> a Gino''s production <<'
      print *,'    --------------------------' 
      print *,''
      print *,'Attension, questo programa fa la griglia uniforme retta'
      print *,'per averla in coord cilindriche modificare icord a 2'
      print *,'e ricompilare'




c -----------------
c  evoco coordmkr
c ---------------------

      call coordmkr
      print *,' Start computing for ', lg,' block(s)'
      print *,''

c ----------------------
c output di tutta la cella su file cell
c -----------------------------------
      do 1 ill=1,lg
      do 1 i=1,iv(ill)
      do 1 j=1,jv(ill)
      do 1 k=1,kv(ill)
c        write(33,*),i,j,k
	write(33,*)(x2(ill,i,j,k,l),l=1,3)
1     continue

      write(11)lg

c  ------------------------------------
c  calcolo dei tensori di metrica
c  scrivo x2 (solo i centri cella) in metric, il solutore li legge così...
c   evoco metric
c  scrivo metric e xyu
c  
c  ------------------------------------

      do 2 ig=1,lg
      
       ncx=(iv(ig)-1)/2
       ncy=(jv(ig)-1)/2
       ncz=(kv(ig)-1)/2

       write(12)((((x2(ig,2*i,2*j,2*k,l),i=1,ncx),j=1,ncy),
     &                                        k=1,ncz),l=1,3)
       
       write(11)ncx,ncy,ncz
       write(11)((((x2(ig,2*i,2*j,2*k,l),i=1,ncx),j=1,ncy),
     &                                        k=1,ncz),l=1,3)
       
       call metrica(ig)
	
c --------------
c  scrittura file metric x il solutore
c   prima era in metrica.f ..
c --------------
	
      write(12)(((xix(i,j,k), i = 1,ncx), j = 1,ncy),
     &                        k = 1,ncz)
      write(12)(((xiy(i,j,k), i = 1,ncx), j = 1,ncy),
     &                        k = 1,ncz)
      write(12)(((xiz(i,j,k), i = 1,ncx), j = 1,ncy),
     &                        k = 1,ncz)

      write(12)(((etx(i,j,k), i = 1,ncx), j = 1,ncy),
     &                        k = 1,ncz)
      write(12)(((ety(i,j,k), i = 1,ncx), j = 1,ncy),
     &                        k = 1,ncz)
      write(12)(((etz(i,j,k), i = 1,ncx), j = 1,ncy),
     &                        k = 1,ncz)

      write(12)(((ztx(i,j,k), i = 1,ncx), j = 1,ncy),
     &                        k = 1,ncz)
      write(12)(((zty(i,j,k), i = 1,ncx), j = 1,ncy),
     &                        k = 1,ncz)
      write(12)(((ztz(i,j,k), i = 1,ncx), j = 1,ncy),
     &                        k = 1,ncz)
   

      write(12)(((jac(i,j,k), i = 1,ncx), j = 1,ncy),
     &                        k = 1,ncz)

      write(12)(((g11(i,j,k), i = 1,ncx), j = 1,ncy),
     &                        k = 1,ncz)
      write(12)(((g22(i,j,k), i = 1,ncx), j = 1,ncy),
     &                        k = 1,ncz)
      write(12)(((g33(i,j,k), i = 1,ncx), j = 1,ncy),
     &                        k = 1,ncz)

c -----------------------------------
c output su file "met", commentato xche nn si sa bene ancora a cosa serve..
c
c ------------------------------------

c
c
c
      write(14)(((jac(i,j,k), i = 1,ncx), j = 1,ncy),
     &          k = 1,ncz)
      write(14)(((xix(i,j,k), i = 1,ncx), j = 1,ncy),
     &          k = 1,ncz)
      write(14)(((xiy(i,j,k), i = 1,ncx), j = 1,ncy),
     &          k = 1,ncz)
      write(14)(((xiz(i,j,k), i = 1,ncx), j = 1,ncy),
     &          k = 1,ncz)
      write(14)(((etx(i,j,k), i = 1,ncx), j = 1,ncy),
     &          k = 1,ncz)
      write(14)(((ety(i,j,k), i = 1,ncx), j = 1,ncy),
     &          k = 1,ncz)
      write(14)(((etz(i,j,k), i = 1,ncx), j = 1,ncy),
     &          k = 1,ncz)
      write(14)(((ztx(i,j,k), i = 1,ncx), j = 1,ncy),
     &          k = 1,ncz)
      write(14)(((zty(i,j,k), i = 1,ncx), j = 1,ncy),
     &          k = 1,ncz)
      write(14)(((ztz(i,j,k), i = 1,ncx), j = 1,ncy),
     &          k = 1,ncz)







2     continue


c ------------------------
c  questa parte commentata soprintende al calcolo dei coeffcienti
c   per l'interpolazione in caso di griglie multiblocco
c -----------------------


c      if(lg.gt.1)then
c        print *,' calculating interboundary coefficients'
c        open(22,file = 'coef.in',form='unformatted')
c        call coefint
c      endif

      stop
         end














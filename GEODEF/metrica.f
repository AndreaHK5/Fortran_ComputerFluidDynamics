c-----------------------------------------------------------------------------
c    
c     subroutine per il calcolo di tensori controvarianti e porcellate varie
c     (che poi sinceramente preferivo la griglia, magari anche
c     con un due salCiccie)
c
c
c  attenzione, c'è anche la subroutine molecule!!!!!!!!!!!
c-----------------------------------------------------
      
      subroutine metrica(ig)
      implicit double precision (a-h,o-z)
      include 'main.param'
      parameter(nx2=2*il+1, ny2=2*jl+1, nz2=2*kl+1 )
      integer iv(ng), jv(ng), kv(ng),lg
      double precision x2(ng,nx2,ny2,nz2,3)
      double precision xix(il,jl,kl), etx(il,jl,kl),
     &                 xiy(il,jl,kl), ety(il,jl,kl),
     &                 xiz(il,jl,kl), etz(il,jl,kl),
     &                 ztx(il,jl,kl), zty(il,jl,kl),
     &                 ztz(il,jl,kl), jac(il,jl,kl),
     &                 g11(il,jl,kl), g22(il,jl,kl),
     &                 g33(il,jl,kl)
      integer ncx,ncy,ncz
      double precision resci(il,jl,kl), rescj(il,jl,kl),
     &                 resck(il,jl,kl)
      
      common /met1/x2,iv,jv,kv,lg
      common /met2/xix,xiy,xiz,etx,ety,etz,ztx,zty,ztz,
     & jac,g11,g22,g33


c----------------------------------------------------------
c ncx, ncy,ncz numero di celle (iv è il numero totale di vertici di x2) 
c----------------------------------------------------------

      ncx=(iv(ig)-1)/2
      ncy=(jv(ig)-1)/2
      ncz=(kv(ig)-1)/2

c --------------------------
c inizializzazione indici per ricerca max e min
c ---------------------------------

      ajacm = 1000000.d0
      iMjac = 0
      jMjac = 0
      kMjac = 0
      imjac = 0
      jmjac = 0
      kmjac = 0
      g12m = 0.d0
      g13m = 0.d0
      g23m = 0.d0
      i12  = 1
      j12  = 1
      k12  = 1
      i13  = 1
      j13  = 1
      k13  = 1
      i23  = 1
      j23  = 1
      k23  = 1
      
c ----------------------
c per ogni cella 
c evochioamo molecule su faccia x, y, z e su centro cella
c 
c troviamo le celle con max volume, max coseni direttori nelle tre direzioni
c ---------------------------------

      do 2,i=1,ncx
      do 2,j=1,ncy
      do 2,k=1,ncz
      call molecule(ig,i,j,k,1,g12t,g13t,g23t)
      call molecule(ig,i,j,k,2,g12t,g13t,g23t)
      call molecule(ig,i,j,k,3,g12t,g13t,g23t)
      call molecule(ig,i,j,k,4,g12t,g13t,g23t)

c -----------------------------
c controllo sui volumi di cella (jac nei centro cella) massimo e minimo
c e check su orthogonality of the grid
c -------------------------------

      if(dabs(jac(i,j,k)).lt.ajacm) then
            ajacm = dabs(jac(i,j,k))
            imjac  = i
            jmjac  = j
            kmjac  = k
       endif
       if(dabs(jac(i,j,k)).lt.1.d-12) then
            iMjac = i
            jMjac = j
            kMjac = k
        endif
c-----------------------------------------
c coseni direttori
c-----------------------------------------
      if(dabs(g12t).gt.g12m) then
            g12m = dabs(g12t)
            i12  = i
            j12  = j
            k12  = k
       endif
      if(dabs(g13t).gt.g13m) then
            g13m = dabs(g13t)
            i13  = i
            j13  = j
            k13  = k
       endif
      if(dabs(g23t).gt.g23m) then
            g23m = dabs(g23t)
            i23  = i
            j23  = j
            k23  = k
       endif

 2     continue

c ----------------
c output su file checkmetric dei massimi valori dei coseni direttori
c ----------------

      write(21,294)
      write(21,*)'MAX',ajacm,imjac,jmjac,kmjac

      write(21,295)g12m,i12,j12,k12
      write(21,296)g13m,i13,j13,k13
      write(21,297)g23m,i23,j23,k23
 294  format('Check the orthogonality of the grid',/)

 295  format(1x,'g12max = ',e10.4,' i = ',i3,
     &' j = ',i3,' k = ',i3,/)
 296  format(1x,'g13max = ',e10.4,' i = ',i3,
     &' j = ',i3,' k = ',i3,/)
 297  format(1x,'g23max = ',e10.4,' i = ',i3,
     &' j = ',i3,' k = ',i3,/)




c------------------------------------------------------
c    scrittura file output di metrica su 12
c 
c    output xix, ..., on faces and jac at p nodes
c------------------------------------------------------

      tm = 0.d0
      im = 1
      jm = 1
      km = 1
      do 400 i = 2,ncx-1
      do 400 j = 2,ncy-1
      do 400 k = 2,ncz-1
         resci(i,j,k) = xix(i,j,k) - xix(i-1,j,k)   +
     &                  etx(i,j,k) - etx(i,j-1,k)   +
     &                  ztx(i,j,k) - ztx(i,j,k-1)
         rescj(i,j,k) = xiy(i,j,k) - xiy(i-1,j,k)   +
     &                  ety(i,j,k) - ety(i,j-1,k)   +
     &                  zty(i,j,k) - zty(i,j,k-1)
         resck(i,j,k) = xiz(i,j,k) - xiz(i-1,j,k)   +
     &                  etz(i,j,k) - etz(i,j-1,k)   +
     &                  ztz(i,j,k) - ztz(i,j,k-1)
            t1 = resci(i,j,k)
            t2 = rescj(i,j,k)
            t3 = resck(i,j,k)
            t  = max(t1,t2,t3)
            if(dabs(t).gt.dabs(tm)) then
               tm = t
               im = i
               jm = j
               km = k
            endif
 400        continue
         
      print *,'max deviation of cells:',tm
      print *,'@ cell i=',im,' j=',jm,' k= ',km
      print *,''         

	 write(21,2111)ig,tm,im,jm,km
 2111 format(/,'Check closedness of cells',/,'ig = ',i3,
     &       /,' max. dev. = ',d10.4,' at i = ',i3,' j = ',i3,
     &      ' k = ',i3)
         if(dabs(tm).gt.1.e-12) then
            print *,'tm = ',tm,' check file mout'
         endif
      return
      end

c --------------------------------
c   output of metrics at p for plotting
c   tenere anche questo output??
c   Questo output è stato spostato al programma principale grilgia.f
c  ------------------------------

c ------------
c  qui c'era  anche un output di x2 che è stato spostato nel programma principale
c --------














         subroutine molecule(ig,ici,icj,ick,icomp,g12t,g13t,g23t)


c    input: ig=numero blocco
c    ici,icy,ick=numero di cella
c    icomp= tipo di conto 


      implicit double precision (a-h,o-z)
      include 'main.param'
      parameter(nx2=2*il+1, ny2=2*jl+1, nz2=2*kl+1 )
      integer iv(ng), jv(ng), kv(ng),lg
      double precision x2(ng,nx2,ny2,nz2,3)
      double precision xix(il,jl,kl), etx(il,jl,kl),
     &                 xiy(il,jl,kl), ety(il,jl,kl),
     &                 xiz(il,jl,kl), etz(il,jl,kl),
     &                 ztx(il,jl,kl), zty(il,jl,kl),
     &                 ztz(il,jl,kl), jac(il,jl,kl),
     &                 g11(il,jl,kl), g22(il,jl,kl),
     &                 g33(il,jl,kl)
      double precision xxi,xet,xzt,yxi,yet,yzt,zxi,zet,zzt
      double precision i1,j1,k1
      
      common /met1/x2,iv,jv,kv,lg
      common /met2/xix,xiy,xiz,etx,ety,etz,ztx,zty,ztz,
     & jac,g11,g22,g33


c---------------------------------------------------
c 
c icomp -> tipo di calcolo da fare:
c icomp=1 nodo u
c icomp=2 nodo v
c icomp=3 nodo w
c icomp=4 centro cella
c 
c aggiornamento relativo degli indici di calcolo per prendere i
c   punti giusti in x2
c
c---------------------------------------------------
         
	 if(icomp.eq.1) then
         i1=2*ici+1
         j1=2*icj
         k1=2*ick
         endif
         if(icomp.eq.2) then
         i1=2*ici
         j1=2*icj+1
         k1=2*ick
         endif
         if(icomp.eq.3) then
         i1=2*ici
         j1=2*icj
         k1=2*ick+1
         endif
         if(icomp.eq.4) then
         i1=2*ici
         j1=2*icj
         k1=2*ick
         endif
   

c calcolo degli elementi di [J^-1] per punto centrale i1,j1,k1
c che viene posto al centro della faccia in direzione x, y o z se icomp=1 2 o 3
c o al centro della cella se icomp=4

         xxi = x2(ig,i1+1,j1,k1,1) - x2(ig,i1-1,j1,k1,1)
         xet = x2(ig,i1,j1+1,k1,1) - x2(ig,i1,j1-1,k1,1)
         xzt = x2(ig,i1,j1,k1+1,1) - x2(ig,i1,j1,k1-1,1)
         yxi = x2(ig,i1+1,j1,k1,2) - x2(ig,i1-1,j1,k1,2)
         yet = x2(ig,i1,j1+1,k1,2) - x2(ig,i1,j1-1,k1,2)
         yzt = x2(ig,i1,j1,k1+1,2) - x2(ig,i1,j1,k1-1,2)
         zxi = x2(ig,i1+1,j1,k1,3) - x2(ig,i1-1,j1,k1,3)
         zet = x2(ig,i1,j1+1,k1,3) - x2(ig,i1,j1-1,k1,3)
         zzt = x2(ig,i1,j1,k1+1,3) - x2(ig,i1,j1,k1-1,3)

c ---------------------
c 
c xi = vettore c v lungo i
c et = vettore c v lungo j
c zt = vettore c v lungo k
c 
c  e adesso ci facciamo il determinante della jacob
c
c ----------------------

         jac(ici,icj,ick) = xxi*yet*zzt + xet*yzt*zxi + xzt*zet*yxi
     &                - xzt*yet*zxi - xet*yxi*zzt - xxi*zet*yzt


         if(abs(jac(ici,icj,ick)).lt.error) then
            temp = 0.
         else
            temp = 1./jac(ici,icj,ick)
         endif

c-----------------------------------------------
c componenti del vettore controvariante faccia lungo i 
c    quindi uso le componenti dei vettori lungo j e k (et, zt)
c inoltre g11 è il quadrato del modulo del vettore controvariante lungo i
c    diviso per il volume della cella
c-----------------------------------------------

         if(icomp.eq.1) then
         xix(ici,icj,ick) =   yet*zzt - yzt*zet
         xiy(ici,icj,ick) = -(xet*zzt - xzt*zet)
         xiz(ici,icj,ick) =   xet*yzt - xzt*yet
         g11(ici,icj,ick) = temp*(xix(ici,icj,ick)*xix(ici,icj,ick)
     &                    + xiy(ici,icj,ick)*xiy(ici,icj,ick)
     &                    + xiz(ici,icj,ick)*xiz(ici,icj,ick))
         return
         endif

c-----------------------------------------------
c componenti del vettore controvariante faccia y normal
c-----------------------------------------------

         if(icomp.eq.2)then
         etx(ici,icj,ick) = -(yxi*zzt - yzt*zxi)
         ety(ici,icj,ick) =   xxi*zzt - xzt*zxi
         etz(ici,icj,ick) = -(xxi*yzt - xzt*yxi)
         g22(ici,icj,ick) = temp*(etx(ici,icj,ick)*etx(ici,icj,ick)
     &                    + ety(ici,icj,ick)*ety(ici,icj,ick)
     &                    + etz(ici,icj,ick)*etz(ici,icj,ick))
         return
         endif

c-----------------------------------------------
c componenti del vettore controvariante faccia z normal
c-----------------------------------------------

         if(icomp.eq.3)then
         ztx(ici,icj,ick) =   yxi*zet - yet*zxi
         zty(ici,icj,ick) = -(xxi*zet - xet*zxi)
         ztz(ici,icj,ick) =   xxi*yet - xet*yxi
         g33(ici,icj,ick) = temp*(ztx(ici,icj,ick)*ztx(ici,icj,ick)
     &                    + zty(ici,icj,ick)*zty(ici,icj,ick)
     &                    + ztz(ici,icj,ick)*ztz(ici,icj,ick))
         return
         endif
         
	 if(icomp.eq.4)then

c ----------------------------
c inversa di [J-1]: colonne =vettori controvarianti
c 
c   compute g12, g13, g23
c
c ---------------------------------

         xixt =   yet*zzt - yzt*zet
         xiyt = -(xet*zzt - xzt*zet)
         xizt =   xet*yzt - xzt*yet
         etxt = -(yxi*zzt - yzt*zxi)
         etyt =   xxi*zzt - xzt*zxi
         etzt = -(xxi*yzt - xzt*yxi)
         ztxt =   yxi*zet - yet*zxi
         ztyt = -(xxi*zet - xet*zxi)
         ztzt =   xxi*yet - xet*yxi

c-----------------------------
c g12 g13 g23 coseni direttori -> grado di distorsione della mesh rispetto
c a mesh ortogonale di partenza
c-----------------------------

         if(g11(ici,icj,ick)*g22(ici,icj,ick)*g33(ici,icj,ick).eq.0)then
         continue
         else
         g12u = temp*(xixt*etxt + xiyt*etyt + xizt*etzt)/
     &           g11(ici,icj,ick)
         g12v = temp*(xixt*etxt + xiyt*etyt + xizt*etzt)/
     &           g22(ici,icj,ick)
         g13t = temp*(xixt*ztxt + xiyt*ztyt + xizt*ztzt)/
     &           g33(ici,icj,ick)
         g23t = temp*(etxt*ztxt + etyt*ztyt + etzt*ztzt)/
     &           g33(ici,icj,ick)
         g12t = max(g12u, g12v)
         endif
         endif
         
	 return
         end   




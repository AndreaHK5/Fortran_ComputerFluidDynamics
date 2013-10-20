       subroutine velxyz(P,V)



      implicit double precision (a-h,o-z)
      include 'main.param'

c----------------------------------------------------------------------
c  gli passi un vettore di posizione e ti torna la velocità del fluido
c   in quel posto (figo eh?)
c
c   Subroutine già debugata e controllata, questa è roba buona mica 
c     coi conservti!!!
c----------------------------------------------------------------------
      double precision P(3),V(3)
      double precision G(il,jl,kl,3)
      double precision A(il,jl,kl,3)     
      double precision V1(3),V2(3),V3(3),V4(3),V5(3),V6(3),V7(3),V8(3)
      double precision V12(3),V34(3),V56(3),V78(3)
      double precision V13(3),V57(3),V24(3),V68(3)
      double precision V1234(3),V5678(3)
      double precision V1357(3),V2468(3)

      integer ilg(ng),jlg(ng),klg(ng)
      integer index_i,index_j,index,k
      common /general/G,A,ilg,jlg,klg

c ----------------------------------------------- 
c  Trovo gli indici della cela a cui apaprtiene il punto in questione
c -----------------------------------------------


      index_i=1
      index_j=1
      index_k=1
      
      do i=1,ilg(1)
       j=1
       k=1

       if ((P(1).GE.G(i,j,k,1)).AND.
     &    (P(1).LT.G(i+1,j,k,1))) then
         index_i=i
        goto 1
        endif
        enddo
 1       continue

        do j=1,jlg(1)
        k=1
        if ((P(2).GE.G(index_i,j,k,2)).AND.
     &     (P(2).LT.G(index_i,j+1,k,2))) then
         index_j=j
         goto 2
         endif
         enddo
 2       continue

         do k=1,klg(1)
         if ((P(3).GE.G(index_i,index_j,k,1)).AND.
     &      (P(3).LT.G(index_i,index_j,k+1,1))) then
         index_k=k
         goto 3
         endif
         enddo
3        continue



c ---------------------------
c Placing velocity on cell points
c --------------------------      

       V1(1)=A(index_i,index_j,index_k,1)
       V1(2)=A(index_i,index_j,index_k,2)
       V1(3)=A(index_i,index_j,index_k,3)
       

       V2(1)=A(index_i,index_j+1,index_k,1)
       V2(2)=A(index_i,index_j+1,index_k,2)
       V2(3)=A(index_i,index_j+1,index_k,3)
  
       V3(1)=A(index_i+1,index_j,index_k,1)
       V3(2)=A(index_i+1,index_j,index_k,2)
       V3(3)=A(index_i+1,index_j,index_k,3)

       V4(1)=A(index_i+1,index_j+1,index_k,1)
       V4(2)=A(index_i+1,index_j+1,index_k,2)
       V4(3)=A(index_i+1,index_j+1,index_k,3)
       

       V5(1)=A(index_i,index_j,index_k+1,1)
       V5(2)=A(index_i,index_j,index_k+1,2)
       V5(3)=A(index_i,index_j,index_k+1,3)

       V6(1)=A(index_i,index_j+1,index_k+1,1)
       V6(2)=A(index_i,index_j+1,index_k+1,2)
       V6(3)=A(index_i,index_j+1,index_k+1,3)

       V7(1)=A(index_i+1,index_j,index_k+1,1)
       V7(2)=A(index_i+1,index_j,index_k+1,2)
       V7(3)=A(index_i+1,index_j,index_k+1,3)

       V8(1)=A(index_i+1,index_j+1,index_k+1,1)
       V8(2)=A(index_i+1,index_j+1,index_k+1,2)
       V8(3)=A(index_i+1,index_j+1,index_k+1,3)



c ------------------------------------
c  Now find the interpolation on y basis
c ------------------------------------
       V12(1)=V1(1)+(V2(1)-V1(1))*((P(2)-
     &   G(index_i,index_j,index_k,2))/
     & (G(index_i,index_j+1,index_k,2)-
     & G(index_i,index_j,index_k,2)))
       
      V12(2)=V1(2)+(V2(2)-V1(2))*((P(2)-
     &  G(index_i,index_j,index_k,2))/
     & (G(index_i,index_j+1,index_k,2)-
     & G(index_i,index_j,index_k,2)))

       V12(3)=V1(3)+(V2(3)-V1(3))*((P(2)-
     &   G(index_i,index_j,index_k,2))/
     & (G(index_i,index_j+1,index_k,2)-
     & G(index_i,index_j,index_k,2)))

     
    
       V34(1)=V3(1)+(V4(1)-V3(1))*((P(2)-
     & G(index_i+1,index_j,index_k,2))/
     & (G(index_i+1,index_j+1,index_k,2)-
     & G(index_i+1,index_j,index_k,2)))
     
      V34(2)=V3(2)+(V4(2)-V3(2))*((P(2)-
     & G(index_i+1,index_j,index_k,2))/
     & (G(index_i+1,index_j+1,index_k,2)-
     & G(index_i+1,index_j,index_k,2)))

       V34(3)=V3(3)+(V4(3)-V3(3))*((P(2)-
     & G(index_i+1,index_j,index_k,2))/
     & (G(index_i+1,index_j+1,index_k,2)-
     & G(index_i+1,index_j,index_k,2)))
     


      V56(1)=V5(1)+(V6(1)-V5(1))*((P(2)-
     & G(index_i,index_j,index_k+1,2))/
     & (G(index_i,index_j+1,index_k+1,2)-
     & G(index_i,index_j,index_k+1,2)))

       V56(2)=V5(2)+(V6(2)-V5(2))*((P(2)-
     & G(index_i,index_j,index_k+1,2))/
     & (G(index_i,index_j+1,index_k+1,2)-
     & G(index_i,index_j,index_k+1,2)))
 
       V56(3)=V5(3)+(V6(3)-V5(3))*((P(2)-
     & G(index_i,index_j,index_k+1,2))/
     & (G(index_i,index_j+1,index_k+1,2)-
     & G(index_i,index_j,index_k+1,2)))



       V78(1)=V7(1)+(V8(1)-V7(1))*((P(2)-
     & G(index_i+1,index_j,index_k+1,2))/
     & (G(index_i+1,index_j+1,index_k+1,2)-
     & G(index_i+1,index_j,index_k+1,2)))

       V78(2)=V7(2)+(V8(2)-V7(2))*((P(2)-
     & G(index_i+1,index_j,index_k+1,2))/
     & (G(index_i+1,index_j+1,index_k+1,2)-
     & G(index_i+1,index_j,index_k+1,2)))

       V78(3)=V7(3)+(V8(3)-V7(3))*((P(2)-
     & G(index_i+1,index_j,index_k+1,2))/
     & (G(index_i+1,index_j+1,index_k+1,2)-
     & G(index_i+1,index_j,index_k+1,2)))

c ------------------------------------
c  Interpolation on x basis
c ------------------------------------


       V1234(1)=V12(1)+(V34(1)-V12(1))*((P(1)-
     & G(index_i,index_j,index_k,1))/
     & (G(index_i+1,index_j,index_k,1)-
     & G(index_i,index_j,index_k,1)))

       V1234(2)=V12(2)+(V34(2)-V12(2))*((P(1)-
     & G(index_i,index_j,index_k,1))/
     & (G(index_i+1,index_j,index_k,1)-
     & G(index_i,index_j,index_k,1)))

       V1234(3)=V12(3)+(V34(3)-V12(3))*((P(1)-
     & G(index_i,index_j,index_k,1))/
     & (G(index_i+1,index_j,index_k,1)-
     & G(index_i,index_j,index_k,1)))



       V5678(1)=V56(1)+(V78(1)-V56(1))*((P(1)-
     & G(index_i,index_j,index_k+1,1))/
     & (G(index_i+1,index_j,index_k+1,1)-
     & G(index_i,index_j,index_k+1,1)))

       V5678(2)=V56(2)+(V78(2)-V56(2))*((P(1)-
     & G(index_i,index_j,index_k+1,1))/
     & (G(index_i+1,index_j,index_k+1,1)-
     & G(index_i,index_j,index_k+1,1)))

       V5678(3)=V56(3)+(V78(3)-V56(3))*((P(1)-
     & G(index_i,index_j,index_k+1,1))/
     & (G(index_i+1,index_j,index_k+1,1)-
     & G(index_i,index_j,index_k+1,1)))

c ---------------------------------
c  Finally get the actual velocity in point P via z interpolation
c ----------------------------------


       V(1)=V1234(1)+(V5678(1)-V1234(1))*((P(3)-
     & G(index_i,index_j,index_k,3))/
     & (G(index_i,index_j,index_k+1,3)-
     & G(index_i,index_j,index_k,3)))

       V(2)=V1234(2)+(V5678(2)-V1234(2))*((P(3)-
     & G(index_i,index_j,index_k,3))/
     & (G(index_i,index_j,index_k+1,3)-
     & G(index_i,index_j,index_k,3)))

       V(3)=V1234(3)+(V5678(3)-V1234(3))*((P(3)-
     & G(index_i,index_j,index_k,3))/
     & (G(index_i,index_j,index_k+1,3)-
     & G(index_i,index_j,index_k,3)))


c      open(unit =27,file='pvp_00500',form='formatted')

c      write(27,*) (P(l),l=1,3)
c      write(27,*)  (G(index_i,index_j,index_k,l),l=1,3) 
c      write(27,*)  (G(index_i,index_j+1,index_k,l),l=1,3)
c      write(27,*)  (G(index_i+1,index_j,index_k,l),l=1,3)
c      write(27,*)  (G(index_i+1,index_j+1,index_k,l),l=1,3)
c      write(27,*)  (G(index_i,index_j,index_k+1,l),l=1,3) 
c      write(27,*)  (G(index_i,index_j+1,index_k+1,l),l=1,3)
c      write(27,*)  (G(index_i+1,index_j,index_k+1,l),l=1,3)
c      write(27,*)  (G(index_i+1,index_j+1,index_k+1,l),l=1,3)
c      write(27,*)  (V1(l),l=1,3)
c      write(27,*)  (V2(l),l=1,3)
c      write(27,*)  (V3(l),l=1,3)
c      write(27,*)  (V4(l),l=1,3)
c      write(27,*)  (V5(l),l=1,3)
c      write(27,*)  (V6(l),l=1,3)
c      write(27,*)  (V7(l),l=1,3)
c      write(27,*)  (V8(l),l=1,3)
c      write(27,*)
c      write(27,*)  (V12(l),l=1,3)
c      write(27,*)  (V34(l),l=1,3)
c      write(27,*)  (V56(l),l=1,3)
c      write(27,*)  (V78(l),l=1,3)
c      write(27,*)
c      write(27,*)  (V1234(l),l=1,3)
c      write(27,*)  (V5678(l),l=1,3)
c      write(27,*)
c      write(27,*) (P(t),t=1,3), (V(h),h=1,3)
c      write(27,*)

      return
      end

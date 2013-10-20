      subroutine pospar(P0,P1,U0,U1)
      
    
      double precision P0(3),P1(3),U0(3),U1(3)
      double precision DD,DT,KVF,K02,K04,GAx,GAy,GAz
c      double precision XMIN,XMAX,YMIN,YMAX,ZMIN,ZMAX
c      integer R1,R2,R3,R4,R5,R6
      
      
      common /param/DD,DT,KVF,K02,K04,GAx,GAy,GAz
c      common /rebo/XMAX,XMIN,YMAX,YMIN,ZMAX,ZMIN
     &                                ,R1,R2,R3,R4,R5,R6
c ----------------------------
c Calculation of particle position at the next time step.
c ----------------------------

       do 10 l=1,3
10           P1(l)=P0(l)+DT*.5*(U0(l)+U1(l))


c --------------------------------
c Chiamo rebound per vedere se la particiella esce
c dal dominio e in caso cosa fà...
c --------------------------------


      call rebound(P0,P1,U0,U1)






      return
       end
c

c

      subroutine rebound(P0,P1,U0,U1)     
      
      double precision P0(3),P1(3),U0(3),U1(3)      
      double precision XMIN,XMAX,YMIN,YMAX,ZMIN,ZMAX
      double precision imptime
      integer R1,R2,R3,R4,R5,R6
      integer l
      common /rebo/XMAX,XMIN,YMAX,YMIN,ZMAX,ZMIN
     &                                ,R1,R2,R3,R4,R5,R6
     


c ---------------------------------
c Controlo sulla nuova posizione della particiella, se è ok esco
c ---------------------------------     
     
      if (P1(1).gt.XMIN.and.P1(1).lt.XMAX.and.
     &      P1(2).gt.YMIN.and.P1(2).lt.YMAX.and.
     &      P1(3).gt.ZMIN.and.P1(3).lt.ZMAX) then
      goto 99
      endif

c ----------------------------------
c Se sono registrate violazioni del dominio vado a calcolare nuova posizione
c  e nuova velocità
c ----------------------------------

c -------------       
c Particella che esce dal fondo del dominio
c Prima il caso di fine dominio (la particiella lascia traccia sul bordo
c -------------

      if (P1(1).gt.XMAX) then

       if (R3.eq.0) then
        print *,'particiella uscita da xmax' 
	P1(1)=XMAX 
        imptime=(XMAX-P0(1))*2./(U0(1)+U1(1))
        P1(2)=P0(2)+(U0(2)+U1(2))*imptime/2.
        P1(3)=P0(3)+(U0(3)+U1(3))*imptime/2.
        do l=1,3
	U1(l)=(U1(l)+U0(l))/2.
	enddo
	goto 99
       endif

c -----------------
c Rimbalzo anaelastico
c -----------------
       if (R3.eq.1) then
        print *,'particiella rimbalzata da xmax'
	U1(1)=U1(1)*(-1.)
	P1(1)=XMAX-(P1(1)-XMAX)
        goto 99
       endif 

      endif


c -------------       
c Particella che esce dall'imbocco del dominio
c -------------      
      
      if (P1(1).lt.XMIN) then
       print *,'Rivedere tutto, una particiella esce dall"ingresso??'
       if (R1.eq.0) then
        P1(1)=XMIN 
        imptime=(XMIN-P0(1))*2./(U0(1)+U1(1))
        P1(2)=P0(2)+(U0(2)+U1(2))*imptime/2.
        P1(3)=P0(3)+(U0(3)+U1(3))*imptime/2.
        do l=1,3
	U1(l)=(U1(l)+U0(l))/2.
	enddo
	goto 99
       endif
       
       if (R1.eq.1) then
        U1(1)=U1(1)*(-1.)
	P1(1)=XMIN-(P1(1)-XMIN)
        goto 99
       endif 

      endif

c -----------------
c Out sul bordo Y negativo
c -----------------
      
      if (P1(2).lt.YMIN) then
       if (R2.eq.0) then
        print *,'particiella uscita da ymin'        
	P1(2)=YMIN 
        imptime=(YMIN-P0(2))*2./(U0(2)+U1(2))
        P1(1)=P0(1)+(U0(1)+U1(1))*imptime/2.
        P1(3)=P0(3)+(U0(3)+U1(3))*imptime/2.
        do l=1,3
	U1(l)=(U1(l)+U0(l))/2.
	enddo
	goto 99
       endif
       
       if (R2.eq.1) then
        print *,'particiella rimbalzata da ymin'
	U1(2)=U1(2)*(-1.)
	P1(2)=YMIN-(P1(2)-YMIN)
        goto 99
       endif 

      endif

c -----------------
c Out sul bordo Y positivo
c -----------------
      
      if (P1(2).gt.YMAX) then
       if (R2.eq.0) then 
        print *,'particiella uscita da ymax'
        P1(2)=YMAX
        imptime=(YMAX-P0(2))*2./(U0(2)+U1(2))
        P1(1)=P0(1)+(U0(1)+U1(1))*imptime/2.
        P1(3)=P0(3)+(U0(3)+U1(3))*imptime/2.
        do l=1,3
	U1(l)=(U1(l)+U0(l))/2.
	enddo
	goto 99
       endif
       
       if (R4.eq.1) then        
       print *,'particiella rimbalzata da ymax'
        U1(2)=U1(2)*(-1.)
	P1(2)=YMAX-(P1(2)-YMAX)
        goto 99
       endif 
 
      endif

c -----------------
c Out sul bordo Z positivo
c -----------------
      
      if (P1(3).gt.ZMAX) then
       if (R5.eq.0) then        
       print *,'particiella uscita da zmax'
        P1(3)=ZMAX
        imptime=(ZMAX-P0(3))*2./(U0(3)+U1(3))
        P1(1)=P0(1)+(U0(1)+U1(1))*imptime/2.
        P1(2)=P0(2)+(U0(2)+U1(2))*imptime/2.
        do l=1,3
	U1(l)=(U1(l)+U0(l))/2.
	enddo
	goto 99
       endif
       
       if (R5.eq.1) then       
       print *,'particiella rimbalzata da zmax'
        U1(3)=U1(2)*(-1.)
	P1(3)=ZMAX-(P1(3)-ZMAX)
        goto 99
       endif 

      endif

      if (P1(3).lt.ZMIN) then
       if (R5.eq.0) then       
       print *,'particiella uscita da zmin'
        P1(3)=ZMIN
        imptime=(ZMIN-P0(3))*2./(U0(3)+U1(3))
        P1(1)=P0(1)+(U0(1)+U1(1))*imptime/2.
        P1(2)=P0(2)+(U0(2)+U1(2))*imptime/2.
        do l=1,3
	U1(l)=(U1(l)+U0(l))/2.
	enddo
	goto 99
       endif
       
       if (R5.eq.1) then       
       print *,'particiella uscita da zmin'
        U1(3)=U1(3)*(-1.)
	P1(3)=ZMIN-(P1(3)-ZMIN)
        goto 99
       endif 

      endif


     




 99    return
      end
      
      

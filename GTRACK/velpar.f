      subroutine velpar(U0,V0,U1) 

c ----------------------------------
c A questa subroutine passi vel+pos particiella,velocità fluido
c  e ti dà la velocità aggiornata della particiella
c ----------------------------------

      double precision U0(3),V0(3),U1(3) 
      double precision mod,RV(3),REF,EU(4)
      double precision K03,K05
      double precision DD,DT,KVF,K02,K04,GAx,GAy,GAz
      
      common /param/DD,DT,KVF,K02,K04,GAx,GAy,GAz
c ----------------------------------
c calcolo del modulo della velocità relativa
c ----------------------------------


 60   mod=0.
       do 7 l=1,3
       RV(l)=U0(l)-V0(l)
7       mod=mod+RV(l)*RV(l)
       mod=SQRT(mod)


       REF=mod*DD/KVF
       if (REF.LT.1) K05=1  
       if ((REF.GT.1).AND.(REF.LT.1E3)) K05=(1+0.15*REF**0.687)
       if ((REF.GT.1E3).AND.(REF.LT.1E5)) then
            K05=0.45*REF/24 
       endif
       if (REF.GT.1E5) then 
          print *,'error in calculate REF'
           stop
       endif
                                                          
	                                                                        
                                                                                
                                                  
        K03=K02*K05

       U1(1)=(U0(1)/DT+K03*V0(1)+GAx*(1-K04))/(1/DT+K03)
       U1(2)=(U0(2)/DT+K03*V0(2)+GAy*(1-K04))/(1/DT+K03)
       U1(3)=(U0(3)/DT+K03*V0(3)+GAz*(1-K04))/(1/DT+K03)

c ------------------------------------------
c Check the convergence of particle velocity in case of KO5>1
c ------------------------------------------

       if (K05.GT.1) then
           do 90 l=1,3 
90         EU(l)=ABS(U1(l)-U0(l))
      
           EU(4)=MAX(EU(1),EU(2),EU(3))
           if (EU(4).gt.EPS) then 
            print *,'fail in velpar'
c            stop
       endif
       endif
      return
      end

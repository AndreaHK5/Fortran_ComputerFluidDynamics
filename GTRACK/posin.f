      subroutine posin(P0,maxR,str,th,nm,DD)
      include 'main.param'
      double precision alpha,th0(maxpart),r0(maxpart),pig
      double precision minR,maxR,str,DD,P0(3,maxpart)
      integer th,nm
   
      minR=maxR-DD/(2.)-str
      pig=4.*atan(1.)
      alpha=th*(2.)*pig/(360.)
          
      do n=1,nm
       r0(n)=RAND(0)*str+minR
       th0(n)=alpha/2+RAND(0)*alpha
       P0(3,n)=0.0
       P0(1,n)=r0(n)*COS(th0(n))*(-1)
       P0(2,n)=r0(n)*SIN(th0(n))
       enddo


       return
       end

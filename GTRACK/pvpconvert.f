      program pvpconvert
      
c ---------------------------
c agguanta con grinta i file pvp_????? non formattati 
c  e li converte in formattati
c ---------------------------
      implicit double precision (a-h,o-z)
      include 'main.param'
      double precision P(3,maxpart),U(3,maxpart)      
      integer first,last,step,nstep,nm,n
      
      
      character *9 pvpfile
      character *5 stepcha
      character *20 outfile

c ------------------------
c Lettura numero di particelle da inputtr
c ------------------------     
       
       
       open (unit=2,file='inputtr',form='formatted')
      

       read (2,*)
       read (2,*)
       read (2,*)
       read (2,*)
       read (2,*)
       read (2,*)
       read (2,*)
       read (2,*)
       read (2,*)
       read (2,*) nm
       read (2,*)
       read (2,*)
       read (2,*)
       read (2,*)
       read (2,*)
       read (2,*)
       read (2,*)
       read (2,*)
       read (2,*)

       close (2)
      write (6,*)'Primo pvp-file da convertire: '
      read (5,*) first
      write (6,*)'Ultimo pvp-file da convertire: '
      read (5,*) last
      write(6,*)'Passo dei file da convertire: '
      read (5,*) step
      
      do 1 nstep=first,last,step
      
       write (stepcha,99) nstep
 99   format (i5.5)    
       pvpfile='pvp_'//stepcha 
       open (unit=10,file=pvpfile,form='unformatted')

c -------------    
c Lettura file pvp unformatted
c -------------

       do n=1,nm
         read (10) (P(l,n),l=1,3),(U(h,n),h=1,3)
       enddo
       close (10)
c --------------
c Scrittura in file formatted
c ---------------
      
      outfile='pvpout.'//stepcha
      open (unit=11,file=outfile,form='formatted')
      print *,'writing on ',outfile
      do n=1,nm
         write(11,*) (P(l,n),l=1,3),(U(h,n),h=1,3)
      enddo
      close (11)
      
      
      
 1    continue     
 
      stop
      end

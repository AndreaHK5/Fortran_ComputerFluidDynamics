             program multiread
c -------------------------------------------------
c estrae campo di moto per plottaggio da file binario
c  e da un pacco di file rst_???? tutti in una volta
c
c Parametri esterni in file formatato inputrd
c
c    i't another Gino's production, 'o duri!
c -------------------------------------------------

      implicit double precision (a-h,o-z)
      include 'main.param'
      character *9 nomefile 
      character *20 fileoutput
      character *20 fileinput
      character *5 timestep
      
      double precision u(il,jl,kl,3),p(il,jl,kl)
      double precision x(il,jl,kl,3)

      integer ilg(ng),jlg(ng),klg(ng),chk
      integer nxx,nyy,nzz
      integer first,last,passo,fc,fisso


c ---------------------
c  fc è file counter
c ---------------------


     
      print *,'   '
      print *,'    --------------------------'
      print *,'    >> a Gino''s production <<'
      print *,'    --------------------------' 
      print *,'     (che credevate, o duri!)'
      print *,''
      print *,'leggo-scrivo-traduco sequenzialmente'
      print *,'coi paramentri in inputrd'
      print *,'(se hai messo i parametri sbagliati sei della Lazio)'
      print *,''
      
      open(15, file='inputrd',form='formatted')
      open(11,file = 'xyu',   form = 'unformatted') 
      
      
      
      read (15,*)
      read (15,*) first,last,passo
      read (15,*)
      read (15,*)
      read (15,*)chk
      read (15,*)
      read (15,*)fisso
      



c ----------------------------
c  read xyu to get the geometry
c ----------------------------

      read(11)lg
      
      if (lg.gt.2) then
         print *,'Cant work with multiple blocks, 
     &  just change file writing in source file(now it overwrites..)!'
      endif
      
      
      do 14 ig = 1,lg
       read(11)ilg(ig), jlg(ig), klg(ig)
       read(11)((((x(i,j,k,l),i = 1,ilg(ig)),
     &          j = 1,jlg(ig)), k = 1,klg(ig)), l = 1,3)

        print *,'block #',ig,' cells=',ilg(ig),jlg(ig),klg(ig)
        print *,''

14    continue
    
    
      do 27 fc=first,last,passo
    
      write (timestep,999) fc
      fileinput='rst_'//timestep

      open(22,file=fileinput, form='unformatted')

      read(22)nt, rtime
       
       do 23 ig=1,lg

c -------------------------------------------
c  start reading the input unformatted file
c -------------------------------------------
 
          read(22) ((((u(i,j,k,l),i=1,ilg(ig)),j=1,jlg(ig)),
     &      k=1,klg(ig)),l=1,3)
         read(22)pp
         read(22)pp
         read(22)pp
         read(22)pp
         read(22) (((p(i,j,k),i=1,ilg(ig)),j=1,jlg(ig)), 
     &      k=1,klg(ig))
         read(22)dummy2o
         read(22)dummy3o
         read(22)dummy2o
         read(22)dummy2o
         
	 
	 if(chk.eq.1) then

	  fileoutput='frontout.'//timestep
	  open(55,file = fileoutput,   form = 'formatted')
           do 11 j=2,jlg(ig)-1
             write(55,*)
	   do 11 k=2,klg(ig)-1
	    write(55,1000)(x(fisso,j,k,l),l=1,3),(u(fisso,j,k,l),l=1,3)
     &                                                    ,p(fisso,j,k) 
 11      continue
         close(55)

	 print *,'Scrivo er file (output)=', fileoutput
 


	endif

c ----------------------
c for side parallel planes
c ----------------------         
	 
	if(chk.eq.2) then
 	  
	  fileoutput='sideout.'//timestep
	  open(55,file =fileoutput,   form = 'formatted')

           do 21 i=2,ilg(ig)-1
              write(55,*)
	   do 21 k=2,klg(ig)-1
            write(55,1000)(x(i,fisso,k,l),l=1,3),(u(i,fisso,k,l),l=1,3)
     &                                                    ,p(i,fisso,k)
 21      continue
         close(55)
	 
	 print *,'Scrivo er file (output)=', fileoutput
 

	 endif
 
c ----------------------
c for side parallel planes
c ----------------------         
	 
	if(chk.eq.3) then
 	  
	   	  
	  fileoutput='upout.'//timestep
	  open(55,file =fileoutput,   form = 'formatted')
	  open(55,file = fileoutput,   form = 'formatted')
	 

           do 31 i=2,ilg(ig)-1
             write(55,*)
	   do 31 j=2,jlg(ig)-1
            write(55,1000)(x(i,j,fisso,l),l=1,3),(u(i,j,fisso,l),l=1,3)
     &                                                    ,p(i,j,fisso)
 31      continue
         close(55)	 
	 
	 print *,'scrivo er file (output)=', fileoutput


	 endif

23     continue

27    continue

 1000 format (6(1x,f12.8),1x,e10.3)
 999  format (i5.5)

 
 
 10   continue 

      stop
      end
  

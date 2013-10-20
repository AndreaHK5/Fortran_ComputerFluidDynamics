             program readoutput
c -------------------------------------------------
c estrae campo di moto per plottaggio da file binario
c    i't another Gino's production, 'o duri!
c -------------------------------------------------

      implicit double precision (a-h,o-z)
      include 'main.param'
      character *9 nomefile 
      character *20 fileoutput
      character *5 extension
      
      double precision u(il,jl,kl,3),p(il,jl,kl)
      double precision x(il,jl,kl,3)

      integer ilg(ng),jlg(ng),klg(ng),chk
      integer nxx,nyy,nzz



     
      print *,'   '
      print *,'    --------------------------'
      print *,'    >> a Gino''s production <<'
      print *,'    --------------------------' 
      print *,'     (che credevate, o duri!)'
      print *,''

      
      
      write (6,*) 'Dime dime, che file te vol convertir,mulo?'
      read (5,*)nomefile
      print *,''
      
      open(22,file=nomefile, form='unformatted')
      open(11,file = 'xyu',   form = 'unformatted')
      
      read(22)nt, rtime
      write (extension,999) nt


c ----------------------------
c  read xyu to get the geometry
c ----------------------------

      read(11)lg
      
      if (lg.gt.2) then
         print *,'Cant work with multiple blocks, 
     &  just change file writing in source file(now it overwrites..)!'
      endif
      
      
      do 10 ig = 1,lg
       read(11)ilg(ig), jlg(ig), klg(ig)
       read(11)((((x(i,j,k,l),i = 1,ilg(ig)),
     &          j = 1,jlg(ig)), k = 1,klg(ig)), l = 1,3)

        print *,'block #',ig,' cells=',ilg(ig),jlg(ig),klg(ig)
        print *,''



c ------------------------------ 
c  start reading the file
c ------------------------------
 
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

c -----------------------
c  usk usr what wana do
c -----------------------


         write(6,*)'What section do u wanna output?'
	 write(6,*)'1=frontal (y-z); 2=lateral (x-z); 3=up (x-y)'
	 read(5,*)chk

c ------------------
c  for frontal parallel plane
c ------------------
         
	 if(chk.eq.1) then
	  
	  fileoutput='frontout.'//extension
	  open(55,file = fileoutput,   form = 'formatted')
	  print *,''
	  write(6,*)'on what vertical plane between 2 and',
     &                                             ilg(ig)-1,'?'
	  read (5,*)nxx
	 

           do 11 j=2,jlg(ig)-1
             write(55,*)
	   do 11 k=2,klg(ig)-1
            write(55,1000)(x(nxx,j,k,l),l=1,3),(u(nxx,j,k,l),l=1,3)
     &                                                    ,p(nxx,j,k)
 11      continue
         close(55)
         print *,''
	 print *,'name file output=', fileoutput
	 print *,'check this out!' 
	 print *,''
	endif

c ----------------------
c for side parallel planes
c ----------------------         
	 
	if(chk.eq.2) then
 	  
	  fileoutput='sideout.'//extension
	  open(55,file =fileoutput,   form = 'formatted')
	  print *,''
	  write(6,*)'on what vertical plane between 2 and',
     &                                             jlg(ig)-1,'?'
	  write(6,*)'(Between',jlg(ig)/2,' and',(jlg(ig)/2)+1,
     &                             ' there will be the symmetry plane)'
	  read (5,*)nyy
	 

           do 21 i=2,ilg(ig)-1
              write(55,*)
	   do 21 k=2,klg(ig)-1
            write(55,1000)(x(i,nyy,k,l),l=1,3),(u(i,nyy,k,l),l=1,3)
     &                                                    ,p(i,nyy,k)
 21      continue
         close(55)
         print *,''
	 print *,'name file output= ', fileoutput
	 print *,'check this out!'
	 print *,''
	 endif
 
c ----------------------
c for side parallel planes
c ----------------------         
	 
	if(chk.eq.3) then
 	  
	   	  
	  fileoutput='upout.'//extension
	  open(55,file =fileoutput,   form = 'formatted')
	  open(55,file = fileoutput,   form = 'formatted')
	  print *,''
	  write(6,*)'on what horizontal plane between 2 and',
     &                                             klg(ig)-1,'?'

	  read (5,*)nzz
	 

           do 31 i=2,ilg(ig)-1
             write(55,*)
	   do 31 j=2,jlg(ig)-1
            write(55,1000)(x(i,j,nzz,l),l=1,3),(u(i,j,nzz,l),l=1,3)
     &                                                    ,p(i,j,nzz)
 31      continue
         close(55)
         print *,''
	 print *,'name file output= ',fileoutput
	 print *,'check this out!'
	 print *,''
	 endif
	 print *,'enjoy'
         print *,''

 1000 format (6(1x,f12.8),1x,e10.3)
 999  format (i5.5)

 
 
 10   continue 

      stop
      end
  

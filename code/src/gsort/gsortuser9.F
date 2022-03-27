! This program is to be used when the Euler angles of DANTE are given with respect to the 
! natural reference frame, i.e. with respect to the beam.


      subroutine usersub9(flag)

      integer flag

#include "gsort.inc"
#include "dante.inc"

      record/dante/dante

      common/dante/dante



      if(flag.eq.0) then      ! from GETINPUT
        call gs0_dante(dante)
      elseif(flag.eq.1) then      ! from LISTPROG
        call gs1_dante(dante)
      elseif(flag.eq.2) then      ! from INIT_RUN
        call gs2_dante(dante)
      elseif(flag.eq.3) then      ! from EVANA
        call gs3_dante(dante)
      elseif(flag.eq.4) then      ! from FINIT_RUN
        call gs4_dante
      endif

      return

      end


**********************************************************
**********************************************************
      subroutine gs0_dante(dante)

#include "gsort.inc"
#include "dante.inc"

      record/dante/dante

      character*72 angfil

      logical*1 getmem
      logical*1 gs_testindpar
      external gs_testindpar

c      synt(1)='ANGLES_DANTE Px Py Angles_file_Dante Pt Pp Pt[beam] Pp[beam] rotation_Dante Grazing'
      synt(1)='USERSUB9 Px Py Angles_file_Dante Pt Pp Pt[beam] Pp[beam] rotation_Dante Grazing'

      if(syntax) then
        call gs_syntax(' ')
        return
      endif

      do nn = 1,2
        call gs_getindpar(dante.ind(nn),dante.par(nn),dante.res(nn))
        if(dante.ind(nn).le.0) call gs_syntax('Meaningless for header parameters')
        if(dante.ind(nn) .ne. dante.ind(1))call gs_syntax('The first and second parameters must be in the same detector')
      enddo

      call istr_ch(comline,angfil)
      if(angfil.EQ.'DANTE') then
        dante.afile_Dante=gsdir(1:lgsdir)//'NUMBERING.DANTE'
        dante.ndet_Dante=10
      else
        dante.afile_Dante=angfil  
        dante.ndet_Dante=10
        if(lengthc(dante.afile_Dante).LT.1) call gs_syntax('Filename expected')
      endif

      do nn = 3,6
        call gs_getindpar(dante.ind(nn),dante.par(nn),dante.res(nn))
        if(dante.ind(nn).le.0) call gs_syntax('Meaningless for header parameters')
        if(dante.ind(nn) .ne. dante.ind(1))call gs_syntax('The third, fourth, fifth and 
     $	sixth parameters must be in the same detector')
      enddo

       
      if( istr_r1(comline, rrotation_Dante) .ne. 1 )call gs_syntax('DANTE rotation angle expected here')
        if( rrotation_Dante .lt. 0.0E0 .or. rrotation_Dante .gt. 360.0E0 ) 
     $	call gs_syntax('Invalid value for DANTE rotation angle [0-360]')
        dante.rrotation_Dante = rrotation_Dante

      if( istr_r1(comline, grazing) .ne. 1 )call gs_syntax('PRISMA grazing angle expected here')
        if( grazing .lt. 0.0E0 .or. grazing .gt. 360.0E0 ) 
     $	call gs_syntax('Invalid value for PRISMA grazing angle [0-360]')
        dante.grazing = grazing
	
      
      isize= 6 +  6*MAXDET + 2*MAXRES
      isize= 4*isize
      isize=isize+16 ! di guardia
      if(.NOT.getmem(isize,com(icmd).addr)) Stop 'Error allocating memory for DANTE'
      com(icmd).size=isize

      return

      end


      subroutine gs1_dante(dante)

#include "gsort.inc"
#include "dante.inc"

      record/dante/dante

      do nn = 1,2
        call gs_putindpar(dante.ind(nn),dante.par(nn))
        llist=llist+3
      enddo
        llist=lengthc(list)+3
        write(list(llist:),'(f5.1)')dante.rrotation_Dante
        write(list(llist:),'(f5.1)')dante.grazing

      return
      end


      subroutine gs2_dante(dante)

#include "gsort.inc"
#include "dante.inc"
      record/dante/dante
      
      real alphaDante(0:MAXDET-1),betaDante(0:MAXDET-1),gammaDante(0:MAXDET-1) 

      integer inlu /-1/

      PARAMETER (PI=3.1415926)       ! some constants
      PARAMETER (Degree180=180.0)
      PARAMETER (DR=PI/Degree180)
      PARAMETER (RD=Degree180/PI)

      if(init_mode.ne.0) return


      call gs_get_eulerangles(alphaDante,betaDante,gammaDante,4,dante.afile_Dante,LU2)
      do ii=0,3
      dante.alphaDante(ii) = alphaDante(ii)
      dante.betaDante(ii) = betaDante(ii)
      dante.gammaDante(ii) = gammaDante(ii)
      write(6,*) dante.alphaDante(ii), dante.betaDante(ii), dante.gammaDante(ii)
      
      dante.cosAlphaDante(ii) = cos(alphaDante(ii)*DR)
      dante.sinAlphaDante(ii) = sin(alphaDante(ii)*DR)
      dante.cosBetaDante(ii)  = cos(betaDante(ii)*DR)
      dante.sinBetaDante(ii)  = sin(betaDante(ii)*DR)
      dante.cosGammaDante(ii) = cos(gammaDante(ii)*DR)
      dante.sinGammaDante(ii) = sin(gammaDante(ii)*DR)

      enddo

      write(LU1,*)

      return

      end


      subroutine gs3_dante(dante)

#include "gsort.inc"
#include "dante.inc"

      record/dante/dante
      real*8 tmp_theta, tmp_phi, xx, yy, zz, xx2, yy2, zz2, theta, phi, modulo, posX, posY, posZ
      real*8 xx3, yy3, zz3, grazing
      
!      PARAMETER (radius=62.5) ! 88-25.5mm (check, radius anello - thickness rivelatore)
!      PARAMETER (radius=64.3) ! Good. Use this one
      PARAMETER (radius=250.0) 
      PARAMETER (PI=3.1415926)       ! some constants
      PARAMETER (Degree180=180.0)
      PARAMETER (DR=PI/Degree180)
      PARAMETER (RD=Degree180/PI)


!      if( ndet( dante.ind(1) ) .ne. 1 )return !multiplicita 1 di DANTE
!      if( ndet( dante.ind(1) ) .ne.  ndet( dante.ind(2) ))return ! Abbiamo la X e la Y nel riv.
      if( ndet( dante.ind(1) ) .le. 0 )return
	do ix = doff(dante.ind(1)),doff(dante.ind(1))+ndet(dante.ind(1))-1
	  if( det(ix).ival(dante.par(1)).lt.0 .or. det(ix).ival(dante.par(2)).lt.0 )then
	      det(ix).ival(dante.par(1)) = 0
	      det(ix).xval(dante.par(1)) = 0.0E0
	      det(ix).ival(dante.par(2)) = 0
	      det(ix).xval(dante.par(2)) = 0.0E0
            else

            xx = 0.0E0
            yy = 0.0E0
            zz = 0.0E0
	    

           ! These matrices are left handed
	   ! Utilizzando gli angoli di Euler rotiamo i rivelatori alla sua posizione
           ! devo togliere prima 1000 per la questione del offset messo nel setup file 
	   ! dopo devo dividere per 10. 
           det(ix).xval(dante.par(1)) = det(ix).xval(dante.par(1)) - 1000
           det(ix).xval(dante.par(2)) = det(ix).xval(dante.par(2)) - 1000

           det(ix).xval(dante.par(1)) = det(ix).xval(dante.par(1))*0.1
           det(ix).xval(dante.par(2)) = det(ix).xval(dante.par(2))*0.1

           posX = det(ix).xval(dante.par(1))
           posY = det(ix).xval(dante.par(2))
           posZ = radius

ccccccccccccc rotazione Rz(psi)

           dante.posX = dante.cosAlphaDante(det(ix).id) * posX - dante.sinAlphaDante(det(ix).id) * posY
           dante.posY = dante.sinAlphaDante(det(ix).id) * posX + dante.cosAlphaDante(det(ix).id) * posY
           dante.posZ = posZ

           posX = dante.posX
           posY = dante.posY
           posZ = dante.posZ
	   

ccccccccccccc rotazione Ry(theta)

           dante.posX = dante.cosBetaDante(det(ix).id) * posX + dante.sinBetaDante(det(ix).id) * posZ
           dante.posY = posY
           dante.posZ = -dante.sinBetaDante(det(ix).id) * posX + dante.cosBetaDante(det(ix).id) * posZ

           posX = dante.posX
           posY = dante.posY
           posZ = dante.posZ
	   

ccccccccccccc rotazione Rz(psi)

           dante.posX = dante.cosGammaDante(det(ix).id) * posX - dante.sinGammaDante(det(ix).id) * posY
           dante.posY = dante.sinGammaDante(det(ix).id) * posX + dante.cosGammaDante(det(ix).id) * posY
           dante.posZ = posZ

           posX = dante.posX
           posY = dante.posY
           posZ = dante.posZ
	   
       
           det(ix).xval(dante.par(3)) = !theta
     $      acos(posZ/(dsqrt(posX*posX+posY*posY+posZ*posZ)))
     
           det(ix).xval(dante.par(4)) = !phi
     $      atan(posY/posX)
     
   

!      write(6,*)'posX',det(ix).xval(dante.par(1)),'posY',det(ix).xval(dante.par(2))
!      write(6,*) 'theta',det(ix).xval(dante.par(3))*RD,'phi',det(ix).xval(dante.par(4))*RD
!      write(6,*) 'theta2',theta*RD,'phi2',phi*RD

       
!      write(6,*) 'riv dante',dante.alphaDante(det(ix).id),dante.betaDante(det(ix).id),dante.gammaDante(det(ix).id),det(ix).id,ix 

       det(ix).xval(dante.par(3)) = det(ix).xval(dante.par(3))*RD

       ! mettiamo la phi al suo posto
       if((posX .ge. 0.00D0) .and. (posY .ge. 0.00D0)) then
         det(ix).xval(dante.par(4)) = det(ix).xval(dante.par(4))*RD
         endif
       if((posX .ge. 0.00D0) .and. (posY .lt. 0.00D0)) then
         det(ix).xval(dante.par(4)) = 360 + det(ix).xval(dante.par(4))*RD
	 endif
       if((posX .lt. 0.00D0) .and. (posY .ge. 0.00D0)) then
         det(ix).xval(dante.par(4)) = 180 + det(ix).xval(dante.par(4))*RD
	 endif
       if((posX .lt. 0.00D0) .and. (posY .lt. 0.00D0)) then
         det(ix).xval(dante.par(4)) = 180 + det(ix).xval(dante.par(4))*RD
	 endif



!We rotate DANTE with respect to the beam reference frame

! These matrices are right handed to be able to put the grazing in it easily, 
!otherwise I had to write 360-grazing angle, so it is less understandable.

       tmp_theta = det(ix).xval(dante.par(3))*DR
       tmp_phi = det(ix).xval(dante.par(4))*DR
       
       xx3 = sin(tmp_theta)*cos(tmp_phi)
       
       yy3 = cos(dante.rrotation_Dante*DR)*sin(tmp_theta)*sin(tmp_phi) + 
     $      sin(dante.rrotation_Dante*DR)*cos(tmp_theta)
       
       zz3 = -sin(dante.rrotation_Dante*DR)*sin(tmp_theta)*sin(tmp_phi) + 
     $      cos(dante.rrotation_Dante*DR)*cos(tmp_theta)


       det(ix).xval(dante.par(5)) = acos(zz3)
       det(ix).xval(dante.par(6)) = atan(yy3/xx3)

!       write(6,*)'y,x',yy3,xx3
!       write(6,*)'atan',atan(yy3,xx3)*RD
       
       det(ix).xval(dante.par(5)) = det(ix).xval(dante.par(5))*RD
       
       ! mettiamo la phi al suo posto
       if((xx3 .ge. 0.00D0) .and. (yy3 .ge. 0.00D0)) then
	 det(ix).xval(dante.par(6)) = det(ix).xval(dante.par(6))*RD
	 endif
       if((xx3 .ge. 0.00D0) .and. (yy3 .lt. 0.00D0)) then
	 det(ix).xval(dante.par(6)) = 360 + det(ix).xval(dante.par(6))*RD
	 endif
       if((xx3 .lt. 0.00D0) .and. (yy3 .ge. 0.00D0)) then
	 det(ix).xval(dante.par(6)) = 180 + det(ix).xval(dante.par(6))*RD
	 endif
       if((xx3 .lt. 0.00D0) .and. (yy3 .lt. 0.00D0)) then
	 det(ix).xval(dante.par(6)) = 180 + det(ix).xval(dante.par(6))*RD
	 endif


      ! in gradi e *10, perche cosi lo capisce il resto del programa
      det(ix).xval(dante.par(5)) = det(ix).xval(dante.par(5))*10 
      det(ix).ival(dante.par(5)) = det(ix).xval(dante.par(5))
      det(ix).xval(dante.par(6)) = det(ix).xval(dante.par(6))*10 
      det(ix).ival(dante.par(6)) = det(ix).xval(dante.par(6))


! We rotate DANTE to be placed at the position reference frame of CLARA-PRISMA 
       
       grazing = dante.grazing - dante.rrotation_Dante
       
       tmp_theta = det(ix).xval(dante.par(3))*DR
       tmp_phi = det(ix).xval(dante.par(4))*DR
       
       xx2 = sin(tmp_theta)*cos(tmp_phi)
       
       yy2 = cos(grazing*DR)*sin(tmp_theta)*sin(tmp_phi) - 
     $      sin(grazing*DR)*cos(tmp_theta)
       
       zz2 = sin(grazing*DR)*sin(tmp_theta)*sin(tmp_phi) + 
     $      cos(grazing*DR)*cos(tmp_theta)


       det(ix).xval(dante.par(3)) = acos(zz2)
       det(ix).xval(dante.par(4)) = atan(yy2/xx2)


       det(ix).xval(dante.par(3)) = det(ix).xval(dante.par(3))*RD
       
       ! mettiamo la phi al suo posto
       if((xx2 .ge. 0.00D0) .and. (yy2 .ge. 0.00D0)) then
	 det(ix).xval(dante.par(4)) = det(ix).xval(dante.par(4))*RD
	 endif
       if((xx2 .ge. 0.00D0) .and. (yy2 .lt. 0.00D0)) then
	 det(ix).xval(dante.par(4)) = 360 + det(ix).xval(dante.par(4))*RD
	 endif
       if((xx2 .lt. 0.00D0) .and. (yy2 .ge. 0.00D0)) then
	 det(ix).xval(dante.par(4)) = 180 + det(ix).xval(dante.par(4))*RD
	 endif
       if((xx2 .lt. 0.00D0) .and. (yy2 .lt. 0.00D0)) then
	 det(ix).xval(dante.par(4)) = 180 + det(ix).xval(dante.par(4))*RD
	 endif


      ! in gradi e *10 e perche cosi lo capisce il resto del programa
      det(ix).xval(dante.par(3)) = det(ix).xval(dante.par(3))*10 
      det(ix).ival(dante.par(3)) = det(ix).xval(dante.par(3))
      det(ix).xval(dante.par(4)) = det(ix).xval(dante.par(4))*10 
      det(ix).ival(dante.par(4)) = det(ix).xval(dante.par(4))

!      write(6,*) 'thetadopo',det(ix).xval(dante.par(3)),'phidopo',det(ix).xval(dante.par(4))
!      write(6,*) 'theta2dopo',theta*RD,'phi2dopo',phi*RD, modulo
!      write(6,*) 'thetatmp',tmp_theta*RD,'phitmp',tmp_phi*RD,dante.rrotation_Dante     
!      write(6,*) 'thetadopo',det(ix).xval(dante.par(3)),'phidopo',det(ix).xval(dante.par(4))
	    
	    
   	   endif
	enddo

      
      return

      end
      
      

      subroutine gs4_dante

      return

      end


	subroutine gs_get_eulerangles(al,be,ga,nn,angfile,lun)

	integer nn,lun
	real al(0:nn-1) 
	real be(0:nn-1)
        real ga(0:nn-1)
	character angfile*(*) 

	character*80 line
	integer inlu /-1/

	if(inlu.le.0) call lib$get_lun(INLU)

	langfile=max(1,lengthc(angfile))

6	open(unit=INLU,file=angfile,status='old',READONLY,ERR=7)
	nga=0
	write(6,*)
	write(6,*) 'Reading angles of the detectors from  '//angfile(1:langfile)
	if(lun.ne.6) write(lun,*) 'Reading angles of the detectors from  '//angfile(1:langfile)
	goto 8

7	write(6,*) 'Error reading  '//angfile(1:langfile)
	if(inp_yes('Retry')) goto 6
	call exit(0)

8	read(INLU,'(A)',end=10) line 
	lline=lengthc(line)
	if(lline.gt.0) then
	  iriv=index(line,'Riv#')
	  if(iriv.gt.0) then
	    read(line,'(6x,i4)')ii
	    if(ii.ge.0 .AND. ii.lt.nn) then
	      read(line,'(6x, i4, 8x, f10.4, 10x, f10.4, 10x, f10.4)') ii,al(ii),be(ii),ga(ii)
	      nga=nga+1
	      write(lun,'(i6,''   Alpha='',f10.4,''     Beta='',f10.4,''     Gamma='',f10.4)') ii,al(ii),be(ii),ga(ii)
	    endif
	  endif
	endif
	goto 8

10	close(unit=inlu)
	call lib$free_lun(INLU)
	inlu=-1
	
	if(nga.lt.nn) then
	               write(  6,*) 'W A R N I N G  only',nga,' detector''s angles found' 
	  if(lun.ne.6) write(lun,*) 'W A R N I N G  only',nga,' detector''s angles found' 
	endif

	return

	end

      subroutine usersub8(flag)

      integer flag

#include "gsort.inc"
#include "mcppar.inc"

      record/mcppar/mcp
      common/mcppar/mcp


      if(flag.eq.0) then      ! from GETINPUT
        call gs0_mcppar(mcp)
      elseif(flag.eq.1) then      ! from LISTPROG
        call gs1_mcppar(mcp)
      elseif(flag.eq.2) then      ! from INIT_RUN
        call gs2_mcppar(mcp)
      elseif(flag.eq.3) then      ! from EVANA
        call gs3_mcppar(mcp)
      elseif(flag.eq.4) then      ! from FINIT_RUN
        call gs4_mcppar
      endif

      return

      end




**********************************************************
**********************************************************
      subroutine gs0_mcppar(mcp)

#include "gsort.inc"
#include "mcppar.inc"
      record/mcppar/mcp

      logical*1 getmem
      logical*1 gs_testindpar
      external gs_testindpar
      synt(1)='USERSUB8 At Ap Ax Ay'
c      synt(1)='MCP_PAR At Ap Ax Ay'

      if(syntax) then
        call gs_syntax(' ')
        return
      endif

      do nn = 1,4
        call gs_getindpar(mcp.ind(nn),mcp.par(nn),mcp.res(nn))
        if(mcp.ind(nn).le.0) call gs_syntax('Meaningless for header parameters')
        if(mcp.ind(nn) .ne. mcp.ind(1))call gs_syntax('All parameters must be in the same detector')
      enddo

      isize= 6 +  6*256 + 2*MAXRES
      isize= 4*isize
      isize=isize+16 ! di guardia
      if(.NOT.getmem(isize,com(icmd).addr)) Stop 'Error allocating memory for MCP_PAR'
      com(icmd).size=isize

      return

      end


**********************************************************
**********************************************************
      subroutine gs1_mcppar(mcp)

#include "gsort.inc"
#include "mcppar.inc"
      record/mcppar/mcp

      do nn = 1,4
        call gs_putindpar(mcp.ind(nn),mcp.par(nn))
        llist=llist+3
      enddo

      return
      end



**********************************************************
**********************************************************
	subroutine gs2_mcppar(mcp)

#include "gsort.inc"
#include "mcppar.inc"
	record/mcppar/mcp
	
	return

	end



**********************************************************
**********************************************************
      subroutine gs3_mcppar(mcp)

#include "gsort.inc"
#include "mcppar.inc"
      record/mcppar/mcp

	real*8 theta, phi, xm, ym, y, dist, angle
	
        
        ! checks multiplicity for MCP and PPAC (should be 1)
	if( ndet( mcp.ind(1) ) .ne. 1 )then
	   killev = .true.
	   return
	endif
	
	theta = 0.1 * det(doff(mcp.ind(1))).xval(mcp.par(1))
	phi   = 0.1 * det(doff(mcp.ind(2))).xval(mcp.par(2))
	dist  = 250.0
	angle = 135.0
	
c	write(*,*) 'Theta, phi ', theta, phi
	
        ! trivial when trajectory is along the vertical plane x=0
        if( ((theta*phi) .eq. 0.) .or. (phi .eq. 180.) ) then
          xm = 0.
          ym = dist * cosd(phi) * tand(theta)
        ! otherwise find the intersection of two lines
        else
          y  = tand(theta)*sind(phi)*dist / ( 1. + tand(angle) * tand(theta) * sind(phi) )
          xm = y / cosd(angle)
          ym = -sind(angle) * xm / tand(phi);
        endif

c	write(*,*) 'X, Y ', xm, ym

	
	det(doff(mcp.ind(3))).xval(mcp.par(3)) = 1000. + 20. * xm
c	det(doff(mcp.ind(3))).xval(mcp.par(3)) = xm
	det(doff(mcp.ind(3))).ival(mcp.par(3)) = det(doff(mcp.ind(3))).xval(mcp.par(3))
	
	det(doff(mcp.ind(4))).xval(mcp.par(4)) = 1000. + 20. * ym
c	det(doff(mcp.ind(4))).xval(mcp.par(4)) = ym
	det(doff(mcp.ind(4))).ival(mcp.par(4)) = det(doff(mcp.ind(4))).xval(mcp.par(4))
	
      return

      end


**********************************************************
**********************************************************
      subroutine gs4_mcppar

      return

      end


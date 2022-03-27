
        structure /xfp_ic/
c PPAC data (x)
	    integer ind_xf
	    integer par_xf
c ionization chamber
	    integer ind_ic
	    integer par_ic(4)
	    real    thr
c parameters of the PRISMA geometry
	    real*8  x_ic(4)     ! "x" coordinate to make linear fit
	    real*8  y_ic(4)     ! "y" coordinate to make linear fit
	    real*8  e_ic(4)     ! weights
	    real*8  ifp_dist	! IC-PPAC distance
	    real*8  ic_depth	! IC section depth
	    real*8  ic_width	! IC section width
	end structure


	subroutine gs0_xfp_ic( pms )
	

#include "gsort.inc"
	record /xfp_ic/pms
	
	integer ind_t, ind_f, ind_l, ind_r, ind_e, ind_p, ind_ra, ii, ibsize
	real*8  dfp_dist,xf,yf,fp_half_length
	logical*1 getmem

	synt(1)='XFP_IC Px (output: PPAC X [mm]) Ia Ib Ic Id [low_thr] (ind. energies in the sections of the IC and lower threshold)'
	

	if(syntax) then
	  call gs_syntax(' ')
	  return
	endif

	call gs_getindpar(pms.ind_xf, pms.par_xf, ii)
	if( pms.ind_xf .lt. 1)call gs_syntax('Meaningless for header parameters')
	
	call gs_getindpar(pms.ind_ic, pms.par_ic(1), ii)
	if( pms.ind_ic .lt. 1)call gs_syntax('Meaningless for header parameters')
	
	call gs_getindpar(ind_r, pms.par_ic(2), ii)
	if(ind_r .ne. pms.ind_ic)call gs_syntax('Must be the same detector')
	if(pms.par_ic(1) .eq. pms.par_ic(2))call gs_syntax('All IC parameters must be different')
	
	call gs_getindpar(ind_r, pms.par_ic(3), ii)
	if(ind_r .ne. pms.ind_ic)call gs_syntax('Must be the same detector')
	if(pms.par_ic(1) .eq. pms.par_ic(3))call gs_syntax('All IC parameters must be different')
	if(pms.par_ic(2) .eq. pms.par_ic(3))call gs_syntax('All IC parameters must be different')

	call gs_getindpar(ind_r, pms.par_ic(4), ii)
	if(ind_r .ne. pms.ind_ic)call gs_syntax('Must be the same detector')
	if(pms.par_ic(1) .eq. pms.par_ic(4))call gs_syntax('All IC parameters must be different')
	if(pms.par_ic(2) .eq. pms.par_ic(4))call gs_syntax('All IC parameters must be different')
	if(pms.par_ic(3) .eq. pms.par_ic(4))call gs_syntax('All IC parameters must be different')
	
	if( istr_r1(comline,pms.thr) .ne. 1 ) pms.thr = 0.0
	
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccc All geometry parameters are given values here 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	pms.ifp_dist = 720.0	! IC-PPAC distance
	pms.ic_depth = 250.0  	! IC section depth
	pms.ic_width = 100.0  	! IC section width
	
	do ii=1,4
	  pms.x_ic(ii) = pms.ifp_dist + ( ii - 0.5 ) * pms.ic_depth
	end do
	
        ibsize  = 4096*8
	if(.NOT.getmem( ibsize ,com(icmd).addr ) ) Stop 'Error allocating memory for XFP_IC'
            
	com(icmd).size=ibsize
	return
	end
	

	subroutine gs1_xfp_ic( pms )
	
#include "gsort.inc"
	record /xfp_ic/pms

	call gs_putindpar(pms.ind_xf, pms.par_xf)
	llist = llist+3

	do ii=1,4
	  call gs_putindpar(pms.ind_ic, pms.par_ic(ii))
	  llist = llist+3
	end do
	write(list(llist+1:),'(F7.2)') pms.thr
	call gs_writeline(LLUN,list,INDENT)
		
	return
	end

c  gs_2_prisma_solver: all geometry parameters have been initialized
c                      in gs0_prisma_solver
	subroutine gs2_xfp_ic
	
	return
	end

	subroutine gs3_xfp_ic( pms )
	
#include "gsort.inc"
	record /xfp_ic/pms

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccc Several variables
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc	
	real*8  sumx, sumx2, sumy, sumxy, intercept
	
	do ii=1,4
	  pms.y_ic(ii) = 0.0
	  pms.e_ic(ii) = 0.0
	end do
	
	
        ! Constructs points for linear regression (weighted average of sections with energies)
	do ii=doff(pms.ind_ic), doff(pms.ind_ic)+ndet(pms.ind_ic)-1
	  do jj=1,4
	    if(det(ii).xval(pms.par_ic(jj)) .gt. pms.thr) then
	      pms.e_ic(jj) = pms.e_ic(jj) + det(ii).xval(pms.par_ic(jj))
	      pms.y_ic(jj) = pms.y_ic(jj) + ( det(ii).id + 0.5 ) * pms.ic_width * det(ii).xval(pms.par_ic(jj))
	    endif
	  end do
	end do
	
	do jj=1,4
	  if( pms.e_ic(jj) .gt. 0.0 ) pms.y_ic(jj) = pms.y_ic(jj) / pms.e_ic(jj)
	end do
	
	sumx  = 0.0
	sumx2 = 0.0
	sumy  = 0.0
	sumxy = 0.0
	
	do jj=1,4
	  sumx   = sumx  + pms.x_ic(jj)
	  sumx2  = sumx2 + pms.x_ic(jj)*pms.x_ic(jj)
	  sumy   = sumy  + pms.y_ic(jj)
	  sumxy  = sumxy + pms.y_ic(jj)*pms.x_ic(jj)
	end do
	
	intercept = ( sumx * sumxy - sumx2 * sumy ) / ( sumx * sumx - 4 * sumx2 )
	
	if( intercept .lt. 0.0 )    goto 9000
	if( intercept .gt. 1000.0 ) goto 9000
	
	
	! Finally store results
	det(pms.ind_xf).xval(pms.par_xf) = intercept
	det(pms.ind_xf).ival(pms.par_xf) = det(pms.ind_xf).xval(pms.par_xf)

	
	return

9000    killev = .true.
	return
	end
	

c  gs_4_prisma_solver: nothing to do
	subroutine gs4_xfp_ic
	
	return
	end

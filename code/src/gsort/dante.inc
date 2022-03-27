      structure/dante/
        INTEGER IND(6)		! Tipo di parametro
        INTEGER PAR(6)		! Parametro
        INTEGER RES(6)		! Risoluzione del parametro
        INTEGER NDET_DANTE	! quanti dante
        CHARACTER*72 AFILE_Dante	! angoli dante
	REAL    alphaDante(10)
	REAL    betaDante(10)
	REAL    gammaDante(10)
	REAL    cosAlphaDante(10)
	REAL    cosBetaDante(10)
	REAL    cosGammaDante(10)
	REAL    sinAlphaDante(10)
	REAL    sinBetaDante(10)
	REAL    sinGammaDante(10)
        REAL*4  grazing !Angolo di grazing di Prisma
        REAL*4  rrotation_Dante ! Angolo di rotazione della camera
	REAL    posX
	REAL    posY
	REAL    posZ
      endstructure

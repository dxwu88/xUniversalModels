#include "stdafx.h"
#include "EOSFlash.h"


EOSFlash::EOSFlash()
{
}


EOSFlash::~EOSFlash()
{
}

/***
C
SUBROUTINE EOSFlash(FEED, LIQ, VAP, PSTEMP, PSPRESS, VAPFRAC,
	&Enthalpy, PFLAG, IERROR, EOS_Type)
	C *** CALCULATE VLE FOR Hydroprocessing USING SRK ***
	C
	C  VARIABLES :
C	FEED(KL)    REAL * 8		MOLES OF MATERIAL TO BE FLASHED
C	LIQ(KL)	    REAL * 8		MOLES OF LIQUID
C	VAP(KL)	    REAL * 8		MOLES OF VAPOR
C	PSTEMP		REAL * 8		FLASH TEMPERATURE, DEGREES F
C	PSPRESS		REAL * 8		FLASH PRESSURE, PSIG
C	VAPFRAC		REAL * 8		MOLAR FRACTION  OF VAPOR TO FEED
C	ENTHALPY	REAL * 8		ENTHALPY OF FEED, BTU / HR
C	PFLAG		INT * 4		PRINT FLAG : 0 = NO OUTPUT FILE
C     IERROR		INT * 4		ERROR FLAG : 0 = NO ERRORS
C
c
c  Expose subroutine
C**********************************************************************
!DEC$ATTRIBUTES DLLEXPORT, alias : 'EOSFlash'    ::EOSFlash
!DEC$ATTRIBUTES STDCALL, REFERENCE::EOSFlash
!DEC$ATTRIBUTES MIXED_STR_LEN_ARG::EOSFlash
C**********************************************************************
C
USE DFPORT
use COMPDATA
IMPLICIT NONE
EXTERNAL ZSRK
integer * 4 IERROR, PFLAG, IER, I, EOS_Type
REAL(KIND = 8) FEED(KL), LIQ(KL), VAP(KL), PSTEMP, PSPRESS,
&VAPFRAC, TOTALV, TOTAL, DHLMX, DHVMX, Z, PATM, TK, ZSRK
REAL(KIND = 8) VAPK(KL), TEMP, PRESS,
&X(KL), Y(KL), TOTALL, Enthalpy
character * 24 systime

C   CONVERT T AND P UNITS TO R AND PSIA
TEMP = PSTEMP + 459.67
TK = TEMP / 1.8
PRESS = PSPRESS + 14.696
PATM = PRESS / 14.696
IER = 0
DHLMX = 0.0
DHVMX = 0.0
LIQ = 0.0
VAP = 0.0
C   CALL THE SRK FLASH CALCULATIONS
CALL UFLSHSRK(FEED, VAP, LIQ, TEMP, PRESS, DHLMX, DHVMX, IER, EOS_Type)
IERROR = IER
C   CALCULATE MOLAR VAPOR FRACTION
TOTALL = sum(LIQ)
TOTALV = sum(VAP)
TOTAL = TOTALL + TOTALV
IF(TOTAL.GT.0.0) THEN
VAPFRAC = TOTALV / TOTAL
ELSE
VAPFRAC = 0.0
ENDIF
VAPK = 0.0
DO I = 1, KL
IF(FEED(I).GT.0.0) THEN
Y(I) = VAP(I) / TOTALV
X(I) = LIQ(I) / TOTALL
VAPK(I) = Y(I) / MAX(X(I), 1.0E-20)
ENDIF
ENDDO
C	CALCULATE STREAM ENTHALPY
CALL IdealGasEnthalpy(FEED, TK, Enthalpy)
Enthalpy = Enthalpy - (DHLMX*TOTALL + DHVMX * TOTALV)

C   PRINT FILE IS DESIRED
IF(PFLAG.NE.0) THEN
OPEN(unit = 66, FILE = 'EOSFlash.TXT', STATUS = 'UNKNOWN')
systime = CTIME(TIME())
WRITE(66, *) ' *** SUBROUTINE EOSFlash ***'
write(66, *) ' Current date and time is ', systime

WRITE(66, '(A,F7.1,A,F7.1,A)') ' *** FLASH RESULTS: AT T = ',
&PSTEMP, ' AND P = ', PSPRESS, ' ***'
write(66, *) ' Enthalpy correction DHLMX = ', DHLMX
write(66, *) ' Enthalpy correction DHVMX = ', DHVMX
write(66, *) ' Enthalpy correction DHLMX * Sum(LIQ) = ',
*DHLMX*sum(liq)
write(66, *) ' Enthalpy correction DHVMX * Sum(VAP) = ',
*DHVMX*sum(VAP)
write(66, *) ' Stream Enthalpy = ', Enthalpy
write(66, *) ' Error Flag (K) = ', ier
write(66, *) ' Vapor fraction = ', VAPFRAC
WRITE(66, '(3A)') ' C#        Feed        Liquid      ',
&'X           Vapor',
&'        Y           MW          K-VAP           Z'

DO I = 1, KL
Z = ZSRK(I, TK, PATM)
write(66, '(i3,6x,8g12.4)') i, feed(i), liq(i), X(I),
&vap(i), Y(I), MW(I), VAPK(I), Z
ENDDO

CLOSE(66)
ENDIF

RETURN
END



SUBROUTINE UFLSHSRK(STREAM, VAPOR, LIQUID, TEMP, PRES, DHLMX, DHVMX,
	*IER, EOS_Type)
	C
	C================================================================================
	C  V / L equilibrium flash calculation.  (SRK THERMO)
	C  Copyright(c) 1986 by R.A.Russell, Deerhaven Technical Software Co.
	C
	C     Mike Caracotsios added support for Enthalpy calculations
	C--------------------------------------------------------
	C     STREAM(KL) - Real * 8 - Component mole array, LB - MOLE
	C     VAPOR(KL) - Real * 8 - Vapor mole array, LB - MOLE
	C     LIQUID(KL) - Real * 8 - Liquid mole array, LB - MOLE
	C     TEMP - Real * 8 - Temperature, Deg R
	C     PRES - Real * 8 - Pressure, PSIA
	C     DHLMX - Real * 8 - LIQUID RESIDUAL ENTHALPY, BTU / LB - MOLE
	C     DHVMX - Real * 8 - VAPOR RESIDUAL ENTHALPY, BTU / LB - MOLE
	C     IER - Int * 4 - Error flag, 0 = no error
	C================================================================================
	C
	USE COMPDATA

	IMPLICIT NONE

	REAL * 8 STREAM(KL), VAPOR(KL), ALPHA(KL), ERR(3), OLDALF(KL), VAR(3)
	REAL * 8 LIQUID(KL), LF, KB, LNKB, KB0, LIQ, FEDLF, DELTA, RMIN
	REAL * 8 RMAX, SUMP, SUMPA, V, X, ATOLER, ERROR, DENOM, VRATIO, B, C
	REAL * 8 SUMFED, FEEDWT, R, TEMP, PRES, DHLMX, DHVMX, AKB, HCLIQ
	INTEGER * 4 RCODE, IER, I, MCOUNT, MVCOUNT, KK, ITERA, ITERR
	INTEGER * 4 MINFLG, MAXFLG, NOTCON, EOS_Type

	REAL * 8 TotalFeed, ZFeed(KL), Bubble, Dew, SUM_Bubble, SUM_Dew

	IER = 0
	SUMFED = 0.0D0
	FEEDWT = 0.0D0
	LIQUID = 0.0D0
	VAPOR = 0.0D0
	VAR = 0.0D0
	ERR = 0.0D0
	DO I = 1, KL
	SUMFED = SUMFED + STREAM(I)
	FEEDWT = FEEDWT + STREAM(I) * MW(I)
	ENDDO

	IF(SUMFED.LE. 0.0D0)THEN
	IER = -101
	RETURN
	ENDIF
	C
	!C component 1 is Hydrogen
	!FEDLF = SUM(STREAM(1:KL))
	!IF(FEDLF.LE.0.0) THEN
	!LF = (SUMFED - STREAM(i_H2)) / SUMFED     !GUESS A LIQUID FRACTION
	!ELSE
	!LF = FEDLF / SUMFED     !GUESS A LIQUID FRACTION
	!ENDIF
	!
	!LF = (SUMFED - STREAM(i_H2)) / SUMFED     !GUESS A LIQUID FRACTION

	LF = 0.651037591       !Intial liquid fraction

	KB = 1.0D0
	LNKB = 0.0D0
	R = 1.0D0 - LF
	MCOUNT = 0
	MVCOUNT = 0
	C  Get Kb model coefficients and volatilities.For flash to temp., no
	C  model is used, alpha = K, and 'R' will be found so calc.Kb = 1.

	KK = 0
	CALL UKBMODL(KK, TEMP, PRES, AKB, KB, STREAM, ALPHA, DHLMX, DHVMX, IER,
		&EOS_Type)

	IF(IER.NE. 0)RETURN
	KB0 = KB
	B = 0.0D0
	C = 0.0D0

	C  General convergence loop limited to 500 iterations
	DO ITERA = 1, 500
	ITERR = 0
	DELTA = 0.0D0
	RMIN = 0.0D0
	RMAX = 1.0D0
	MINFLG = 0
	MAXFLG = 0

	C  Evaluate current value of 'R', the inner loop variable.Get the 'p'
	C  factors(use the liquid vector).
	RCODE = 0
	DO WHILE(RCODE.EQ.0)
	SUMP = 0.0D0
	SUMPA = 0.0D0
	DO I = 1, KL
	IF(((1.0D0 - R) + R * KB0*ALPHA(I)).LE. 0.0D0)THEN
	IER = -102
	RETURN
	ENDIF
	LIQUID(I) = STREAM(I) / ((1.0D0 - R) + R * KB0*ALPHA(I))
	SUMP = SUMP + LIQUID(I)
	SUMPA = SUMPA + LIQUID(I) * ALPHA(I)
	ENDDO

	C  Calculate total liquid flow and ln(kb) for the current 'R'.

	IF(R.GT.0.0D0)THEN
	LIQ = (1.0D0 - R) * SUMP
	ELSE
	LIQ = SUMFED
	END IF
	V = SUMFED - LIQ
	IF(SUMPA.LE. 0.0D0.OR.SUMP.LE. 0.0D0)THEN
	IER = -103
	RETURN
	ENDIF
	LNKB = DLOG(SUMP / SUMPA)
	ERR(1) = LNKB

	C  If 'R' is at upper limit and error is still negative, no hydrocarbon
	C  liquid phase.If 'R' is at zero and error is still positive, feed is
	C  subcooled.

	IF((R.EQ.1.0D0.AND.ERR(1).LT.0.0D0).OR.
		+ (R.EQ.0.0D0.AND.ERR(1).GT.0.0D0)) EXIT

	C  Check inner loop convergence.

	ITERR = ITERR + 1
	IF(DABS(ERR(1)).LE. 1.0D - 05)EXIT
	VAR(1) = R
	CALL UCURVEF(VAR, ERR, ITERR, B, C, RCODE, DELTA, RMIN, RMAX,
		*MINFLG, MAXFLG, IER)
	IF(IER.NE. 0)RETURN
	R = VAR(1)
	DELTA = 0.3D0
	END DO

	C  Get new volatilities, check.
	20		NOTCON = 0
	OLDALF = ALPHA
	KK = 1
	CALL UKBMODL(KK, TEMP, PRES, AKB, KB, LIQUID, ALPHA, DHLMX, DHVMX, IER,
		&EOS_Type)

	IF(IER.NE. 0)RETURN

	IF(R.EQ. 0.0D0)THEN           !WHEN "R" IS ZERO, STREAM IS ALL LIQUID
	MCOUNT = MCOUNT + 1         !THIS SOMETIMES CREATES A PROBLEM IN
	ELSE                           !ALPHA CONVERGENCE.COUNT HOW MANY
	MCOUNT = 0                  !CONSECUTIVE TIMES "R" IS ZERO
	ENDIF
	IF(R.EQ.RMAX)THEN           !WHEN "R" IS RMAX, STREAM IS ALL VAPOR
	MVCOUNT = MVCOUNT + 1       !THIS SOMETIMES CREATES A PROBLEM IN
	ELSE                           !ALPHA CONVERGENCE.COUNT HOW MANY
	MVCOUNT = 0                 !CONSECUTIVE TIMES "R" IS RMAX
	ENDIF
	DO I = 1, KL
	X = STREAM(I) / SUMFED
	IF(X.GE. 1.0D - 10)THEN
	C  Make tolerance looser for trace components.
	IF(X.LT. 1.0D - 05)THEN
	ATOLER = 1.0D - 03
	c					 ATOLER = 1.0D - 06
	ELSE IF(X.LT. 1.0D - 03)THEN
	ATOLER = 5.0D - 04
	c					 ATOLER = 5.0D - 07
	ELSE
	ATOLER = 1.0D - 06
	c					 ATOLER = 1.0D - 09
	END IF
	IF(OLDALF(I).EQ. 0.0D0)THEN
	IER = -104
	RETURN
	ENDIF
	ERROR = DABS(ALPHA(I) / OLDALF(I) - 1.0D0)
	IF(ERROR.GT.ATOLER)NOTCON = 1
	ENDIF
	END DO

	C  Get new Kb model coefficients if needed.
	IF(NOTCON.EQ. 0)THEN        !DONE LOAD UP VAPOR AND LIQUID MOLE ARRAYS
	HCLIQ = (1.0D0 - R)*SUMP
	DENOM = (1.0D0 - R)*SUMPA
	IF(DENOM.GT.0.0D0)THEN
	VRATIO = V / DENOM
	ELSE
	VRATIO = 1.0D0
	END IF
	DO I = 1, KL
	IF(HCLIQ.EQ.0.0D0.OR.ALPHA(I)*VRATIO.GT.1.0D0)THEN
	LIQUID(I) = LIQUID(I)*(1.0D0 - R)
	IF(LIQUID(I).LT. 1.0D - 10)LIQUID(I) = 0.0D0
	VAPOR(I) = STREAM(I) - LIQUID(I)
	IF(VAPOR(I).LT.1.0D - 10)THEN
	VAPOR(I) = 0.0D0
	LIQUID(I) = STREAM(I)
	ENDIF
	ELSE
	VAPOR(I) = ALPHA(I)*LIQUID(I)*V / SUMPA
	IF(VAPOR(I).LT.1.0D - 10)VAPOR(I) = 0.0D0
	LIQUID(I) = STREAM(I) - VAPOR(I)
	IF(LIQUID(I).LT.1.0D - 10)THEN
	LIQUID(I) = 0.0D0
	VAPOR(I) = STREAM(I)
	ENDIF
	END IF
	ENDDO
	RETURN
	ELSEIF(MCOUNT.GE. 10)THEN        !MUST BE ALL LIQUID
	LIQUID = STREAM
	VAPOR = 0.0D0
	RETURN
	ELSEIF(MVCOUNT.GE. 10)THEN       !MUST BE ALL VAPOR
	LIQUID = 0.0D0
	VAPOR = STREAM
	RETURN
	ELSE                    !NOT CONVERGED..ADJUST ALPHA..TRY AGAIN
	DO I = 1, KL
	IF(ALPHA(I) / OLDALF(I).GT. 5.0D0)THEN
	ALPHA(I) = OLDALF(I) * 5.0D0
	ELSE IF(ALPHA(I) / OLDALF(I).LT. 0.2D0)THEN
	ALPHA(I) = OLDALF(I) * 0.2D0
	ENDIF
	ENDDO
	ENDIF
	END DO
	IER = -105  !OUT OF ITERATIONS
	RETURN
	END



	SUBROUTINE UKBMODL(KCODE, TEMP, PRES, A, KB, X, ALPHA, DHLMX, DHVMX, IER,
		&EOS_Type)

	C  Initialize or update Kb model coefficients, relative volatilities.
	C  Copyright(c) 1984, 87 by R.A.Russell, Deerhaven Technical Software.

	USE COMPDATA

	IMPLICIT REAL * 8 (A - H, O - Z), INTEGER * 4 (I - N)
	DIMENSION X(KL), Y(KL), ALPHA(KL)
	REAL * 8 KB, K(KL)
	INTEGER * 4 EOS_Type
	C  Get K's at current T,P.

	IER = 0
	IF(KB.LE. 0.0D0)THEN
	IER = -131
	RETURN
	ENDIF
	IF(KCODE.NE.0)THEN
	DO 5 I = 1, KL
	5    Y(I) = ALPHA(I)*X(I)
	END IF
	CALL UKVAL(TEMP, PRES, X, Y, KCODE, K, DHLMX, DHVMX, IER, EOS_Type)
	IF(IER.NE. 0)RETURN

	C  Set volatilities to K / Kb.

	DO 20 I = 1, KL
	20 ALPHA(I) = K(I) / KB
	A = DLOG(KB)

	RETURN
	END


	SUBROUTINE UKVAL(TEMP, PRES, X, Y, XYCODE, K, DHLMX, DHVMX, IER, EOS_Type)
	C  Calculate K values.
	C  Copyright(c) 1984, 87 by R.A.Russell, Deerhaven Technical Software.
	C
	USE COMPDATA

	IMPLICIT REAL * 8 (A - H, O - Z), INTEGER * 4 (I - N)

	INTEGER * 4 XYCODE, EOS_Type
	REAL * 8 K(KL), LNP
	DIMENSION X(KL), Y(KL), XNORM(KL)
	IER = 0
	IF(PRES.LE. 0.0D0)THEN
	IER = -141
	RETURN
	ENDIF
	LNP = DLOG(PRES)

	C  For each component that has tabular data or coefficients, interpolate
	C  for vapor pressure at the temperature.

	K = 0.0D0
	LNNU = 0.0D0
	LNGAM = 0.0D0
	LNPHI = 0.0D0

	C  liquid fugacity.

	IF(XYCODE.EQ.0)THEN
	DO I = 1, KL
	LNNU(I) = DLOG(PC(I)*14.696) - LNP + 5.42D0 *
	&                (1.0D0 + ACEN(I)) * (1.0D0 - (TC(I)*1.8) / TEMP)
	ENDDO
	ELSE
	SUM = 0.0D0
	DO I = 1, KL                  !Get normalized composition
	SUM = SUM + X(I)
	ENDDO
	IF(SUM.LE.0.0D0)THEN
	DO I = 1, KL
	LNNU(I) = DLOG(PC(I)*14.696) - LNP + 5.42D0 * (1.0D0 +
		&ACEN(I)) *(1.0D0 - (TC(I)*1.8) / TEMP)
	ENDDO
	ELSE
	DO I = 1, KL
	XNORM(I) = 0.0D0
	IF(X(I).GT.1.0D - 25)XNORM(I) = X(I) / SUM
	ENDDO

	C  Liquid activity coefficient.

	IPHASE = 1
	If(EOS_Type.eq.0) Then
	CALL USOAVE(IPHASE, TEMP, PRES, XNORM, K, DHLMX, IER, EOS_TYPE)
	Else
	CALL UPREOS(IPHASE, TEMP, PRES, XNORM, K, DHLMX, IER, EOS_TYPE)
	End If

	IF(IER.NE. 0)RETURN
	LNNU = K
	ENDIF
	C  If any vapor fugacity model, get normalized composition.If none,
	C  fugacity coefficients default to unity.

	60    SUM = 0.0D0
	DO 65 I = 1, KL
	65    SUM = SUM + Y(I)
	IF(SUM.LE. 0.0D0)GOTO 80
	DO I = 1, KL
	XNORM(I) = 0.0D0
	IF(Y(I).GT.1.0D - 25)XNORM(I) = Y(I) / SUM
	ENDDO

	C  vapor fugacity.

	IPHASE = 2
	If(EOS_Type.eq.0) Then
	CALL USOAVE(IPHASE, TEMP, PRES, XNORM, K, DHVMX, IER, EOS_TYPE)
	Else
	CALL UPREOS(IPHASE, TEMP, PRES, XNORM, K, DHVMX, IER, EOS_TYPE)
	End If
	IF(IER.NE. 0)RETURN
	LNPHI = K
	ENDIF

	C  Get ln(K) = ln(liq.fug.) + ln(liq.act.) - ln(vap.fug.), trace if called for

	80 DO I = 1, KL
	K(I) = LNNU(I) + LNGAM(I) - LNPHI(I)
	IF(K(I).LT. - 35.0D0)K(I) = -35.0D0
	IF(K(I).GT.35.0D0)K(I) = 35.0D0
	K(I) = DEXP(K(I))
	ENDDO

	RETURN
	END


	SUBROUTINE USOAVE(IPHASE, TEMP, PRES, X, FUG, DH, IER, EOS_Type)
	C  Soave - Redlich - Kwong fugacity coefficients
	C  Copyright(c) 1984, 88, by R.A.Russell, Deerhaven Technical Software.

	USE COMPDATA
	IMPLICIT REAL * 8 (A - H, O - Z), INTEGER * 4 (I - N)

	DIMENSION X(KL), FUG(KL), ACALF(KL), ACDADT(KL)
	INTEGER * 4 EOS_Type
	C  Get the 'b' and start of 'a' summations for the mixture.
	BM = 0.0D0
	S = 0.0D0
	SH = 0.0D0
	DO I = 1, KL
	BM = BM + X(I) * SRK(KL + I)
	IF(I.EQ.i_H2)THEN       !Hydrogen is treated specially.
	ACALF(I) = SRK(I)*1.0963576D0*DEXP(-0.15114D0*TEMP / 75.0D0)
	S = S + X(I)*ACALF(I)
	ELSE
	IF(TEMP.LT. 0.0D0)THEN
	IER = -151
	RETURN
	ENDIF
	SQRTR = DSQRT(TEMP / (TC(I)*1.8))
	ACALF(I) = SRK(I)*(1. + SRK((KL * 2) + I)*(1.0D0 - SQRTR))
	S = S + X(I)*ACALF(I)
	ACDADT(I) = -SRK(I)*SRK(KL * 2 + I)*0.5*SQRTR / TEMP
	SH = SH + X(I)*ACDADT(I)
	ENDIF
	ENDDO
	SH = SH * S
	C  Get the 'a' summation and initialize the 'Q' terms for fugacities.

	AM = S * *2.0D0
	DO 10 I = 1, KL
	10   FUG(I) = ACALF(I)*S

	C  Calc.the 'A' and 'B' terms in the cubic for 'Z' and solve for Z.
	25 AA = AM * PRES / TEMP * *2.0D0
	BB = BM * PRES / TEMP
	CALL UCUBIC2(AA, BB, Z, IPHASE, IER, EOS_Type)
	IF(IER.NE. 0)RETURN
	C
	C  Calculate Fugacities and Enthalpy Departure
	C
	IF(Z.EQ. 0.0D0.OR.BB.EQ. 0.0D0)THEN
	IER = -152
	RETURN
	ENDIF
	C0 = AA * DLOG(1.0D0 + BB / DABS(Z)) / BB
	IF(BM.EQ. 0.0D0)THEN
	IER = -153
	RETURN
	ENDIF
	C1 = (DABS(Z) - 1.0D0 + C0) / BM
	DH = -(ABS(Z) - 1.0 - (1.0 - 2.0*TEMP*SH / AM)*AA*LOG(1.0 + BB / ABS(Z)) / BB)*
	&        1.987*TEMP

	C  Watch out for Z less than B, due to fudged Z from CUBIC.
	c(try different handling of z<b)

	C2 = DLOG(DABS(DABS(Z) - BB))
	c     C2 = -4.6051702D0
	c     IF(DABS(Z) - BB.GT.0.)C2 = DLOG(DABS(Z) - BB)
	IF(AM.EQ. 0.0D0)THEN
	IER = -154
	RETURN
	ENDIF
	C3 = C0 * 2.0D0 / AM
	DO 30 I = 1, KL
	30    FUG(I) = C1 * SRK(KL + I) - C2 - C3 * FUG(I)
	RETURN
	END

	SUBROUTINE UPREOS(IPHASE, T, P, X, FUG, DH, IER, EOS_Type)
	!PR EOS fugacity coefficients

	USE COMPDATA
	IMPLICIT REAL * 8 (A - H, O - Z), INTEGER * 4 (I - N)

	DIMENSION X(KL), FUG(KL), ACALF(KL), ACDADT(KL)
	INTEGER * 4 EOS_Type
	!Variables
	Real(8), Intent(In) ::T, P

	Integer i, j
	Real(8) Ln_Temp, Fug_Temp1, Fug_Temp2, Fug_Temp3, Fug_Temp4, Ln_Phi
	Real(8) SUMA, SUMB, Omega_A, Omega_B, m0, m1, m2, MKi, MKj

	Real(8) A(KL), B(KL), AR(KL), BR(KL), AA, BB, AAR, BBR
	Real(8) AlphaC(KL), MK(KL), AT(KL), TR(KL), PR(KL)
	!Local Variables
	Real(8) SUMAA, SUMBB, SUMAR, SUMBR, SUMAT, TK, PATM


	TK = T / 1.8d0
	PATM = P / 14.696d0

	!RGAS = 8.3144598          !Gas constant in KJ / Kmol / K

	!EOS A and B definitions
	!Peng - Robinson EOS parameters
	Omega_A = 0.45724d0
	Omega_B = 0.0778d0
	m0 = 0.37464d0
	m1 = 1.54226d0
	m2 = -0.26992d0

	Do i = 1, KL
	!EOS A and B definitions
	TR(i) = TK * TCinv(i)
	PR(i) = PATM * PCinv(i)

	MK(i) = m0 + m1 * ACEN(i) + m2 * ACEN(i) ** 2d0
	AT(i) = (1d0 + MK(i) * (1d0 - TR(i) ** 0.5)) ** 2d0

	AlphaC(i) = Omega_A * RGAS ** 2d0 * Tc(i) ** 2d0 * PCinv(i)
	A(i) = AlphaC(i) * AT(i)
	B(i) = Omega_B * RGAS * Tc(i) * PCinv(i)
	AR(i) = A(i) * PATM / (RGAS * TK) ** 2d0
	BR(i) = B(i) * PATM / (RGAS * TK)
	End do

	!One fluid mixing rule
	SUMAR = 0d0
	SUMBR = 0d0
	SUMAA = 0d0
	SUMBB = 0d0
	SUMAT = 0d0
	Do i = 1, KL
	SUMBB = SUMBB + x(i) * B(i)
	SUMBR = SUMBR + x(i) * BR(i)
	Do j = 1, KL
	SUMAA = SUMAA + x(j) * x(i) * (A(j) * A(i)) ** 0.5d0
	SUMAR = SUMAR + x(j) * x(i) * (AR(j) * AR(i)) ** 0.5d0

	!SUMBB = SUMBB + x(j) * x(i) * (B(j) + B(i)) / 2
	!SUMBR = SUMBR + x(j) * x(i) * (BR(j) + BR(i)) / 2
	!
	!Fugacity dA / dT term calculation
	MKi = MK(i) * (A(j) * Tc(i) * PCinv(i)) ** 0.5d0
	MKj = MK(j) * (A(i) * Tc(j) * PCinv(j)) ** 0.5d0
	SUMAT = SUMAT + x(j) * x(i) * (MKi + MKj)
	End do
	End do

	AAR = SUMAR
	BBR = SUMBR
	AA = SUMAA
	BB = SUMBB

	DADT = -(RGAS / 2d0) * (Omega_A / TK) ** 0.5d0 * SUMAT

	!'A' and 'B' terms in the cubic for 'Z' and solve for Z.
	CALL UCUBIC2(AAR, BBR, Z, IPHASE, IER, EOS_Type)
	IF(IER.NE. 0)RETURN

	Do i = 1, KL
	SUMA = 0d0
	SUMB = 0d0
	Do j = 1, KL
	SUMA = SUMA + x(j) * (A(i) * A(j)) ** 0.5d0
	SUMB = SUMB + x(j) * (B(i) + B(j)) / 2d0
	End do

	!Fug_Temp1 = ((2 * SUMB - BB) / BB) * (Z - 1)
	Fug_Temp1 = (B(i) / BB) * (Z - 1d0)
	If(Z - BBR < 0d0) Then
	Fug_Temp2 = Log(1E-30)
	Else
	Fug_Temp2 = Log(Z - BBR)
	End If

	Fug_Temp3 = (2d0 * SUMA / AA) - (B(i) / BB)
	!Fug_Temp3 = (2 * SUMA / AA) - ((2 * SUMB - BB) / BB)
	Ln_Temp = (Z + BBR * (1d0 + 2d0**0.5d0)) / (Z + BBR * (1d0 - 2d0**0.5d0))
	If(Ln_Temp < 0d0) Ln_Temp = 1E-30
	Fug_Temp4 = Log(Ln_Temp) * AAR / (BBR * 8d0 ** 0.5d0)

	FUG(i) = Fug_Temp1 - Fug_Temp2 - Fug_Temp3 * Fug_Temp4
	End do

	Ln_Temp = (Z + BBR * (1d0 + 2d0**0.5d0)) / (Z + BBR * (1d0 - 2d0**0.5d0))
	If(Ln_Temp < 0) Ln_Temp = 1E-30
	DH = Log(Ln_Temp)*AA / (BB * 8d0**0.5d0)*(1d0 - TK / AA * DADT)
	&      -RGAS * TK*(Z - 1d0)
	DH = DH * 101.325d0 / 4.184d0 *1.8d0    !atm to kpa 101.325 kcal to kj 4.184

	END

	!SUBROUTINE UCUBIC(A, B, Z, PHASE, IER)
	!C  Solve cubic equation for compressibility, using Gundersen's method.
	!C  Copyright(c) 1984 by R.A.Russell, Deerhaven Technical Software.
	!
	!IMPLICIT REAL * 8 (A - H, O - Z), INTEGER * 4 (I - N)
	!INTEGER * 4 PHASE
	!F(ARG) = (ARG*ARG*ARG) - (ARG*ARG) + (ARG*Q) - R
	!
	!IER = 0
	!C  Get the coefficients 'Q' and 'R' and see if only one root with no
	!C  maximum or minimum.
	!Q = A - B - B * *2.0D0
	!R = A * B
	!IF(Q.GT.1.0D0 / 3.0D0)THEN
	!IF(F(1.0D0 / 3.0D0).GT.0.0D0)THEN
	!Z = 0.0D0
	!ELSE
	!Z = 1.0D0
	!END IF
	!GO TO 10
	!END IF
	!C  For the vapor phase, 'R' greater than 1 / 27 means one root, greater
	!C  than 1 / 3.  Otherwise find where the minimum and use it if no vapor -
	!C  like root.
	!IF((1.0D0 - 3.0D0*Q).LE. 0.0D0)THEN
	!IER = -171
	!RETURN
	!ENDIF
	!IF(PHASE.NE.1)THEN
	!C      vapor phase
	!IF(R.GT.1.0D0 / 27.0D0)THEN
	!Z = 1.0D0
	!GO TO 10
	!ELSE
	!ZM = (1.0D0 + DSQRT(1.0D0 - 3.0D0*Q)) / 3.0D0
	!IF(F(ZM).GT.0.0D0)GO TO 15
	!Z = 1.0D0
	!GO TO 10
	!END IF
	!ELSE
	!C  For the liquid phase, if the maximum is negative, use Z there.
	!ZM = (1.0D0 - DSQRT(1.0D0 - 3.0D0*Q)) / 3.0D0
	!IF(F(ZM).LT.0.0D0)GO TO 15
	!Z = 0.0D0
	!END IF
	!C
	!C  Solve for Z using Newton's method.  Iterate until no change in Z.
	!10 Z1 = 0.0D0
	!Z2 = 0.0D0
	!11 IF((Q + (3.0D0*Z - 2.0D0)*Z).EQ. 0.0D0)THEN
	!IER = -172
	!RETURN
	!ENDIF
	!Z = Z - F(Z) / (Q + (3.0D0*Z - 2.0D0)*Z)
	!IF(DABS(Z - Z1).GT. 1.0D - 10.OR.DABS(Z - Z2).GT. 1.0D - 10)THEN
	!Z2 = Z1
	!Z1 = Z
	!GO TO 11
	!END IF
	!C
	!C  Restrict vapor Z to 0.25 min., liquid Z to .75 max.
	!IF(PHASE.NE.1)THEN
	!C      vapor phase
	!IF(Z.LT.0.25D0)Z = -0.25D0
	!ELSE
	!C      liquid phase
	!IF(Z.GT.0.75D0)Z = -0.75D0
	!END IF
	!RETURN
	!
	!C  Correct 'B' if max or min was chosen as solution.
	!15 Z = -ZM
	!c  Suspend correction for now; doesn't always work well.  See SOAVE.
	!IF(PHASE.NE.1)THEN
	!C      vapor phase
	!B = B * (1. + F(ZM) / R)
	!ELSE
	!C      liquid phase
	!B = B * (1. + F(ZM) / R)
	!END IF
	!RETURN
	!END

	SUBROUTINE UCUBIC2(A, B, Z, PHASE, IER, EOS_Type)
	!7 / 12 / 2015
	!Solve cubic equation for compressibility analytically
	!By Solomon Gebreyohannes
	!
	IMPLICIT REAL * 8 (A - H, O - Z), INTEGER * 4 (I - N)
	INTEGER * 4 PHASE
	F(ARG) = (ARG*ARG*ARG) - (ARG*ARG) + (ARG*Q) - R
	REAL * 8 k1, k2, k3, r1, r2, r3
	INTEGER * 4 RT, EOS_Type
	IER = 0
	C  Get the coefficients 'Q' and 'R' and see if only one root with no
	C  maximum or minimum.
	Q = A - B - B * *2.0D0
	R = A * B

	!k1 = -1
	!k2 = Q
	!k3 = -R

	If(EOS_Type == 1) Then
	!Peng - Robinson EOS
	K1 = B - 1d0
	K2 = A - (2d0 * B) - (3d0 * B ** 2d0)
	K3 = (B ** 3d0) + (B ** 2d0) - (A * B)
	ElseIf(EOS_Type == 0) Then
	!Soave - Redlich - Kwong EOS
	K1 = -1d0
	K2 = A - B - B * * 2d0
	K3 = -A * B
	End If

	Call Cubic_Solver(k1, k2, k3, r1, r2, r3, RT)

	If(RT.EQ.1) Then
	Z = r1
	ELSEIF(RT.EQ.2) THEN
	IF(PHASE.EQ.1) THEN
	Z = MAX(r1, r2)
	IF(r1.GT.0D0.AND.r1.LT.Z) Z = r1
	IF(r2.GT.0D0.AND.r2.LT.Z) Z = r2
	END IF

	IF(PHASE.EQ.2)  Z = MAX(r1, r2)
	ELSE
	IF(PHASE.EQ.1) THEN
	Z = MAX(r1, r2, r3)
	IF(r1.GT.0D0.AND.r1.LT.Z) Z = r1
	IF(r2.GT.0D0.AND.r2.LT.Z) Z = r2
	IF(r3.GT.0D0.AND.r3.LT.Z) Z = r3
	END IF

	IF(PHASE.EQ.2.D0) Z = MAX(r1, r2, r3)
	End If

	!Incase there is a loss of numerical precision, Use Newton's method to solve for Z. Iterate until no change in Z.
	!If(abs(f(Z)).GT.1.0D - 10) Then
	!10         Z1 = 0.0D0
	!Z2 = 0.0D0
	!DO 11 i = 1, 50
	!IF((Q + (3.0D0*Z - 2.0D0)*Z).EQ. 0.0D0)THEN
	!IER = -172
	!Exit
	!END IF
	!
	!Z = Z - F(Z) / (Q + (3.0D0*Z - 2.0D0)*Z)
	!IF(DABS(Z - Z1).GT. 1.0D - 10.OR.DABS(Z - Z2).GT. 1.0D - 10)THEN
	!Z2 = Z1
	!Z1 = Z
	!ELSE
	!Exit
	!END IF
	!11       CONTINUE
	!ENDIF

	IF(abs(f(Z)).GT.0.1)THEN
	IER = -172
	END IF

	C  Restrict liquid Z to B min.
	IF(PHASE.EQ.1.and.Z.LT.0) Z = B

	C  Restrict vapor Z to 0.25 min.
	IF(PHASE.EQ.2.and.Z.LT.0.25) Z = 0.25

	RETURN
	END

	SUBROUTINE Cubic_Solver(k1, k2, k3, r1, r2, r3, RT)
	!7 / 12 / 2015
	!Cubic equation solver(Analytically)
	!By Solomon Gebreyohannes
	!
	!Equation x ^ 3 + k1 * x ^ 2 + k2 * x + k3 = 0
	!Inputs
	!Coefficients: k1, k2, k3
	!Outputs
	!Root Type(RT)
	!RT = 1   One real[r1]; Two complex[real part, r1; imaginary part, r3]
	!RT = 2   Two real[r1, r2]
	!RT = 3   Three real[r1, r2, r3]
	REAL * 8 k1, k2, k3, r1, r2, r3
	INTEGER * 4 RT
	REAL * 8 A, B, D, M, N, Phi, M_T, N_T, Pi_Val

	A = 1.D0 / 3.D0*(3.D0*k2 - k1 * *2.D0)
	B = 1.D0 / 27.D0*(2.D0*k1**3 - 9.D0*k1 * k2 + 27.D0*k3)
	D = A * *3.D0 / 27.D0 + B * *2.D0 / 4.D0        !discriminant

	If(D.GT.0.D0) Then
	!Root type 1
	RT = 1
	M_T = -B / 2.D0 + Sqrt(D)
	If(M_T.GE.0.D0) Then
	M = M_T * * (1.D0 / 3.D0)
	Else
	M = -Abs(M_T) ** (1.D0 / 3.D0)
	End If

	N_T = -B / 2.D0 - Sqrt(D)
	If(N_T.GE.0.D0) Then
	N = N_T * * (1.D0 / 3.D0)
	Else
	N = -Abs(N_T) ** (1.D0 / 3.D0)
	End If

	!Real root
	r1 = M + N - k1 / 3.D0
	!Complex roots
	r2 = -1.D0 / 2.D0 * (M + N) - k1 / 3.D0      !Real
	r3 = Sqrt(3.D0) / 2.D0 * (M - N)             !Imaginary

	ElseIf(D.EQ.0.D0) Then
	!Root type 2
	RT = 2
	M_T = -B / 2.D0 + Sqrt(D)
	If(M_T.GE.0.D0) Then
	M = M_T * * (1.D0 / 3.D0)
	Else
	M = -Abs(M_T) ** (1.D0 / 3.D0)
	End If

	N_T = -B / 2.D0 - Sqrt(D)
	If(N_T.GE.0.D0) Then
	N = N_T * * (1.D0 / 3.D0)
	Else
	N = -Abs(N_T) ** (1.D0 / 3.D0)
	End If

	!Real roots
	r1 = M + N - k1 / 3.D0
	r2 = -0.5D0 * (M + N) - k1 / 3.D0
	r3 = -0.5D0 * (M + N) - k1 / 3.D0

	ElseIf(D.LT.0.D0) Then
	!Root type 3
	RT = 3
	Phi = 1.D0 / 3.D0 * Acos(-0.5D0*B / Sqrt(-A * *3.D0 / 27.D0))
	Pi_Val = 4.D0*DATAN(1.D0)

	!Real roots
	r1 = 2.D0*Sqrt(-A / 3.D0)*Cos(Phi - 2.D0*Pi_Val*0.D0 / 3.D0) - k1 / 3.D0
	r2 = 2.D0*Sqrt(-A / 3.D0)*Cos(Phi - 2.D0*Pi_Val*1.D0 / 3.D0) - k1 / 3.D0
	r3 = 2.D0*Sqrt(-A / 3.D0)*Cos(Phi - 2.D0*Pi_Val*2.D0 / 3.D0) - k1 / 3.D0
	End If

	END SUBROUTINE

	SUBROUTINE UCURVEF(V, E, NP, B, C, RCODE, DELTA, VMIN, VMAX, MINFLG,
		+MAXFLG, IER)
	C  Curve fit of 3 points with prediction of root.
	C  Copyright(c) 1984 by R.A.Russell, Deerhaven Technical Software Co.

	IMPLICIT REAL * 8 (A - H, O - Z), INTEGER * 4 (I - N)
	INTEGER * 4 RCODE
	DIMENSION V(3), E(3)
	IER = 0
	RCODE = 0
	C  First, note if current 'V' is new upper or lower limit.
	IF(E(1).GT.0.0D0)THEN
	IF(V(1).LE.VMIN)THEN
	RCODE = 1
	RETURN
	END IF
	VMAX = V(1)
	MAXFLG = 1
	ELSE
	IF(V(1).GE.VMAX)THEN
	RCODE = 1
	RETURN
	END IF
	VMIN = V(1)
	MINFLG = 1
	END IF
	C
	C  If the current point was predicted from three points, and the error
	C  was not reduced to half the previous lowest, and points above and
	C  below the answer have been used, then go halfway between bounds.
	IF(NP.GT.3.AND.B.NE.0.0D0)THEN
	IF(DABS(E(1)).GT.DMIN1(DABS(E(2)), DABS(E(3)))*0.5D0.AND.
		1	  MINFLG.NE.0.AND.MAXFLG.NE.0)THEN
	VNEW = (VMIN + VMAX)*0.5D0
	C = 0.0D0
	GO TO 55
	END IF
	END IF
	C
	C  For one point, use 'B' (= 1 / slope) if given.Otherwise use the given
	C  change size if given.As a last resort, use a 1 percent change.
	IF(NP.GT.1.AND.E(1).NE.E(2))GO TO 25
	IF(B.NE.0.0D0)GO TO 40
	10 IF(DELTA.NE.0.0D0)THEN
	DV = DELTA
	GO TO 20
	ELSE IF(V(1).EQ.0.0D0)THEN
	DV = 0.01D0*(VMAX - VMIN)
	GO TO 20
	ELSE
	DV = 0.01D0*V(1)
	END IF
	20 DV = DSIGN(DV, -E(1))
	GO TO 50
	C
	C  For two or more points, calc.the 'alpha' and 'beta' terms for the
	C  first two points.
	25 IF(E(1).EQ.E(2))THEN
	IER = -181
	RETURN
	ENDIF
	ALPHA2 = (V(1) - V(2)) / (E(1) - E(2))
	BETA2 = (V(1)*E(1) - V(2)*E(2)) / (E(1) - E(2))
	IF(NP.LT.3.OR.E(1).EQ.E(3))GO TO 30
	C
	C  For 3 points, calc.the 'alpha' and 'beta' terms for the first and
	C  third points and from there the 'C' coefficient.
	IF(E(1).NE.E(3))THEN
	ALPHA3 = (V(1) - V(3)) / (E(1) - E(3))
	BETA3 = (V(1)*E(1) - V(3)*E(3)) / (E(1) - E(3))
	IF(BETA2.NE.BETA3)C = (ALPHA2 - ALPHA3) / (BETA3 - BETA2)
	END IF
	C
	30 B = ALPHA2 + BETA2 * C
	C
	40 DV = E(1)*(V(1)*C - B)
	C
	C  If the change is in the wrong direction, use default change.
	IF(DV*E(1).GE.0.0D0)THEN
	B = 0.0D0
	C = 0.0D0
	GO TO 10
	END IF
	C
	C  Restrict size of predicted change.
	IF(DELTA.NE.0.0D0)DV = DSIGN(DMIN1(DABS(DV), DABS(DELTA)), DV)
	C
	C  Get the new variable and check against bounds.
	50 VNEW = V(1) + DV
	IF(VNEW.LE.VMIN)THEN
	IF(MINFLG.EQ.0)THEN
	VNEW = VMIN
	GO TO 60
	ELSE
	VNEW = (V(1) + VMIN)*0.5D0
	END IF
	ELSE IF(VNEW.GE.VMAX)THEN
	IF(MAXFLG.EQ.0)THEN
	VNEW = VMAX
	GO TO 60
	ELSE
	VNEW = (V(1) + VMAX)*0.5D0
	END IF
	END IF
	C
	C  If the new value is still at a bound, there is a precision problem.
	55 IF(VNEW.LE.VMIN.OR.VNEW.GE.VMAX.OR.VNEW.EQ.V(1))THEN
	RCODE = 1
	RETURN
	END IF
	C  Save the most recent two points and install the new value.
	60 V(3) = V(2)
	V(2) = V(1)
	E(3) = E(2)
	E(2) = E(1)
	V(1) = VNEW
	RETURN
	END


		SUBROUTINE NRSolve(xin, x, n, check, UserFunction, JacobianFlag,
			*UserJacobian, fnorm, xbounds, TOLF, ierr)
		implicit none
		C	USES fdjac, fmin, lnsrch, lubksb, ludcmp
		C	Given an initial guess x(1:n) for a root in n dimensions, find the root by a globally
		C	convergent Newton's method. The vector of functions to be zeroed, called fvec(1:n)
		C	in the routine below, is returned by a user - supplied subroutine that must be called funcv
		C	and have the declaration subroutine funcv(n, x, fvec).The output quantity check
		C	is false on a normal return and true if the routine has converged to a local minimum of the
		C	function fmin defined below.In this case try restarting from a different initial guess.
		C	Parameters : NP is the maximum expected value of n; MAXITS is the maximum number of
		C	iterations; TOLF sets the convergence criterion on function values; TOLMIN sets the criterion
		C	for deciding whether spurious convergence to a minimum of fmin has occurred; TOLX is
		C	the convergence criterion on.x; STPMX is the scaled maximum step length allowed in line
		C	searches.
		INTEGER n, nn, MAXITS, ierr
		LOGICAL check, JacobianFlag
		REAL * 8 x(n), xin(n), fvec(n), TOLF, TOLMIN, TOLX, STPMX, xbounds(2, n)
		C	PARAMETER(MAXITS = 100, TOLF = 1.0d - 4, TOLMIN = 1.0d - 6, TOLX = 1.0d - 7,
			C     *			STPMX = 100.0d0)
		C	PARAMETER(MAXITS = 100, TOLF = 1.0d - 8, TOLMIN = 1.0d - 10, TOLX = 1.0d - 20,
			C     *			STPMX = 50.0d0)
		INTEGER i, its, j, indx(n)
		REAL * 8 d, den, f, fold, stpmax, sum, temp, test, fjac(n, n),
		*g(n), p(n), xold(n), fmin, fnorm, fjac2(n, n)

		EXTERNAL fmin, UserFunction, UserJacobian
		MAXITS = 100
		TOLMIN = 1.0D - 10
		TOLX = 1.0D - 20
		STPMX = 50.0d0
		c	check that the number of equations is less than the max number
		ierr = 0
		fjac = 0.0d0
		g = 0.0d0
		p = 0.0d0
		fvec = 0.0d0
		nn = n
		x = xin
		xold = 0.0d0
		C	The vector fvec is also computed by this call.
		f = fmin(x, fvec, n, UserFunction)
		fnorm = 0.0d0
		test = 0.0d0
		C	Test for initial guess being a root.Use more stringent
		C	test than simply TOLF.
		do i = 1, n
			if (abs(fvec(i)).gt.test)test = abs(fvec(i))
				enddo
				fnorm = test
				if (test.lt..01*TOLF)then
					check = .false.
					return
					endif
					sum = 0. !Calculate stpmax for line searches.
					do i = 1, n
						sum = sum + x(i)**2
						enddo
						sum = dsqrt(sum)
						stpmax = STPMX * sum / 100.0d0
						c	stpmax = STPMX * max(sum, float(n))

						c	if (n.eq.3) then
						c      	WRITE(66, *) ' *** Convergence in NR method ***'
						c		write(66, *) ' its,x(3),fvec(3)'
						c	endif

						do its = 1, MAXITS !Start of iteration loop.


							c		if (n.eq.3) write(66, '(i3,6x,6g12.4)') its, x, fvec


							c		If analytic Jacobian is available, you can replace the routine
							c		fdjac below with your own routine.
							if (JacobianFlag) then
								call UserJacobian(n, x, fjac)
								c			call fdjac(n, x, fvec, fjac2, UserFunction)
								c
								c			if (its.eq.1) then
								c      	WRITE(66, *) ' *** Jacobian comparison ***'
								c		write(66, *) ' j,i,analytic,delta'
								c				do i = 1, n !Compute rf for the line search.
								c					do j = 1, n
								c		   write(66, '(i3,6x,I3,6X,2g12.4)') J, I, fjac(J, i), fjac2(J, I)
								c					enddo
								c				enddo
								c			endif

								else
								call fdjac(n, x, fvec, fjac, UserFunction)
								endif
								do i = 1, n !Compute rf for the line search.
									sum = 0.
									do j = 1, n
										sum = sum + fjac(j, i)*fvec(j)
										enddo
										g(i) = sum
										enddo
										c		Store x and f
										xold = x
										fold = f
										c		Right - hand side for linear equations.
										p = -fvec
										c		Solve linear equations by LU decomposition.
										call ludcmp(fjac, n, n, indx, d, ierr)
										if (ierr.ne.0) return
											call lubksb(fjac, n, n, indx, p)
											call lnsrch(n, xold, fold, g, p, x, f, fvec, stpmax, check,
												&fmin, UserFunction, xbounds, ierr)
											if (ierr.ne.0) return
												c		lnsrch returns new x and f.It also calculates fvec at the new x when it calls fmin.


												c		if (n.eq.3) write(66, '(a11,6g12.4)') ' after ls: ', x, fvec


												test = 0.
												fnorm = 0.0d0
												c		Test for convergence on function values.
												do i = 1, n
													if (abs(fvec(i)).gt.test)test = abs(fvec(i))
														enddo
														fnorm = test


														c		if (n.eq.3) write(66, '(a8,6g12.4)') ' test = ', test


														if (test.lt.TOLF)then
															check = .false.
															return
															endif
															if (check)then
																c			Check for gradient of f zero, i.e., spurious convergence.
																test = 0.
																den = max(f, .5*n)
																do i = 1, n
																	temp = abs(g(i))*max(abs(x(i)), 1.) / den
																	if (temp.gt.test)test = temp
																		enddo
																		if (test.lt.TOLMIN)then
																			check = .true.
																		else
																			check = .false.
																			endif
																			return
																			endif
																			c		Test for convergence on del - x
																			test = 0.
																			do i = 1, n
																				temp = (abs(x(i) - xold(i))) / max(abs(x(i)), 1.)
																				if (temp.gt.test)test = temp
																					enddo
																					if (test.lt.TOLX) then
																						return
																						endif
																						enddo
																						c	MAXITS exceeded in newt
																						ierr = 200
																						return
																						END


																						SUBROUTINE ludcmp(a, n, np, indx, d, ierror)
																						IMPLICIT NONE
																						INTEGER n, np, indx(n), NMAX, ierror
																						REAL * 8 d, a(np, np), TINY
																						C	Largest expected n, and a small number.
																						C	PARAMETER(NMAX = 500, TINY = 1.0d - 20)
																						C	Given a matrix a(1:n, 1 : n), with physical dimension np by np, this routine replaces it by
																						C	the LU decomposition of a rowwise permutation of itself.a and n are input.a is output,
																						C	arranged as in equation(2.3.14) above; indx(1:n) is an output vector that records the
																						C	row permutation effected by the partial pivoting; d is output as ±1 depending on whether
																						C	the number of row interchanges was even or odd, respectively.This routine is used in
																						C	combination with lubksb to solve linear equations or invert a matrix.
																						INTEGER i, imax, j, k
																						REAL * 8 aamax, dum, sum, vv(500)
																						NMAX = 500
																						TINY = 1.0d - 20
																						C	vv stores the implicit scaling of each row.
																						d = 1. !No row interchanges yet.
																						C	Loop over rows to get the implicit scaling information.
																						do i = 1, n
																							aamax = 0.
																							do j = 1, n
																								if (abs(a(i, j)).gt.aamax) aamax = abs(a(i, j))
																									enddo
																									if (aamax.eq.0.) then
																										C			'singular matrix in ludcmp' No nonzero largest element.
																										ierror = 1
																										return
																										endif
																										vv(i) = 1. / aamax
																										C		Save the scaling.
																										enddo
																										do j = 1, n
																											C		This is the loop over columns of Crout’s method.
																											do i = 1, j - 1
																												C			This is equation(2.3.12) except for i = j.
																												sum = a(i, j)
																												do k = 1, i - 1
																													sum = sum - a(i, k)*a(k, j)
																													enddo
																													a(i, j) = sum
																													enddo
																													C		Initialize for the search for largest pivot element.
																													aamax = 0.
																													do i = j, n
																														C			This is i = j of equation(2.3.12)
																														C			and i = j + 1...N of equation(2.3.13).
																														sum = a(i, j)
																														do k = 1, j - 1
																															sum = sum - a(i, k)*a(k, j)
																															enddo
																															a(i, j) = sum
																															C			Figure of merit for the pivot.
																															dum = vv(i)*abs(sum)
																															if (dum.ge.aamax) then
																																C				Is it better than the best so far ?
																																imax = i
																																aamax = dum
																																endif
																																enddo
																																C		Do we need to interchange rows ?
																																if ((j.ne.imax).and.(imax.gt.0)) then
																																	C			Yes, do so...
																																	do k = 1, n
																																		dum = a(imax, k)
																																		a(imax, k) = a(j, k)
																																		a(j, k) = dum
																																		enddo
																																		C			 ...and change the parity of d.
																																		d = -d
																																		C			Also interchange the scale factor.
																																		vv(imax) = vv(j)
																																		endif
																																		indx(j) = imax
																																		if (a(j, j).eq.0.)a(j, j) = TINY
																																			C		If the pivot element is zero the matrix is singular
																																			C(at least to the precision of the algorithm).
																																			C		For some applications on singular matrices, it is
																																			C		desirable to substitute TINY for zero.
																																			if (j.ne.n)then
																																				C			Now, finally, divide by the pivot element.
																																				dum = 1. / a(j, j)
																																				do i = j + 1, n
																																					a(i, j) = a(i, j)*dum
																																					enddo
																																					endif
																																					C		Go back for the next column in the reduction.
																																					enddo
																																					return
																																					END


																																					SUBROUTINE lubksb(a, n, np, indx, b)
																																					IMPLICIT NONE
																																					INTEGER n, np, indx(n)
																																					REAL * 8 a(np, np), b(n)
																																					C	Solves the set of n linear equations A · X = B.Here a is input, not as the matrix A but
																																					C	rather as its LU decomposition, determined by the routine ludcmp.indx is input as the
																																					C	permutation vector returned by ludcmp.b(1:n) is input as the right - hand side vector B,
																																					C	and returns with the solution vector X.a, n, np, and indx are not modified by this routine
																																					C	and can be left in place for successive calls with different right - hand sides b.This routine
																																					C	takes into account the possibility that b will begin with many zero elements, so it is efficient
																																					C	for use in matrix inversion.
																																					INTEGER i, ii, j, ll
																																					REAL * 8 sum
																																					ii = 0
																																					C	When ii is set to a positive value, it will become the index
																																					C	of the first nonvanishing element of b.We now do
																																					C	the forward substitution, equation(2.3.6).The only new
																																					C	wrinkle is to unscramble the permutation as we go.
																																					do i = 1, n
																																						ll = indx(i)
																																						sum = b(ll)
																																						b(ll) = b(i)
																																						if (ii.ne.0)then
																																							do j = ii, i - 1
																																								sum = sum - a(i, j)*b(j)
																																								enddo
																																						else if (sum.ne.0.) then
																																								C			A nonzero element was encountered, so from now on we will
																																							C			have to do the sums in the loop above.
																																							ii = i
																																							endif
																																							b(i) = sum
																																							enddo
																																							do i = n, 1, -1
																																								C		Now we do the backsubstitution, equation(2.3.7).
																																								sum = b(i)
																																								do j = i + 1, n
																																									sum = sum - a(i, j)*b(j)
																																									enddo
																																									C		Store a component of the solution vector X.
																																									b(i) = sum / a(i, i)
																																									enddo
																																									C	All done!
																																									return
																																									END



																																									SUBROUTINE lnsrch(n, xold, fold, g, p, x, f, fvec, stpmax,
																																										&check, func, UserFunction, xbounds, ierr)
																																									implicit none
																																									INTEGER n
																																									LOGICAL check
																																									REAL * 8 f, fold, stpmax, g(n), p(n), x(n), xold(n), func, ALF, TOLX, fvec(n)
																																									REAL * 8 xbounds(2, n), CheckBounds
																																									C	PARAMETER(ALF = 1.0d - 4, TOLX = 1.0d - 7)
																																									EXTERNAL func, UserFunction, CheckBounds
																																									C	USES func
																																									c	Given an n - dimensional point xold(1:n), the value of the function and gradient there,
																																									c	fold and g(1:n), and a direction p(1:n), finds a new point x(1:n) along the direction
																																									c	p from xold where the function func has decreased "sufficiently".The new function value
																																									c	is returned in f.stpmax is an input quantity that limits the length of the steps so that you
																																									c	do not try to evaluate the function in regions where it is undefined or subject to overflow.
																																									c	p is usually the Newton direction.The output quantity check is false on a normal exit.
																																									c	It is true when x is too close to xold.In a minimization algorithm, this usually signals
																																									c	convergence and can be ignored.However, in a zero - finding algorithm the calling program
																																									c	should check whether the convergence is spurious.
																																									c	Parameters : ALF ensures sufficient decrease in function value; TOLX is the convergence
																																									c	criterion on.x.
																																									c
																																									c	Modify to not allow steps beyond bounds.Bounds for x in xbounds(2, n) array.Elements(1, n)
																																									c	are lower bounds and elements(2, n) are upper bounds.A value of - 100.0d0 means no bound.
																																									c
																																									INTEGER i, ierr, k
																																									REAL * 8 a, alam, alam2, alamin, b, disc, f2, rhs1, rhs2, slope,
																																									*sum, temp, test, tmplam

																																									ALF = 1.0D - 4
																																									TOLX = 1.0D - 07
																																									ierr = 0
																																									check = .false.
																																									sum = 0.0d0
																																									do  i = 1, n
																																										sum = sum + p(i)*p(i)
																																										enddo
																																										sum = sqrt(sum)
																																										if (sum.gt.stpmax)then !Scale if attempted step is too big.
																																											do i = 1, n
																																												p(i) = p(i)*stpmax / sum
																																												enddo
																																												endif
																																												slope = 0.0d0
																																												do i = 1, n
																																													slope = slope + g(i)*p(i)
																																													enddo
																																													if (slope.ge.0.0d0) then
																																														c		roundoff problem in lnsrch
																																														ierr = 100
																																														return
																																														endif
																																														test = 0.0d0 !Compute.min.
																																														do i = 1, n
																																															temp = abs(p(i)) / max(abs(xold(i)), 1.0d0)
																																															if (temp.gt.test)test = temp
																																																enddo
																																																alamin = TOLX / test
																																																alam = 1.0d0 !Always try full Newton step first.
																																																do k = 1, 100
																																																	do i = 1, n
																																																		x(i) = xold(i) + alam * p(i)
																																																		x(i) = CheckBounds(x(i), xold(i), xbounds(1, i), xbounds(2, i))
																																																		enddo
																																																		c		evaluate function at this point
																																																		f = func(x, fvec, n, UserFunction)
																																																		if (isnan(f)) then
																																																			ierr = 100
																																																			return
																																																			endif
																																																			if (alam.lt.alamin)then
																																																				c			Convergence on delta - x.For zero finding,
																																																				c			the calling program should verify the convergence.
																																																				c			do i = 1, n
																																																				c				x(i) = xold(i)
																																																				c			enddo
																																																				c			check = .true.
																																																				return
																																																			else if (f.le.fold + ALF * alam*slope)then !Sufficient function decrease.
																																																				return
																																																			else  !Backtrack.
																																																				if (k.eq.1)then !First time.
																																																					tmplam = -slope / (2.0d0*(f - fold - slope))
																																																				else  !Subsequent backtracks.
																																																					rhs1 = f - fold - alam * slope
																																																					rhs2 = f2 - fold - alam2 * slope
																																																					a = (rhs1 / alam * *2 - rhs2 / alam2 * *2) / (alam - alam2)
																																																					b = (-alam2 * rhs1 / alam * *2 + alam * rhs2 / alam2 * *2) /
																																																					*(alam - alam2)
																																																					if (a.eq.0.0d0)then
																																																						tmplam = -slope / (2.0d0*b)
																																																					else
																																																						disc = b * b - 3.*a*slope
																																																						if (disc.lt.0.0d0)then
																																																							tmplam = .5*alam
																																																						else if (b.le.0.0d0)then
																																																							tmplam = (-b + dsqrt(disc)) / (3.0d0*a)
																																																						else
																																																							tmplam = -slope / (b + dsqrt(disc))
																																																							endif
																																																							endif
																																																							if (tmplam.gt.0.5d0*alam)tmplam = 0.5d0*alam
																																																								endif
																																																								endif
																																																								alam2 = alam
																																																								f2 = f
																																																								alam = max(tmplam, 0.1d0*alam)
																																																								enddo
																																																								c	if this point is reached there is a problem with line search
																																																								ierr = 100
																																																								return
																																																								END



																																																								SUBROUTINE fdjac(n, x, fvec, df, UserFunction)
																																																								implicit none
																																																								C	USES funcv
																																																								c	Computes forward - difference approximation to Jacobian.On input, x(1:n) is the point
																																																								c	at which the Jacobian is to be evaluated, fvec(1:n) is the vector of function values at
																																																								c	the point, and np is the physical dimension of the Jacobian array df(1:n, 1 : n) which is
																																																								c	output.subroutine funcv(n, x, f) is a fixed - name, user - supplied routine that returns
																																																								c	the vector of functions at x.
																																																								c	Parameters : EPS is the approximate square root of the machine precision.
																																																								INTEGER n, i, j
																																																								REAL * 8 df(n, n), fvec(n), x(n), EPS
																																																								PARAMETER(EPS = 1.0d - 4)
																																																								REAL * 8 h, temp, f(n)
																																																								external UserFunction
																																																								do j = 1, n
																																																									temp = x(j)
																																																									h = EPS * abs(temp)
																																																									if (h.eq.0.0d0)h = EPS
																																																										x(j) = temp + h !Trick to reduce finite precision error.
																																																										h = x(j) - temp
																																																										c		call funcv(n, x, f)
																																																										CALL UserFunction(X, f, N)
																																																										x(j) = temp
																																																										do i = 1, n !Forward difference formula.
																																																											df(i, j) = (f(i) - fvec(i)) / h
																																																											enddo
																																																											enddo
																																																											return
																																																											END


																																																											REAL * 8 FUNCTION fmin(x, fvec, n, UserFunction)
																																																											implicit none
																																																											C	USES funcv
																																																											c	Returns f = 1 / 2F*F at x subroutine funcv(n, x, f) is a fixed - name, user - supplied
																																																											c	routine that returns the vector of functions at x.
																																																											INTEGER n, i
																																																											REAL * 8 x(n), fvec(n), sum
																																																											external UserFunction

																																																											CALL UserFunction(x, fvec, n)
																																																											sum = 0.
																																																											do i = 1, n
																																																												sum = sum + fvec(i)**2
																																																												enddo
																																																												fmin = 0.5*sum
																																																												return
																																																												END


																																																												REAL * 8 FUNCTION CheckBounds(x, xold, xmin, xmax)
																																																												implicit none
																																																												C	Checks new x vs the bounds.If the new x is past a boundary then the value returned is
																																																												c	one - half the distance between the boundary and the old x value.Otherwise the x value
																																																												c	is returned.A boundary value of - 100. means no boundary condition set.
																																																												REAL * 8 x, xold, xmin, xmax
																																																												c	check for transgression of lower bound
																																																												if ((x.lt.xmin).and.(abs(xmin + 100.0d0).gt.0.0001d0))then
																																																													if (xold.le.xmin)then
																																																														CheckBounds = xmin
																																																													else
																																																														CheckBounds = xmin + (xold - xmin)*0.5
																																																														endif
																																																														return
																																																														endif
																																																														c	check for transgression of lower bound
																																																														if ((x.gt.xmax).and.(abs(xmax + 100.0d0).gt.0.0001d0))then
																																																															if (xold.ge.xmax)then
																																																																CheckBounds = xmax
																																																															else
																																																																CheckBounds = xold + (xmax - xold)*0.5
																																																																endif
																																																																return
																																																																endif
																																																																c	just return x, no bound transgressed
																																																																CheckBounds = x
																																																																return
																																																																END



																																																																SUBROUTINE zbrent(func, x1, x2, tol, zroot, ierror)
																																																																INTEGER ITMAX, ierror
																																																																REAL * 8 zroot, tol, x1, x2, func, EPS
																																																																EXTERNAL func
																																																																PARAMETER(ITMAX = 100, EPS = 3.e-8)
																																																																c	FROM Numerical Recipes
																																																																c	func - reference to function to find root of
																																																																c	x1, x2 - bounds on interval for function
																																																																c	tol - tolerance for solution
																																																																c	ierror - error flag, 0 = no errors
																																																																c
																																																																c	Using Brent’s method, find the root of a function func known to lie between x1 and x2.
																																																																c	The root, returned as zbrent, will be refined until its accuracy is tol.
																																																																c	Parameters : Maximum allowed number of iterations, and machine floating - point precision.
																																																																INTEGER iter
																																																																REAL * 8 a, b, c, d, e, fa, fb, fc, p, q, r, s, tol1, xm
																																																																ierror = 0
																																																																a = x1
																																																																b = x2
																																																																fa = func(a)
																																																																fb = func(b)
																																																																zroot = 0.0d0
																																																																if ((fa.gt.0..and.fb.gt.0.). or .(fa.lt.0..and.fb.lt.0.)) then
																																																																	c		root must be bracketed for zbrent
																																																																	ierror = 100
																																																																	return
																																																																	endif

																																																																	c = b
																																																																	fc = fb
																																																																	do iter = 1, ITMAX
																																																																		if ((fb.gt.0..and.fc.gt.0.). or .(fb.lt.0..and.fc.lt.0.))then
																																																																			c			Rename a, b, c and adjust bounding interval d.
																																																																			c = a
																																																																			fc = fa
																																																																			d = b - a
																																																																			e = d
																																																																			endif
																																																																			if (abs(fc).lt.abs(fb)) then
																																																																				a = b
																																																																				b = c
																																																																				c = a
																																																																				fa = fb
																																																																				fb = fc
																																																																				fc = fa
																																																																				endif
																																																																				c		Convergence check.
																																																																				tol1 = 2.*EPS*abs(b) + 0.5*tol
																																																																				xm = .5*(c - b)
																																																																				if (abs(xm).le.tol1 . or .fb.eq.0.)then
																																																																					zroot = b
																																																																					return
																																																																					endif
																																																																					c		Attempt inverse quadratic interpolation.
																																																																					if (abs(e).ge.tol1 .and.abs(fa).gt.abs(fb)) then
																																																																						s = fb / fa
																																																																						if (a.eq.c) then
																																																																							p = 2.0d0*xm*s
																																																																							q = 1.0d0 - s
																																																																						else
																																																																							q = fa / fc
																																																																							r = fb / fc
																																																																							p = s * (2.0d0*xm*q*(q - r) - (b - a)*(r - 1.0d0))
																																																																							q = (q - 1.0d0)*(r - 1.0d0)*(s - 1.0d0)
																																																																							endif
																																																																							c			Check whether in bounds.
																																																																							if (p.gt.0.0d0) q = -q
																																																																								p = abs(p)
																																																																								if (2.*p.lt.min(3.*xm*q - abs(tol1*q), abs(e*q))) then
																																																																									c				Accept interpolation.
																																																																									e = d
																																																																									d = p / q
																																																																								else
																																																																									c				Interpolation failed, use bisection.
																																																																									d = xm
																																																																									e = d
																																																																									endif
																																																																							else
																																																																									c			Bounds decreasing too slowly, use bisection.
																																																																								d = xm
																																																																								e = d
																																																																								endif
																																																																								c		Move last best guess to a
																																																																								a = b
																																																																								fa = fb
																																																																								if (abs(d).gt.tol1) then
																																																																									c			Evaluate new trial root.
																																																																									b = b + d
																																																																								else
																																																																									b = b + sign(tol1, xm)
																																																																									endif
																																																																									fb = func(b)
																																																																									enddo
																																																																									c	zbrent exceeding maximum iterations
																																																																									ierror = 101
																																																																									zroot = b
																																																																									return
																																																																									END
***/

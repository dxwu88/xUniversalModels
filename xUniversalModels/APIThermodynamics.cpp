#include "stdafx.h"
#include "Util.h"
#include "APIThermodynamics.h"

APIThermodynamics::APIThermodynamics()
{
}

APIThermodynamics::~APIThermodynamics()
{
}

//C     This subroutine calculates the temperature needed to bring the feed
//C     to InEnthalpy by performing a series of PT flash calculations.
//C
//C  VARIABLES :
//C	FeedMoles(KL)     REAL * 8		LB - MOLES / HR OF MATERIAL
//C	InEnthalpy		REAL * 8		Target ENTHALPY OF STREAM, BTU / HR
//C	temperature				REAL * 8		Stream TEMPERATURE to match enthalpy, K
//C	PRESSatm		    REAL * 8		Stream PRESSURE, ATM
//C	StartT			REAL * 8		Search starting TEMPERATURE, K
//C	VF    		    REAL * 8		Molar vapor fraction
//C     Error			INT * 4		Error flag, 0 = no errors
int APIThermodynamics::PHFlash(double streamMoles[], double InEnthalpy, double temperature, double PRESSatm, double StartT, double VF,	int& Error)
{
	double EA, EB;
	double EC = 0.0;

	Error = 0;
	double errrel = 1.0E-5;

	if (StartT == 0.0) StartT = 450.0;

	//streamMoles = FeedMoles
	double PAtm = PRESSatm;
	double TargetEnthalpy = InEnthalpy;
	double Tol = 1.0E-2;
	double A = 0.0;
	double B = 0.0;
	bool FoundLower = false;
	bool FoundUpper = false;
	double C = StartT;

	// Try at the guess temperature
	Error = CalcEnthalpy(streamMoles, C, PAtm, EC, VF);
	if (Error) return Error;

	if (abs((EC - InEnthalpy) / InEnthalpy) < errrel)
	{
		temperature = C;
		return SUCCESS;
	}
	else
	{
		if (EC > TargetEnthalpy)
		{
			EB = EC;
			B = temperature;
			FoundUpper = true;
		}
		else
		{
			EA = EC;
			A = temperature;
			FoundLower = true;
		}
	}

	// Search the interval for zero
	// Find A and B to bound the interval
	if (FoundLower)
	{
		// FIND THE UPPER BOUND OF THE INTERVAL, B
		if (B == 0.0) {
			for (int i = 0; i < 15; ++i)
			{
				B = StartT + Flash::DELTAS[i];
				Error = CalcEnthalpy(streamMoles, B, PAtm, EB, VF);
				if (Error) return Error;
						
				if (EB >= TargetEnthalpy)
				{
					FoundUpper = true;
					break;
				}
			}
		}
	}
	else
	{
		// FIND THE LOWER BOUND OF THE INTERVAL, A
		if (A == 0.0)
		{
			for (int i = 0; i < 15; ++i)
			{
				A = StartT - Flash::DELTAS[i];
				Error = CalcEnthalpy(streamMoles, A, PAtm, EA, VF);
				if (Error) return Error;

				if (EA <= TargetEnthalpy)
						{
					FoundLower = true;
					break;
				}
			}
		}
	}

	// check if problem is not bounded
	if (!FoundLower | !FoundUpper)
	{
		Error = 91;
		return Error;
	}

	// check if the interval is small(< 10°K) and if so just interpolate
	// this is faster than DZBREN
	if ((B - A) < 10.0)
	{
		temperature = A + (B - A)*(TargetEnthalpy - EA) / (EB - EA);
		Error = CalcEnthalpy(streamMoles, temperature, PAtm, EC, VF);
		if (Error) return Error;

		// do a second interpolation based on the centerpoint
		// check result and rebound interval
		if (EC < EB &&  EC > EA)
		{
			if (abs((EC - InEnthalpy) / InEnthalpy) <= errrel)
			{
				return SUCCESS;
			}
			else
			{
				if (EC > TargetEnthalpy)
				{
					EB = EC;
					B = temperature;
				}
				else
				{
					EA = EC;
					A = temperature;
				}
						
				temperature = A + (B - A)*(TargetEnthalpy - EA) / (EB - EA);
				Error = CalcEnthalpy(streamMoles, temperature, PAtm, EA, VF);
				return Error;
			}
		}
		else
		{
			// the interpolated enthalpy is outside the bounds, some
			// problem has occurred
			Error = 95;
			return Error;
		}
	}
	else
	{
		// find solution with Brent's method
		//zbrent(FunctionE, A, B, Tol, temperature, Error);
		Util::brents_fun(FunctionEnthalpy2Target, A, B, Tol, 500);

		if (Error) return Error;
		
		Error = CalcEnthalpy(streamMoles, temperature, PAtm, EC, VF);
	}

	return Error;
}

// This function is used in the search for the temperature to match the
// stream enthalpy. The TargetEnthalpy is in HeatCom
double APIThermodynamics::FunctionEnthalpy2Target(double temperature)
{
	double Enthalpy = 0.0;
	double Patm, VF;
	double TargetEnthalpy = 0.0;

	//int Error = CalcEnthalpy(streamMoles, temperature, Patm, Enthalpy, VF);
	//return Enthalpy - TargetEnthalpy;
	
	double x = temperature;
	return (x*x*x - 4 * x - 9);
}

/***
				Subroutine IdealGasEnthalpy(FeedMoles, TK, Enthalpy)
				C     This subroutine calculates the ideal gas enthalpy of the stream
				C     in units of BTU.
				C
				C  VARIABLES :
C	FeedMoles(KL)   REAL * 8		LB - MOLES / HR OF MATERIAL
C	TK				REAL * 8		TEMPERATURE, K
C	Enthalpy		REAL * 8		IDEAL GAS ENTHALPY OF STREAM, BTU / HR
C
C**********************************************************************
use CompData
IMPLICIT NONE

Real * 8  FeedMoles(KL), TK, Enthalpy,
&EMass, T2, T3, T4, T5, HFORM, EMassT, EMassT0, T0
Integer * 4  I

C   Calculate Ideal gas enthalpy(KJ / KG)
T2 = TK * TK
T3 = TK * T2
T4 = TK * T3
T5 = TK * T4
Enthalpy = 0.0
HFORM = 0.0
T0 = 298.15         !Refernce T(K)

do i = 1, KL
IF(FeedMoles(I).GT.0.0) Then
!Emass, KJ / KG
select case (HEquationType(i))
case ('POLY5')
!poly5 form of ethalpy equation(non - DIPPR components)
EMassT = AAE(i) + TK * BE(i) + T2 * CE(i) + T3 * DE(i) +
&T4*EE(i) + T5 * FE(i)

EMassT0 = AAE(i) + T0 * BE(i) + T0 * *2.0D0*CE(i) +
&T0**3.0D0*DE(i) + T0 * *4.0D0*EE(i) + T0 * *5.0D0*FE(i)

case ('DIPPR117')
!DIPPR 117 form of enthalpy equation
EMassT = AAE(I)*TK + BE(I)*CE(I)*
&                  (COSH(CE(I) / TK) / SINH(CE(I) / TK)) -
&DE(I)*EE(I)*TANH(EE(I) / TK) + FE(I)

EMassT0 = AAE(I)*T0 + BE(I)*CE(I)*
&                  (COSH(CE(I) / T0) / SINH(CE(I) / T0)) -
&DE(I)*EE(I)*TANH(EE(I) / T0) + FE(I)
end select

EMass = EMassT - EMassT0

C		  HFORM, KJ   HF25, KJ / KG - Mole
HFORM = HFORM + FEEDMOLES(I)*0.4535924*HF25(I)
Enthalpy = Enthalpy + FEEDMOLES(I)*0.4535924*MW(I)*EMass
ENDIF
ENDDO
Enthalpy = Enthalpy + HFORM
c  Convert to BTU from KJ
Enthalpy = Enthalpy / 1.055056

RETURN

END



Subroutine IGEnthalpy(FeedMoles, TK, Enthalpy)
C     This subroutine calculates the ideal gas enthalpy of the stream
C     in units of BTU. (A wrapper for IdealGasEnthalpy)
C
C  VARIABLES :
C	FeedMoles(KL)   REAL * 8		LB - MOLES / HR OF MATERIAL
C	TK			  REAL * 8		TEMPERATURE, K
C	Enthalpy		  REAL * 8		IDEAL GAS ENTHALPY OF STREAM, BTU / HR
C	HFORM 		  REAL * 8		Heat of Formation OF STREAM, BTU / HR
C
C**********************************************************************
!DEC$ATTRIBUTES DLLEXPORT, alias:'IGENTHALPY'  ::IGEnthalpy
!DEC$ATTRIBUTES STDCALL, REFERENCE::IGEnthalpy
!DEC$ATTRIBUTES MIXED_STR_LEN_ARG::IGEnthalpy
C**********************************************************************
use CompData
IMPLICIT NONE

Real * 8  FeedMoles(KL), TK, Enthalpy,
&EMass, T2, T3, T4, T5, HFORM
Integer * 4  I

CALL IdealGasEnthalpy(FeedMoles, TK, Enthalpy)
RETURN

END
***/

// This subroutine calculates the enthalpy of the stream
// in units of BTU, non - ideal enthalpy corrections from RapidFlash.
int APIThermodynamics::CalcEnthalpy(double streamMoles[], double temperature, double PRatm, double Enthalpy, double VF)
{
	int RFlag = 0;
	double LIQ = 0.0;
	double VAP = 0.0;
	return 0; // SRKPTFlash(streamMoles, LIQ, VAP, temperature, PRatm, VF, Enthalpy, RFlag);
}

/***

	function ZPR(ICOMP, T, PATM)
	C ************************************************************************
	C  THIS FUNCTION CALCULATES THE COMPRESSIBILITY FACTOR, Z, AT TEMP = T(K)
	C   AND PRESSURE = P(ATM) FOR COMPONENT "ICOMP".
	C
	C  VARIABLES :
C
C       ICOMP    INTEGER * 4   THE LUMPED FEED COMPONENT NUMBER
C       T        REAL * 8      TEMPERATURE, DEG K
C       P        REAL * 8      PRESSURE, ATM
C
C  METHOD :
C
C       CALCULATE THE FACTOR Z USING FUNCTION ZP
C
C ************************************************************************
c  Expose FUNCTION
C**********************************************************************
!DEC$ATTRIBUTES DLLEXPORT, alias:'ZPR'         ::ZPR
!DEC$ATTRIBUTES STDCALL, REFERENCE::ZPR
!DEC$ATTRIBUTES MIXED_STR_LEN_ARG::ZPR
C**********************************************************************
REAL * 8  T, ZPR, PATM, ZP
INTEGER * 4 ICOMP

zpr = zp(icomp, t, patm)
return
end


FUNCTION ZP(ICOMP, T, PATM)
C ************************************************************************
C  THIS FUNCTION CALCULATES THE COMPRESSIBILITY FACTOR, Z, AT TEMP = T(K)
C   AND PRESSURE = P(ATM) FOR COMPONENT "ICOMP".
C
C  VARIABLES :
C
C       ICOMP    INTEGER * 4   THE LUMPED FEED COMPONENT NUMBER
C       T        REAL * 8      TEMPERATURE, DEG K
C       P        REAL * 8      PRESSURE, ATM
C
C  METHOD :
C
C       CALCULATE THE FACTOR Z USING REDUCED TEMPERATURE AND PRESSURE WITH
C       THE CORRELATION OF PITZER ET AL.THIS IS IN PERRY'S CHEMICAL ENGINEERS
C       HANDOOK 5TH EDITION, PAGE 4 - 61. THE PRINCIPLE OF CORRESPONDING STATES IS
C       AUGMENTED WITH A THIRD PARAMETER, THE ACENTRIC FACTOR.
C
C       THIS DOES NOT WORK FOR HYDROGEN SO WHEN ICOMP = 4 IT RETURNS 1.0.
C
C ************************************************************************

USE COMPDATA

IMPLICIT NONE


REAL * 8  T, ZP, PATM, TR, PR, B1, B0
INTEGER * 4 ICOMP

C   CHECK FOR HYDROGEN
IF(ICOMP.EQ.1) THEN
ZP = 1.0
RETURN
ENDIF

C   CALCULATE THE REDUCED PRESSURE AND TEMPERATURE
TR = T * TCINV(ICOMP)
PR = PATM * PCINV(ICOMP)
C
B0 = 0.1445 - (0.33 / Tr) - (0.1385 / (Tr*Tr)) -
&(0.0121 / (Tr*Tr*Tr))
B1 = 0.073 + (0.46 / Tr) - (0.5 / (Tr*Tr)) - (0.097 / (Tr*Tr*Tr)) -
&(0.0073 / (Tr*Tr*Tr*Tr*Tr*Tr*Tr*Tr))
ZP = 1.0 + (Pr / Tr) * (B0 + ACEN(ICOMP) * B1)

IF(ZP.LE.0.00001) ZP = 0.00001

END FUNCTION




FUNCTION ZSRKEXT(ICOMP, T, PATM)
C ************************************************************************
C  THIS FUNCTION CALCULATES THE COMPRESSIBILITY FACTOR, Z, AT TEMP = T(K)
C   AND PRESSURE = P(ATM) FOR COMPONENT "ICOMP".
C
C  VARIABLES :
C
C       ICOMP    INTEGER * 4   THE LUMPED FEED COMPONENT NUMBER
C       T        REAL * 8      TEMPERATURE, DEG K
C       P        REAL * 8      PRESSURE, ATM
C
C  METHOD :
C
C       CALCULATE THE FACTOR Z USING REDUCED TEMPERATURE AND PRESSURE WITH
C       THE SRK EQUATION OF STATE.SEE REID, PRAUSNITX & POLING,
C       'THE PROPERTIES OF GASES AND LIQUIDS, 4TH ED.', P43
C
C       THIS DOES NOT WORK FOR HYDROGEN SO WHEN ICOMP = 4 IT RETURNS 1.0.
C
C ************************************************************************
C
c  Expose FUNCTION
C**********************************************************************
!DEC$ATTRIBUTES DLLEXPORT, alias:'ZSRK1'       ::ZSRKEXT
!DEC$ATTRIBUTES STDCALL, REFERENCE::ZSRKEXT
!DEC$ATTRIBUTES MIXED_STR_LEN_ARG::ZSRKEXT
C**********************************************************************
EXTERNAL ZSRK

REAL * 8  T, ZSRK, PATM, ZSRKEXT
INTEGER * 4 ICOMP
ZSRKEXT = ZSRK(ICOMP, T, PATM)
END FUNCTION



FUNCTION ZSRK(ICOMP, T, PATM)
C ************************************************************************
C  THIS FUNCTION CALCULATES THE COMPRESSIBILITY FACTOR, Z, AT TEMP = T(K)
C   AND PRESSURE = P(ATM) FOR COMPONENT "ICOMP".
C
C  VARIABLES :
C
C       ICOMP    INTEGER * 4   THE LUMPED FEED COMPONENT NUMBER
C       T        REAL * 8      TEMPERATURE, DEG K
C       P        REAL * 8      PRESSURE, ATM
C
C  METHOD :
C
C       CALCULATE THE FACTOR Z USING REDUCED TEMPERATURE AND PRESSURE WITH
C       THE SRK EQUATION OF STATE.SEE REID, PRAUSNITX & POLING,
C       'THE PROPERTIES OF GASES AND LIQUIDS, 4TH ED.', P43
C
C       THIS DOES NOT WORK FOR HYDROGEN SO WHEN ICOMP = 4 IT RETURNS 1.0.
C
C ************************************************************************

USE COMPDATA

IMPLICIT NONE

EXTERNAL ZP

REAL * 8  T, ZSRK, PATM, TR, PR, Z(3), TCr, PCr, A, B, W
REAL * 8  AX, BX, CX, AA, BB, Q, R, PI, ZP
INTEGER * 4 ICOMP

PI = 3.141593

C   CHECK FOR HYDROGEN
IF(ICOMP.EQ.1) THEN
ZSRK = 1.0
RETURN
ENDIF

C   CALCULATE THE REDUCED PRESSURE AND TEMPERATURE
TR = T * TCINV(ICOMP)
PR = PATM * PCINV(ICOMP)
TCr = TC(ICOMP)
PCr = PC(ICOMP)
W = Acen(ICOMP)

A = 0.42748 * (RGAS*RGAS) * ((TCr*TCr) / PCr) *
&    (1.0 + (0.48 + 1.574*W - 0.176*W*W) * (1.0 - SQRT(Tr)))**2
B = 0.08664 * TCr * RGAS / PCr

C   SOLVE FOR Z AS ROOT OF CUBIC EQUATION(SEE REID, PRAUSNITX & POLING, P42)
C   Z**3 + AX * Z**2 + BX * Z + CX = 0
C	AX = -1
C     BX = AA - BB - BB * *2
C     CX = -AA * BB
C
AA = A * PATM / (RGAS*RGAS*T*T)
BB = B * PATM / (RGAS*T)
AX = -1.0
BX = AA - BB - BB * *2
CX = -AA * BB

Q = (AX*AX - 3.0*BX) / 9.0
R = (2.0*AX*AX*AX - 9.0*AX*BX + 27.0*CX) / 54.0

IF((Q*Q*Q - R * R).GE.0.0) THEN
C       THE EQUATION HAS THREE REAL ROOTS
W = DACOS(R / SQRT(Q**3))
Z(1) = -2.0 * SQRT(Q) * DCOS(W / 3.0) - AX / 3.0
Z(2) = -2.0 * SQRT(Q) * DCOS((W + 2.0*PI) / 3.0) - AX / 3.0
Z(3) = -2.0 * SQRT(Q) * DCOS((W - 2.0*PI) / 3.0) - AX / 3.0
C	  SELECT THE LARGEST VALUE
W = Z(1)
IF(Z(2).GT.W) W = Z(2)
IF(Z(3).GT.W) W = Z(3)
ZSRK = W
ELSE
C       THE EQUATION HAS ONE REAL ROOT
A = -SIGN(1.0D0, R)*(SQRT(R*R - Q * Q*Q) + DABS(R))**(1.0 / 3.0)
IF(A.NE.0.0) THEN
B = Q / A
ELSE
B = 0.0
ENDIF
Z(1) = (A + B) - AX / 3.0
ZSRK = Z(1)
ENDIF

C   CHECK IF <= 0, TRY PITZER CORRELATION
IF(ZSRK.LE.0.0) ZSRK = ZP(ICOMP, T, PATM)

END FUNCTION


FUNCTION PHISRK(IComp, TK, Patm)
C**********************************************************************
!DEC$ATTRIBUTES DLLEXPORT, alias:'PHISRK'      ::PHISRK
!DEC$ATTRIBUTES STDCALL, REFERENCE::PHISRK
!DEC$ATTRIBUTES MIXED_STR_LEN_ARG::PHISRK
C**********************************************************************
C  This function returns the fugacity coefficient for an input component
C(ICOMP) at an input temperature(TK) and pressure(PATM).The SRK EOS
C  is the basis for the calculation.
C
C  VARIABLES :
C
C       ICOMP    INTEGER * 4   THE LUMPED FEED COMPONENT NUMBER
C       T        REAL * 8      TEMPERATURE, DEG K
C       P        REAL * 8      PRESSURE, ATM
C
C  METHOD :
C
C       Calculate the SRK compressibility factor, Z, and then find the
C       fugacity coefficent.See Soave, Chem Eng Sci, 27, pp.1197 - 1203
C[1972]
C
C       THIS DOES NOT WORK FOR HYDROGEN SO WHEN ICOMP = 4 IT RETURNS 1.0.
C
C**********************************************************************

USE COMPDATA


IMPLICIT NONE
EXTERNAL ZSRK
C     ******************************************************************
INTEGER * 4 ICOMP
REAL(KIND = 8) A, B, Z, PHISRK
REAL(KIND = 8) TR, PR, ALPHA, W, ZSRK, TK, PATM
C     ******************************************************************
C   CHECK FOR HYDROGEN
IF(ICOMP.EQ.1) THEN
PHISRK = 1.0
RETURN
ENDIF
C   INITIALIZE
TR = TK * TCINV(ICOMP)
PR = PATM * PCINV(ICOMP)
W = Acen(ICOMP)
Z = ZSRK(ICOMP, TK, PATM)
C   CALCULATE ALPHA
ALPHA = (1.0 + (0.48 + 1.574*W - 0.176*W*W)*(1.0 - DSQRT(TR)))**2
C   CALCULATE A AND B
A = 0.42747 * ALPHA * PR / (TR*TR)
B = 0.08664 * PR / TR
PHISRK = EXP(Z - 1.0 - LOG(Z - B) - A * LOG((Z + B) / Z) / B)

RETURN
END

***/
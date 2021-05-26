import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


## Read in different data sets of interest:
data_temp = pd.read_csv('./clean_inputs_for_greenroof_tests/AirTemp_clean_Fond_Du_lac_May_to_October_2019.csv',parse_dates=True,index_col=0)
data_longwave = pd.read_csv('./clean_inputs_for_greenroof_tests/Merra_2_longwave_Fond_du_lac_May_to_October_2019.csv',parse_dates=True,index_col=0)
data_pressure = pd.read_csv('./clean_inputs_for_greenroof_tests/Pressure_clean_mb_Fond_Du_lac_May_to_October_2019.csv',parse_dates=True,index_col=0)
data_RH = pd.read_csv('./clean_inputs_for_greenroof_tests/RH_clean_Fond_Du_lac_May_to_October_2019.csv',parse_dates=True,index_col=0)
data_shortwave = pd.read_csv('./clean_inputs_for_greenroof_tests/Shortwave_clean_Fond_Du_lac_May_to_October_2019.csv',parse_dates=True,index_col=0)
data_WindSpeed = pd.read_csv('./clean_inputs_for_greenroof_tests/WindSpeed_clean_Fond_Du_lac_May_to_October_2019.csv',parse_dates=True,index_col=0)

data_precip = pd.read_csv('./clean_inputs_for_greenroof_tests/Precip_clean_Fond_Du_lac_May_to_October_2019.csv',parse_dates=True,index_col=0)

print(data_temp)

## create values to catch values
SMC_total = np.zeros(52849)
CMC_total = np.zeros(52849)
DETENTION_STORAGE_IN_total = np.zeros(52849)
RUNOFF_total = np.zeros(52849)

TGRP_total = np.zeros(52849)
TGRL_total = np.zeros(52849)
TBOT_total = np.zeros(52849)

H_total = np.zeros(52849)
LH_total = np.zeros(52849)
G_total = np.zeros(52849)
NETRAD_total = np.zeros(52849)

SMC_total[0] = 0.4
DETENTION_STORAGE_IN_total[0] = 0
CMC_total[0] = 0.0
RUNOFF_total[0] = 0.0

TGRL_total[0] = 292 # kelvin
TGRP_total[0] = 292
TBOT_total[0] = 292

iter_total = 0

#####
# define constants
global ZA  # level that atmospheric observations were taken
global Z0R  # roughness length of the roof
global AKANDA_URBAN  # coefficient modifying the kanda approach to computin surface layer coefficients (1.29)
global CP  # heat cpacity of dry air [cgs units]
global CPP  # heat capacity in J/K/kg

# Pulled from Green Roof Paper
global DT #time step of interest in seconds
global FVEG #fractional vegitation
global SAI_LAI_storage_factor #from NOAHMP TABLES. How much per LAI can be stored in a plant
global TOTAL_VEG #the amount of vegitation is present (LAI and SAI together)
global STORAGE_MAX #the amount of storage allowed in retention 1 (m)
global a #parameter to help switch on and off the water transfer from retention 1 to 2
global DETENTION_storage_MAX #amount of storage in the retention 2
global k #parameter to go ahead of the head for routing
global n #nonlinear routing parameter
global SMAX


## Values that are needed to run this model with description

## These are needed for GLOBAL 
FVEG = 0.8 #fraction of the grid square that is covered by vegetation (unitless) (literture suggests 0.6 to 1.0)

DT = 300.0 #time step  in seconds 
SAI_LAI_storage_factor = 0.1 # a factor that transforms the LAI and SAI values into a depth (unitless)
TOTAL_VEG = 1.5 # the total LAI and SAI values together (m^2 per m^2) Literature suggests 1 - 3
STORAGE_MAX = 76.2/1000 # max storage of the growing medium in dpeth (m) Based on 6 inch depth 
## in addition, storage max is used to back calculate the depth of sil in the green roof model
a = 0.75 # ratio to determine *when raining* water will go to detention and runoff. Bascially above this value, water 
# begins to drain
DETENTION_storage_MAX = 10.0/1000 # maximum amount of detention storage (mm)
k = (0.002/60)*DT # the green roof model uses a non-linear formula. K is the multiplier before the value of Q
n = 1.67 # non-linear exponent 

for i in range(0,52848):
	PS = data_pressure['BP'].values[i] # Hpa
	TA = data_temp[' Value'].values[i] + 273.15 #degrees Kelvin
	SSG = data_shortwave[' Value'].values[i] # watts per meter squared
	LLG = data_longwave['Value'].values[i] # watts per meter squared


	Sat_press = 6.1078*10**(7.5*(TA-273.15)/((TA-273.15)+237.3)) #uses the tetens equation (hPa)
	Pressure_vapor = Sat_press*data_RH[' Value'].values[i]/100 #vapor pressure is in (hPa)
	RHOO = ((PS - Pressure_vapor)*100)/(287.058*TA) + Pressure_vapor*100/(461.495*TA) #air density in kg/m^3
	RHO = RHOO * 0.001 # air density in cgs units

	w = Pressure_vapor*287.058/(461.495*(PS - Pressure_vapor))
	QA = w/(w+1) #specific humidity in kg/kg


	UA = data_WindSpeed['Value'].values[i] # wind speed in m/s
	TH2V = TA*(1000/PS)**.286 #potential temperature in K
	
	TGRP = TGRP_total[i] #surface temperature in K

	RAIN = (data_precip['Value'].values[i]/(DT))*(3600) # conver to MM/hr

	RUNOFF = RUNOFF_total[0] #mm

	## Atmospheric Values needed to run the model 
	# TGRP = 301 # Green Roof surface temperature (Kelvin)
	# QA = 0.010  # Specific Humidity of the atmosphere (kg/kg)
	ZA = 10  # Level at which the atmospheric values are being measured (meters)
	Z0R = 0.25  # Surface Roughness (meters)
	# TA = 295.73517 # Atmospheric Air Temperature (Kelvin)
	# TH2V = 22.690049+273.15  # Atmospheric Potential Air Temperature (Kelvin)
	# UA = 0.842624   # Atmospheric Wind Speed (m/s)
	# PS = 983.9154 # Surface Pressure hPa
	# RAIN = 0.0 # rainfall depth in current time step (mm/hr)
	# RHOO = 1.225 #air density kg m^3
	# RHO = RHOO * 0.001  # air density cgs units


	## radiation information
	# SSG = 674.7253 # watts per meter squared
	# LLG = 379.1754 # watts per meter squared
	ALBEDO = 0.07 #shortwave Based on Grassland Values in NOAH-MP Table
	EPSV = 0.93 # longwave emissivity (current value, hard coded in, but is consistent with cambell an norton)


	SX = (SSG)/697.7/60 # DOWNWARD SHORTWAVE RADIATION Now in Langlys per minute
	RX = LLG/697.7/60 #downward longwave radiation in langlys per minute
	SGR = SX*(1-ALBEDO) # net amount of shortwave (second value is currently the bulk albedo)
	SIG = 8.17*10**-11 # the stefan boltzman constant in calories

	##water balance and soil values
	SMC = SMC_total[i] # soil moisture in m^3 per m^3
	CMC = CMC_total[i] # canopy storage depth at current time step
	DETENTION_STORAGE_IN = DETENTION_STORAGE_IN_total[i] # amount of water that is being stored in outflow area of green roof
	ET = 0.0 # evapotranspiration values
	TBOT = TBOT_total[i] # bottom of the soil temperature (Kelvin)
	TGRL = TGRL_total[i] # soil temperature at central node (kelvin)



	## These parameters all controll the 'soil'
	SMCMAX = 0.5 # M3 per M3 #needs to be set
	SMCDRY = 0.1 # M3 per M3 # needs to be set
	SMCREF = 0.25# M3 per M3 # needs to be set

	DF = 1 # this is the diffusivity of the soil, and will be overwritten (just initilized)
	QUARTZ = 0.4 # not very sensitive to our solutions, so keeping 0.4
	ZSOILR = STORAGE_MAX/(SMCMAX) #based on a 6 inch deep green roof (extensive)
	DZGR = STORAGE_MAX/SMCMAX #based on a 6 inch deep green roof (extensive)
	CSOIL = 1.125*10**6 #Specific heat (1500 J/kg K ) * Density (750 kg/m^3)


	## Transfer coefficients from M-O Theory(Unitless) (just initilization)
	CMGR_URB = 0.1 # transfer coefficient for momentum. Will be overwritten
	CHGR_URB = 0.1 # transfer coefficient for heat and water, will be overwritten
	RLMO_URB = 1 # stability paramters (1 / Obukov length)
	CDGR = 1.0 # This is a drag coffeicient for the surface layer

	## paramters that are specified in the header files of fortran and locations of the 
	## URBAN.TBL. Should not need to be touched
	AKANDA_URBAN = 1.29  # unitles parameter
	CP = 0.24 # specific heat capacity of dry air cgs units
	CPP = 1004.5 # specific heat capacity of dry air [J/K/kg]
	SHDFAC = FVEG # Veg Fraction
	FXEXP = 2.0 # unitless exponent that is used to control the ET
	CFACTR = 0.5 # value used to translate between LAI 
	CMCMAX = TOTAL_VEG*CFACTR/1000 #from Precip_heat calculation of the max amount (mm)
	LAI = TOTAL_VEG # controled by the total veg parameter for continutint
	RSMIN = 100 # minimum resistacne of the vegitation [s m^-1]
	RSMAX = 5000 # maximim resistance of teh vegitation [s m^-1]
	RGL = 100 # paramter used in radiation stress function 
	HS = 36 # paramter used in the vapor pressure deficit function
	EL = 583 # latent heat of vaporization in CGS [calories per gram]
	ELL = 2.442*10**6 # latent heat of vaporization [ joules per kg ]




	######
	# define functions
	def SFCDIF_URB(z, z0, ground_theta, air_theta, wind_speed, akanda, akms, akhs, RLMO, cdgr):
		# based on the urban model, we are going to use the Paulson surface functions

		zu = z0  # local roughness
		RDZ = 1/z  # inverse of where we made the measurements

		CXCH = 0.001*RDZ  # the 0.001 2as hardcoded it. I believe this is the ratio of the heights (on surface and atmosphere)
		DTHV = air_theta - ground_theta  # delta theta for temperatures

		# beljars correction of USTAR

		DU2 = max(wind_speed**2, 1.000*10**-4)  # minimum values
		# add if statements to avoid the tangent linear problems near zero
		BTGH = (1 / 270) * 9.81 * 1000  # based on paramters in the beginning of the function. This seems to be hard coded, and the 1000
		# is particulalry worrysome (assumed PBL height). Lookeing at Stuhl et al. should be okay (PBLH scale ~1000)

		if (BTGH*akhs*DTHV) != 0:
			WSTAR2 = 1.2**2 * np.abs(BTGH*akhs*DTHV)**(2/3)
		else:
			WSTAR2 = 0.0


		#now apply the Zilitinkevitch Approach for ZT
		USTAR = max(np.sqrt(akms*np.sqrt(DU2+WSTAR2)),0.07)


		ZT = np.exp(2.0-akanda*(258.2**2 * USTAR*z0)**0.25)*z0 #thermally varrying roughness length
		ZLSU = z + zu #zeta urban

		ZLST = z + ZT #zeta thermal

		RLOGU = np.log(ZLSU/zu) #log of the zeta

		RLOGT = np.log(ZLST/ZT) #also log of the zeta

		RLMO = (0.4*1/270 * 9.8)*akhs*DTHV/(USTAR**3) #1/L (one over the obuhov length)

		for i in range(0,6):

			ZETALT = max(ZLST*RLMO,-5)

			RLMO = ZETALT/ZLST
			ZETALU = ZLSU*RLMO
			ZETAU = zu*RLMO
			ZETAT = ZT*RLMO

			if RLMO < 0.0:
				XLU4 = 1.0 - 16*ZETALU
				XLT4 = 1.0 - 16*ZETALT
				XU4 = 1.0 - 16*ZETAU
				XT4 = 1.0 -16*ZETAT

				XLU = np.sqrt(np.sqrt(XLU4))
				XLT = np.sqrt(np.sqrt(XLT4))
				XU = np.sqrt(np.sqrt(XU4))
				XT = np.sqrt(np.sqrt(XT4))

				PSMZ = PSPMU(XU)
				SIMM = PSPMU(XLU) - PSMZ + RLOGU
				PSHZ = PSPHU(XT)
				SIMH = PSPHU(XLT) - PSHZ + RLOGT
			else:
				ZETALU = min(ZETALU,1.0)
				ZETALT = min(ZETALT,1.0)

				PSMZ = PSPMS(ZETAU)
				SIMM = PSPMS(ZETALU) - PSMZ + RLOGU
				PSHZ = PSPHS(ZETAT)
				SIMH = PSPHS(ZETALT) - PSHZ + RLOGT

			USTAR = max(np.sqrt(akms*np.sqrt(DU2 + WSTAR2)),1.0*10**-4)
			ZT = np.exp(2.0 - akanda*(258.2**2 * USTAR * z0)**0.25)*z0

			ZLST = z + ZT #zeta thermal

			RLOGT = np.log(ZLST/ZT) #also log of the zeta

			USTARK = USTAR*0.4

			akms = max(USTARK/SIMM,CXCH)
			akhs = max(USTARK/SIMH,CXCH)

			if (BTGH*akhs*DTHV) != 0:
				WSTAR2 = 1.2**2 * np.abs(BTGH*akhs*DTHV)**(2/3)
			else:
				WSTAR2 = 0.0

			RLMN  = (0.4*1/270 * 9.8)*akhs*DTHV/(USTAR**3) #1/L (one over the obuhov length)
			RLMA = RLMO * .15 + RLMN * (1 - .15)

			RLMO = RLMA

		CD = USTAR*USTAR/wind_speed**2

		return(akhs, akms, RLMO, CD)

	def PSPHS(ZZ):
		#pauslons surface exchange functions based on SFCDIF_URB in NOAHMP
		#Heat Stable
		return 5.0*ZZ

	def PSPMS(ZZ):
		#pauslons surface exchange functions based on SFCDIF_URB in NOAHMP
		#momentum Stable
		return 5.0 *ZZ

	def PSPMU(XX):
		#pauslons surface exchange functions based on SFCDIF_URB in NOAHMP
		#Momentum UnStable

		return -2.0* np.log( (XX + 1)*0.5) - np.log( (XX**2 + 1)*0.5) + 2.0*np.arctan(XX) - 3.14159265/2.

	def PSPHU(XX):
		#pauslons surface exchange functions based on SFCDIF_URB in NOAHMP
		#Momentum UnStable
		return -2.0 * np.log( (XX**2 +1) * 0.5)

	def DIREVAP(ETP,SMC,SHDFAC,SMCMAX,SMCDRY,FXEXP):
		#based on the direct soil evaporation in the WRF URban model
		#one underying assumption with how we are doing this at the moment (green roof implementation)
		#is that the entire moisture column is going to be partcipitating in the direct evaporation, intead
		#of the top soil layer (~0.05 m in WRF)

		SRATIO = (SMC - SMCDRY)/(SMCMAX - SMCDRY)

		if SRATIO > 0:
			# FX > 1 represents the demand control
			# FX < 1 represents flux control
			FX = SRATIO**FXEXP
			FX = max(min(FX,1),0)
		else:
			FX = 0

		EDIR = FX*(1-SHDFAC)*ETP*0.001

		return EDIR

	def TRANSP(ET,SHDFAC,ETP1,CMC,CFACTR,CMCMAX,LAI,RSMIN,RSMAX,RGL,SX,TS,TA,QA,SMC,SMCWLT,SMCREF,CPP,PS,CH,EPSV,DELT,HS):
		#Based on the transpiration function that is in WRF Urban

		SLV = 2.501*10**6
		SIGMA = 5.67*10**-8
		ETT = 0.0
		ET = 0

		#initialize the canopy resistance multipliers:
		RCS = 0.0
		RCT = 0.0
		RCQ = 0.0
		RCSOIL = 0.0

		#Contripution due to the incoming solar radiation
		#canopy resistance due to solar

		FF = 0.55*2*SX*697.7*(60/(RGL*LAI))
		RCS = (FF + RSMIN/RSMAX) / (1.0+FF)
		RCS = max(RCS,0.0001)

		#Contribution due to the air temperature
		#RCT from Noilhan and Planton (1989, MWR)

		RCT = 1.0 - 0.0016*( (298 - TA)**2.0)
		RCT = max(RCT,0.0001)

		#contribution due to the vapor pressure deficit at the first model level
		#RCQ expression from SSIB (Niyogi and Raman, 1997)

		EA = 6.11*np.exp( (2.5*10**6 /461.51) * (TA - 273.15)/(273.15*TA)) #Full Clausies Claperon Equation
		WS = 0.622*EA/1013 #saturation vapor pressure
		RCQ = 1.0/(1.0 + HS*(WS - QA))
		RCQ = max(RCQ,0.01)

		#contribution due to soil moisture avaiability (only need to do a single value because we have a SINGLE)
		#value

		GX = (SMC - SMCWLT) / (SMCREF - SMCWLT)
		if GX > 1:
			GX = 1
		elif GX < 0:
			GX = 0

		RCSOIL =  GX
		RCSOIL = max(RCSOIL,0.0001)

		RC = RSMIN/ (LAI *RCS * RCT * RCQ * RCSOIL) #the total resistance (invert to get a conductance!)

		DESDT = 0.622*SLV*EA/(461.51*TA*TA*1013)
		DELTA = (SLV/CPP)*DESDT
		RR = (4*EPSV*SIGMA*287.04/CPP) * TA**4/(TS*CH) + 1 #so this is from an old PBL scheme. It is used to relate a 'Plant Coefficient' to the resistance of the soil

		PC = (RR + DELTA)/(RR * (1 + RC * CH) + DELTA) #plant coefficient based on the total resistance that was calculated before
		#thie Plant coefficient is from an old scheme from OSU!

		if CMC != 0.0:
			ETT1 = SHDFAC * PC * ETP1 * (1.0 - (CMC/CMCMAX)**CFACTR) * 0.001
		else:
			ETT1 = SHDFAC*PC*ETP1*0.001

		#evapotransipriaton total
		ETT += ETT1

		if CMC > 0.0:
			EC = SHDFAC* ((CMC/CMCMAX)**CFACTR) * ETP1 * 0.001
		else:
			EC = 0.0

		CMC2MS = CMC/DELT
		EC = min(CMC2MS,EC)

		return(EC,ETT)

	def TDFCND(DF,SMC,QZ,SMCMAX):
		#Calculate the Thermal COnductivity of the Green Roof Growing Medium
		#Based on the Peters-Lidard Approach - Peters-Lidard et al. 1998

		#Saturation Ratio
		SATRATIO = SMC/SMCMAX
		#Thermal Conductivity of Water
		THKW = 0.57
		#"Other Soil Components"
		THKO = 2.0
		#Quartz Conductivity
		THKQTZ = 7.7
		#Solids conductivity
		THKS = (THKQTZ**QZ)*(THKO **(1 - QZ))

		#Saturated Thermal Conductivity
		THKSAT = THKS**(1-SMCMAX) * THKW**(SMCMAX)

		#DRY THERMAL CONDUCTIVITY
		GAMMD = (1 - SMCMAX)*2700

		THKDRY = (0.135*GAMMD + 64.7)/(2700 - .947*GAMMD)

		#Kersten Number for Fine soils, wiht at least 5% of particles having
		#Diameters less than 2.0 * 10 ** -6

		if SATRATIO > 0.1:
			AKE = np.log10(SATRATIO) + 1.0
		else:
			AKE = 0.0

		DF = AKE * (THKSAT - THKDRY) + THKDRY
		return DF

	def SHFLX(TGRL, SMCMAX, SMC, TGRP, DF, STORAGE_MAX, TBOT, YY, ZZ1, CSOIL):

		# This is an adjusted version of the heat flow equation that is used in
		# in the WRF Urban Scheme. The formulation solve the heat diffusion Equation
		# using a backward Euler scheme and a centered difference
		# will need to use the CPAR if we use the neuman boundary conditions

		# First Calculate the Heat Capacity
		CH2O = 4.2*10**6
		CAIR = 1004.0

		HCPCT = SMC*CH2O + (1.0 - SMCMAX)*CSOIL + (SMCMAX - SMC)*CAIR

		# create soil heat beta (AKA Diffusivity) as thermal conductivity over the
		# volumetric heat capacity
		SOIL_HEAT_BETA = DF/HCPCT

		# Calculate the DX, based on the STORAGE_MAX and the SMCMAX
		DX = (STORAGE_MAX/SMCMAX)/2
		# now creat the s coeffcient (s = beta * DT/(DX**2)) from the
		# numerical differentiation scheme

		A = np.zeros((3,3))
		B_vals = np.zeros((3,1))

		S = SOIL_HEAT_BETA*DT/(DX**2)

		A[0,0] = 1
		A[1,0] = -S
		A[1,1] = 1+2*S
		A[1,2] = -S
		A[2,1] = -2*S
		A[2,2] = 1+2*S

		B_vals[0,0] = TGRP
		B_vals[1,0] = TGRL
		B_vals[2,0] = TBOT

		# print(A)
		# print(np.linalg.inv(A))
		# print(np.dot(np.linalg.inv(A),B_vals))
		

		# Now solve for the New temperature:
		TGRP,TGRL,TBOT = np.dot(np.linalg.inv(A),B_vals)
		# print(TBOT)
		# ss
			#(TGRL +S*TBOT + S*TGRP)/(1+2*S)

		# we now adjust the Surface temperature ?

		TGRP = (YY + (ZZ1 - 1)*TGRL)/ZZ1

		return(TGRP,TGRL,TBOT)

	def green_roof_soil(CANLIQ, PRECIP_input, MAXLIQ, EC_IN, ET_IN, EDIR_IN,SMC,SMCMAX,SMCDRY,DETENTION_STORAGE):

		# Used to call the interception_throughfall
		# The retention_1
		# The retention_2 functions
		# The inputs needed for this program are the canopy liquid, the incoming
		# precip, the maximum amount of liquid able to be held by the canopy, and
		# the amount of canopy ET, the amount of 'storage' in, which is the a mapping
		# between the soil moisture and a reservoir capacity, the amount of direct
		# ET and the amount transpired, and the amount of storage in the dentention
		# package. It is important to note that the water that enters into the retention_2
		# areas is thought to be 'runoff', and will be added to the apropriate runof
		# terms.

		# This is based on the locatelli et al 2014 paper for modeling green roofs
		# Key parameters that should be tuned:
		# global DT #time step of interest in seconds
		# global FVEG #fractional vegitation
		# global SAI_LAI_storage_factor #from NOAHMP TABLES. How much per LAI can be stored in a plant
		# global TOTAL_VEG #the amount of vegitation is present (LAI and SAI together)
		# global STORAGE_MAX #the amount of storage allowed in retention 1 (m)
		# global a #parameter to help switch on and off the water transfer from retention 1 to 2
		# global DETENTION_storage_MAX #amount of storage in the retention 2 (m)
		# global k #parameter to go ahead of the head for routing per time step
		# global n #nonlinear routing parameter per time step
		# global SMAX

		# Adust precip input (Input is in mm/h, but the ET and EC are in m/s)
		PRECIP_input = (PRECIP_input/3600)/1000

		# First Linearly Map from soil moisture to STORAGE
		STORAGE = ((STORAGE_MAX)/(SMCMAX - SMCDRY))*(SMC-SMCDRY) #units of STORAGE_MAX (m)


		# Calculate the amount that is intercepted and the amount of water
		# that is throughfall:

		Precip_dripped_to_ground, Precip_throughfall_to_ground, CANLIQ = interception_throughfall(CANLIQ, PRECIP_input, MAXLIQ, EC_IN)
		STORAGE, out_box2_to_box3 = retention_1(Precip_dripped_to_ground+Precip_throughfall_to_ground, STORAGE, ET_IN, EDIR_IN)
		DETENTION_STORAGE, runoff_final = retention_2(out_box2_to_box3,DETENTION_STORAGE)
	

		# print(DETENTION_STORAGE)
		# Now we map from the storage back to soil Moisture
		SMC = ((SMCMAX - SMCDRY)/(STORAGE_MAX))*STORAGE + SMCDRY #volumetric soil Moisture


		return(SMC, CANLIQ, DETENTION_STORAGE, runoff_final)

	def interception_throughfall(CANLIQ,PRECIP,MAXLIQ,EC):
		QINTER = FVEG*PRECIP #line 1231 in precip heat
		QINTER = min(QINTER, (MAXLIQ - CANLIQ)/DT * (1-np.exp(-PRECIP*DT/MAXLIQ)))
		QINTER = np.max(QINTER,0)

		QDRIPR = FVEG*PRECIP - QINTER
		QTHROR = (1-FVEG)*PRECIP
		CANLIQ = max(0.0,CANLIQ+QINTER*DT - EC*DT)

		return(QDRIPR,QTHROR,CANLIQ)

	def retention_1(PRECIP_IN,STORAGE,ET,EDIR):

		#based on a single simple reservoir. This is pulled directly from equation 4 from Locatelli et al 2014

		if STORAGE/STORAGE_MAX > a:
			if PRECIP_IN > 0:

				outflow = (PRECIP_IN) * ((STORAGE/STORAGE_MAX) - a)/(1-a)

			else:
				outflow = 0.0


		else:
			outflow = 0



		storage_new = max(0.0, STORAGE + PRECIP_IN*DT - outflow*DT - ET*DT - EDIR*DT)

		return(storage_new,outflow)

	def retention_2(INFLOW,STORAGE_DETENTION):

		runoff= k*((STORAGE_DETENTION)**n)

		if runoff > STORAGE_DETENTION:
			runoff = STORAGE_DETENTION

		STORAGE_DETENTION_new = INFLOW*DT + STORAGE_DETENTION - runoff*DT

		if STORAGE_DETENTION_new > DETENTION_storage_MAX:
			runoff += (STORAGE_DETENTION_new - DETENTION_storage_MAX)/DT
			STORAGE_DETENTION_new = DETENTION_storage_MAX
		
		if STORAGE_DETENTION_new < 0.0:
			STORAGE_DETENTION_new = 0


		return(STORAGE_DETENTION_new, runoff)


	## Main CoDE!!


	T1VGR = TGRP*(1.0+.61*QA) #convert to virtual temperature of green roof layer

	CMGR_URB, CHGR_URB, RLMO_URB, CDGR = SFCDIF_URB(ZA,Z0R,T1VGR,TH2V,UA,AKANDA_URBAN,CMGR_URB,CHGR_URB,RLMO_URB,CDGR) #probably works?
	ALPHAGR = RHO * CP * CHGR_URB
	CHGR = ALPHAGR/RHO/CP/UA


	# we are now okay to move on

	for it in range(0, 200):
		ES = 6.11*np.exp( (2.5*10**6 /461.51) * (TGRP - 273.15)/(273.15*TGRP)) #Full Clausies Claperon Equation
		DESDT = (2.5*10**6 / 461.51)*ES/(TGRP**2)


		QS0GR = 0.622*ES/(PS -0.378*ES) # Slope dependent on pressure

		DQS0GRDTGR = DESDT*0.622*PS/((PS - 0.378*ES)**2)

		EPGR = RHOO*CHGR*UA*(QS0GR-QA) # potential evaporation from the atmosphere side KG/m2/s

		if EPGR >= 0.0:

			#direct evaporation
			EDIR = DIREVAP(EPGR,SMC,SHDFAC,SMCMAX,SMCDRY,FXEXP)
			#evapotranspiration
			ECR, ETTR = TRANSP(ET,SHDFAC,EPGR,CMC,CFACTR,CMCMAX,LAI,RSMIN,RSMAX,RGL,SX,TGRP,TA,QA,SMC,SMCDRY,SMCREF,CPP,PS,CHGR,EPSV,DT,HS)
			#update soil layers
			SMC_new, CMC_new, DETENTION_STORAGE_IN_new, RUNOFF_new = green_roof_soil(CMC, RAIN, CMCMAX, ECR, ETTR, EDIR, SMC, SMCMAX, SMCDRY, DETENTION_STORAGE_IN)

		else:
	
			DEW = -EPGR
			RAINDR = RAIN+DEW*3600
			EDIR = 0.0
			ECR = 0.0
			ETTR = 0.0

			SMC_new, CMC_new, DETENTION_STORAGE_IN_new, RUNOFF_new = green_roof_soil(CMC, RAINDR, CMCMAX, ECR, ETTR, EDIR, SMC, SMCMAX, SMCDRY, DETENTION_STORAGE_IN)

		# print(RUNOFF_new)
		#Convert Modeled ET from m to s to kg m^2 per second
		# if np.isnan(DETENTION_STORAGE_IN):
		# 	print(RUNOFF)
		# 	break
		EDIR = EDIR * 1000.0
		ETTR = ETTR * 1000.0
		ECR = ECR * 1000.0
		ETAR = EDIR + ETTR + ECR

		if ETAR < 1*10**-20:
			BETGR = 0.0
		else:
			BETGR = ETAR/EPGR

		ELEGR = ETAR*RHO*EL/RHOO * 100

		DF = TDFCND(DF,SMC,QUARTZ,SMCMAX)

		DF1 = DF*np.exp(-2.0*SHDFAC)

		#Longwave?
		RGR = EPSV*(RX - (SIG*(TGRP**4)/60))
		RGRR = (RGR+SGR) * 697.7 * 60
		RCH = RHOO*CPP*CHGR
		RR1 = EPSV*(TA**4) * 6.48*10**-8/(PS*CHGR) + 1.0

		if RAIN > 0.0:
			RR2 = RR1 + RAIN/3600 * 4.218*10**3 / RCH
		else:
			RR2 = RR1

		YY = TA + (RGRR/RCH - BETGR *EPGR * ELL/RCH)/RR2
		ZZ1 = DF1/(-0.5 * ZSOILR * RCH * RR2) + 1.0

		#calculate the Sensible heat flux
		HGR = RHO * CP * CHGR * UA * (TGRP - TA) * 100
		RUNOFF = RUNOFF/DT

		#calculate the ground heat flux
		G0GR = DF1 * (TGRP - TGRL)/(DZGR/2)/697.7/60

		#Residual of energy
		FV = SGR + RGR - HGR - ELEGR - G0GR


		#calculating the slope paramters of each term to track of we are in a stable solution
		DRRDTGR = (-4.0 * EPSV*SIG*TGRP**3)/60
		DHRDTGR = RHO*CP*CHGR*UA*100
		DELERDTGR = RHO*EL*CHGR*UA*BETGR*DQS0GRDTGR*100

		DG0RDTGR = 2*(DF1/DZGR) * (1.0/4.1868)*1*10**-4

		DFDVT = DRRDTGR - DHRDTGR - DELERDTGR - DG0RDTGR
		DTGR = FV/DFDVT/6

		TGR = TGRP - DTGR
		TGRP = TGR

		if (np.abs(FV) < 0.0001 and np.abs(DTGR) < 0.0001):
			SMC = SMC_new
			DETENTION_STORAGE_IN = DETENTION_STORAGE_IN_new
			CMC = CMC_new
			RUNOFF = RUNOFF_new
			break

	# print(TGRP-273.15)	
	TGRP_, TGRL ,TBOT= SHFLX(TGRL, SMCMAX, SMC, TGRP, DF1, STORAGE_MAX, TBOT, YY, ZZ1, CSOIL)
	# print(TGRP)


	FLXTHGR = HGR/RHO/CP/100
	FLXHUMGR = ELEGR/RHO/EL/100
	# print('WOOT')
	# print(FLXTHGR*CPP*RHOO) #watts per meter squared
	# print(FLXHUMGR*RHOO*ELL)
	# print((RGR+SGR)*697.7*60)
	# print(G0GR*697.7*60)

	H_total[i] = FLXTHGR*CPP*RHOO
	LH_total[i] = FLXHUMGR*RHOO*ELL
	G_total[i] = G0GR*697.7*60
	NETRAD_total[i] = (RGR+SGR)*697.7*60

	SMC_total[i+1] = SMC
	TGRP_total[i+1] = TGRP
	TGRL_total[i+1] = TGRL
	TBOT_total[i+1] = TBOT


	DETENTION_STORAGE_IN_total[i+1] = DETENTION_STORAGE_IN
	CMC_total[i+1] = CMC
	RUNOFF_total[i+1] = RUNOFF

	iter_total += it
print(iter_total)

final_dataframe = pd.DataFrame({'H':H_total,'LH':LH_total,'G':G_total,'RN':NETRAD_total,'SMC':SMC_total,'DETENTION_STORAGE':DETENTION_STORAGE_IN_total,'CMC':CMC_total,'RUNOFF':RUNOFF_total,'TGRP':TGRP_total,'TGRL':TGRL_total,'TBOT':TBOT_total})
final_dataframe.to_csv('Fond_du_lac_2019_growing_season.csv')
# plt.plot(SMC_total)
plt.figure(figsize=(15,9))
# plt.subplot(311)
# plt.xlim([dates_subset[0],dates_subset[-1]])
# plt.bar(dates_subset, data_subset.P.values,width=0.01,linewidth=0.001)
# plt.ylabel(r'Rainfall $[mm]$',fontsize=18)
# plt.legend(['Rainfall'],frameon=False,fontsize=18)
# plt.xticks([])
# plt.yticks(fontsize=14)

# # plt.subplot(312)
# # plt.plot(dates_subset, data_subset.TA_1_1_1.values,color='#BD4D24',linewidth=2.5)
# # plt.xlim([dates_subset[0],dates_subset[-1]])
# # plt.ylabel(r'Temperature $[C]$',fontsize=16)
# # plt.xticks([])
# # # plt.xticks([])

plt.subplot(311)
plt.plot(SMC_total,color='k',label = 'Soil Moisture',linewidth=2.5)
plt.ylabel(r'SWC $[m^3 m^{-3}]$',fontsize=18)
# plt.xticks([])
plt.yticks(fontsize=14)
plt.legend(fontsize=18,frameon=False)


plt.subplot(312)
plt.plot(LH_total)
plt.show()

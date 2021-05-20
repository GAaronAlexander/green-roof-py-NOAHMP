import numpy as np
import matplotlib.pyplot as plt



#this is a toy model that based off of the WRF simulations
#it contains three 'boxes' that are throught to be a) the vegitation canopy
#b) the storage of the soil water and C) the detention portion of the green roof



global DT #time step of interest in seconds
global FVEG #fractional vegitation
global SAI_LAI_storage_factor #from NOAHMP TABLES. How much per LAI can be stored in a plant
global TOTAL_VEG #the amount of vegitation is present (LAI and SAI together)
global STORAGE_MAX #the amount of storage allowed in retention 1 (mm)
global a #parameter to help switch on and off the water transfer from retention 1 to 2
global DETENTION_storage_MAX #amount of storage in the retention 2
global k #parameter to go ahead of the head for routing
global n #nonlinear routing parameter
global SMAX

FVEG = 0.8
DT = 60.0
SAI_LAI_storage_factor = 0.1
TOTAL_VEG = 5.0
STORAGE_MAX = 10.0/1000
a = 0.65
DETENTION_storage_MAX = 2.0/1000
k = 0.0038
k1=  0.215/60  #per time step
n = 1.67
alpha = 1



#this is based on the PRECIP_HEAT information
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


	STORAGE_DETENTION_new = INFLOW*DT + STORAGE_DETENTION - runoff*DT

	if STORAGE_DETENTION_new > DETENTION_storage_MAX:
		runoff += (STORAGE_DETENTION_new - DETENTION_storage_MAX)/DT
		STORAGE_DETENTION_new = DETENTION_storage_MAX

	return(STORAGE_DETENTION_new,runoff*DT)

#This is based on the green roof system shown in the


MAX_amount_intercepted = TOTAL_VEG*SAI_LAI_storage_factor/1000 #from Precip_heat calculation of the max amount
CANLIQ  = 0.0
stored_water = 1.0
out_box2_to_box3 = 0
detention_storage = 0
PRECIP_input = 0.006 #mm/s

x = np.linspace(0,np.pi,10000)
x = np.pad(x,(0,24000),'constant')

EC = 0.000000001
ET = 0.00000001
EDIR = 0.000000005

PRECIP_input = .0000001*np.sin(x)
stored_water = np.zeros_like(PRECIP_input)
detention_storage = np.zeros_like(PRECIP_input)
CANLIQ = np.zeros_like(PRECIP_input)
runoff_final = np.zeros_like(PRECIP_input)
out_box2_to_box3 = np.zeros_like(PRECIP_input)

##set the intial conditions
CANLIQ[0] = 0.0
stored_water[0] = 0.0
detention_storage[0] = 0.0
runoff_final[0] = 0.0
top_storage = 0.0
out_box2_to_box3[0] = 0.0




for i in range (1,len(x)):

	Precip_dripped_to_ground, Precip_throughfall_to_ground, CANLIQ[i] = interception_throughfall(CANLIQ[i-1],PRECIP_input[i],MAX_amount_intercepted,EC)

	stored_water[i], out_box2_to_box3[i] = retention_1(Precip_dripped_to_ground+Precip_throughfall_to_ground,stored_water[i-1],ET,EDIR)

	detention_storage[i], runoff_final[i] = retention_2(out_box2_to_box3[i],detention_storage[i-1])


x = np.linspace(0,len(x),len(x))
plt.subplot(511)
plt.plot(x,PRECIP_input,linewidth=3.0,color='K',label='Desgin Rain Storm: ' +str(round(np.sum(PRECIP_input*DT*1000),3)) + ' mm',alpha=alpha)
plt.ylabel(r'Rain Rate $[mm \enspace s^{-1}]$')
plt.legend()
plt.subplot(512)
plt.plot(x,CANLIQ,linewidth=2.0,color='g',label='Canopy Liquid Content: ' +  str(round(CANLIQ[-1],3)) + ' mm',alpha=alpha)
plt.ylabel(r'Canopy Storage $[mm]$')
plt.legend()
plt.subplot(513)
plt.plot(x,stored_water,linewidth=2.0, color='b',label='Growing Media Moisture Content: ' + str(round(stored_water[-1],3)) + ' mm' ,alpha=alpha)
plt.ylabel(r'Non-Drainable Storage $[mm]$')
plt.legend()
plt.subplot(514)
plt.plot(x,detention_storage,linewidth=2.0,color='darkorchid',label='Detention Pond Storage: ' + str(round(detention_storage[-1],3)) + ' mm',alpha=alpha)
plt.ylabel(r'Drainable Storage $[mm]$')
plt.legend()
plt.subplot(515)
plt.plot(x,runoff_final,linewidth=2.0,color='r',label='System Outflow: ' + str(round(np.sum(runoff_final),3)),alpha=alpha)
plt.ylabel(r'Outflow $[mm]$')
plt.legend()

print(np.sum(PRECIP_input*DT))
print(np.sum(runoff_final))
print(detention_storage[-1])
print(stored_water[-1])
print(CANLIQ[-1])
# plt.tight_layout()
plt.show()

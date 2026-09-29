import os
import numpy as np
import ttide
import netCDF4 as NC
import datetime
from datetime import datetime
import matplotlib.dates as mdates
import pandas as pd

workpath='/work/oda/ag15419/tmp/ttide_obs/'
flag_modorobs=2 # 0= mod; 1=emodnet obs; 2=ispra obs

######################################
if flag_modorobs == 2:
 # 2) ISPRA OBS:
 # All ISPRA TG - read list:
 stations_obs=['ancona','carloforte','catania','imperia','lampedusa','livorno','messina','ortona','palinuro','portoempedocle','portotorres','reggiocalabria','trieste','venezia','vieste']
 #
 # Lit comparison ISPRA
 #stations_obs=['ancona','carloforte','catania','lampedusa','livorno','ortona','portoempedocle','reggiocalabria','trieste','venezia','vieste']

elif flag_modorobs == 1:
 # 1) EMODNET OBS:
 # All EMODNET:
 stations_obs = ['Ajaccio','Algeciras','Almeria','Barcelona','Centuri','Ginostra','Ibiza','IleRousse','iske','LaFigueirette','Malaga','Marseille','Melilla','Monaco','Motril','PalmadeMallorca','PortLaNouvelle','PortVendres','Sagunto','Sete','Solenzara','Tarifa','Valencia','zygi1']
 #
 # Lit comparison EMODNET
 #stations_obs = ['Almeria','Malaga','Marseille','Tarifa']

elif flag_modorobs == 0:
 # 0) MOD
 # All TG - plot list 4 MOD extraction :
 stations_mod=['AjaccioTG','AlmeriaTG','BarcelonaTG','CenturiTG','GinostraTG','IbizaTG','IleRousseTG','LaFigueiretteTG','MarseilleTG','MonacoTG','PalmadeMallorcaTG','PortLaNouvelleTG','PortVendresTG','SaguntoTG','SeteTG','SolenzaraTG','ValenciaTG','carloforte','imperia','lampedusa','livorno','palinuro','portoempedocle','portotorres','catania','messina','reggiocalabria','AlgecirasTG','MalagaTG','MelillaTG','MotrilTG','TarifaTG','ancona','ortona','trieste','venezia','vieste','iskeTG','zygi1TG']
 #
 # 2018 TG only:
 #stations_mod=['iskeTG','zygi1TG']
 #
 #Lit comparison MOD
 #stations_mod=['ancona','carloforte','catania','lampedusa','livorno','ortona','portoempedocle','reggiocalabria','trieste','venezia','vieste','Almeria','Malaga','Marseille','Tarifa']


tidal_comp=['M2','S2','K1','O1','N2','P1','Q1','K2']

# MOD
if flag_modorobs == 0:
 # Loop on TG
 for stn in range (0,len(stations_mod)):
   # Open time_series file and get values
   fT1_mod = NC.Dataset(workpath+stations_mod[stn]+'_mod_Tides8_v31.nc','r')
   xin = fT1_mod.variables['sossheig'][:,0,0] *100.0 # want cm not meters
   latitudes=fT1_mod.variables['lat'][0]
   # Set the time
   if stations_mod[stn] == 'zygi1TG' or stations_mod[stn] == 'iskeTG' or stations_mod[stn] == 'GinostraTG':
      time=datetime(2018, 7, 1, 0, 30, 0) 
   else:
      time=datetime(2017, 7, 1, 0, 30, 0)
   print ('Extracting.. STZ LAT TIME', stations_mod[stn] , latitudes, time)
   # Extract Amp and Pha
   harmanOUT = ttide.t_tide(np.array(xin), dt=1., stime=time, lat=latitudes, constitnames=tidal_comp, out_style='classic', outfile =workpath+stations_mod[stn]+'_mod_tt.txt')

# EMODNET OBS
if flag_modorobs == 1:
 # Loop on TG
 for stn in range (0,len(stations_obs)):
   # Open file and get values
   fT1_obs = NC.Dataset(workpath+stations_obs[stn]+'TG_obs.nc_h.nc','r')
   xin = fT1_obs.variables["sossheig"][:,0]*100.0 # want cm not meters
   # Read latitude from model file
   fT1_mod = NC.Dataset(workpath+stations_obs[stn]+'TG_mod_Tides8_v31.nc','r')
   latitudes=fT1_mod.variables['lat'][0]
   # Set the time
   if stations_obs[stn] == 'zygi1' or stations_obs[stn] == 'iske':
      time=datetime(2018, 7, 1, 0, 30, 0)
   elif stations_obs[stn] == 'Ginostra':
      time=datetime(2018, 7, 4, 8, 40, 0)
   elif stations_obs[stn] == 'Toulon':
      time=datetime(2017, 7, 5, 12, 30, 0)
   else:
      time=datetime(2017, 7, 1, 0, 30, 0)
   print ('Extracting.. STZ LAT TIME', stations_obs[stn] , latitudes, time)   
   # Extract Amp and Pha
   harmanOUT = ttide.t_tide(np.array(xin), dt=1, stime=time, lat=latitudes, constitnames=tidal_comp, out_style='classic', outfile =workpath+stations_obs[stn]+'_obs_tt.txt')

# ISPRA OBS
if flag_modorobs == 2:
 # Loop on TG
 for stn in range (0,len(stations_obs)):
   # Open file and get values
   fT1_obs = pd.read_csv(workpath+'obs_'+stations_obs[stn]+'.csv',sep=';',usecols=['idNum', 'station', 'year', 'month', 'day', 'hour', 'minu', 'sec', 'value'])
   xin = fT1_obs.value[:]
   print(xin)
   # Read latitude from model file
   fT1_mod = NC.Dataset(workpath+stations_obs[stn]+'_mod_Tides8_v31.nc','r')
   latitudes=fT1_mod.variables['lat'][0]
   # Set the time
   time=datetime(2017, 6, 30, 10, 30, 0)
   print ('Extracting.. STZ LAT TIME', stations_obs[stn] , latitudes, time)
   # Extract Amp and Pha
   harmanOUT = ttide.t_tide(np.array(xin), dt=1., stime=time, lat=latitudes, constitnames=tidal_comp, out_style='classic', outfile =workpath+stations_obs[stn]+'_obs_tt.txt')
#

# Interpolation to hourly data
# for TOBERED in $(ls *_obs.nc ); do echo $TOBERED ; cdo inttime,2017-07-01,00:30:00,1hour $TOBERED ${TOBERED}_h.nc ; done
# cdo inttime,2018-07-04,08:40:00,1hour GinostraTG_obs.nc GinostraTG_obs.nc_h.nc 
# cdo inttime,2018-07-01,00:30:00,1hour iskeTG_obs.nc iskeTG_obs.nc_h.nc
# cdo inttime,2018-07-01,00:30:00,1hour zygi1TG_obs.nc zygi1TG_obs.nc_h.nc
# cdo inttime,2017-07-05,12:30:00,1hour ToulonTG_obs.nc ToulonTG_obs.nc_h.nc 


# Built arrays from file

# ISPRA:
#for t_comp in 'M2' 'S2' 'K1' 'O1' 'N2' 'P1' 'Q1' 'K2'; do echo -en "\nI_${t_comp}_amp_mod=["; for stz in 'ancona' 'carloforte' 'catania' 'imperia' 'lampedusa' 'livorno' 'messina' 'ortona' 'palinuro' 'portoempedocle' 'portotorres' 'reggiocalabria' 'trieste' 'venezia' 'vieste'; do AMP=$( grep $t_comp /work/oda/ag15419/tmp/ttide_obs/${stz}_mod_tt.txt | cut -f 2,3 -d"." | cut -f 5,6 -d" " ); echo -en "${AMP},"; done; echo -en "]"; done
#for t_comp in 'M2' 'S2' 'K1' 'O1' 'N2' 'P1' 'Q1' 'K2'; do echo -en "\nI_${t_comp}_amp_obs=["; for stz in 'ancona' 'carloforte' 'catania' 'imperia' 'lampedusa' 'livorno' 'messina' 'ortona' 'palinuro' 'portoempedocle' 'portotorres' 'reggiocalabria' 'trieste' 'venezia' 'vieste'; do AMP=$( grep $t_comp /work/oda/ag15419/tmp/ttide_obs/${stz}_obs_tt.txt | cut -f 2,3 -d"." | cut -f 5,6 -d" " ); echo -en "${AMP},"; done; echo -en "]"; done

#for t_comp in 'M2' 'S2' 'K1' 'O1' 'N2' 'P1' 'Q1' 'K2'; do echo -en "\nI_${t_comp}_pha_mod=["; for stz in 'ancona' 'carloforte' 'catania' 'imperia' 'lampedusa' 'livorno' 'messina' 'ortona' 'palinuro' 'portoempedocle' 'portotorres' 'reggiocalabria' 'trieste' 'venezia' 'vieste'; do PHA=$( grep $t_comp /work/oda/ag15419/tmp/ttide_obs/${stz}_mod_tt.txt | cut -f 4,5 -d"." | cut -f 3,4,5,6,7 -d" "); echo -en "${PHA},"; done; echo -en "]"; done
#for t_comp in 'M2' 'S2' 'K1' 'O1' 'N2' 'P1' 'Q1' 'K2'; do echo -en "\nI_${t_comp}_pha_obs=["; for stz in 'ancona' 'carloforte' 'catania' 'imperia' 'lampedusa' 'livorno' 'messina' 'ortona' 'palinuro' 'portoempedocle' 'portotorres' 'reggiocalabria' 'trieste' 'venezia' 'vieste'; do PHA=$( grep $t_comp /work/oda/ag15419/tmp/ttide_obs/${stz}_obs_tt.txt | cut -f 4,5 -d"." | cut -f 3,4,5,6,7 -d" "); echo -en "${PHA},"; done; echo -en "]"; done

# ISPRA for comp wrt literature:
#for t_comp in 'M2' 'S2' 'K1' 'O1' 'N2' 'P1' 'Q1' 'K2'; do echo -en "\nI_${t_comp}_amp_mod=["; for stz in 'ancona' 'carloforte' 'catania'  'lampedusa' 'livorno'  'ortona'  'portoempedocle'  'reggiocalabria' 'trieste' 'venezia' 'vieste'; do AMP=$( grep $t_comp /work/oda/ag15419/tmp/ttide_obs/${stz}_mod_tt.txt | cut -f 2,3 -d"." | cut -f 5,6 -d" " ); echo -en "${AMP},"; done; echo -en "]"; done
#for t_comp in 'M2' 'S2' 'K1' 'O1' 'N2' 'P1' 'Q1' 'K2'; do echo -en "\nI_${t_comp}_amp_obs=["; for stz in 'ancona' 'carloforte' 'catania'  'lampedusa' 'livorno'  'ortona'  'portoempedocle'  'reggiocalabria' 'trieste' 'venezia' 'vieste'; do AMP=$( grep $t_comp /work/oda/ag15419/tmp/ttide_obs/${stz}_obs_tt.txt | cut -f 2,3 -d"." | cut -f 5,6 -d" " ); echo -en "${AMP},"; done; echo -en "]"; done

#for t_comp in 'M2' 'S2' 'K1' 'O1' 'N2' 'P1' 'Q1' 'K2'; do echo -en "\nI_${t_comp}_pha_mod=["; for stz in 'ancona' 'carloforte' 'catania'  'lampedusa' 'livorno'  'ortona'  'portoempedocle'  'reggiocalabria' 'trieste' 'venezia' 'vieste'; do PHA=$( grep $t_comp /work/oda/ag15419/tmp/ttide_obs/${stz}_mod_tt.txt | cut -f 4,5 -d"." | cut -f 3,4,5,6,7 -d" "); echo -en "${PHA},"; done; echo -en "]"; done
#for t_comp in 'M2' 'S2' 'K1' 'O1' 'N2' 'P1' 'Q1' 'K2'; do echo -en "\nI_${t_comp}_pha_obs=["; for stz in 'ancona' 'carloforte' 'catania'  'lampedusa' 'livorno'  'ortona'  'portoempedocle'  'reggiocalabria' 'trieste' 'venezia' 'vieste'; do PHA=$( grep $t_comp /work/oda/ag15419/tmp/ttide_obs/${stz}_obs_tt.txt | cut -f 4,5 -d"." | cut -f 3,4,5,6,7 -d" "); echo -en "${PHA},"; done; echo -en "]"; done

# EMODNET:
#for t_comp in 'M2' 'S2' 'K1' 'O1' 'N2' 'P1' 'Q1' 'K2'; do echo -en "\n${t_comp}_amp_mod=["; for stz in 'Ajaccio' 'Algeciras' 'Almeria' 'Barcelona' 'Centuri' 'Ginostra' 'Ibiza' 'IleRousse' 'iske' 'LaFigueirette' 'Malaga' 'Marseille' 'Melilla' 'Monaco' 'Motril' 'PalmadeMallorca' 'PortLaNouvelle' 'PortVendres' 'Sagunto' 'Sete' 'Solenzara' 'Tarifa' 'Valencia' 'zygi1' ; do AMP=$( grep $t_comp /work/oda/ag15419/tmp/ttide_obs/${stz}TG_mod_tt.txt | cut -f 2,3 -d"." | cut -f 5,6 -d" " ); echo -en "${AMP},"; done; echo -en "]"; done
#for t_comp in 'M2' 'S2' 'K1' 'O1' 'N2' 'P1' 'Q1' 'K2'; do echo -en "\n${t_comp}_amp_obs=["; for stz in 'Ajaccio' 'Algeciras' 'Almeria' 'Barcelona' 'Centuri' 'Ginostra' 'Ibiza' 'IleRousse' 'iske' 'LaFigueirette' 'Malaga' 'Marseille' 'Melilla' 'Monaco' 'Motril' 'PalmadeMallorca' 'PortLaNouvelle' 'PortVendres' 'Sagunto' 'Sete' 'Solenzara' 'Tarifa' 'Valencia' 'zygi1' ; do AMP=$( grep $t_comp /work/oda/ag15419/tmp/ttide_obs/${stz}_obs_tt.txt | cut -f 2,3 -d"." | cut -f 5,6 -d" " ); echo -en "${AMP},"; done; echo -en "]"; done

#for t_comp in 'M2' 'S2' 'K1' 'O1' 'N2' 'P1' 'Q1' 'K2'; do echo -en "\n${t_comp}_pha_mod=["; for stz in 'Ajaccio' 'Algeciras' 'Almeria' 'Barcelona' 'Centuri' 'Ginostra' 'Ibiza' 'IleRousse' 'iske' 'LaFigueirette' 'Malaga' 'Marseille' 'Melilla' 'Monaco' 'Motril' 'PalmadeMallorca' 'PortLaNouvelle' 'PortVendres' 'Sagunto' 'Sete' 'Solenzara' 'Tarifa' 'Valencia' 'zygi1' ; do PHA=$( grep $t_comp /work/oda/ag15419/tmp/ttide_obs/${stz}TG_mod_tt.txt | cut -f 4,5 -d"." | cut -f 3,4,5,6,7 -d" "); echo -en "${PHA},"; done; echo -en "]"; done
#for t_comp in 'M2' 'S2' 'K1' 'O1' 'N2' 'P1' 'Q1' 'K2'; do echo -en "\n${t_comp}_pha_obs=["; for stz in 'Ajaccio' 'Algeciras' 'Almeria' 'Barcelona' 'Centuri' 'Ginostra' 'Ibiza' 'IleRousse' 'iske' 'LaFigueirette' 'Malaga' 'Marseille' 'Melilla' 'Monaco' 'Motril' 'PalmadeMallorca' 'PortLaNouvelle' 'PortVendres' 'Sagunto' 'Sete' 'Solenzara' 'Tarifa' 'Valencia' 'zygi1' ; do PHA=$( grep $t_comp /work/oda/ag15419/tmp/ttide_obs/${stz}TG_obs_tt.txt | cut -f 4,5 -d"." | cut -f 3,4,5,6,7 -d" "); echo -en "${PHA},"; done; echo -en "]"; done



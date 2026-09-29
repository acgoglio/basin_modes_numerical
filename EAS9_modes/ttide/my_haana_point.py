#
# Script for HARMONIC ANALYSIS AND POST_PROC
#
# imports
import matplotlib.pyplot as plt
import numpy as np
import netCDF4 as NC
from scipy.optimize import curve_fit
from scipy import stats
import collections
import pandas as pd
import csv
import math
import datetime
from operator import itemgetter 
import plotly
from plotly import graph_objects as go # for bar plot
#
#######################################################
where_box='Med' # Med or AtlBox
path = '/work/oda/ag15419/tmp/tides8_v31/HA_v31/tt_tide/' # workdir path
mod_file_template='Tides8_v31'
iniend_dates= '01/07/2017 - 31/12/2017'
dates_lab='20170701_20171231'
name = path
tpxo_flag = 1 # to compare also wrt TPXO data
cos_pha = 0 # to compare cosine values of phases instead of abs values (to avoid +-360*n corrections..)
flag_15stats = 1 # to compare results with literature
flag_tab = 1 # For name strings length (if =1 the name is complete else only the first 3 characters are considered) 
ttide_flag = 2 # To use ttide results intsead of fit from Salish Sea script

########################################################
# DO NOT CHANGE THE CODE BELOW THIS LINE!!!
########################################################
# Check on flags..
if where_box=='AtlBox':
   tpxo_flag = 0
   flag_15stats = 0
   ttide_flag = 0
if where_box=='Med' and flag_15stats==1 :
   tpxo_flag == 1
   flag_tab==1
#if where_box=='Med' and flag_tab==0 :
#   flag_15stats = 0

#

if tpxo_flag == 1 or tpxo_flag == 0:
   TPXO_M2=[0.075,0.098,0.041,0.082,0.117,0.011,0.080,0.073,0.064,0.075,0.025,0.057,0.055,0.017,0.061,0.103,0.016,0.060,0.078,0.064,0.091,0.120,0.040,0.072,0.061,0.039,0.046,0.285,0.188,0.117,0.155,0.405,0.082,0.067,0.278,0.231,0.096,0.124,0.099]
   TPXO_P_M2=[-137.2,47.5,-150.9,-138.9,-126.7,-149.4,-138.5,-141.7,-145.8,-141.0,-158.6,-150.1,-149.2,-162.6,-149.2,-132.2,-159.4,-131.0,-140.2,40.3,-133.2,-128.2,68.1,-136.1,57.1,60.8,67.7,50.4,50.5,66.2,49.4,52.4,-56.9,68.0,-108.3,-99.9,72.3,-127.1,-120.6]
   TPXO_S2=[0.030,0.038,0.015,0.033,0.045,0.004,0.032,0.028,0.025,0.029,0.009,0.022,0.021,0.004,0.024,0.038,0.004,0.026,0.031,0.048,0.035,0.046,0.033,0.030,0.034,0.024,0.026,0.105,0.071,0.046,0.059,0.147,0.031,0.042,0.160,0.128,0.059,0.070,0.056]
   TPXO_P_S2=[-117.4,73.0,-126.9,-118.9,-105.6,-119.1,-118.2,-122.0,-123.7,-121.0,-130.1,-125.8,-125.2,-135.0,-126.4,-107.0,-146.3,-107.9,-121.6,54.3,-115.1,-107.7,68.8,-115.7,67.2,65.6,73.3,75.6,76.2,92.5,75.2,78.9,-51.6,75.0,-107.4,-94.2,81.2,-116.2,-109.7]
   TPXO_N2=[0.016,0.019,0.010,0.018,0.024,0.004,0.017,0.016,0.014,0.017,0.006,0.013,0.013,0.004,0.014,0.022,0.004,0.014,0.017,0.009,0.020,0.025,0.008,0.016,0.011,0.007,0.008,0.059,0.038,0.025,0.032,0.087,0.007,0.008,0.042,0.036,0.012,0.021,0.016]    
   TPXO_P_N2=[-142.4,33.0,-150.9,-141.8,-138.4,-139.2,-144.5,-145.3,-146.3,-147.3,-149.0,-153.4,-151.4,-153.4,-149.7,-144.2,-153.4,-135.2,-144.5,40.6,-139.1,-139.9,97.6,-140.2,63.4,74.1,76.0,34.0,35.0,56.3,34.7,36.0,-90.0,113.2,-94.1,-75.1,94.8,-129.1,-124.4]
   TPXO_K2=[0.009,0.010,0.004,0.009,0.012,0.002,0.009,0.008,0.007,0.008,0.003,0.006,0.006,0.001,0.007,0.010,0.001,0.007,0.009,0.015,0.010,0.012,0.009,0.009,0.009,0.006,0.007,0.030,0.021,0.013,0.016,0.041,0.006,0.013,0.031,0.024,0.016,0.023,0.019]
   TPXO_P_K2=[-116.6,73.3,-116.6,-116.6,-104.0,-104.4,-116.6,-119.7,-123.7,-119.7,-108.4,-121.0,-121.0,-135.0,-123.7,-106.7,-135.0,-105.9,-116.6,42.3,-114.0,-104.0,58.0,-110.6,54.5,51.3,56.3,76.4,76.0,94.4,76.0,75.6,-38.7,67.4,-97.4,-82.4,71.6,-110.0,-105.5]
   TPXO_K1=[0.037,0.039,0.037,0.037,0.031,0.038,0.037,0.036,0.036,0.036,0.038,0.037,0.037,0.039,0.036,0.029,0.039,0.037,0.036,0.003,0.036,0.031,0.018,0.037,0.019,0.016,0.018,0.030,0.037,0.045,0.037,0.019,0.129,0.091,0.173,0.162,0.057,0.030,0.024]
   TPXO_P_K1=[-176.9,154.1,166.0,-175.3,-165.1,164.9,-175.4,178.4,173.7,-180.0,167.8,169.0,169.0,160.6,170.5,-178.0,160.6,-180.0,-178.4,-28.0,-178.0,-165.1,86.8,-176.9,55.5,68.2,58.4,135.0,145.9,146.0,149.3,104.5,72.4,67.4,68.6,72.3,71.5,-86.2,-85.2]
   TPXO_O1=[0.016,0.022,0.023,0.015,0.011,0.022,0.016,0.018,0.020,0.017,0.021,0.022,0.022,0.025,0.022,0.015,0.025,0.017,0.017,0.009,0.014,0.011,0.009,0.016,0.009,0.008,0.009,0.011,0.019,0.023,0.020,0.008,0.035,0.026,0.044,0.044,0.017,0.021,0.017]
   TPXO_P_O1=[111.8,114.2,105.3,113.2,131.2,108.9,111.8,106.4,104.7,107.4,106.7,105.9,103.4,106.3,103.4,109.7,107.8,110.6,107.4,110.6,114.8,131.2,83.7,111.8,54.5,60.3,54.5,131.0,117.9,97.4,116.6,-109.6,61.3,57.5,58.7,64.5,59.0,-109.3,-110.6]
   TPXO_P1=[0.012,0.013,0.013,0.012,0.010,0.013,0.012,0.012,0.012,0.012,0.013,0.013,0.013,0.013,0.013,0.009,0.013,0.012,0.012,0.001,0.012,0.010,0.007,0.012,0.006,0.005,0.006,0.008,0.011,0.014,0.012,0.006,0.037,0.027,0.050,0.048,0.017,0.012,0.009]
   TPXO_P_P1=[175.2,151.4,157.4,175.2,-174.3,157.4,175.2,170.5,166.0,170.5,161.6,161.6,161.6,153.4,161.6,173.7,153.4,170.5,170.5,45.0,171.8,-174.3,74.1,175.2,45.0,53.1,45.0,129.8,142.1,140.7,145.0,90.0,71.1,66.3,72.6,79.5,69.4,-90.0,-83.8]
   TPXO_Q1=[0.003,0.003,0.004,0.003,0.002,0.004,0.004,0.004,0.004,0.004,0.004,0.004,0.004,0.004,0.004,0.004,0.004,0.003,0.004,0.003,0.002,0.002,0.002,0.003,0.002,0.001,0.001,0.002,0.002,0.002,0.002,0.005,0.005,0.004,0.008,0.007,0.003,0.003,0.002]
   TPXO_P_Q1=[45.0,90.0,56.3,18.4,26.6,56.3,33.7,33.7,33.7,33.7,56.3,56.3,56.3,56.3,45.0,33.7,63.4,45.0,33.7,71.6,26.6,26.6,63.4,45.0,63.4,45.0,45.0,-180.0,116.6,26.6,90.0,-142.5,111.8,90.0,113.2,122.3,90.0,-108.4,-116.6]

# Values from literature (Almeria,Marsiglia,Carloforte,Lampedusa,Livorno,P.Empedocle,Catania,R.Calabria,Malaga,Tarifa/Gibraltar,Ancona,Ortona,Trieste,Venezia,Vieste)
# [Palma et al.; 2020]
PALMA_M2=['',8.0,'','','','','','',7.1,'','','','','','','','',6.9,'',6.6,9.4,'',3.7,'',6.8,'',4.4,'',16.5,'','',29.9,6.1,5.8,20.8,18.5,8.0,'','']
PALMA_P_M2=['',46.6,'','','','','','',223.6,'','','','','','','','',232.3,'',33.2,229.4,'',65.4,'',61.6,'',67.7,'',48.4,'','',50.3,302.3,54.6,247.3,257.7,65.8,'','']
PALMA_d_M2=['',1.2,'','','','','','',0.8,'','','','','','','','',0.4,'',1.4,0.9,'',1.2,'',0.4,'',1.9,'',2.7,'','',2.3,0.5,1.2,5.5,4.9,0.7,'','']
PALMA_S2=['',3.4,'','','','','','',2.5,'','','','','','','','',2.7,'',3.9,3.5,'',2.6,'',4.2,'',2.4,'',6.6,'','',11.1,3.4,4.2,13.0,11.5,5.3,'','']
PALMA_P_S2=['',68.8,'','','','','','',240.7,'','','','','','','','',251.0,'',58.3,244.4,'',74.7,'',65.8,'',66.3,'',71.0,'','',74.9,310.9,61.3,247.6,258.0,69.1,'','']
PALMA_d_S2=['',0.9,'','','','','','',0.5,'','','','','','','','',0.2,'',0.3,0.1,'',0.6,'',0.8,'',0.7,'',1.0,'','',0.7,0.3,1.1,2.9,3.0,1.3,'','']
PALMA_K1=['',3.1,'','','','','','',3.5,'','','','','','','','',3.8,'',0.9,4.1,'',1.7,'',3.6,'',2.0,'',2.6,'','',1.8,12.9,8.1,17.0,16.8,5.0,'','']
PALMA_P_K1=['',185.0,'','','','','','',201.4,'','','','','','','','',204.7,'',15.2,205.8,'',99.8,'',45.6,'',62.8,'',178.9,'','',169.9,65.7,64.2,50.2,55.4,72.4,'','']
PALMA_d_K1=['',1.5,'','','','','','',1.5,'','','','','','','','',1.7,'',0.2,1.7,'',0.3,'',2.1,'',0.8,'',1.2,'','',1.3,1.4,1.7,3.9,2.0,1.0,'','']
PALMA_O1=['',2.3,'','','','','','',1.9,'','','','','','','','',1.7,'',0.4,1.8,'',1.1,'',1.2,'',1.0,'',2.1,'','',1.7,4.2,2.0,5.2,5.2,1.8,'','']
PALMA_P_O1=['',132.4,'','','','','','',128.1,'','','','','','','','',135.9,'',48.4,136.4,'',63.7,'',44.9,'',42.2,'',137.4,'','',153.8,44.0,127.8,28.9,34.6,47.4,'','']
PALMA_d_O1=['',0.6,'','','','','','',0.9,'','','','','','','','',1.3,'',1.0,0.9,'',0.4,'',0.1,'',0.1,'',0.4,'','',0.8,1.2,0.8,1.6,1.5,0.3,'','']
# [Tsimplis et al.; 1995]
TSIMPLIS_M2=['',9.2,'','','','','','',7.6,'','','','','','','','',7.3,'',6.4,9.8,'',3.9,'',5.5,'',3.2,'',20.1,'','',29.6,6.6,6.4,24.8,21.5,8.2,'','']
TSIMPLIS_P_M2=['',43.0,'','','','','','',221.0,'','','','','','','','',229.0,'',31.0,226.0,'',75.0,'',55.0,'',65.0,'',46.0,'','',45.0,288.0,54.0,242.0,255.0,64.0,'','']
TSIMPLIS_d_M2=['',1.5,'','','','','','',0.8,'','','','','','','','',0.9,'',1.6,1.4,'',0.7,'',1.1,'',3.0,'',4.0,'','',0.7,1.6,1.1,3.2,2.7,0.6,'','']
TSIMPLIS_S2=['',3.4,'','','','','','',2.6,'','','','','','','','',2.7,'',4.9,3.5,'',3.3,'',3.1,'',2.1,'',6.9,'','',9.9,3.8,4.9,16.6,14.4,5.9,'','']
TSIMPLIS_P_S2=['',73.0,'','','','','','',240,'','','','','','','','',251.0,'',45.0,243,'',62,'',55.0,'',54.0,'',77.0,'','',77,293.0,59.0,242.0,254.0,66.0,'','']
TSIMPLIS_d_S2=['',0.7,'','','','','','',0.6,'','','','','','','','',0.2,'',1.2,0.1,'',0.9,'',0.8,'',1.2,'',0.3,'','',1.2,1.5,1.5,3.9,2.6,1.8,'','']
TSIMPLIS_K1=['',4.2,'','','','','','',4.3,'','','','','','','','',4.4,'',0.7,4.6,'',2.3,'',1.9,'',1.5,'',3.5,'','',2.6,14.3,9.7,19.0,18.6,6.7,'','']
TSIMPLIS_P_K1=['',169.0,'','','','','','',173,'','','','','','','','',178.0,'',1.0,181.0,'',90.0,'',49.0,'',65,'',171.0,'','',180,58.0,57.0,45.0,51.0,59.0,'','']
TSIMPLIS_d_K1=['',1.4,'','','','','','',1.4,'','','','','','','','',1.2,'',0.2,0.6,'',0.5,'',0.5,'',0.5,'',1.0,'','',2.0,3.5,2.7,4.5,3.2,2.5,'','']
TSIMPLIS_O1=['',2.5,'','','','','','',2.2,'','','','','','','','',1.9,'',0.7,1.8,'',1.3,'',1.1,'',1.0,'',2.1,'','',1.8,4.8,3.3,6.3,6.2,2.3,'','']
TSIMPLIS_P_O1=['',122.0,'','','','','','',107,'','','','','','','','',113.0,'',73.0,115.0,'',65.0,'',45.0,'',50,'',133.0,'','',156.0,1.0,49.0,39.0,45,49,'','']
TSIMPLIS_d_O1=['',0.5,'','','','','','',0.3,'','','','','','','','',0.6,'','',0.3,'',0.3,'',0.0,'',0.2,'',0.3,'','',0.9,1.0,0.3,0.5,0.8,0.8,'','']

# ============ ISPRA ===============

# STZ EX

if where_box == 'Med':

   # For comparison with Tsimplis et al. 1995
   if flag_15stats == 1: 
      print ('Lit comparison..')
      stations_obs = ['ancona','carloforte','catania','lampedusa','livorno','ortona','portoempedocle','reggiocalabria','trieste','venezia','vieste']
      stations_mod = stations_obs
      stations_lab = ['Ancona','Carloforte','Catania','Lampedusa','Livorno','Ortona','P.Empedocle','R.Calabria','Trieste','Venezia','Vieste']
      stations_col=[ 'blue','green','orange','green','green','blue','green','orange','blue','blue','blue']
   else:
   # x 2017 NO GIRNE
      stations_obs = ['ancona','carloforte','catania','imperia','lampedusa','livorno','messina','ortona','palinuro','portoempedocle','portotorres','reggiocalabria','trieste','venezia','vieste']
      stations_mod = ['ancona','carloforte','catania','imperia','lampedusa','livorno','messina','ortona','palinuro','portoempedocle','portotorres','reggiocalabria','trieste','venezia','vieste']
      stations_lab = ['Ancona','Carloforte','Catania','Imperia','Lampedusa','Livorno','Messina','Ortona','Palinuro','P.Empedocle','P.Torres','R.Calabria','Trieste','Venezia','Vieste']
      stations_col=[ 'blue','green','orange','green','green','green','orange','blue','green','green','green','orange','blue','blue','blue' ]

elif where_box == 'AtlBox':
   stations_obs = []
   stations_mod = stations_obs
   stations_lab = stations_obs
   stations_col = stations_obs

##
numsta_obs=len(stations_obs)
numsta_mod=len(stations_mod)
print("ISPRA Stations num:",len(stations_obs),len(stations_mod) )


# TTIDE METHODOLOGY
if ttide_flag == 1:

 I_stations_lab=stations_lab
 I_numsta=numsta_obs
 I_stations_col=stations_col

 #I_%comp%_amp/pha_mod/obs=[]

 if flag_15stats == 0:

  I_M2_amp_mod=[5.8198,6.1121,5.9116,7.4113,5.8557,8.5442,3.0175,6.3895,11.5321,4.0252,6.8800,5.5154,24.8351,21.9406,8.7600]
  I_S2_amp_mod=[3.1273,2.3397,2.9233,2.6831,3.2933,3.1613,1.7678,4.5953,4.5990,2.2132,2.5750,2.7195,15.4861,13.6757,5.6851]
  I_K1_amp_mod=[12.6234,3.6305,1.8898,3.4995,0.7324,3.5033,1.0566,8.7689,2.7299,1.7052,3.5556,1.6563,17.3597,16.7605,5.0634]
  I_O1_amp_mod=[3.9300,1.6296,1.0153,1.6909,0.6475,1.5429,0.5940,2.8863,0.8668,1.1074,1.5922,0.9394,5.0544,4.9113,1.8294]
  I_N2_amp_mod=[1.0118,1.2407,1.0770,1.5522,0.8134,1.7737,0.5206,1.0132,2.3305,0.7642,1.3942,1.0052,4.2636,3.7479,1.4630]
  I_P1_amp_mod=[3.7570,1.1749,0.6204,1.0377,0.1673,1.0540,0.2552,2.5885,0.9421,0.3316,1.1248,0.5453,4.7057,4.8115,1.6146]
  I_Q1_amp_mod=[1.0211,0.2793,0.1945,0.2725,0.1953,0.2372,0.1892,0.6233,0.1484,0.1488,0.2778,0.2200,1.4955,1.4418,0.2975]
  I_K2_amp_mod=[0.9966,0.6661,0.8556,0.7728,1.0068,0.9107,0.4295,1.2103,1.2491,0.6993,0.7838,0.7660,3.8812,3.5495,1.5315]
  
  I_M2_pha_mod=[296.66,230.04,56.01,222.38,42.10,227.54,51.33,63.30,231.18,78.78,225.69,61.68,249.22,257.45,69.35]
  I_S2_pha_mod=[312.14,250.12,62.88,238.39,50.05,244.39,55.09,75.71,250.10,72.74,244.79,67.79,258.94,266.32,79.96]
  I_K1_pha_mod=[68.54,181.78,37.17,180.76,305.15,179.72,15.15,65.77,198.99,83.37,183.49,30.13,57.12,61.56,70.35]
  I_O1_pha_mod=[56.85,113.04,41.56,108.40,94.59,113.95,22.94,54.43,122.05,67.71,110.07,40.74,48.52,51.46,57.72]
  I_N2_pha_mod=[291.95,215.86,53.85,210.33,49.26,215.66,60.74,61.94,217.62,96.32,212.43,61.14,248.05,257.49,67.49]
  I_P1_pha_mod=[61.54,183.16,23.38,187.43,343.03,192.12,27.68,59.06,207.00,71.38,188.65,39.82,44.41,51.82,64.61]
  I_Q1_pha_mod=[94.98,50.81,42.47,41.28,35.24,60.41,18.65,92.57,1.21,57.55,42.06,37.04,76.18,80.88,88.37]
  I_K2_pha_mod=[302.50,244.50,56.00,235.01,48.92,239.71,47.84,59.19,246.27,66.21,240.80,54.06,243.76,253.96,66.90]
  
  I_M2_amp_obs=[1.4066,1.6148,0.3021,1.5422,0.5292,0.9806,0.4420,0.5887,1.0634,0.4765,0.4532,0.3625,9.0177,4.8987,0.3291]
  I_S2_amp_obs=[2.7891,2.5459,1.0957,2.4888,0.8004,1.3132,1.1680,0.7857,1.6257,1.3499,0.6870,0.8180,16.7153,12.1345,0.3884]
  I_K1_amp_obs=[11.9488,4.0061,0.3439,3.5316,0.1420,2.6768,0.2647,3.9043,0.9101,0.5544,1.0086,0.3104,17.0437,16.5540,0.2956]
  I_O1_amp_obs=[0.8322,0.5271,0.2279,0.4238,0.1217,0.1366,0.2789,0.3328,0.1719,0.2229,0.2382,0.0991,1.8671,0.6870,0.0854]
  I_N2_amp_obs=[0.4597,0.3368,0.1468,0.3379,0.0778,0.3797,0.2722,0.4686,0.3149,0.5629,0.2366,0.2609,4.0261,0.8332,0.1675]
  I_P1_amp_obs=[1.0455,0.7843,0.3329,0.7798,0.1825,0.7578,0.0769,2.9403,0.4059,0.4655,0.6899,0.2981,4.5756,3.3568,0.1971]
  I_Q1_amp_obs=[0.3882,0.1716,0.0900,0.0778,0.0690,0.3045,0.1233,0.2893,0.2303,0.3045,0.0534,0.0910,1.0502,0.5161,0.2628]
  I_K2_amp_obs=[1.1329,0.5483,2.8830,0.7843,4.2492,2.7491,2.5282,4.8959,4.9595,4.4329,3.2774,2.8670,2.0959,2.7229,0.4593]
  
  I_M2_pha_obs=[23.43,323.13,189.02,317.86,196.23,236.41,132.17,191.26,26.13,334.22,28.75,152.17,312.35,339.54,113.51]
  I_S2_pha_obs=[305.06,238.20,78.82,235.95,353.56,196.35,27.87,324.70,123.37,87.27,234.29,12.23,230.60,241.62,228.69]
  I_K1_pha_obs=[78.12,196.77,59.46,198.78,356.78,257.86,253.98,71.59,206.46,100.24,200.93,37.41,65.45,70.75,147.80]
  I_O1_pha_obs=[105.35,177.67,172.26,178.74,224.00,133.56,155.57,105.37,135.27,172.27,149.03,112.31,75.35,102.29,228.41]
  I_N2_pha_obs=[124.47,22.93,306.23,17.12,257.03,37.81,204.54,281.02,124.18,324.64,137.99,327.25,21.24,46.18,352.70]
  I_P1_pha_obs=[15.48,126.21,343.83,139.86,227.73,12.93,88.96,23.55,94.43,37.80,116.35,6.09,323.45,351.80,86.61]
  I_Q1_pha_obs=[358.67,125.89,244.16,264.49,272.82,291.24,189.90,279.64,112.45,255.56,90.66,305.84,65.12,66.74,133.92]
  I_K2_pha_obs=[59.40,1.46,192.94,4.20,169.67,212.47,129.58,182.35,7.89,183.75,20.63,186.45,314.17,17.54,216.39]

 elif flag_15stats == 1:

  I_M2_amp_mod=[5.8198,6.1121,5.9116,5.8557,8.5442,6.3895,4.0252,5.5154,24.8351,21.9406,8.7600]
  I_S2_amp_mod=[3.1273,2.3397,2.9233,3.2933,3.1613,4.5953,2.2132,2.7195,15.4861,13.6757,5.6851]
  I_K1_amp_mod=[12.6234,3.6305,1.8898,0.7324,3.5033,8.7689,1.7052,1.6563,17.3597,16.7605,5.0634]
  I_O1_amp_mod=[3.9300,1.6296,1.0153,0.6475,1.5429,2.8863,1.1074,0.9394,5.0544,4.9113,1.8294]
  I_N2_amp_mod=[1.0118,1.2407,1.0770,0.8134,1.7737,1.0132,0.7642,1.0052,4.2636,3.7479,1.4630]
  I_P1_amp_mod=[3.7570,1.1749,0.6204,0.1673,1.0540,2.5885,0.3316,0.5453,4.7057,4.8115,1.6146]
  I_Q1_amp_mod=[1.0211,0.2793,0.1945,0.1953,0.2372,0.6233,0.1488,0.2200,1.4955,1.4418,0.2975]
  I_K2_amp_mod=[0.9966,0.6661,0.8556,1.0068,0.9107,1.2103,0.6993,0.7660,3.8812,3.5495,1.5315]
  
  I_M2_pha_mod=[296.66,230.04,56.01,42.10,227.54,63.30,78.78,61.68,249.22,257.45,69.35]
  I_S2_pha_mod=[312.14,250.12,62.88,50.05,244.39,75.71,72.74,67.79,258.94,266.32,79.96]
  I_K1_pha_mod=[68.54,181.78,37.17,305.15,179.72,65.77,83.37,30.13,57.12,61.56,70.35]
  I_O1_pha_mod=[56.85,113.04,41.56,94.59,113.95,54.43,67.71,40.74,48.52,51.46,57.72]
  I_N2_pha_mod=[291.95,215.86,53.85,49.26,215.66,61.94,96.32,61.14,248.05,257.49,67.49]
  I_P1_pha_mod=[61.54,183.16,23.38,343.03,192.12,59.06,71.38,39.82,44.41,51.82,64.61]
  I_Q1_pha_mod=[94.98,50.81,42.47,35.24,60.41,92.57,57.55,37.04,76.18,80.88,88.37]
  I_K2_pha_mod=[302.50,244.50,56.00,48.92,239.71,59.19,66.21,54.06,243.76,253.96,66.90]
  
  I_M2_amp_obs=[1.4066,1.6148,0.3021,0.5292,0.9806,0.5887,0.4765,0.3625,9.0177,4.8987,0.3291]
  I_S2_amp_obs=[2.7891,2.5459,1.0957,0.8004,1.3132,0.7857,1.3499,0.8180,16.7153,12.1345,0.3884]
  I_K1_amp_obs=[11.9488,4.0061,0.3439,0.1420,2.6768,3.9043,0.5544,0.3104,17.0437,16.5540,0.2956]
  I_O1_amp_obs=[0.8322,0.5271,0.2279,0.1217,0.1366,0.3328,0.2229,0.0991,1.8671,0.6870,0.0854]
  I_N2_amp_obs=[0.4597,0.3368,0.1468,0.0778,0.3797,0.4686,0.5629,0.2609,4.0261,0.8332,0.1675]
  I_P1_amp_obs=[1.0455,0.7843,0.3329,0.1825,0.7578,2.9403,0.4655,0.2981,4.5756,3.3568,0.1971]
  I_Q1_amp_obs=[0.3882,0.1716,0.0900,0.0690,0.3045,0.2893,0.3045,0.0910,1.0502,0.5161,0.2628]
  I_K2_amp_obs=[1.1329,0.5483,2.8830,4.2492,2.7491,4.8959,4.4329,2.8670,2.0959,2.7229,0.4593]
  
  I_M2_pha_obs=[23.43,323.13,189.02,196.23,236.41,191.26,334.22,152.17,312.35,339.54,113.51]
  I_S2_pha_obs=[305.06,238.20,78.82,353.56,196.35,324.70,87.27,12.23,230.60,241.62,228.69]
  I_K1_pha_obs=[78.12,196.77,59.46,356.78,257.86,71.59,100.24,37.41,65.45,70.75,147.80]
  I_O1_pha_obs=[105.35,177.67,172.26,224.00,133.56,105.37,172.27,112.31,75.35,102.29,228.41]
  I_N2_pha_obs=[124.47,22.93,306.23,257.03,37.81,281.02,324.64,327.25,21.24,46.18,352.70]
  I_P1_pha_obs=[15.48,126.21,343.83,227.73,12.93,23.55,37.80,6.09,323.45,351.80,86.61]
  I_Q1_pha_obs=[358.67,125.89,244.16,272.82,291.24,279.64,255.56,305.84,65.12,66.74,133.92]
  I_K2_pha_obs=[59.40,1.46,192.94,169.67,212.47,182.35,183.75,186.45,314.17,17.54,216.39]

# SALISH SEA METHODOLOGY or MIX
else:


   # Nodal correction and phase 0 referred to 01/07/2017 computed by means of OTPSnc modified script
   #
   # Phase at t=0 referred to 01/07/2017:
   P_shift=[1.731557,0.0000000,0.1730037,1.558554,6.050721,6.110182,5.877717,3.487600] # 'M2','S2','K1','O1','N2','P1','Q1','K2'
   # Phase nodal correction at t=0 referred to 01/07/2017:
   P_nc=[-0.02035177632486352,0.0,-0.09736291220980081,0.127925145499746,-0.02035177632486352,0.0,0.122784747704008,-0.180230370789512]
   # Tot phase offset
   P_shift=np.array(P_shift)+np.array(P_nc)
   # Hand shift..(??)
   P_shift_hand=[0,0,np.pi,np.pi,np.pi/2.0,np.pi,-np.pi/2.0,0]
   P_shift=P_shift+P_shift_hand
   #
   # Amplitude nodal correction at t=0 referred to 01/07/2017:
   A_nc=[1.03158251188290,1.00000000000000,0.906424071424833,0.846920536108457,1.03158251188290,1.00000000000000,0.849287063114331,0.787346337862811]
   # 
   
   
   ##
   ###constants and fitting
   ### M2
   M2freq = 28.984106 # degrees per hour
   M2freq = M2freq*np.pi/180. # radians per hour
   ###K1
   K1freq = 15.041069*np.pi/180.
   ###O1
   O1freq = 13.943036*np.pi/180.
   ###S2
   S2freq = 30.000002*np.pi/180.
   ###P1
   P1freq = 14.958932*np.pi/180.
   ###N2
   N2freq = 28.439730*np.pi/180.
   ###Q1
   Q1freq = 13.398661*np.pi/180.
   ###K2
   K2freq = 30.082138*np.pi/180.
   ##
   
   ### initial phase calculation
   
   # M2
   M2ft = 1.03669420605307
   M2uvt = -12.5753366566839
   
   # S2
   S2ft = 1.000000
   S2uvt = 0.000000
   
   # N2
   N2ft = 1.0
   N2uvt = 0.0
   
   # K1
   K1ft = 0.886345136935456
   K1uvt = -0.4407793529972482
   
   # O1
   O1ft = 0.813447400642605
   O1uvt = -12.5067505117739
   
   # Q1
   Q1ft = 1.0
   Q1uvt = 0.0
   
   # K2
   K2ft = 1.0
   K2uvt = 0.0
   
   # P1
   P1ft = 1.0
   P1uvt = 0.0
   
   ##
   
   # Add nodal corrections and phase at t=0 (01072017)
   def octuple_obs(x, M2amp_obs, M2pha_obs, K1amp_obs, K1pha_obs, O1amp_obs, O1pha_obs, S2amp_obs, S2pha_obs, P1amp_obs, P1pha_obs, N2amp_obs, N2pha_obs, Q1amp_obs, Q1pha_obs, K2amp_obs, K2pha_obs):
       return (A_nc[0]*M2amp_obs*np.cos(M2freq*x-M2pha_obs*np.pi/180.+P_shift[0])+
               A_nc[2]*K1amp_obs*np.cos(K1freq*x-K1pha_obs*np.pi/180.+P_shift[2])+
               A_nc[3]*O1amp_obs*np.cos(O1freq*x-O1pha_obs*np.pi/180.+P_shift[3])+
               A_nc[1]*S2amp_obs*np.cos(S2freq*x-S2pha_obs*np.pi/180.+P_shift[1])+
               A_nc[5]*P1amp_obs*np.cos(P1freq*x-P1pha_obs*np.pi/180.+P_shift[5])+
               A_nc[4]*N2amp_obs*np.cos(N2freq*x-N2pha_obs*np.pi/180.+P_shift[4])+
               A_nc[6]*Q1amp_obs*np.cos(Q1freq*x-Q1pha_obs*np.pi/180.+P_shift[6])+
               A_nc[7]*K2amp_obs*np.cos(K2freq*x-K2pha_obs*np.pi/180.+P_shift[7]))
   
   def octuple_mod(x, M2amp_mod, M2pha_mod, K1amp_mod, K1pha_mod, O1amp_mod, O1pha_mod, S2amp_mod, S2pha_mod, P1amp_mod, P1pha_mod, N2amp_mod, N2pha_mod, Q1amp_mod, Q1pha_mod, K2amp_mod, K2pha_mod):
       return (A_nc[0]*M2amp_mod*np.cos(M2freq*x-M2pha_mod*np.pi/180.+P_shift[0])+
               A_nc[2]*K1amp_mod*np.cos(K1freq*x-K1pha_mod*np.pi/180.+P_shift[2])+
               A_nc[3]*O1amp_mod*np.cos(O1freq*x-O1pha_mod*np.pi/180.+P_shift[3])+
               A_nc[1]*S2amp_mod*np.cos(S2freq*x-S2pha_mod*np.pi/180.+P_shift[1])+
               A_nc[5]*P1amp_mod*np.cos(P1freq*x-P1pha_mod*np.pi/180.+P_shift[5])+
               A_nc[4]*N2amp_mod*np.cos(N2freq*x-N2pha_mod*np.pi/180.+P_shift[4])+
               A_nc[6]*Q1amp_mod*np.cos(Q1freq*x-Q1pha_mod*np.pi/180.+P_shift[6])+
               A_nc[7]*K2amp_mod*np.cos(K2freq*x-K2pha_mod*np.pi/180.+P_shift[7]))
   
   
   #allocate space for our arrays
   M2_amp_obs=[]; M2_pha_obs=[]; K1_amp_obs=[]; K1_pha_obs=[]
   O1_amp_obs=[]; O1_pha_obs=[]; S2_amp_obs=[]; S2_pha_obs=[]
   P1_amp_obs=[]; P1_pha_obs=[]; N2_amp_obs=[]; N2_pha_obs=[]
   Q1_amp_obs=[]; Q1_pha_obs=[]; K2_amp_obs=[]; K2_pha_obs=[]
   
   M2_amp_mod=[]; M2_pha_mod=[]; K1_amp_mod=[]; K1_pha_mod=[]
   O1_amp_mod=[]; O1_pha_mod=[]; S2_amp_mod=[]; S2_pha_mod=[]
   P1_amp_mod=[]; P1_pha_mod=[]; N2_amp_mod=[]; N2_pha_mod=[]
   Q1_amp_mod=[]; Q1_pha_mod=[]; K2_amp_mod=[]; K2_pha_mod=[]
   
   
   # ISPRA
   for stn in range (0,len(stations_mod)): 
       print('STZ:',stations_mod[stn])
   
       fT1_obs = pd.read_csv(name+'obs_'+stations_obs[stn]+'.csv',sep=';',usecols=['idNum', 'station', 'year', 'month', 'day', 'hour', 'minu', 'sec', 'value'])
       fT1_mod = NC.Dataset(name+stations_mod[stn]+'_mod_'+mod_file_template+'.nc','r') 
       
       # ISPRA
       ssh_obs_arr= fT1_obs.value[:] 
       # NEMO:
       time_mod = (fT1_mod.variables["time_counter"][:]/3600.)-582912.0   # want hours (not seconds) from 1/1/2015 (not 1/1/1950)
       print('time mod:',time_mod)
       ssh_mod = fT1_mod.variables["sossheig"][:,0,0]*100.0 # want cm not meters
       # times
       te_obs= ssh_obs_arr.shape[0]
       te_mod = ssh_mod.shape[0]
       #
       ts_mod=0
       ts_obs=0
   
       # ISPRA
       time_obs=[]
       ssh_obs=[]
       for otime in range (0,te_obs):
           time_obs_el = (datetime.datetime(fT1_obs.year[otime],fT1_obs.month[otime],fT1_obs.day[otime],fT1_obs.hour[otime],fT1_obs.minu[otime])-datetime.datetime(2016, 7, 1, 0, 0, 0)).total_seconds()/3600 
           # 
           num_of_nan=0
           if str(ssh_obs_arr[otime]) != 'nan' : # rm nan values from arrays
              time_obs.append(time_obs_el)
              ssh_obs.append(float(ssh_obs_arr[otime]))
           else:
              num_of_nan=num_of_nan+1
   
       print ('# of nan = ',num_of_nan)
       print ("Fit mod...")
       fitted_mod, cov_mod = curve_fit(octuple_mod,time_mod[ts_mod:te_mod],ssh_mod[ts_mod:te_mod])
       print ("Fit obs...")
       fitted_obs, cov_obs = curve_fit(octuple_obs,time_obs[ts_obs:te_obs],ssh_obs[ts_obs:te_obs])
   
   
       txtname=path+"pre_imm"+stations_mod[stn]+'_'+dates_lab+".txt"
   
       # M2 mod
       if fitted_mod[0] < 0:
           print('Am2_mod<0')
           fitted_mod[0] = -fitted_mod[0]
           fitted_mod[1] = fitted_mod[1]+180
       M2_amp_mod.append(fitted_mod[0]*M2ft)
   
       pha_mod = fitted_mod[1]+M2uvt
       if  pha_mod > 360:
           print('Pm2_mod>360')
           pha_mod=pha_mod-360
   
       if  pha_mod > 360:
           print('Pm2_mod>360')
           pha_mod=pha_mod-360
   
       elif pha_mod < -360:
           pha_mod=pha_mod+360
   
       if pha_mod < 0:
            print('Pm2_mod<0 ==> mod: ')
            pha_mod = pha_mod+360
            print (pha_mod)
   
       fitted_mod[1] = pha_mod
       M2_pha_mod.append(pha_mod)
   
       # K1 mod    
       if fitted_mod[2] < 0:
           print('Ak1_mod<0')
           fitted_mod[2] = - fitted_mod[2]
           fitted_mod[3] = fitted_mod[3] + 180
       K1_amp_mod.append(fitted_mod[2]*K1ft)
       pha_mod = fitted_mod[3] + K1uvt
       if  pha_mod > 360:
           print('Pk1_mod>360')
           pha_mod = pha_mod-360
   
       if pha_mod < 0:
           print('Pk1_mod<0')
           pha_mod = pha_mod+360
   
       if pha_mod < 0:
           print('Pk1_mod<0')
           pha_mod = pha_mod+360
   
       fitted_mod[3] = pha_mod
       K1_pha_mod.append(pha_mod)
       
       # O1 mod
       if fitted_mod[4] < 0:
           print('Ao1_mod<0')
           fitted_mod[4] = -fitted_mod[4]
           fitted_mod[5] = fitted_mod[5]+180
       O1_amp_mod.append(fitted_mod[4]*O1ft)
       pha_mod= fitted_mod[5]+O1uvt
       if  pha_mod > 360:
           print('Po1_mod>360')
           pha_mod=pha_mod-360
       elif pha_mod < 0:
           print('Po1_mod<0')
           pha_mod = pha_mod+360
   
       fitted_mod[5] = pha_mod
       O1_pha_mod.append(pha_mod)
       
       # S2 mod
       if fitted_mod[6] < 0:
           print('As2_mod<0')
           fitted_mod[6] = -fitted_mod[6]
           fitted_mod[7] = fitted_mod[7]+180
       S2_amp_mod.append(fitted_mod[6]*S2ft)
       pha_mod= fitted_mod[7]+S2uvt
   
       if  pha_mod > 360:
           print('Ps2_mod>360')
           pha_mod=pha_mod-360
       elif pha_mod < -360:
           pha_mod=pha_mod+360
   
       if pha_mod < 0:
           print('Ps2_mod<0')
           pha_mod = pha_mod+360
   
       fitted_mod[7] = pha_mod
       S2_pha_mod.append(pha_mod)
       
       
       # P1 mod
       if fitted_mod[8] < 0:
               fitted_mod[8] = -fitted_mod[8]
               fitted_mod[9] = fitted_mod[9]+180
       P1_amp_mod.append(fitted_mod[8]*P1ft)
       pha_mod= fitted_mod[9]+P1uvt
       if  pha_mod > 360:
           pha_mod=pha_mod-360
       elif pha_mod < 0:
           pha_mod = pha_mod+360
       if  pha_mod > 320:
           pha_mod=pha_mod-360
   
       fitted_mod[9] = pha_mod
       P1_pha_mod.append(pha_mod)
   
       # N2 mod
       if fitted_mod[10] < 0:
               fitted_mod[10] = -fitted_mod[10]
               fitted_mod[11] = fitted_mod[11]+180
       N2_amp_mod.append(fitted_mod[10]*N2ft)
       pha_mod= fitted_mod[11]+N2uvt
   
       if pha_mod > 720:
          pha_mod = pha_mod-720
       elif pha_mod > 360:
           pha_mod = pha_mod-360
       elif pha_mod < -360:
           pha_mod = pha_mod+360
   
       if pha_mod > 360:
           pha_mod = pha_mod-360
   
       if pha_mod < 0:
          pha_mod = pha_mod+360
   
       fitted_mod[11] = pha_mod
       N2_pha_mod.append(pha_mod)
       
        # Q1 mod
       if fitted_mod[12] < 0:
               fitted_mod[12] = -fitted_mod[12]
               fitted_mod[13] = fitted_mod[13]+180
       Q1_amp_mod.append(fitted_mod[12]*Q1ft)
       pha_mod= fitted_mod[13]+Q1uvt
       if  pha_mod > 360:
           pha_mod=pha_mod-360
       elif pha_mod < -360:
          pha_mod=pha_mod+360
   
       if pha_mod < 0:
           pha_mod = pha_mod+360
       fitted_mod[13] = pha_mod
       Q1_pha_mod.append(pha_mod)
       
       # K2 mod
       if fitted_mod[14] < 0:
               fitted_mod[14] = -fitted_mod[14]
               fitted_mod[15] = fitted_mod[15]+180
       K2_amp_mod.append(fitted_mod[14]*K2ft)
       pha_mod= fitted_mod[15]+K2uvt
       if pha_mod > 720:
          pha_mod = pha_mod-720
       elif pha_mod > 360:
           pha_mod = pha_mod-360
       elif pha_mod < -360:
           pha_mod = pha_mod+360
   
       if pha_mod < 0:
           pha_mod = pha_mod+360
       fitted_mod[15] = pha_mod
       K2_pha_mod.append(pha_mod)
       
       print('##########FIT MOD')
       print(fitted_mod[:])
   
   ##################
   ### OBS ###
   ##################
       if fitted_obs[0] < 0:
           print('Am2_obs<0')
           fitted_obs[0] = -fitted_obs[0]
           fitted_obs[1] = fitted_obs[1]+180
       M2_amp_obs.append(fitted_obs[0]*M2ft)
   
       pha_obs = fitted_obs[1]+M2uvt
       if  pha_obs > 360:
           print('Pm2_obs>360')
           pha_obs=pha_obs-360
   
       if  pha_obs > 360:
           print('Pm2_obs>360')
           pha_obs=pha_obs-360
   
       elif pha_obs < -360:
           pha_obs=pha_obs+360
   
       if pha_obs < 0:
            print('Pm2_obs<0 ==> mod: ')
            pha_obs = pha_obs+360
            print (pha_obs)
   
       fitted_obs[1] = pha_obs
       M2_pha_obs.append(pha_obs)
   
       if fitted_obs[2] < 0:
           print('Ak1_obs<0')
           fitted_obs[2] = - fitted_obs[2]
           fitted_obs[3] = fitted_obs[3] + 180
       K1_amp_obs.append(fitted_obs[2]*K1ft)
       pha_obs = fitted_obs[3] + K1uvt
       if  pha_obs > 360:
           print('Pk1_obs>360')
           pha_obs = pha_obs-360
   
       if pha_obs < 0:
           print('Pk1_obs<0')
           pha_obs = pha_obs+360
   
       if pha_obs < 0:
           print('Pk1_obs<0')
           pha_obs = pha_obs+360
   
       fitted_obs[3] = pha_obs
       K1_pha_obs.append(pha_obs)
   
       if fitted_obs[4] < 0:
           print('Ao1_obs<0')
           fitted_obs[4] = -fitted_obs[4]
           fitted_obs[5] = fitted_obs[5]+180
       O1_amp_obs.append(fitted_obs[4]*O1ft)
       pha_obs= fitted_obs[5]+O1uvt
       if  pha_obs > 360:
           print('Po1_obs>360')
           pha_obs=pha_obs-360
       elif pha_obs < 0:
           print('Po1_obs<0')
           pha_obs = pha_obs+360
   
       fitted_obs[5] = pha_obs
       O1_pha_obs.append(pha_obs)
   
       if fitted_obs[6] < 0:
           print('As2_obs<0')
           fitted_obs[6] = -fitted_obs[6]
           fitted_obs[7] = fitted_obs[7]+180
       S2_amp_obs.append(fitted_obs[6]*S2ft)
       pha_obs= fitted_obs[7]+S2uvt
   
       if  pha_obs > 360:
           print('Ps2_obs>360')
           pha_obs=pha_obs-360
       if pha_obs < -360:
           pha_obs=pha_obs+360
       if pha_obs < 0:
           print('Ps2_obs<0')
           pha_obs = pha_obs+360
       if pha_obs < 0:
           print('Ps2_obs<0')
           pha_obs = pha_obs+360
       fitted_obs[7] = pha_obs
       S2_pha_obs.append(pha_obs)
   
   
       #####
       if fitted_obs[8] < 0:
               fitted_obs[8] = -fitted_obs[8]
               fitted_obs[9] = fitted_obs[9]+180
       P1_amp_obs.append(fitted_obs[8]*P1ft)
       pha_obs= fitted_obs[9]+P1uvt
       if  pha_obs > 360:
           pha_obs=pha_obs-360
       elif pha_obs < 0:
           pha_obs = pha_obs+360
       fitted_obs[9] = pha_obs
       P1_pha_obs.append(pha_obs)
   
       if fitted_obs[10] < 0:
               fitted_obs[10] = -fitted_obs[10]
               fitted_obs[11] = fitted_obs[11]+180
       N2_amp_obs.append(fitted_obs[10]*N2ft)
       pha_obs= fitted_obs[11]+N2uvt
   
       if pha_obs > 720:
          pha_obs = pha_obs-720
       elif pha_obs > 360:
           pha_obs = pha_obs-360
       elif pha_obs < -360:
           pha_obs = pha_obs+360
   
       if pha_obs > 360:
           pha_obs = pha_obs-360
   
       if pha_obs < 0:
          pha_obs = pha_obs+360
       
   
       fitted_obs[11] = pha_obs
       N2_pha_obs.append(pha_obs)
   
       if fitted_obs[12] < 0:
               fitted_obs[12] = -fitted_obs[12]
               fitted_obs[13] = fitted_obs[13]+180
       Q1_amp_obs.append(fitted_obs[12]*Q1ft)
       pha_obs= fitted_obs[13]+Q1uvt
       if  pha_obs > 360:
           pha_obs=pha_obs-360
       elif pha_obs < -360:
          pha_obs=pha_obs+360
   
       if pha_obs < 0:
           pha_obs = pha_obs+360
       fitted_obs[13] = pha_obs
       Q1_pha_obs.append(pha_obs)
   
       if fitted_obs[14] < 0:
               fitted_obs[14] = -fitted_obs[14]
               fitted_obs[15] = fitted_obs[15]+180
       K2_amp_obs.append(fitted_obs[14]*K2ft)
       pha_obs= fitted_obs[15]+K2uvt
       if pha_obs > 720:
          pha_obs = pha_obs-720
       elif pha_obs > 360:
           pha_obs = pha_obs-360
       elif pha_obs < -360:
           pha_obs = pha_obs+360
   
       if pha_obs < 0:
           pha_obs = pha_obs+360
       fitted_obs[15] = pha_obs
       K2_pha_obs.append(pha_obs)
   
       print('##########FIT OBS')
       print(fitted_obs[:])
   
   
       txtname=path+"imm"+stations_mod[stn]+'_'+dates_lab+".txt"
   
   tidal_c=['M2','N2','K2','S2','K1','O1','Q1','P1']
   
   for comp in ('M2','N2','K2','S2','K1','O1','Q1','P1'):
   
       ###
       nameA_obs=comp+'_amp_obs'
       nameA_mod=comp+'_amp_mod'
   
       nameP_obs=comp+'_pha_obs'
       nameP_mod=comp+'_pha_mod'
   
   
   
       # Save ISPRA arrays
      
       IAmod_name='I_'+nameA_mod
       IAobs_name='I_'+nameA_obs
   
       globals()[IAmod_name]=np.array(globals()[nameA_mod])
       globals()[IAobs_name]=np.array(globals()[nameA_obs])
   
       I_stations_lab=stations_lab
       I_numsta=numsta_obs
       I_stations_col=stations_col
   
       IPmod_name='I_'+nameP_mod
       IPobs_name='I_'+nameP_obs
   
   
       globals()[IPmod_name]=np.array(globals()[nameP_mod])
       globals()[IPobs_name]=np.array(globals()[nameP_obs])
  
       # MIX case
       if ttide_flag == 2:
         if flag_15stats == 0:
          I_M2_amp_mod=[5.8198,6.1121,5.9116,7.4113,5.8557,8.5442,3.0175,6.3895,11.5321,4.0252,6.8800,5.5154,24.8351,21.9406,8.7600]
          I_S2_amp_mod=[3.1273,2.3397,2.9233,2.6831,3.2933,3.1613,1.7678,4.5953,4.5990,2.2132,2.5750,2.7195,15.4861,13.6757,5.6851]
          I_K1_amp_mod=[12.6234,3.6305,1.8898,3.4995,0.7324,3.5033,1.0566,8.7689,2.7299,1.7052,3.5556,1.6563,17.3597,16.7605,5.0634]
          I_O1_amp_mod=[3.9300,1.6296,1.0153,1.6909,0.6475,1.5429,0.5940,2.8863,0.8668,1.1074,1.5922,0.9394,5.0544,4.9113,1.8294]
          I_N2_amp_mod=[1.0118,1.2407,1.0770,1.5522,0.8134,1.7737,0.5206,1.0132,2.3305,0.7642,1.3942,1.0052,4.2636,3.7479,1.4630]
          I_P1_amp_mod=[3.7570,1.1749,0.6204,1.0377,0.1673,1.0540,0.2552,2.5885,0.9421,0.3316,1.1248,0.5453,4.7057,4.8115,1.6146]
          I_Q1_amp_mod=[1.0211,0.2793,0.1945,0.2725,0.1953,0.2372,0.1892,0.6233,0.1484,0.1488,0.2778,0.2200,1.4955,1.4418,0.2975]
          I_K2_amp_mod=[0.9966,0.6661,0.8556,0.7728,1.0068,0.9107,0.4295,1.2103,1.2491,0.6993,0.7838,0.7660,3.8812,3.5495,1.5315]
        
          I_M2_pha_mod=[296.66,230.04,56.01,222.38,42.10,227.54,51.33,63.30,231.18,78.78,225.69,61.68,249.22,257.45,69.35]
          I_S2_pha_mod=[312.14,250.12,62.88,238.39,50.05,244.39,55.09,75.71,250.10,72.74,244.79,67.79,258.94,266.32,79.96]
          I_K1_pha_mod=[68.54,181.78,37.17,180.76,305.15,179.72,15.15,65.77,198.99,83.37,183.49,30.13,57.12,61.56,70.35]
          I_O1_pha_mod=[56.85,113.04,41.56,108.40,94.59,113.95,22.94,54.43,122.05,67.71,110.07,40.74,48.52,51.46,57.72]
          I_N2_pha_mod=[291.95,215.86,53.85,210.33,49.26,215.66,60.74,61.94,217.62,96.32,212.43,61.14,248.05,257.49,67.49]
          I_P1_pha_mod=[61.54,183.16,23.38,187.43,343.03,192.12,27.68,59.06,207.00,71.38,188.65,39.82,44.41,51.82,64.61]
          I_Q1_pha_mod=[94.98,50.81,42.47,41.28,35.24,60.41,18.65,92.57,1.21,57.55,42.06,37.04,76.18,80.88,88.37]
          I_K2_pha_mod=[302.50,244.50,56.00,235.01,48.92,239.71,47.84,59.19,246.27,66.21,240.80,54.06,243.76,253.96,66.90]

         elif flag_15stats == 1:
 
          I_M2_amp_mod=[5.8198,6.1121,5.9116,5.8557,8.5442,6.3895,4.0252,5.5154,24.8351,21.9406,8.7600]
          I_S2_amp_mod=[3.1273,2.3397,2.9233,3.2933,3.1613,4.5953,2.2132,2.7195,15.4861,13.6757,5.6851]
          I_K1_amp_mod=[12.6234,3.6305,1.8898,0.7324,3.5033,8.7689,1.7052,1.6563,17.3597,16.7605,5.0634]
          I_O1_amp_mod=[3.9300,1.6296,1.0153,0.6475,1.5429,2.8863,1.1074,0.9394,5.0544,4.9113,1.8294]
          I_N2_amp_mod=[1.0118,1.2407,1.0770,0.8134,1.7737,1.0132,0.7642,1.0052,4.2636,3.7479,1.4630]
          I_P1_amp_mod=[3.7570,1.1749,0.6204,0.1673,1.0540,2.5885,0.3316,0.5453,4.7057,4.8115,1.6146]
          I_Q1_amp_mod=[1.0211,0.2793,0.1945,0.1953,0.2372,0.6233,0.1488,0.2200,1.4955,1.4418,0.2975]
          I_K2_amp_mod=[0.9966,0.6661,0.8556,1.0068,0.9107,1.2103,0.6993,0.7660,3.8812,3.5495,1.5315]
         
          I_M2_pha_mod=[296.66,230.04,56.01,42.10,227.54,63.30,78.78,61.68,249.22,257.45,69.35]
          I_S2_pha_mod=[312.14,250.12,62.88,50.05,244.39,75.71,72.74,67.79,258.94,266.32,79.96]
          I_K1_pha_mod=[68.54,181.78,37.17,305.15,179.72,65.77,83.37,30.13,57.12,61.56,70.35]
          I_O1_pha_mod=[56.85,113.04,41.56,94.59,113.95,54.43,67.71,40.74,48.52,51.46,57.72]
          I_N2_pha_mod=[291.95,215.86,53.85,49.26,215.66,61.94,96.32,61.14,248.05,257.49,67.49]
          I_P1_pha_mod=[61.54,183.16,23.38,343.03,192.12,59.06,71.38,39.82,44.41,51.82,64.61]
          I_Q1_pha_mod=[94.98,50.81,42.47,35.24,60.41,92.57,57.55,37.04,76.18,80.88,88.37]
          I_K2_pha_mod=[302.50,244.50,56.00,48.92,239.71,59.19,66.21,54.06,243.76,253.96,66.90]

# ======== EMODnet ========== 

if where_box == 'Med':
   
   if flag_15stats == 1:
      print ('Lit comparison..')
      # For comparison with Tsimplis et al. 1995 and Palma et al 2020
      stations_mod = ['Almeria','Malaga','Marseille','Tarifa']
      stations_obs = stations_mod
      stations_lab = ['Almeria','Malaga','Marseille','Tarifa']
      stations_col = ['green','red','green','green']
   else:
      # STZ EMODNET OK 2017 HOURLY and RENAMED + EAST MED STZ + NEW MESS
      stations_mod = ['Ajaccio','Algeciras','Almeria','Barcelona','Centuri','Ginostra','Ibiza','IleRousse','iske','LaFigueirette','Malaga','Marseille','Melilla','Monaco','Motril','PalmadeMallorca','PortLaNouvelle','PortVendres','Sagunto','Sete','Solenzara','Tarifa','Valencia','zygi1']
      stations_obs = ['Ajaccio','Algeciras','Almeria','Barcelona','Centuri','Ginostra','Ibiza','IleRousse','iske','LaFigueirette','Malaga','Marseille','Melilla','Monaco','Motril','PalmadeMallorca','PortLaNouvelle','PortVendres','Sagunto','Sete','Solenzara','Tarifa','Valencia','zygi1']
      stations_lab = ['Ajaccio','Algeciras','Almeria','Barcelona','Centuri','Ginostra','Ibiza','IleRousse','iske','LaFigueirette','Malaga','Marseille','Melilla','Monaco','Motril','P.deMallorca','P.LaNouvelle','P.Vendres','Sagunto','Sete','Solenzara','Tarifa','Valencia','zygi1']
      stations_col = ['green','red','green','green','green','green','green','green','magenta','green','red','green','red','green','red','green','green','green','green','green','green','red','green','magenta']
 
   ##
if where_box == 'AtlBox':

   # ATLANTIC BOX OK - ARC - B.P
   stations_mod = ['AngletConvergent','BayonneBoucau','BayonneQuaiDeLesseps','Bilbao','Bonanza','Ciboure','Coruna','Ferrol2','Ferrol','Gijon','Huelva','Langosteira','Leixoes','Marin','Nazare','Peniche','Santander','Sines','Socoa','Vigo','Villagarcia']
   stations_obs = ['AngletConvergent','BayonneBoucau','BayonneQuaiDeLesseps','Bilbao','Bonanza','Ciboure','Coruna','Ferrol2','Ferrol','Gijon','Huelva','Langosteira','Leixoes','Marin','Nazare','Peniche','Santander','Sines','Socoa','Vigo','Villagarcia']
   stations_lab = ['AngletConvergent','B.Boucau','B.QuaiDeLesseps','Bilbao','Bonanza','Ciboure','Coruna','Fe2rrol','Ferrol','Gijon','Huelva','Langosteira','Leixoes','Marin','Nazare','Peniche','Santander','Sines','Socoa','Vigo','V.garcia']
   stations_col = ['cyan','cyan','cyan','cyan','deeppink','cyan','tab:olive','tab:olive','tab:olive','cyan','deeppink','tab:olive','tab:olive','tab:olive','tab:olive','tab:olive','cyan','tab:olive','cyan','tab:olive','tab:olive']


numsta_obs=len(stations_obs)
numsta_mod=len(stations_mod)

print("EmodNET Stations num:",len(stations_obs),len(stations_mod) )


##
##
# TTIDE METHODOLOGY
if ttide_flag == 1:

 numsta=numsta_obs

 #%comp%_amp/pha_mod/obs=[]

 if flag_15stats == 0:

  M2_amp_mod=[7.0103,28.9092,8.6594,4.3750,7.7169,11.3976,1.7285,7.4710,10.4221,6.9434,16.9520,6.2713,11.0818,7.0897,14.0424,2.6202,5.8046,5.5151,1.9854,6.1165,9.8293,35.7517,1.9118,8.2126]
  S2_amp_mod=[2.5883,11.0820,3.6346,1.3329,2.8740,4.5656,0.3778,2.7701,5.2174,2.4843,6.7628,2.1442,4.4477,2.5680,5.7389,0.7514,1.8585,1.7417,0.3261,1.9916,3.6701,13.7243,0.3147,4.1841]
  K1_amp_mod=[3.4281,2.6905,3.7932,3.5498,3.5426,2.8028,3.6944,3.4931,2.9552,3.5328,3.3829,3.5554,4.0055,3.5134,3.4551,3.5835,3.5095,3.5583,4.0131,3.5682,2.4663,2.3121,4.0217,2.4203]
  O1_amp_mod=[1.5264,1.0258,2.2624,2.1916,1.5482,0.8286,2.1499,1.5743,2.2651,1.7817,2.0658,1.8895,2.1781,1.7547,2.0377,2.0676,2.0541,2.0518,2.4673,1.9678,1.3360,0.5640,2.4617,1.8698]
  N2_amp_mod=[1.4298,6.0143,1.7139,0.9657,1.5907,2.2883,0.4147,1.5241,1.7481,1.4446,3.3983,1.3409,2.2455,1.4851,2.8533,0.5797,1.2677,1.2086,0.5220,1.3354,2.0463,7.6738,0.5109,1.3800]
  P1_amp_mod=[1.1125,0.7606,0.9898,1.0428,1.1360,1.0159,1.1241,1.1336,1.2444,1.0502,0.9461,1.0612,1.1095,1.0411,1.0197,1.0272,0.7452,0.9095,0.9163,0.8181,0.8803,0.7404,0.9422,0.9659]
  Q1_amp_mod=[0.2648,0.4300,0.2762,0.3006,0.2045,0.1701,0.2570,0.2397,0.3682,0.3012,0.3013,0.3198,0.1608,0.3067,0.2486,0.3141,0.3528,0.3319,0.3836,0.3546,0.2631,0.5689,0.3554,0.2596]
  K2_amp_mod=[0.7690,3.0959,1.0769,0.4023,0.8607,1.1749,0.1129,0.8018,1.6381,0.6961,1.9309,0.6185,1.1977,0.7098,1.6231,0.2102,0.5451,0.4984,0.0962,0.5836,1.0170,3.9425,0.0892,1.2774]
  
  M2_pha_mod=[224.48,46.85,50.99,218.98,223.84,333.68,225.47,223.26,232.33,221.26,49.04,219.45,68.36,221.38,46.80,214.39,217.84,218.76,215.43,217.38,227.84,45.35,215.04,239.98]
  S2_pha_mod=[243.07,72.20,73.03,239.51,240.65,252.65,249.40,240.69,240.35,238.99,72.93,237.56,86.98,237.34,71.89,235.04,237.50,239.90,234.21,236.44,249.58,72.16,229.31,246.72]
  K1_pha_mod=[183.43,129.84,157.05,166.88,187.17,198.10,167.05,185.36,275.60,178.35,147.80,173.09,148.78,180.33,151.87,166.86,165.21,166.15,162.00,164.59,182.69,118.14,162.88,280.23]
  O1_pha_mod=[110.82,137.09,114.90,104.65,113.23,227.57,107.06,112.10,258.97,105.93,118.71,103.32,98.81,107.59,115.61,104.11,103.64,103.44,106.22,102.86,107.10,175.89,106.16,257.60]
  N2_pha_mod=[212.02,30.39,34.51,204.11,211.15,230.85,205.95,210.36,231.98,209.45,32.74,207.23,55.16,210.13,29.81,200.02,205.50,206.13,199.62,205.90,213.15,28.92,199.28,240.55]
  P1_pha_mod=[190.01,125.93,148.35,171.02,187.79,201.20,162.22,188.35,282.93,186.03,144.81,178.01,146.10,183.50,145.92,169.10,171.55,172.76,158.89,180.06,190.69,125.34,158.23,281.80]
  Q1_pha_mod=[38.52,181.41,92.14,47.77,46.37,26.63,54.85,36.46,244.69,52.07,130.23,41.02,72.48,47.60,117.41,47.02,38.17,39.84,78.32,31.83,48.82,197.68,77.97,227.69]
  K2_pha_mod=[237.85,69.95,69.52,234.34,237.23,244.29,252.86,236.33,235.85,234.48,68.84,232.80,82.69,231.78,68.08,234.23,233.84,233.14,233.50,232.17,244.94,67.83,237.59,240.99]
  
  M2_amp_obs=[7.4715,30.3436,9.0690,4.3331,8.3190,10.8910,1.7189,7.8536,12.3054,7.5119,17.9133,5.2311,11.6687,7.7943,15.3551,2.5371,6.1381,5.7004,1.6583,6.2630,10.2843,39.4856,1.6562,8.4823]
  S2_amp_obs=[2.9420,10.3727,3.6266,1.5767,3.2268,4.3456,0.6943,3.1823,7.1780,2.9008,6.9659,1.9811,4.3057,3.0450,5.8578,0.9566,2.2019,2.0782,0.5125,2.3178,3.8503,15.1493,0.5282,4.7878]
  K1_amp_obs=[3.5971,2.3674,3.4971,3.6526,3.3993,2.4437,3.8159,3.3989,2.2235,3.5997,3.1557,3.0829,2.9969,3.7168,3.0848,3.7255,3.7955,3.7787,3.8894,3.6675,2.5581,2.5552,3.8908,2.1071]
  O1_amp_obs=[1.4417,0.6732,1.8663,2.0219,1.3613,0.8672,2.0264,1.5679,1.5454,1.7053,1.9252,1.3299,1.7885,1.7053,1.7004,1.8435,2.0104,1.9205,2.2537,1.7521,1.2808,0.8039,2.2029,2.1458]
  N2_amp_obs=[1.6153,6.1808,1.7053,0.9125,1.7089,2.7568,0.4058,1.7619,2.3380,1.5555,3.4832,0.8766,2.2295,1.5930,3.1253,0.4236,1.2477,1.2114,0.3412,1.4094,2.2743,8.6489,0.3984,1.2875]
  P1_amp_obs=[1.0621,0.7671,1.0994,1.2512,0.9408,0.6735,1.3297,1.0771,0.7758,0.9756,1.0421,1.0179,0.7732,0.9755,1.1975,1.1974,0.6872,0.6792,0.7284,0.7189,0.7260,1.3385,0.8917,0.6880]
  Q1_amp_obs=[0.2833,0.2814,0.2922,0.2913,0.3450,0.3207,0.3751,0.2534,0.6675,0.3554,0.4505,0.2582,0.2494,0.2947,0.2277,0.2077,0.3606,0.4035,0.2791,0.3752,0.4061,0.7376,0.2097,0.3954]
  K2_amp_obs=[0.7878,2.6082,1.3188,0.4938,1.0635,1.1266,0.1722,1.1446,2.6180,0.8707,1.8456,1.1672,1.6494,0.9065,1.8958,0.1887,0.6125,0.5263,0.1398,0.6817,1.0555,4.1837,0.1458,1.6327]
  
  M2_pha_obs=[224.30,47.41,51.59,211.75,222.15,233.80,212.34,221.26,229.06,218.45,48.97,219.63,65.38,220.71,45.65,207.68,212.50,214.73,196.44,213.55,227.78,41.31,195.19,237.88]
  S2_pha_obs=[240.16,71.36,72.47,226.50,240.43,252.43,228.48,241.80,242.54,236.90,71.09,237.88,88.15,235.87,71.68,221.32,229.21,234.59,199.67,230.09,248.53,65.33,193.72,248.78]
  K1_pha_obs=[181.30,120.56,159.87,170.68,184.82,201.80,171.76,187.63,270.66,175.81,156.29,168.83,143.02,181.27,158.45,168.96,167.22,165.57,162.63,166.41,186.05,117.61,167.28,264.81]
  O1_pha_obs=[105.30,155.17,119.16,98.35,109.43,123.08,104.01,109.37,254.48,100.88,120.25,96.27,128.15,100.96,117.30,100.22,95.96,97.19,105.51,94.68,107.68,24.24,102.78,245.40]
  N2_pha_obs=[213.60,21.98,35.82,200.95,212.77,212.77,197.46,217.50,239.00,205.47,34.34,218.43,52.25,209.97,33.18,189.66,200.52,204.14,163.96,201.86,213.90,27.43,181.14,241.64]
  P1_pha_obs=[193.70,128.56,135.47,160.14,177.95,189.16,154.82,198.99,280.28,192.07,129.20,216.56,158.36,181.07,131.90,169.45,174.89,179.50,171.11,177.88,186.24,128.81,145.38,314.95]
  Q1_pha_obs=[56.77,174.82,95.18,67.87,92.84,28.87,80.53,80.05,286.84,67.59,146.16,100.60,100.62,77.21,129.58,56.75,64.81,66.90,103.41,39.89,66.82,194.19,74.89,250.60]
  K2_pha_obs=[232.14,49.15,69.25,228.32,234.32,285.89,206.41,218.33,242.74,234.21,73.42,273.10,109.24,235.96,65.04,215.07,228.19,227.29,208.70,235.03,248.56,62.08,223.07,257.59]

 elif flag_15stats == 1:

  M2_amp_mod=[8.6594,16.9520,6.2713,35.7517]
  S2_amp_mod=[3.6346,6.7628,2.1442,13.7243]
  K1_amp_mod=[3.7932,3.3829,3.5554,2.3121]
  O1_amp_mod=[2.2624,2.0658,1.8895,0.5640]
  N2_amp_mod=[1.7139,3.3983,1.3409,7.6738]
  P1_amp_mod=[0.9898,0.9461,1.0612,0.7404]
  Q1_amp_mod=[0.2762,0.3013,0.3198,0.5689]
  K2_amp_mod=[1.0769,1.9309,0.6185,3.9425]
  
  M2_pha_mod=[50.99,49.04,219.45,45.35]
  S2_pha_mod=[73.03,72.93,237.56,72.16]
  K1_pha_mod=[157.05,147.80,173.09,118.14]
  O1_pha_mod=[114.90,118.71,103.32,175.89]
  N2_pha_mod=[34.51,32.74,207.23,28.92]
  P1_pha_mod=[148.35,144.81,178.01,125.34]
  Q1_pha_mod=[92.14,130.23,41.02,197.68]
  K2_pha_mod=[69.52,68.84,232.80,67.83]
  
  M2_amp_obs=[9.0690,17.9133,5.2311,39.4856]
  S2_amp_obs=[3.6266,6.9659,1.9811,15.1493]
  K1_amp_obs=[3.4971,3.1557,3.0829,2.5552]
  O1_amp_obs=[1.8663,1.9252,1.3299,0.8039]
  N2_amp_obs=[1.7053,3.4832,0.8766,8.6489]
  P1_amp_obs=[1.0994,1.0421,1.0179,1.3385]
  Q1_amp_obs=[0.2922,0.4505,0.2582,0.7376]
  K2_amp_obs=[1.3188,1.8456,1.1672,4.1837]
  
  M2_pha_obs=[50.99,49.04,219.45,45.35]
  S2_pha_obs=[73.03,72.93,237.56,72.16]
  K1_pha_obs=[157.05,147.80,173.09,118.14]
  O1_pha_obs=[114.90,118.71,103.32,175.89]
  N2_pha_obs=[34.51,32.74,207.23,28.92]
  P1_pha_obs=[148.35,144.81,178.01,125.34]
  Q1_pha_obs=[92.14,130.23,41.02,197.68]
  K2_pha_obs=[69.52,68.84,232.80,67.83]

# SALISH SEA METHODOLOGY or MIX case
else:
   
   #constants and fitting
   # M2
   M2freq = 28.984106 # degrees per hour
   M2freq = M2freq*np.pi/180. # radians per hour
   ###K1
   K1freq = 15.041069*np.pi/180.
   ###O1
   O1freq = 13.943036*np.pi/180.
   ###S2
   S2freq = 30.000002*np.pi/180.
   ###P1
   P1freq = 14.958932*np.pi/180.
   ###N2
   N2freq = 28.439730*np.pi/180.
   ###Q1
   Q1freq = 13.398661*np.pi/180.
   ###K2
   K2freq = 30.082138*np.pi/180.
   ##
   
   ### initial phase calculation
   
   # M2
   M2ft = 1.03669420605307
   M2uvt = -12.5753366566839
   
   # S2
   S2ft = 1.000000
   S2uvt = 0.000000
   
   # N2
   N2ft = 1.0
   N2uvt = 0.0
   
   # K1
   K1ft = 0.886345136935456
   K1uvt = -0.4407793529972482
   
   # O1
   O1ft = 0.813447400642605
   O1uvt = -12.5067505117739
   
   # Q1
   Q1ft = 1.0
   Q1uvt = 0.0
   
   # K2
   K2ft = 1.0
   K2uvt = 0.0
   
   # P1
   P1ft = 1.0
   P1uvt = 0.0
   
   # Add nodal corrections and phase at t=0 (01072017)
   def octuple_obs(x, M2amp_obs, M2pha_obs, K1amp_obs, K1pha_obs, O1amp_obs, O1pha_obs, S2amp_obs, S2pha_obs, P1amp_obs, P1pha_obs, N2amp_obs, N2pha_obs, Q1amp_obs, Q1pha_obs, K2amp_obs, K2pha_obs):
       return (A_nc[0]*M2amp_obs*np.cos(M2freq*x-M2pha_obs*np.pi/180.+P_shift[0])+
               A_nc[2]*K1amp_obs*np.cos(K1freq*x-K1pha_obs*np.pi/180.+P_shift[2])+
               A_nc[3]*O1amp_obs*np.cos(O1freq*x-O1pha_obs*np.pi/180.+P_shift[3])+
               A_nc[1]*S2amp_obs*np.cos(S2freq*x-S2pha_obs*np.pi/180.+P_shift[1])+
               A_nc[5]*P1amp_obs*np.cos(P1freq*x-P1pha_obs*np.pi/180.+P_shift[5])+
               A_nc[4]*N2amp_obs*np.cos(N2freq*x-N2pha_obs*np.pi/180.+P_shift[4])+
               A_nc[6]*Q1amp_obs*np.cos(Q1freq*x-Q1pha_obs*np.pi/180.+P_shift[6])+
               A_nc[7]*K2amp_obs*np.cos(K2freq*x-K2pha_obs*np.pi/180.+P_shift[7]))
   
   def octuple_mod(x, M2amp_mod, M2pha_mod, K1amp_mod, K1pha_mod, O1amp_mod, O1pha_mod, S2amp_mod, S2pha_mod, P1amp_mod, P1pha_mod, N2amp_mod, N2pha_mod, Q1amp_mod, Q1pha_mod, K2amp_mod, K2pha_mod):
       return (A_nc[0]*M2amp_mod*np.cos(M2freq*x-M2pha_mod*np.pi/180.+P_shift[0])+
               A_nc[2]*K1amp_mod*np.cos(K1freq*x-K1pha_mod*np.pi/180.+P_shift[2])+
               A_nc[3]*O1amp_mod*np.cos(O1freq*x-O1pha_mod*np.pi/180.+P_shift[3])+
               A_nc[1]*S2amp_mod*np.cos(S2freq*x-S2pha_mod*np.pi/180.+P_shift[1])+
               A_nc[5]*P1amp_mod*np.cos(P1freq*x-P1pha_mod*np.pi/180.+P_shift[5])+
               A_nc[4]*N2amp_mod*np.cos(N2freq*x-N2pha_mod*np.pi/180.+P_shift[4])+
               A_nc[6]*Q1amp_mod*np.cos(Q1freq*x-Q1pha_mod*np.pi/180.+P_shift[6])+
               A_nc[7]*K2amp_mod*np.cos(K2freq*x-K2pha_mod*np.pi/180.+P_shift[7]))
   
   ##
   
   #allocate space for our arrays
   
   M2_amp_mod=[]; M2_pha_mod=[]; K1_amp_mod=[]; K1_pha_mod=[]
   O1_amp_mod=[]; O1_pha_mod=[]; S2_amp_mod=[]; S2_pha_mod=[]
   P1_amp_mod=[]; P1_pha_mod=[]; N2_amp_mod=[]; N2_pha_mod=[]
   Q1_amp_mod=[]; Q1_pha_mod=[]; K2_amp_mod=[]; K2_pha_mod=[]
   
   M2_amp_obs=[]; M2_pha_obs=[]; K1_amp_obs=[]; K1_pha_obs=[]
   O1_amp_obs=[]; O1_pha_obs=[]; S2_amp_obs=[]; S2_pha_obs=[]
   P1_amp_obs=[]; P1_pha_obs=[]; N2_amp_obs=[]; N2_pha_obs=[]
   Q1_amp_obs=[]; Q1_pha_obs=[]; K2_amp_obs=[]; K2_pha_obs=[]
   
   M2_Diff_A=[]
   M2_Diff_P=[]
   
   # EMODNET
   for stn in range (0,len(stations_mod)): 
       print(stations_mod[stn])
   
       fT1_mod = NC.Dataset(name+stations_mod[stn]+'TG_mod_'+mod_file_template+'.nc','r') # 'TG_mod_Tides8.nc'
       fT1_obs = NC.Dataset(name+stations_obs[stn]+'TG_obs.nc','r')
       # mod:
       time_mod = (fT1_mod.variables["time_counter"][:]/3600.)-582912.0   # want hours (not seconds) from 1/1/2015 (not 1/1/1950)
       print('time mod:',time_mod)
       ssh_mod = fT1_mod.variables["sossheig"][:,0,0]*100.0 # want cm not meters
       # TG obs
       time_obs = (fT1_obs.variables["TIME"][:]*24)-582912.0 # want hours (not days)
       print ('time obs:',time_obs)
       ssh_obs = fT1_obs.variables["sossheig"][:,0]*100.0
       print(ssh_obs)
   
   
       # times
       te_mod = ssh_mod.shape[0]
       print (te_mod)
       te_obs = ssh_obs.shape[0]
       print (te_obs)
   
       ts_mod=0
       ts_obs=0
       # Fit
       print ("Fit mod...")
       fitted_mod, cov_mod = curve_fit(octuple_mod,time_mod[ts_mod:te_mod],ssh_mod[ts_mod:te_mod])
       print ("Fit obs...")
       fitted_obs, cov_obs = curve_fit(octuple_obs,time_obs[ts_obs:te_obs],ssh_obs[ts_obs:te_obs])
   
   
   
       txtname=path+"pre_omm"+stations_mod[stn]+'_'+dates_lab+".txt"
   
   
       if fitted_mod[0] < 0:
           print('Am2_mod<0')
           fitted_mod[0] = -fitted_mod[0]
           fitted_mod[1] = fitted_mod[1]+180
       M2_amp_mod.append(fitted_mod[0]*M2ft)
   
       pha_mod = fitted_mod[1]+M2uvt
   
       while pha_mod > 360:
             pha_mod=pha_mod-360
       while pha_mod < 0:
             pha_mod = pha_mod+360
       if pha_mod <= 40:
          pha_mod = pha_mod+360
   
       fitted_mod[1] = pha_mod
       M2_pha_mod.append(pha_mod)
   
       if fitted_mod[2] < 0:
           print('Ak1_mod<0')
           fitted_mod[2] = - fitted_mod[2]
           fitted_mod[3] = fitted_mod[3] + 180
       K1_amp_mod.append(fitted_mod[2]*K1ft)
       pha_mod = fitted_mod[3] + K1uvt
   
       while pha_mod > 360:
             pha_mod=pha_mod-360
       while pha_mod < 0:
             pha_mod = pha_mod+360
       if pha_mod <= 40:
          pha_mod = pha_mod+360
   
       fitted_mod[3] = pha_mod
       K1_pha_mod.append(pha_mod)
   
       # O1 tide 8
       if fitted_mod[4] < 0:
           print('Ao1_mod<0')
           fitted_mod[4] = -fitted_mod[4]
           fitted_mod[5] = fitted_mod[5]+180
       O1_amp_mod.append(fitted_mod[4]*O1ft)
       pha_mod= fitted_mod[5]+O1uvt
   
       while pha_mod > 360:
             pha_mod=pha_mod-360
       while pha_mod < 0:
             pha_mod = pha_mod+360
       #if pha_mod <= 40:
       #   pha_mod = pha_mod+360
       if pha_mod > 360:
             pha_mod=pha_mod-360
   
       fitted_mod[5] = pha_mod
       O1_pha_mod.append(pha_mod)
   
       if fitted_mod[6] < 0:
           print('As2_mod<0')
           fitted_mod[6] = -fitted_mod[6]
           fitted_mod[7] = fitted_mod[7]+180
       S2_amp_mod.append(fitted_mod[6]*S2ft)
       pha_mod= fitted_mod[7]+S2uvt
   
       while pha_mod > 360:
             pha_mod=pha_mod-360
       while pha_mod < 0:
             pha_mod = pha_mod+360
       if pha_mod <= 40:
          pha_mod = pha_mod+360
       fitted_mod[7] = pha_mod
       S2_pha_mod.append(pha_mod)
   
   
       #####
       # P1 tide 8
       if fitted_mod[8] < 0:
               fitted_mod[8] = -fitted_mod[8]
               fitted_mod[9] = fitted_mod[9]+180
       P1_amp_mod.append(fitted_mod[8]*P1ft)
       pha_mod= fitted_mod[9]+P1uvt
   
       while pha_mod > 360:
             pha_mod=pha_mod-360
       while pha_mod < 0:
             pha_mod = pha_mod+360
       if pha_mod <= 40:
          pha_mod = pha_mod+360
   
       fitted_mod[9] = pha_mod
       P1_pha_mod.append(pha_mod)
   
       if fitted_mod[10] < 0:
               fitted_mod[10] = -fitted_mod[10]
               fitted_mod[11] = fitted_mod[11]+180
       N2_amp_mod.append(fitted_mod[10]*N2ft)
       pha_mod= fitted_mod[11]+N2uvt
   
       while pha_mod > 360:
             pha_mod=pha_mod-360
       while pha_mod < 0:
             pha_mod = pha_mod+360
       if pha_mod <= 40:
          pha_mod = pha_mod+360
   
       fitted_mod[11] = pha_mod
       N2_pha_mod.append(pha_mod)
   
       # Q1 tide 8
       if fitted_mod[12] < 0:
               fitted_mod[12] = -fitted_mod[12]
               fitted_mod[13] = fitted_mod[13]+180
       Q1_amp_mod.append(fitted_mod[12]*Q1ft)
       pha_mod= fitted_mod[13]+Q1uvt
       
       while pha_mod > 360:
             pha_mod=pha_mod-360
       while pha_mod < 0:
             pha_mod = pha_mod+360
       #if pha_mod <= 40:
       #   pha_mod = pha_mod+360
       if pha_mod > 360:
             pha_mod=pha_mod-360
   
       fitted_mod[13] = pha_mod
       Q1_pha_mod.append(pha_mod)
   
       # K2 tide 8
       if fitted_mod[14] < 0:
               fitted_mod[14] = -fitted_mod[14]
               fitted_mod[15] = fitted_mod[15]+180
       K2_amp_mod.append(fitted_mod[14]*K2ft)
       pha_mod= fitted_mod[15]+K2uvt
   
       while pha_mod > 360:
             pha_mod=pha_mod-360
       while pha_mod < 0:
             pha_mod = pha_mod+360
       if pha_mod <= 40:
          pha_mod = pha_mod+360
       if pha_mod > 360:
             pha_mod=pha_mod-360
   
       fitted_mod[15] = pha_mod
       K2_pha_mod.append(pha_mod)
   
       print('##########FIT MOD')
       print(fitted_mod[:])
   
   
       # M2 obs
       if fitted_obs[0] < 0:
           print('Am2_obs<0')
           fitted_obs[0] = -fitted_obs[0]
           fitted_obs[1] = fitted_obs[1]+180
       M2_amp_obs.append(fitted_obs[0]*M2ft)
       pha_obs = fitted_obs[1]+M2uvt
   
       while pha_obs > 360:
             pha_obs=pha_obs-360
       while pha_obs < 0:
             pha_obs = pha_obs+360
       if pha_obs <= 40:
          pha_obs = pha_obs+360
   
       fitted_obs[1] = pha_obs
       M2_pha_obs.append(pha_obs)
   
       # K1 obs
       if fitted_obs[2] < 0:
           print('Ak1_obs<0')
           fitted_obs[2] = - fitted_obs[2]
           fitted_obs[3] = fitted_obs[3] + 180
       K1_amp_obs.append(fitted_obs[2]*K1ft)
       pha_obs = fitted_obs[3] + K1uvt
   
       while pha_obs > 360:
             pha_obs=pha_obs-360
       while pha_obs < 0:
             pha_obs = pha_obs+360
       if pha_obs <= 40:
          pha_obs = pha_obs+360
   
       fitted_obs[3] = pha_obs
       K1_pha_obs.append(pha_obs)
   
       # O1 obs
       if fitted_obs[4] < 0:
           print('Ao1_obs<0')
           fitted_obs[4] = -fitted_obs[4]
           fitted_obs[5] = fitted_obs[5]+180
       O1_amp_obs.append(fitted_obs[4]*O1ft)
       pha_obs= fitted_obs[5]+O1uvt
   
       while pha_obs > 360:
             pha_obs=pha_obs-360
       while pha_obs < 0:
             pha_obs = pha_obs+360
       #if pha_obs <= 40:
       #   pha_obs = pha_obs+360
       if pha_obs > 360:
          pha_obs=pha_obs-360
   
       fitted_obs[5] = pha_obs
       O1_pha_obs.append(pha_obs)
   
       # S2 obs
       if fitted_obs[6] < 0:
           print('As2_obs<0')
           fitted_obs[6] = -fitted_obs[6]
           fitted_obs[7] = fitted_obs[7]+180
       S2_amp_obs.append(fitted_obs[6]*S2ft)
       pha_obs= fitted_obs[7]+S2uvt
   
       while pha_obs > 360:
             pha_obs=pha_obs-360
       while pha_obs < 0:
             pha_obs = pha_obs+360
       if pha_obs <= 40:
          pha_obs = pha_obs+360
   
       fitted_obs[7] = pha_obs
       S2_pha_obs.append(pha_obs)
   
   
       #####
   
       # P1 obs
       if fitted_obs[8] < 0:
               fitted_obs[8] = -fitted_obs[8]
               fitted_obs[9] = fitted_obs[9]+180
       P1_amp_obs.append(fitted_obs[8]*P1ft)
       pha_obs= fitted_obs[9]+P1uvt
   
       while pha_obs > 360:
             pha_obs=pha_obs-360
       while pha_obs < 0:
             pha_obs = pha_obs+360
       if pha_obs <= 40:
          pha_obs = pha_obs+360
   
       fitted_obs[9] = pha_obs
       P1_pha_obs.append(pha_obs)
   
       # N2 obs
       if fitted_obs[10] < 0:
               fitted_obs[10] = -fitted_obs[10]
               fitted_obs[11] = fitted_obs[11]+180
       N2_amp_obs.append(fitted_obs[10]*N2ft)
       pha_obs= fitted_obs[11]+N2uvt
   
       while pha_obs > 360:
             pha_obs=pha_obs-360
       while pha_obs < 0:
             pha_obs = pha_obs+360
       if pha_obs <= 40:
          pha_obs = pha_obs+360
   
       fitted_obs[11] = pha_obs
       N2_pha_obs.append(pha_obs)
   
   
       # Q1 obs
       if fitted_obs[12] < 0:
               fitted_obs[12] = -fitted_obs[12]
               fitted_obs[13] = fitted_obs[13]+180
       Q1_amp_obs.append(fitted_obs[12]*Q1ft)
       pha_obs= fitted_obs[13]+Q1uvt
   
       while pha_obs > 360:
             pha_obs=pha_obs-360
       while pha_obs < 0:
             pha_obs = pha_obs+360
       #if pha_obs <= 40:
       #   pha_obs = pha_obs+360
       if pha_obs > 360:
             pha_obs=pha_obs-360
   
       fitted_obs[13] = pha_obs
       Q1_pha_obs.append(pha_obs)
   
       # K2 obs
       if fitted_obs[14] < 0:
               fitted_obs[14] = -fitted_obs[14]
               fitted_obs[15] = fitted_obs[15]+180
       K2_amp_obs.append(fitted_obs[14]*K2ft)
       pha_obs= fitted_obs[15]+K2uvt
   
       while pha_obs > 360:
             pha_obs=pha_obs-360
       while pha_obs < 0:
             pha_obs = pha_obs+360
       if pha_obs <= 40:
          pha_obs = pha_obs+360
       if pha_obs > 360:
             pha_obs=pha_obs-360
   
       fitted_obs[15] = pha_obs
       K2_pha_obs.append(pha_obs)
   
       print('##########FIT OBS')
       print(fitted_obs[:])
   
   
   
       txtname=path+"omm_"+stations_mod[stn]+'_'+dates_lab+"_mod.txt"
       
       # MIX case
       if flag_15stats == 2:
         if flag_15stats == 0:
           M2_amp_mod=[7.0103,28.9092,8.6594,4.3750,7.7169,11.3976,1.7285,7.4710,10.4221,6.9434,16.9520,6.2713,11.0818,7.0897,14.0424,2.6202,5.8046,5.5151,1.9854,6.1165,9.8293,35.7517,1.9118,8.2126]
           S2_amp_mod=[2.5883,11.0820,3.6346,1.3329,2.8740,4.5656,0.3778,2.7701,5.2174,2.4843,6.7628,2.1442,4.4477,2.5680,5.7389,0.7514,1.8585,1.7417,0.3261,1.9916,3.6701,13.7243,0.3147,4.1841]
           K1_amp_mod=[3.4281,2.6905,3.7932,3.5498,3.5426,2.8028,3.6944,3.4931,2.9552,3.5328,3.3829,3.5554,4.0055,3.5134,3.4551,3.5835,3.5095,3.5583,4.0131,3.5682,2.4663,2.3121,4.0217,2.4203]
           O1_amp_mod=[1.5264,1.0258,2.2624,2.1916,1.5482,0.8286,2.1499,1.5743,2.2651,1.7817,2.0658,1.8895,2.1781,1.7547,2.0377,2.0676,2.0541,2.0518,2.4673,1.9678,1.3360,0.5640,2.4617,1.8698]
           N2_amp_mod=[1.4298,6.0143,1.7139,0.9657,1.5907,2.2883,0.4147,1.5241,1.7481,1.4446,3.3983,1.3409,2.2455,1.4851,2.8533,0.5797,1.2677,1.2086,0.5220,1.3354,2.0463,7.6738,0.5109,1.3800]
           P1_amp_mod=[1.1125,0.7606,0.9898,1.0428,1.1360,1.0159,1.1241,1.1336,1.2444,1.0502,0.9461,1.0612,1.1095,1.0411,1.0197,1.0272,0.7452,0.9095,0.9163,0.8181,0.8803,0.7404,0.9422,0.9659]
           Q1_amp_mod=[0.2648,0.4300,0.2762,0.3006,0.2045,0.1701,0.2570,0.2397,0.3682,0.3012,0.3013,0.3198,0.1608,0.3067,0.2486,0.3141,0.3528,0.3319,0.3836,0.3546,0.2631,0.5689,0.3554,0.2596]
           K2_amp_mod=[0.7690,3.0959,1.0769,0.4023,0.8607,1.1749,0.1129,0.8018,1.6381,0.6961,1.9309,0.6185,1.1977,0.7098,1.6231,0.2102,0.5451,0.4984,0.0962,0.5836,1.0170,3.9425,0.0892,1.2774]
          
           M2_pha_mod=[224.48,46.85,50.99,218.98,223.84,333.68,225.47,223.26,232.33,221.26,49.04,219.45,68.36,221.38,46.80,214.39,217.84,218.76,215.43,217.38,227.84,45.35,215.04,239.98]
           S2_pha_mod=[243.07,72.20,73.03,239.51,240.65,252.65,249.40,240.69,240.35,238.99,72.93,237.56,86.98,237.34,71.89,235.04,237.50,239.90,234.21,236.44,249.58,72.16,229.31,246.72]
           K1_pha_mod=[183.43,129.84,157.05,166.88,187.17,198.10,167.05,185.36,275.60,178.35,147.80,173.09,148.78,180.33,151.87,166.86,165.21,166.15,162.00,164.59,182.69,118.14,162.88,280.23]
           O1_pha_mod=[110.82,137.09,114.90,104.65,113.23,227.57,107.06,112.10,258.97,105.93,118.71,103.32,98.81,107.59,115.61,104.11,103.64,103.44,106.22,102.86,107.10,175.89,106.16,257.60]
           N2_pha_mod=[212.02,30.39,34.51,204.11,211.15,230.85,205.95,210.36,231.98,209.45,32.74,207.23,55.16,210.13,29.81,200.02,205.50,206.13,199.62,205.90,213.15,28.92,199.28,240.55]
           P1_pha_mod=[190.01,125.93,148.35,171.02,187.79,201.20,162.22,188.35,282.93,186.03,144.81,178.01,146.10,183.50,145.92,169.10,171.55,172.76,158.89,180.06,190.69,125.34,158.23,281.80]
           Q1_pha_mod=[38.52,181.41,92.14,47.77,46.37,26.63,54.85,36.46,244.69,52.07,130.23,41.02,72.48,47.60,117.41,47.02,38.17,39.84,78.32,31.83,48.82,197.68,77.97,227.69]
           K2_pha_mod=[237.85,69.95,69.52,234.34,237.23,244.29,252.86,236.33,235.85,234.48,68.84,232.80,82.69,231.78,68.08,234.23,233.84,233.14,233.50,232.17,244.94,67.83,237.59,240.99]


         elif flag_15stats == 1:

           M2_amp_mod=[8.6594,16.9520,6.2713,35.7517]
           S2_amp_mod=[3.6346,6.7628,2.1442,13.7243]
           K1_amp_mod=[3.7932,3.3829,3.5554,2.3121]
           O1_amp_mod=[2.2624,2.0658,1.8895,0.5640]
           N2_amp_mod=[1.7139,3.3983,1.3409,7.6738]
           P1_amp_mod=[0.9898,0.9461,1.0612,0.7404]
           Q1_amp_mod=[0.2762,0.3013,0.3198,0.5689]
           K2_amp_mod=[1.0769,1.9309,0.6185,3.9425]
          
           M2_pha_mod=[50.99,49.04,219.45,45.35]
           S2_pha_mod=[73.03,72.93,237.56,72.16]
           K1_pha_mod=[157.05,147.80,173.09,118.14]
           O1_pha_mod=[114.90,118.71,103.32,175.89]
           N2_pha_mod=[34.51,32.74,207.23,28.92]
           P1_pha_mod=[148.35,144.81,178.01,125.34]
           Q1_pha_mod=[92.14,130.23,41.02,197.68]
           K2_pha_mod=[69.52,68.84,232.80,67.83]

##################################
# ============ ISPRA + EMODnet


TOT_stations_lab_Garea=[]
TOT_A_obs_Garea=[]
TOT_P_obs_Garea=[]
TOT_A_mod_Garea=[]
TOT_P_mod_Garea=[]

TOT_stations_lab_Oarea=[]
TOT_A_obs_Oarea=[]
TOT_P_obs_Oarea=[]
TOT_A_mod_Oarea=[]
TOT_P_mod_Oarea=[]

TOT_stations_lab_Aarea=[]
TOT_A_obs_Aarea=[]
TOT_P_obs_Aarea=[]
TOT_A_mod_Aarea=[]
TOT_P_mod_Aarea=[]

TOT_stations_lab_Tarea=[]
TOT_A_obs_Tarea=[]
TOT_P_obs_Tarea=[]
TOT_A_mod_Tarea=[]
TOT_P_mod_Tarea=[]

TOT_stations_lab_TAarea=[]
TOT_A_obs_TAarea=[]
TOT_P_obs_TAarea=[]
TOT_A_mod_TAarea=[]
TOT_P_mod_TAarea=[]

TOT_stations_lab_TGarea=[]
TOT_A_obs_TGarea=[]
TOT_P_obs_TGarea=[]
TOT_A_mod_TGarea=[]
TOT_P_mod_TGarea=[]

TOT_stations_lab_Earea=[]
TOT_A_obs_Earea=[]
TOT_P_obs_Earea=[]
TOT_A_mod_Earea=[]
TOT_P_mod_Earea=[]

TOT_stations_lab_ord=[]
TOT_A_obs_ord=[]
TOT_P_obs_ord=[]
TOT_A_mod_ord=[]
TOT_P_mod_ord=[]

# Initialize the tables for Amp and Pha statistics
    # Table for TEX Amp
if where_box=='Med':
   Amp_file = open(path+"amp_stats.txt","w")
elif where_box=='AtlBox':
   Amp_file = open(path+"amp_stats_AB.txt","w")
Amp_file.write('\\begin{table}')
Amp_file.write('\\footnotesize')
Amp_file.write('\hspace{-2.5cm}')
Amp_file.write('\\begin{tabular}{||c||c|c|c||c|c|c||}')
Amp_file.write('     \hline')
Amp_file.write('     \hline')
Amp_file.write('     & \multicolumn{3}{c||}{Absolute bias $\|Mod-Obs\|$} & \multicolumn{3}{c||}{Relative bias $\|\\frac{Mod-Obs}{Obs}\|\cdot 100$ } \\\\')
Amp_file.write('      & \multicolumn{3}{c||}{}& \multicolumn{3}{c||}{} \\\\')
Amp_file.write('     \hline')
Amp_file.write('     Tidal & Max & 95th percentile& Mean & Max & 95th percentile& Mean \\\\')
Amp_file.write('     components & & & & & & \\\\')
Amp_file.write('     \hline')
Amp_file.write('     \hline')



    # Table for TEX Pha
if where_box=='Med':
   Pha_file = open(path+"pha_stats.txt","w")
elif where_box=='AtlBox':
   Pha_file = open(path+"pha_stats_AB.txt","w")
Pha_file.write('\\begin{table}')
Pha_file.write('\\begin{tabular}{||c||c|c|c||}')
Pha_file.write('     \hline')
Pha_file.write('     \hline')
Pha_file.write('     & \multicolumn{3}{c||}{Absolute bias $\|Mod-Obs\|$} \\\\')
Pha_file.write('     \hline')
Pha_file.write('     Tidal & Max & 95th percentile & Mean \\\\')
Pha_file.write('     components & & & \\\\')
Pha_file.write('     \hline')
Pha_file.write('     \hline')

# Initialize the tables for lin reg statistics
if where_box=='Med':
   LinReg_file = open(path+"linreg_stats.txt","w")
elif where_box=='AtlBox':
   LinReg_file = open(path+"linreg_stats_AB.txt","w")
LinReg_file.write('\\begin{table}')
LinReg_file.write('\\begin{tabular}{||c||c|c||c|c||}')
LinReg_file.write('     \hline')
LinReg_file.write('     \hline')
LinReg_file.write('     & \multicolumn{2}{c||}{Amplitudes} & \multicolumn{2}{c||}{Phases} \\\\')
LinReg_file.write('     Tidal components & Slope & $R^{2}$ & Slope & $R^{2}$ \\\\')
LinReg_file.write('     \hline')
LinReg_file.write('     \hline')


# Initialize global lists of arrays (with all the constituents) 
N_comp=8
N_stz=len(I_stations_lab+stations_lab)

GLOB_A_mod=[[ 0 for i in range(N_comp+1) ] for j in range(N_stz+1) ]
GLOB_P_mod=[[ 0 for i in range(N_comp+1) ] for j in range(N_stz+1) ] 
GLOB_A_obs=[[ 0 for i in range(N_comp+1) ] for j in range(N_stz+1) ] 
GLOB_P_obs=[[ 0 for i in range(N_comp+1) ] for j in range(N_stz+1) ] 

# Initialize global matrix for Foreman distances and Root Mean Square distance
d_foreman=[[ 0 for i in range(N_comp+1) ] for j in range(N_stz+1) ]
RMSd=[0 for i in range(N_comp)]

##############################################################
# Loop on the 8 components 
comp_idx=0
for comp in ('M2','S2','K1','O1','N2','P1','Q1','K2'):
  if 1 == 1:

    TOT_stations_lab_Garea=[]
    TOT_A_obs_Garea=[]
    TOT_P_obs_Garea=[]
    TOT_A_mod_Garea=[]
    TOT_P_mod_Garea=[]
    
    TOT_stations_lab_Oarea=[]
    TOT_A_obs_Oarea=[]
    TOT_P_obs_Oarea=[]
    TOT_A_mod_Oarea=[]
    TOT_P_mod_Oarea=[]
    
    TOT_stations_lab_Aarea=[]
    TOT_A_obs_Aarea=[]
    TOT_P_obs_Aarea=[]
    TOT_A_mod_Aarea=[]
    TOT_P_mod_Aarea=[]
    
    TOT_stations_lab_Marea=[]
    TOT_A_obs_Marea=[]
    TOT_P_obs_Marea=[]
    TOT_A_mod_Marea=[]
    TOT_P_mod_Marea=[]

    TOT_stations_lab_Tarea=[] # Atlantic area
    TOT_A_obs_Tarea=[]
    TOT_P_obs_Tarea=[]
    TOT_A_mod_Tarea=[]
    TOT_P_mod_Tarea=[]

    TOT_stations_lab_TAarea=[] # Atlantic area
    TOT_A_obs_TAarea=[]
    TOT_P_obs_TAarea=[]
    TOT_A_mod_TAarea=[]
    TOT_P_mod_TAarea=[]

    TOT_stations_lab_TGarea=[] # Atlantic area
    TOT_A_obs_TGarea=[]
    TOT_P_obs_TGarea=[]
    TOT_A_mod_TGarea=[]
    TOT_P_mod_TGarea=[]

    TOT_stations_lab_Earea=[] # East Med Area
    TOT_A_obs_Earea=[]
    TOT_P_obs_Earea=[]
    TOT_A_mod_Earea=[]
    TOT_P_mod_Earea=[]


    TOT_stations_lab_ord=[]
    TOT_A_obs_ord=[]
    TOT_P_obs_ord=[]
    TOT_A_mod_ord=[]
    TOT_P_mod_ord=[]


    nameA_obs=comp+'_amp_obs'
    nameA_mod=comp+'_amp_mod'

    nameP_obs=comp+'_pha_obs'
    nameP_mod=comp+'_pha_mod'
 
    I_nameA_obs='I_'+comp+'_amp_obs'
    I_nameA_mod='I_'+comp+'_amp_mod'

    I_nameP_obs='I_'+comp+'_pha_obs'
    I_nameP_mod='I_'+comp+'_pha_mod'

    TOT_A_obs=np.append(globals()[nameA_obs],globals()[I_nameA_obs])
    TOT_A_obs=np.append(np.array(globals()[nameA_obs]),np.array(globals()[I_nameA_obs]))    
    TOT_A_mod=np.append(np.array(globals()[nameA_mod]),np.array(globals()[I_nameA_mod]))    

    TOT_P_obs=np.append(np.array(globals()[nameP_obs]),np.array(globals()[I_nameP_obs]))    
    TOT_P_mod=np.append(np.array(globals()[nameP_mod]),np.array(globals()[I_nameP_mod]))

    TOT_P_obs=TOT_P_obs[:]
    TOT_P_mod=TOT_P_mod[:]
    #
    for array_idx in range (0,len(TOT_P_obs)):
       while TOT_P_obs[array_idx] > 360:
          TOT_P_obs[array_idx]=TOT_P_obs[array_idx]-360
       while TOT_P_obs[array_idx] < 0:
          TOT_P_obs[array_idx]=TOT_P_obs[array_idx]+360
    for array_idx in range (0,len(TOT_P_mod)):
       while TOT_P_mod[array_idx] > 360:
          TOT_P_mod[array_idx]=TOT_P_mod[array_idx]-360
       while TOT_P_mod[array_idx] < 0:
          TOT_P_mod[array_idx]=TOT_P_mod[array_idx]+360

    # Adjustment for linear regressions
    for idx_diff_mo in range(0,len(TOT_P_mod)):
        #
        if TOT_P_mod[idx_diff_mo]-TOT_P_obs[idx_diff_mo] > 340:
           TOT_P_obs[idx_diff_mo]=TOT_P_obs[idx_diff_mo]+360
        elif TOT_P_mod[idx_diff_mo]-TOT_P_obs[idx_diff_mo] < -340 :
           TOT_P_mod[idx_diff_mo]=TOT_P_mod[idx_diff_mo]+360
           if TOT_P_mod[idx_diff_mo]-TOT_P_obs[idx_diff_mo] < -340 :
              TOT_P_mod[idx_diff_mo]=TOT_P_mod[idx_diff_mo]+360
        #
        #if TOT_P_mod[idx_diff_mo] < 40:
        #   TOT_P_mod[idx_diff_mo]=TOT_P_mod[idx_diff_mo]+360
        #if TOT_P_obs[idx_diff_mo] < 40 :
        #   TOT_P_obs[idx_diff_mo]=TOT_P_obs[idx_diff_mo]+360


    # Adjustment for Amp and Pha plots wrt TPXO9 
    if where_box == 'Med' and len(TOT_P_mod)!= 15: #and tpxo_flag == 1:
     # For simu_tides8_v31 and eas6_v2:
     # tmp messina and lampedusa (and palinuro)
      if comp == 'M2':
         TOT_P_mod[30]=TOT_P_mod[30]+360
      #if comp == 'O1':
      #   TOT_P_mod[30]=TOT_P_mod[30]-360
      if comp == 'N2':
         TOT_P_obs[30]=TOT_P_obs[30]-360
      if comp == 'P1':
         #TOT_P_mod[30]=TOT_P_mod[30]-360
         TOT_P_mod[28]=TOT_P_mod[28]-360
      if comp == 'K1':
         TOT_P_mod[30]=TOT_P_mod[30]+360
      #if comp == 'K2':
      #   TOT_P_mod[30]=TOT_P_mod[30]+360
      #   TOT_P_obs[30]=TOT_P_obs[30]+360
      # For eas6_v2
      if comp == 'Q1':
         TOT_P_mod[32]=TOT_P_mod[32]-360

    # Cut 3 characters for stns names (if flag_tab != 1)
    TOT_stations_lab=np.append(stations_lab,I_stations_lab)
    TOT_stations_lab_long=TOT_stations_lab
    if flag_tab != 1:
       lab=0
       for stz_nm in TOT_stations_lab:
           TOT_stations_lab[lab]=stz_nm[:3]
           lab=lab+1


    TOT_stations_col=np.append(stations_col,I_stations_col)
    # Split in different datasets: Gibraltar area, Adriatic area and others


    # Count stz per area (color based)
    howmany_Garea=0
    howmany_Aarea=0
    howmany_Oarea=0
    howmany_Marea=0
    howmany_Tarea=0
    howmany_TAarea=0
    howmany_TGarea=0
    howmany_Earea=0


    for idx_area in range(0,len(TOT_stations_col)):
        if TOT_stations_col[idx_area] == 'red': # Gibraltar area
           TOT_stations_lab_Garea.append(TOT_stations_lab[idx_area])
           TOT_A_obs_Garea.append(TOT_A_obs[idx_area])
           TOT_P_obs_Garea.append(TOT_P_obs[idx_area])   
           TOT_A_mod_Garea.append(TOT_A_mod[idx_area])
           TOT_P_mod_Garea.append(TOT_P_mod[idx_area])
           howmany_Garea=howmany_Garea+1
        elif TOT_stations_col[idx_area] == 'blue': # Adriatic area
           TOT_stations_lab_Aarea.append(TOT_stations_lab[idx_area])
           TOT_A_obs_Aarea.append(TOT_A_obs[idx_area])
           TOT_P_obs_Aarea.append(TOT_P_obs[idx_area])
           TOT_A_mod_Aarea.append(TOT_A_mod[idx_area])
           TOT_P_mod_Aarea.append(TOT_P_mod[idx_area])
           howmany_Aarea=howmany_Aarea+1
        elif TOT_stations_col[idx_area] == 'green': # Other areas
           TOT_stations_lab_Oarea.append(TOT_stations_lab[idx_area])
           TOT_A_obs_Oarea.append(TOT_A_obs[idx_area])
           TOT_P_obs_Oarea.append(TOT_P_obs[idx_area])
           TOT_A_mod_Oarea.append(TOT_A_mod[idx_area])
           TOT_P_mod_Oarea.append(TOT_P_mod[idx_area])
           howmany_Oarea=howmany_Oarea+1
        elif TOT_stations_col[idx_area] == 'orange': # Messina area
           TOT_stations_lab_Marea.append(TOT_stations_lab[idx_area])
           TOT_A_obs_Marea.append(TOT_A_obs[idx_area])
           TOT_P_obs_Marea.append(TOT_P_obs[idx_area])
           TOT_A_mod_Marea.append(TOT_A_mod[idx_area])
           TOT_P_mod_Marea.append(TOT_P_mod[idx_area])
           howmany_Marea=howmany_Marea+1
        elif TOT_stations_col[idx_area] == 'cyan': # Atlantic Box Biscay Bay area
           TOT_stations_lab_Tarea.append(TOT_stations_lab[idx_area])
           TOT_A_obs_Tarea.append(TOT_A_obs[idx_area])
           TOT_P_obs_Tarea.append(TOT_P_obs[idx_area])
           TOT_A_mod_Tarea.append(TOT_A_mod[idx_area])
           TOT_P_mod_Tarea.append(TOT_P_mod[idx_area])
           howmany_Tarea=howmany_Tarea+1
        elif TOT_stations_col[idx_area] == 'tab:olive': # Atlantic Box Atlantic coast area
           TOT_stations_lab_TAarea.append(TOT_stations_lab[idx_area])
           TOT_A_obs_TAarea.append(TOT_A_obs[idx_area])
           TOT_P_obs_TAarea.append(TOT_P_obs[idx_area])
           TOT_A_mod_TAarea.append(TOT_A_mod[idx_area])
           TOT_P_mod_TAarea.append(TOT_P_mod[idx_area])
           howmany_TAarea=howmany_TAarea+1
        elif TOT_stations_col[idx_area] == 'deeppink': # Atlantic Box Gibrltar area
           TOT_stations_lab_TGarea.append(TOT_stations_lab[idx_area])
           TOT_A_obs_TGarea.append(TOT_A_obs[idx_area])
           TOT_P_obs_TGarea.append(TOT_P_obs[idx_area])
           TOT_A_mod_TGarea.append(TOT_A_mod[idx_area])
           TOT_P_mod_TGarea.append(TOT_P_mod[idx_area])
           howmany_TGarea=howmany_TGarea+1
        elif TOT_stations_col[idx_area] == 'magenta': # East Med area
           TOT_stations_lab_Earea.append(TOT_stations_lab[idx_area])
           TOT_A_obs_Earea.append(TOT_A_obs[idx_area])
           TOT_P_obs_Earea.append(TOT_P_obs[idx_area])
           TOT_A_mod_Earea.append(TOT_A_mod[idx_area])
           TOT_P_mod_Earea.append(TOT_P_mod[idx_area])
           howmany_Earea=howmany_Earea+1
    #print ('PROVA Bay Tarea:',TOT_stations_lab_Tarea,howmany_Tarea,TOT_A_mod_Tarea,TOT_A_obs_Tarea)
    #print ('PROVA girne Earea:',TOT_stations_lab_Earea,howmany_Earea,TOT_A_mod_Earea,TOT_A_obs_Earea)

    TOT_A_obs_ord=TOT_A_obs_Oarea
    TOT_A_obs_ord.extend(TOT_A_obs_Marea)
    TOT_A_obs_ord.extend(TOT_A_obs_Garea)
    TOT_A_obs_ord.extend(TOT_A_obs_Aarea)
    TOT_A_obs_ord.extend(TOT_A_obs_Tarea)
    TOT_A_obs_ord.extend(TOT_A_obs_TAarea)
    TOT_A_obs_ord.extend(TOT_A_obs_TGarea)
    TOT_A_obs_ord.extend(TOT_A_obs_Earea)

    TOT_P_obs_ord=TOT_P_obs_Oarea
    TOT_P_obs_ord.extend(TOT_P_obs_Marea)
    TOT_P_obs_ord.extend(TOT_P_obs_Garea)
    TOT_P_obs_ord.extend(TOT_P_obs_Aarea)
    TOT_P_obs_ord.extend(TOT_P_obs_Tarea)
    TOT_P_obs_ord.extend(TOT_P_obs_TAarea)
    TOT_P_obs_ord.extend(TOT_P_obs_TGarea)
    TOT_P_obs_ord.extend(TOT_P_obs_Earea)

    TOT_A_mod_ord=TOT_A_mod_Oarea
    TOT_A_mod_ord.extend(TOT_A_mod_Marea)
    TOT_A_mod_ord.extend(TOT_A_mod_Garea) 
    TOT_A_mod_ord.extend(TOT_A_mod_Aarea)
    TOT_A_mod_ord.extend(TOT_A_mod_Tarea)
    TOT_A_mod_ord.extend(TOT_A_mod_TAarea)
    TOT_A_mod_ord.extend(TOT_A_mod_TGarea)
    TOT_A_mod_ord.extend(TOT_A_mod_Earea)

    TOT_P_mod_ord=TOT_P_mod_Oarea
    TOT_P_mod_ord.extend(TOT_P_mod_Marea)
    TOT_P_mod_ord.extend(TOT_P_mod_Garea)
    TOT_P_mod_ord.extend(TOT_P_mod_Aarea)
    TOT_P_mod_ord.extend(TOT_P_mod_Tarea)
    TOT_P_mod_ord.extend(TOT_P_mod_TAarea)
    TOT_P_mod_ord.extend(TOT_P_mod_TGarea)
    TOT_P_mod_ord.extend(TOT_P_mod_Earea)

    TOT_stations_lab_ord=TOT_stations_lab_Oarea
    TOT_stations_lab_ord.extend(TOT_stations_lab_Marea)
    TOT_stations_lab_ord.extend(TOT_stations_lab_Garea)
    TOT_stations_lab_ord.extend(TOT_stations_lab_Aarea) 
    TOT_stations_lab_ord.extend(TOT_stations_lab_Tarea)
    TOT_stations_lab_ord.extend(TOT_stations_lab_TAarea)
    TOT_stations_lab_ord.extend(TOT_stations_lab_TGarea)
    TOT_stations_lab_ord.extend(TOT_stations_lab_Earea)

    #print ('PROVA TOT stz lab, TOT A mod, TOT A obs:',TOT_stations_lab_ord,TOT_A_mod_ord,TOT_A_obs_ord)

    TOT_A_obs=[]
    TOT_P_obs=[]
    TOT_A_mod=[]
    TOT_P_mod=[]
    TOT_stations_lab=[]

    TOT_A_obs=np.array(TOT_A_obs_ord)
    TOT_P_obs=np.array(TOT_P_obs_ord)
    #for i in range (0,len(TOT_P_obs_ord)):
    # if TOT_P_obs[i]<40:
    #    TOT_P_obs[i]=TOT_P_obs[i]+360
    # if TOT_P_obs[i]>400:
    #    TOT_P_obs[i]=TOT_P_obs[i]-360
    TOT_A_mod=np.array(TOT_A_mod_ord)
    TOT_P_mod=np.array(TOT_P_mod_ord)
    #for i in range (0,len(TOT_P_mod_ord)):
    # if TOT_P_mod[i]<40:
    #    TOT_P_mod[i]=TOT_P_mod[i]+360
    # if TOT_P_mod[i]>400:
    #    TOT_P_mod[i]=TOT_P_mod[i]-360
    TOT_stations_lab=TOT_stations_lab_ord

    #print ('PROVA TOT stz lab, TOT A mod, TOT A obs:',TOT_stations_lab,TOT_A_mod,TOT_A_obs)
    
    if ttide_flag != 1 and flag_15stats == 0:
       # Messina corr
       if comp == 'M2':
          TOT_P_obs[25]=TOT_P_obs[25]+360
       if comp == 'K1':
          TOT_P_obs=np.where(TOT_P_obs<200,TOT_P_obs+180,TOT_P_obs-180)
       elif comp == 'Q1' or comp == 'P1' or comp == 'O1':
          TOT_P_obs=np.where(TOT_P_obs<150,TOT_P_obs+180,TOT_P_obs-180)
       elif comp == 'N2':
          TOT_P_obs=np.where(TOT_P_obs<350,TOT_P_obs,TOT_P_obs-360)

    # Diff and diff_mean computation
    # diff Amp
    x_textA=[]
    y_textA=[]
    x_textA=TOT_A_obs
    y_textA=TOT_A_mod
    diffA_mo=y_textA-x_textA


    diffA_mo_Oarea=diffA_mo[0:howmany_Oarea]
    diffA_mo_Marea=diffA_mo[howmany_Oarea:howmany_Oarea+howmany_Marea]
    diffA_mo_Garea=diffA_mo[howmany_Oarea+howmany_Marea:howmany_Oarea+howmany_Marea+howmany_Garea]
    diffA_mo_Aarea=diffA_mo[howmany_Oarea+howmany_Marea+howmany_Garea:howmany_Oarea+howmany_Marea+howmany_Garea+howmany_Aarea]
    diffA_mo_Tarea=diffA_mo[howmany_Oarea+howmany_Marea+howmany_Garea+howmany_Aarea:howmany_Oarea+howmany_Marea+howmany_Garea+howmany_Aarea+howmany_Tarea]
    diffA_mo_TAarea=diffA_mo[howmany_Oarea+howmany_Marea+howmany_Garea+howmany_Aarea+howmany_Tarea:howmany_Oarea+howmany_Marea+howmany_Garea+howmany_Aarea+howmany_Tarea+howmany_TAarea]
    diffA_mo_TGarea=diffA_mo[howmany_Oarea+howmany_Marea+howmany_Garea+howmany_Aarea+howmany_Tarea+howmany_TAarea:howmany_Oarea+howmany_Marea+howmany_Garea+howmany_Aarea+howmany_Tarea+howmany_TAarea+howmany_TGarea]
    diffA_mo_Earea=diffA_mo[howmany_Oarea+howmany_Marea+howmany_Garea+howmany_Aarea+howmany_Tarea+howmany_TAarea+howmany_TGarea:howmany_Oarea+howmany_Marea+howmany_Garea+howmany_Aarea+howmany_Tarea+howmany_TAarea+howmany_TGarea+howmany_Earea]


    pdiffA_mo=(y_textA-x_textA)*100.0/x_textA
    mean_diffA=np.mean(diffA_mo)

    # STATISTICS x TAB AMP
    meanabs_diffA=np.mean(abs(diffA_mo))
    max_diffA=np.max(abs(diffA_mo))
    wheremax_AmpAbs=TOT_stations_lab[np.argmax(abs(diffA_mo))]
    perc95_diffA=np.percentile(abs(diffA_mo),95)
    mean_pdiffA=np.mean(pdiffA_mo)
    meanabs_pdiffA=np.mean(abs(pdiffA_mo))
    max_pdiffA=np.max(abs(pdiffA_mo))
    wheremax_AmpPerc=TOT_stations_lab[np.argmax(abs(pdiffA_mo))]
    perc95_pdiffA=np.percentile(abs(pdiffA_mo),95)

    # Table for TEX Amp
    #print(comp,' &',str(round(max_diffA,2)),' cm ({\color{red}{',wheremax_AmpAbs,'}})&',str(round(perc95_diffA,2)),' cm &',str(round(meanabs_diffA,1)),' cm &',str(round(max_pdiffA,2)),' \% ({\color{orange}{',wheremax_AmpPerc,'}})&',str(round(perc95_pdiffA,2)),' \% &',str(round(meanabs_pdiffA,1)),'\% \\\\ ')
    print(comp,' &',str(round(max_diffA,2)),' cm ({\color{red}{',wheremax_AmpAbs,'}})&',str(round(perc95_diffA,2)),' cm &',str(round(meanabs_diffA,1)),' cm &',str(round(max_pdiffA,2)),' \% ({\color{orange}{',wheremax_AmpPerc,'}})&',str(round(perc95_pdiffA,2)),' \% &',str(round(meanabs_pdiffA,1)),'\% \\\\ ', file=Amp_file)
    #print(comp,' &',str(round(max_diffA,2)),' cm ({\color{red}{Algeciras}})&',str(round(perc95_diffA,2)),' cm &',str(round(meanabs_diffA,1)),' cm &',str(round(max_pdiffA,2)),' \% ({\color{orange}{Messina}})&',str(round(perc95_pdiffA,2)),' \% &',str(round(meanabs_pdiffA,1)),'\% \\\\ ', file=Amp_file)
    print('\hline', file=Amp_file)


    # diff Pha
    x_textP=[]
    y_textP=[]
    x_textP=TOT_P_obs
    y_textP=TOT_P_mod
    if cos_pha == 0:
       diffP_mo=y_textP-x_textP
       pdiffP_mo=(y_textP-x_textP)*100.0/x_textP
    elif cos_pha == 1:
       diffP_mo=np.cos(np.array(y_textP)*np.pi/180)-np.cos(np.array(x_textP)*np.pi/180)
       pdiffP_mo=(np.cos(np.array(y_textP)*np.pi/180)-np.cos(np.array(x_textP))*np.pi/180)*100.0/np.cos(np.array(x_textP)*np.pi/180)

    diffP_mo_Oarea=diffP_mo[0:howmany_Oarea]
    diffP_mo_Marea=diffP_mo[howmany_Oarea:howmany_Oarea+howmany_Marea]
    diffP_mo_Garea=diffP_mo[howmany_Oarea+howmany_Marea:howmany_Oarea+howmany_Marea+howmany_Garea]
    diffP_mo_Aarea=diffP_mo[howmany_Oarea+howmany_Marea+howmany_Garea:howmany_Oarea+howmany_Marea+howmany_Garea+howmany_Aarea]
    diffP_mo_Tarea=diffP_mo[howmany_Oarea+howmany_Marea+howmany_Garea+howmany_Aarea:howmany_Oarea+howmany_Marea+howmany_Garea+howmany_Aarea+howmany_Tarea]
    diffP_mo_TAarea=diffP_mo[howmany_Oarea+howmany_Marea+howmany_Garea+howmany_Aarea+howmany_Tarea:howmany_Oarea+howmany_Marea+howmany_Garea+howmany_Aarea+howmany_Tarea+howmany_TAarea]
    diffP_mo_TGarea=diffP_mo[howmany_Oarea+howmany_Marea+howmany_Garea+howmany_Aarea+howmany_Tarea+howmany_TAarea:howmany_Oarea+howmany_Marea+howmany_Garea+howmany_Aarea+howmany_Tarea+howmany_TAarea+howmany_TGarea]
    diffP_mo_Earea=diffP_mo[howmany_Oarea+howmany_Marea+howmany_Garea+howmany_Aarea+howmany_Tarea+howmany_TAarea+howmany_TGarea:howmany_Oarea+howmany_Marea+howmany_Garea+howmany_Aarea+howmany_Tarea+howmany_TAarea+howmany_TGarea+howmany_Earea]

       
    mean_diffP=np.mean(diffP_mo)

    # STATISTICS x TAB PHA
    max_diffP=np.max(abs(diffP_mo))
    wheremax_PhaAbs=TOT_stations_lab[np.argmax(abs(diffP_mo))]
    perc95_diffP=np.percentile(abs(diffP_mo),95)
    mean_pdiffP=np.mean(pdiffP_mo)
    meanabs_pdiffP=np.mean(abs(pdiffP_mo))
    max_pdiffP=np.max(abs(pdiffP_mo))
    perc95_pdiffP=np.percentile(abs(pdiffP_mo),95)

    # Table for TEX Pha
    print(comp,' & ',round(max_diffP,2),'$^{\circ}$ ({\color{orange}{',wheremax_PhaAbs,'}})& ',round(perc95_diffP,2),'$^{\circ}$ &',round(meanabs_pdiffP,1),'$^{\circ}$ \\\\ ', file=Pha_file)
    #print(comp,' & ',round(max_diffP,2),'$^{\circ}$ ({\color{orange}{Messina}})& ',round(perc95_diffP,2),'$^{\circ}$ &',round(meanabs_pdiffP,1),'$^{\circ}$ \\\\ ', file=Pha_file)
    print('\hline', file=Pha_file)

    # Plot 1 ( A and P x stz )
    plt.figure(figsize=(24,10))
    plt.rc('font', size=12)
    # Amp 
    plt.subplot(2,1,1)
    plt.xticks(fontsize=12)
    #plt.plot(TOT_stations_lab, np.array(TOT_A_mod), 'r-o', label = 'Tides 8')
    plt.plot(TOT_A_obs, '-s', color = 'black' ,label = 'Obs')
    
    if tpxo_flag == 1:
       TPXO_AMP=globals()['TPXO_'+comp]
       TPXO_PHA=globals()['TPXO_P_'+comp]
       TPXO_AMP=np.multiply(TPXO_AMP,100) # Want cm not m!
       TPXO_PHA=TPXO_PHA #+P_shift[comp_idx]#+P_offset[comp_idx] # rm the tpxo zero offset 
       for i in range (0,len(TPXO_PHA)):
           while TPXO_PHA[i]<0:
              TPXO_PHA[i]=TPXO_PHA[i]+360
              #print ('TPXO_PHA < 0 ', comp)
           while TPXO_PHA[i]>360:
              TPXO_PHA[i]=TPXO_PHA[i]-360
           #print ('TPXO_PHA: ',TPXO_PHA[i])
       plt.plot(TPXO_AMP, '--v', color = 'black' ,label = 'TPXO')
       #print ('Ven',comp,TPXO_AMP[35])

       if where_box == 'Med' and ttide_flag != 1:
       # tmp messina and lampedusa
          if comp == 'M2' or comp == 'K1':
          #if comp == 'M2' or comp == 'K2' or comp == 'K1':
             TPXO_PHA[25]=TPXO_PHA[25]+360
    #
    if howmany_Oarea != 0:
       plt.plot(TOT_stations_lab_Oarea[0:howmany_Oarea], np.array(TOT_A_mod_Oarea[0:howmany_Oarea]), 'g-o', label = 'Mod (Other areas)')
    if howmany_Marea != 0:
       plt.plot(TOT_stations_lab_Marea, np.array(TOT_A_mod_Marea), '-o', color='orange', label = 'Mod (Messina area)')
    if howmany_Garea != 0:
       plt.plot(TOT_stations_lab_Garea, np.array(TOT_A_mod_Garea), 'r-o', label = 'Mod (Gibraltar area)')
    if howmany_Aarea != 0:
       plt.plot(TOT_stations_lab_Aarea, np.array(TOT_A_mod_Aarea), 'b-o', label = 'Mod (Adriatic area)')
    if howmany_Tarea != 0:
       plt.plot(TOT_stations_lab_Tarea, np.array(TOT_A_mod_Tarea), 'c-o', label = 'Mod (Atlantic Box - Biscay Bay)')
    if howmany_TAarea != 0:
       plt.plot(TOT_stations_lab_TAarea, np.array(TOT_A_mod_TAarea), '-o', color='tab:olive', label = 'Mod (Atlantic Box - Atlantic coast)')
    if howmany_TGarea != 0:
       plt.plot(TOT_stations_lab_TGarea, np.array(TOT_A_mod_TGarea), '-o', color='deeppink', label = 'Mod (Atlantic Box - Gibraltar strait)')
    if howmany_Earea != 0:
       plt.plot(TOT_stations_lab_Earea, np.array(TOT_A_mod_Earea), 'm-o', label = 'Mod (East Med area)')
    plt.title(comp+' Amplitude [cm] '+iniend_dates)
    plt.legend( loc='upper left',fontsize = 'large' )
    plt.grid ()
    plt.ylabel ('Amplitude [cm]')
    plt.ylim(bottom=0.0)
    # Amp diff
    plt.subplot(2,1,2)
    plt.xticks(fontsize=12)
    #plt.plot(TOT_stations_lab, diffA_mo, 'r-o', label = 'Tides_8 - Obs')    
    #plt.plot(TOT_stations_lab, diffA_mo, '-o', color = 'black' , label = 'Tides_8 - Obs')
    if howmany_Oarea != 0:
       plt.plot(TOT_stations_lab_Oarea[0:howmany_Oarea], diffA_mo_Oarea, 'g-o', label = 'Mod - Obs (Other areas)')
    if howmany_Marea != 0:
       plt.plot(TOT_stations_lab_Marea, diffA_mo_Marea, '-o', color='orange', label = 'Mod - Obs (Messina area)')
    if howmany_Garea != 0:
       plt.plot(TOT_stations_lab_Garea, diffA_mo_Garea, 'r-o', label = 'Mod - Obs (Gibraltar area)')
    if howmany_Aarea != 0:
       plt.plot(TOT_stations_lab_Aarea, diffA_mo_Aarea, 'b-o', label = 'Mod - Obs (Adriatic area)')
    if howmany_Tarea != 0:
       plt.plot(TOT_stations_lab_Tarea, diffA_mo_Tarea, '-o', color='tab:cyan', label = 'Mod - Obs (Atlantic Box - Biscay Bay)')
    if howmany_TAarea != 0:
       plt.plot(TOT_stations_lab_TAarea, diffA_mo_TAarea, '-o', color='tab:olive', label = 'Mod - Obs (Atlantic Box - Atlantic coast)')
    if howmany_TGarea != 0:
       plt.plot(TOT_stations_lab_TGarea, diffA_mo_TGarea, '-o', color='deeppink', label = 'Mod - Obs (Atlantic Box  - Gibraltar strait)')
    if howmany_Earea != 0:
       plt.plot(TOT_stations_lab_Earea, diffA_mo_Earea, 'm-o', label = 'Mod - Obs (East Med area)')
    if tpxo_flag == 1 and flag_15stats==0 :
       plt.plot(TPXO_AMP-TOT_A_obs, '--v', color = 'black' ,label = 'TPXO-OBS')
    plt.title(comp+' Amplitude diff (Mod - Obs) [cm] ')
    plt.grid ()
    plt.axhline(y=0, color='black')
    App_mean_diffA=round(mean_diffA,2)
    plt.axhline(y=mean_diffA, color='grey', linestyle='dashed' ,label = '<Amp Diff>='+str(App_mean_diffA)+'[cm]')
    plt.ylabel ('Mod - Obs [cm]')
    plt.legend( loc='upper left',fontsize = 'large' )
    # Add %
    i=0
    for word in pdiffA_mo:
         word_approx=round(word,1)
         string_to=str(word_approx)+'%'
         plt.text(TOT_stations_lab[i],diffA_mo[i]+.03,string_to,fontsize=12,color = 'black')
         i=i+1
    
    #plt.savefig(path+comp+'_'+dates_lab+'_Aallmm.jpg')
    if where_box=='Med' and tpxo_flag == 1 :
       plt.savefig(path+comp+'_'+dates_lab+'_Atpxo.jpg')
    elif where_box=='Med' and tpxo_flag == 0 :
       plt.savefig(path+comp+'_'+dates_lab+'_A.jpg')
    elif where_box=='AtlBox':
       plt.savefig(path+comp+'_'+dates_lab+'_A_AB.jpg')
    plt.clf()
       ###
    # Pha
    if cos_pha == 0:

     plt.figure(figsize=(24,10))
     plt.rc('font', size=12)
     plt.subplot(2,1,1)
     plt.xticks(fontsize=12)
     #plt.plot(TOT_stations_lab, TOT_P_mod, 'r-o', label = 'Tides 8')
     plt.plot(TOT_P_obs, '-s', color = 'black' , label = 'Obs')
     if tpxo_flag == 1:
        plt.plot(TPXO_PHA, '--v', color = 'black' ,label = 'TPXO')
 
     if howmany_Oarea != 0:
        plt.plot(TOT_stations_lab_Oarea[0:howmany_Oarea], np.array(TOT_P_mod_Oarea[0:howmany_Oarea]), 'g-o', label = 'Mod (Other areas)')
     if howmany_Marea != 0:
        plt.plot(TOT_stations_lab_Marea, np.array(TOT_P_mod_Marea), '-o',color='orange', label = 'Mod (Messina area)')
     if howmany_Garea != 0:
        plt.plot(TOT_stations_lab_Garea, np.array(TOT_P_mod_Garea), 'r-o', label = 'Mod (Gibraltar area)')
     if howmany_Aarea != 0:
        plt.plot(TOT_stations_lab_Aarea, np.array(TOT_P_mod_Aarea), 'b-o', label = 'Mod (Adriatic area)')
     if howmany_Tarea != 0:
        plt.plot(TOT_stations_lab_Tarea, np.array(TOT_P_mod_Tarea), 'c-o', label = 'Mod (Atlantic Box - Biscay Bay)')
     if howmany_TAarea != 0:
        plt.plot(TOT_stations_lab_TAarea, np.array(TOT_P_mod_TAarea), '-o', color='tab:olive', label = 'Mod (Atlantic Box - Atlantic coast)')
     if howmany_TGarea != 0:
        plt.plot(TOT_stations_lab_TGarea, np.array(TOT_P_mod_TGarea), '-o', color='deeppink', label = 'Mod (Atlantic Box - Gibraltar strait)')
     if howmany_Earea != 0:
        plt.plot(TOT_stations_lab_Earea, np.array(TOT_P_mod_Earea), 'm-o', label = 'Mod (East Med area)')
     plt.title(comp+' Phase [deg] '+iniend_dates)
     plt.grid ()
     plt.ylabel ('Phase [deg]')
     plt.ylim(-50.0, 450.0)
     plt.legend( loc='upper left',fontsize = 'large')
     # Pha diff
     plt.subplot(2,1,2)
     plt.xticks(fontsize=12)
     #plt.plot(TOT_stations_lab,diffP_mo, '-o',color = 'black', label = 'Tides_8 - Obs')
     if howmany_Oarea != 0:
        plt.plot(TOT_stations_lab_Oarea[0:howmany_Oarea], diffP_mo_Oarea, 'g-o', label = 'Mod - Obs (Other areas)')
     if howmany_Marea != 0:
        plt.plot(TOT_stations_lab_Marea, diffP_mo_Marea, '-o', color='orange', label = 'Mod - Obs (Messina area)')
     if howmany_Garea != 0:
        plt.plot(TOT_stations_lab_Garea, diffP_mo_Garea, 'r-o', label = 'Mod - Obs (Gibraltar area)')
     if howmany_Aarea != 0:
        plt.plot(TOT_stations_lab_Aarea, diffP_mo_Aarea, 'b-o', label = 'Mod - Obs (Adriatic area)')
     if howmany_Tarea != 0:
        plt.plot(TOT_stations_lab_Tarea, diffP_mo_Tarea, 'c-o', label = 'Mod - Obs (Atlantic Box - Biscay Bay)')
     if howmany_TAarea != 0:
        plt.plot(TOT_stations_lab_TAarea, diffP_mo_TAarea, '-o', color='tab:olive', label = 'Mod - Obs (Atlantic Box - Atlantic coast)')
     if howmany_TGarea != 0:
        plt.plot(TOT_stations_lab_TGarea, diffP_mo_TGarea, '-o', color='deeppink', label = 'Mod - Obs (Atlantic Box - Gibraltar strait)')
     if howmany_Earea != 0:
        plt.plot(TOT_stations_lab_Earea, diffP_mo_Earea, 'm-o', label = 'Mod - Obs (East Med area)')
     if tpxo_flag == 1 and flag_15stats==0 :
        plt.plot(TPXO_PHA-TOT_P_obs, '--v', color = 'black' ,label = 'TPXO-OBS')
     plt.title(comp+' Phase diff (Mod - Obs) [deg] ')
     plt.grid ()
     plt.ylabel ('Mod - Obs [deg]')
     App_mean_diffP=round(mean_diffP,2)
     plt.axhline(y=mean_diffP, color='gray' , linestyle='dashed' ,label = '<Pha Diff>='+str(App_mean_diffP)+'[deg]')
     plt.axhline(y=0, color='black')
     plt.legend( loc='upper left',fontsize = 'large' )
     # Add %
     i=0
     for word in pdiffP_mo:
          word_approx=round(word,1)
          string_to=str(word_approx)+'%'
          plt.text(TOT_stations_lab[i],diffP_mo[i]+.03,string_to,fontsize=12,color = 'black')
          i=i+1
     #
     #plt.savefig(path+comp+'_'+dates_lab+'_Pallmm.jpg')
     if where_box=='Med' and tpxo_flag == 1 :
       plt.savefig(path+comp+'_'+dates_lab+'_Ptpxo.jpg')
     elif where_box=='Med' and tpxo_flag == 0 :
       plt.savefig(path+comp+'_'+dates_lab+'_P.jpg')
     elif where_box=='AtlBox':
       plt.savefig(path+comp+'_'+dates_lab+'_P_AB.jpg')
     plt.clf()

    elif cos_pha == 1:

     plt.figure(figsize=(24,10))
     plt.rc('font', size=12)
     plt.subplot(2,1,1)
     plt.xticks(fontsize=12)
     #plt.plot(TOT_stations_lab, TOT_P_mod, 'r-o', label = 'Tides 8')
     plt.plot(np.cos(TOT_P_obs*np.pi/180), '-s', color = 'black' , label = 'Obs')
     if tpxo_flag == 1:
        plt.plot(np.cos(TPXO_PHA), '--v', color = 'black' ,label = 'TPXO')
 
     if howmany_Oarea != 0:
        plt.plot(TOT_stations_lab_Oarea[0:howmany_Oarea], np.cos(np.array(TOT_P_mod_Oarea[0:howmany_Oarea])*np.pi/180), 'g-o', label = 'Mod (Other areas)')
     if howmany_Marea != 0:
        plt.plot(TOT_stations_lab_Marea, np.cos(np.array(TOT_P_mod_Marea)*np.pi/180), '-o',color='orange', label = 'Mod (Messina area)')
     if howmany_Garea != 0:
        plt.plot(TOT_stations_lab_Garea, np.cos(np.array(TOT_P_mod_Garea)*np.pi/180), 'r-o', label = 'Mod (Gibraltar area)')
     if howmany_Aarea != 0:
        plt.plot(TOT_stations_lab_Aarea, np.cos(np.array(TOT_P_mod_Aarea)*np.pi/180), 'b-o', label = 'Mod (Adriatic area)')
     if howmany_Tarea != 0:
        plt.plot(TOT_stations_lab_Tarea, np.cos(np.array(TOT_P_mod_Tarea)*np.pi/180), 'c-o', label = 'Mod (Atlantic Box - Biscay Bay)')
     if howmany_TAarea != 0:
        plt.plot(TOT_stations_lab_TAarea, np.cos(np.array(TOT_P_mod_TAarea)*np.pi/180), '-o', color='tab:olive', label = 'Mod (Atlantic Box - Atlantic coast)')
     if howmany_TGarea != 0:
        plt.plot(TOT_stations_lab_TGarea, np.cos(np.array(TOT_P_mod_TGarea)*np.pi/180), '-o', color='deeppink', label = 'Mod (Atlantic Box - Gibraltar strait)')
     if howmany_Earea != 0:
        plt.plot(TOT_stations_lab_Earea, np.cos(np.array(TOT_P_mod_Earea)*np.pi/180), 'm-o', label = 'Mod (East Med area)')
     plt.title(comp+' cos(Phase) '+iniend_dates)
     plt.grid ()
     plt.ylabel ('cos(Phase)')
     plt.ylim(-1, 1)
     plt.legend( loc='upper left',fontsize = 'large')
     # Pha diff
     plt.subplot(2,1,2)
     plt.xticks(fontsize=12)
     #plt.plot(TOT_stations_lab,diffP_mo, '-o',color = 'black', label = 'Tides_8 - Obs')
     if howmany_Oarea != 0:
        plt.plot(TOT_stations_lab_Oarea[0:howmany_Oarea], diffP_mo_Oarea, 'g-o', label = 'Mod - Obs (Other areas)')
     if howmany_Marea != 0:
        plt.plot(TOT_stations_lab_Marea, diffP_mo_Marea, '-o', color='orange', label = 'Mod - Obs (Messina area)')
     if howmany_Garea != 0:
        plt.plot(TOT_stations_lab_Garea, diffP_mo_Garea, 'r-o', label = 'Mod - Obs (Gibraltar area)')
     if howmany_Aarea != 0:
        plt.plot(TOT_stations_lab_Aarea, diffP_mo_Aarea, 'b-o', label = 'Mod - Obs (Adriatic area)')
     if howmany_Tarea != 0:
        plt.plot(TOT_stations_lab_Tarea, diffP_mo_Tarea, 'c-o', label = 'Mod - Obs (Atlantic Box - Biscay Bay)')
     if howmany_TAarea != 0:
        plt.plot(TOT_stations_lab_TAarea, diffP_mo_TAarea, '-o', color='tab:olive', label = 'Mod - Obs (Atlantic Box - Atlantic coast)')
     if howmany_TGarea != 0:
        plt.plot(TOT_stations_lab_TGarea, diffP_mo_TGarea, '-o', color='deeppink', label = 'Mod - Obs (Atlantic Box - Gibraltar strait)')
     if howmany_Earea != 0:
        plt.plot(TOT_stations_lab_Earea, diffP_mo_Earea, 'm-o', label = 'Mod - Obs (East Med area)')
     if tpxo_flag == 1:
        plt.plot(np.cos(np.array(TPXO_PHA)*np.pi/180)-np.cos(np.array(TOT_P_obs)*np.pi/180), '--v', color = 'black' ,label = 'TPXO-OBS')
     plt.title(comp+' cos(Phase) diff (Mod - Obs) ')
     plt.grid ()
     plt.ylabel ('Mod - Obs')
     App_mean_diffP=round(mean_diffP,2)
     plt.axhline(y=mean_diffP, color='gray' , linestyle='dashed' ,label = '<Pha Diff>='+str(App_mean_diffP)+'[deg]')
     plt.axhline(y=0, color='black')
     plt.legend( loc='upper left',fontsize = 'large' )
     # Add %
     i=0
     for word in pdiffP_mo:
          word_approx=round(word,1)
          string_to=str(word_approx)+'%'
          #plt.text(TOT_stations_lab[i],diffP_mo[i]+.03,string_to,fontsize=12,color = 'black')
          i=i+1
     #
     #plt.savefig(path+comp+'_'+dates_lab+'_Pallmm.jpg')
     if where_box=='Med' and tpxo_flag == 1 :
       plt.savefig(path+comp+'_'+dates_lab+'_cosPtpxo.jpg')
     elif where_box=='Med' and tpxo_flag == 0 :
       plt.savefig(path+comp+'_'+dates_lab+'_cosP.jpg')
     elif where_box=='AtlBox':
       plt.savefig(path+comp+'_'+dates_lab+'_cosP_AB.jpg')
     plt.clf()



######### LIN REG


   # Plot 2 ( Lin reg mod Vs obs x A and P x stz )
    if cos_pha == 0:
       plt.figure(figsize=(6,12))
    elif cos_pha == 1:
       plt.figure(figsize=(7,12))
    plt.rc('font', size=12)
    #
    plt.subplot(2,1,1)
    plt.title(comp+' Amplitude [cm] '+'('+iniend_dates+')')
    plt.grid ()
    # Arrays defn
    x_text=[]
    y_text=[]
    y_text2=[]
    x_text=np.array(TOT_A_obs)
    y_text=np.array(TOT_A_mod)
    top=np.maximum(x_text,y_text)
    top=max(top[:])
    # Linear regression
    # Mod T8
    slope, intercept, r_value, p_value, std_err = stats.linregress(x_text,y_text)
    m_A=[]
    q_A=[]
    fitted_A=[]
    cov_A=[]
    def line_A(x, m_A, q_A):
        return (m_A*x+q_A)
    fitted_A, cov_A = curve_fit(line_A,x_text[:],y_text[:])
    m_A=fitted_A[0]
    q_A=fitted_A[1]
    rx=np.linspace(0.0,top)
    retta_A=m_A*rx+q_A
    m_A_approx=round(m_A,2)
    print (cov_A)
    perr = np.abs(np.diag(cov_A))
    print(perr)
    m_Ae_approx=round(perr[0],2)
    r_A=round(r_value,2)
    lr_leg_str='( Slope='+str(m_A_approx)+'; R2='+str(r_A)+')'

    # Arrays defn
    x_text=[]
    x_text=np.array(TOT_A_obs)
    # Plot points
    #plt.plot(np.array(TOT_A_obs), np.array(TOT_A_mod),'ro',label = 'Tides 8 '+lr_leg_str)
    if howmany_Oarea != 0:
      plt.plot(np.array(TOT_A_obs_Oarea), np.array(TOT_A_mod_Oarea), 'go', label = 'Other areas')
    if howmany_Marea != 0:
      plt.plot(np.array(TOT_A_obs_Marea), np.array(TOT_A_mod_Marea), 'o',color='orange', label = 'Messina area')
    if howmany_Garea != 0:
      plt.plot(np.array(TOT_A_obs_Garea), np.array(TOT_A_mod_Garea), 'ro', label = 'Gibraltar area')
    if howmany_Aarea != 0:
      plt.plot(np.array(TOT_A_obs_Aarea), np.array(TOT_A_mod_Aarea), 'bo', label = 'Adriatic area')
    if howmany_Tarea != 0:
      plt.plot(np.array(TOT_A_obs_Tarea), np.array(TOT_A_mod_Tarea), 'co', label = 'Atlantic Box - Biscay Bay')
    if howmany_TAarea != 0:
      plt.plot(np.array(TOT_A_obs_TAarea), np.array(TOT_A_mod_TAarea), 'o',color='tab:olive', label = 'Atlantic Box - Atlantic coast')
    if howmany_TGarea != 0:
      plt.plot(np.array(TOT_A_obs_TGarea), np.array(TOT_A_mod_TGarea), 'o', color='deeppink', label = 'Atlantic Box - Gibraltar strait')
    if howmany_Earea != 0:
      plt.plot(np.array(TOT_A_obs_Earea), np.array(TOT_A_mod_Earea), 'mo', label = 'East Med area')
    # Ax settings
    plt.xlabel ('OBS Amplitude [cm]')
    plt.ylabel ('MOD Amplitude [cm]')
    bottom_ax, top_ax = plt.xlim()
    bottom_ay, top_ay = plt.ylim()
    top_a=np.maximum(top_ax,top_ay)
    bottom_a=np.minimum(bottom_ax,bottom_ay)
    plt.ylim(bottom_a, top_a)
    plt.xlim(bottom_a, top_a)
    # Plot Lines
    plt.plot(rx,retta_A,color = 'red',label=lr_leg_str)
    plt.plot([bottom_a,top_a], [bottom_a,top_a], 'k-', color = 'black')
    # Legend
    plt.legend( loc='upper left' )
    # point label (stn names)
    i=0
    for word in TOT_stations_lab:
        plt.text(x_text[i]+.03,y_text[i]+.03,word,fontsize=12,color = 'black')
        i=i+1

    plt.subplot(2,1,2)
    ### Pha linear reg
    if cos_pha == 0:
     plt.title(comp+' Phase [deg] '+'('+iniend_dates+')')
     plt.grid ()
     plt.ylim(0.0, 450.0)
     plt.xlim(0.0, 450.0)
     plt.plot([0.0, 450.0], [0.0, 450.0], 'k-', color = 'black')
     x_text=[]
     y_text=[]
     x_text=np.array(TOT_P_obs)
     y_text=np.array(TOT_P_mod)
     slopeP, interceptP, r_valueP, p_valueP, std_errP = stats.linregress(x_text,y_text)
     rx=np.linspace(0.0,450.0)
     retta_P=slopeP*rx+interceptP
     #plt.plot(rx,retta_P,color = 'red')
     lr_leg_str='( Slope='+str(round(slopeP,2))+'; R2='+str(round(r_valueP,2))+')'
     plt.plot(rx,retta_P,color = 'red',label=lr_leg_str)
     #plt.plot(np.array(TOT_P_obs), np.array(TOT_P_mod),'ro',label = 'Tides 8 '+lr_leg_str)
     if howmany_Oarea != 0:
        plt.plot(np.array(TOT_P_obs_Oarea), np.array(TOT_P_mod_Oarea), 'go', label = 'Other areas')
     if howmany_Marea != 0:
        plt.plot(np.array(TOT_P_obs_Marea), np.array(TOT_P_mod_Marea), 'o',color='orange', label = 'Messina area')
     if howmany_Garea != 0:
        plt.plot(np.array(TOT_P_obs_Garea), np.array(TOT_P_mod_Garea), 'ro', label = 'Gibraltar area')
     if howmany_Aarea != 0:
        plt.plot(np.array(TOT_P_obs_Aarea), np.array(TOT_P_mod_Aarea), 'bo', label = 'Adriatic area')
     if howmany_Tarea != 0:
        plt.plot(np.array(TOT_P_obs_Tarea), np.array(TOT_P_mod_Tarea), 'co', label = 'Atlantic Box - Biscay Bay')
     if howmany_TAarea != 0:
        plt.plot(np.array(TOT_P_obs_TAarea), np.array(TOT_P_mod_TAarea), 'o',color='tab:olive', label = 'Atlantic Box - Atlantic coast')
     if howmany_TGarea != 0:
        plt.plot(np.array(TOT_P_obs_TGarea), np.array(TOT_P_mod_TGarea), 'o', color='deeppink', label = 'Atlantic Box - Gibraltar strait')
     if howmany_Earea != 0:
        plt.plot(np.array(TOT_P_obs_Earea), np.array(TOT_P_mod_Earea), 'mo', label = 'East Med area') 
 
     # Axes
     plt.xlabel ('OBS Phase [deg]')
     plt.ylabel ('MOD Phase [deg]')
     # Legend
     plt.legend( loc='lower right' )
     i=0
     for word in TOT_stations_lab:
         plt.text(x_text[i]+.03,y_text[i]+.03,word,fontsize=12,color = 'black')
         i=i+1

     if where_box=='Med':
       plt.savefig(path+comp+'_'+dates_lab+'_all_lr.jpg')
     elif where_box=='AtlBox':
       plt.savefig(path+comp+'_'+dates_lab+'_all_lr_AB.jpg')
     plt.clf()

    elif cos_pha == 1:
     plt.title(comp+' cos(Phase)'+'('+iniend_dates+')')
     plt.grid ()
     plt.ylim(-1, 1)
     plt.xlim(-1, 1)
     plt.plot([-1, 1], [-1, 1], 'k-', color = 'black')
     x_text=[]
     y_text=[]
     x_text=np.cos(np.array(TOT_P_obs)*np.pi/180)
     y_text=np.cos(np.array(TOT_P_mod)*np.pi/180)
     slopeP, interceptP, r_valueP, p_valueP, std_errP = stats.linregress(x_text,y_text)
     rx=np.linspace(-1,1)
     retta_P=slopeP*rx+interceptP
     #plt.plot(rx,retta_P,color = 'red')
     lr_leg_str='( Slope='+str(round(slopeP,2))+'; R2='+str(round(r_valueP,2))+')'
     plt.plot(rx,retta_P,color = 'red',label=lr_leg_str)
     #plt.plot(np.array(TOT_P_obs), np.array(TOT_P_mod),'ro',label = 'Tides 8 '+lr_leg_str)
     if howmany_Oarea != 0:
        plt.plot(np.cos(np.array(TOT_P_obs_Oarea)*np.pi/180), np.cos(np.array(TOT_P_mod_Oarea)*np.pi/180), 'go', label = 'Other areas')
     if howmany_Marea != 0:
        plt.plot(np.cos(np.array(TOT_P_obs_Marea)*np.pi/180), np.cos(np.array(TOT_P_mod_Marea)*np.pi/180), 'o',color='orange', label = 'Messina area')
     if howmany_Garea != 0:
        plt.plot(np.cos(np.array(TOT_P_obs_Garea)*np.pi/180), np.cos(np.array(TOT_P_mod_Garea)*np.pi/180), 'ro', label = 'Gibraltar area')
     if howmany_Aarea != 0:
        plt.plot(np.cos(np.array(TOT_P_obs_Aarea)*np.pi/180), np.cos(np.array(TOT_P_mod_Aarea)*np.pi/180), 'bo', label = 'Adriatic area')
     if howmany_Tarea != 0:
        plt.plot(np.cos(np.array(TOT_P_obs_Tarea)*np.pi/180), np.cos(np.array(TOT_P_mod_Tarea)*np.pi/180), 'co', label = 'Atlantic Box - Biscay Bay')
     if howmany_TAarea != 0:
        plt.plot(np.cos(np.array(TOT_P_obs_TAarea)*np.pi/180), np.cos(np.array(TOT_P_mod_TAarea)*np.pi/180), 'o',color='tab:olive', label = 'Atlantic Box - Atlantic coast')
     if howmany_TGarea != 0:
        plt.plot(np.cos(np.array(TOT_P_obs_TGarea)*np.pi/180), np.cos(np.array(TOT_P_mod_TGarea)*np.pi/180), 'o', color='deeppink', label = 'Atlantic Box - Gibraltar strait')
     if howmany_Earea != 0:
        plt.plot(np.cos(np.array(TOT_P_obs_Earea)*np.pi/180), np.cos(np.array(TOT_P_mod_Earea)*np.pi/180), 'mo', label = 'East Med area')
 
     # Axes
     plt.xlabel ('OBS cos(Phase)')
     plt.ylabel ('MOD cos(Phase)')
     # Legend
     plt.legend( loc='lower right' )
     i=0
     for word in TOT_stations_lab:
         plt.text(x_text[i]+.03,y_text[i]+.03,word,fontsize=12,color = 'black')
         i=i+1

     if where_box=='Med':
       plt.savefig(path+comp+'_'+dates_lab+'_all_lr_cos.jpg')
     elif where_box=='AtlBox':
       plt.savefig(path+comp+'_'+dates_lab+'_all_lr_cos_AB.jpg')
     plt.clf()


    print(comp,' & ',str(m_A_approx),' & ',str(r_A),' & ',str(round(slopeP,2)),' & ',str(round(r_valueP,2)),' \\\\ ', file=LinReg_file)



  # Save val in GLOBAL arrays

  # Components names
  GLOB_A_mod[0][comp_idx+1]=comp
  GLOB_P_mod[0][comp_idx+1]=comp
  GLOB_A_obs[0][comp_idx+1]=comp
  GLOB_P_obs[0][comp_idx+1]=comp

  d_foreman[0][comp_idx+1]=comp 

  #if comp_idx == 0:
  # Stations names
  for nnn_stz in range (0,len(TOT_stations_lab_ord)):
      GLOB_A_mod[nnn_stz+1][0]=TOT_stations_lab_ord[nnn_stz]
      GLOB_P_mod[nnn_stz+1][0]=TOT_stations_lab_ord[nnn_stz]
      GLOB_A_obs[nnn_stz+1][0]=TOT_stations_lab_ord[nnn_stz]
      GLOB_P_obs[nnn_stz+1][0]=TOT_stations_lab_ord[nnn_stz]

      d_foreman[nnn_stz+1][0]=TOT_stations_lab_ord[nnn_stz]

  # Put right numbers in the matrix!
  for nnn_AP in range (0,len(TOT_stations_lab_ord)):
      GLOB_A_mod[nnn_AP+1][comp_idx+1]=TOT_A_mod_ord[nnn_AP]
      GLOB_P_mod[nnn_AP+1][comp_idx+1]=TOT_P_mod_ord[nnn_AP]
      GLOB_A_obs[nnn_AP+1][comp_idx+1]=TOT_A_obs_ord[nnn_AP]
      GLOB_P_obs[nnn_AP+1][comp_idx+1]=TOT_P_obs_ord[nnn_AP]


      # Compute Distances in the complex plane [Foreman et al. 93]
      d_foreman[nnn_AP+1][comp_idx+1]=np.sqrt((TOT_A_obs_ord[nnn_AP]*np.cos((np.pi/180.0)*TOT_P_obs_ord[nnn_AP])-(TOT_A_mod_ord[nnn_AP]*np.cos((np.pi/180.0)*TOT_P_mod_ord[nnn_AP])))**2+((TOT_A_obs_ord[nnn_AP]*np.sin((np.pi/180.0)*TOT_P_obs_ord[nnn_AP])-(TOT_A_mod_ord[nnn_AP]*np.sin((np.pi/180.0)*TOT_P_mod_ord[nnn_AP])))**2))

      # Root Mean Square difference
      RMSd[comp_idx]=RMSd[comp_idx]+(TOT_A_obs_ord[nnn_AP]*np.cos((np.pi/180.0)*TOT_P_obs_ord[nnn_AP])-(TOT_A_mod_ord[nnn_AP]*np.cos((np.pi/180.0)*TOT_P_mod_ord[nnn_AP])))**2+((TOT_A_obs_ord[nnn_AP]*np.sin((np.pi/180.0)*TOT_P_obs_ord[nnn_AP])-(TOT_A_mod_ord[nnn_AP]*np.sin((np.pi/180.0)*TOT_P_mod_ord[nnn_AP])))**2)

  RMSd[comp_idx]=np.sqrt((1/(2*len(TOT_stations_lab_ord)))*RMSd[comp_idx])

  comp_idx=comp_idx+1


###################### END LOOP ON COMPONENTS
# Close tables
Amp_file.close()
Pha_file.close() 
LinReg_file.close()

comp=['M2','S2','K1','O1','N2','P1','Q1','K2']

for j in range (1,9): # Loop on components
 if tpxo_flag != 0:
  TPXO_AMP=globals()['TPXO_'+comp[j-1]]
  TPXO_PHA=globals()['TPXO_P_'+comp[j-1]]#-P_shift[j-1]#-P_offset[j-1] # rm the tpxo zero offset 
  for i in range (0,len(TPXO_PHA)):
           if TPXO_PHA[i]<0:
              TPXO_PHA[i]=TPXO_PHA[i]+360
           if TPXO_PHA[i]<0:
              TPXO_PHA[i]=TPXO_PHA[i]+360
 if where_box=='Med':
    Tab_filename=path+'tab_'+comp[j-1]+'.txt'
 elif where_box=='AtlBox':
    Tab_filename=path+'tab_'+comp[j-1]+'_AB.txt'

 # TABLE A_mod/A_obs P_mod/P_obs alpha(ratio Amp) dvec Per stz
 Tab_file = open(Tab_filename,"w")
 Tab_file.write('\\begin{table}')
 Tab_file.write('\\footnotesize')
 Tab_file.write('\hspace{-2.5cm}')
 Tab_file.write('\\begin{tabular}{||c||c|c|c|c||c|c|c||c||}')
 Tab_file.write('     \hline')
 Tab_file.write('     \hline')
 Tab_file.write('     TG Name & \multicolumn{4}{||c||}{Amplitudes} & \multicolumn{3}{||c||}{Phases} & Vectorial \\\\')
 Tab_file.write('     \hline')
 Tab_file.write('      & $A_{mod}$ [cm] & $A_{obs}$ [cm] & $A_{tpxo}$ [cm] & $\\alpha$ & P_{mod} [deg] & P_{obs} [deg] &P_{tpxo} [deg]& Vectorial \\\\')
 Tab_file.write('      & & & & & & & & & distance [cm] \\\\')
 Tab_file.write('     \hline')
 Tab_file.write('     \hline')
 for i in range(1,N_stz+1): # Loop on stz
   ALPHA=GLOB_A_mod[i][j]/GLOB_A_obs[i][j]
   print (GLOB_A_mod[i][0],'&',end =" ",file=Tab_file)
   #d_foreman=np.array(d_foreman)
   if tpxo_flag == 0: # TO BE checked!
    print(np.round(np.array(GLOB_A_mod[i][j]),1),'&',np.round(np.array(GLOB_A_obs[i][j]),1),'&',np.round(ALPHA,2),'&',np.round(np.array(GLOB_P_mod[i][j]),1),'&',np.round(np.array(GLOB_P_obs[i][j]),1),'&',np.round(np.array(d_foreman[i][j]),1),'\\\\',file=Tab_file)
   else: # OK
    if j > 4:
       print(np.round(np.array(GLOB_A_mod[i][j]),1),'&',np.round(np.array(GLOB_A_obs[i][j]),1),'&',np.round(np.array(TPXO_AMP[i-1])*100,1),'&',np.round(ALPHA,2),'&',np.round(np.array(GLOB_P_mod[i][j]),1),'&',np.round(np.array(GLOB_P_obs[i][j]),1),'&',np.round(np.array(TPXO_PHA[i-1]),1),'&',np.round(np.array(d_foreman[i][j]),1),'\\\\',file=Tab_file)
    else:
       A_lit1=globals()['TSIMPLIS_'+str(comp[j-1])]
       P_lit1=globals()['TSIMPLIS_P_'+str(comp[j-1])]
       d_lit1=globals()['TSIMPLIS_d_'+str(comp[j-1])]
       A_lit2=globals()['PALMA_'+str(comp[j-1])]
       P_lit2=globals()['PALMA_P_'+str(comp[j-1])]
       d_lit2=globals()['PALMA_d_'+str(comp[j-1])]
       #
       #print ('Prova literature in tab:',comp[j-1],i,GLOB_A_mod[i][0])
       
       #print (A_lit1[i-1],len(A_lit1))
       #print (P_lit1[i-1],len(P_lit1))
       #print (d_lit1[i-1],len(d_lit1))
       #print (A_lit2[i-1],len(A_lit2))
       #print (P_lit2[i-1],len(P_lit2))
       #print (d_lit2[i-1],len(d_lit2))
       #
       print(np.round(np.array(GLOB_A_mod[i][j]),1),'&',np.round(np.array(GLOB_A_obs[i][j]),1),'&',A_lit1[i-1],'/',A_lit2[i-1],'&',np.round(np.array(TPXO_AMP[i-1])*100,1),'&',np.round(ALPHA,2),'&',np.round(np.array(GLOB_P_mod[i][j]),1),'&',np.round(np.array(GLOB_P_obs[i][j]),1),'&',P_lit1[i-1],'/',P_lit2[i-1],'&',np.round(np.array(TPXO_PHA[i-1]),1),'&',np.round(np.array(d_foreman[i][j]),1),'&',d_lit1[i-1],'/',d_lit2[i-1],'\\\\',file=Tab_file)
 Tab_file.write('\hline')
 print ('RMSd &&&&&&&&',RMSd[j-1],'\\\\',file=Tab_file)
 Tab_file.close() 


# Colors array
TOT_color_stz=[]
for idx_color in range (0,howmany_Oarea):
    TOT_color_stz.append('green')
for idx_color in range (howmany_Oarea,howmany_Oarea+howmany_Marea):
    TOT_color_stz.append('orange')
for idx_color in range (howmany_Oarea+howmany_Marea,howmany_Oarea+howmany_Marea+howmany_Garea):
    TOT_color_stz.append('red')
for idx_color in range (howmany_Oarea+howmany_Marea+howmany_Garea,howmany_Oarea+howmany_Marea+howmany_Garea+howmany_Aarea):
    TOT_color_stz.append('blue')
for idx_color in range (howmany_Oarea+howmany_Marea+howmany_Garea+howmany_Aarea,howmany_Oarea+howmany_Marea+howmany_Garea+howmany_Aarea+howmany_Tarea):
    TOT_color_stz.append('tab:cyan')
for idx_color in range (howmany_Oarea+howmany_Marea+howmany_Garea+howmany_Aarea+howmany_Tarea,howmany_Oarea+howmany_Marea+howmany_Garea+howmany_Aarea+howmany_Tarea+howmany_TAarea):
    TOT_color_stz.append('tab:olive')
for idx_color in range (howmany_Oarea+howmany_Marea+howmany_Garea+howmany_Aarea+howmany_Tarea+howmany_TAarea,howmany_Oarea+howmany_Marea+howmany_Garea+howmany_Aarea+howmany_Tarea+howmany_TAarea+howmany_TGarea):
    TOT_color_stz.append('deeppink')
for idx_color in range (howmany_Oarea+howmany_Marea+howmany_Garea+howmany_Aarea+howmany_Tarea+howmany_TAarea+howmany_TGarea,howmany_Oarea+howmany_Marea+howmany_Garea+howmany_Aarea+howmany_Tarea+howmany_TAarea+howmany_TGarea+howmany_Earea):
    TOT_color_stz.append('magenta')


# AMPLITUDES Global plots

plt.figure(figsize=(80,24))
plt.rc('font', size=40)

labels=[ GLOB_A_mod[i][0] for i in range(1,N_stz+1) ]
mod_M2 = [ GLOB_A_mod[j][1] for j in range (1,N_stz+1)]
obs_M2 = [ GLOB_A_obs[j][1] for j in range (1,N_stz+1)]
mod_S2 = [ GLOB_A_mod[j][2] for j in range (1,N_stz+1)]
obs_S2 = [ GLOB_A_obs[j][2] for j in range (1,N_stz+1)]
mod_K1 = [ GLOB_A_mod[j][3] for j in range (1,N_stz+1)]
obs_K1 = [ GLOB_A_obs[j][3] for j in range (1,N_stz+1)]
mod_O1 = [ GLOB_A_mod[j][4] for j in range (1,N_stz+1)]
obs_O1 = [ GLOB_A_obs[j][4] for j in range (1,N_stz+1)]
#
mod_N2 = [ GLOB_A_mod[j][5] for j in range (1,N_stz+1)]
obs_N2 = [ GLOB_A_obs[j][5] for j in range (1,N_stz+1)]
mod_P1 = [ GLOB_A_mod[j][6] for j in range (1,N_stz+1)]
obs_P1 = [ GLOB_A_obs[j][6] for j in range (1,N_stz+1)]
mod_Q1 = [ GLOB_A_mod[j][7] for j in range (1,N_stz+1)]
obs_Q1 = [ GLOB_A_obs[j][7] for j in range (1,N_stz+1)]
mod_K2 = [ GLOB_A_mod[j][8] for j in range (1,N_stz+1)]
obs_K2 = [ GLOB_A_obs[j][8] for j in range (1,N_stz+1)]

x = np.arange(len(labels))  # the label locations
width = 0.45  # the width of the bars

fig,ax=plt.subplots( figsize=(80,24))
#
rects1_M2 = ax.bar(x - width/2, mod_M2, width-0.05, color='#1f77b4', label='M2')
rects2_M2 = ax.bar(x + width/2, obs_M2, width-0.05, color='#1f77b4')

topbottom_M=mod_M2
topbottom_O=obs_M2
rects1_S2 = ax.bar(x - width/2, mod_S2, width-0.05, bottom=topbottom_M, color='#ff7f03',label='S2')
rects2_S2 = ax.bar(x + width/2, obs_S2, width-0.05, bottom=topbottom_O, color='#ff7f03')

zipped_lists_M = zip(topbottom_M,mod_S2)
zipped_lists_O = zip(topbottom_O,obs_S2)
topbottom_M=[x + y for (x, y) in zipped_lists_M] 
topbottom_O=[x + y for (x, y) in zipped_lists_O] 
rects1_K1 = ax.bar(x - width/2, mod_K1, width-0.05, bottom=topbottom_M, color='#2ca02c', label='K1')
rects2_K1 = ax.bar(x + width/2, obs_K1, width-0.05, bottom=topbottom_O, color='#2ca02c')

zipped_lists_M = zip(topbottom_M,mod_K1)
zipped_lists_O = zip(topbottom_O,obs_K1)
topbottom_M=[x + y for (x, y) in zipped_lists_M]
topbottom_O=[x + y for (x, y) in zipped_lists_O]
rects1_O1 = ax.bar(x - width/2, mod_O1, width-0.05, bottom=topbottom_M, color='#d62728', label='O1')
rects2_O1 = ax.bar(x + width/2, obs_O1, width-0.05, bottom=topbottom_O, color='#d62728')
#
zipped_lists_M = zip(topbottom_M,mod_O1)
zipped_lists_O = zip(topbottom_O,obs_O1)
topbottom_M=[x + y for (x, y) in zipped_lists_M]
topbottom_O=[x + y for (x, y) in zipped_lists_O]
rects1_N2 = ax.bar(x - width/2, mod_N2, width-0.05, bottom=topbottom_M, color='#bcdb22', label='N2')
rects2_N2 = ax.bar(x + width/2, obs_N2, width-0.05, bottom=topbottom_O, color='#bcdb22')

zipped_lists_M = zip(topbottom_M,mod_N2)
zipped_lists_O = zip(topbottom_O,obs_N2)
topbottom_M=[x + y for (x, y) in zipped_lists_M]
topbottom_O=[x + y for (x, y) in zipped_lists_O]
rects1_P1 = ax.bar(x - width/2, mod_P1, width-0.05, bottom=topbottom_M, color='#17becf', label='P1')
rects2_P1 = ax.bar(x + width/2, obs_P1, width-0.05, bottom=topbottom_O, color='#17becf')

zipped_lists_M = zip(topbottom_M,mod_P1)
zipped_lists_O = zip(topbottom_O,obs_P1)
topbottom_M=[x + y for (x, y) in zipped_lists_M]
topbottom_O=[x + y for (x, y) in zipped_lists_O]
rects1_Q1 = ax.bar(x - width/2, mod_Q1, width-0.05, bottom=topbottom_M, color='#9467bd', label='Q1')
rects2_Q1 = ax.bar(x + width/2, obs_Q1, width-0.05, bottom=topbottom_O, color='#9467bd')

zipped_lists_M = zip(topbottom_M,mod_Q1)
zipped_lists_O = zip(topbottom_O,obs_Q1)
topbottom_M=[x + y for (x, y) in zipped_lists_M]
topbottom_O=[x + y for (x, y) in zipped_lists_O]
rects1_K2 = ax.bar(x - width/2, mod_K2, width-0.05, bottom=topbottom_M, color='#e377c2', label='K2')
rects2_K2 = ax.bar(x + width/2, obs_K2, width-0.05, bottom=topbottom_O, color='#e377c2')

zipped_lists_M = zip(topbottom_M,mod_K2)
zipped_lists_O = zip(topbottom_O,obs_K2)
topbottom_M=[x + y for (x, y) in zipped_lists_M]
topbottom_O=[x + y for (x, y) in zipped_lists_O]

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('Tidal Amplitudes [cm]')
ax.set_title('Tidal Amplitudes')
ax.set_xticks(x)
ax.set_xticklabels(labels,fontweight='bold')
for xtick, color in zip(ax.get_xticklabels(), TOT_color_stz):
    xtick.set_color(color)
ax.legend()

plt.grid ()

def autolabel1(rects,g_height):
    bin_idx=0
    for rect in rects:
        #height = rect.get_height()
        height = g_height
        ax.annotate('M',
                    xy=(rect.get_x() + rect.get_width() / 2, height[bin_idx]),
                    xytext=(0, 3),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom')
        bin_idx=bin_idx+1

def autolabel2(rects,g_height):
    bin_idx=0
    for rect in rects:
        #height = rect.get_height()
        height = g_height
        ax.annotate('O',
                    xy=(rect.get_x() + rect.get_width() / 2, height[bin_idx]),
                    xytext=(0, 3),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom')
        bin_idx=bin_idx+1

autolabel1(rects1_K2,topbottom_M)
autolabel2(rects2_K2,topbottom_O)

if where_box=='Med':
       plt.savefig(path+'GLOBAL_A_'+dates_lab+'.jpg')
elif where_box=='AtlBox':
       plt.savefig(path+'GLOBAL_A_'+dates_lab+'_AB.jpg')
plt.clf()

# PHASES Global plots

plt.figure(figsize=(80,24))
plt.rc('font', size=40)

labels=[ GLOB_P_mod[i][0] for i in range(1,N_stz+1) ]
mod_M2 = [ GLOB_P_mod[j][1] for j in range (1,N_stz+1)]
obs_M2 = [ GLOB_P_obs[j][1] for j in range (1,N_stz+1)]
mod_S2 = [ GLOB_P_mod[j][2] for j in range (1,N_stz+1)]
obs_S2 = [ GLOB_P_obs[j][2] for j in range (1,N_stz+1)]
mod_K1 = [ GLOB_P_mod[j][3] for j in range (1,N_stz+1)]
obs_K1 = [ GLOB_P_obs[j][3] for j in range (1,N_stz+1)]
mod_O1 = [ GLOB_P_mod[j][4] for j in range (1,N_stz+1)]
obs_O1 = [ GLOB_P_obs[j][4] for j in range (1,N_stz+1)]
#
mod_N2 = [ GLOB_P_mod[j][5] for j in range (1,N_stz+1)]
obs_N2 = [ GLOB_P_obs[j][5] for j in range (1,N_stz+1)]
mod_P1 = [ GLOB_P_mod[j][6] for j in range (1,N_stz+1)]
obs_P1 = [ GLOB_P_obs[j][6] for j in range (1,N_stz+1)]
mod_Q1 = [ GLOB_P_mod[j][7] for j in range (1,N_stz+1)]
obs_Q1 = [ GLOB_P_obs[j][7] for j in range (1,N_stz+1)]
mod_K2 = [ GLOB_P_mod[j][8] for j in range (1,N_stz+1)]
obs_K2 = [ GLOB_P_obs[j][8] for j in range (1,N_stz+1)]

x = np.arange(len(labels))  # the label locations
width = 0.45  # the width of the bars

fig,ax=plt.subplots( figsize=(80,24))
#
rects1_M2 = ax.bar(x - width/2, mod_M2, width-0.05, color='#1f77b4', label='M2')
rects2_M2 = ax.bar(x + width/2, obs_M2, width-0.05, color='#1f77b4')

topbottom_M=mod_M2
topbottom_O=obs_M2
rects1_S2 = ax.bar(x - width/2, mod_S2, width-0.05, bottom=topbottom_M, color='#ff7f03',label='S2')
rects2_S2 = ax.bar(x + width/2, obs_S2, width-0.05, bottom=topbottom_O, color='#ff7f03')

zipped_lists_M = zip(topbottom_M,mod_S2)
zipped_lists_O = zip(topbottom_O,obs_S2)
topbottom_M=[x + y for (x, y) in zipped_lists_M]
topbottom_O=[x + y for (x, y) in zipped_lists_O]
rects1_K1 = ax.bar(x - width/2, mod_K1, width-0.05, bottom=topbottom_M, color='#2ca02c', label='K1')
rects2_K1 = ax.bar(x + width/2, obs_K1, width-0.05, bottom=topbottom_O, color='#2ca02c')

zipped_lists_M = zip(topbottom_M,mod_K1)
zipped_lists_O = zip(topbottom_O,obs_K1)
topbottom_M=[x + y for (x, y) in zipped_lists_M]
topbottom_O=[x + y for (x, y) in zipped_lists_O]
rects1_O1 = ax.bar(x - width/2, mod_O1, width-0.05, bottom=topbottom_M, color='#d62728', label='O1')
rects2_O1 = ax.bar(x + width/2, obs_O1, width-0.05, bottom=topbottom_O, color='#d62728')
#
zipped_lists_M = zip(topbottom_M,mod_O1)
zipped_lists_O = zip(topbottom_O,obs_O1)
topbottom_M=[x + y for (x, y) in zipped_lists_M]
topbottom_O=[x + y for (x, y) in zipped_lists_O]
rects1_N2 = ax.bar(x - width/2, mod_N2, width-0.05, bottom=topbottom_M, color='#bcdb22', label='N2')
rects2_N2 = ax.bar(x + width/2, obs_N2, width-0.05, bottom=topbottom_O, color='#bcdb22')

zipped_lists_M = zip(topbottom_M,mod_N2)
zipped_lists_O = zip(topbottom_O,obs_N2)
topbottom_M=[x + y for (x, y) in zipped_lists_M]
topbottom_O=[x + y for (x, y) in zipped_lists_O]
rects1_P1 = ax.bar(x - width/2, mod_P1, width-0.05, bottom=topbottom_M, color='#17becf', label='P1')
rects2_P1 = ax.bar(x + width/2, obs_P1, width-0.05, bottom=topbottom_O, color='#17becf')

zipped_lists_M = zip(topbottom_M,mod_P1)
zipped_lists_O = zip(topbottom_O,obs_P1)
topbottom_M=[x + y for (x, y) in zipped_lists_M]
topbottom_O=[x + y for (x, y) in zipped_lists_O]
rects1_Q1 = ax.bar(x - width/2, mod_Q1, width-0.05, bottom=topbottom_M, color='#9467bd', label='Q1')
rects2_Q1 = ax.bar(x + width/2, obs_Q1, width-0.05, bottom=topbottom_O, color='#9467bd')

zipped_lists_M = zip(topbottom_M,mod_Q1)
zipped_lists_O = zip(topbottom_O,obs_Q1)
topbottom_M=[x + y for (x, y) in zipped_lists_M]
topbottom_O=[x + y for (x, y) in zipped_lists_O]
rects1_K2 = ax.bar(x - width/2, mod_K2, width-0.05, bottom=topbottom_M, color='#e377c2', label='K2')
rects2_K2 = ax.bar(x + width/2, obs_K2, width-0.05, bottom=topbottom_O, color='#e377c2')

zipped_lists_M = zip(topbottom_M,mod_K2)
zipped_lists_O = zip(topbottom_O,obs_K2)
topbottom_M=[x + y for (x, y) in zipped_lists_M]
topbottom_O=[x + y for (x, y) in zipped_lists_O]

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('Tidal Phases [deg]')
ax.set_title('Tidal Phases')
ax.set_xticks(x)
ax.set_xticklabels(labels,fontweight='bold')
for xtick, color in zip(ax.get_xticklabels(), TOT_color_stz):
    xtick.set_color(color)
ax.legend()

plt.grid ()

def autolabel1(rects,g_height):
    bin_idx=0
    for rect in rects:
        #height = rect.get_height()
        height = g_height
        ax.annotate('M',
                    xy=(rect.get_x() + rect.get_width() / 2, height[bin_idx]),
                    xytext=(0, 3),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom')
        bin_idx=bin_idx+1

def autolabel2(rects,g_height):
    bin_idx=0
    for rect in rects:
        #height = rect.get_height()
        height = g_height
        ax.annotate('O',
                    xy=(rect.get_x() + rect.get_width() / 2, height[bin_idx]),
                    xytext=(0, 3),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom')
        bin_idx=bin_idx+1

autolabel1(rects1_K2,topbottom_M)
autolabel2(rects2_K2,topbottom_O)
if where_box=='Med':
       plt.savefig(path+'GLOBAL_P_'+dates_lab+'.jpg')
elif where_box=='AtlBox':
       plt.savefig(path+'GLOBAL_P_'+dates_lab+'_AB.jpg')
plt.clf()



# GLOBAL DISTANCES PLOT

plt.figure(figsize=(80,24))
plt.rc('font', size=40)

labels=[ d_foreman[i][0] for i in range(1,N_stz+1) ]
mod_M2 = [ d_foreman[j][1] for j in range (1,N_stz+1)]
mod_S2 = [ d_foreman[j][2] for j in range (1,N_stz+1)]
mod_K1 = [ d_foreman[j][3] for j in range (1,N_stz+1)]
mod_O1 = [ d_foreman[j][4] for j in range (1,N_stz+1)]
#
mod_N2 = [ d_foreman[j][5] for j in range (1,N_stz+1)]
mod_P1 = [ d_foreman[j][6] for j in range (1,N_stz+1)]
mod_Q1 = [ d_foreman[j][7] for j in range (1,N_stz+1)]
mod_K2 = [ d_foreman[j][8] for j in range (1,N_stz+1)]

# Compute and print statistics on vectorial distances
if flag_15stats == 1:
   # Table for Comp with literature on 15 TGs
   Lit_file = open(path+"stats_15lit.txt","w")
   print ('PROVA:', Lit_file)
   print (labels,file=Lit_file)
   print ('M2: ', np.mean(mod_M2), np.max(mod_M2),file=Lit_file)
   print (mod_M2,file=Lit_file)
   print ('S2: ', np.mean(mod_S2), np.max(mod_S2),file=Lit_file)
   print (mod_S2,file=Lit_file)
   print ('K1: ', np.mean(mod_K1), np.max(mod_K1),file=Lit_file)
   print (mod_K1,file=Lit_file)
   print ('O1: ', np.mean(mod_O1), np.max(mod_O1),file=Lit_file)
   print (mod_O1,file=Lit_file)

# Plot vectorial distances
x = np.arange(len(labels))  # the label locations
width = 0.45*2  # the width of the bars

fig,ax=plt.subplots( figsize=(80,24))
#
rects1_M2 = ax.bar(x, mod_M2, width-0.05, color='#1f77b4', label='M2')

topbottom_M=mod_M2
rects1_S2 = ax.bar(x, mod_S2, width-0.05, bottom=topbottom_M, color='#ff7f03',label='S2')

zipped_lists_M = zip(topbottom_M,mod_S2)
topbottom_M=[x + y for (x, y) in zipped_lists_M]
rects1_K1 = ax.bar(x, mod_K1, width-0.05, bottom=topbottom_M, color='#2ca02c', label='K1')

zipped_lists_M = zip(topbottom_M,mod_K1)
topbottom_M=[x + y for (x, y) in zipped_lists_M]
rects1_O1 = ax.bar(x, mod_O1, width-0.05, bottom=topbottom_M, color='#d62728', label='O1')
#
zipped_lists_M = zip(topbottom_M,mod_O1)
topbottom_M=[x + y for (x, y) in zipped_lists_M]
rects1_N2 = ax.bar(x, mod_N2, width-0.05, bottom=topbottom_M, color='#bcdb22', label='N2')

zipped_lists_M = zip(topbottom_M,mod_N2)
topbottom_M=[x + y for (x, y) in zipped_lists_M]
rects1_P1 = ax.bar(x, mod_P1, width-0.05, bottom=topbottom_M, color='#17becf', label='P1')

zipped_lists_M = zip(topbottom_M,mod_P1)
topbottom_M=[x + y for (x, y) in zipped_lists_M]
rects1_Q1 = ax.bar(x, mod_Q1, width-0.05, bottom=topbottom_M, color='#9467bd', label='Q1')

zipped_lists_M = zip(topbottom_M,mod_Q1)
topbottom_M=[x + y for (x, y) in zipped_lists_M]
rects1_K2 = ax.bar(x, mod_K2, width-0.05, bottom=topbottom_M, color='#e377c2', label='K2')

zipped_lists_M = zip(topbottom_M,mod_K2)
topbottom_M=[x + y for (x, y) in zipped_lists_M]

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('Vectorial Differences [cm]')
ax.set_title('Vectorial Differences per TG per tidal component [Foreman et al. 1993]')
ax.set_xticks(x)
ax.set_xticklabels(labels,fontweight='bold')
for xtick, color in zip(ax.get_xticklabels(), TOT_color_stz):
    xtick.set_color(color)
ax.legend()

plt.grid ()

if where_box=='Med':
       plt.savefig(path+'GLOBAL_d_Foreman_'+dates_lab+'.jpg')
elif where_box=='AtlBox':
       plt.savefig(path+'GLOBAL_d_Foreman_'+dates_lab+'_AB.jpg')
plt.clf()


# GLOBAL RMS DISTANCES PLOT

plt.figure(figsize=(20,10))
plt.rc('font', size=20)

labels=[ d_foreman[0][j] for j in range(1,N_comp+1) ]

x = np.arange(len(labels))  # the label locations
width = 0.45*2  # the width of the bars

fig,ax=plt.subplots( figsize=(20,10))
#
comp_color=['#1f77b4','#ff7f03','#2ca02c','#d62728','#bcdb22','#17becf','#9467bd','#e377c2']
rects1 = ax.bar(x, RMSd, width-0.05,color=comp_color)


# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('RMSd [cm]')
ax.set_title('Root Mean Square differences - '+where_box)
ax.set_xticks(x)
ax.set_xticklabels(labels,fontweight='bold')
#for xtick, color in zip(ax.get_xticklabels(), TOT_color_stz):
#    xtick.set_color(color)
#ax.legend()

plt.grid ()


def autolabel(rects):
    bin_idx=0
    for rect in rects:
        height = rect.get_height()
        #height = g_height
        ax.annotate(round(RMSd[bin_idx],1),
                    xy=(rect.get_x() + rect.get_width() / 2, height ),
                    xytext=(0, 3),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom')
        bin_idx=bin_idx+1

autolabel(rects1)
if where_box=='Med':
       plt.savefig(path+'GLOBAL_RMSd_'+dates_lab+'.jpg')
elif where_box=='AtlBox':
       plt.savefig(path+'GLOBAL_RMSd_'+dates_lab+'_AB.jpg')
plt.clf()


print ('Output path: ',path)

import glob
import sys
import xarray as xr
import netCDF4 as nc
import numpy as np
from scipy.signal import periodogram
import heapq
from functools import partial
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from scipy.signal import detrend
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks
mpl.use('Agg')

########
# Inputs and outputs

start_date = "20150101" #"20200101" #"20150201" #"20160101"
end_date = "20170104" #"20250101" #"20160901" #"20240101"

#all_files = sorted(glob.glob("/work/cmcc/ag15419/exp/fix_mfseas9_longrun_surge_2NT_AB/EXP00/20*/wind/out/ecmwf_y20*.nc"))
#all_files = sorted(glob.glob("/work/cmcc/ag15419/exp/fix_mfseas9_longrun_hmslp_2NT_AB/EXP00/20*/wind/out/ecmwf_y20*.nc"))
all_files = sorted(glob.glob("/work/cmcc/ag15419/exp/fix_mfseas9_longrun_barotropic_final22/EXP00/20*/wind/out/ecmwf_y20*.nc"))

# Exp tag
Med_reg=str(sys.argv[3])
exp='bt_ecmwf_'+Med_reg

# Lat and lon indexes
lat_idx = int(sys.argv[2]) #72 #138 #358 #360
lon_idx = int(sys.argv[1]) #1127 #331 #744 #746

# Model time step in seconds
dt = 3600*6 #450 #90 #60*60 #50 #10 #3*60

# Number of modes to analyze
n_modes = 12

# Minimum peak amplitude ; min,max width ; min distance between peaks to detect peaks in Amp plots (meters,hours, points respectively)
#amp_peak_height=0.000001 #0.0001
#amp_peak_width=(0, 20)
#amp_peak_distance=71 # 31 11 3
#tolerance = 1  # Minimum peakdistance between peaks in hours

# Flag and threshold [h] for filtering the spectrum the threshold is also used as plot minimum 
flag_filter='true'
th_filter=40 #*240

# Flag for spectrum detrending:
flag_detrend='false'

# Flag for Gaussian smoothing of the spectrum: true, false or plot (to use the original spt but add the plot of the smoothed spt)
flag_smooth='plot'
sigma=15 #4 11
#def moving_average(data, window_size):
#    return np.convolve(data, np.ones(window_size) / window_size, mode='same')
#window_size=11

###################
# Select the period

infile = []
visited_dates = set() 
for f in all_files:
    parts = f.split("/")
    date_folder = parts[7]
    if date_folder not in visited_dates:
        infile.append(f)  
        visited_dates.add(date_folder)

# Initialize SSH time series
ssh_ts_all = []

# Read data from NetCDF files
grid_info = False
for nc2open in infile:
    print('Processing:', nc2open)
    model = nc.Dataset(nc2open, 'r')
    ssh_ts = np.array(model.variables['msl'][:, lat_idx, lon_idx])
    ssh_ts_all = np.concatenate((ssh_ts_all, ssh_ts))

    if not grid_info:
        lats = 0 #np.round(model.variables['nav_lat'][lat_idx, lon_idx],2)
        lons = 0 #np.round(model.variables['nav_lon'][lat_idx, lon_idx],2)
        grid_info = True
        print ('I am working on',lats,lons)

    #print(f"Total number of points in the SSH time series: {len(ssh_ts_all)}")
    model.close()

print ("I am workign on period",start_date,"-",end_date," Freq:",dt,"s Num of inputs:",len(ssh_ts_all))

# Convert SSH time series to NumPy array and remove NaNs
time_series_point = np.array(ssh_ts_all)
#print ('Check input time series:',time_series_point.shape,time_series_point)
valid_indices = np.logical_not(np.isnan(time_series_point))
time_series_clean = time_series_point[valid_indices]
# Conversion Pa->hPa
time_series_clean=time_series_clean/100.0
#print ('Check clean time series:',time_series_clean.shape,time_series_clean)

#### SPECTRUM ANALYSIS

# Compute FFT
spt_len = len(time_series_clean)
print('Time series values:', spt_len)
fft = np.fft.fft(time_series_clean)
freq = np.fft.fftfreq(spt_len, d=dt)

# Select only positive frequencies (first half of FFT, excluding Nyquist frequency)
half_spt_len = spt_len // 2
freq_positive = freq[:half_spt_len]  # Only positive frequencies
fft_positive = fft[:half_spt_len]  # Only corresponding FFT values

# Ensure we are not including zero frequency
mask = freq_positive > 0
freq_positive = freq_positive[mask]
fft_positive = fft_positive[mask]
#freq_positive = freq_positive[freq_positive > 0]
#fft_positive = fft_positive[:len(freq_positive)]  # Match the length of fft_positive with freq_positive
#print(f"Lunghezza freq_positive: {len(freq_positive)}")
#print(f"Lunghezza fft_positive: {len(fft_positive)}")

# Compute Power Spectrum Density (normalized)
spt = (np.abs(fft_positive) ** 2) / spt_len #**2

# Compute Periods in hours
periods = 1 / freq_positive / 3600  # Convert periods to hours
# Check that periods span the expected range
#print(f"Periods (hours) range from {periods[0]:.2f} h to {periods[-1]:.2f} h.")

# Compute Amplitudes (correct scaling)
amplitudes = (2 / spt_len) * np.abs(fft_positive)

if flag_filter=='true':
   print ('Filter = true')
   # Apply high-pass filter: Set frequencies below the threshold to zero
   high_pass_threshold = 1 / (th_filter * 3600)  # Corresponding to 48 hours in Hz
   fft_positive[freq_positive < high_pass_threshold] = 0  # Filter out frequencies below the threshold

   # Recompute the power spectrum and amplitudes after filtering
   spt = (np.abs(fft_positive) ** 2) / spt_len #**2
   amplitudes = (2 / spt_len) * np.abs(fft_positive)

# Smooth
if flag_smooth == 'true':
   print ('Smooth = true')
   spt_smooth = gaussian_filter1d(spt, sigma=sigma)
   #spt_smooth = moving_average(spt, window_size)
   amp_smooth = gaussian_filter1d(amplitudes, sigma=sigma)
else:
   spt_smooth = spt
   amp_smooth = amplitudes
   if flag_smooth == 'plot':
      spt_smooth_2plot = gaussian_filter1d(spt, sigma=sigma)
      #spt_smooth_2plot = moving_average(spt, window_size)
      amp_smooth_2plot = gaussian_filter1d(amplitudes, sigma=sigma)

# Detrend
if flag_detrend == 'true':
   print ('Detrend = true')
   spt_det = detrend(spt_smooth)
   #detrend_window_size = int(spt_len/2)
   #for i in range(0, len(spt_smooth), detrend_window_size):
   #    spt_det[i:i+window_size] = detrend(spt_smooth[i:i+window_size])
else:
   spt_det = spt

# Found peaks in spt and in amplitude:
#peaks, _ = find_peaks(spt_smooth,prominence=amp_peak_height,width=amp_peak_width,distance=amp_peak_distance)
#peak_frequencies = freq_positive[peaks]
#peak_amplitudes = spt_smoothed[peaks]

#amp_peaks, _ = find_peaks(amp_smooth,prominence=amp_peak_height,width=amp_peak_width,distance=amp_peak_distance)
#amp_peaks, _ = find_peaks(amp_smooth,distance=amp_peak_distance)
amp_peaks, _ = find_peaks(amp_smooth)
amp_peak_frequencies = freq_positive[amp_peaks]
amp_peak_amplitudes = amp_smooth[amp_peaks]

# Order based on amplitudes
sorted_indices_peak_amp = np.argsort(amp_peak_amplitudes)[::-1]
amp_peak_amplitudes_sorted = amp_peak_amplitudes[sorted_indices_peak_amp]
amp_peak_frequencies_sorted = amp_peak_frequencies[sorted_indices_peak_amp]

# Remove close modes:
period_sorted = 1/amp_peak_frequencies_sorted/3600  # periodi dei modi
amplitude_sorted = amp_peak_amplitudes_sorted  # ampiezze dei modi

print("Ini Periods:", period_sorted)

final_periods = []
final_amplitudes = []

for i in range(len(period_sorted)):
    if period_sorted[i] > 30:
        tolerance = 10  
    elif 25 <= period_sorted[i] <= 30:
        tolerance = 4  
    elif 12 <= period_sorted[i] < 25:
        tolerance = 2  
    elif 6 <= period_sorted[i] < 12:
        tolerance = 1  
    else:
        tolerance = 0.1  

    keep_mode = True  # Presumo di tenere il modo, salvo verifica successiva

    for j in range(len(final_periods)):
        if abs(period_sorted[i] - final_periods[j]) <= tolerance:
            if amplitude_sorted[i] <= final_amplitudes[j]:
                keep_mode = False  # Non aggiungere il modo
                break
            else:
                final_periods[j] = period_sorted[i]
                final_amplitudes[j] = amplitude_sorted[i]
                keep_mode = False  # Non aggiungere di nuovo lo stesso modo
                break
    
    if keep_mode:
        final_periods.append(period_sorted[i])
        final_amplitudes.append(amplitude_sorted[i])


#for i in range(len(period_sorted)):
#    keep_mode = True
#
#    for j in range(len(final_periods)):
#        if abs(period_sorted[i] - final_periods[j]) <= tolerance:
#            if amplitude_sorted[i] <= final_amplitudes[j]:
#                keep_mode = False  
#                break
#            else:
#                final_periods[j] = period_sorted[i]
#                final_amplitudes[j] = amplitude_sorted[i]
#                keep_mode = False  
#                break
#    
#    if keep_mode:
#        final_periods.append(period_sorted[i])
#        final_amplitudes.append(amplitude_sorted[i])

print("Final Periods:", final_periods)
#print("Final Amplitudes:", final_amplitudes)
amp_peak_amplitudes_sorted=np.array(final_amplitudes)
#amp_peak_frequencies_sorted=1/(np.array(final_periods)*3600)
amp_peak_period_sorted=np.array(final_periods)


# Time array for SSH plot (convert to hours)
ssh_time = np.arange(0, spt_len * dt, dt) / 3600  

# Select the main spectral peaks based on spectral density
top_indices_spt_det = np.argsort(spt_det)[-n_modes:]  # Get indices of top n_modes power values
sorted_indices_spt_det = np.argsort(spt_det[top_indices_spt_det])[::-1]  # Sort by descending power

# Extract the corresponding frequencies, periods, and amplitudes for density-based selection
top_freq_positive_spt_det = freq_positive[top_indices_spt_det][sorted_indices_spt_det]
top_periods_spt_det = periods[top_indices_spt_det][sorted_indices_spt_det]
top_amplitudes_spt_det = amplitudes[top_indices_spt_det][sorted_indices_spt_det]

# Print the main spectral modes based on spectral density
#print("Main Spectral Modes (based on power spectrum density):")
#for i in range(n_modes):
#    print(f"Mode {i+1}: Period = {top_periods_spt_det[i]:.2f} h, Freq = {top_freq_positive_spt_det[i]:.6f} Hz, Amp = {top_amplitudes_spt_det[i]:.3f} m")

# Now select the main modes based on amplitude
n_valid = min(len(amplitudes), len(freq_positive))  # Ensure valid index range
top_indices_amp = np.argpartition(amplitudes[:n_valid], -n_modes)[-n_modes:]  # Select indices of top amplitudes
sorted_indices_amp = np.argsort(amplitudes[top_indices_amp])[::-1]  # Sort by descending amplitude

# Extract the corresponding frequencies, periods, and amplitudes for amplitude-based selection
top_freq_positive_amp = freq_positive[top_indices_amp][sorted_indices_amp]
top_periods_amp = periods[top_indices_amp][sorted_indices_amp]
top_amplitudes_amp = amplitudes[top_indices_amp][sorted_indices_amp]

# Print the main spectral modes based on amplitude
#print("Main Spectral Modes (based on amplitude):")
#for i in range(n_modes):
#    print(f"Mode {i+1}: Period = {top_periods_amp[i]:.2f} h, Freq = {top_freq_positive_amp[i]:.6f} Hz, Amp = {top_amplitudes_amp[i]:.3f} m")

#######################
# PLOT SSH
# Plot the whole period
plt.figure(figsize=(18, 8))
plt.rc('font', size=20)
plt.title(f'MSL Pressure at lat='+str(lats)+' lon='+str(lons)+' '+Med_reg)
plt.plot(ssh_time, time_series_clean, '-',linewidth=2, label=f'MSLP at lat='+str(lats)+' lon='+str(lons))
plt.xlabel('Time (h)')
plt.ylabel('MSLP (hPa)')
plt.grid()
plt.legend(loc='upper right')
plt.savefig(f'mslp_{lat_idx}_{lon_idx}_{exp}.png')

# Plot the last 24h
#plt.figure(figsize=(18, 8))
#plt.rc('font', size=20)
#plt.title(f'SSH at lat='+str(lats)+' lon='+str(lons)+' '+Med_reg)
#plt.plot(ssh_time, time_series_clean, '-',linewidth=4, label=f'SSH at lat='+str(lats)+' lon='+str(lons))
#plt.xlim(0,96)
#plt.xlabel('Time (h)')
#plt.ylabel('SSH (m)')
#plt.grid()
#plt.legend(loc='upper right')
#plt.savefig(f'ssh_24h_{lat_idx}_{lon_idx}_{exp}.png') 

# PLOT POWER SPECTRUM
#plt.figure(figsize=(15, 9))
#plt.title(f'Power Spectrum at lat='+str(lats)+' lon='+str(lons)+' '+Med_reg)
#plt.loglog(periods, spt, marker='o', linestyle='-', label='Power Spectrum')
#if flag_smooth == 'true':
#   plt.loglog(periods, spt_smooth, marker='o', linestyle='-', label='Smoothed Power Spectrum')
#plt.xlabel('Period (h)')
#plt.ylabel('Power Spectrum')
#plt.axvline(24, color='black', linestyle='-')
#plt.axvline(12, color='black', linestyle='-')
#plt.axvline(6, color='black', linestyle='-')
#
## Mark the main modes based on spectral density
#for i in range(n_modes):
#    plt.axvline(top_periods_spt_det[i], color='red', linestyle='--', label=f'SPT Mode {i+1} (T = {top_periods_spt_det[i]:.2f} h, Amp = {top_amplitudes_spt_det[i]:.3f} m)')
#
## Mark the main modes based on amplitude
#for i in range(n_modes):
#    plt.axvline(top_periods_amp[i], color='blue', linestyle='--', label=f'Amp Mode {i+1} (T = {top_periods_amp[i]:.2f} h, Amp = {top_amplitudes_amp[i]:.3f} m)')
#
## Mark the main modes based on peak finder
#for i in range(0,len(peak_frequencies)):
#    plt.axvline(1/peak_frequencies[i]/3600, color='green',linestyle='--')
#
#plt.xlim(th_filter,0.5)
#plt.grid()
#plt.legend()
#plt.savefig(f'spt_{lat_idx}_{lon_idx}_{exp}.png')
#
## Non-log spt plot
#plt.figure(figsize=(15, 9))
#plt.title(f'Power Spectrum at lat={lats} lon={lons}'+' '+Med_reg)
#plt.loglog(periods, spt, marker='o', linestyle='-', label='Power Spectrum')
#if flag_smooth == 'true':
#   plt.loglog(periods, spt_smooth, marker='o', linestyle='-', label='Smoothed Power Spectrum')
#plt.xlabel('Period (h)')
#plt.ylabel('Power Spectrum')
#plt.axvline(24, color='black', linestyle='-')
#plt.axvline(12, color='black', linestyle='-')
#plt.axvline(6, color='black', linestyle='-')
#
## Mark the main modes based on spectral density
#for i in range(n_modes):
#    plt.axvline(top_periods_spt_det[i], color='red', linestyle='--', label=f'SPT Mode {i+1} (T = {top_periods_spt_det[i]:.2f} h, Amp = {top_amplitudes_spt_det[i]:.3f} m)')
#
## Mark the main modes based on amplitude
#for i in range(n_modes):
#    plt.axvline(top_periods_amp[i], color='blue', linestyle='--', label=f'Amp Mode {i+1} (T = {top_periods_amp[i]:.2f} h, Amp = {top_amplitudes_amp[i]:.3f} m)')
#
#plt.xlim(th_filter,0.5)
#plt.yscale('linear')
#plt.grid()
#plt.legend()
#plt.savefig(f'spt_nolog_{lat_idx}_{lon_idx}_{exp}.png')
#
## Peaks plot
plt.figure(figsize=(18, 9))
plt.rc('font', size=20)
plt.title(f'Power Spectrum at lat='+str(lats)+' lon='+str(lons)+' '+Med_reg)
plt.loglog(periods, spt_det, marker='o', linestyle='-', label='Power Spectrum')
if flag_smooth == 'true':
   plt.loglog(periods, spt_smooth, marker='o', linestyle='-', label='Smoothed Power Spectrum')
elif flag_smooth == 'plot':
   plt.loglog(periods, spt_smooth_2plot, marker='o', linestyle='-', label='Smoothed Power Spectrum')
plt.xlabel('Period (h)')
plt.ylabel('Power Spectrum')

plt.axvline(24, color='black', linestyle=':',linewidth=4)
plt.axvline(12, color='black', linestyle=':',linewidth=4)
plt.axvline(6, color='black', linestyle=':',linewidth=4)

#
## Mark the main modes based on spectral density
#for i in range(n_modes):
#    plt.axvline(top_periods_spt_det[i], color='red', linestyle='--', label=f'SPT Mode {i+1} (T = {top_periods_spt_det[i]:.2f} h, Amp = {top_amplitudes_spt_det[i]:.3f} m)')
#
## Mark the main modes based on amplitude
#for i in range(n_modes):
#    plt.axvline(top_periods_amp[i], color='blue', linestyle='--', label=f'Amp Mode {i+1} (T = {top_periods_amp[i]:.2f} h, Amp = {top_amplitudes_amp[i]:.3f} m)')
#
plt.xlim(th_filter,0.5)
plt.text(24,plt.ylim()[0],'24', ha='center', va='top')
plt.text(12,plt.ylim()[0],'12', ha='center', va='top')
plt.text(6,plt.ylim()[0],'6', ha='center', va='top')

plt.grid()
plt.legend(loc='upper right')
plt.savefig(f'mslp_spt_det_{lat_idx}_{lon_idx}_{exp}.png')

# Amplitude plots and tables
plt.figure(figsize=(27, 18))
ax = plt.subplot(111)
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.75, box.height])
plt.title(f'Modes amplitudes at lat='+str(lats)+' lon='+str(lons)+' '+Med_reg)

#plt.axvline(24, color='black', linestyle=':',linewidth=4)
#plt.axvline(12, color='black', linestyle=':',linewidth=4)
#plt.axvline(6, color='black', linestyle=':',linewidth=4)

# Mark the main modes based on peak finder
mode_colors = plt.cm.rainbow(np.linspace(0, 1, n_modes)) #len(amp_peak_frequencies)))
for i in range(0,n_modes): #len(amp_peak_frequencies)):
    try:
       plt.axvline(amp_peak_period_sorted[i], color=mode_colors[i],linestyle='--',linewidth=4,label=f'Mode {i} (T={amp_peak_period_sorted[i]:.2f} h, Amp={amp_peak_amplitudes_sorted[i]:.3f} hPa)')
       #plt.text(1/amp_peak_frequencies_sorted[i]/3600, plt.ylim()[0] - 0.1, f'{1/amp_peak_frequencies_sorted[i]/3600}', ha='center', va='top')
    except:
       print ('Nan')
f_Nyq=dt*2/3600
plt.loglog(periods[periods>f_Nyq], amplitudes[periods>f_Nyq], marker='o', linestyle='-',linewidth=4,color='black', label='Modes Amplitudes')
if flag_smooth == 'true':
   plt.loglog(periods[periods>f_Nyq], amp_smooth[periods>f_Nyq], marker='o', linestyle='-',linewidth=4,color='tab:green', label='Smoothed Modes Amplitudes')
elif flag_smooth == 'plot':
   plt.loglog(periods[periods>f_Nyq], amp_smooth_2plot[periods>f_Nyq], marker='o', linestyle='-',linewidth=4,color='tab:green', label='Smoothed Modes Amplitudes')
plt.xlabel('Period (h)')
plt.ylabel('Mode Amplitude (hPa)')
plt.xlim(th_filter-1,dt*1/3600)
#plt.ylim(0.0000001,0.1)
plt.text(24,plt.ylim()[0],'24', ha='center', va='top')
plt.text(12,plt.ylim()[0],'12', ha='center', va='top')
plt.text(6,plt.ylim()[0],'6', ha='center', va='top')
plt.grid()
plt.legend(loc='upper right') #'center left', bbox_to_anchor=(1, 0.5))

# Aggiunta della tabella
table_data = []
for i in range(n_modes):
    try:
       table_data.append([f'Mode {i}', f'{amp_peak_period_sorted[i]:.2f} h', f'{amp_peak_amplitudes_sorted[i]:.3f} hPa'])
    except:
       print ('No mode',i)

col_labels = ['Mode', 'Period (h)', 'Amplitude (hPa)']

table_ax = plt.table(cellText=table_data,
                     colLabels=col_labels,
                     loc='bottom',
                     cellLoc='center',
                     bbox=[0.15, -0.4, 0.7, 0.3])  # bbox controlla la posizione e la dimensione della tabella

table_ax.auto_set_font_size(False)
table_ax.scale(1.2, 1.2)

for i, label in enumerate(col_labels):
    cell = table_ax.get_celld()[(0, i)]  # Ottieni la cella della riga 0, colonna i
    cell.set_facecolor('#d3d3d3')  # Grigio chiaro
    cell.set_text_props(weight='bold')  # Grassetto per le intestazioni
    
for key, cell in table_ax.get_celld().items():
    cell.set_height(0.1)  # Altezza personalizzata
    cell.set_fontsize(16)

plt.tight_layout()
plt.savefig(f'mslp_amp_{lat_idx}_{lon_idx}_{exp}.png')

# Amp no log plot
plt.figure(figsize=(27, 9))
ax = plt.subplot(111)
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.75, box.height])
plt.title(f'Modes amplitudes at lat='+str(lats)+' lon='+str(lons))

#plt.axvline(24, color='black', linestyle=':',linewidth=4)
#plt.axvline(12, color='black', linestyle=':',linewidth=4)
#plt.axvline(6, color='black', linestyle=':',linewidth=4)

# Mark the main modes based on peak finder
mode_colors = plt.cm.rainbow(np.linspace(0, 1, n_modes)) #len(amp_peak_frequencies)))
for i in range(0,n_modes): #len(amp_peak_frequencies)):
    try:
       plt.axvline(amp_peak_period_sorted[i], color=mode_colors[i],linestyle='--',linewidth=4,label=f'Mode {i} (T={amp_peak_period_sorted[i]:.2f} h, Amp={amp_peak_amplitudes_sorted[i]:.3f} hPa)')
       #plt.text(1/amp_peak_frequencies_sorted[i]/3600, plt.ylim()[0] - 0.1, f'{1/amp_peak_frequencies_sorted[i]/3600}', ha='center', va='top')
    except:
       print ('Nan')
f_Nyq=dt*2/3600
plt.loglog(periods[periods>f_Nyq], amplitudes[periods>f_Nyq], marker='o', linestyle='-',linewidth=4,color='black', label='Modes Amplitudes')
if flag_smooth == 'true':
   plt.loglog(periods[periods>f_Nyq], amp_smooth[periods>f_Nyq], marker='o', linestyle='-',linewidth=4,color='tab:green', label='Smoothed Modes Amplitudes')
elif flag_smooth == 'plot':
   plt.loglog(periods[periods>f_Nyq], amp_smooth_2plot[periods>f_Nyq], marker='o', linestyle='-',linewidth=4,color='tab:green', label='Smoothed Modes Amplitudes')
plt.xlabel('Period (h)')
plt.ylabel('Mode Amplitude (hPa)')
plt.xlim(th_filter-1,dt*1/3600)
##plt.ylim(0.0,0.1)
#plt.ylim(0.0000001,0.05)

plt.text(24,plt.ylim()[0],'24', ha='center', va='top')
plt.text(12,plt.ylim()[0],'12', ha='center', va='top')
plt.text(6,plt.ylim()[0],'6', ha='center', va='top')

plt.grid()
plt.yscale('linear')
#plt.xscale('linear')
plt.legend(loc='upper right') #'center left', bbox_to_anchor=(1, 0.5))
plt.tight_layout()
plt.savefig(f'mslp_amp_nolog_{lat_idx}_{lon_idx}_{exp}.png')


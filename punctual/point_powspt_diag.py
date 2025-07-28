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
from point_ini import *
mpl.use('Agg')

########

# Exp tag
Med_reg=str(sys.argv[3])
exp=tag+Med_reg

# Lat and lon indexes
lat_idx = int(sys.argv[2]) 
lon_idx = int(sys.argv[1]) 

all_files=sorted(glob.glob(file_template))
print ('all_files',all_files)

###################
# Select the period
infile = []
for f in all_files:
    print ('file',f)
    parts = f.split("/")
    file_date = parts[7] # 7 6  
    if start_date <= file_date <= end_date: 
            infile.append(f)

# Initialize SSH time series
ssh_ts_all = []

# Read data from NetCDF files
grid_info = False
for nc2open in infile:
    print('Processing:', nc2open)
    model = nc.Dataset(nc2open, 'r')
    ssh_ts = np.array(model.variables['sossheig'][:, lat_idx, lon_idx])
    ssh_ts_all = np.concatenate((ssh_ts_all, ssh_ts))

    if not grid_info:
        lats = np.round(model.variables['nav_lat'][lat_idx, lon_idx],2)
        lons = np.round(model.variables['nav_lon'][lat_idx, lon_idx],2)
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
#print ('Check clean time series:',time_series_clean.shape,time_series_clean)

#### SPECTRUM ANALYSIS

spt_len = len(time_series_clean)
print('Time series values:', spt_len)

# Compute FFT
if flag_segmented_spectrum:

    print(f"Segmented spectrum: averaging over {segment_len_days}-day segments")

    segment_len = int((segment_len_days * 86400) / dt)
    num_segments = len(time_series_clean) // segment_len
    print(f"Segment length (in steps): {segment_len}, Total segments: {num_segments}")

    all_amplitudes = []

    for i in range(num_segments):
        segment = time_series_clean[i * segment_len : (i + 1) * segment_len]
        fft_segment = np.fft.fft(segment)
        freq_segment = np.fft.fftfreq(len(segment), d=dt)
        mask = freq_segment > 0
        amp_segment = (np.abs(fft_segment[mask]) ** 2) / len(segment)
        all_amplitudes.append(amp_segment)

    # Average the spectra
    amplitudes = np.mean(all_amplitudes, axis=0)
    # Select only positive frequencies 
    freq_positive = freq_segment[mask]  # Same for all segments
    # Compute Periods in hours
    periods = 1 / freq_positive / 3600

else:
    # Classical spectrum from full series
    fft = np.fft.fft(time_series_clean)
    freq = np.fft.fftfreq(len(time_series_clean), d=dt)

    # Select only positive frequencies (excluding zero and Nyquist)
    half_len = len(time_series_clean) // 2
    freq_positive = freq[:half_len]
    fft_positive = fft[:half_len]
    mask = freq_positive > 0
    freq_positive = freq_positive[mask]
    fft_positive = fft_positive[mask]

    # Compute Power Spectrum Density (normalized)
    amplitudes = (np.abs(fft_positive) ** 2) / len(time_series_clean) #**2
    # Compute Periods in hours
    periods = 1 / freq_positive / 3600

if flag_filter == 'true':
    print('Filter = true')
    high_pass_threshold = 1 / (th_filter * 3600)  # Hz

    # Mask to apply on freq_positive
    freq_mask = freq_positive >= high_pass_threshold

    # Apply the mask to filter frequencies and corresponding amplitudes
    amplitudes = amplitudes[freq_mask]
    freq_positive = freq_positive[freq_mask]
    periods = 1 / freq_positive / 3600

# Smooth
if flag_smooth == 'true':
   print ('Smooth = true')
   amp_smooth = gaussian_filter1d(amplitudes, sigma=sigma)
else:
   amp_smooth = amplitudes
   if flag_smooth == 'plot':
      amp_smooth_2plot = gaussian_filter1d(amplitudes, sigma=sigma)

if energy_threshold_ratio == 0: 

   # Found peaks in spt and in amplitude:
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

else:

   # Find peaks in the smoothed power spectrum (amp_smooth = PSD)
   amp_peaks, _ = find_peaks(amp_smooth)
   amp_peak_frequencies = freq_positive[amp_peaks]
   amp_peak_amplitudes = amp_smooth[amp_peaks]  # these are energies (m²)
   
   # Compute total spectral energy (sum of all PSD values)
   total_energy = np.sum(amp_smooth)
   
   # Define energy threshold: keep only peaks contributing ≥ 10% of total
   energy_threshold = energy_threshold_ratio * total_energy
   
   # Filter peaks by energy contribution
   mask_significant = amp_peak_amplitudes >= energy_threshold
   amp_peak_amplitudes = amp_peak_amplitudes[mask_significant]
   amp_peak_frequencies = amp_peak_frequencies[mask_significant]
   
   # Sort by descending amplitude
   sorted_indices_peak_amp = np.argsort(amp_peak_amplitudes)[::-1]
   amp_peak_amplitudes_sorted = amp_peak_amplitudes[sorted_indices_peak_amp]
   amp_peak_frequencies_sorted = amp_peak_frequencies[sorted_indices_peak_amp]

   period_sorted = 1 / amp_peak_frequencies_sorted / 3600
   amplitude_sorted = amp_peak_amplitudes_sorted

   print(f"Energy threshold: {energy_threshold:.4e}, Significant peaks: {len(amp_peak_amplitudes)}")

final_periods = []
final_amplitudes = []

for i in range(len(period_sorted)):
    if period_sorted[i] > 28: #30:
        tolerance = 10  
    elif 25 <= period_sorted[i] <= 28: #30:
        tolerance = 5 #10  
    elif 12 <= period_sorted[i] < 25:
        tolerance = 2  
    elif 6 <= period_sorted[i] < 12:
        tolerance = 0.5   
    else:
        tolerance = 0.1  

    keep_mode = True  

    for j in range(len(final_periods)):
        if abs(period_sorted[i] - final_periods[j]) <= tolerance:
            if amplitude_sorted[i] <= final_amplitudes[j]:
                keep_mode = False  
                break
            else:
                final_periods[j] = period_sorted[i]
                final_amplitudes[j] = amplitude_sorted[i]
                keep_mode = False  
                break
    
    if keep_mode:
        final_periods.append(period_sorted[i])
        final_amplitudes.append(amplitude_sorted[i])

#print("Final Periods:", final_periods)
#print("Final Amplitudes:", final_amplitudes)
amp_peak_amplitudes_sorted=np.array(final_amplitudes)
amp_peak_period_sorted=np.array(final_periods)

# Count the modes:
if n_modes == 'auto':
   n_modes=len(amp_peak_period_sorted)
   print (n_modes,'modes found')

# Order by Period
if flag_T_order == 1:
   sorted_by_period = np.argsort(amp_peak_period_sorted)[::-1]
   amp_peak_period_sorted = amp_peak_period_sorted[sorted_by_period]
   amp_peak_amplitudes_sorted = amp_peak_amplitudes_sorted[sorted_by_period]

# Time array for SSH plot (convert to hours)
ssh_time = np.arange(0, spt_len * dt, dt) / 3600  

# Select the main modes based on amplitude
n_valid = min(len(amplitudes), len(freq_positive))  # Ensure valid index range
top_indices_amp = np.argpartition(amplitudes[:n_valid], -n_modes)[-n_modes:]  # Select indices of top amplitudes
sorted_indices_amp = np.argsort(amplitudes[top_indices_amp])[::-1]  # Sort by descending amplitude

# Extract the corresponding frequencies, periods, and amplitudes for amplitude-based selection
top_freq_positive_amp = freq_positive[top_indices_amp][sorted_indices_amp]
top_periods_amp = periods[top_indices_amp][sorted_indices_amp]
top_amplitudes_amp = amplitudes[top_indices_amp][sorted_indices_amp]


#######################
# PLOT SSH
# Plot the whole period
plt.figure(figsize=(18, 8))
plt.rc('font', size=20)
plt.title(f'SSH at lat='+str(lats)+' lon='+str(lons)+' '+Med_reg)
plt.plot(ssh_time, time_series_clean, '-',linewidth=2, label=f'SSH at lat='+str(lats)+' lon='+str(lons))
plt.xlabel('Time (h)')
plt.ylabel('SSH (m)')
plt.grid()
plt.legend(loc='upper right')
plt.savefig(work_dir+f'ssh_{lat_idx}_{lon_idx}_{exp}.png')

# PLOT POWER SPECTRUM
plt.figure(figsize=(27, 14))
ax = plt.subplot(111)
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.75, box.height])
plt.title(f'Power Spectrum at lat='+str(lats)+' lon='+str(lons)+' '+Med_reg)

# Mark the main modes based on peak finder
mode_colors = plt.cm.rainbow(np.linspace(0, 1, n_modes)) 
for i in range(0,n_modes): 
    try:
       plt.axvline(amp_peak_period_sorted[i], color=mode_colors[i],linestyle='--',linewidth=4,label=f'Mode {i} (T={amp_peak_period_sorted[i]:.2f} h, E={amp_peak_amplitudes_sorted[i]:.3f} m2*s2)')
       #plt.text(1/amp_peak_frequencies_sorted[i]/3600, plt.ylim()[0] - 0.1, f'{1/amp_peak_frequencies_sorted[i]/3600}', ha='center', va='top')
    except:
       print ('Nan')
f_Nyq=dt*2/3600
plt.loglog(periods[periods>f_Nyq], amplitudes[periods>f_Nyq], marker='o', linestyle='-',linewidth=4,color='navy', label='Modes Energy')
if flag_smooth == 'true':
   plt.loglog(periods[periods>f_Nyq], amp_smooth[periods>f_Nyq], marker='o', linestyle='-',linewidth=4,color='tab:green', label='Smoothed Modes Energy')
elif flag_smooth == 'plot':
   plt.loglog(periods[periods>f_Nyq], amp_smooth_2plot[periods>f_Nyq], marker='o', linestyle='-',linewidth=4,color='tab:green', label='Smoothed Modes Energy')
plt.xlabel('Period (h)')
plt.ylabel('Power Spectrum (m2*s2)')
plt.xlim(th_filter-1,dt*1/3600)
plt.ylim(0.0000001,0.5)
plt.text(24,plt.ylim()[0],'24', ha='center', va='top')
plt.text(12,plt.ylim()[0],'12', ha='center', va='top')
plt.text(6,plt.ylim()[0],'6', ha='center', va='top')
plt.grid()
plt.legend(loc='upper right') 

# Add the table
table_data = []
for i in range(n_modes):
    table_data.append([f'Mode {i}', f'{amp_peak_period_sorted[i]:.2f} h', f'{amp_peak_amplitudes_sorted[i]:.3f} m2*s2'])

col_labels = ['Mode', 'Period (h)', 'Energy (m2*s2)']

table_ax = plt.table(cellText=table_data,
                     colLabels=col_labels,
                     loc='bottom',
                     cellLoc='center',
                     bbox=[0.15, -0.4, 0.7, 0.3])  

#table_ax.auto_set_font_size(False)
table_ax.scale(1.2, 1.2)

for i, label in enumerate(col_labels):
    cell = table_ax.get_celld()[(0, i)]  
    cell.set_facecolor('#d3d3d3')  
    cell.set_text_props(weight='bold')
    
for key, cell in table_ax.get_celld().items():
    cell.set_height(0.2)  
    #cell.set_fontsize(16)

plt.tight_layout()
plt.savefig(work_dir+f'pow_{lat_idx}_{lon_idx}_{exp}.png')

# Amp no log plot
plt.figure(figsize=(27, 14))
ax = plt.subplot(111)
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.75, box.height])
plt.title(f'Power Spectrum at lat='+str(lats)+' lon='+str(lons)+' '+Med_reg)

# Mark the main modes based on peak finder
mode_colors = plt.cm.rainbow(np.linspace(0, 1, n_modes)) 
for i in range(0,n_modes): 
    try:
       plt.axvline(amp_peak_period_sorted[i], color=mode_colors[i],linestyle='--',linewidth=4,label=f'Mode {i} (T={amp_peak_period_sorted[i]:.2f} h, E={amp_peak_amplitudes_sorted[i]:.3f} m2*s2)')
       #plt.text(1/amp_peak_frequencies_sorted[i]/3600, plt.ylim()[0] - 0.1, f'{1/amp_peak_frequencies_sorted[i]/3600}', ha='center', va='top')
    except:
       print ('Nan')
f_Nyq=dt*2/3600
plt.loglog(periods[periods>f_Nyq], amplitudes[periods>f_Nyq], marker='o', linestyle='-',linewidth=4,color='navy', label='Modes Energy')
if flag_smooth == 'true':
   plt.loglog(periods[periods>f_Nyq], amp_smooth[periods>f_Nyq], marker='o', linestyle='-',linewidth=4,color='tab:green', label='Smoothed Modes Energy')
elif flag_smooth == 'plot':
   plt.loglog(periods[periods>f_Nyq], amp_smooth_2plot[periods>f_Nyq], marker='o', linestyle='-',linewidth=4,color='tab:green', label='Smoothed Modes Energy')
plt.xlabel('Period (h)')
plt.ylabel('Power Spectrum (m2*s2)')
plt.xlim(th_filter-1,dt*1/3600)
plt.ylim(0.00000001,0.04)

plt.text(24,plt.ylim()[0],'24', ha='center', va='top')
plt.text(12,plt.ylim()[0],'12', ha='center', va='top')
plt.text(6,plt.ylim()[0],'6', ha='center', va='top')

plt.grid()
plt.yscale('linear')
plt.legend(loc='upper right') 

# Add the table
table_data = []
for i in range(n_modes):
    table_data.append([f'Mode {i}', f'{amp_peak_period_sorted[i]:.2f} h', f'{amp_peak_amplitudes_sorted[i]:.3f} m2*s2'])

col_labels = ['Mode', 'Period (h)', 'Energy (m2*s2)']

table_ax = plt.table(cellText=table_data,
                     colLabels=col_labels,
                     loc='bottom',
                     cellLoc='center',
                     bbox=[0.15, -0.4, 0.7, 0.3])

#table_ax.auto_set_font_size(False)
table_ax.scale(1.2, 1.2)

for i, label in enumerate(col_labels):
    cell = table_ax.get_celld()[(0, i)]
    cell.set_facecolor('#d3d3d3')
    cell.set_text_props(weight='bold')

for key, cell in table_ax.get_celld().items():
    cell.set_height(0.2)
    #cell.set_fontsize(16)

plt.tight_layout()
plt.savefig(work_dir+f'pow_nolog_{lat_idx}_{lon_idx}_{exp}.png')


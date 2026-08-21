import glob
import sys
import re
import os
import xarray as xr
import netCDF4 as nc
import numpy as np
from scipy.signal import periodogram
import heapq
from functools import partial
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter, NullFormatter
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

###################
# Select the period
# Cylc archive structure:
#infile = []
#for f in all_files:
#    print ('file',f)
#    parts = f.split("/")
#    file_date = parts[7] # 7 6
#    if start_date <= file_date <= end_date:
#            infile.append(f)
# General archive structure:
print("Num of in files:", len(all_files))
if len(all_files) == 0:
    print("WARNING: template file not found ", file_template)
date_re = re.compile(r'(\d{8})')
infile = []
no_date_files = []
for f in all_files:
    basename = os.path.basename(f)
    m = date_re.search(basename)
    if m:
        file_date = m.group(1)
        if start_date <= file_date <= end_date:
            infile.append(f)
    else:
        no_date_files.append(f)
# Debug
print("Selected files:", len(infile))
if len(no_date_files) > 0:
    print("This files do not have a date in the name..:")
    for ff in no_date_files[:5]:
        print("  ", ff)
if len(infile) == 0:
    print("WARNING: NO files range", start_date, "-", end_date)

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
        lats = float(np.round(model.variables['nav_lat'][lat_idx, lon_idx], 2))
        lons = float(np.round(model.variables['nav_lon'][lat_idx, lon_idx], 2))
        grid_info = True
        print(f'I am working on {lats:.2f} {lons:.2f}')

    model.close()

print ("I am workign on period",start_date,"-",end_date," Freq:",dt,"s Num of inputs:",len(ssh_ts_all))

# Convert SSH time series to NumPy array and remove NaNs
time_series_point = np.array(ssh_ts_all)
valid_indices = np.logical_not(np.isnan(time_series_point))
time_series_clean = time_series_point[valid_indices]

#### SPECTRUM ANALYSIS
# Compute FFT

spt_len = len(time_series_clean)
print('Time series values:', spt_len)

if flag_segmented_spectrum:

    # ------------------------------------------------------------------
    # Metodo di Welch
    # ------------------------------------------------------------------
    segment_len = int((segment_len_days * 86400) / dt)
    step        = int((segment_step_days * 86400) / dt)

    if len(time_series_clean) < segment_len:
        raise ValueError(
            f"Time series too short: {len(time_series_clean)} steps vs "
            f"segment length {segment_len}. Reduce segment_len_days or set "
            "flag_segmented_spectrum=False."
        )

    starts = list(range(0, len(time_series_clean) - segment_len + 1, step))
    print(f"Segmented spectrum: {segment_len_days}-day windows shifted by "
          f"{segment_step_days} day(s) -> {len(starts)} windows "
          f"(segment length: {segment_len} steps)")

    spt_segments = []
    amp_segments = []

    for i0 in starts:
        #segment = time_series_clean[i0 : i0 + segment_len]
        segment = detrend(time_series_clean[i0 : i0 + segment_len], type='linear')

        if flag_hanning != 0:
            hanning_window = np.hanning(len(segment))
            segment_windowed = segment * hanning_window
            segment_windowed /= hanning_window.mean()  # normalizzazione
        else:
            segment_windowed = segment

        N_real = len(segment_windowed)

        if flag_nfft != 0:
            fft_segment = np.fft.fft(segment_windowed, n=N_fft)
            freq_segment = np.fft.fftfreq(N_fft, d=dt)
        else:
            fft_segment = np.fft.fft(segment_windowed)
            freq_segment = np.fft.fftfreq(len(segment_windowed), d=dt)

        # Select only positive frequencies
        N_used = len(fft_segment)
        half_len = N_used // 2
        freq_positive = freq_segment[:half_len]
        fft_positive = fft_segment[:half_len]
        mask = freq_positive > 0
        freq_positive = freq_positive[mask]
        fft_positive = fft_positive[mask]

        # Power Spectral Density
        spt_seg = (np.abs(fft_positive)**2) / N_real

        amp_seg = (2 / N_real) * np.abs(fft_positive)

        spt_segments.append(spt_seg)
        amp_segments.append(amp_seg)

    # Average the spectra over the windows
    amplitudes = np.mean(amp_segments, axis=0)
    # Select only positive frequencies
    freq_positive = freq_positive  # Same for all segments
    # Compute Periods in hours
    periods = 1 / freq_positive / 3600

else:
    # Classical spectrum from full series
    fft = np.fft.fft(time_series_clean)
    freq = np.fft.fftfreq(spt_len, d=dt)

    # Select only positive frequencies (excluding zero and Nyquist)
    half_spt_len = spt_len // 2
    freq_positive = freq[:half_spt_len]
    fft_positive = fft[:half_spt_len]

    mask = freq_positive > 0
    freq_positive = freq_positive[mask]
    fft_positive = fft_positive[mask]

    amplitudes = (2 / spt_len) * np.abs(fft_positive)
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

if amplitude_threshold_ratio == 0:
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

   # Find peaks in the smoothed amplitude spectrum
   amp_peaks, _ = find_peaks(amp_smooth)
   amp_peak_frequencies = freq_positive[amp_peaks]
   amp_peak_amplitudes = amp_smooth[amp_peaks]

   # Compute total spectral amplitude
   total_amplitude = np.sum(amp_smooth)

   # Define threshold: keep only peaks contributing >= amplitude_threshold_ratio
   amplitude_threshold = amplitude_threshold_ratio * total_amplitude

   # Filter peaks by amplitude contribution
   mask_significant = amp_peak_amplitudes >= amplitude_threshold
   amp_peak_amplitudes = amp_peak_amplitudes[mask_significant]
   amp_peak_frequencies = amp_peak_frequencies[mask_significant]

   # Sort by descending amplitude
   sorted_indices_peak_amp = np.argsort(amp_peak_amplitudes)[::-1]
   amp_peak_amplitudes_sorted = amp_peak_amplitudes[sorted_indices_peak_amp]
   amp_peak_frequencies_sorted = amp_peak_frequencies[sorted_indices_peak_amp]

   period_sorted = 1 / amp_peak_frequencies_sorted / 3600
   amplitude_sorted = amp_peak_amplitudes_sorted

   print(f"Amplitude threshold: {amplitude_threshold:.4e}, Significant peaks: {len(amp_peak_amplitudes)}")


# ----------------------------------------------------------------------
# Calcolo della tolleranza per l'accorpamento dei picchi vicini
#   tolerance = (T^2 * delta_f) * (1 + extra_unc), arrotondata a 1 cifra
#   significativa, con un valore minimo min_unc.
# ----------------------------------------------------------------------
Ttot_h  = segment_len_days * 24
delta_f = 1.0 / Ttot_h        # risoluzione spettrale [cph]

def round_to_1_sigfig(x):
    if x == 0:
        return 0
    return round(x, -int(np.floor(np.log10(abs(x)))))

final_periods = []
final_amplitudes = []

for i in range(len(period_sorted)):

    tolerance = (period_sorted[i] ** 2) * delta_f
    tolerance = tolerance + extra_unc * tolerance
    tolerance = round_to_1_sigfig(tolerance)
    if tolerance < min_unc:
        tolerance = min_unc

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

print("Final periods (h):", np.round(amp_peak_period_sorted, 2))
print("Final amplitudes (m):", np.round(amp_peak_amplitudes_sorted, 4))

# Time array for SSH plot (convert to hours)
ssh_time = np.arange(0, spt_len * dt, dt) / 3600

# Now select the main modes based on amplitude
n_valid = min(len(amplitudes), len(freq_positive))  # Ensure valid index range
top_indices_amp = np.argpartition(amplitudes[:n_valid], -n_modes)[-n_modes:]  # Select indices of top amplitudes
sorted_indices_amp = np.argsort(amplitudes[top_indices_amp])[::-1]  # Sort by descending amplitude

# Extract the corresponding frequencies, periods, and amplitudes
top_freq_positive_amp = freq_positive[top_indices_amp][sorted_indices_amp]
top_periods_amp = periods[top_indices_amp][sorted_indices_amp]
top_amplitudes_amp = amplitudes[top_indices_amp][sorted_indices_amp]


# ======================================================================
# Funzioni di plot
# ======================================================================

f_Nyq = dt * 2 / 3600
mode_colors = plt.cm.rainbow(np.linspace(0, 1, max(n_modes, 1)))

# Tick dell'asse x (asse logaritmico, decrescente da th_filter a f_Nyq)
xticks_all = [36, 24, 18, 12, 9, 6, 4, 3, 2]


def set_xticks_periods(ax):
    """Tick espliciti in ore sull'asse dei periodi (log, decrescente)."""
    xmin, xmax = sorted([th_filter - 1, dt * 1 / 3600])
    ticks = [t for t in xticks_all if xmin <= t <= xmax]
    ax.set_xticks(ticks)
    ax.set_xticklabels([str(t) for t in ticks])
    ax.xaxis.set_minor_formatter(NullFormatter())


def plot_spectrum(ax, log_scale):
    """Disegna lo spettro di ampiezza con i modi identificati."""

    # Linee verticali sui modi identificati
    for i in range(0, n_modes):
        try:
            ax.axvline(amp_peak_period_sorted[i], color=mode_colors[i],
                       linestyle='--', linewidth=4,
                       label=f'T={amp_peak_period_sorted[i]:.2f} h, '
                             f'Amp={amp_peak_amplitudes_sorted[i]:.3f} m')
        except Exception:
            print('Nan')

    sel = periods > f_Nyq
    ax.plot(periods[sel], amplitudes[sel], marker='o', linestyle='-',
            linewidth=4, color='black', label='Amplitude Spectrum')

    if flag_smooth == 'true':
        ax.plot(periods[sel], amp_smooth[sel], marker='o', linestyle='-',
                linewidth=4, color='tab:green', label='Smoothed Amplitude Spectrum')
    elif flag_smooth == 'plot':
        ax.plot(periods[sel], amp_smooth_2plot[sel], marker='o', linestyle='-',
                linewidth=4, color='tab:green', label='Smoothed Amplitude Spectrum')

    ax.set_xscale('log')
    if log_scale:
        ax.set_yscale('log')
        ax.set_ylim(0.0000001, 0.5)
    else:
        ax.set_yscale('linear')
        ax.set_ylim(0.0, 0.08)

    ax.set_xlim(th_filter - 1, dt * 1 / 3600)
    ax.set_xlabel('Period (h)')
    ax.set_ylabel('Mode Amplitude (m)')
    set_xticks_periods(ax)
    ax.grid(True)


def plot_ssh(ax):
    """Disegna la serie temporale di SSH."""
    ax.plot(ssh_time, time_series_clean, '-', linewidth=2,
            label=f'SSH at lat={lats:.2f} lon={lons:.2f}')
    ax.set_xlabel('Time (h)')
    ax.set_ylabel('SSH (m)')
    ax.grid(True)
    ax.axhline(y=0, color='k', linewidth=1.8)
    ax.legend(loc='upper right', fontsize=24)


#######################
# PLOT SSH (invariato)
plt.figure(figsize=(18, 8))
plt.rc('font', size=20)
plt.title(f'SSH at lat={lats:.2f} lon={lons:.2f} {Med_reg}')
plt.plot(ssh_time, time_series_clean, '-', linewidth=2,
         label=f'SSH at lat={lats:.2f} lon={lons:.2f}')
plt.xlabel('Time (h)')
plt.ylabel('SSH (m)')
plt.grid()
plt.legend(loc='upper right')
plt.savefig(work_dir+f'ssh_{lat_idx}_{lon_idx}_{exp}.png')
plt.close()


#######################
# PLOT AMPLITUDE SPECTRUM - LOG
plt.rc('font', size=24)
fig = plt.figure(figsize=(27, 14))
ax = plt.subplot(111)
ax.set_title(f'Modes amplitudes at lat={lats:.2f} lon={lons:.2f} {Med_reg}')
plot_spectrum(ax, log_scale=True)
ax.legend(loc='upper right', fontsize=26)
plt.tight_layout()
plt.savefig(work_dir+f'amp_{lat_idx}_{lon_idx}_{exp}.png')
plt.close()


#######################
# PLOT AMPLITUDE SPECTRUM - NO LOG
plt.rc('font', size=24)
fig = plt.figure(figsize=(27, 14))
ax = plt.subplot(111)
ax.set_title(f'Modes amplitudes at lat={lats:.2f} lon={lons:.2f}')
plot_spectrum(ax, log_scale=False)
ax.legend(loc='upper right', fontsize=26)
plt.tight_layout()
plt.savefig(work_dir+f'amp_nolog_{lat_idx}_{lon_idx}_{exp}.png')
plt.close()


#######################
# PLOT COMBINATO: SSH (sinistra) + spettro lineare (destra)
plt.rc('font', size=24)
fig = plt.figure(figsize=(38, 11))
gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.5], wspace=0.18)

ax_ssh = fig.add_subplot(gs[0, 0])
ax_ssh.set_title(f'Sea Level - {Med_reg} (lat={lats:.2f} lon={lons:.2f})')
plot_ssh(ax_ssh)

ax_spec = fig.add_subplot(gs[0, 1])
ax_spec.set_title(f'Amplitude Spectrum - {Med_reg} (lat={lats:.2f} lon={lons:.2f})')
plot_spectrum(ax_spec, log_scale=False)
ax_spec.legend(loc='upper right', fontsize=24)

plt.savefig(work_dir+f'ssh_amp_{lat_idx}_{lon_idx}_{exp}.png', bbox_inches='tight')
plt.close()

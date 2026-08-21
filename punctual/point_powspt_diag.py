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

spt_len = len(time_series_clean)
print('Time series values:', spt_len)

# Compute FFT
if flag_segmented_spectrum:

    # ------------------------------------------------------------------
    # MODIFICA 1: finestre SOVRAPPOSTE (metodo di Welch), come nello
    # script di ampiezza. Prima i segmenti erano consecutivi e non
    # sovrapposti: con 744 punti e segment_len=480 si otteneva UNA sola
    # finestra (nessuna media) e si scartavano le ultime 264 ore.
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

    all_powers = []

    for i0 in starts:
        segment = time_series_clean[i0 : i0 + segment_len]

        # Hanning
        if flag_hanning != 0:
            hanning_window = np.hanning(len(segment))
            segment_windowed = segment * hanning_window
            segment_windowed /= hanning_window.mean()  # normalizzazione
        else:
            segment_windowed = segment

        # MODIFICA 2: numero di campioni REALI per la normalizzazione.
        # Lo zero-padding (N_fft) interpola lo spettro ma non aggiunge
        # segnale, quindi non deve entrare nel fattore di normalizzazione.
        N_real = len(segment_windowed)

        # FFT
        if flag_nfft != 0:
            fft_segment = np.fft.fft(segment_windowed, n=N_fft)
            freq_segment = np.fft.fftfreq(N_fft, d=dt)
        else:
            fft_segment = np.fft.fft(segment_windowed)
            freq_segment = np.fft.fftfreq(len(segment_windowed), d=dt)

        # Frequenze positive
        N_used = len(fft_segment)
        half_len = N_used // 2
        freq_positive = freq_segment[:half_len]
        fft_positive = fft_segment[:half_len]
        mask = freq_positive > 0
        freq_positive = freq_positive[mask]
        fft_positive = fft_positive[mask]

        # ---------------------------------------------------------------
        # MODIFICA CHIAVE: qui veniva calcolata l'AMPIEZZA
        #     amp_segment = (2 / N_used) * np.abs(fft_positive)
        # cioe' esattamente la stessa quantita' dello script di ampiezza,
        # pur essendo poi etichettata come "Power Spectrum" / "Energy".
        # Con flag_segmented_spectrum=True i due script producevano quindi
        # lo stesso spettro, e il confronto ampiezza-vs-potenza descritto
        # nel capitolo non poteva evidenziare alcuna differenza.
        # Ora si calcola la potenza come quadrato dell'ampiezza spettrale
        # (proporzionale a |FFT|^2, come richiesto), mediata sulle finestre.
        # ---------------------------------------------------------------
        amp_segment = (2 / N_real) * np.abs(fft_positive)
        pow_segment = amp_segment ** 2
        all_powers.append(pow_segment)

    # Media dei segmenti
    amplitudes = np.mean(all_powers, axis=0)

    # Frequenze positive (uguali per tutti i segmenti)
    freq_positive = freq_positive

    # Periodi in ore (solo se dt in secondi)
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

    # Power spectrum, coerente con il ramo segmentato (ampiezza al quadrato)
    amplitudes = ((2 / len(time_series_clean)) * np.abs(fft_positive)) ** 2
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

   # Found peaks in the power spectrum:
   amp_peaks, _ = find_peaks(amp_smooth)
   amp_peak_frequencies = freq_positive[amp_peaks]
   amp_peak_amplitudes = amp_smooth[amp_peaks]

   # Order based on power
   sorted_indices_peak_amp = np.argsort(amp_peak_amplitudes)[::-1]
   amp_peak_amplitudes_sorted = amp_peak_amplitudes[sorted_indices_peak_amp]
   amp_peak_frequencies_sorted = amp_peak_frequencies[sorted_indices_peak_amp]

   # Remove close modes:
   period_sorted = 1/amp_peak_frequencies_sorted/3600  # periodi dei modi
   amplitude_sorted = amp_peak_amplitudes_sorted  # energie dei modi

   print("Ini Periods:", period_sorted)

else:

   # Find peaks in the power spectrum
   amp_peaks, _ = find_peaks(amp_smooth)
   amp_peak_frequencies = freq_positive[amp_peaks]
   amp_peak_amplitudes = amp_smooth[amp_peaks]  # these are energies

   # Compute total spectral energy
   total_energy = np.sum(amp_smooth)

   # Define energy threshold
   energy_threshold = energy_threshold_ratio * total_energy

   # Filter peaks by energy contribution
   mask_significant = amp_peak_amplitudes >= energy_threshold
   amp_peak_amplitudes = amp_peak_amplitudes[mask_significant]
   amp_peak_frequencies = amp_peak_frequencies[mask_significant]

   # Sort by descending energy
   sorted_indices_peak_amp = np.argsort(amp_peak_amplitudes)[::-1]
   amp_peak_amplitudes_sorted = amp_peak_amplitudes[sorted_indices_peak_amp]
   amp_peak_frequencies_sorted = amp_peak_frequencies[sorted_indices_peak_amp]

   period_sorted = 1 / amp_peak_frequencies_sorted / 3600
   amplitude_sorted = amp_peak_amplitudes_sorted

   print(f"Energy threshold: {energy_threshold:.4e}, Significant peaks: {len(amp_peak_amplitudes)}")


# ----------------------------------------------------------------------
# MODIFICA 3: tolleranza per l'accorpamento dei picchi vicini.
# Prima era una scala fissa (10 / 5 / 2 / 0.5 / 0.1 h) non derivata da
# nulla: con tolleranza 0.5 h nella fascia 6-12 h, i modi a 7.21 e 6.71 h
# (distanti esattamente 0.50 h) venivano fusi e il piu' debole scartato.
# Ora si usa la stessa catena dello script areale (mode_period_tab_amp.py):
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
print("Final power (m2):", amp_peak_amplitudes_sorted)

# Time array for SSH plot (convert to hours)
ssh_time = np.arange(0, spt_len * dt, dt) / 3600

# Select the main modes based on power
n_valid = min(len(amplitudes), len(freq_positive))  # Ensure valid index range
top_indices_amp = np.argpartition(amplitudes[:n_valid], -n_modes)[-n_modes:]
sorted_indices_amp = np.argsort(amplitudes[top_indices_amp])[::-1]

top_freq_positive_amp = freq_positive[top_indices_amp][sorted_indices_amp]
top_periods_amp = periods[top_indices_amp][sorted_indices_amp]
top_amplitudes_amp = amplitudes[top_indices_amp][sorted_indices_amp]

# Compute the inertial freq.
Omega = 7.292115e-5  # rad/s
phi = np.deg2rad(lats) # lat in rad
f_c = 2 * Omega * np.sin(phi) # inertial freq.
T_f = (2 * np.pi / f_c) / 3600 # inertial period


# ======================================================================
# MODIFICA 4: funzioni di plot
#  - via la tabella sotto il grafico (duplicava la legenda)
#  - etichette senza "Mode N": solo periodo ed energia
#  - legenda dentro gli assi, in alto a destra, piu' grande
#  - tick dell'asse x espliciti
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


def plot_power(ax, log_scale, show_inertial=True):
    """Disegna lo spettro di potenza con i modi identificati."""

    if show_inertial:
        ax.axvline(T_f, color='black', linestyle='--', linewidth=3,
                   label=f'Inertial period = {T_f:.1f} h')

    # Linee verticali sui modi identificati
    for i in range(0, n_modes):
        try:
            ax.axvline(amp_peak_period_sorted[i], color=mode_colors[i],
                       linestyle='--', linewidth=4,
                       label=f'T={amp_peak_period_sorted[i]:.2f} h, '
                             f'E={amp_peak_amplitudes_sorted[i]:.2e} m$^2$')
        except Exception:
            print('Nan')

    sel = periods > f_Nyq
    ax.plot(periods[sel], amplitudes[sel], marker='o', linestyle='-',
            linewidth=4, color='navy', label='Power Spectrum')

    if flag_smooth == 'true':
        ax.plot(periods[sel], amp_smooth[sel], marker='o', linestyle='-',
                linewidth=4, color='tab:green', label='Smoothed Power Spectrum')
    elif flag_smooth == 'plot':
        ax.plot(periods[sel], amp_smooth_2plot[sel], marker='o', linestyle='-',
                linewidth=4, color='tab:green', label='Smoothed Power Spectrum')

    ax.set_xscale('log')
    if log_scale:
        ax.set_yscale('log')
        ax.set_ylim(1e-12, 1.0)
    else:
        ax.set_yscale('linear')
        # NOTA: con la potenza definita come ampiezza^2 i valori cambiano di
        # ordine di grandezza rispetto a prima. Scala automatica al primo
        # lancio: se preferisci un limite fisso, sostituisci questa riga.
        ymax = np.nanmax(amplitudes[sel]) if np.any(sel) else 1.0
        ax.set_ylim(0.0, 1.15 * ymax)

    ax.set_xlim(th_filter - 1, dt * 1 / 3600)
    ax.set_xlabel('Period (h)')
    ax.set_ylabel('Power Spectrum (m$^2$)')
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
# PLOT SSH
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
# PLOT POWER SPECTRUM - LOG
plt.rc('font', size=24)
fig = plt.figure(figsize=(27, 14))
ax = plt.subplot(111)
ax.set_title(f'Power Spectrum at lat={lats:.2f} lon={lons:.2f} {Med_reg}')
plot_power(ax, log_scale=True, show_inertial=True)
ax.legend(loc='upper right', fontsize=26)
plt.tight_layout()
plt.savefig(work_dir+f'pow_{lat_idx}_{lon_idx}_{exp}.png')
plt.close()


#######################
# PLOT POWER SPECTRUM - NO LOG
plt.rc('font', size=24)
fig = plt.figure(figsize=(27, 14))
ax = plt.subplot(111)
ax.set_title(f'Power Spectrum at lat={lats:.2f} lon={lons:.2f} {Med_reg}')
plot_power(ax, log_scale=False, show_inertial=False)
ax.legend(loc='upper right', fontsize=26)
plt.tight_layout()
plt.savefig(work_dir+f'pow_nolog_{lat_idx}_{lon_idx}_{exp}.png')
plt.close()


#######################
# PLOT COMBINATO: SSH (sinistra) + spettro di potenza lineare (destra)
plt.rc('font', size=24)
fig = plt.figure(figsize=(38, 11))
gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.5], wspace=0.18)

ax_ssh = fig.add_subplot(gs[0, 0])
ax_ssh.set_title(f'Sea Level - {Med_reg} (lat={lats:.2f} lon={lons:.2f})')
plot_ssh(ax_ssh)

ax_pow = fig.add_subplot(gs[0, 1])
ax_pow.set_title(f'Power Spectrum - {Med_reg} (lat={lats:.2f} lon={lons:.2f})')
plot_power(ax_pow, log_scale=False, show_inertial=False)
ax_pow.legend(loc='upper right', fontsize=22)

plt.savefig(work_dir+f'ssh_pow_{lat_idx}_{lon_idx}_{exp}.png', bbox_inches='tight')
plt.close()


# -------------------------------------------
# Extract the main modes at each grid point
# -------------------------------------------

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
from area_ini import *
mpl.use('Agg')


def round_to_1_sigfig(x):
   if x == 0:
      return 0
   return round(x, -int(np.floor(np.log10(abs(x)))))


def amp_main_modes(lat_idx,lon_idx,ssh_ts_all,dt):

   global n_modes, flag_filter, th_filter
   global energy_threshold_ratio, amplitude_threshold_ratio
   global flag_segmented_spectrum
   global segment_len_days, segment_step_days, flag_T_order
   global extra_unc, min_unc

   # Convert SSH time series to NumPy array and remove NaNs
   time_series_point = np.array(ssh_ts_all)
   valid_indices = np.logical_not(np.isnan(time_series_point))
   time_series_clean = time_series_point[valid_indices]

   if len(time_series_clean) == 0:
        return np.full(8, np.nan), np.full(8, np.nan)

   #### SPECTRUM ANALYSIS

   spt_len = len(time_series_clean)
   #print('Time series values:', spt_len)

   if flag_segmented_spectrum:

      # ------------------------------------------------------------------
      # MODIFICA 1 (= Mod. 1 di point_ampspt_diag.py): finestre SOVRAPPOSTE
      # (metodo di Welch). Prima i segmenti erano consecutivi e non
      # sovrapposti:
      #     num_segments = len(serie) // segment_len
      # che con 720 punti e segment_len=480 dava UNA sola finestra (quindi
      # nessuna media) e scartava le ultime 240 ore. Ora le finestre sono
      # traslate di segment_step_days, come descritto nel Cap. 3.
      # NB: la risoluzione spettrale resta 1/segment_len_days, quindi le
      # incertezze sui periodi NON cambiano; migliora il rapporto S/N.
      # ------------------------------------------------------------------
      segment_len = int((segment_len_days * 86400) / dt)
      step        = int((segment_step_days * 86400) / dt)

      if len(time_series_clean) < segment_len:
         raise ValueError(
            f"Time series too short: {len(time_series_clean)} steps vs "
            f"segment length {segment_len}. Reduce segment_len_days or set "
            "flag_segmented_spectrum=False.")

      starts = list(range(0, len(time_series_clean) - segment_len + 1, step))
      #print(f"Segmented spectrum: {segment_len_days}-day windows shifted by "
      #      f"{segment_step_days} day(s) -> {len(starts)} windows")

      spt_segments = []
      amp_segments = []

      for i0 in starts:

         # MODIFICA 4 (= detrend di point_ampspt_diag.py): rimozione del
         # trend lineare su ogni finestra, prima della finestratura.
         segment = detrend(time_series_clean[i0 : i0 + segment_len], type='linear')

         # Apply Hanning window if requested
         if flag_hanning != 0:
                window = np.hanning(len(segment))
                segment_windowed = segment * window
                segment_windowed /= window.mean()  # normalization
         else:
                segment_windowed = segment

         # MODIFICA 2 (= Mod. 2 di point_ampspt_diag.py): numero di campioni
         # REALI, usato per la normalizzazione. Lo zero-padding (N_fft)
         # interpola lo spettro ma non aggiunge segnale, quindi non deve
         # entrare nel fattore 2/N. Prima era N_used = N_fft = 512 contro 480
         # campioni reali -> ampiezze sottostimate del 6%.
         N_real = len(segment_windowed)

         # FFT with or without zero-padding
         if flag_nfft != 0:
                N_used = N_fft
                fft_segment = np.fft.fft(segment_windowed, n=N_fft)
                freq_segment = np.fft.fftfreq(N_fft, d=dt)
         else:
                N_used = len(segment_windowed)
                fft_segment = np.fft.fft(segment_windowed)
                freq_segment = np.fft.fftfreq(N_used, d=dt)

         # Half spectrum
         half_len = N_used // 2
         freq_half = freq_segment[:half_len]
         fft_half  = fft_segment[:half_len]

         # Select positive frequencies
         mask = freq_half > 0
         freq_positive = freq_half[mask]
         fft_positive = fft_half[mask]

         # Compute spectra (normalized on the REAL number of samples)
         spt_segment = (np.abs(fft_positive) ** 2) / N_real
         amp_segment = (2 / N_real) * np.abs(fft_positive)

         spt_segments.append(spt_segment)
         amp_segments.append(amp_segment)

      # Average the spectra over segments
      spt = np.mean(spt_segments, axis=0)
      amplitudes = np.mean(amp_segments, axis=0)

      # freq_positive is the same for every window
      # Compute periods in hours
      periods = 1 / freq_positive / 3600

   else:

      # Classical spectrum from full series
      spt_len = len(time_series_clean)
      fft = np.fft.fft(time_series_clean)
      freq = np.fft.fftfreq(spt_len, d=dt)

      # Select only positive frequencies (excluding zero and Nyquist)
      half_spt_len = spt_len // 2
      freq_positive = freq[:half_spt_len]
      fft_positive = fft[:half_spt_len]
      mask = freq_positive > 0
      freq_positive = freq_positive[mask]
      fft_positive = fft_positive[mask]

      spt = (np.abs(fft_positive) ** 2) / spt_len
      amplitudes = (2 / spt_len) * np.abs(fft_positive)
      periods = 1 / freq_positive / 3600

   # ---------------------------------------------------------------------
   # MODIFICA 3: filtro passa-alto applicato agli spettri GIA' MEDIATI.
   # Prima il blocco stava fuori dal loop sui segmenti e ricalcolava
   # spt/amplitudes da `fft_positive`, che al termine del loop contiene solo
   # l'ULTIMO segmento: la media sulle finestre veniva quindi buttata via.
   # Inoltre normalizzava su spt_len (serie intera) invece che su N_real.
   # Il bug era invisibile finche' num_segments valeva 1 (vedi Mod. 1).
   # Ora si maschera, come gia' fa f_point_powspt.py e come fa il puntuale.
   # ---------------------------------------------------------------------
   if flag_filter == 'true':
      high_pass_threshold = 1 / (th_filter * 3600)  # Hz
      freq_mask = freq_positive >= high_pass_threshold

      amplitudes    = amplitudes[freq_mask]
      spt           = spt[freq_mask]
      freq_positive = freq_positive[freq_mask]
      periods       = 1 / freq_positive / 3600

   # NB: nessuno smoothing, coerentemente con flag_smooth='false' nel puntuale
   spt_smooth = spt
   amp_smooth = amplitudes

   if amplitude_threshold_ratio == 0:
      # Found peaks in amplitude spectrum:
      amp_peaks, _ = find_peaks(amp_smooth)
      amp_peak_frequencies = freq_positive[amp_peaks]
      amp_peak_amplitudes = amp_smooth[amp_peaks]

      # Order based on amplitudes
      sorted_indices_peak_amp = np.argsort(amp_peak_amplitudes)[::-1]
      amp_peak_amplitudes_sorted = amp_peak_amplitudes[sorted_indices_peak_amp]
      amp_peak_frequencies_sorted = amp_peak_frequencies[sorted_indices_peak_amp]

      period_sorted = 1/amp_peak_frequencies_sorted/3600
      amplitude_sorted = amp_peak_amplitudes_sorted

   else:

      # Find peaks in the amplitude spectrum
      amp_peaks, _ = find_peaks(amp_smooth)
      amp_peak_frequencies = freq_positive[amp_peaks]
      amp_peak_amplitudes = amp_smooth[amp_peaks]

      # Compute total spectral amplitude
      total_amplitude = np.sum(amp_smooth)

      # Keep only peaks contributing >= amplitude_threshold_ratio
      amplitude_threshold = amplitude_threshold_ratio * total_amplitude

      mask_significant = amp_peak_amplitudes >= amplitude_threshold
      amp_peak_amplitudes = amp_peak_amplitudes[mask_significant]
      amp_peak_frequencies = amp_peak_frequencies[mask_significant]

      # Sort by descending amplitude
      sorted_indices_peak_amp = np.argsort(amp_peak_amplitudes)[::-1]
      amp_peak_amplitudes_sorted = amp_peak_amplitudes[sorted_indices_peak_amp]
      amp_peak_frequencies_sorted = amp_peak_frequencies[sorted_indices_peak_amp]

      period_sorted = 1 / amp_peak_frequencies_sorted / 3600
      amplitude_sorted = amp_peak_amplitudes_sorted

   # ---------------------------------------------------------------------
   # MODIFICA 5 (= Mod. 3 di point_ampspt_diag.py): accorpamento dei picchi
   # vicini. Qui il blocco era del tutto ASSENTE (c'era solo il commento
   # "# Remove close modes:"), quindi ogni massimo locale dello spettro
   # finiva nel netCDF e la de-duplicazione avveniva solo a valle, nel
   # raggruppamento globale di mode_period_tab_amp.py. Conseguenza: due bin
   # adiacenti entrambi massimi locali potevano finire in due gruppi diversi.
   # Ora si usa la stessa catena del puntuale:
   #   tolerance = (T^2 * delta_f) * (1 + extra_unc), arrotondata a 1 cifra
   #   significativa, con valore minimo min_unc.
   # I picchi sono ordinati per ampiezza decrescente, quindi in caso di
   # conflitto sopravvive il piu' forte.
   # ---------------------------------------------------------------------
   if flag_segmented_spectrum:
      Ttot_h = segment_len_days * 24
   else:
      Ttot_h = spt_len * dt / 3600.
   delta_f = 1.0 / Ttot_h        # risoluzione spettrale [cph]

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

   # Final arrays
   amp_peak_amplitudes_sorted = np.array(final_amplitudes)
   amp_peak_period_sorted     = np.array(final_periods)

   # Count the modes:
   # MODIFICA 6: n_modes_all era assegnato solo dentro il ramo 'auto' ma
   # usato sempre -> NameError latente se n_modes != 'auto'.
   if n_modes == 'auto':
      n_modes_all = len(amp_peak_period_sorted)
   else:
      n_modes_all = min(int(n_modes), len(amp_peak_period_sorted))

   # Order by Period
   if flag_T_order == 1:
      sorted_by_period = np.argsort(amp_peak_period_sorted)[::-1]
      amp_peak_period_sorted = amp_peak_period_sorted[sorted_by_period]
      amp_peak_amplitudes_sorted = amp_peak_amplitudes_sorted[sorted_by_period]

   #############
   # Return the values

   amp_peak_periods_main = []
   amp_peak_amplitudes_main = []

   for i in range(0, n_modes_all):
     try:
       amp_peak_periods_main.append(amp_peak_period_sorted[i])
       amp_peak_amplitudes_main.append(amp_peak_amplitudes_sorted[i])
     except:
       amp_peak_periods_main.append(np.nan)
       amp_peak_amplitudes_main.append(np.nan)

   return np.array(amp_peak_periods_main), np.array(amp_peak_amplitudes_main)

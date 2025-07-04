
# -------------------------------------------
# Extract the main 8 modes at each grid point
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
mpl.use('Agg')

def pow_main_modes(lat_idx,lon_idx,ssh_ts_all,dt):

   ########
   
   # Flag and threshold [h] for filtering the spectrum the threshold is also used as plot minimum 
   flag_filter='true'
   th_filter=40

   # Filter out modes with low energy, e.g. 0.10 means that all the modes with energy<10% of the total energy are rm (to avoid the filtering set energy_threshold_ratio = 0)
   energy_threshold_ratio = 0.002

   # Number of modes to analyze (n_modes = 'auto' if you want to analyze all the modes)
   n_modes = 'auto'

   # Flag: use segmented (averaged) spectrum or full time series
   flag_segmented_spectrum = True  # False to disable and use full spectrum
   segment_len_days = 5  # length of each segment (in days) if segmented spectrum is used

   # To order by period instead of by amplitude set flag_T_order = 1
   flag_T_order = 1

   ###################
   
   # Convert SSH time series to NumPy array and remove NaNs
   time_series_point = np.array(ssh_ts_all)
   valid_indices = np.logical_not(np.isnan(time_series_point))
   time_series_clean = time_series_point[valid_indices]
   
   if len(time_series_clean) == 0:
        return np.full(8, np.nan), np.full(8, np.nan)

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

       # Compute Power Spectrum Density (normalized)
       amplitudes = (np.abs(fft_positive) ** 2) / spt_len
       # Compute Periods in hours
       periods = 1 / freq_positive / 3600


   if flag_filter == 'true':
      high_pass_threshold = 1 / (th_filter * 3600)  # threshold Hz

      if flag_segmented_spectrum:
        freq_mask = freq_positive >= high_pass_threshold
        amplitudes = amplitudes[freq_mask]
        freq_positive = freq_positive[freq_mask]
        periods = 1 / freq_positive / 3600
      else:
        # Apply high-pass filter: Set frequencies below the threshold to zero
        fft_positive[freq_positive < high_pass_threshold] = 0
        # Recompute the power spectrum and amplitudes after filtering
        amplitudes = (np.abs(fft_positive) ** 2) / spt_len

   amp_smooth = amplitudes

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

   final_periods = period_sorted
   final_amplitudes = amplitude_sorted

#   final_periods = []
#   final_amplitudes = []
#
#   for i in range(len(period_sorted)):
#       if period_sorted[i] > 28: #30:
#           tolerance = 10  
#       elif 25 <= period_sorted[i] <= 28: #30:
#           tolerance = 5 #10  
#       elif 12 <= period_sorted[i] < 25:
#           tolerance = 2  
#       elif 6 <= period_sorted[i] < 12:
#           tolerance = 0.5   
#       else:
#           tolerance = 0.1  
# 
#   keep_mode = True  
#   for j in range(len(final_periods)):
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
#   if keep_mode:
#        final_periods.append(period_sorted[i])
#        final_amplitudes.append(amplitude_sorted[i])

   # Final arrays
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

   #############
   # Return the values

   amp_peak_periods_main=[]
   amp_peak_amplitudes_main=[]

   for i in range(0,n_modes):
     try:
       amp_peak_periods_main.append(amp_peak_period_sorted[i])
       amp_peak_amplitudes_main.append(amp_peak_amplitudes_sorted[i])
     except:
       amp_peak_periods_main.append(np.nan)
       amp_peak_amplitudes_main.append(np.nan)
   
   return np.array(amp_peak_periods_main),np.array(amp_peak_amplitudes_main)


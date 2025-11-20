
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
from area_ini import *
mpl.use('Agg')

def amp_main_modes(lat_idx,lon_idx,ssh_ts_all,dt):

   global n_modes, flag_filter, th_filter
   global energy_threshold_ratio, flag_segmented_spectrum
   global segment_len_days, flag_T_order

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

      #print(f"Segmented spectrum: averaging over {segment_len_days}-day segments")

      segment_len = int((segment_len_days * 86400) / dt)
      num_segments = len(time_series_clean) // segment_len
      #print(f"Segment length (in steps): {segment_len}, Total segments: {num_segments}")

      spt_segments = []
      amp_segments = []

      for i in range(num_segments):
         segment = time_series_clean[i * segment_len : (i + 1) * segment_len]

         # Apply Hanning window if requested
         if flag_hanning != 0:
                window = np.hanning(len(segment))
                segment_windowed = segment * window
                segment_windowed /= window.mean()  # normalization
         else:
                segment_windowed = segment

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

         # Compute spectra
         spt_segment = (np.abs(fft_positive) ** 2) / N_used
         amp_segment = (2 / N_used) * np.abs(fft_positive)

         spt_segments.append(spt_segment)
         amp_segments.append(amp_segment)

      # Average the spectra over segments
      spt = np.mean(spt_segments, axis=0)
      amplitudes = np.mean(amp_segments, axis=0)

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
      periods = 1 / freq_positive / 3600
      amplitudes = (2 / spt_len) * np.abs(fft_positive)

   if flag_filter=='true':
      #print ('Filter = true')
      # Apply high-pass filter: Set frequencies below the threshold to zero
      high_pass_threshold = 1 / (th_filter * 3600)  # Corresponding to 48 hours in Hz
      fft_positive[freq_positive < high_pass_threshold] = 0  # Filter out frequencies below the threshold

      # Recompute the power spectrum and amplitudes after filtering
      spt = (np.abs(fft_positive) ** 2) / spt_len #**2
      amplitudes = (2 / spt_len) * np.abs(fft_positive)

   spt_smooth = spt
   amp_smooth = amplitudes

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
   
      #print("Ini Periods:", period_sorted)
   
   else:
   
      # Find peaks in the smoothed power spectrum (amp_smooth = PSD)
      amp_peaks, _ = find_peaks(amp_smooth)
      amp_peak_frequencies = freq_positive[amp_peaks]
      amp_peak_amplitudes = amp_smooth[amp_peaks]
   
      # Compute total spectral amplitude (sum of all PSD values)
      total_amplitude = np.sum(amp_smooth)
   
      # Define energy threshold: keep only peaks contributing ≥ 10% of total
      amplitude_threshold = amplitude_threshold_ratio * total_amplitude
   
      # Filter peaks by energy contribution
      mask_significant = amp_peak_amplitudes >= amplitude_threshold
      amp_peak_amplitudes = amp_peak_amplitudes[mask_significant]
      amp_peak_frequencies = amp_peak_frequencies[mask_significant]
   
      # Sort by descending amplitude
      sorted_indices_peak_amp = np.argsort(amp_peak_amplitudes)[::-1]
      amp_peak_amplitudes_sorted = amp_peak_amplitudes[sorted_indices_peak_amp]
      amp_peak_frequencies_sorted = amp_peak_frequencies[sorted_indices_peak_amp]
   
      period_sorted = 1 / amp_peak_frequencies_sorted / 3600
      amplitude_sorted = amp_peak_amplitudes_sorted
   
      #print(f"Amplitude threshold: {amplitude_threshold:.4e}, Significant peaks: {len(amp_peak_amplitudes)}")

   final_periods = period_sorted
   final_amplitudes = amplitude_sorted

   # Final arrays
   amp_peak_amplitudes_sorted=np.array(final_amplitudes) 
   amp_peak_period_sorted=np.array(final_periods)

   # Count the modes:
   if n_modes == 'auto':
      n_modes_all=len(amp_peak_period_sorted)
      #print (n_modes_all,'modes found')

   # Order by Period
   if flag_T_order == 1:
      sorted_by_period = np.argsort(amp_peak_period_sorted)[::-1]
      amp_peak_period_sorted = amp_peak_period_sorted[sorted_by_period]
      amp_peak_amplitudes_sorted = amp_peak_amplitudes_sorted[sorted_by_period]

   #############
   # Return the values

   amp_peak_periods_main=[]
   amp_peak_amplitudes_main=[]

   for i in range(0,n_modes_all):
     try:
       amp_peak_periods_main.append(amp_peak_period_sorted[i])
       amp_peak_amplitudes_main.append(amp_peak_amplitudes_sorted[i])
     except:
       amp_peak_periods_main.append(np.nan)
       amp_peak_amplitudes_main.append(np.nan)
   
   return np.array(amp_peak_periods_main),np.array(amp_peak_amplitudes_main)


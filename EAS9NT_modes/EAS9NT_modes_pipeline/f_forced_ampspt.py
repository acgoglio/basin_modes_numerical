"""
Funzione core per l'analisi spettrale su un singolo punto griglia,
per la run FORZATA (§3.3.2).

Differenze rispetto a f_point_ampspt.py (pipeline free-oscillation):

1. BUG FIX: il filtro passa-alto viene applicato come MASCHERA
   sull'array gia' mediato sui segmenti (spt/amplitudes), esattamente
   come fa correttamente f_point_powspt.py - non ricalcolato da
   fft_positive (che dopo il loop resta agganciato all'ULTIMO
   segmento soltanto) con normalizzazione su spt_len invece che sulla
   lunghezza del segmento.

2. NIENTE discovery cieca / grouping libero: dopo il peak-finding,
   i picchi vengono usati per DUE scopi distinti:
     a) diagnostico: tutti i picchi trovati, per l'istogramma di
        copertura di validazione (stile Fig. 3.5)
     b) quantitativo: per ciascuno dei 10 modi di riferimento
        (forced_windows_ini.REFERENCE_MODES), si cerca il picco piu'
        vicino dentro la banda T_i +/- tol_i; se non c'e' nessun picco
        in banda, fallback sul valore dello spettro al bin piu'
        vicino a T_i (nessun punto scartato silenziosamente).

Non gestisce piu' n_modes='auto'/flag_T_order (non ha senso quando i
modi sono fissati a priori) - se servono ancora altrove vanno tenuti
in f_point_ampspt.py invariato.
"""

import numpy as np
from scipy.signal import detrend, find_peaks


def compute_point_spectrum(ssh_ts_all, dt, *, flag_hanning, flag_nfft, N_fft,
                            flag_segmented_spectrum, segment_len_days,
                            flag_filter, th_filter):
    """
    Calcola lo spettro di ampiezza per un singolo punto, con la
    normalizzazione corretta indipendentemente dalla combinazione di
    flag (bug fix rispetto a f_point_ampspt.py).

    Ritorna: freq_positive [Hz], periods [h], amplitudes (ampiezza,
    non densita' di potenza), oppure (None, None, None) se la serie e'
    tutta NaN o troppo corta per il segment_len richiesto.
    """
    time_series_point = np.array(ssh_ts_all)
    valid = np.logical_not(np.isnan(time_series_point))
    time_series_clean = time_series_point[valid]

    if len(time_series_clean) == 0:
        return None, None, None

    if flag_segmented_spectrum:
        segment_len = int((segment_len_days * 86400) / dt)
        num_segments = len(time_series_clean) // segment_len

        if num_segments == 0:
            # Finestra piu' corta del segmento richiesto (es. evento di
            # 10gg con segment_len_days=20): niente fallback silenzioso,
            # segnalato esplicitamente al chiamante.
            return None, None, None

        spt_segments = []
        amp_segments = []
        freq_positive = None

        for i in range(num_segments):
            segment = time_series_clean[i * segment_len:(i + 1) * segment_len]
            segment = detrend(segment)  # detrend per finestra (era mancante)

            if flag_hanning != 0:
                window = np.hanning(len(segment))
                segment_w = segment * window
                segment_w /= window.mean()
            else:
                segment_w = segment

            if flag_nfft != 0:
                N_used = N_fft
                fft_segment = np.fft.fft(segment_w, n=N_fft)
                freq_segment = np.fft.fftfreq(N_fft, d=dt)
            else:
                N_used = len(segment_w)
                fft_segment = np.fft.fft(segment_w)
                freq_segment = np.fft.fftfreq(N_used, d=dt)

            half_len = N_used // 2
            freq_half = freq_segment[:half_len]
            fft_half = fft_segment[:half_len]

            mask = freq_half > 0
            freq_pos_seg = freq_half[mask]
            fft_pos_seg = fft_half[mask]

            spt_segments.append((np.abs(fft_pos_seg) ** 2) / N_used)
            amp_segments.append((2 / N_used) * np.abs(fft_pos_seg))
            freq_positive = freq_pos_seg  # identica per ogni segmento

        amplitudes = np.mean(amp_segments, axis=0)
        spt = np.mean(spt_segments, axis=0)
        norm_len = None  # non serve piu': gia' mediato correttamente

    else:
        spt_len = len(time_series_clean)
        clean = detrend(time_series_clean)
        fft = np.fft.fft(clean)
        freq = np.fft.fftfreq(spt_len, d=dt)

        half = spt_len // 2
        freq_positive = freq[:half]
        fft_positive = fft[:half]
        mask = freq_positive > 0
        freq_positive = freq_positive[mask]
        fft_positive = fft_positive[mask]

        spt = (np.abs(fft_positive) ** 2) / spt_len
        amplitudes = (2 / spt_len) * np.abs(fft_positive)

    # ---- BUG FIX: filtro come maschera sull'array gia' mediato ----
    if flag_filter == "true":
        high_pass_threshold = 1 / (th_filter * 3600)
        keep = freq_positive >= high_pass_threshold
        # Maschero (non azzero-e-basta) per restare coerenti con la
        # sostituzione del vecchio comportamento; i bin sotto soglia
        # non parteciperanno al peak-finding perche' find_peaks lavora
        # sull'array intero - li mettiamo a 0 in amplitudes/spt cosi'
        # non vengono trovati come picchi ne' contano nel 99th perc.
        amplitudes = np.where(keep, amplitudes, 0.0)
        spt = np.where(keep, spt, 0.0)

    periods = 1 / freq_positive / 3600  # ore

    return freq_positive, periods, amplitudes


def find_all_peaks(periods, amplitudes):
    """
    Peak-finding sull'intero spettro (per uso sia diagnostico che
    quantitativo). Ritorna array di (periodo, ampiezza) dei picchi.
    """
    if periods is None:
        return np.array([]), np.array([])
    peak_idx, _ = find_peaks(amplitudes)
    return periods[peak_idx], amplitudes[peak_idx]


def extract_reference_mode_amplitudes(periods, amplitudes, peak_periods,
                                       peak_amplitudes, reference_modes,
                                       flag_use_fallback=False):
    """
    Per ciascun modo di riferimento (label, T_ref, tol), ritorna
    l'ampiezza nel punto corrente:

      - se esiste un picco (da find_all_peaks) dentro T_ref +/- tol,
        si prende il picco piu' vicino a T_ref tra quelli in banda
        (massimo locale reale, non il centro banda fisso)
      - altrimenti:
          * flag_use_fallback=True  -> fallback: valore dello spettro
            al bin piu' vicino a T_ref (nessun NaN, ma valore non
            necessariamente un vero massimo locale)
          * flag_use_fallback=False (default) -> NaN: il punto non
            contribuisce a quel modo/finestra (nessun dato "debole"
            mescolato con i veri picchi)

    Ritorna: array shape (n_modes,) di ampiezze, e array bool
    shape (n_modes,) 'used_fallback' per diagnosticare quanto spesso
    si sarebbe dovuto ricorrere al fallback (utile per QA anche
    quando flag_use_fallback=False, per capire quanti punti sono
    diventati NaN e perche').
    """
    n_modes = len(reference_modes)
    out_amp = np.full(n_modes, np.nan)
    used_fallback = np.zeros(n_modes, dtype=bool)

    if periods is None:
        return out_amp, used_fallback

    for i, (label, T_ref, tol) in enumerate(reference_modes):
        if len(peak_periods) > 0:
            in_band = np.abs(peak_periods - T_ref) <= tol
        else:
            in_band = np.array([], dtype=bool)

        if np.any(in_band):
            cand_periods = peak_periods[in_band]
            cand_amps = peak_amplitudes[in_band]
            # picco piu' vicino a T_ref tra quelli in banda
            closest = np.argmin(np.abs(cand_periods - T_ref))
            out_amp[i] = cand_amps[closest]
        else:
            used_fallback[i] = True
            if flag_use_fallback:
                # fallback: bin piu' vicino a T_ref nello spettro completo
                closest_bin = np.argmin(np.abs(periods - T_ref))
                out_amp[i] = amplitudes[closest_bin]
            else:
                # NaN: punto escluso per questo modo/finestra
                out_amp[i] = np.nan

    return out_amp, used_fallback

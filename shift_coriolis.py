import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.use("Agg")  # For non-interactive backend

# Diff con e senza Coriolis:

# Periodi in secondi
T_hours = np.linspace(2, 40, 1000)
T_sec = T_hours * 3600
omega = 2 * np.pi / T_sec  # frequenza in rad/s

# Valori di f per il Mediterraneo (circa)
f_low = 7e-5
f_high = 1e-4

# Calcolo della variazione percentuale: (omega_f - omega)/omega * 100
# Dove omega_f = sqrt(omega^2 + f^2)
delta_low = (np.sqrt(omega**2 + f_low**2) - omega) / omega * 100
delta_high = (np.sqrt(omega**2 + f_high**2) - omega) / omega * 100

# Plot
plt.figure(figsize=(8, 5))
plt.plot(T_hours, delta_low, color='tab:orange', label=r'South Med')
plt.plot(T_hours, delta_high, color='red', label=r'North Med')
plt.axhline(50, color='gray', linestyle='--', label='50% threshold')
plt.xlim(0, 40)
plt.ylim(0, np.max(delta_high[T_hours <= 40])*1.1)
plt.xlabel('Period (hours)')
plt.ylabel('Relative frequency shift [%]')
plt.title('Effect of Coriolis on Mode Frequency')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig('coriolis_shift_perc.png')
plt.close()

# Calcolo della differenza in periodo
# omega = 2*pi / T (in rad/s)
# f_low, f_high dati

T = T_hours * 3600  # Converti T da ore a secondi
omega = 2 * np.pi / T

T_f_low = 2 * np.pi / np.sqrt(omega**2 + f_low**2)
T_f_high = 2 * np.pi / np.sqrt(omega**2 + f_high**2)

delta_T_low = T_f_low - T  # shift assoluto in secondi
delta_T_high = T_f_high - T

plt.figure(figsize=(8, 5))
plt.plot(T_hours, delta_T_low / 3600, color='tab:orange', label='South Med')
plt.plot(T_hours, delta_T_high / 3600, color='red', label='North Med')
plt.xlim(0, 40)
plt.xlabel('Modes Period (hours)')
plt.ylabel('Absolute period shift [hours]')
plt.title('Effect of Coriolis on Mode Period')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig('coriolis_shift_abs.png')
plt.close()

T_hours_bk = T_hours

# Incertezza dovuta alla risoluzione spaziale

# Parametri
dt = 3600  # 1h in secondi
N_1 = 240    # 10 giorni di dati orari
N_2 = 480    # 20 giorni di dati orari
N_3 = 120    # 5 giorni di dati orari
df_1 = 1 / (N_1 * dt)  # risoluzione spettrale (Hz)
df_2 = 1 / (N_2 * dt)  # risoluzione spettrale (Hz)
df_3 = 1 / (N_3 * dt)  # risoluzione spettrale (Hz)
extra_unc = 0.1

# Periodi in ore da 2 a 40
T_hours = np.linspace(2, 40, 500)
T_seconds = T_hours * 3600

# Calcolo incertezza sul periodo
delta_T_seconds_1 = df_1 * T_seconds**2
delta_T_hours_1 = delta_T_seconds_1 / 3600  # in ore
delta_T_hours_1 = delta_T_hours_1 + extra_unc*delta_T_hours_1

delta_T_seconds_2 = df_2 * T_seconds**2
delta_T_hours_2 = delta_T_seconds_2 / 3600  # in ore
delta_T_hours_2 = delta_T_hours_2 + extra_unc*delta_T_hours_2

delta_T_seconds_3 = df_3 * T_seconds**2
delta_T_hours_3 = delta_T_seconds_3 / 3600  # in ore
delta_T_hours_3 = delta_T_hours_3 + extra_unc*delta_T_hours_3

# Plot
plt.figure(figsize=(8,5))
plt.plot(T_hours, delta_T_hours_1, color='tab:blue', label='10-days spt windows - hourly time-series')
#plt.plot(T_hours, delta_T_hours_2, color='navy', label='20-days spt windows - hourly time-series')
#plt.plot(T_hours, delta_T_hours_3, color='cyan', label='5-days spt windows - hourly time-series')
plt.xlabel('Modes Period [h]')
plt.ylabel('Period uncertainty [h]')
plt.title('Spectrum Modes Uncertainty')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig('T_uncertainty.png')
plt.close()

# Plot in %
plt.figure(figsize=(8,5))
plt.plot(T_hours, delta_T_hours_1*100/T_hours, color='tab:blue', label='10-days spt windows - hourly time-series')
#plt.plot(T_hours, delta_T_hours_2, color='navy', label='20-days spt windows - hourly time-series')
#plt.plot(T_hours, delta_T_hours_3, color='cyan', label='5-days spt windows - hourly time-series')
plt.xlabel('Modes Period [h]')
plt.ylabel('Period uncertainty [%]')
plt.title('Spectrum Modes % Uncertainty')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig('T_uncertainty_perc.png')
plt.close()


### Plot incertezza + effetto aggiunta Coriolis
plt.figure(figsize=(8, 5))
plt.plot(T_hours, delta_T_hours_1, color='tab:blue', label='10-days spt windows - hourly time-series')
plt.plot(T_hours_bk, delta_T_low / 3600, color='tab:orange', label='South Med')
plt.plot(T_hours_bk, delta_T_high / 3600, color='red', label='North Med')
#plt.plot(T_hours_bk, (delta_T_high-delta_T_low)/3600,'k--', label='North Med - South Med') 
plt.axhline(y=0, color='k', linestyle='--', label='')
plt.xlim(0, 40)
plt.xlabel('Modes Period (hours)')
plt.ylabel('Absolute period shift / Period Uncertainty [hours]')
plt.title('Effect of Coriolis on Mode Period Vs Spt uncertainty')
plt.plot(T_hours, delta_T_hours_1, color='tab:blue', label='10-days spt windows - hourly time-series')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig('coriolis_shift_spt_uncertainty.png')
plt.close()

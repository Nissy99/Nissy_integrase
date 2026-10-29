import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Parameters
promoter_strength = 24  # µM min^-1
initial_triggers = [3.6, 1.8, 0.72, 0.072, 0.0072, 0.00072, 0]  # µM
K_m = 3.0e3  # uM
initial_dna_int = 4.25e-3  # µM  Dv = 4.25e-4  B = 2.12e-3  Sc = 4.25e-3
initial_dna_luc = 3.3e-3  # µM
translation_rate_luciferase = 1.0e5  # au
mrna_degradation_rate = 6.25e-2  # min^-1
flipping_rate = 0.001666  # min^-1
translation_rate_integrase = 10  # au
leakage_rate = 0.001  # min^-1 switch leakage
K_ta = 10  # min^-1 trigger activation rate

total_time = 200  # min

# List DataFramesfor each simulation
all_data = []

# Initial conditions
for initial_trigger in initial_triggers:
    mrna_integrase = 0
    protein_integrase = 0
    mrna_luciferase = 0
    protein_luciferase = 0
    time = 0
    time_points = [time]
    mrna_integrase_levels = [mrna_integrase]
    protein_integrase_levels = [protein_integrase]
    mrna_luciferase_levels = [mrna_luciferase]
    protein_luciferase_levels = [protein_luciferase]


    # Calculate events
    def calculate_rates(mrna_integrase, protein_integrase, mrna_luciferase, protein_luciferase):
        transcription_rate_integrase = promoter_strength * initial_dna_int

        if mrna_integrase <= initial_trigger:
            translation_integrase = translation_rate_integrase * K_ta * mrna_integrase * initial_trigger
            translation_leakage = translation_rate_integrase * leakage_rate * mrna_integrase
            translation_rate_eff_integrase = translation_leakage + translation_integrase
        else:
            translation_integrase = translation_rate_integrase * K_ta * initial_trigger * initial_trigger
            translation_leakage = translation_rate_integrase * leakage_rate * mrna_integrase
            translation_rate_eff_integrase = translation_leakage + translation_integrase

        degradation_rates_integrase = mrna_degradation_rate * mrna_integrase

        transcription_rate_luciferase = promoter_strength * initial_dna_luc
        mrna_luciferase_on = mrna_luciferase * flipping_rate * (protein_integrase / 4)
        translation_rate_eff_luciferase = translation_rate_luciferase * mrna_luciferase_on
        degradation_rates_luciferase = mrna_degradation_rate * mrna_luciferase

        return (transcription_rate_integrase, translation_rate_eff_integrase, degradation_rates_integrase,
                transcription_rate_luciferase, translation_rate_eff_luciferase, degradation_rates_luciferase)


    # Stochastic simulation - Gillespie
    while time < total_time:
        (transcription_rate_integrase, translation_rate_eff_integrase, degradation_rates_integrase,
         transcription_rate_luciferase, translation_rate_eff_luciferase,
         degradation_rates_luciferase) = calculate_rates(mrna_integrase, protein_integrase, mrna_luciferase,
                                                         protein_luciferase)

        total_rate = (transcription_rate_integrase + translation_rate_eff_integrase + degradation_rates_integrase +
                      transcription_rate_luciferase + translation_rate_eff_luciferase + degradation_rates_luciferase)

        if total_rate == 0:
            break

        time += 1

        mrna_integrase += (transcription_rate_integrase * K_m) / (transcription_rate_integrase + K_m)
        protein_integrase += (translation_rate_eff_integrase * K_m) / (translation_rate_eff_integrase + K_m)
        mrna_integrase -= degradation_rates_integrase
        mrna_luciferase += (transcription_rate_luciferase * K_m) / (transcription_rate_luciferase + K_m)
        protein_luciferase += (translation_rate_eff_luciferase * K_m) / (translation_rate_eff_luciferase + K_m)
        mrna_luciferase -= degradation_rates_luciferase

        time_points.append(time)
        mrna_integrase_levels.append(mrna_integrase)
        protein_integrase_levels.append(protein_integrase)
        mrna_luciferase_levels.append(mrna_luciferase)
        protein_luciferase_levels.append(protein_luciferase)

    data = {
        'Time (min)': time_points,
        'Integrase mRNA Count': mrna_integrase_levels,
        'Integrase Protein Count': protein_integrase_levels,
        'luciferase mRNA Count': mrna_luciferase_levels,
        'luciferase Protein Count': protein_luciferase_levels,
        'Initial Trigger (uM)': [initial_trigger] * len(time_points)
    }
    df = pd.DataFrame(data)
    all_data.append(df)

# Combine DataFrames
combined_df = pd.concat(all_data, ignore_index=True)

# Save results as CSV
combined_df.to_csv('expression_data_combined.csv', index=False)

# Plot combine results
plt.figure(figsize=(6, 5))
for initial_trigger in initial_triggers:
    subset = combined_df[combined_df['Initial Trigger (uM)'] == initial_trigger]
    plt.step(subset['Time (min)'], subset['luciferase Protein Count'], where='post', label=f'{initial_trigger} µM')
plt.xlabel('Time (min)')
plt.ylabel('Luciferase (au)')
plt.title('Switch-Integrase')
plt.legend()
plt.xlim(0, 200)
plt.savefig('expression_plot_combined_integrase.pdf')
plt.show()

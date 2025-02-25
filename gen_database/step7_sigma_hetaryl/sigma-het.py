# calculate Pka
import pandas as pd

# Constants 
R = 0.001987204  # Gas constant in kcal/mol/K
T = 298.15  # Temperature in Kelvin
conversion_factor = 627.5094741  # Hartree to kcal/mol
ln_conversion = 2.30258509  # Natural log to log base 10 conversion

df = pd.read_csv('energies.csv')

if 'G(T)_SPC_c' in df.columns and 'G(T)_SPC_b' in df.columns:
    # new 'pKa' column
    df['pKa'] = (df['G(T)_SPC_c'] * conversion_factor - 270.29 - df['G(T)_SPC_b'] * conversion_factor) / (ln_conversion * R * T)
else:
    print("Error: Required columns are missing")

if "Structure" in df.columns:
    reference_pka = df[df['Structure'] == '24_1_0_1']['pKa'].iloc[0]  # pKa of reference 
    df['sigma_het'] = reference_pka - df['pKa']
    df.to_csv('sigma_het.csv', index=False)
    print("sigma_het.csv")
else:
    print("Error: 'Structure' column is missing")
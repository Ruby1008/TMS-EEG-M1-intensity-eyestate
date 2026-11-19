"""Fig 4 bandwidth per IMF (decision 2026-10-16: IMF4 definition = peak and FWHM; 10-90th percentile in the Fig 4 caption).
Input: data/Fig4_IMF_frequency_distribution_main17.csv (probability per 1/8-octave bin; imf_frequency_distribution.m) and
data/Fig4_IMF_peak_mean_Hz_main17.csv. FWHM = lowest and highest bin centre with probability >= half the maximum.
Output: data/Fig4_IMF_bandwidth_main17.csv"""
import numpy as np, pandas as pd
from pathlib import Path
D = Path(__file__).resolve().parents[1] / "data"
H = pd.read_csv(D / "Fig4_IMF_frequency_distribution_main17.csv"); PK = pd.read_csv(D / "Fig4_IMF_peak_mean_Hz_main17.csv").set_index("IMF")
f = H.freq_Hz.values; rows = []
for j in range(1, 11):
    c = H[f"IMF{j}"].values; cs = np.cumsum(c) / c.sum(); q = lambda x: f[np.searchsorted(cs, x)]
    rows.append(dict(IMF=j, peak_Hz=PK.loc[j, "peak_Hz"], mean_Hz=PK.loc[j, "mean_Hz"], fwhm_low_Hz=f[c >= c.max() / 2].min(), fwhm_high_Hz=f[c >= c.max() / 2].max(),
                     p10_Hz=q(.10), p90_Hz=q(.90), p025_Hz=q(.025), p975_Hz=q(.975)))
B = pd.DataFrame(rows); B.to_csv(D / "Fig4_IMF_bandwidth_main17.csv", index=False); print(B.round(1).to_string())

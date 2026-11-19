"""N1 (user 2026-10-28): eye-state tests on the conventional TEP (0.1-30 Hz, _F_TEP), ROI, 17 participants, real stimulation.
(a) 2 x 3 rmANOVA (eye x intensity) per component, conventional-TEP windows (P30 28-48, N45 41-61, P60 48-68, N100 120-140 ms; same
    windows and same averaging as code/tep_stats.py, which these values reproduce: data/TEP_component_anova.csv).
(b) EC - EO per component x intensity: mean, SD, dz, paired t, uncorrected p (12 rows).
(c) Sliding-window rmANOVA (21-sample windows starting at 0, 10, ..., 390 ms, as tep_stats.py) eye and interaction, BH-FDR over 40 windows.
(d) 120% RMT N100 EC - EO: mean, SD, the 17 individual values.
Outputs: N1/N1_conventional_TEP_eye.csv, N1/N1_anova.csv, N1/N1_slidingwindow.csv, N1/N1_N100_120_individual.csv, N1/N1_text.md"""
import numpy as np, pandas as pd
import r1110_lib as L

OUT = L.RV / "N1"
rows_a, rows_b = [], []
for comp, (a0, a1) in L.WIN_TEP.items():
    df = pd.DataFrame([dict(sub=s, eye=e, intensity=lv, y=L.wmean(*L.tep_wave(s, e, lv), a0, a1)) for s in L.S17 for e in L.EYES for lv in L.LEV])
    A = L.anova2(df); A.insert(0, "measure", f"{comp} {a0}-{a1} ms"); rows_a.append(A)
    for lv in L.LEV:
        P = df[df.intensity == lv].pivot(index="sub", columns="eye", values="y")
        r = L.paired(P["EC"], P["EO"])
        rows_b.append(dict(measure=f"{comp} {a0}-{a1} ms", intensity=lv, n=r["n"], mean_EO=P["EO"].mean(), mean_EC=P["EC"].mean(),
                           mean_EC_minus_EO=r["mean_diff"], SD_diff=r["SD_diff"], dz=r["dz"], t=r["t"], df=r["df"], p=r["p"]))
AN = pd.concat(rows_a, ignore_index=True)
AN[["measure", "effect", "F", "df1", "df2", "p", "eta_p2", "eps"]].to_csv(OUT / "N1_anova.csv", index=False)
B = pd.DataFrame(rows_b); B.to_csv(OUT / "N1_conventional_TEP_eye.csv", index=False)
# reproduction check against the main-analysis file
OLD = pd.read_csv(L.DATA / "TEP_component_anova.csv"); OLD = OLD[OLD.data == "TEP"].set_index("component")
for comp in L.WIN_TEP:
    a = AN[(AN.measure.str.startswith(comp + " ")) & (AN.effect == "intensity")].iloc[0]
    assert abs(a.F - OLD.loc[comp, "F_intensity"]) < 1e-9 and abs(a.p - OLD.loc[comp, "p_intensity"]) < 1e-9, comp
    e = AN[(AN.measure.str.startswith(comp + " ")) & (AN.effect == "eye")].iloc[0]; assert abs(e.F - OLD.loc[comp, "F_eye"]) < 1e-9
# (c) sliding window
starts = np.arange(0, 400, 10); sw = []
for a in starts:
    df = pd.DataFrame([dict(sub=s, eye=e, intensity=lv, y=L.winmean21(*L.tep_wave(s, e, lv), a)) for s in L.S17 for e in L.EYES for lv in L.LEV])
    A = L.anova2(df).set_index("effect")
    sw.append(dict(win_start_ms=a, win_end_ms=a + 20, **{f"{c}_{k}": A.loc[eff, c] for k, eff in (("eye", "eye"), ("intensity", "intensity"), ("interaction", "eye x intensity")) for c in ("F", "df1", "df2", "p")}))
SW = pd.DataFrame(sw)
for k in ("eye", "intensity", "interaction"): SW[f"q_{k}"] = L.bh(SW[f"p_{k}"])
SW.to_csv(OUT / "N1_slidingwindow.csv", index=False)
OLDSW = pd.read_csv(L.DATA / "TEP_slidingwindow_anova_TEP.csv")
assert np.allclose(SW.p_eye, OLDSW.p_eye) and np.allclose(SW.p_interaction, OLDSW.p_interaction)
# (d) individual EC - EO at 120% N100
df = pd.DataFrame([dict(sub=s, eye=e, y=L.wmean(*L.tep_wave(s, e, 120), 120, 140)) for s in L.S17 for e in L.EYES]).pivot(index="sub", columns="eye", values="y")
df["EC_minus_EO"] = df["EC"] - df["EO"]; df = df.reset_index().rename(columns={"sub": "participant"}); df["participant"] = df.participant.map(lambda s: f"{s:03d}")
df.to_csv(OUT / "N1_N100_120_individual.csv", index=False)
d = df.EC_minus_EO
print(AN.round(3).to_string()); print(B.round(3).to_string())
print("sliding: min q_eye", SW.q_eye.min().round(3), "min q_int", SW.q_interaction.min().round(3), "| uncorrected eye p<.05:", list(SW.win_start_ms[SW.p_eye < .05]),
      "int p<.05:", list(SW.win_start_ms[SW.p_interaction < .05]))
print("120 N100 EC-EO mean", d.mean().round(3), "SD", d.std(ddof=1).round(3), "median", d.median().round(3), "n pos", (d > 0).sum(), "range", d.min().round(2), d.max().round(2))
print("EO mean", df.EO.mean().round(3), "EC mean", df.EC.mean().round(3))
print(df.round(3).to_string())

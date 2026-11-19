"""Round 8 analyses (user 2026-10-20 evening, section 1), part A (Python). All pre-specified, all reported; facts only.
1-1  Pre-stimulus ITPC without baseline correction (-100..-10 ms; ROI; IMF6, IMF7): 2 x 3 rmANOVA (GG), 17 participants; uncorrected ITPC = baseline-corrected ROI ITPC
     (data/ITPC_roi_real.mat) + ROI mean baseline ITPC (data/ITPC_roi_real_bl.mat, itpc_bl_roi_export.m). Finite-trial bias of ITPC under random phase: sqrt(pi / (4 n)).
3a   Manipulation check: pre-stimulus 8-12 Hz log10 power (Welch 500-sample Hann, -1000..-10 ms; data/Alpha_prestim_values.csv, intensities pooled); EO vs EC paired t + dz per
     region (occipital O1/O2/POz, sensorimotor FC3/FC1/C3/C1), Holm over the 2 regions per stimulation type; region x state rmANOVA (interaction); real and sham; IMF6 version listed.
3b   Volume-conduction check: surface Laplacian (spherical-spline CSD, Perrin et al. 1989; m = 4, 50 Legendre terms, lambda = 1e-5; 63 EEG channels, positions from
     waveguard64_woRef.ced) applied to -1000..-10 ms of every trial of _F_ERM; mu power (8-12 Hz, same Welch) at C3 (main) and the ROI mean; EO vs EC before vs after the Laplacian.
3c   (i) existing rmANOVAs (main 17 / MEP 22); (iii) whole-scalp cluster-based permutation EO vs EC on the ERM (0-400 ms, 61 channels, intensities pooled; paired t, two-tailed;
     5000 sign-flip permutations, seed 20261021; cluster-forming p < .05; neighbours = 3-D distance < 1.6 x median nearest-neighbour distance, adjacent samples; mass = sum t);
     (iv) sensitivity: participants without any 20-50 uV residual condition (n = 10); (ii) and (v) in R (eyestate_null_round8.R) on the tables written here.
3d   TOST bounds: the original TOST (dz = +-0.5; data/TOST_EO_EC.csv) and, per measure, the smallest equivalence bound that would pass TOST at alpha = .05 (= |mean| + t(.95, df) x SE),
     in raw units and dz, also at alpha = .05/21 (Bonferroni over the 21 tests).
3f   Individual differences (exploratory): per participant EC - EO change (intensities pooled) of Laplacian mu at C3 vs the change of log MEP, ERM N100 and IMF6 ITPC (longest window);
     Spearman, Holm over 3.
Outputs: data/EyeState_manipulation.csv, data/EyeState_null_robustness.csv (rows for 1-1 and 3c i/iii/iv; R appends ii/v/e), data/EyeState_TOST_bounds.csv,
         data/EyeState_individual_differences.csv, data/EyeState_cluster_perm.csv, data/EyeState_cluster_perm_tmap.npz, data/ITPC_prestim_uncorrected.csv,
         data/EyeState_allavailable_long.csv, data/EyeState_block_order.csv, data/EyeState_laplacian_values.csv"""
data_dir = r"path/to/data"   # release: set this path
import re, numpy as np, pandas as pd, pingouin as pg, scipy.io as sio
from scipy import stats, signal
from scipy.special import eval_legendre
from scipy.sparse import coo_matrix, csr_matrix, kron, identity
from scipy.sparse.csgraph import connected_components
import tep_io as io

D = io.DATA; S17 = io.MAIN17; ROWS = []
def add(**k): ROWS.append(k)
def rm(df, dv, within=("eye", "intensity")):
    a = pg.rm_anova(data=df, dv=dv, within=list(within), subject="sub", correction=True, effsize="np2").set_index("Source"); return a
def gg(a, src): r = a.loc[src]; return r["F"], r["ddof1"], r["ddof2"], (r["p_GG_corr"] if "p_GG_corr" in r and pd.notna(r.get("p_GG_corr")) else r["p_unc"]), r["np2"]
# ======================================================= 1-1 pre-stimulus ITPC without baseline correction
M = sio.loadmat(D / "ITPC_roi_real.mat", squeeze_me=True); B = sio.loadmat(D / "ITPC_roi_real_bl.mat", squeeze_me=True)
ids = [str(i) for i in M["ids"]]; idsB = [str(i) for i in B["ids"]]; assert ids == idsB; t = np.round(M["t"].astype(float)); idx = {i: k for k, i in enumerate(ids)}
pre = (t >= -100) & (t <= -10); rows = []
for s in S17:
    for e in io.EYES:
        k = idx[f"spTMS_{s:03d}_Real_{e}"]
        for li, lv in enumerate(io.LEV):
            for imf in (6, 7):
                corr = M["X"][k, imf - 1][pre, li].mean(); raw = corr + B["BL"][k, imf - 1, li]; n = B["nT"][k, li]
                rows.append(dict(sub=s, eye=e, intensity=lv, IMF=imf, itpc_corrected=corr, itpc_uncorrected=raw, n_trials=n, bias_random_phase=np.sqrt(np.pi / (4 * n))))
PI = pd.DataFrame(rows); PI.to_csv(D / "ITPC_prestim_uncorrected.csv", index=False)
for imf in (6, 7):
    g = PI[PI.IMF == imf]; a = rm(g, "itpc_uncorrected")
    for src in ("eye", "intensity", "eye * intensity"):
        F, d1, d2, p, np2 = gg(a, src)
        add(section="1-1", analysis="pre-stimulus ITPC (-100..-10 ms), no baseline correction, 2 x 3 rmANOVA (GG)", measure=f"IMF{imf}", effect=src, n=17, F=F, df1=d1, df2=d2, p=p, np2=np2)
    for (e, lv), h in g.groupby(["eye", "intensity"]):
        add(section="1-1", analysis="pre-stimulus ITPC, condition mean (no baseline correction) and random-phase bias sqrt(pi/4n)", measure=f"IMF{imf}", effect=f"{e} {lv}%", n=len(h),
            mean=h.itpc_uncorrected.mean(), se=h.itpc_uncorrected.std(ddof=1) / np.sqrt(len(h)), bias=h.bias_random_phase.mean(), n_trials_median=h.n_trials.median())
# ======================================================= 3a manipulation check (Welch / IMF6)
AV = pd.read_csv(D / "Alpha_prestim_values.csv"); AV = AV[AV.intensity.astype(str) == "all"]
MAN = []
for meth in ("Welch", "IMF6"):
    for stim in ("Real", "Sham"):
        g = AV[(AV.method == meth) & (AV.stim == stim)]; tests = []
        for reg in ("occipital", "sensorimotor"):
            w = g[g.region == reg].pivot(index="sub", columns="eye", values="log10_power").dropna(); d = w.EC - w.EO; r = stats.ttest_rel(w.EC, w.EO)
            tests.append(dict(check="3a", method=meth, stim=stim, region=reg, n=len(d), mean_EO=w.EO.mean(), mean_EC=w.EC.mean(), diff_EC_minus_EO=d.mean(), t=r.statistic, df=len(d) - 1, p=r.pvalue, dz=d.mean() / d.std(ddof=1)))
        ph = pg.multicomp([x["p"] for x in tests], method="holm")[1]
        for x, h in zip(tests, ph): x["p_holm_2"] = h; MAN.append(x)
        w2 = g[g.region.isin(["occipital", "sensorimotor"])]; subs_ok = w2.groupby("sub").size(); w2 = w2[w2["sub"].isin(subs_ok[subs_ok == 4].index)]
        a = pg.rm_anova(data=w2, dv="log10_power", within=["region", "eye"], subject="sub").set_index("Source"); r = a.loc["region * eye"]
        MAN.append(dict(check="3a", method=meth, stim=stim, region="region x state interaction (rmANOVA)", n=w2["sub"].nunique(), F=r["F"], df=f"{int(r['ddof1'])}, {int(r['ddof2'])}", p=r["p_unc"], np2=r["np2"] if "np2" in r else np.nan))
# ======================================================= 3b surface Laplacian
CED = pd.read_csv(data_dir + r"\waveguard64_woRef.ced", sep="\t", skipinitialspace=True); CED.columns = [c.strip() for c in CED.columns]; CED["labels"] = CED["labels"].str.strip()
def csd_matrix(labs, m=4, nterms=50, lam=1e-5):
    P = np.array([[float(CED[CED.labels == l].X.iloc[0]), float(CED[CED.labels == l].Y.iloc[0]), float(CED[CED.labels == l].Z.iloc[0])] for l in labs]); P /= np.linalg.norm(P, axis=1, keepdims=True)
    cosang = np.clip(P @ P.T, -1, 1); n = np.arange(1, nterms + 1)
    Ln = np.stack([eval_legendre(k, cosang) for k in n])                                 # nterms x ch x ch
    G = np.tensordot((2 * n + 1) / (n * (n + 1)) ** m, Ln, 1) / (4 * np.pi); H = -np.tensordot((2 * n + 1) / (n * (n + 1)) ** (m - 1), Ln, 1) / (4 * np.pi)
    Gi = np.linalg.inv(G + lam * np.eye(len(labs))); s = Gi.sum(0); C = Gi - np.outer(Gi.sum(1), s) / s.sum()      # c = Gi (V - c0), c0 = sum(Gi V) / sum(Gi)
    return H @ C
LAP = []; Tcache = {}
for stim in ("Real", "Sham"):
    for s in S17:
        for e in io.EYES:
            id_ = f"spTMS_{s:03d}_{stim}_{e}"
            if not (io.FINAL / f"{id_}_F_ERM.set").exists(): continue
            x, labs, tt, inten, _ = io.load_set(io.FINAL, id_ + "_F_ERM"); tt = np.round(tt).astype(int); w = (tt >= -1000) & (tt <= -10)
            key = tuple(labs)
            if key not in Tcache: Tcache[key] = csd_matrix(labs)
            seg = x[:, w, :].astype(float); csd = np.einsum("ij,jtk->itk", Tcache[key], seg)
            for ver, dat in (("before Laplacian", seg), ("after Laplacian (CSD)", csd)):
                f, Pw = signal.welch(dat, fs=1000, window="hann", nperseg=500, noverlap=250, axis=1); fa = (f >= 8) & (f <= 12)
                for reg, chs in (("C3", ["C3"]), ("ROI", io.ROI), ("occipital", ["O1", "O2", "POz"])):
                    ci = [labs.index(c) for c in chs]; pw = Pw[ci][:, fa].mean((0, 1)).mean()
                    LAP.append(dict(stim=stim, sub=s, eye=e, version=ver, region=reg, log10_power=np.log10(pw)))
            print(id_, flush=True)
LAPd = pd.DataFrame(LAP); LAPd.to_csv(D / "EyeState_laplacian_values.csv", index=False)
for stim in ("Real", "Sham"):
    for ver in ("before Laplacian", "after Laplacian (CSD)"):
        for reg in ("C3", "ROI", "occipital"):
            w = LAPd[(LAPd.stim == stim) & (LAPd.version == ver) & (LAPd.region == reg)].pivot(index="sub", columns="eye", values="log10_power").dropna(); d = w.EC - w.EO; r = stats.ttest_rel(w.EC, w.EO)
            MAN.append(dict(check="3b", method=f"Welch, {ver}", stim=stim, region=reg, n=len(d), mean_EO=w.EO.mean(), mean_EC=w.EC.mean(), diff_EC_minus_EO=d.mean(), t=r.statistic, df=len(d) - 1, p=r.pvalue, dz=d.mean() / d.std(ddof=1)))
pd.DataFrame(MAN).to_csv(D / "EyeState_manipulation.csv", index=False)
# ======================================================= data tables: main 17, all available, sensitivity set
MEPm = pd.read_csv(D / "MEP_matrix_new.csv").set_index("sub")
WIN = {"P30": (27, 47), "N45": (39, 59), "P60": (50, 70), "N100": (120, 140)}
WN = pd.read_csv(D / "Fig8_ITPC_windows.csv").set_index("IMF"); IW = {imf: (WN.loc[imf, "longest_start_ms"], WN.loc[imf, "longest_end_ms"]) for imf in (6, 7)}
long = []
for e in io.EYES:
    for s in sorted(set(io.AVAIL["EO"]) | set(io.AVAIL["EC"])):
        if s not in io.AVAIL[e]: continue
        R = io.load_erm(s, e); te = np.round(R["time"])
        k = idx[f"spTMS_{s:03d}_Real_{e}"]
        for li, lv in enumerate(io.LEV):
            for c, (a0, a1) in WIN.items(): long.append(dict(sub=s, eye=e, intensity=lv, measure=f"ERM {c}", value=R["erm_roi"][li][(te >= a0) & (te <= a1)].mean()))
            for imf in (6, 7): long.append(dict(sub=s, eye=e, intensity=lv, measure=f"ITPC IMF{imf}", value=M["X"][k, imf - 1][(t >= IW[imf][0]) & (t <= IW[imf][1]), li].mean()))
for s in MEPm.index:
    for e in io.EYES:
        for lv in io.LEV:
            v = MEPm.loc[s, f"{e}_{lv}"]
            if pd.notna(v): long.append(dict(sub=s, eye=e, intensity=lv, measure="log MEP", value=np.log(v)))
L = pd.DataFrame(long); L["main17"] = L["sub"].isin(S17); L.to_csv(D / "EyeState_allavailable_long.csv", index=False)
_q = io.Q[(io.Q.status == "processed") & (io.Q.stim == "Real") & (io.Q.stageD == "sensitivity")]; SENS = [s for s in S17 if s not in set(_q["sub"])]
MEP22 = list(MEPm.dropna().index)
# block order from the recording start times in the .cnt file names
DL = pd.read_csv(io.REV / "Preproc" / "datalist_Real.csv", encoding="utf-8-sig"); ords = []
for s in sorted(set(DL["sub"])):
    tm = {}
    for e in io.EYES:
        r = DL[(DL["sub"] == s) & (DL.cond == f"Real_{e}")]
        if len(r): m_ = re.search(r"(\d{4}-\d{2}-\d{2}_\d{2}-\d{2}-\d{2})", str(r.cnt.iloc[0])); tm[e] = m_.group(1) if m_ else None
    if tm.get("EO") and tm.get("EC"): ords.append(dict(sub=s, EO_start=tm["EO"], EC_start=tm["EC"], order="EO first" if tm["EO"] < tm["EC"] else "EC first"))
ORD = pd.DataFrame(ords); ORD.to_csv(D / "EyeState_block_order.csv", index=False)
# ======================================================= 3c (i) and (iv): rmANOVA of eye state, main and sensitivity sets
for label, subs_eeg, subs_mep in (("(i) main analysis", S17, MEP22), ("(iv) sensitivity, no 20-50 uV residual", SENS, [s for s in MEP22 if s in SENS])):
    for meas in ["log MEP", "ERM P30", "ERM N45", "ERM P60", "ERM N100", "ITPC IMF6", "ITPC IMF7"]:
        subs = subs_mep if meas == "log MEP" else subs_eeg
        g = L[(L.measure == meas) & L["sub"].isin(subs)]; full = g.groupby("sub").size(); g = g[g["sub"].isin(full[full == 6].index)]
        a = rm(g, "value")
        for src in ("eye", "eye * intensity"):
            F, d1, d2, p, np2 = gg(a, src); add(section="3c", analysis=label, measure=meas, effect=src, n=g["sub"].nunique(), F=F, df1=d1, df2=d2, p=p, np2=np2)
# ======================================================= 3c (iii): cluster-based permutation EO vs EC, ERM 0-400 ms, 61 channels
E0 = io.load_erm(S17[0], "EO"); te = np.round(E0["time"]).astype(int); labs61 = [str(x).strip() for x in E0["chanlabels"]]; tw = (te >= 0) & (te <= 400)
Dd = np.stack([io.load_erm(s, "EO")["erm_allch"].mean(1)[:, tw] - io.load_erm(s, "EC")["erm_allch"].mean(1)[:, tw] for s in S17])   # n x ch x time (EO - EC, intensities pooled)
n, nC, nT = Dd.shape
P3 = np.array([[float(CED[CED.labels == l].X.iloc[0]), float(CED[CED.labels == l].Y.iloc[0]), float(CED[CED.labels == l].Z.iloc[0])] for l in labs61])
Dm = np.linalg.norm(P3[:, None] - P3[None], axis=-1); np.fill_diagonal(Dm, np.inf); adjC = (Dm < 1.6 * np.median(Dm.min(1))).astype(float)
adjT = coo_matrix((np.ones(2 * (nT - 1)), (np.r_[np.arange(nT - 1), np.arange(1, nT)], np.r_[np.arange(1, nT), np.arange(nT - 1)])), shape=(nT, nT))
A = (kron(csr_matrix(adjC), identity(nT)) + kron(identity(nC), adjT)).tocsr()                          # node = ch * nT + time
tcrit = stats.t.isf(.025, n - 1)
def tmap(X): m = X.mean(0); sd = X.std(0, ddof=1); return m / (sd / np.sqrt(n))
def clus(T):
    out = []
    for sgn in (1, -1):
        mask = (sgn * T).ravel() > tcrit
        if not mask.any(): continue
        nodes = np.where(mask)[0]; ncomp, labl = connected_components(A[nodes][:, nodes], directed=False)
        mass = np.bincount(labl, weights=T.ravel()[nodes], minlength=ncomp)
        for c in range(ncomp): out.append((mass[c], nodes[labl == c]))
    return out
Tobs = tmap(Dd); obs = clus(Tobs); rng = np.random.default_rng(20261021); null = np.zeros(5000)
for k in range(5000):
    fl = rng.choice([-1.0, 1.0], n)[:, None, None]; cl = clus(tmap(Dd * fl)); null[k] = max([abs(m_) for m_, _ in cl], default=0)
CL = []
for m_, nodes in sorted(obs, key=lambda z: -abs(z[0])):
    ch = sorted({labs61[i // nT] for i in nodes}); tms = te[tw][nodes % nT]
    CL.append(dict(sign="EO > EC" if m_ > 0 else "EO < EC", mass_t=m_, p=(1 + (null >= abs(m_)).sum()) / 5001, start_ms=int(tms.min()), end_ms=int(tms.max()), n_channels=len(ch), channels=", ".join(ch), n_nodes=len(nodes)))
CLd = pd.DataFrame(CL, columns=["sign", "mass_t", "p", "start_ms", "end_ms", "n_channels", "channels", "n_nodes"]); CLd.to_csv(D / "EyeState_cluster_perm.csv", index=False)
np.savez_compressed(D / "EyeState_cluster_perm_tmap.npz", t=Tobs, labels=np.array(labs61), time=te[tw])
add(section="3c", analysis="(iii) cluster-based permutation EO vs EC, ERM 0-400 ms, 61 channels, intensities pooled (5000)", measure="ERM, all channels", effect="eye", n=n,
    n_clusters=len(CLd), min_cluster_p=CLd.p.min() if len(CLd) else np.nan, note="; ".join(f"{r.sign} {r.start_ms}-{r.end_ms} ms, {r.n_channels} ch, p = {r.p:.3f}" for r in CLd.head(4).itertuples()) or "no cluster formed")
print("clusters:", CLd.head(6).to_string())
# ======================================================= 3d TOST bounds
TO = pd.read_csv(D / "TOST_EO_EC.csv"); tb = []
for r in TO.itertuples():
    se = r.sd_diff / np.sqrt(r.n); df = r.n - 1
    for alpha, lab in ((.05, ".05"), (.05 / 21, ".05/21")):
        b = abs(r.mean_diff_EO_minus_EC) + stats.t.ppf(1 - alpha, df) * se
        tb.append(dict(measure=r.measure, intensity=r.intensity, unit=r.unit, n=r.n, mean_diff_EO_minus_EC=r.mean_diff_EO_minus_EC, sd_diff=r.sd_diff, dz=r.dz, alpha=lab,
                       original_bound_dz=0.5, original_bound_raw=r.bound_high, p_tost_original=r.p_tost, p_tost_holm_original=r.p_tost_holm, smallest_passing_bound_raw=b, smallest_passing_bound_dz=b / r.sd_diff))
TB = pd.DataFrame(tb); TB.to_csv(D / "EyeState_TOST_bounds.csv", index=False)
# ======================================================= 3f individual differences
lap = LAPd[(LAPd.stim == "Real") & (LAPd.version == "after Laplacian (CSD)") & (LAPd.region == "C3")].pivot(index="sub", columns="eye", values="log10_power"); dmu = lap.EC - lap.EO
def dmeas(meas):
    g = L[(L.measure == meas) & L["sub"].isin(S17)].pivot_table(index=["sub", "intensity"], columns="eye", values="value"); return (g.EC - g.EO).groupby("sub").mean()
ID = []
for meas in ("log MEP", "ERM N100", "ITPC IMF6"):
    y = dmeas(meas); j = dmu.index.intersection(y.dropna().index); r = stats.spearmanr(dmu.loc[j], y.loc[j]); ID.append(dict(x="EC - EO Laplacian mu (C3)", y=f"EC - EO {meas}", n=len(j), rho=r.statistic, p=r.pvalue))
ID = pd.DataFrame(ID); ID["p_holm"] = pg.multicomp(ID.p.values, method="holm")[1]; ID.to_csv(D / "EyeState_individual_differences.csv", index=False)
for r in ID.itertuples(): add(section="3f", analysis="individual differences (exploratory), Spearman, Holm over 3", measure=r.y, effect="rho", n=r.n, rho=r.rho, p=r.p, p_holm=r.p_holm)
pd.DataFrame(ROWS).to_csv(D / "EyeState_null_robustness.csv", index=False)
pd.set_option("display.width", 230); print(pd.DataFrame(MAN).round(4).to_string(index=False)); print(pd.DataFrame(ROWS).round(4).to_string(index=False)); print(ORD.order.value_counts())
print(TB[TB.alpha == ".05"][["measure", "intensity", "dz", "smallest_passing_bound_dz", "smallest_passing_bound_raw"]].round(3).to_string(index=False))

"""N2 (user 2026-10-28): ITPC of the sham recordings, with the main-analysis ITPC rule v2 (code/itpc_v2.py): same layer-1 CEEMDAN phase
(HHSA/itpc_<id>.mat from itpc_extract), ROI FC3/FC1/C3/C1 (itpc_roi_export_1110.m -> _check/ITPC_roi_sham.mat; the same export of the real data
reproduces data/ITPC_roi_real.mat exactly), baseline -500..-100 ms (subtracted from the exported -500..-10 baseline-corrected ITPC, as itpc_v2.py),
IMF3-7, tests at 10..500 ms, BH-FDR q = .05 with one family per effect over IMF3-7 x 10..500 ms (rm2/bh/segs_of taken from itpc_round15.py as
itpc_v2.py does).
(a) sham 2 x 3 rmANOVA maps (participants with both sham eye states among the 17), intensity FDR mask, Jaccard real vs sham (IMF6, IMF7).
(b) summary ITPC per participant x condition: IMF6 27-456 ms, IMF7 33-471 ms (the real-data longest intensity segments, data/ITPC_v2_windows.csv).
(c) stimulation x intensity x eye rmANOVA (GG; participants with all 12 cells) on (b).
(d) real vs sham paired t per eye x intensity (Holm over the six tests, per IMF).
(e) sham intensity post hoc (eyes pooled; Holm over three contrasts, per IMF).
(f) FigureS_shamITPC_1110: layout of Figure 5B (IMF6, IMF7), real solid, sham dashed, 80/100/120% colours of Figure 5B, mean +- SE.
Outputs (N2/): N2_sham_ITPC_summary.csv, N2_sham_ITPC_tests.csv, N2_jaccard.csv, N2_sham_segments.csv, FigureS_shamITPC_1110.png/.pdf"""
import itertools
import numpy as np, pandas as pd, scipy.io as sio
from scipy import stats
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import r1110_lib as L

MFC = L.MF / "code"; OUT = L.RV / "N2"
_src = (MFC / "itpc_round15.py").read_text(encoding="utf-8"); ns = {"__file__": str(MFC / "itpc_round15.py")}; exec(_src.split("OLD = analyse(")[0], ns)
rm2, bh, segs_of, subs17 = ns["rm2"], ns["bh"], ns["segs_of"], ns["subs"]
assert list(subs17) == list(L.S17)
IMFS = [3, 4, 5, 6, 7]; T0 = 10; WIN = {6: (27, 456), 7: (33, 471)}

def load(path):
    M = sio.loadmat(path, squeeze_me=True); X = M["X"]; ids = [str(i) for i in M["ids"]]; t = np.round(M["t"].astype(float))
    return X - X[:, :, (t >= -500) & (t <= -100), :].mean(2, keepdims=True), {i: k for k, i in enumerate(ids)}, t
XR, iR, t = load(L.DATA / "ITPC_roi_real.mat"); XS, iS, tS = load(L.CK / "ITPC_roi_sham.mat"); assert np.array_equal(t, tS)
QQ = L.io.Q; shamok = QQ[(QQ.status == "processed") & (QQ.stim == "Sham") & (QQ.stageD != "exclude")]
S_EO = sorted(set(shamok[shamok.eyes == "EO"]["sub"]) & set(L.S17)); S_EC = sorted(set(shamok[shamok.eyes == "EC"]["sub"]) & set(L.S17))
S16 = sorted(set(S_EO) & set(S_EC)); print("sham EO", len(S_EO), "EC", len(S_EC), "both", len(S16), "missing EO:", sorted(set(L.S17) - set(S_EO)))
tw = (t >= -100) & (t <= 500); tt = t[tw]; test = tt >= T0

def cube(X, idx, stim, subjects):   # subjects x eye(EO, EC) x intensity x IMF x time
    Y = np.stack([np.stack([X[idx[f"spTMS_{s:03d}_{stim}_{e}"]][np.array(IMFS) - 1][:, tw, :] for e in ("EO", "EC")], 0) for s in subjects], 0)
    return np.moveaxis(Y, 4, 2)
def maps(Y):
    R = rm2(Y)
    for k in R:
        q = np.full_like(R[k]["p"], np.nan); q[:, test] = bh(R[k]["p"][:, test]); R[k]["q"] = q; R[k]["sig"] = np.nan_to_num(q, nan=1) < .05
    return R
YR17 = cube(XR, iR, "Real", L.S17); YR16 = cube(XR, iR, "Real", S16); YS16 = cube(XS, iS, "Sham", S16)
RR17, RR16, RS = maps(YR17), maps(YR16), maps(YS16)
# check: the real 17 intensity mask equals the main analysis (data/ITPC_v2_real_significant_segments.csv)
segR = pd.read_csv(L.DATA / "ITPC_v2_real_significant_segments.csv")
for i, imf in enumerate(IMFS):
    mine = [(tt[a], tt[b]) for a, b in segs_of(RR17["intensity"]["sig"][i])]
    ref = [(r.start_ms, r.end_ms) for r in segR[(segR.effect == "intensity") & (segR.IMF == imf)].itertuples()]
    assert mine == ref, (imf, mine[:3], ref[:3])
seg = []
for tag, R in (("sham_16", RS), ("real_16", RR16), ("real_17", RR17)):
    for k, v in R.items():
        for i, imf in enumerate(IMFS):
            for a, b in segs_of(v["sig"][i]): seg.append(dict(data=tag, effect=k, IMF=imf, start_ms=tt[a], end_ms=tt[b], n_points=b - a + 1, F_max=v["F"][i, a:b + 1].max()))
SEG = pd.DataFrame(seg); SEG.to_csv(OUT / "N2_sham_segments.csv", index=False)
J = []
for ref_tag, RRx in (("real_17 (main analysis)", RR17), ("real_16 (same participants)", RR16)):
    for imf in (6, 7):
        a = RRx["intensity"]["sig"][IMFS.index(imf)]; b = RS["intensity"]["sig"][IMFS.index(imf)]
        J.append(dict(real=ref_tag, IMF=imf, n_points_real=int(a.sum()), n_points_sham=int(b.sum()), n_both=int((a & b).sum()), jaccard=(a & b).sum() / (a | b).sum() if (a | b).any() else np.nan))
JD = pd.DataFrame(J); JD.to_csv(OUT / "N2_jaccard.csv", index=False)
SC = {k: {imf: int(RS[k]["sig"][IMFS.index(imf)].sum()) for imf in IMFS} for k in RS}

# (b) summary values
rows = []
for stim, X, idx, subjects in (("real", XR, iR, L.S17), ("sham", XS, iS, sorted(set(S_EO) | set(S_EC)))):
    for s in subjects:
        for e in ("EO", "EC"):
            key = f"spTMS_{s:03d}_{stim.capitalize()}_{e}"
            if key not in idx or (stim == "sham" and s not in (S_EO if e == "EO" else S_EC)): continue
            for k, lv in enumerate(L.LEV):
                for imf, (a0, a1) in WIN.items():
                    w = (tt >= a0) & (tt <= a1); rows.append(dict(participant=f"{s:03d}", stim=stim, eye=e, intensity=lv, IMF=imf, window_ms=f"{a0}-{a1}",
                                                                  ITPC_mean=X[idx[key], imf - 1][tw][w, k].mean(), in_16=s in S16))
SUM = pd.DataFrame(rows); SUM.to_csv(OUT / "N2_sham_ITPC_summary.csv", index=False)

# (c) 3-way rmANOVA with GG (orthonormal contrasts)
def helmert(k):
    H = np.zeros((k, k - 1))
    for j in range(1, k): H[:j, j - 1] = 1; H[j, j - 1] = -j
    return H / np.linalg.norm(H, axis=0)
def rm_anova_n(Y, names):
    """Y: n x k1 x k2 x ... ; returns F, GG df, p, eta_p2, eps per effect (all subsets of factors)."""
    n = Y.shape[0]; ks = Y.shape[1:]; out = []
    for r in range(1, len(ks) + 1):
        for S in itertools.combinations(range(len(ks)), r):
            Z = Y.reshape(n, -1)
            Ms = [helmert(k) if f in S else np.ones((k, 1)) / np.sqrt(k) for f, k in enumerate(ks)]
            M = Ms[0]
            for m in Ms[1:]: M = np.kron(M, m)
            Z = Z @ M; df = Z.shape[1]; mz = Z.mean(0); Zc = Z - mz
            ss_eff = n * (mz ** 2).sum(); ss_err = (Zc ** 2).sum(); dfe = df * (n - 1)
            F = (ss_eff / df) / (ss_err / dfe); C = Zc.T @ Zc / (n - 1)
            eps = 1.0 if df == 1 else min(1.0, np.trace(C) ** 2 / (df * np.trace(C @ C)))
            out.append(dict(effect=" x ".join(names[f] for f in S), F=F, df1=df * eps, df2=dfe * eps, p=stats.f.sf(F, df * eps, dfe * eps), eta_p2=ss_eff / (ss_eff + ss_err), eps=eps))
    return pd.DataFrame(out)
# validation of rm_anova_n against pingouin (2-way, real ITPC IMF6 summary)
_d = SUM[(SUM.stim == "real") & (SUM.IMF == 6)].rename(columns={"participant": "sub"})
_pg = L.anova2(_d, dv="ITPC_mean").set_index("effect")
_Y = np.stack([[[_d[(_d["sub"] == s) & (_d.eye == e) & (_d.intensity == lv)].ITPC_mean.iloc[0] for lv in L.LEV] for e in ("EO", "EC")] for s in sorted(_d["sub"].unique())])
_me = rm_anova_n(_Y, ["eye", "intensity"]).set_index("effect")
for a_, b_ in (("eye", "eye"), ("intensity", "intensity"), ("eye x intensity", "eye x intensity")):
    assert abs(_pg.loc[a_, "F"] - _me.loc[b_, "F"]) < 1e-8 and abs(_pg.loc[a_, "p"] - _me.loc[b_, "p"]) < 1e-6, (a_, _pg.loc[a_].to_dict(), _me.loc[b_].to_dict())
T = []
for imf in (6, 7):
    d = SUM[(SUM.IMF == imf) & SUM.in_16]
    Y = np.stack([[[[d[(d.participant == f"{s:03d}") & (d.stim == st) & (d.eye == e) & (d.intensity == lv)].ITPC_mean.iloc[0] for e in ("EO", "EC")] for lv in L.LEV]
                   for st in ("real", "sham")] for s in S16])            # n x stim x intensity x eye
    A3 = rm_anova_n(Y, ["stimulation", "intensity", "eye"])
    for r in A3.itertuples(): T.append(dict(test="3-way rmANOVA (16)", IMF=imf, effect_or_contrast=r.effect, F_or_t=r.F, df1=r.df1, df2=r.df2, p=r.p, p_holm=np.nan, dz=np.nan, eta_p2=r.eta_p2, mean_diff=np.nan))
    for st in ("real", "sham"):   # 2 x 3 per stimulation (16)
        A2 = rm_anova_n(Y[:, ("real", "sham").index(st)].transpose(0, 2, 1), ["eye", "intensity"])
        for r in A2.itertuples(): T.append(dict(test=f"2 x 3 rmANOVA {st} (16)", IMF=imf, effect_or_contrast=r.effect, F_or_t=r.F, df1=r.df1, df2=r.df2, p=r.p, p_holm=np.nan, dz=np.nan, eta_p2=r.eta_p2, mean_diff=np.nan))
    # (d) real vs sham per condition (16), Holm over 6
    tmp = []
    for e in ("EO", "EC"):
        for lv in L.LEV:
            r = L.paired(Y[:, 0, L.LEV.index(lv), ("EO", "EC").index(e)], Y[:, 1, L.LEV.index(lv), ("EO", "EC").index(e)])
            tmp.append(dict(test="real vs sham paired t (16)", IMF=imf, effect_or_contrast=f"{e} {lv}% RMT, real - sham", F_or_t=r["t"], df1=r["df"], df2=np.nan, p=r["p"], dz=r["dz"], eta_p2=np.nan, mean_diff=r["mean_diff"]))
    hp = L.holm([x["p"] for x in tmp]); [x.update(p_holm=h) for x, h in zip(tmp, hp)]; T += tmp
    # (e) sham intensity post hoc (eyes pooled, 16), Holm over 3; real for reference
    for st in ("sham", "real"):
        P = Y[:, ("real", "sham").index(st)].mean(2); tmp = []
        for a_, b_ in ((120, 100), (120, 80), (100, 80)):
            r = L.paired(P[:, L.LEV.index(a_)], P[:, L.LEV.index(b_)])
            tmp.append(dict(test=f"{st} intensity post hoc (16, eyes pooled)", IMF=imf, effect_or_contrast=f"{a_} vs {b_}", F_or_t=r["t"], df1=r["df"], df2=np.nan, p=r["p"], dz=r["dz"], eta_p2=np.nan, mean_diff=r["mean_diff"]))
        hp = L.holm([x["p"] for x in tmp]); [x.update(p_holm=h) for x, h in zip(tmp, hp)]; T += tmp
    # stimulation x intensity follow-up: (real - sham) intensity slope contrasts
    D_ = (Y[:, 0] - Y[:, 1]).mean(2); tmp = []
    for a_, b_ in ((120, 100), (120, 80), (100, 80)):
        r = L.paired(D_[:, L.LEV.index(a_)], D_[:, L.LEV.index(b_)])
        tmp.append(dict(test="(real - sham) intensity contrast (16)", IMF=imf, effect_or_contrast=f"{a_} vs {b_}", F_or_t=r["t"], df1=r["df"], df2=np.nan, p=r["p"], dz=r["dz"], eta_p2=np.nan, mean_diff=r["mean_diff"]))
    hp = L.holm([x["p"] for x in tmp]); [x.update(p_holm=h) for x, h in zip(tmp, hp)]; T += tmp
TD = pd.DataFrame(T); TD.to_csv(OUT / "N2_sham_ITPC_tests.csv", index=False)
DESC = SUM[SUM.in_16].groupby(["IMF", "stim", "intensity"]).ITPC_mean.agg(["mean", "std"]).reset_index(); DESC.to_csv(L.CK / "N2_summary_descriptives.csv", index=False)

# (f) figure (layout of Figure 5B, Fig5_v2 colours)
plt.rcParams.update({"font.family": "Times New Roman", "font.size": 9, "axes.linewidth": .8})
COL = {80: "#9BB7B4", 100: "#2E6E6A", 120: "#C5550A"}
fig = plt.figure(figsize=(7.2, 3.3))
for j, imf in enumerate((6, 7)):
    ax = fig.add_axes([0.07 + j * 0.49, 0.25, 0.41, 0.63])
    for stim, Y, ls in (("Real", YR16, "-"), ("Sham", YS16, "--")):
        Ym = Y.mean(1)
        for li, lv in enumerate(L.LEV):
            x = Ym[:, li, IMFS.index(imf)]; m = x.mean(0); se = x.std(0, ddof=1) / np.sqrt(len(x))
            ax.fill_between(tt, m - se, m + se, color=COL[lv], alpha=.18 if stim == "Real" else .12, lw=0)
            ax.plot(tt, m, color=COL[lv], lw=1.4, ls=ls, label=f"{stim}, {lv}% RMT")
    ax.axvline(0, color="k", lw=.6, ls=":"); ax.axhline(0, color="0.6", lw=.5); ax.set_xlim(-100, 500); ax.set_xlabel("Time (ms)")
    ax.set_ylabel("ITPC (baseline-corrected)" if j == 0 else ""); ax.set_title(f"IMF{imf} ({'alpha' if imf == 6 else 'theta'})", fontsize=9)
    if j == 0: ax.text(-0.16, 1.04, "", transform=ax.transAxes); hl = ax.get_legend_handles_labels()
order = [0, 3, 1, 4, 2, 5]
fig.legend([hl[0][i] for i in order], [hl[1][i] for i in order], loc="lower center", bbox_to_anchor=(0.5, -0.01), ncol=3, frameon=False, fontsize=7.5)
for ext in ("png", "pdf"): fig.savefig(OUT / f"FigureS_shamITPC_1110.{ext}", dpi=600, bbox_inches="tight")
plt.close(fig)
pd.set_option("display.width", 250)
print("sham intensity sig points per IMF:", SC["intensity"], "eye:", SC["eye"], "interaction:", SC["interaction"])
print(SEG[(SEG.data == "sham_16")].to_string()); print(JD.to_string())
print(TD.round(4).to_string()); print(DESC.round(3).to_string())

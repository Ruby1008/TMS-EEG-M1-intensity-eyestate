# -*- coding: utf-8 -*-
"""MEP reprocessing with one fixed rule set for every session.

Input  : the two CSV files written by mep_extract.m (raw EMG, -100 to +150 ms
         around the trigger edge, 1000 Hz, mV).
Output : trial table, cell means, participant x condition matrix, MEP
         probability, QC flags and plots, statistics, and a comparison with
         the old VETA values.

Rules (pre-specified; change only through the command-line options):
  pulse        onset of the EMG stimulus artifact in each trial (the trigger channel
               lags the pulse by a variable 5-12 ms); session median if no artifact
  pre-stim RMS -54 to -5 ms before the pulse (50 samples); excluded if > 0.05 mV
  MEP          peak-to-peak amplitude 15-50 ms after the pulse, measured on every
               trial (no zeroing, no change-point detection). The window starts at
               15 ms because the TMS artifact and its slow recovery in the EMG can
               reach about 13 ms at 120% RMT (sensitivity analysis: 10-60 ms).
               Before the measurement the artifact spike (-2 to +6 ms) is replaced
               by linear interpolation and a 2nd-order zero-phase 20 Hz high-pass
               removes the slow artifact recovery (sensitivity analysis: 60 Hz)
  MEP present  peak-to-peak >= 0.05 mV (reported as MEP probability)
  cell mean    mean over all accepted trials (present or not); cells with fewer
               than 10 accepted trials are set to missing
  statistics   ln(cell mean), 2 x 3 repeated-measures ANOVA with Greenhouse-Geisser
               correction, Holm-corrected paired t-tests of intensity (eye states
               averaged) with Cohen's dz, EO vs EC two one-sided tests (dz = +/-0.5)
               and the smallest supported equivalence bound

Two passes:
  1) run without --qc-log  -> qc_log_template.csv and qc_plots/*.png
  2) a rater inspects the plots and, for flagged trials only, enters a 'decision'
     and a 'reason' (an empty decision means accept):
       reject  technical failure (pre-activation, coil displacement, saturated
               signal); the trial is excluded. "No MEP" is never a reason.
       absent  the window holds only a smooth artifact recovery without an MEP;
               the trial is kept as "no MEP" at the session noise level
     The script is then rerun with --qc-log <file>.
"""
import argparse
import json
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

INTENSITY = {1: 80, 2: 100, 3: 120}
CELLS = [(e, i) for e in ("EO", "EC") for i in (80, 100, 120)]


# --------------------------------------------------------------------------- measures
def artifact_onsets(X, t, d, a):
    """TMS pulse time per trial (ms re trigger edge): first sample of the EMG stimulus artifact.

    The trigger channel lags the pulse by a variable 5-12 ms, so the pulse is taken
    from the artifact: the first sample within -25 to +5 ms of the trigger at which
    the absolute sample-to-sample change exceeds max(8 x its 95th percentile from
    -90 to -30 ms, 0.01 mV). Trials without a detectable artifact take the median
    onset of their session.
    """
    D = np.abs(np.diff(X, axis=1))
    td = t[1:]
    base = np.nanpercentile(D[:, (td >= -90) & (td <= -30)], 95, axis=1)
    sw = (td >= -25) & (td <= 5)
    on = np.full(len(X), np.nan)
    for i in range(len(X)):
        k = np.where(D[i, sw] > max(a.art_k * base[i], a.art_min))[0]
        if len(k):
            on[i] = td[sw][k[0]] - 1
    d["pulse_ms_re_trigger"] = on
    d["pulse_source"] = np.where(np.isnan(on), "session median", "EMG artifact")
    med = d.groupby(["subject", "eye"]).pulse_ms_re_trigger.transform("median")
    d["pulse_ms_re_trigger"] = d.pulse_ms_re_trigger.fillna(med).round().astype(int)
    return d


def clean_trace(x, t, onset, a):
    """blank the artifact spike (linear interpolation) and high-pass (zero phase) if requested"""
    if not a.hp_hz:
        return x
    from scipy.signal import butter, filtfilt
    y = x.copy()
    i0 = int(np.searchsorted(t, onset - 2))
    i1 = int(np.searchsorted(t, onset + a.blank_ms))
    y[i0:i1 + 1] = np.linspace(y[i0], y[i1], i1 - i0 + 1)
    b, c = butter(2, a.hp_hz / 500.0, "high")
    return filtfilt(b, c, y, padtype="odd", padlen=60)


def measure(meta, wave, a):
    t = np.array([int(c[1:]) for c in wave.columns if c != "row_id"])
    X = wave.drop(columns="row_id").to_numpy(dtype=float)
    d = artifact_onsets(X, t, meta.copy(), a)
    on = d.pulse_ms_re_trigger.to_numpy()
    Xc = np.vstack([clean_trace(X[i], t, on[i], a) for i in range(len(X))])

    def seg(M, lo, hi):
        """values from lo to hi ms after the pulse, per trial (NaN outside the epoch)"""
        out = np.full((len(M), hi - lo + 1), np.nan)
        for i in range(len(M)):
            idx = np.searchsorted(t, on[i] + np.arange(lo, hi + 1))
            ok = (idx >= 0) & (idx < len(t))
            out[i, ok] = M[i, idx[ok]]
        return out

    pre = seg(X, a.rms_start, a.rms_end)
    win = seg(Xc, a.win_start, a.win_end)
    sens = seg(Xc, a.sens_start, a.sens_end)
    tail = seg(Xc, a.win_start - 3, a.win_start)
    d["rms_pre"] = np.sqrt(np.nanmean(pre ** 2, axis=1))
    d["p2p"] = np.nanmax(win, axis=1) - np.nanmin(win, axis=1)
    d["p2p_sens_window"] = np.nanmax(sens, axis=1) - np.nanmin(sens, axis=1)
    tw = np.arange(a.win_start, a.win_end + 1)
    d["t_max"] = tw[np.nanargmax(win, axis=1)]
    d["t_min"] = tw[np.nanargmin(win, axis=1)]
    d["artifact_tail"] = np.nanmax(np.abs(tail), axis=1)
    d["present"] = d["p2p"] >= a.present
    return d, t, X, Xc


def qc_flags(d, a):
    f = pd.Series("", index=d.index)
    edge = d[["t_max", "t_min"]].isin([a.win_start, a.win_start + 1, a.win_end - 1, a.win_end]).any(axis=1)
    add = lambda mask, tag: f.where(~mask, f + tag + ";")
    f = add((d.rms_pre > a.rms_sens) & (d.rms_pre <= a.rms_max), "borderline_preactivation")
    f = add((d.p2p >= a.present * 0.6) & (d.p2p < a.present * 1.4), "borderline_MEP")
    f = add(edge & d.present, "extremum_at_window_edge")
    f = add(d.artifact_tail > a.present, "artifact_decay_at_window_start")
    f = add(d.pulse_source.ne("EMG artifact") & d.present, "pulse_time_from_session_median")
    return f.str.rstrip(";")


# ------------------------------------------------------------------------ statistics
def rm_anova_2x3(Y):
    """Y: n x 6 array, columns ordered as CELLS (eye x intensity). Returns DataFrame."""
    n = Y.shape[0]
    ce = np.array([[1], [-1]]) / np.sqrt(2)
    ci = np.array([[1, 0], [-1, 1], [0, -1]], dtype=float)
    ci = np.linalg.qr(ci)[0]
    me = np.ones((2, 1)) / np.sqrt(2)
    mi = np.ones((3, 1)) / np.sqrt(3)
    effects = {"eye state": np.kron(ce, mi), "intensity": np.kron(me, ci), "eye state x intensity": np.kron(ce, ci)}
    rows = []
    for name, C in effects.items():
        Z = Y @ C
        k = C.shape[1]
        zbar = Z.mean(0)
        ss_eff = n * np.sum(zbar ** 2)
        ss_err = np.sum((Z - zbar) ** 2)
        df1, df2 = k, (n - 1) * k
        F = (ss_eff / df1) / (ss_err / df2)
        S = np.cov(Z, rowvar=False).reshape(k, k)
        eps = 1.0 if k == 1 else np.trace(S) ** 2 / (k * np.trace(S @ S))
        p_gg = stats.f.sf(F, df1 * eps, df2 * eps)
        rows.append(dict(effect=name, F=F, df1=df1, df2=df2, eps_GG=eps, df1_GG=df1 * eps, df2_GG=df2 * eps,
                         p_GG=p_gg, eta_p2=ss_eff / (ss_eff + ss_err)))
    return pd.DataFrame(rows)


def holm(p):
    p = np.asarray(p)
    order = np.argsort(p)
    adj = np.empty_like(p)
    run = 0
    for r, j in enumerate(order):
        run = max(run, (len(p) - r) * p[j])
        adj[j] = min(1, run)
    return adj


def paired(x, y):
    dlt = x - y
    n = len(dlt)
    t_, p_ = stats.ttest_rel(x, y)
    return dict(n=n, mean_diff=dlt.mean(), t=t_, p=p_, dz=dlt.mean() / dlt.std(ddof=1))


def tost(x, y, bound_dz=0.5, alpha=0.05):
    dlt = x - y
    n = len(dlt)
    sd = dlt.std(ddof=1)
    se = sd / np.sqrt(n)
    b = bound_dz * sd
    p_lo = stats.t.sf((dlt.mean() + b) / se, n - 1)
    p_hi = stats.t.cdf((dlt.mean() - b) / se, n - 1)
    dz = dlt.mean() / sd
    smallest = abs(dz) + stats.t.ppf(1 - alpha, n - 1) / np.sqrt(n)
    return dict(n=n, dz=dz, p_TOST=max(p_lo, p_hi), equivalent=max(p_lo, p_hi) < alpha, smallest_bound_dz=smallest)


def run_stats(wide, label):
    cols = ["%s_%d" % c for c in CELLS]
    cc = wide.dropna(subset=cols)
    Y = np.log(cc[cols].to_numpy())
    out = {"label": label, "n_complete": len(cc), "subjects": list(cc.index)}
    out["anova"] = rm_anova_2x3(Y).to_dict("records")
    ln = pd.DataFrame(Y, columns=cols, index=cc.index)
    inten = {i: ln[["EO_%d" % i, "EC_%d" % i]].mean(1) for i in (80, 100, 120)}
    pairs = [(100, 80), (120, 80), (120, 100)]
    res = [dict(contrast="%d vs %d" % pr, **paired(inten[pr[0]].values, inten[pr[1]].values)) for pr in pairs]
    for r, ph in zip(res, holm([r["p"] for r in res])):
        r["p_Holm"] = ph
    out["intensity_posthoc"] = res
    out["eye_TOST"] = [dict(intensity=i, **tost(ln["EO_%d" % i].values, ln["EC_%d" % i].values)) for i in (80, 100, 120)]
    return out


# ------------------------------------------------------------------------------ main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("extract_dir", help="folder with mep_trials_meta.csv and mep_waveforms.csv")
    ap.add_argument("out_dir")
    ap.add_argument("--no-pulse", help="CSV with columns subject, eye, trial_number (trials without a TMS pulse)")
    ap.add_argument("--intensity-overrides", help="CSV with columns subject, eye, trial_number, intensity (e.g. 022 EO from the EEG codes)")
    ap.add_argument("--exclude-subjects", default="", help="comma-separated subject codes, e.g. 014,025,040")
    ap.add_argument("--qc-log", help="filled qc_log_template.csv (pass 2)")
    ap.add_argument("--stim", default="Real")
    ap.add_argument("--rms-start", type=int, default=-54)
    ap.add_argument("--rms-end", type=int, default=-5)
    ap.add_argument("--rms-max", type=float, default=0.05)
    ap.add_argument("--rms-sens", type=float, default=0.02, help="stricter RMS limit for the sensitivity analysis")
    ap.add_argument("--win-start", type=int, default=15)
    ap.add_argument("--win-end", type=int, default=50)
    ap.add_argument("--sens-start", type=int, default=10)
    ap.add_argument("--sens-end", type=int, default=60)
    ap.add_argument("--present", type=float, default=0.05)
    ap.add_argument("--art-k", type=float, default=8.0, help="artifact threshold: multiple of the baseline 95th percentile change")
    ap.add_argument("--art-min", type=float, default=0.01, help="artifact threshold floor (mV per sample)")
    ap.add_argument("--hp-hz", type=float, default=20.0, help="zero-phase 2nd-order high-pass after spike blanking (0 = off)")
    ap.add_argument("--sens-hp", type=float, default=60.0, help="high-pass for the sensitivity analysis (0 = skip)")
    ap.add_argument("--blank-ms", type=int, default=6, help="spike blanking from -2 to this many ms after the pulse")
    ap.add_argument("--min-trials", type=int, default=10)
    ap.add_argument("--no-plots", action="store_true")
    ap.add_argument("--plot-subjects", default="", help="comma-separated subject codes to plot (default: all)")
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)

    meta = pd.read_csv(os.path.join(a.extract_dir, "mep_trials_meta.csv"), dtype={"subject": str})
    wave = pd.read_csv(os.path.join(a.extract_dir, "mep_waveforms.csv"))
    assert (meta.row_id.values == wave.row_id.values).all()
    keep = meta.stim.eq(a.stim)
    meta, wave = meta[keep].reset_index(drop=True), wave[keep].reset_index(drop=True)
    dup = meta.groupby(["subject", "eye"]).file.nunique()
    if (dup > 1).any():
        sys.exit("several files for the same subject x eye (choose one batch):\n%s" % dup[dup > 1])

    d, t, X, Xc = measure(meta, wave, a)
    if a.sens_hp:
        b = argparse.Namespace(**vars(a))
        b.hp_hz = a.sens_hp
        d["p2p_sens_hp"] = measure(meta, wave, b)[0].p2p.values
    else:
        d["p2p_sens_hp"] = np.nan
    d["intensity"] = d.TMS_Condition.map(INTENSITY)
    if a.intensity_overrides:
        ov = pd.read_csv(a.intensity_overrides, dtype={"subject": str})
        d = d.merge(ov.rename(columns={"intensity": "intensity_override"}), on=["subject", "eye", "trial_number"], how="left")
        changed = d.intensity_override.notna() & (d.intensity_override != d.intensity)
        d.loc[d.intensity_override.notna(), "intensity"] = d.intensity_override
        print("intensity labels overridden:", int(changed.sum()))
        d = d.drop(columns="intensity_override")

    reason = pd.Series("", index=d.index)
    reason[d.pulse_sample.isna()] = "no trigger"
    if a.no_pulse:
        npul = pd.read_csv(a.no_pulse, dtype={"subject": str})
        key = d.set_index(["subject", "eye", "trial_number"]).index
        mask = key.isin(npul.set_index(["subject", "eye", "trial_number"]).index)
        reason[mask & reason.eq("")] = "no TMS pulse"
    reason[(d.rms_pre > a.rms_max) & reason.eq("")] = "pre-stimulus RMS > %.2f mV" % a.rms_max
    d["qc_flags"] = qc_flags(d, a)

    d["qc_decision"] = ""
    if a.qc_log:
        log = pd.read_csv(a.qc_log, dtype={"subject": str})
        dec = log.decision.fillna("").astype(str).str.strip().str.lower()
        if (~dec.isin(["", "accept", "reject", "absent"])).any():
            sys.exit("decision must be empty, accept, reject or absent")
        if log.loc[dec.ne("") & dec.ne("accept"), "reason"].isna().any():
            sys.exit("every reject/absent decision in the QC log needs a reason")
        log = log.assign(decision=dec).set_index(["subject", "eye", "trial_number"])
        key = pd.MultiIndex.from_frame(d[["subject", "eye", "trial_number"]])
        d["qc_decision"] = pd.Series(key.map(log.decision), index=d.index).fillna("")
        rej = d.qc_decision.eq("reject")
        reason[rej & reason.eq("")] = "visual QC: " + pd.Series(key.map(log.reason), index=d.index)[rej]
    excl_subj = [s.strip() for s in a.exclude_subjects.split(",") if s.strip()]
    reason[d.subject.isin(excl_subj) & reason.eq("")] = "subject excluded"
    d["exclusion"] = reason
    d["accepted"] = reason.eq("")
    # 'absent': the window holds only artifact recovery. The trial is kept, coded as
    # no MEP, and its amplitudes are set to the session noise level (median of the
    # accepted trials without an MEP in that session).
    ab = d.accepted & d.qc_decision.eq("absent")
    for col in ("p2p", "p2p_sens_window", "p2p_sens_hp"):
        if not ab.any() or d[col].isna().all():
            continue
        ok = d.accepted & ~ab & (d[col] < a.present)
        floor = d[ok].groupby(["subject", "eye"])[col].median()
        fill = d[ab].set_index(["subject", "eye"]).index.map(floor)
        d.loc[ab, col] = np.where(pd.isna(fill), d.loc[ok, col].median(), fill).astype(float)
    d.loc[ab, "present"] = False
    d["ln_p2p"] = np.log(d.p2p)

    # ---------------- cells
    acc = d[d.accepted]
    cells = (acc.groupby(["subject", "eye", "intensity"])
             .agg(n_accepted=("p2p", "size"), mean_p2p=("p2p", "mean"), mean_p2p_sens_window=("p2p_sens_window", "mean"),
                  mean_p2p_sens_hp=("p2p_sens_hp", "mean"),
                  mep_probability=("present", "mean"), mean_p2p_present_only=("p2p", lambda s: s[acc.loc[s.index, "present"]].mean()))
             .reset_index())
    rec = d[~d.subject.isin(excl_subj)].groupby(["subject", "eye", "intensity"]).size().rename("n_recorded").reset_index()
    cells = rec.merge(cells, how="left", on=["subject", "eye", "intensity"])
    cells["n_accepted"] = cells.n_accepted.fillna(0).astype(int)
    cells["below_min_trials"] = cells.n_accepted < a.min_trials

    def wide_of(col, rms_strict=False):
        src = cells
        if rms_strict:
            s = acc[acc.rms_pre <= a.rms_sens].groupby(["subject", "eye", "intensity"]).p2p.agg(["size", "mean"]).reset_index()
            s.loc[s["size"] < a.min_trials, "mean"] = np.nan
            s = s.rename(columns={"mean": col})
            src = s
        else:
            src = cells.copy()
            src.loc[src.below_min_trials, col] = np.nan
        w = src.pivot_table(index="subject", columns=["eye", "intensity"], values=col)
        w.columns = ["%s_%d" % c for c in w.columns]
        return w.reindex(columns=["%s_%d" % c for c in CELLS])

    W = wide_of("mean_p2p")
    W.to_csv(os.path.join(a.out_dir, "MEP_matrix_reprocessed.csv"))
    d.drop(columns=[]).to_csv(os.path.join(a.out_dir, "mep_trials_reprocessed.csv"), index=False)
    cells.to_csv(os.path.join(a.out_dir, "mep_cells_reprocessed.csv"), index=False)

    # ---------------- counts per condition
    cnt = d[~d.subject.isin(excl_subj)].assign(present_acc=lambda x: x.present & x.accepted).groupby(["eye", "intensity"]).agg(
        recorded=("row_id", "size"),
        excl_no_pulse=("exclusion", lambda s: s.str.startswith("no ").sum()),
        excl_rms=("exclusion", lambda s: s.str.startswith("pre-stimulus").sum()),
        excl_visual=("exclusion", lambda s: s.str.startswith("visual").sum()),
        accepted=("accepted", "sum"), mep_present=("present_acc", "sum")).reset_index()
    cnt["mep_probability_pct"] = 100 * cnt.mep_present / cnt.accepted
    cnt.to_csv(os.path.join(a.out_dir, "mep_counts_by_condition.csv"), index=False)

    # ---------------- statistics: primary and sensitivity
    results = [run_stats(W, "primary: window %d-%d ms, RMS <= %.2f mV" % (a.win_start, a.win_end, a.rms_max)),
               run_stats(wide_of("mean_p2p_sens_window"), "sensitivity: window %d-%d ms" % (a.sens_start, a.sens_end)),
               run_stats(wide_of("p2p", rms_strict=True), "sensitivity: RMS <= %.2f mV" % a.rms_sens),
               run_stats(wide_of("mean_p2p_sens_hp"), "sensitivity: high-pass %g Hz" % a.sens_hp),
               run_stats(wide_of("mean_p2p_present_only"), "sensitivity: trials with an MEP only")]
    Wp = wide_of("mep_probability")
    Wp.to_csv(os.path.join(a.out_dir, "MEP_probability_matrix.csv"))
    with open(os.path.join(a.out_dir, "mep_stats_reprocessed.json"), "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2, default=float)
    with open(os.path.join(a.out_dir, "mep_stats_reprocessed.txt"), "w", encoding="utf-8") as f:
        f.write("pass: %s\n\n" % ("2 (QC log applied)" if a.qc_log else "1 (no QC log, provisional)"))
        f.write(cnt.round(1).to_string(index=False) + "\n\n")
        for r in results:
            f.write("== %s (n = %d)\n" % (r["label"], r["n_complete"]))
            f.write(pd.DataFrame(r["anova"]).round(4).to_string(index=False) + "\n")
            f.write(pd.DataFrame(r["intensity_posthoc"]).round(4).to_string(index=False) + "\n")
            f.write(pd.DataFrame(r["eye_TOST"]).round(4).to_string(index=False) + "\n\n")

    # ---------------- comparison with the old VETA values (accepted in both)
    both = d.accepted & d.veta_trial_accept.eq(1)
    cmp_ = dict(n_trials=int(both.sum()),
                r_new_vs_old_amplitude=float(np.corrcoef(d.loc[both & (d.veta_amplitude > 0), "p2p"], d.loc[both & (d.veta_amplitude > 0), "veta_amplitude"])[0, 1]),
                old_zero_amplitude_now_present=int((both & d.veta_amplitude.eq(0) & d.present).sum()),
                old_excluded_now_accepted=int((d.accepted & d.veta_trial_accept.eq(0)).sum()),
                old_accepted_now_excluded=int((~d.accepted & d.veta_trial_accept.eq(1)).sum()))
    json.dump(cmp_, open(os.path.join(a.out_dir, "comparison_old_vs_new.json"), "w"), indent=2)

    # ---------------- QC template and plots
    fl = d[d.qc_flags.ne("") & d.exclusion.eq("")]
    tmpl = fl[["row_id", "file", "subject", "eye", "trial_number", "intensity", "pulse_ms_re_trigger", "pulse_source",
               "rms_pre", "p2p", "t_max", "t_min", "qc_flags"]].copy()
    tmpl["decision"] = ""
    tmpl["reason"] = ""
    tmpl["rater"] = ""
    tmpl.to_csv(os.path.join(a.out_dir, "qc_log_template.csv"), index=False)
    if not a.no_plots:
        plot_sessions(d, t, X, Xc, a)
    print(cnt.round(1).to_string(index=False))
    print("flagged for visual QC:", len(tmpl), "| outputs in", a.out_dir)


def plot_sessions(d, t, X, Xc, a):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    od = os.path.join(a.out_dir, "qc_plots")
    os.makedirs(od, exist_ok=True)
    only = [s.strip() for s in a.plot_subjects.split(",") if s.strip()]
    for (subj, eye), g in d.groupby(["subject", "eye"]):
        if only and subj not in only:
            continue
        n = len(g)
        ncol = 14
        nrow = int(np.ceil(n / ncol))
        fig, axs = plt.subplots(nrow, ncol, figsize=(ncol * 1.6, nrow * 1.25), sharex=True)
        for ax, (i, r) in zip(axs.flat, g.iterrows()):
            ax.axvspan(a.win_start, a.win_end, color="0.9")
            tr = t - r.pulse_ms_re_trigger
            show = (tr >= -20) & (tr <= 90)
            ax.plot(tr[show], X[i, show], lw=0.6, color="k")
            if a.hp_hz:
                ax.plot(tr[show], Xc[i, show], lw=0.6, color="tab:blue")
            col = "red" if r.qc_flags else ("0.6" if r.exclusion else "k")
            ax.set_title("#%d %d%% %.2f" % (r.trial_number, r.intensity, r.p2p), fontsize=6, color=col)
            ax.tick_params(labelsize=4)
        for ax in list(axs.flat)[n:]:
            ax.axis("off")
        fig.suptitle("%s %s  (time re pulse; red title: flagged for QC, gray: excluded; title: trial, intensity, p2p mV; "
                     "shaded: %d-%d ms; black raw%s)" % (subj, eye, a.win_start, a.win_end, ", blue cleaned" if a.hp_hz else ""), fontsize=9)
        fig.subplots_adjust(left=0.02, right=0.99, bottom=0.03, top=0.94, wspace=0.35, hspace=0.6)
        fig.savefig(os.path.join(od, "%s_%s.png" % (subj, eye)), dpi=110)
        plt.close(fig)


if __name__ == "__main__":
    main()

"""QC summary table for the Supplementary (README 6-3).

Reads Logs/preproc_summary.json (export_preproc_summary.m) and the data lists, writes
  QC/QC_summary_wide.csv   one row per dataset (incl. excluded datasets, with reason)
  QC/QC_summary_long.csv   one row per dataset x intensity
Flags (README 6-3): residual 11-25 ms > 20 uV, rho < 0.3 (Real only; Sham pulse size does not
scale with intensity), any intensity < 25 trials, round-2 ICA removed > 20 components.
Usage: python build_qc_table.py [--dry-run]
"""
import csv, json, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
TH = {"resid": 20, "rho": 0.3, "trials": 25, "r2": 20}
STAGE_D = {"exclude": 50, "sens": 20}      # stage D (user 2026-10-10)
LEVELS = [80, 100, 120]
DRY = "--dry-run" in sys.argv

D = json.load(open(ROOT / "Logs" / "preproc_summary.json", encoding="utf-8"))
datalist = []
for f in ("datalist_Real.csv", "datalist_Sham.csv"):
    if (ROOT / f).exists():
        datalist += list(csv.DictReader(open(ROOT / f, encoding="utf-8-sig")))


def count(ics, pred):
    return sum(1 for c in ics if c["removed"] and pred(c["reason"]))


wide, long = [], []
for x in sorted(D, key=lambda x: (x["sub"], x["cond"])):
    i1, i2 = x["ica1"], x["ica2"]
    qc = {q["intensity"]: q for q in x["qc"]}
    nPer = dict(zip(x["pair"]["levels"], x["trialsPerIntensityKept"]))
    resid = [qc[v]["residual11to25_uV"] if v in qc else None for v in LEVELS]
    rmax = max(r for r in resid if r is not None)
    stageD = "exclude" if rmax > STAGE_D["exclude"] else "sensitivity" if rmax > STAGE_D["sens"] else "include"
    r2n = sum(c["removed"] for c in i2)
    rho = x["pair"]["rho"]
    flags = []
    if rmax > TH["resid"]: flags.append(f"residual>{TH['resid']}uV")
    if "Real" in x["cond"] and (rho is None or rho < TH["rho"]): flags.append(f"rho<{TH['rho']}")
    if any(nPer.get(v, 0) < TH["trials"] for v in LEVELS): flags.append(f"trials<{TH['trials']}")
    if r2n > TH["r2"]: flags.append(f"R2>{TH['r2']}")
    ovr = list(x.get("ica1Overrides") or []) + list(x.get("ica2Overrides") or [])
    row = {
        "id": x["id"], "sub": x["sub"], "stim": x["cond"].split("_")[0], "eyes": x["cond"].split("_")[1],
        "status": "processed", "labelSource": x["pair"]["source"], "rho": None if rho is None else round(rho, 3),
        "nBadChan": len(x["badChannels"] or []), "badChannels": " ".join(x["badChannels"] or []),
        "trialsEpoched": x["trialsEpoched"], "trialsKept": x["trialsKept"],
        **{f"trials_{v}": nPer.get(v) for v in LEVELS},
        "R1_decay": count(i1, lambda r: "decay" in r), "R1_tmsMuscle_only": count(i1, lambda r: r == "tmsMuscle"),
        "R1_other": count(i1, lambda r: "decay" not in r and r != "tmsMuscle"), "R1_total": sum(c["removed"] for c in i1),
        "R1_varPct": round(sum(c["varPct"] for c in i1 if c["removed"]), 2),
        **{f"R2_{k}": count(i2, lambda r, k=k: r == k) for k in ("blink", "eyeMove", "muscle", "elecNoise")},
        "R2_other": count(i2, lambda r: r not in ("blink", "eyeMove", "muscle", "elecNoise")), "R2_total": r2n,
        "R2_varPct": round(sum(c["varPct"] for c in i2 if c["removed"]), 2),
        "nManualOverrides": len(ovr),
        **{f"preShift_{v}_uV": round(qc[v]["preShift_uV"], 3) if v in qc else None for v in LEVELS},
        **{f"preSlope_{v}_uVper10ms": round(qc[v]["preSlope_uVper10ms"], 4) if v in qc else None for v in LEVELS},
        **{f"residual11to25_{v}_uV": None if r is None else round(r, 2) for v, r in zip(LEVELS, resid)},
        "stageD": stageD, "flags": "; ".join(flags), "excludeReason": "",
    }
    wide.append(row)
    for v in LEVELS:
        if v in qc:
            q = qc[v]
            long.append({"id": x["id"], "sub": x["sub"], "stim": row["stim"], "eyes": row["eyes"], "intensity": v,
                         "nTrials": q["nTrials"], "preShift_uV": q["preShift_uV"], "preSlope_uVper10ms": q["preSlope_uVper10ms"],
                         "residual11to25_uV": q["residual11to25_uV"], "p30_uV": q["p30_uV"], "stageD": stageD})

done = {r["id"] for r in wide}
for d in datalist:
    if d["id"] not in done and str(d.get("exclude", "0")).strip() not in ("0", ""):
        s, c = d["id"].split("_")[1], d["id"].split("_", 2)[2]
        wide.append({"id": d["id"], "sub": int(s), "stim": c.split("_")[0], "eyes": c.split("_")[1], "status": "excluded",
                     "excludeReason": d.get("excludeReason") or d.get("note") or d.get("reason") or ""})
cols = list(next(r for r in wide if r["status"] == "processed").keys())
wide.sort(key=lambda r: (r["sub"], r["stim"], r["eyes"]))

P = [r for r in wide if r["status"] == "processed"]
print(f"processed {len(P)}, excluded {len(wide) - len(P)}; long rows {len(long)}")
print("stage D:", {k: sum(r['stageD'] == k for r in P) for k in ("include", "sensitivity", "exclude")})
print("flagged:", sum(bool(r["flags"]) for r in P))
for r in P:
    if r["flags"]: print(f"  {r['id']}: {r['flags']}")
print("excluded:", [(r["id"], r["excludeReason"]) for r in wide if r["status"] == "excluded"])
if DRY:
    print("--dry-run: nothing written"); print("columns:", cols)
    for r in P[:2]: print(r)
else:
    for name, rows, cs in (("QC_summary_wide.csv", wide, cols), ("QC_summary_long.csv", long, list(long[0].keys()))):
        with open(ROOT / "QC" / name, "w", newline="", encoding="utf-8-sig") as f:
            w = csv.DictWriter(f, fieldnames=cs, restval=""); w.writeheader(); w.writerows(rows)
        print("wrote", ROOT / "QC" / name)

    # Excel for the Supplementary: notes + wide + long, flagged rows highlighted
    from openpyxl import Workbook
    from openpyxl.styles import Font, PatternFill, Alignment
    from openpyxl.utils import get_column_letter
    fill = {"exclude": "F4CCCC", "sensitivity": "FCE5CD", "flag": "FFF2CC", "excluded": "D9D9D9"}
    wb = Workbook(); ws = wb.active; ws.title = "說明"
    notes = [
        ("TEP 預處理 QC 彙整表（README 6-3）", ""),
        ("來源", "Logs/preproc_summary.json（export_preproc_summary.m）；程式 code/build_qc_table.py"),
        ("QC_wide", "一個條件一列；status = excluded 的列為未處理的資料（excludeReason）"),
        ("QC_long", "一個條件 × 強度一列"),
        ("trials_*", "各強度保留的 trial 數"),
        ("R1_*", "第一輪 ICA 移除數：decay（規則 A ∪ C，含同時被 tmsMuscle 抓到）、tmsMuscle_only、other；R1_varPct = 移除的解釋變異量 %"),
        ("R2_*", "第二輪 ICA 移除數：blink、eyeMove、muscle、elecNoise；R2_varPct = 移除的解釋變異量 %"),
        ("nManualOverrides", "ica_overrides.csv 中套用到此條件的人工決定數"),
        ("preShift / preSlope", "ROI = FC3/FC1/C3/C1 平均（_F_ERM，30 Hz low-pass 之前）。偏移 = mean(−20～−6 ms) − mean(−100～−50 ms)（µV）；斜率 = −50～−6 ms 線性斜率（µV/10 ms）"),
        ("residual11to25", "11–25 ms ROI 平均波形的最大絕對值（µV）"),
        ("p30_uV（QC_long）", "15–50 ms ROI 平均振幅（µV）"),
        ("stageD", f"任一強度殘留 > {STAGE_D['exclude']} µV → exclude；{STAGE_D['sens']}–{STAGE_D['exclude']} µV → sensitivity；其他 include"),
        ("flags", f"殘留 > {TH['resid']} µV、ρ < {TH['rho']}（只看 Real）、任一強度 < {TH['trials']} trials、第二輪移除 > {TH['r2']} 個 component"),
        ("顏色", "紅 = stageD exclude；橘 = sensitivity；黃 = 其他警示；灰 = 未處理（排除）"),
    ]
    for r in notes: ws.append(list(r))
    ws["A1"].font = Font(bold=True, size=13); ws.column_dimensions["A"].width = 22; ws.column_dimensions["B"].width = 110
    for name, rows, cs in (("QC_wide", wide, cols), ("QC_long", long, list(long[0].keys()))):
        ws = wb.create_sheet(name); ws.append(cs)
        for c in ws[1]: c.font = Font(bold=True); c.alignment = Alignment(wrap_text=True, vertical="top")
        for r in rows:
            ws.append([r.get(c) for c in cs])
            key = ("excluded" if r.get("status") == "excluded" else r.get("stageD") if r.get("stageD") in ("exclude", "sensitivity")
                   else "flag" if r.get("flags") else None)
            if key:
                for c in ws[ws.max_row]: c.fill = PatternFill("solid", fgColor=fill[key])
        for j, c in enumerate(cs, 1):
            ws.column_dimensions[get_column_letter(j)].width = max(8, min(40, max(len(str(r.get(c) or "")) for r in rows) + 2, len(c) * 0.8 + 2))
        ws.freeze_panes = "B2"; ws.auto_filter.ref = ws.dimensions
    wb.save(ROOT / "QC" / "QC_summary.xlsx"); print("wrote", ROOT / "QC" / "QC_summary.xlsx")

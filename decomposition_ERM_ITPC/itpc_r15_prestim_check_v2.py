"""Round 15 merge, part 1 step 0, v2 rule (user 2026-10-22): baseline -500..-100 ms, IMF3-7 (IMF8 excluded), main test range 10..500 ms.
Verification: the same rmANOVA (GG; functions of itpc_round15.py) tested at IMF3-7 x every 1-ms point from -500 to 0 ms; BH-FDR q = .05 with one
family per effect over IMF3-7 x -500..0 ms. Required: 0 significant points for eye state, intensity and interaction.
Main sample (n = 17). Output: data/ITPC_r15_prestim_check_v2.csv."""
import sys
from pathlib import Path
import numpy as np, pandas as pd, scipy.io as sio
sys.argv = [sys.argv[0]]
ROOT = Path(__file__).resolve().parent.parent; D = ROOT / "data"
src = (Path(__file__).parent / "itpc_round15.py").read_text(encoding="utf-8")
ns = {"__file__": str(Path(__file__).parent / "itpc_round15.py")}; exec(src.split("OLD = analyse(")[0], ns)   # definitions only
rm2, bh, subs = ns["rm2"], ns["bh"], ns["subs"]
IMFS = [3, 4, 5, 6, 7]
M = sio.loadmat(D / "ITPC_roi_real.mat", squeeze_me=True); X = M["X"]; ids = [str(i) for i in M["ids"]]; t = np.round(M["t"].astype(float)); idx = {i: k for k, i in enumerate(ids)}
XN = X - X[:, :, (t >= -500) & (t <= -100), :].mean(2, keepdims=True)
w = (t >= -500) & (t <= 0); tt = t[w]
Y = np.stack([np.stack([XN[idx[f"spTMS_{s:03d}_Real_{e}"]][np.array(IMFS) - 1][:, w, :] for e in ("EO", "EC")], 0) for s in subs], 0)
R = rm2(np.moveaxis(Y, 4, 2)); rows = []
for k, v in R.items():
    q = bh(v["p"])
    for a, imf in enumerate(IMFS):
        sig = q[a] < .05
        rows.append(dict(effect=k, IMF=imf, n_points_tested=len(tt), n_sig_q05=int(sig.sum()), min_q=float(q[a].min()),
                         sig_range_ms=f"{tt[sig].min():.0f}..{tt[sig].max():.0f}" if sig.any() else ""))
T = pd.DataFrame(rows); T.to_csv(D / "ITPC_r15_prestim_check_v2.csv", index=False)
print(T.to_string()); print(T.groupby("effect").n_sig_q05.sum().to_string())

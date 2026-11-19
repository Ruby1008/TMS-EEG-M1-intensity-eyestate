"""E (user 2026-10-28): delivered and accepted EMG trials per participant x eye state x intensity, 23 participants, from the final EMG trial table
MEP_reprocess_2026-10-22/main_qc2/mep_trials_reprocessed.csv (one row per recorded trial, intensity after the label corrections, 'accepted' = trial kept
after the exclusions: no TMS pulse, pre-stimulus RMS > 0.05 mV, visual QC technical failure). main_qc2_008 has the identical table (checked).
Output: Figures_2026-11-08n/trial_counts_1108n.csv; summary printed and written to _check/trial_counts_summary.json. Reads only."""
data_dir = r"path/to/data"   # release: set this path
import json
from pathlib import Path
import pandas as pd
MF = Path(data_dir + r"\Revision_2026\Manuscript_Figures"); N = MF / "Figures_2026-11-08n"; M = MF / "MEP_reprocess_2026-10-22"
T = pd.read_csv(M / "main_qc2" / "mep_trials_reprocessed.csv", dtype={"subject": str}); T2 = pd.read_csv(M / "main_qc2_008" / "mep_trials_reprocessed.csv", dtype={"subject": str})
assert T.equals(T2) and set(T.stim) == {"Real"}
G = T.groupby(["subject", "eye", "intensity"]).agg(delivered=("accepted", "size"), accepted=("accepted", "sum")).reset_index().rename(columns={"subject": "participant"})
G["accepted_pct"] = (100 * G.accepted / G.delivered).round(1); G = G[["participant", "eye", "intensity", "delivered", "accepted", "accepted_pct"]]
G.to_csv(N / "trial_counts_1108n.csv", index=False)
S = dict(participants=int(G.participant.nunique()), cells=len(G), delivered_min=int(G.delivered.min()), delivered_max=int(G.delivered.max()),
         accepted_min=int(G.accepted.min()), accepted_max=int(G.accepted.max()), pct_min=float(G.accepted_pct.min()), pct_max=float(G.accepted_pct.max()),
         cells_delivered_over_42=G[G.delivered > 42].to_dict("records"), cells_delivered_under_42=G[G.delivered < 42].to_dict("records"),
         cells_accepted_min=G[G.accepted == G.accepted.min()].to_dict("records"), cells_accepted_max=G[G.accepted == G.accepted.max()].to_dict("records"),
         cells_pct_min=G[G.accepted_pct == G.accepted_pct.min()].to_dict("records"),
         per_condition=G.groupby(["eye", "intensity"]).agg(delivered_total=("delivered", "sum"), accepted_total=("accepted", "sum"), accepted_min=("accepted", "min"), accepted_max=("accepted", "max")).reset_index().to_dict("records"),
         exclusions=T.exclusion.value_counts().to_dict())
G42 = G[G.delivered == 42]; S["among_42_delivered_accepted_range"] = [int(G42.accepted.min()), int(G42.accepted.max())]
json.dump(S, open(N / "_check" / "trial_counts_summary.json", "w"), indent=1, default=str); print(json.dumps(S, indent=1, default=str))

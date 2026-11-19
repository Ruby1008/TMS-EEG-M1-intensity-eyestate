"""Part 4 A2 (user 2026-10-22): blinded visual-QC page for the flagged trials of the reprocessed MEP (main after A1).
Trials: main/mep_trials_reprocessed.csv, accepted == True and qc_flags containing artifact_decay_at_window_start.
Per trial: EMG from extract/mep_waveforms.csv, time re pulse (pulse_ms_re_trigger of main), -20 to 90 ms; black = raw (median of -20 to -3 ms removed),
blue = processed with mep_reprocess.clean_trace (20 Hz high-pass, blanking -2 to 6 ms, as main); gray band 15-50 ms. Vertical axis scaled per panel,
no values (blinding). Order shuffled with seed 20261029; qc_id = Q1, Q2, ... in shuffled order, skipping numbers that contain the digit 0.
Blinding: subject, eye state, intensity, amplitude and file name are not written anywhere in the HTML. The PNGs are stored with a custom 64-character
alphabet (no digits, no 'E', no '%', no '_') and decoded in the browser, so that the raw HTML text can be searched literally.
Outputs: qc/MEP_QC_2026-10-29.html, qc/qc_id_map.csv (qc_id -> row_id, subject, eye, trial_number, intensity; not embedded), qc/blinding_check.csv.
Frozen after creation (MD5 in HANDOFF section 38): do not regenerate."""
data_dir = r"path/to/data"   # release: set this path
import base64, io, json, re, sys
from argparse import Namespace
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from PIL import Image
R = Path(data_dir + r"\Revision_2026\Manuscript_Figures\MEP_reprocess_2026-10-22"); Q = R / "qc"
sys.path.insert(0, str(R.parent / "code" / "mep_reprocess")); import mep_reprocess as MR
HTML = Q / "MEP_QC_2026-10-29.html"; MAP = Q / "qc_id_map.csv"
assert not HTML.exists() and not MAP.exists(), "frozen files exist: do not regenerate"
plt.rcParams.update({"font.family": "Times New Roman", "font.size": 8})
T = pd.read_csv(R / "main" / "mep_trials_reprocessed.csv", dtype={"subject": str})
S = T[T.accepted & T.qc_flags.fillna("").str.contains("artifact_decay_at_window_start")].copy()
S = S.sample(frac=1, random_state=20261029).reset_index(drop=True); ZF = [n for n in range(1, 100000) if "0" not in str(n)][:len(S)]; S["qc_id"] = [f"Q{n}" for n in ZF]   # serial numbers without the digit 0 (cannot contain a subject code 0xx)
W = pd.read_csv(R / "extract" / "mep_waveforms.csv").set_index("row_id"); t = np.array([int(c[1:]) for c in W.columns])
A = Namespace(hp_hz=20.0, blank_ms=6)                       # = main (mep_reprocess.py defaults used for main)
STD = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/"
CUS = "ABCDFGHIJKLMNOPQRSTUVWXYZ" + "abcdefghijklmnopqrstuvwxyz" + "!#$()*+,-./:;"
assert len(CUS) == 64 and len(set(CUS)) == 64 and not re.search(r"[0-9E%_\"'\\<>&]", CUS)
TR = str.maketrans(STD, CUS)
def panel(x, on):
    xc = MR.clean_trace(x, t, on, A); tr = t - on; k = (tr >= -20) & (tr <= 90)
    raw = x - np.median(x[(tr >= -20) & (tr <= -3)]); pr = xc
    lim = 1.15 * max(np.abs(raw[(tr >= 10) & (tr <= 90)]).max(), np.abs(pr[(tr >= 5) & (tr <= 90)]).max(), 1e-6)
    fig, ax = plt.subplots(figsize=(3.6, 2.0), dpi=100)
    ax.axvspan(15, 50, color="0.88", lw=0); ax.axvline(0, color="0.6", lw=0.6, ls=":")
    ax.plot(tr[k], raw[k], color="black", lw=0.8); ax.plot(tr[k], pr[k], color="#1f5fbf", lw=1.0)
    ax.set_xlim(-20, 90); ax.set_ylim(-lim, lim); ax.set_yticks([]); ax.set_xticks([-20, 0, 15, 50, 90]); ax.tick_params(length=2, pad=1)
    for s in ("top", "right", "left"): ax.spines[s].set_visible(False)
    fig.subplots_adjust(left=0.02, right=0.98, top=0.97, bottom=0.14)
    b = io.BytesIO(); fig.savefig(b, format="png"); plt.close(fig)
    im = Image.open(io.BytesIO(b.getvalue())).convert("RGB").quantize(colors=16, method=Image.Quantize.MEDIANCUT)
    o = io.BytesIO(); im.save(o, format="PNG", optimize=True); return base64.b64encode(o.getvalue()).decode().translate(TR)
IMG = {}
for i, r in enumerate(S.itertuples()):
    IMG[r.qc_id] = panel(W.loc[r.row_id].values.astype(float), r.pulse_ms_re_trigger)
    if i % 200 == 0: print(i, len(S))
S[["qc_id", "row_id", "subject", "eye", "trial_number", "intensity"]].to_csv(MAP, index=False)
ids = list(S.qc_id)
HEAD = """<!DOCTYPE html><html lang="zh-Hant"><head><meta charset="utf-8"><title>MEP 目視 QC</title>
<style>
:root{--bg:#fafafa;--fg:#111;--card:#fff;--line:#ccc;--absent:#fff3c4;--reject:#ffd6d6}
body{margin:0;font-family:"Noto Sans TC","Microsoft JhengHei",sans-serif;background:var(--bg);color:var(--fg)}
#bar{position:sticky;top:0;z-index:5;background:#1d3a5f;color:#fff;padding:6px 12px;display:flex;flex-wrap:wrap;gap:12px;align-items:center}
#bar button{font-size:14px;padding:3px 10px}
#crit{margin:12px 16px;padding:10px 14px;background:#fff;border:1px solid var(--line);line-height:1.6}
#grid{display:grid;grid-template-columns:repeat(auto-fill,minmax(380px,1fr));gap:10px;margin:12px 16px}
.cell{background:var(--card);border:1px solid var(--line);padding:6px}
.cell.absent{background:var(--absent)}.cell.reject{background:var(--reject)}
.cell h3{margin:0 0 4px;font-size:15px}.cell img{display:block;max-width:360px;width:auto}
.cell label{margin-right:12px;font-size:14px}.cell textarea{box-sizing:border-box;display:block;width:360px;height:36px;margin-top:4px}
.need{outline:2px solid #c00}
</style></head><body>
<div id="bar"><span id="pg"></span><button id="prev">上一頁</button><button id="next">下一頁</button>
<span>跳到 QC 編號 <input id="jump" size="7" placeholder="Q1"><button id="go">前往</button></span>
<span id="cnt"></span><button id="exp">匯出 CSV</button></div>
<div id="crit"><b>判準</b>（每格一個 trial；時間以 TMS 脈衝為 0，−20～90 ms；黑線＝原始 EMG，藍線＝處理後（內插 −2～6 ms＋20 Hz 高通）；灰底＝測量窗 15–50 ms；縱軸每格自動縮放、不顯示數值）<br>
<b>accept</b>：灰底內有 MEP 的快速多相波形，或看不出問題。<br>
<b>absent</b>：灰底內只有平滑的偽跡回復、沒有 MEP 波形（trial 保留，記為無 MEP）。<br>
<b>reject</b>：技術失敗（預先收縮、訊號飽和、波形斷裂），trial 排除。<br>
<b>不要以「看起來沒有 MEP」為由選 reject。</b> absent 與 reject 必須填理由。進度自動存在這台電腦的瀏覽器（localStorage），關掉重開不會遺失；完成後按「匯出 CSV」。</div>
<div id="grid"></div>
<script>
"""
JS = r"""
const CUS="ABCDFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz!#$()*+,-./:;", STD=(function(){let u="",l="",d="";for(let i=0;i<26;i++){u+=String.fromCharCode(65+i);l+=String.fromCharCode(97+i)}for(let i=0;i<10;i++)d+=String.fromCharCode(48+i);return u+l+d+"+/"})();
const MAPC={}; for(let i=0;i<CUS.length;i++) MAPC[CUS[i]]=STD[i];
function dec(s){let o="";for(const ch of s) o+=(MAPC[ch]!==undefined?MAPC[ch]:ch);return "data:image/png;base64,"+o;}
const KEY="mepqc_v1_store"; let st={}; try{st=JSON.parse(localStorage.getItem(KEY)||"{}")}catch(e){st={}}
function save(){try{localStorage.setItem(KEY,JSON.stringify(st))}catch(e){alert("無法寫入 localStorage："+e)}}
const PER=Math.floor(PAGE_N), NP=Math.ceil(IDS.length/PER); let page=0;
function get(id){return st[id]||{d:"accept",r:""}}
function counts(){let a=0,j=0;for(const id of IDS){const d=get(id).d;if(d==="absent")a++;if(d==="reject")j++}
 document.getElementById("cnt").textContent="absent "+a+" ／ reject "+j+" ／ 共 "+IDS.length+" 格";}
function render(){const g=document.getElementById("grid");g.innerHTML="";
 const s=page*PER, e=Math.min(IDS.length,s+PER);
 for(let i=s;i<e;i++){const id=IDS[i],v=get(id);const c=document.createElement("div");c.className="cell "+(v.d!=="accept"?v.d:"");c.id="c_"+id;
  c.innerHTML="<h3>"+id+"</h3><img alt='"+id+"'><div>"+["accept","absent","reject"].map(x=>"<label><input type='radio' name='r_"+id+"' value='"+x+"'"+(v.d===x?" checked":"")+"> "+x+"</label>").join("")+
   "</div><textarea placeholder='理由（absent／reject 必填）'></textarea>";
  c.querySelector("img").src=dec(IMG[id]); const ta=c.querySelector("textarea"); ta.value=v.r;
  c.querySelectorAll("input").forEach(r=>r.addEventListener("change",()=>{const o=get(id);o.d=r.value;st[id]=o;save();c.className="cell "+(o.d!=="accept"?o.d:"");counts();}));
  ta.addEventListener("input",()=>{const o=get(id);o.r=ta.value;st[id]=o;save();});
  g.appendChild(c);}
 document.getElementById("pg").textContent="第 "+(page+1)+" ／ "+NP+" 頁";counts();window.scrollTo(0,0);}
document.getElementById("prev").onclick=()=>{if(page>0){page--;render()}};
document.getElementById("next").onclick=()=>{if(page<NP-1){page++;render()}};
document.getElementById("go").onclick=()=>{const id=document.getElementById("jump").value.trim().toUpperCase();const k=IDS.indexOf(id);
 if(k<0){alert("找不到 "+id);return} page=Math.floor(k/PER);render();const el=document.getElementById("c_"+id);if(el){el.scrollIntoView();el.classList.add("need");setTimeout(()=>el.classList.remove("need"),1500)}};
function q(s){s=String(s).replace(/\r?\n/g," ");return /[",]/.test(s)?'"'+s.replace(/"/g,'""')+'"':s}
document.getElementById("exp").onclick=()=>{const miss=IDS.filter(id=>{const v=get(id);return v.d!=="accept"&&!v.r.trim()});
 if(miss.length){alert("以下 "+miss.length+" 格選了 absent／reject 但沒有理由，未匯出：\n"+miss.slice(0,40).join(", ")+(miss.length>40?" …":""));return}
 const rows=["qc_id,decision,reason"].concat(IDS.map(id=>{const v=get(id);return [id,v.d,q(v.r.trim())].join(",")}));
 const b=new Blob(["\ufeff"+rows.join("\r\n")+"\r\n"],{type:"text/csv"});const a=document.createElement("a");a.href=URL.createObjectURL(b);a.download="MEP_QC_decisions.csv";document.body.appendChild(a);a.click();a.remove();};
render();
</script></body></html>
"""
body = HEAD + "const IDS=" + json.dumps(ids) + ";\nconst PAGE_N=" + "100" + ";\nconst IMG={" + ",\n".join(f'"{k}":"{v}"' for k, v in IMG.items()) + "};\n" + JS
HTML.write_text(body, encoding="utf-8")
# ---- blinding check on the raw HTML text
txt = HTML.read_text(encoding="utf-8"); subs = sorted(T.subject.unique())
rows = []
for s in subs: rows.append(dict(term=f"subject code {s} (standalone 3-digit token)", hits=len(re.findall(rf"(?<![0-9]){s}(?![0-9])", txt))))
for s in subs: rows.append(dict(term=f"subject code {s} (any position)", hits=txt.count(s)))
for term in ("EO", "EC", "80%", "100%", "120%", "TMS_MEP", "spTMS", "Veta", "VETA", "preprocessed", ".mat", "subject", "intensity", "RMT"): rows.append(dict(term=term, hits=txt.count(term)))
for col in ("p2p", "rms_pre"):
    vals = {f"{v:.4f}" for v in S[col]}; rows.append(dict(term=f"any {col} value of the QC trials (4 decimals)", hits=sum(txt.count(v) for v in vals)))
B = pd.DataFrame(rows); B.to_csv(Q / "blinding_check.csv", index=False)
print("trials", len(S), "pages", -(-len(S) // 100), "HTML MB", round(HTML.stat().st_size / 1e6, 1)); print(B[B.hits > 0].to_string(index=False) if (B.hits > 0).any() else "blinding check: 0 hits for every term")

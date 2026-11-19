"""Part 5 E3 (user 2026-10-23): blinded visual-QC page v2 for the same 1,111 trials as v1 (qc/qc_id_map.csv, same qc_id), new order (seed 20261031).
Drawing as the diagnostic page (qc/diag/make_mep_diagnostic.py): per trial two stacked panels, top = raw EMG (black, minus the median of -20..-3 ms),
bottom = processed EMG (blue, mep_reprocess.clean_trace as main), -20..90 ms re pulse (pulse time of main), gray band 15-50 ms.
Vertical axis fixed per session: +-1.2 x 95th percentile of |processed EMG| at 15-50 ms over all accepted trials of the session (main); the scale is
written in the top-left corner of the image ("+-0.40 mV"). Reference: light-gray session mean of the raw EMG of the unflagged accepted trials with an MEP,
overlaid in the raw panel (fewer than 5 -> other eye state of the participant; fewer than 5 there too -> "no reference" in the image).
No separate reference cell (would reveal the session grouping). Blinding as v1: only the qc_id is shown; images stored with the custom alphabet of v1.
Decisions: accept (default), absent, reject, unsure; absent / reject / unsure need a reason. Export qc_id, decision, reason.
Outputs: qc/v2/MEP_QC_v2_2026-10-31.html, qc/v2/qc_v2_order.csv (qc_id, page position), qc/v2/blinding_check_v2.csv. Frozen after creation."""
import base64, io, json, re, sys
from argparse import Namespace
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from PIL import Image
V = Path(__file__).parent; Q = V.parent; R = Q.parent
sys.path.insert(0, str(R.parent / "code" / "mep_reprocess")); import mep_reprocess as MR
HTML = V / "MEP_QC_v2_2026-10-31.html"
assert not HTML.exists(), "frozen file exists: do not regenerate"
plt.rcParams.update({"font.family": "Times New Roman", "font.size": 8})
T = pd.read_csv(R / "main" / "mep_trials_reprocessed.csv", dtype={"subject": str})
M = pd.read_csv(Q / "qc_id_map.csv", dtype={"subject": str}); assert len(M) == 1111
S = M.merge(T[["row_id", "pulse_ms_re_trigger"]], on="row_id", how="left"); assert S.pulse_ms_re_trigger.notna().all()
S = S.sample(frac=1, random_state=20261031).reset_index(drop=True)
W = pd.read_csv(R / "extract" / "mep_waveforms.csv").set_index("row_id"); t = np.array([int(c[1:]) for c in W.columns]); A = Namespace(hp_hz=20.0, blank_ms=6)
STD = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/"
CUS = "ABCDFGHIJKLMNOPQRSTUVWXYZ" + "abcdefghijklmnopqrstuvwxyz" + "!#$()*+,-./:;"
assert len(CUS) == 64 and len(set(CUS)) == 64 and not re.search(r"[0-9E%_\"'\\<>&]", CUS)
TR = str.maketrans(STD, CUS)
def sig(row_id, on):
    x = W.loc[row_id].values.astype(float); xc = MR.clean_trace(x, t, on, A); tr = t - on; raw = x - np.median(x[(tr >= -20) & (tr <= -3)]); return tr, raw, xc
GRID = np.arange(-20, 91)
SCALE, REF = {}, {}
for (s, e), g in T[T.accepted].groupby(["subject", "eye"]):
    v = [np.abs(xc[(tr >= 15) & (tr <= 50)]) for tr, raw, xc in (sig(r.row_id, r.pulse_ms_re_trigger) for r in g.itertuples())]
    SCALE[(s, e)] = 1.2 * float(np.percentile(np.concatenate(v), 95))
def ref(s, e):
    for ee in (e, "EC" if e == "EO" else "EO"):
        g = T[(T.subject == s) & (T.eye == ee) & T.accepted & T.qc_flags.isna() & T.present]
        if len(g) >= 5:
            return np.mean([np.interp(GRID, *sig(r.row_id, r.pulse_ms_re_trigger)[:2]) for r in g.itertuples()], axis=0)
    return None
for k in SCALE: REF[k] = ref(*k)
def panel(r):
    tr, raw, xc = sig(r.row_id, r.pulse_ms_re_trigger); lim = SCALE[(r.subject, r.eye)]; rf = REF[(r.subject, r.eye)]; k = (tr >= -20) & (tr <= 90)
    fig, ax = plt.subplots(2, 1, figsize=(3.6, 2.9), dpi=100, sharex=True, gridspec_kw={"hspace": .06})
    for a in ax: a.axvspan(15, 50, color="0.88", lw=0); a.axvline(0, color="0.6", lw=.6, ls=":"); a.axhline(0, color="0.85", lw=.4); a.set_ylim(-lim, lim); a.set_yticks([]); a.set_xlim(-20, 90)
    if rf is not None: ax[0].plot(GRID, np.clip(rf, -lim * 1.02, lim * 1.02), color="0.72", lw=2.2)
    else: ax[0].text(89, -lim * .85, "no reference", ha="right", fontsize=7, color="0.4")
    ax[0].plot(tr[k], np.clip(raw[k], -lim * 1.02, lim * 1.02), color="black", lw=.8)
    ax[1].plot(tr[k], np.clip(xc[k], -lim * 1.02, lim * 1.02), color="#1f5fbf", lw=1.0)
    ax[0].text(-19, lim * .78, f"±{lim:.2f} mV", fontsize=7.5); ax[0].text(89, lim * .78, "raw", ha="right", fontsize=7); ax[1].text(89, lim * .78, "processed", ha="right", fontsize=7, color="#1f5fbf")
    ax[1].set_xticks([-20, 0, 15, 50, 90]); ax[1].tick_params(length=2, pad=1)
    fig.subplots_adjust(left=0.02, right=0.98, top=0.98, bottom=0.09)
    b = io.BytesIO(); fig.savefig(b, format="png"); plt.close(fig)
    im = Image.open(io.BytesIO(b.getvalue())).convert("RGB").quantize(colors=16, method=Image.Quantize.MEDIANCUT)
    o = io.BytesIO(); im.save(o, format="PNG", optimize=True); return base64.b64encode(o.getvalue()).decode().translate(TR)
IMG = {}
for i, r in enumerate(S.itertuples()):
    IMG[r.qc_id] = panel(r)
    if i % 200 == 0: print(i, len(S))
S.assign(position=np.arange(1, len(S) + 1))[["qc_id", "position"]].to_csv(V / "qc_v2_order.csv", index=False)
ids = list(S.qc_id)
HEAD = """<!DOCTYPE html><html lang="zh-Hant"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>MEP 目視 QC v2</title>
<style>
:root{--bg:#fafafa;--fg:#111;--card:#fff;--line:#ccc;--absent:#fff3c4;--reject:#ffd6d6;--unsure:#e3ecff}
body{margin:0;font-family:"Noto Sans TC","Microsoft JhengHei",sans-serif;background:var(--bg);color:var(--fg)}
#bar{position:sticky;top:0;z-index:5;background:#1d3a5f;color:#fff;padding:6px 12px;display:flex;flex-wrap:wrap;gap:12px;align-items:center}
#bar button{font-size:14px;padding:3px 10px}
#crit{margin:12px 16px;padding:10px 14px;background:#fff;border:1px solid var(--line);line-height:1.6}
#grid{display:grid;grid-template-columns:repeat(auto-fill,minmax(380px,1fr));gap:10px;margin:12px 16px}
.cell{background:var(--card);border:1px solid var(--line);padding:6px}
.cell.absent{background:var(--absent)}.cell.reject{background:var(--reject)}.cell.unsure{background:var(--unsure)}
.cell h3{margin:0 0 4px;font-size:15px}.cell img{display:block;max-width:360px;width:auto}
.cell label{margin-right:10px;font-size:14px}.cell textarea{box-sizing:border-box;display:block;width:auto;max-width:360px;height:36px;margin-top:4px}
.need{outline:2px solid #c00}
</style></head><body>
<div id="bar"><span id="pg"></span><button id="prev">上一頁</button><button id="next">下一頁</button>
<span>跳到 QC 編號 <input id="jump" size="7" placeholder="Q1"><button id="go">前往</button></span>
<span id="cnt"></span><button id="exp">匯出 CSV</button></div>
<div id="crit"><b>畫法</b>：每格一個 trial，上下兩小圖。上＝原始 EMG（黑），下＝處理後（藍，內插 −2～6 ms＋高通濾波）；時間以 TMS 脈衝為 0，−20～90 ms；灰底＝測量窗 15–50 ms。
縱軸是固定尺度（同一個 session 的所有格子相同，不自動縮放），尺度寫在每格左上角（例如「±0.40 mV」）；超出尺度的部分被截在邊緣。
上圖的<b>淡灰色粗線</b>是同一個 session 中無旗標且有 MEP 的 trial 的平均（參考 MEP）；寫「no reference」的格子沒有參考波形。<br>
<b>對照淡灰色參考波形：MEP 是快速、多相的尖波；偽跡回復是單一平滑弧線，通常從視窗起點就在走。</b><br>
<b>accept</b>（預設）：灰底內有 MEP 波形，或看不出問題。<br>
<b>absent</b>：灰底內只有平滑的偽跡回復、沒有 MEP 波形（trial 保留，記為無 MEP）。<br>
<b>reject</b>：技術失敗（預先收縮、訊號飽和、波形斷裂），trial 排除。<b>不要以「看起來沒有 MEP」為由選 reject。</b><br>
<b>unsure</b>：無法判斷。<br>
absent、reject、unsure 必須填理由。進度自動存在這台電腦的瀏覽器（localStorage）；完成後按「匯出 CSV」，缺理由時不會匯出。</div>
<div id="grid"></div>
<script>
"""
JS = r"""
const CUS="ABCDFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz!#$()*+,-./:;", STD=(function(){let u="",l="",d="";for(let i=0;i<26;i++){u+=String.fromCharCode(65+i);l+=String.fromCharCode(97+i)}for(let i=0;i<10;i++)d+=String.fromCharCode(48+i);return u+l+d+"+/"})();
const MAPC={}; for(let i=0;i<CUS.length;i++) MAPC[CUS[i]]=STD[i];
function dec(s){let o="";for(const ch of s) o+=(MAPC[ch]!==undefined?MAPC[ch]:ch);return "data:image/png;base64,"+o;}
const KEY="mepqc_v2_store"; let st={}; try{st=JSON.parse(localStorage.getItem(KEY)||"{}")}catch(e){st={}}
function save(){try{localStorage.setItem(KEY,JSON.stringify(st))}catch(e){alert("無法寫入 localStorage："+e)}}
const PER=Math.floor(PAGE_N), NP=Math.ceil(IDS.length/PER); let page=0; const OPTS=["accept","absent","reject","unsure"];
function get(id){return st[id]||{d:"accept",r:""}}
function counts(){const c={accept:0,absent:0,reject:0,unsure:0};for(const id of IDS)c[get(id).d]++;
 document.getElementById("cnt").textContent="accept "+c.accept+" ／ absent "+c.absent+" ／ reject "+c.reject+" ／ unsure "+c.unsure+" ／ 共 "+IDS.length+" 格";}
function render(){const g=document.getElementById("grid");g.innerHTML="";
 const s=page*PER, e=Math.min(IDS.length,s+PER);
 for(let i=s;i<e;i++){const id=IDS[i],v=get(id);const c=document.createElement("div");c.className="cell "+(v.d!=="accept"?v.d:"");c.id="c_"+id;
  c.innerHTML="<h3>"+id+"</h3><img alt='"+id+"'><div>"+OPTS.map(x=>"<label><input type='radio' name='r_"+id+"' value='"+x+"'"+(v.d===x?" checked":"")+"> "+x+"</label>").join("")+
   "</div><textarea placeholder='理由（absent／reject／unsure 必填）'></textarea>";
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
function missing(){return IDS.filter(id=>{const v=get(id);return v.d!=="accept"&&!v.r.trim()})}
function buildCSV(){return ["qc_id,decision,reason"].concat(IDS.map(id=>{const v=get(id);return [id,v.d,q(v.r.trim())].join(",")})).join("\r\n")+"\r\n"}
document.getElementById("exp").onclick=()=>{const miss=missing();
 if(miss.length){alert("以下 "+miss.length+" 格選了 absent／reject／unsure 但沒有理由，未匯出：\n"+miss.slice(0,40).join(", ")+(miss.length>40?" …":""));return}
 const b=new Blob(["﻿"+buildCSV()],{type:"text/csv"});const a=document.createElement("a");a.href=URL.createObjectURL(b);a.download="MEP_QC_v2_decisions.csv";document.body.appendChild(a);a.click();a.remove();};
render();
</script></body></html>
"""
body = HEAD + "const IDS=" + json.dumps(ids) + ";\nconst PAGE_N=" + "100" + ";\nconst IMG={" + ",\n".join(f'"{k}":"{v}"' for k, v in IMG.items()) + "};\n" + JS
HTML.write_text(body, encoding="utf-8")
# ---- blinding check on the raw HTML text (same list as v1)
txt = HTML.read_text(encoding="utf-8"); subs = sorted(T.subject.unique()); rows = []
for s in subs: rows.append(dict(term=f"subject code {s} (standalone 3-digit token)", hits=len(re.findall(rf"(?<![0-9]){s}(?![0-9])", txt))))
for s in subs: rows.append(dict(term=f"subject code {s} (any position)", hits=txt.count(s)))
for term in ("EO", "EC", "80%", "100%", "120%", "TMS_MEP", "spTMS", "Veta", "VETA", "preprocessed", ".mat", "subject", "intensity", "RMT"): rows.append(dict(term=term, hits=txt.count(term)))
T2 = T.set_index("row_id").loc[S.row_id]
for col in ("p2p", "rms_pre"):
    vals = {f"{v:.4f}" for v in T2[col]}; rows.append(dict(term=f"any {col} value of the QC trials (4 decimals)", hits=sum(txt.count(v) for v in vals)))
vals = {f"{v:.2f}" for v in SCALE.values()}; rows.append(dict(term="any session scale value (2 decimals, outside the images)", hits=sum(len(re.findall(rf"(?<![0-9.]){re.escape(v)}(?![0-9])", txt)) for v in vals)))
B = pd.DataFrame(rows); B.to_csv(V / "blinding_check_v2.csv", index=False)
nref = sum(v is None for v in REF.values())
print("trials", len(S), "pages", -(-len(S) // 100), "HTML MB", round(HTML.stat().st_size / 1e6, 1), "| sessions", len(SCALE), "without reference", nref)
print(B[B.hits > 0].to_string(index=False) if (B.hits > 0).any() else "blinding check: 0 hits for every term")

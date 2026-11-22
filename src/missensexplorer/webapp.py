"""MissensExplorer web interface (Streamlit).

Run with:  missensexplorer web     (or)     streamlit run streamlit_app.py
"""

from __future__ import annotations

import html
import json
import numbers
from urllib.parse import quote

import altair as alt
import httpx
import pandas as pd
import streamlit as st
import streamlit.components.v1 as components

from missensexplorer import __version__
from missensexplorer.aa import AA_NAMES
from missensexplorer.cache import Cache
from missensexplorer.inputs import vcf_samples
from missensexplorer.predict import PredictOptions, is_damaging, predict
from missensexplorer.report import to_bytes, to_fasta

st.set_page_config(page_title="MissensExplorer", page_icon="🧬", layout="wide")

EXAMPLE = """BRAF V600E
TP53: R175H, R248Q, R273H, P72R
KRAS:p.G12D
BRCA1:c.5096G>A
chr17:7675088 C>T
rs113488022
CFTR p.Phe508del"""

FORMATS_HELP = """
| You write | Meaning |
|---|---|
| `BRAF V600E`, `BRAF:p.Val600Glu` | protein change in a gene |
| `TP53: R175H, R248Q` | several changes in one gene |
| `BRCA1:c.5096G>A`, `NM_004333.6:c.1799T>A` | cDNA (HGVS c.) |
| `chr7:140753336 A>T`, `7-140753336-A-T` | genomic position + alleles |
| `NC_000007.14:g.140753336A>T` | genomic HGVS |
| `rs113488022` | dbSNP id (all alleles) |
| `V600E` | bare change, when *Default gene* is set |
"""

# Status colours (fixed palette); always paired with a shape, a label and the y-position.
STATUS = {
    "Damaging": ("#d03b3b", "triangle-up", "🔴"),
    "Uncertain": ("#fab219", "diamond", "🟡"),
    "Benign / silent": ("#0ca30c", "circle", "🟢"),
    "Not assessed": ("#8a8986", "square", "⚪"),
}


def status_of(label: str | None) -> str:
    label = label or ""
    if is_damaging(label):
        return "Damaging"
    low = label.lower()
    if "benign" in low or label.startswith("Silent"):
        return "Benign / silent"
    if label.startswith(("Uncertain", "Possible", "Protein-altering", "Splice region")):
        return "Uncertain"
    return "Not assessed"


def _isnum(v) -> bool:
    return isinstance(v, numbers.Real) and not isinstance(v, bool) and v == v


def _val(v, default=None):
    if v is None or (isinstance(v, float) and v != v) or v == "":
        return default
    return v


# --------------------------------------------------------------------------- sidebar

with st.sidebar:
    st.header("Settings")
    assembly = st.selectbox("Genome build", ["auto", "GRCh38", "GRCh37"],
                            help="Used for genomic coordinates and VCF files. 'auto' reads the VCF header "
                                 "and checks REF alleles against both builds.")
    predictors = st.toggle("Pathogenicity predictors", True,
                           help="AlphaMissense, REVEL, CADD, EVE/popEVE, SpliceAI, SIFT and PolyPhen-2 "
                                "via the Ensembl VEP REST plugins.")
    context = st.toggle("Protein context (UniProt + AlphaFold)", True)
    clinvar_on = st.toggle("ClinVar classification", True)
    noncoding = st.toggle("Include non-coding variants", False,
                          help="Also list intronic, UTR, upstream/downstream and intergenic consequences.")
    workers = st.slider("Parallel requests", 1, 8, 4)
    if st.button("Clear cache", help="Cached API answers make repeated analyses instant (30 days)."):
        Cache().clear()
        st.toast("Cache cleared")
    st.divider()
    st.caption(
        f"MissensExplorer {__version__} · data: Ensembl VEP, AlphaMissense, REVEL, CADD, EVE, SpliceAI, "
        "UniProt, AlphaFold DB, ClinVar, gnomAD.\n\n"
        "**Research use only** — computational predictions are not a clinical diagnosis."
    )

# --------------------------------------------------------------------------- input

st.title("🧬 MissensExplorer")
st.markdown("Type mutations or upload VCF / CSV files and see **how each variant changes the protein**.")

# Shareable links: ?q=BRAF+V600E;TP53+R175H&gene=TP53 fills the box and runs once.
if "q" in st.query_params and not st.session_state.get("_q_done"):
    st.session_state["text"] = st.query_params["q"].replace(";", "\n")
    st.session_state["_auto_run"] = True
    st.session_state["_q_done"] = True

tab_text, tab_files = st.tabs(["✍️  Type mutations", "📁  Upload files"])
with tab_text:
    c1, c2 = st.columns([3, 1])
    with c2:
        gene = st.text_input("Default gene (optional)", placeholder="e.g. TP53", value=st.query_params.get("gene", ""),
                             help="Lets you type bare changes such as 'R175H' or 'c.524G>A'.")
        st.button("Load example", on_click=lambda: st.session_state.update(text=EXAMPLE), width="stretch")
        with st.popover("Accepted formats", width="stretch"):
            st.markdown(FORMATS_HELP)
    with c1:
        text = st.text_area("Variants — one or more per line", key="text", height=190, placeholder=EXAMPLE)

with tab_files:
    uploads = st.file_uploader("VCF, VCF.gz, CSV, TSV, XLSX or TXT", accept_multiple_files=True,
                               type=["vcf", "gz", "bgz", "csv", "tsv", "tab", "txt", "xlsx"])
    f1, f2, f3 = st.columns([2, 1.3, 1])
    genes_txt = f1.text_input("Only these genes (optional)", placeholder="TP53, BRCA1, BRCA2",
                              help="Keeps only variants inside these genes — much faster for exomes/genomes.")
    samples: list[str] = []
    for up in uploads or []:
        try:
            for s in vcf_samples(up.name, up.getvalue()):
                if s not in samples:
                    samples.append(s)
        except Exception:  # not a VCF
            pass
    sample = None
    if samples:
        choice = f2.selectbox("Sample", ["All carriers"] + samples, index=1 if len(samples) == 1 else 0)
        sample = None if choice == "All carriers" else choice
    pass_only = f3.checkbox("FILTER = PASS only")

run = st.button("Predict protein effect", type="primary", icon="🔬", width="stretch")
run = run or st.session_state.pop("_auto_run", False)

if run:
    files = [(u.name, u.getvalue()) for u in uploads or []]
    if not (text or "").strip() and not files:
        st.warning("Type at least one variant or upload a file.")
    else:
        bar = st.progress(0.0, text="Starting…")

        def cb(stage, done, total):
            bar.progress(min(1.0, done / max(total, 1)), text=f"{stage} · {done}/{total}")

        opts = PredictOptions(
            assembly=assembly, default_gene=gene or None,
            genes=[g for g in genes_txt.replace(";", ",").split(",") if g.strip()] or None,
            sample=sample, pass_only=pass_only, include_noncoding=noncoding, predictors=predictors,
            protein_context=context, clinvar=clinvar_on, workers=workers,
        )
        try:
            st.session_state["res"] = predict(text=text or "", files=files, options=opts, progress=cb)
        except Exception as e:  # show, don't crash the app
            st.error(f"Analysis failed: {e}")
        bar.empty()


# --------------------------------------------------------------------------- helpers for results

@st.cache_data(show_spinner=False, max_entries=64, ttl=86400)
def fetch_text(url: str) -> str | None:
    try:
        r = httpx.get(url, timeout=60, follow_redirects=True)
        return r.text if r.status_code == 200 else None
    except httpx.HTTPError:
        return None


def structure_html(pdb: str, pos: int, label: str, height: int) -> str:
    return f"""
<div id="mol" style="width:100%;height:{height}px;position:relative;border-radius:8px;overflow:hidden"></div>
<script src="https://cdn.jsdelivr.net/npm/3dmol@2.4.2/build/3Dmol-min.js"></script>
<script>
  const v = $3Dmol.createViewer(document.getElementById("mol"), {{backgroundColor: "white", backgroundAlpha: 0}});
  v.addModel({json.dumps(pdb)}, "pdb");
  const plddt = a => a.b >= 90 ? "#0053D6" : a.b >= 70 ? "#65CBF3" : a.b >= 50 ? "#FFDB13" : "#FF7D45";
  v.setStyle({{}}, {{cartoon: {{colorfunc: plddt}}}});
  v.addStyle({{resi: {pos}}}, {{stick: {{colorscheme: "magentaCarbon", radius: 0.28}}}});
  v.addStyle({{resi: {pos}, atom: "CA"}}, {{sphere: {{color: "#d6008f", radius: 1.1}}}});
  v.addLabel({json.dumps(label)}, {{backgroundColor: "#222", backgroundOpacity: 0.85, fontColor: "white",
             fontSize: 13, inFront: true, showBackground: true}}, {{resi: {pos}, atom: "CA"}});
  v.zoomTo({{within: {{distance: 18, sel: {{resi: {pos}}}}}}});
  v.render();
</script>"""


def lollipop(gdf: pd.DataFrame, length: int | None, domains: list[dict], selected: int | None):
    """Variants along the protein: x = residue, y = effect score, colour + shape = prediction."""
    d = gdf.copy()
    d["status"] = d["prediction"].map(status_of)
    d["score"] = d["effect_score"].where(d["effect_score"].notna(), 0.5)
    d["zero"] = 0.0
    d["label"] = d["protein_change"].fillna("")
    d["sel"] = d.index == selected
    names = list(STATUS)
    color = alt.Color("status:N", title=None,
                      scale=alt.Scale(domain=names, range=[STATUS[k][0] for k in names]),
                      legend=alt.Legend(orient="top", symbolSize=90))
    shape = alt.Shape("status:N", scale=alt.Scale(domain=names, range=[STATUS[k][1] for k in names]), legend=None)
    xmax = length or int(d["protein_position"].max() * 1.05) + 1
    xs = alt.Scale(domain=[0, xmax], nice=False)
    ys = alt.Scale(domain=[-0.18, 1.15], nice=False)
    x = alt.X("protein_position:Q", title="Residue", scale=xs)
    y = alt.Y("score:Q", title="Effect score", scale=ys, axis=alt.Axis(values=[0, 0.5, 1], grid=True))
    tooltip = [alt.Tooltip("label:N", title="Change"), alt.Tooltip("prediction:N", title="Prediction"),
               alt.Tooltip("effect_score:Q", title="Effect score", format=".2f"),
               alt.Tooltip("alphamissense:Q", title="AlphaMissense", format=".2f"),
               alt.Tooltip("clinvar:N", title="ClinVar"), alt.Tooltip("input:N", title="Input")]
    base = alt.Chart(d)
    stems = base.mark_rule(strokeWidth=2, opacity=0.6).encode(x=x, y=alt.Y("zero:Q", scale=ys), y2="score:Q", color=color)
    dark = getattr(getattr(st.context, "theme", None), "type", None) == "dark"
    surface, ink = ("#0e1117", "#fafafa") if dark else ("#ffffff", "#262730")
    pts = base.mark_point(filled=True, size=110, stroke=surface, strokeWidth=2, opacity=1).encode(
        x=x, y=y, color=color, shape=shape, tooltip=tooltip)
    labels = base.transform_filter(alt.datum.sel).mark_text(dy=-14, fontSize=12, fontWeight="bold", color=ink).encode(
        x=x, y=y, text="label:N")
    layers = [stems, pts, labels]
    if domains:
        dd = pd.DataFrame(domains)
        dd["y0"], dd["y1"], dd["ym"] = -0.15, -0.05, -0.10
        band = alt.Chart(dd).mark_rect(cornerRadius=4, opacity=0.3, color="#7f7f7a").encode(
            x=alt.X("start:Q", scale=xs), x2="end:Q", y=alt.Y("y0:Q", scale=ys), y2="y1:Q",
            tooltip=[alt.Tooltip("name:N", title="Domain"), alt.Tooltip("start:Q", title="From"),
                     alt.Tooltip("end:Q", title="To")])
        txt = alt.Chart(dd).mark_text(fontSize=10, color=ink).encode(
            x=alt.X("mid:Q", scale=xs), y=alt.Y("ym:Q", scale=ys), text="short:N")
        layers = [band, txt] + layers
    return alt.layer(*layers).properties(height=230).configure_view(strokeWidth=0)


def link(label, url):
    return f"[{label}]({url})"


# --------------------------------------------------------------------------- results

res = st.session_state.get("res")
if res is not None:
    df = res.table
    s = res.stats
    st.divider()
    m = st.columns(6)
    m[0].metric("Variants analysed", s.get("inputs", 0))
    m[1].metric("Protein-level results", s.get("rows", 0))
    m[2].metric("Likely damaging", s.get("likely_damaging", 0))
    m[3].metric("Not analysed", s.get("errors", 0))
    m[4].metric("Genome build", res.assembly)
    m[5].metric("Time", f"{s.get('seconds', 0)} s")
    for w in res.warnings:
        st.warning(w, icon="⚠️")
    skipped = [f"{s[k]} {lbl}" for k, lbl in [("no_protein_effect", "without coding/splice consequence"),
                                                 ("file_skipped_not_carried", "not carried by the sample(s)"),
                                                 ("file_skipped_outside_genes", "outside the selected genes"),
                                                 ("file_skipped_symbolic", "symbolic/structural alleles")] if s.get(k)]
    if skipped:
        st.caption("Skipped: " + " · ".join(skipped))

    if df.empty:
        st.info("No protein-level consequences found. Turn on *Include non-coding variants* in the sidebar to "
                "see everything VEP reports.")
    else:
        df = df.copy()
        df["status"] = df["prediction"].map(status_of)
        fc1, fc2 = st.columns([2, 1])
        present = [k for k in STATUS if k in set(df["status"])]
        show = fc1.segmented_control("Show", present, selection_mode="multi", default=present,
                                     format_func=lambda k: f"{STATUS[k][2]} {k}")
        q = fc2.text_input("Search", placeholder="gene, change or input", label_visibility="collapsed")
        view = df[df["status"].isin(show or present)]
        if q:
            ql = q.lower()
            view = view[view[["gene", "protein_change", "input", "hgvsp"]].astype(str)
                        .apply(lambda c: c.str.lower().str.contains(ql, regex=False)).any(axis=1)]

        disp = pd.DataFrame({
            "Input": view["input"],
            "Gene": view["gene"],
            "Change": view["protein_change"],
            "Consequence": view["consequence"].fillna("").str.replace("_variant", "").str.replace("_", " "),
            "Prediction": view["status"].map(lambda k: STATUS[k][2]) + " " + view["prediction"].fillna(""),
            "Effect": view["effect_score"],
            "AlphaMissense": view["alphamissense"],
            "REVEL": view["revel"],
            "CADD": view["cadd_phred"],
            "pLDDT": view["plddt"],
            "ClinVar": view["clinvar"].fillna("") + view["clinvar_stars"].map(
                lambda x: " " + "★" * int(x) if _isnum(x) and x > 0 else ""),
            "gnomAD AF": view["gnomad_af"],
            "Zygosity": view["zygosity"] if view["zygosity"].notna().any() else view["carriers"],
        })
        disp = disp.loc[:, [c for c in disp.columns if disp[c].notna().any() and (disp[c].astype(str) != "").any()]]
        sel = st.dataframe(
            disp, hide_index=True, on_select="rerun", selection_mode="single-row", key="grid",
            column_config={
                "Effect": st.column_config.ProgressColumn("Effect", min_value=0, max_value=1, format="%.2f",
                                                          help="Weighted fraction of predictors calling the variant damaging (1 = loss of function)"),
                "AlphaMissense": st.column_config.ProgressColumn("AlphaMissense", min_value=0, max_value=1, format="%.2f"),
                "REVEL": st.column_config.ProgressColumn("REVEL", min_value=0, max_value=1, format="%.2f"),
                "CADD": st.column_config.NumberColumn("CADD", format="%.1f"),
                "pLDDT": st.column_config.NumberColumn("pLDDT", format="%.0f", help="AlphaFold confidence at the residue"),
                "gnomAD AF": st.column_config.NumberColumn("gnomAD AF", format="%.2e"),
            },
        )
        st.caption("Select a row to see the protein-level details.")

        d1, d2, d3, d4, _ = st.columns([1, 1, 1, 1.3, 2])
        d1.download_button("CSV", to_bytes(res, "csv"), "missensexplorer.csv", "text/csv", icon="⬇️", on_click="ignore")
        d2.download_button("Excel", to_bytes(res, "xlsx"), "missensexplorer.xlsx",
                           "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet", icon="⬇️", on_click="ignore")
        d3.download_button("JSON", to_bytes(res, "json"), "missensexplorer.json", "application/json", icon="⬇️", on_click="ignore")
        if res.sequences:
            d4.download_button("Protein FASTA", to_fasta(res), "missensexplorer_proteins.fasta", "text/plain",
                               icon="🧬", on_click="ignore", help="Wild-type and mutant protein sequences")

        rows = sel.selection.rows if sel and sel.selection else []
        idx = view.index[rows[0]] if rows and rows[0] < len(view) else (view.index[0] if len(view) else None)

        if idx is not None:
            r = df.loc[idx]
            st.divider()
            stt = status_of(r["prediction"])
            title = f"{_val(r['gene'], '')} {_val(r['hgvsp'], '') or ''}"
            st.subheader(f"{title}  ({_val(r['protein_change'], '—')})")
            st.markdown(f"{STATUS[stt][2]} **{r['prediction']}** · "
                        f"{str(_val(r['consequence'], '')).replace('_', ' ')} · `{r['input']}`")
            left, right = st.columns([1.15, 1])

            with left:
                if _val(r["prediction_basis"]):
                    st.markdown(f"**Why:** {r['prediction_basis']}")
                preds = []
                for name, val, extra, rng in [
                    ("AlphaMissense", r["alphamissense"], r["alphamissense_class"], "0–1 · >0.56 likely pathogenic"),
                    ("REVEL", r["revel"], None, "0–1 · ≥0.644 supports pathogenic"),
                    ("CADD (phred)", r["cadd_phred"], None, "≥25.3 supports pathogenic"),
                    ("EVE", r["eve"], r["eve_class"], "0–1 · evolutionary model"),
                    ("popEVE", r["popeve"], None, "more negative = more deleterious"),
                    ("SIFT", r["sift_score"], r["sift"], "≤0.05 deleterious"),
                    ("PolyPhen-2", r["polyphen_score"], r["polyphen"], "≥0.85 probably damaging"),
                    ("SpliceAI (max Δ)", r["spliceai_max"], r["spliceai_effect"], "≥0.2 possible · ≥0.5 likely"),
                ]:
                    if _isnum(val) or _val(extra):
                        preds.append({"Predictor": name, "Score": val if _isnum(val) else None,
                                      "Call": str(_val(extra, "")).replace("_", " "), "Scale": rng})
                if preds:
                    st.dataframe(pd.DataFrame(preds), hide_index=True,
                                 column_config={"Score": st.column_config.NumberColumn(format="%.3f")})

                ra, aa = _val(r["ref_aa"]), _val(r["alt_aa"])
                if _isnum(r["grantham"]) and ra and aa:
                    st.markdown(
                        f"**{AA_NAMES.get(ra, ra)} → {AA_NAMES.get(aa, aa)}** ({r['aa_class_change']}) · "
                        f"Grantham {int(r['grantham'])} ({r['grantham_class']}) · BLOSUM62 {int(r['blosum62'])} · "
                        f"Δhydropathy {r['hydropathy_change']:+.1f} · Δvolume {r['volume_change']:+.0f} Å³")
                    if _val(r["physchem_notes"]):
                        st.caption(r["physchem_notes"])
                if _val(r["protein_window"]):
                    w = html.escape(r["protein_window"])
                    w = w.replace("[", "<span style='background:#d6008f;color:white;border-radius:4px;padding:0 4px'>").replace("]", "</span>")
                    st.markdown(f"<div style='font-family:monospace;font-size:1.05rem;letter-spacing:.06em'>{w}</div>",
                                unsafe_allow_html=True)
                    note = _val(r["mutant_protein"])
                    if note and note != "ok":
                        st.caption(note)
                if _val(r["nmd"]):
                    st.caption(r["nmd"])

                st.markdown("##### Protein context")
                facts = []
                if _val(r["protein_name"]):
                    facts.append(f"**{r['protein_name']}** · residue {int(r['protein_position'])}"
                                 + (f" of {int(r['protein_length'])}" if _isnum(r["protein_length"]) else ""))
                if _val(r["key_site"]):
                    facts.append(f"🎯 **Functional site:** {r['key_site']}")
                if _val(r["protein_features"]):
                    facts.append(f"**Region:** {r['protein_features']}")
                if _isnum(r["plddt"]):
                    facts.append(f"**AlphaFold pLDDT:** {r['plddt']:.0f} ({r['plddt_class']})")
                st.markdown("\n\n".join(facts) or "_No UniProt context for this transcript._")
                if _val(r["uniprot_variants_here"]):
                    with st.expander("Known UniProt variants at this residue"):
                        for item in str(r["uniprot_variants_here"]).split(" | "):
                            st.markdown(f"- {item}")
                if _val(r["uniprot_mutagenesis_here"]):
                    with st.expander("Experimental mutagenesis at this residue"):
                        for item in str(r["uniprot_mutagenesis_here"]).split(" | "):
                            st.markdown(f"- {item}")
                if _val(r["protein_function"]):
                    with st.expander("Protein function (UniProt)"):
                        st.write(r["protein_function"])

                st.markdown("##### Clinical & population")
                cl = []
                if _val(r["clinvar"]):
                    stars = "★" * int(r["clinvar_stars"]) if _isnum(r["clinvar_stars"]) else ""
                    txt = f"**ClinVar:** {r['clinvar']} {stars}"
                    if _val(r["clinvar_review"]):
                        txt += f" — _{r['clinvar_review']}_"
                    cl.append(txt)
                if _val(r["clinvar_oncogenicity"]):
                    cl.append(f"**Somatic / oncogenicity:** {r['clinvar_oncogenicity']}")
                if _val(r["clinvar_conditions"]):
                    cl.append(f"**Conditions:** {r['clinvar_conditions']}")
                if _isnum(r["gnomad_af"]):
                    cl.append(f"**gnomAD allele frequency:** {r['gnomad_af']:.2e}"
                              + (f" (highest population {r['gnomad_popmax_af']:.2e})" if _isnum(r["gnomad_popmax_af"]) else ""))
                else:
                    cl.append("**gnomAD:** not observed")
                if _val(r["zygosity"]):
                    cl.append(f"**Zygosity:** {r['zygosity']} ({r['genotype']}, {r['sample']})")
                elif _val(r["carriers"]):
                    cl.append(f"**Carriers:** {r['carriers']}")
                st.markdown("\n\n".join(cl))

                links = []
                if _val(r["uniprot"]):
                    links += [link("UniProt", f"https://www.uniprot.org/uniprotkb/{r['uniprot']}/entry"),
                              link("AlphaFold", f"https://alphafold.ebi.ac.uk/entry/{r['uniprot']}")]
                if _val(r["clinvar_ids"]):
                    links.append(link("ClinVar", f"https://www.ncbi.nlm.nih.gov/clinvar/variation/{int(str(r['clinvar_ids'])[3:])}/"))
                if _val(r["rsid"]):
                    rs = str(r["rsid"]).split(",")[0].strip()
                    links.append(link("dbSNP", f"https://www.ncbi.nlm.nih.gov/snp/{rs}"))
                ref, altn = _val(r["ref"]), _val(r["alt"])
                if ref and altn and "-" not in (ref + altn) and res.assembly == "GRCh38":
                    links.append(link("gnomAD", f"https://gnomad.broadinstitute.org/variant/{r['chrom']}-{int(r['pos'])}-{ref}-{altn}?dataset=gnomad_r4"))
                if _val(r["transcript"]):
                    links.append(link("Ensembl", f"https://www.ensembl.org/Homo_sapiens/Transcript/Summary?t={r['transcript']}"))
                if _val(r["pubmed_ids"]):
                    ids = str(r["pubmed_ids"]).replace(" ", "").split(",")
                    links.append(link(f"PubMed ({int(r['pubmed_count'])})",
                                      "https://pubmed.ncbi.nlm.nih.gov/?term=" + quote(" OR ".join(ids[:25]))))
                if links:
                    st.markdown(" · ".join(links))

            with right:
                af = res.alphafold.get(_val(r["uniprot"], "")) if _val(r["uniprot"]) else None
                pos = int(r["protein_position"]) if _isnum(r["protein_position"]) else None
                if af and af.get("available") and pos and _isnum(r["plddt"]):
                    with st.spinner("Loading AlphaFold structure…"):
                        pdb = fetch_text(af["pdb_url"])
                    if pdb:
                        components.html(structure_html(pdb, pos, str(_val(r["protein_change"], pos)), 430), height=440)
                        st.caption(f"AlphaFold model {af.get('entry_id')} (wild type), coloured by confidence: "
                                   "🟦 very high · 🩵 confident · 🟨 low · 🟧 very low. Magenta = mutated residue.")
                else:
                    st.info("No AlphaFold structure mapped to this residue.")

                gene_rows = df[(df["gene"] == r["gene"]) & df["protein_position"].notna()]
                if len(gene_rows):
                    entry = res.uniprot.get(_val(r["uniprot"], "")) or {}
                    length = int(r["protein_length"]) if _isnum(r["protein_length"]) else entry.get("length")
                    doms = []
                    if not _val(r["notes"]) or "skipped" not in str(r["notes"]):
                        for f in entry.get("features", []):
                            if f["type"] in ("Domain", "Zinc finger", "DNA binding", "Transmembrane"):
                                name = f["description"] or f["type"]
                                doms.append({"name": name, "start": f["start"], "end": f["end"],
                                             "mid": (f["start"] + f["end"]) / 2,
                                             "short": name if (f["end"] - f["start"]) > (length or 1000) * 0.12 else ""})
                    st.markdown(f"##### Variants along {r['gene']}")
                    st.altair_chart(lollipop(gene_rows, length, doms, idx), width="stretch")

    if len(res.errors):
        with st.expander(f"⚠️ {len(res.errors)} input(s) could not be analysed"):
            st.dataframe(res.errors, hide_index=True)

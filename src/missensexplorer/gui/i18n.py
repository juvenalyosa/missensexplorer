"""Spanish / English strings for the desktop app."""

from __future__ import annotations

import re

LANG = "es"

S = {
    # sidebar
    "brand_sub": ("Del cambio en el ADN al efecto en la proteína", "From DNA change to protein effect"),
    "variants": ("VARIANTES", "VARIANTS"),
    "variants_hint": ("Una o varias por línea · BRAF V600E · TP53: R175H, R248Q · rs113488022 · chr7:140753336 A>T",
                      "One or more per line · BRAF V600E · TP53: R175H, R248Q · rs113488022 · chr7:140753336 A>T"),
    "default_gene": ("Gen por defecto", "Default gene"),
    "default_gene_tip": ("Permite escribir cambios sin gen, como «R175H» o «c.524G>A».",
                         "Lets you type bare changes such as 'R175H' or 'c.524G>A'."),
    "example": ("Ejemplo", "Example"),
    "formats": ("Formatos", "Formats"),
    "files": ("ARCHIVOS", "FILES"),
    "drop": ("Arrastra un VCF, CSV, Excel o TXT\no haz clic para elegir", "Drop a VCF, CSV, Excel or TXT file\nor click to browse"),
    "only_genes": ("Solo estos genes (TP53, BRCA1…)", "Only these genes (TP53, BRCA1…)"),
    "only_genes_tip": ("Filtra el VCF localmente antes de consultar: mucho más rápido para exomas.",
                       "Filters the VCF locally before any request: much faster for exomes."),
    "sample": ("Muestra", "Sample"),
    "all_carriers": ("Todos los portadores", "All carriers"),
    "pass_only": ("Solo FILTER = PASS", "FILTER = PASS only"),
    "options": ("OPCIONES", "OPTIONS"),
    "options_summary": ("Genoma {build} · {n} de 4 fuentes activas", "Build {build} · {n} of 4 sources on"),
    "build": ("Genoma", "Genome build"),
    "build_auto": ("Automático", "Automatic"),
    "opt_predictors": ("Predictores de patogenicidad", "Pathogenicity predictors"),
    "opt_predictors_tip": ("AlphaMissense, REVEL, CADD, EVE/popEVE, SpliceAI, SIFT y PolyPhen-2.",
                           "AlphaMissense, REVEL, CADD, EVE/popEVE, SpliceAI, SIFT and PolyPhen-2."),
    "opt_context": ("Contexto proteico (UniProt, AlphaFold)", "Protein context (UniProt, AlphaFold)"),
    "opt_clinvar": ("Clasificación ClinVar", "ClinVar classification"),
    "opt_noncoding": ("Incluir variantes no codificantes", "Include non-coding variants"),
    "workers": ("Peticiones en paralelo", "Parallel requests"),
    "predict": ("Predecir efecto en la proteína", "Predict protein effect"),
    "cancel": ("Cancelar", "Cancel"),
    "cancelling": ("Cancelando…", "Cancelling…"),
    "research_only": ("Solo para investigación. Las predicciones computacionales no son un diagnóstico clínico.",
                      "Research use only. Computational predictions are not a clinical diagnosis."),
    "need_input": ("Escribe al menos una variante o añade un archivo.", "Type at least one variant or add a file."),
    # welcome
    "welcome_title": ("¿Cómo cambia tu variante la proteína?", "How does your variant change the protein?"),
    "welcome_lead": ("Escribe mutaciones o arrastra un VCF. MissensExplorer localiza el transcrito, calcula el cambio "
                     "de aminoácido y reúne predictores, estructura y datos clínicos en una sola vista.",
                     "Type mutations or drop a VCF. MissensExplorer finds the transcript, works out the amino-acid "
                     "change and brings predictors, structure and clinical data into one view."),
    "step1_t": ("1 · Entrada", "1 · Input"),
    "step1_b": ("Texto libre, HGVS, rsIDs, coordenadas o archivos VCF, CSV y Excel. El genoma se detecta solo.",
                "Free text, HGVS, rsIDs, coordinates or VCF, CSV and Excel files. The build is detected for you."),
    "step2_t": ("2 · Predicción", "2 · Prediction"),
    "step2_b": ("Votación ponderada de AlphaMissense, REVEL, CADD, EVE, SIFT y PolyPhen-2; SpliceAI para splicing.",
                "Weighted vote of AlphaMissense, REVEL, CADD, EVE, SIFT and PolyPhen-2; SpliceAI for splicing."),
    "step3_t": ("3 · Contexto", "3 · Context"),
    "step3_b": ("Dominio UniProt, estructura AlphaFold, ClinVar del alelo exacto, gnomAD y la proteína mutante.",
                "UniProt domain, AlphaFold structure, exact-allele ClinVar, gnomAD and the mutant protein."),
    "try": ("Prueba con:", "Try:"),
    # header / results
    "results": ("Resultados", "Results"),
    "subtitle": ("{inputs} entradas · {rows} resultados · {build} · {secs} s",
                 "{inputs} inputs · {rows} results · {build} · {secs} s"),
    "export": ("Exportar", "Export"),
    "not_analysed": ("{n} sin analizar", "{n} not analysed"),
    "search": ("Buscar gen, cambio o entrada", "Search gene, change or input"),
    "kpi_inputs": ("Variantes analizadas", "Variants analysed"),
    "kpi_rows": ("Cambios en proteína", "Protein-level results"),
    "kpi_damaging": ("Probablemente dañinas", "Likely damaging"),
    "kpi_clinvar": ("Patogénicas en ClinVar", "Pathogenic in ClinVar"),
    "skipped": ("Omitidas: {items}", "Skipped: {items}"),
    "skip_noprot": ("{n} sin consecuencia codificante", "{n} without coding/splice consequence"),
    "skip_carried": ("{n} no portadas por la muestra", "{n} not carried by the sample(s)"),
    "skip_genes": ("{n} fuera de los genes elegidos", "{n} outside the selected genes"),
    "skip_symbolic": ("{n} alelos simbólicos/estructurales", "{n} symbolic/structural alleles"),
    "no_rows": ("No se encontraron consecuencias a nivel de proteína. Activa «Incluir variantes no codificantes» "
                "para ver todo lo que reporta VEP.",
                "No protein-level consequences found. Turn on 'Include non-coding variants' to see everything VEP reports."),
    "errors_title": ("Entradas sin analizar", "Inputs not analysed"),
    "col_input": ("Entrada", "Input"), "col_gene": ("Gen", "Gene"), "col_change": ("Cambio", "Change"),
    "col_consequence": ("Consecuencia", "Consequence"), "col_prediction": ("Predicción", "Prediction"),
    "col_effect": ("Efecto", "Effect"), "col_clinvar": ("ClinVar", "ClinVar"), "col_gnomad": ("gnomAD AF", "gnomAD AF"),
    "col_zyg": ("Cigosidad", "Zygosity"), "col_source": ("Origen", "Source"), "col_error": ("Motivo", "Reason"),
    # detail
    "why": ("POR QUÉ", "WHY"),
    "predictors": ("PREDICTORES", "PREDICTORS"),
    "physchem": ("CAMBIO FISICOQUÍMICO", "PHYSICO-CHEMISTRY"),
    "sequence": ("SECUENCIA", "SEQUENCE"),
    "context": ("CONTEXTO PROTEICO", "PROTEIN CONTEXT"),
    "clinical": ("CLÍNICA Y POBLACIÓN", "CLINICAL & POPULATION"),
    "links": ("ENLACES", "LINKS"),
    "residue_of": ("residuo {pos} de {len}", "residue {pos} of {len}"),
    "residue": ("residuo {pos}", "residue {pos}"),
    "key_site": ("Sitio funcional", "Functional site"),
    "region": ("Región", "Region"),
    "plddt": ("Confianza AlphaFold (pLDDT)", "AlphaFold confidence (pLDDT)"),
    "no_context": ("Sin contexto de UniProt para este transcrito.", "No UniProt context for this transcript."),
    "uniprot_variants": ("Variantes conocidas en UniProt en este residuo", "Known UniProt variants at this residue"),
    "mutagenesis": ("Mutagénesis experimental en este residuo", "Experimental mutagenesis at this residue"),
    "function": ("Función de la proteína", "Protein function"),
    "clinvar": ("ClinVar", "ClinVar"),
    "oncogenicity": ("Somático / oncogenicidad", "Somatic / oncogenicity"),
    "conditions": ("Condiciones", "Conditions"),
    "gnomad": ("Frecuencia en gnomAD", "gnomAD allele frequency"),
    "gnomad_none": ("No observada en gnomAD", "Not observed in gnomAD"),
    "popmax": ("máx. por población {v}", "highest population {v}"),
    "zygosity": ("Cigosidad", "Zygosity"),
    "carriers": ("Portadores", "Carriers"),
    "grantham": ("Grantham", "Grantham"), "blosum": ("BLOSUM62", "BLOSUM62"),
    "hydropathy": ("Δ hidropatía", "Δ hydropathy"), "volume": ("Δ volumen", "Δ volume"), "charge": ("Δ carga", "Δ charge"),
    "select_row": ("Selecciona una fila para ver el detalle de la variante.", "Select a row to see the variant details."),
    "tab_structure": ("Estructura 3D", "3D structure"),
    "tab_lollipop": ("Variantes en la proteína", "Variants along the protein"),
    "tab_sequence": ("Proteína mutante", "Mutant protein"),
    "no_structure": ("No hay un modelo de AlphaFold asociado a este residuo.", "No AlphaFold model mapped to this residue."),
    "loading_structure": ("Cargando estructura de AlphaFold…", "Loading AlphaFold structure…"),
    "structure_caption": ("Modelo AlphaFold {id} (tipo silvestre), coloreado por confianza. Magenta: residuo mutado.",
                          "AlphaFold model {id} (wild type), coloured by confidence. Magenta: mutated residue."),
    "no_webengine": ("El visor 3D necesita QtWebEngine. Abre el modelo en el navegador:",
                     "The 3D viewer needs QtWebEngine. Open the model in your browser:"),
    "open_alphafold": ("Abrir en AlphaFold DB", "Open in AlphaFold DB"),
    "residue_axis": ("Residuo", "Residue"), "score_axis": ("Puntuación de efecto", "Effect score"),
    "no_sequence": ("No hay secuencia mutante para esta variante.", "No mutant sequence for this variant."),
    "copy_fasta": ("Copiar FASTA", "Copy FASTA"),
    "copied": ("Copiado al portapapeles", "Copied to the clipboard"),
    "wild_type": ("Silvestre", "Wild type"), "mutant": ("Mutante", "Mutant"),
    # menus / dialogs
    "m_file": ("Archivo", "File"), "m_open": ("Añadir archivos…", "Add files…"),
    "m_export_csv": ("Exportar CSV…", "Export CSV…"), "m_export_xlsx": ("Exportar Excel…", "Export Excel…"),
    "m_export_json": ("Exportar JSON…", "Export JSON…"), "m_export_fasta": ("Exportar FASTA de proteínas…", "Export protein FASTA…"),
    "m_quit": ("Salir", "Quit"), "m_view": ("Ver", "View"), "m_language": ("Idioma", "Language"),
    "m_help": ("Ayuda", "Help"), "m_formats": ("Formatos aceptados", "Accepted formats"),
    "m_about": ("Acerca de MissensExplorer", "About MissensExplorer"), "m_clear_cache": ("Vaciar caché", "Clear cache"),
    "cache_cleared": ("Caché vaciada", "Cache cleared"),
    "saved": ("Guardado: {path}", "Saved: {path}"),
    "lang_restart": ("El idioma cambiará al reiniciar la aplicación.", "The language changes when you restart the app."),
    "failed": ("El análisis falló", "Analysis failed"),
    "about": ("<b>MissensExplorer {v}</b><br>Predice cómo las variantes de ADN cambian la proteína humana.<br><br>"
              "Datos: Ensembl VEP, AlphaMissense, REVEL, CADD, EVE/popEVE, SpliceAI, SIFT, PolyPhen-2, UniProt, "
              "AlphaFold DB, ClinVar y gnomAD.<br><br>Universidad Simón Bolívar · Escuela Nacional del Deporte<br><br>"
              "<i>Solo para investigación.</i>",
              "<b>MissensExplorer {v}</b><br>Predicts how DNA variants change the human protein.<br><br>"
              "Data: Ensembl VEP, AlphaMissense, REVEL, CADD, EVE/popEVE, SpliceAI, SIFT, PolyPhen-2, UniProt, "
              "AlphaFold DB, ClinVar and gnomAD.<br><br>Universidad Simón Bolívar · Escuela Nacional del Deporte<br><br>"
              "<i>Research use only.</i>"),
    "formats_html": (
        "<table cellpadding='5'>"
        "<tr><td><code>BRAF V600E</code>, <code>BRAF:p.Val600Glu</code></td><td>cambio proteico en un gen</td></tr>"
        "<tr><td><code>TP53: R175H, R248Q</code></td><td>varios cambios en el mismo gen</td></tr>"
        "<tr><td><code>BRCA1:c.5096G&gt;A</code>, <code>NM_004333.6:c.1799T&gt;A</code></td><td>HGVS de cDNA</td></tr>"
        "<tr><td><code>chr7:140753336 A&gt;T</code>, <code>7-140753336-A-T</code></td><td>posición genómica y alelos</td></tr>"
        "<tr><td><code>NC_000007.14:g.140753336A&gt;T</code></td><td>HGVS genómico</td></tr>"
        "<tr><td><code>rs113488022</code></td><td>identificador dbSNP (todos los alelos)</td></tr>"
        "<tr><td><code>V600E</code></td><td>cambio suelto, con «Gen por defecto»</td></tr>"
        "<tr><td>VCF / VCF.gz, CSV, TSV, Excel, TXT</td><td>archivos (arrastrar y soltar)</td></tr></table>",
        "<table cellpadding='5'>"
        "<tr><td><code>BRAF V600E</code>, <code>BRAF:p.Val600Glu</code></td><td>protein change in a gene</td></tr>"
        "<tr><td><code>TP53: R175H, R248Q</code></td><td>several changes in one gene</td></tr>"
        "<tr><td><code>BRCA1:c.5096G&gt;A</code>, <code>NM_004333.6:c.1799T&gt;A</code></td><td>cDNA HGVS</td></tr>"
        "<tr><td><code>chr7:140753336 A&gt;T</code>, <code>7-140753336-A-T</code></td><td>genomic position and alleles</td></tr>"
        "<tr><td><code>NC_000007.14:g.140753336A&gt;T</code></td><td>genomic HGVS</td></tr>"
        "<tr><td><code>rs113488022</code></td><td>dbSNP id (all alleles)</td></tr>"
        "<tr><td><code>V600E</code></td><td>bare change, with 'Default gene'</td></tr>"
        "<tr><td>VCF / VCF.gz, CSV, TSV, Excel, TXT</td><td>files (drag and drop)</td></tr></table>"),
    # status names
    "Damaging": ("Dañina", "Damaging"), "Uncertain": ("Incierta", "Uncertain"),
    "Benign / silent": ("Benigna / silenciosa", "Benign / silent"), "Not assessed": ("Sin evaluar", "Not assessed"),
    "all": ("Todas", "All"),
    # predictor calls
    "call_damaging": ("dañino", "damaging"), "call_benign": ("benigno", "benign"), "call_neutral": ("intermedio", "intermediate"),
}

# Prediction labels produced by predict.classify()
PREDICTION_ES = {
    "Likely damaging": "Probablemente dañina",
    "Likely damaging (single predictor)": "Probablemente dañina (un solo predictor)",
    "Likely benign": "Probablemente benigna",
    "Likely benign (single predictor)": "Probablemente benigna (un solo predictor)",
    "Uncertain (conflicting predictors)": "Incierta (predictores en conflicto)",
    "Missense (no predictor scores)": "Missense (sin puntuaciones)",
    "Loss of function (truncating)": "Pérdida de función (truncante)",
    "Loss of function (canonical splice site)": "Pérdida de función (sitio de splicing canónico)",
    "Loss of function (start codon lost)": "Pérdida de función (codón de inicio perdido)",
    "Protein-altering in-frame (likely damaging)": "Cambio en marco (probablemente dañino)",
    "Protein-altering in-frame": "Cambio en marco de lectura",
    "Likely splice-altering": "Probable alteración del splicing",
    "Possible splice effect": "Posible efecto en el splicing",
    "Silent (no amino-acid change)": "Silenciosa (sin cambio de aminoácido)",
    "Splice region (no strong splice signal)": "Región de splicing (sin señal fuerte)",
    "Non-coding / no protein change": "No codificante / sin cambio proteico",
}

PHRASES_ES = [
    (r"NMD escape likely \(single-exon gene\)", "Probable escape del NMD (gen de un solo exón)"),
    (r"NMD escape likely \(last exon\)", "Probable escape del NMD (último exón)"),
    (r"NMD uncertain \(penultimate exon; depends on distance to last junction\)",
     "NMD incierto (penúltimo exón; depende de la distancia a la última unión)"),
    (r"NMD likely \(premature stop upstream of last exon\)", "NMD probable (parada prematura antes del último exón)"),
    (r"truncated at residue (\d+): (\d+) of (\d+) residues lost \((\d+)%\)",
     r"truncada en el residuo \1: se pierden \2 de \3 residuos (\4 %)"),
    (r"(\d+) damaging / (\d+) benign predictors \(weighted (\d+)% damaging\)",
     r"\1 predictores dañinos / \2 benignos (ponderado: \3 % dañino)"),
    (r"frameshift: downstream sequence changes, not reconstructed",
     "cambio del marco de lectura: la secuencia posterior cambia y no se reconstruye"),
    (r"HIGH impact consequence", "consecuencia de impacto ALTO"),
    (r"no scores available", "sin puntuaciones disponibles"),
    (r"\(may also affect splicing\)", "(también puede afectar el splicing)"),
    (r"large hydropathy shift \(more hydrophobic\)", "gran cambio de hidropatía (más hidrofóbico)"),
    (r"large hydropathy shift \(more polar\)", "gran cambio de hidropatía (más polar)"),
    (r"large size change \(bulkier\)", "gran cambio de tamaño (más voluminoso)"),
    (r"large size change \(smaller\)", "gran cambio de tamaño (más pequeño)"),
    (r"proline introduced \(breaks helices/strands\)", "se introduce prolina (rompe hélices/láminas)"),
    (r"proline lost", "se pierde una prolina"),
    (r"glycine lost \(backbone flexibility\)", "se pierde una glicina (flexibilidad del esqueleto)"),
    (r"glycine introduced", "se introduce glicina"),
    (r"cysteine lost \(possible disulfide/metal site\)", "se pierde una cisteína (posible puente disulfuro o sitio metálico)"),
    (r"cysteine introduced \(unpaired thiol\)", "se introduce cisteína (tiol libre)"),
    (r"\bcharge ([+-]\d) -> ([+-]\d)", r"carga \1 → \2"),
    (r"\bmoderately conservative\b", "moderadamente conservador"),
    (r"\bmoderately radical\b", "moderadamente radical"),
    (r"\bconservative\b", "conservador"),
    (r"\bradical\b", "radical"),
    (r"\bvery low \(likely disordered\)", "muy baja (probablemente desordenada)"),
    (r"\bvery high\b", "muy alta"),
    (r"\bconfident\b", "alta"),
    (r"^low$", "baja"),
    (r"\bheterozygous\b", "heterocigota"),
    (r"\bhomozygous\b", "homocigota"),
    (r"\bhemizygous\b", "hemicigota"),
    (r"UniProt canonical isoform differs from transcript; features mapped by local sequence match",
     "La isoforma canónica de UniProt difiere del transcrito; las anotaciones se asignaron por coincidencia local"),
    (r"UniProt canonical numbering differs from this transcript; UniProt/AlphaFold features skipped",
     "La numeración de UniProt difiere de este transcrito; se omitieron UniProt/AlphaFold"),
    (r"^Gene\(s\) not found in Ensembl: ", "Gen(es) no encontrados en Ensembl: "),
    (r"unknown gene or transcript symbol", "gen o transcrito desconocido"),
    (r"reference allele does not match the genome/transcript", "el alelo de referencia no coincide con el genoma/transcrito"),
    (r"needs a gene, e\.g\. 'BRAF V600E' \(or set the default gene\)",
     "falta el gen, p. ej. «BRAF V600E» (o define el gen por defecto)"),
    (r"unrecognised variant format", "formato de variante no reconocido"),
    (r"REF alleles match (\S+) \((\d+%)\) better than (\S+) \((\d+%)\); using (\S+)\.",
     r"Los alelos REF coinciden mejor con \1 (\2) que con \3 (\4); se usa \5."),
    (r"Only (\d+%) of checked REF alleles match (\S+)(.*?)\. Check the genome build and coordinates\.",
     r"Solo el \1 de los alelos REF revisados coincide con \2\3. Revisa el genoma y las coordenadas."),
    (r"Multi-sample VCF \((\d+) samples\): carriers listed per variant\. Choose one sample to get its zygosity\.",
     r"VCF multimuestra (\1 muestras): se listan los portadores por variante. Elige una muestra para ver su cigosidad."),
    (r"could not parse", "no se pudo interpretar"),
]

AA_CLASS_ES = {"hydrophobic": "hidrofóbico", "aromatic": "aromático", "polar": "polar", "positive": "positivo",
               "negative": "negativo", "glycine": "glicina", "proline": "prolina"}

AA_NAMES_ES = {
    "A": "Alanina", "R": "Arginina", "N": "Asparagina", "D": "Aspartato", "C": "Cisteína", "Q": "Glutamina",
    "E": "Glutamato", "G": "Glicina", "H": "Histidina", "I": "Isoleucina", "L": "Leucina", "K": "Lisina",
    "M": "Metionina", "F": "Fenilalanina", "P": "Prolina", "S": "Serina", "T": "Treonina", "W": "Triptófano",
    "Y": "Tirosina", "V": "Valina", "*": "Stop",
}
_PHRASES = [(re.compile(a), b) for a, b in PHRASES_ES]


STAGES_ES = {
    "Reading input": "Leyendo la entrada",
    "Checking genome build": "Comprobando el genoma",
    "Annotating variants (VEP)": "Anotando variantes (VEP)",
    "Annotating + predictors (VEP)": "Anotando y puntuando (VEP)",
    "Annotating variants (VEP, fast pass)": "Anotando variantes (VEP, pasada rápida)",
    "Scoring coding variants (predictors)": "Puntuando variantes codificantes",
    "ClinVar classification": "Clasificación ClinVar",
    "Protein context (UniProt, AlphaFold)": "Contexto proteico (UniProt, AlphaFold)",
    "Done": "Listo",
}


CALLS_ES = {
    "likely_pathogenic": "probablemente patogénica", "likely_benign": "probablemente benigna", "ambiguous": "ambigua",
    "pathogenic": "patogénica", "benign": "benigna", "uncertain": "incierta",
    "deleterious": "deletérea", "deleterious_low_confidence": "deletérea (baja confianza)",
    "tolerated": "tolerada", "tolerated_low_confidence": "tolerada (baja confianza)",
    "probably_damaging": "probablemente dañina", "possibly_damaging": "posiblemente dañina", "unknown": "desconocida",
    "likely splice-altering": "probable alteración del splicing", "possible": "posible",
}


def tr_call(call: str | None) -> str:
    if not call:
        return ""
    if LANG == "es":
        return CALLS_ES.get(call.lower(), call.replace("_", " "))
    return call.replace("_", " ")


def tr_stage(stage: str) -> str:
    return STAGES_ES.get(stage, stage) if LANG == "es" else stage


def set_language(lang: str) -> None:
    global LANG
    LANG = "es" if lang.lower().startswith("es") else "en"


def t(key: str, **kw) -> str:
    pair = S.get(key)
    if pair is None:
        return key.format(**kw) if kw else key
    s = pair[0] if LANG == "es" else pair[1]
    return s.format(**kw) if kw else s


def tr_prediction(label: str | None) -> str:
    if not label:
        return ""
    return PREDICTION_ES.get(label, label) if LANG == "es" else label


def aa_name(code: str) -> str:
    from ..aa import AA_NAMES

    return (AA_NAMES_ES if LANG == "es" else AA_NAMES).get(code, code)


def aa_class(name: str | None) -> str:
    if not name:
        return ""
    return AA_CLASS_ES.get(name, name) if LANG == "es" else name


def tr_text(text: str | None) -> str:
    """Translate the fixed English phrases the pipeline produces (Spanish UI only)."""
    if not text:
        return ""
    if LANG != "es":
        return text
    out = text
    for rx, rep in _PHRASES:
        out = rx.sub(rep, out)
    return out

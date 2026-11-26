"""Spanish / English strings for the desktop app."""

from __future__ import annotations

import re
from collections.abc import Callable
from typing import TypeAlias

_Replacement: TypeAlias = str | Callable[[re.Match[str]], str]
"""What a matched phrase becomes: a ``re.sub`` template, or a function of the match."""

LANG = "es"

S: dict[str, tuple[str, str]] = {
    # sidebar
    "brand_sub": ("Del cambio en el ADN al efecto en la proteína", "From DNA change to protein effect"),
    "variants": ("VARIANTES", "VARIANTS"),
    "variants_hint": (
        "Una o varias por línea · BRAF V600E · TP53: R175H, R248Q · rs113488022 · chr7:140753336 A>T",
        "One or more per line · BRAF V600E · TP53: R175H, R248Q · rs113488022 · chr7:140753336 A>T",
    ),
    "default_gene": ("Gen por defecto", "Default gene"),
    "default_gene_tip": (
        "Permite escribir cambios sin gen, como «R175H» o «c.524G>A».",
        "Lets you type bare changes such as 'R175H' or 'c.524G>A'.",
    ),
    "example": ("Ejemplo", "Example"),
    "formats": ("Formatos", "Formats"),
    "files": ("ARCHIVOS", "FILES"),
    "drop": (
        "Arrastra un VCF, CSV, Excel o TXT\no haz clic para elegir",
        "Drop a VCF, CSV, Excel or TXT file\nor click to browse",
    ),
    "only_genes": ("Solo estos genes (TP53, BRCA1…)", "Only these genes (TP53, BRCA1…)"),
    "only_genes_tip": (
        "Filtra el VCF localmente antes de consultar: mucho más rápido para exomas.",
        "Filters the VCF locally before any request: much faster for exomes.",
    ),
    "sample": ("Muestra", "Sample"),
    "all_carriers": ("Todos los portadores", "All carriers"),
    "pass_only": ("Solo FILTER = PASS", "FILTER = PASS only"),
    "options": ("OPCIONES", "OPTIONS"),
    "options_summary": ("Genoma {build} · {n} de 4 fuentes activas", "Build {build} · {n} of 4 sources on"),
    "build": ("Genoma", "Genome build"),
    "build_auto": ("Automático", "Automatic"),
    "opt_predictors": ("Predictores de patogenicidad", "Pathogenicity predictors"),
    "opt_predictors_tip": (
        "AlphaMissense, REVEL, CADD, EVE/popEVE, SpliceAI, SIFT y PolyPhen-2.",
        "AlphaMissense, REVEL, CADD, EVE/popEVE, SpliceAI, SIFT and PolyPhen-2.",
    ),
    "opt_context": ("Contexto proteico (UniProt, AlphaFold)", "Protein context (UniProt, AlphaFold)"),
    "opt_clinvar": ("Clasificación ClinVar", "ClinVar classification"),
    "opt_noncoding": ("Incluir variantes no codificantes", "Include non-coding variants"),
    "opt_oncogenicity": ("Oncogenicidad somática (tumores)", "Somatic oncogenicity (tumours)"),
    "opt_oncogenicity_tip": (
        "Puntos de oncogenicidad ClinGen/CGC/VICC para variantes de tumores: hotspots de cancerhotspots.org, el "
        "papel del gen (IntOGen, o UniProt cuando IntOGen no decide), clasificaciones de oncogenicidad de ClinVar, "
        "estudios funcionales de CIViC, gnomAD y la evidencia computacional calibrada.",
        "ClinGen/CGC/VICC oncogenicity points for tumour variants: cancerhotspots.org hotspots, the gene's role "
        "(IntOGen, or UniProt when IntOGen does not decide), ClinVar oncogenicity classifications, CIViC functional "
        "studies, gnomAD and the calibrated computational evidence.",
    ),
    "workers": ("Peticiones en paralelo", "Parallel requests"),
    "freq_rule": ("Frecuencia", "Frequency"),
    "freq_ba1": ("≥ 5 % en gnomAD → benigna (BA1)", "≥ 5% in gnomAD → benign (BA1)"),
    "freq_bs1": ("≥ 1 % en gnomAD → benigna (BS1)", "≥ 1% in gnomAD → benign (BS1)"),
    "freq_off": ("No usar la frecuencia", "Do not use frequency"),
    "freq_tip": (
        "Una variante común en gnomAD (grupo de ancestría con mayor frecuencia) se marca como benigna aunque "
        "los predictores digan otra cosa. 5 % es el criterio BA1 de ACMG; 1 % es más estricto, pero confunde "
        "variantes patogénicas de enfermedades recesivas o de efecto fundador.",
        "A variant that is common in gnomAD (highest genetic ancestry group) is called benign whatever the "
        "predictors say. 5% is the ACMG BA1 criterion; 1% is stricter but mislabels pathogenic variants of "
        "recessive or founder conditions.",
    ),
    "acmg_src": ("PP3/BP4", "PP3/BP4"),
    "popeve_note": ("más negativo = más deletéreo", "more negative = more deleterious"),
    "acmg_src_model": ("Puntaje del modelo (calibración propia)", "Model score (own calibration)"),
    "acmg_src_revel": ("REVEL (umbrales de ClinGen)", "REVEL (ClinGen thresholds)"),
    "acmg_src_tip": (
        "Puntaje usado para la evidencia computacional ACMG PP3/BP4. El del modelo se calibró con el método "
        "de ClinGen (Pejaver et al. 2022) sobre los mismos datos de ClinVar que REVEL; la calibración es "
        "de este proyecto, no publicada por ClinGen. REVEL usa los umbrales que publicó ClinGen "
        "(AlphaMissense si falta REVEL).",
        "Score used for the ACMG computational evidence PP3/BP4. The model's was calibrated with the "
        "ClinGen method (Pejaver et al. 2022) on the same ClinVar data as REVEL's; the calibration is this "
        "project's, not one published by ClinGen. REVEL uses the thresholds ClinGen published "
        "(AlphaMissense when REVEL is missing).",
    ),
    "engine": ("Anotación", "Annotation"),
    "engine_rest": ("Ensembl en línea", "Ensembl online"),
    "engine_local": ("VEP local", "Local VEP"),
    "engine_tip": (
        "En línea usa la API pública de Ensembl. VEP local usa una instalación de VEP en este equipo "
        "(nativa o Docker): sin límites de uso y los datos no salen del equipo.",
        "Online uses the public Ensembl API. Local VEP uses a VEP install on this computer (native or "
        "Docker): no usage limits and the data never leave the computer.",
    ),
    "vep_cmd": ("Comando de VEP (p. ej. vep --cache --offline)", "VEP command (e.g. vep --cache --offline)"),
    "vep_args": ("Argumentos extra (plugins…)", "Extra arguments (plugins…)"),
    "offline": ("Sin conexión (modo privado)", "Offline (private mode)"),
    "offline_tip": (
        "No usa internet: solo VEP local. Omite ClinVar, UniProt y AlphaFold.",
        "Never uses the internet: local VEP only. Skips ClinVar, UniProt and AlphaFold.",
    ),
    "retry": ("Reintentar", "Retry"),
    "m_run_info": ("Información del análisis", "Run information"),
    "m_open_log": ("Abrir registro de errores", "Open error log"),
    "run_info_title": ("Información del análisis", "Run information"),
    "unexpected": ("Error inesperado", "Unexpected error"),
    "unexpected_body": (
        "Ocurrió un error inesperado. Quedó guardado en el registro:\n{path}",
        "An unexpected error occurred. It was saved to the log:\n{path}",
    ),
    "local_label": ("VEP local", "local VEP"),
    "predict": ("Predecir efecto en la proteína", "Predict protein effect"),
    "cancel": ("Cancelar", "Cancel"),
    "cancelling": ("Cancelando…", "Cancelling…"),
    "research_only": (
        "Solo para investigación. Las predicciones computacionales no son un diagnóstico clínico.",
        "Research use only. Computational predictions are not a clinical diagnosis.",
    ),
    "need_input": ("Escribe al menos una variante o añade un archivo.", "Type at least one variant or add a file."),
    # welcome
    "welcome_title": ("¿Cómo cambia tu variante la proteína?", "How does your variant change the protein?"),
    "welcome_lead": (
        "Escribe mutaciones o arrastra un VCF. MissensExplorer localiza el transcrito, calcula el cambio "
        "de aminoácido y reúne predictores, estructura y datos clínicos en una sola vista.",
        "Type mutations or drop a VCF. MissensExplorer finds the transcript, works out the amino-acid "
        "change and brings predictors, structure and clinical data into one view.",
    ),
    "step1_t": ("1 · Entrada", "1 · Input"),
    "step1_b": (
        "Texto libre, HGVS, rsIDs, coordenadas o archivos VCF, CSV y Excel. El genoma se detecta solo.",
        "Free text, HGVS, rsIDs, coordinates or VCF, CSV and Excel files. The build is detected for you.",
    ),
    "step2_t": ("2 · Predicción", "2 · Prediction"),
    "step2_b": (
        "Un modelo entrenado con ClinVar combina AlphaMissense, REVEL, CADD y popEVE; SpliceAI para "
        "splicing; criterios ACMG.",
        "A model trained on ClinVar combines AlphaMissense, REVEL, CADD and popEVE; SpliceAI for splicing; "
        "ACMG criteria.",
    ),
    "step3_t": ("3 · Contexto", "3 · Context"),
    "step3_b": (
        "Dominio UniProt, estructura AlphaFold, ClinVar del alelo exacto, gnomAD y la proteína mutante.",
        "UniProt domain, AlphaFold structure, exact-allele ClinVar, gnomAD and the mutant protein.",
    ),
    "try": ("Prueba con:", "Try:"),
    # header / results
    "results": ("Resultados", "Results"),
    "subtitle": (
        "{inputs} entradas · {rows} resultados · {build} · {engine} · {secs} s",
        "{inputs} inputs · {rows} results · {build} · {engine} · {secs} s",
    ),
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
    "skip_filter": ("{n} sin FILTER = PASS", "{n} not passing FILTER"),
    "no_rows": (
        "No se encontraron consecuencias a nivel de proteína. Activa «Incluir variantes no codificantes» "
        "para ver todo lo que reporta VEP.",
        "No protein-level consequences found. Turn on 'Include non-coding variants' to see everything VEP reports.",
    ),
    "errors_title": ("Entradas sin analizar", "Inputs not analysed"),
    "col_input": ("Entrada", "Input"),
    "col_gene": ("Gen", "Gene"),
    "col_change": ("Cambio", "Change"),
    "col_consequence": ("Consecuencia", "Consequence"),
    "col_prediction": ("Predicción", "Prediction"),
    "col_effect": ("Efecto", "Effect"),
    "col_with_freq": ("Con frecuencia", "With frequency"),
    "col_clinvar": ("ClinVar", "ClinVar"),
    "col_oncogenicity": ("Oncogenicidad", "Oncogenicity"),
    "col_gnomad": ("gnomAD AF", "gnomAD AF"),
    "col_zyg": ("Cigosidad", "Zygosity"),
    "col_source": ("Origen", "Source"),
    "col_error": ("Motivo", "Reason"),
    # detail
    "why": ("POR QUÉ", "WHY"),
    "contrib_damaging": ("dañina", "damaging"),
    "contrib_benign": ("benigna", "benign"),
    "model_note": (
        "Aporte de cada predictor al puntaje del modelo (log-odds). El modelo es una regresión logística "
        "entrenada con variantes de ClinVar revisadas por expertos.",
        "Each predictor's contribution to the model score (log-odds). The model is a logistic regression "
        "trained on expert-reviewed ClinVar variants.",
    ),
    "acmg": ("EVIDENCIA ACMG (ORIENTATIVA)", "ACMG EVIDENCE (INDICATIVE)"),
    "acmg_pp3": ("Computacional PP3/BP4", "Computational PP3/BP4"),
    "acmg_freq": ("Frecuencia poblacional", "Population frequency"),
    "acmg_codes": ("Criterios", "Criteria"),
    "acmg_vcep": ("Especificación del panel de ClinGen", "ClinGen expert panel specification"),
    "acmg_notes": ("Avisos sobre estos criterios", "Notes on these criteria"),
    "with_freq_tip": (
        "Puntaje del modelo combinado con la frecuencia en gnomAD, para ordenar las variantes de una muestra. "
        "Una variante ausente de gnomAD puntúa alto casi siempre. No es evidencia ACMG: la frecuencia ya "
        "cuenta en BA1, BS1 y PM2.",
        "Model score combined with the gnomAD frequency, for ranking a sample's variants. A variant absent "
        "from gnomAD almost always scores high. Not ACMG evidence: frequency already counts in BA1, BS1 "
        "and PM2.",
    ),
    "acmg_ps1pm5": ("Mismo residuo PS1/PM5", "Same residue PS1/PM5"),
    "acmg_pm2": ("Ausencia PM2", "Absence PM2"),
    "acmg_note_residue": (
        "PS1/PM5: variantes missense patogénicas de ClinVar en el mismo residuo (no se comprueba si su "
        "clasificación es independiente de PS1/PM5). PM2: ausencia confirmada en gnomAD en una posición "
        "con lecturas suficientes; si la variante está presente, depende de la enfermedad.",
        "PS1/PM5: pathogenic ClinVar missense variants at the same residue (whether their "
        "classification is independent of PS1/PM5 is not checked). PM2: absence confirmed in gnomAD "
        "at a position with enough reads; if the variant is present, it depends on the disease.",
    ),
    "points": ("{n} puntos", "{n} points"),
    "acmg_note": (
        "Umbrales calibrados por ClinGen (REVEL; AlphaMissense si falta REVEL). Es solo una parte de la "
        "clasificación ACMG: no reemplaza la revisión de un especialista.",
        "ClinGen-calibrated thresholds (REVEL; AlphaMissense when REVEL is missing). This is only part of "
        "an ACMG classification: it does not replace review by a specialist.",
    ),
    "acmg_note_model": (
        "Puntaje del modelo con umbrales calibrados por el método de ClinGen (Pejaver et al. 2022) "
        "sobre los mismos datos que REVEL; calibración propia, no publicada por ClinGen. Es solo una "
        "parte de la clasificación ACMG: no reemplaza la revisión de un especialista.",
        "Model score with thresholds calibrated by the ClinGen method (Pejaver et al. 2022) on the same "
        "data as REVEL's; this project's calibration, not one published by ClinGen. This is only part "
        "of an ACMG classification: it does not replace review by a specialist.",
    ),
    # Rows without PP3/BP4 (it applies to missense variants only), so no note on its thresholds.
    "acmg_note_general": (
        "Es solo una parte de la clasificación ACMG: no reemplaza la revisión de un especialista.",
        "This is only part of an ACMG classification: it does not replace review by a specialist.",
    ),
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
    "somatic": ("Oncogenicidad somática (automática)", "Somatic oncogenicity (automatic)"),
    "somatic_class": ("Clase", "Class"),
    "somatic_codes": ("Criterios", "Criteria"),
    "somatic_evidence": ("Evidencia", "Evidence"),
    "somatic_clinvar": ("ClinVar, esta misma variante", "ClinVar, this very variant"),
    "somatic_clinical": ("Evidencia clínica (CIViC)", "Clinical evidence (CIViC)"),
    "somatic_note": (
        "Norma ClinGen/CGC/VICC (Horak et al. 2022) con los criterios que se calculan solos: hotspots de "
        "cancerhotspots.org, papel del gen (IntOGen, o UniProt cuando IntOGen no decide), otras variantes del "
        "residuo en ClinVar, estudios funcionales curados en CIViC (OS2/SBS2), gnomAD (OP4 solo con ausencia "
        "confirmada) y la evidencia computacional calibrada (OP1/SBP1, en todos los genes). CIViC cubre pocos "
        "estudios funcionales: revisa la literatura. La evidencia clínica de CIViC (fármacos, pronóstico) no suma "
        "puntos, y la clasificación de ClinVar de esta misma variante se muestra aparte y no cuenta como evidencia.",
        "ClinGen/CGC/VICC standard (Horak et al. 2022) with the criteria that can be computed: hotspots from "
        "cancerhotspots.org, the gene's role (IntOGen, or UniProt when IntOGen does not decide), other ClinVar "
        "variants at the residue, functional studies curated in CIViC (OS2/SBS2), gnomAD (OP4 only on confirmed "
        "absence) and the calibrated computational evidence (OP1/SBP1, in every gene). CIViC covers few functional "
        "studies: check the literature. CIViC's clinical evidence (drugs, prognosis) adds no points, and ClinVar's "
        "classification of this very variant is shown apart and does not count as evidence.",
    ),
    "conditions": ("Condiciones", "Conditions"),
    "gnomad": ("Frecuencia en gnomAD", "gnomAD allele frequency"),
    "gnomad_none": ("No observada en gnomAD", "Not observed in gnomAD"),
    "popmax": ("máx. por población {v}", "highest population {v}"),
    "zygosity": ("Cigosidad", "Zygosity"),
    "carriers": ("Portadores", "Carriers"),
    "grantham": ("Grantham", "Grantham"),
    "blosum": ("BLOSUM62", "BLOSUM62"),
    "hydropathy": ("Δ hidropatía", "Δ hydropathy"),
    "volume": ("Δ volumen", "Δ volume"),
    "charge": ("Δ carga", "Δ charge"),
    "select_row": (
        "Selecciona una fila para ver el detalle de la variante.",
        "Select a row to see the variant details.",
    ),
    "tab_structure": ("Estructura 3D", "3D structure"),
    "tab_lollipop": ("Variantes en la proteína", "Variants along the protein"),
    "tab_sequence": ("Proteína mutante", "Mutant protein"),
    "no_structure": (
        "No hay un modelo de AlphaFold asociado a este residuo.",
        "No AlphaFold model mapped to this residue.",
    ),
    "loading_structure": ("Cargando estructura de AlphaFold…", "Loading AlphaFold structure…"),
    "structure_caption": (
        "Modelo AlphaFold {id} (tipo silvestre), coloreado por confianza. Magenta: residuo mutado.",
        "AlphaFold model {id} (wild type), coloured by confidence. Magenta: mutated residue.",
    ),
    "no_webengine": (
        "El visor 3D necesita QtWebEngine. Abre el modelo en el navegador:",
        "The 3D viewer needs QtWebEngine. Open the model in your browser:",
    ),
    "open_alphafold": ("Abrir en AlphaFold DB", "Open in AlphaFold DB"),
    "residue_axis": ("Residuo", "Residue"),
    "score_axis": ("Puntuación de efecto", "Effect score"),
    "no_sequence": ("No hay secuencia mutante para esta variante.", "No mutant sequence for this variant."),
    "copy_fasta": ("Copiar FASTA", "Copy FASTA"),
    "copied": ("Copiado al portapapeles", "Copied to the clipboard"),
    "wild_type": ("Silvestre", "Wild type"),
    "mutant": ("Mutante", "Mutant"),
    # menus / dialogs
    "m_file": ("Archivo", "File"),
    "m_open": ("Añadir archivos…", "Add files…"),
    "m_export_csv": ("Exportar CSV…", "Export CSV…"),
    "m_export_xlsx": ("Exportar Excel…", "Export Excel…"),
    "m_export_json": ("Exportar JSON…", "Export JSON…"),
    "m_export_fasta": ("Exportar FASTA de proteínas…", "Export protein FASTA…"),
    "m_quit": ("Salir", "Quit"),
    "m_view": ("Ver", "View"),
    "m_language": ("Idioma", "Language"),
    "m_help": ("Ayuda", "Help"),
    "m_formats": ("Formatos aceptados", "Accepted formats"),
    "m_about": ("Acerca de MissensExplorer", "About MissensExplorer"),
    "m_clear_cache": ("Vaciar caché", "Clear cache"),
    "cache_cleared": ("Caché vaciada", "Cache cleared"),
    "saved": ("Guardado: {path}", "Saved: {path}"),
    "lang_restart": (
        "El idioma cambiará al reiniciar la aplicación.",
        "The language changes when you restart the app.",
    ),
    "failed": ("El análisis falló", "Analysis failed"),
    "about": (
        "<b>MissensExplorer {v}</b><br>Predice cómo las variantes de ADN cambian la proteína humana.<br><br>"
        "Datos: Ensembl VEP, AlphaMissense, REVEL, CADD, EVE/popEVE, SpliceAI, SIFT, PolyPhen-2, UniProt, "
        "AlphaFold DB, ClinVar y gnomAD.<br><br>Universidad Simón Bolívar · Escuela Nacional del Deporte<br><br>"
        "<i>Solo para investigación.</i>",
        "<b>MissensExplorer {v}</b><br>Predicts how DNA variants change the human protein.<br><br>"
        "Data: Ensembl VEP, AlphaMissense, REVEL, CADD, EVE/popEVE, SpliceAI, SIFT, PolyPhen-2, UniProt, "
        "AlphaFold DB, ClinVar and gnomAD.<br><br>Universidad Simón Bolívar · Escuela Nacional del Deporte<br><br>"
        "<i>Research use only.</i>",
    ),
    "formats_html": (
        "<table cellpadding='5'>"
        "<tr><td><code>BRAF V600E</code>, <code>BRAF:p.Val600Glu</code></td><td>cambio proteico en un gen</td></tr>"
        "<tr><td><code>TP53: R175H, R248Q</code></td><td>varios cambios en el mismo gen</td></tr>"
        "<tr><td><code>BRCA1:c.5096G&gt;A</code>, <code>NM_004333.6:c.1799T&gt;A</code></td><td>HGVS de cDNA</td></tr>"
        "<tr><td><code>chr7:140753336 A&gt;T</code>, <code>7-140753336-A-T</code></td>"
        "<td>posición genómica y alelos</td></tr>"
        "<tr><td><code>NC_000007.14:g.140753336A&gt;T</code></td><td>HGVS genómico</td></tr>"
        "<tr><td><code>rs113488022</code></td><td>identificador dbSNP (todos los alelos)</td></tr>"
        "<tr><td><code>V600E</code></td><td>cambio suelto, con «Gen por defecto»</td></tr>"
        "<tr><td>VCF / VCF.gz, CSV, TSV, Excel, TXT</td><td>archivos (arrastrar y soltar)</td></tr></table>",
        "<table cellpadding='5'>"
        "<tr><td><code>BRAF V600E</code>, <code>BRAF:p.Val600Glu</code></td><td>protein change in a gene</td></tr>"
        "<tr><td><code>TP53: R175H, R248Q</code></td><td>several changes in one gene</td></tr>"
        "<tr><td><code>BRCA1:c.5096G&gt;A</code>, <code>NM_004333.6:c.1799T&gt;A</code></td><td>cDNA HGVS</td></tr>"
        "<tr><td><code>chr7:140753336 A&gt;T</code>, <code>7-140753336-A-T</code></td>"
        "<td>genomic position and alleles</td></tr>"
        "<tr><td><code>NC_000007.14:g.140753336A&gt;T</code></td><td>genomic HGVS</td></tr>"
        "<tr><td><code>rs113488022</code></td><td>dbSNP id (all alleles)</td></tr>"
        "<tr><td><code>V600E</code></td><td>bare change, with 'Default gene'</td></tr>"
        "<tr><td>VCF / VCF.gz, CSV, TSV, Excel, TXT</td><td>files (drag and drop)</td></tr></table>",
    ),
    # status names
    "Damaging": ("Dañina", "Damaging"),
    "Uncertain": ("Incierta", "Uncertain"),
    "Benign / silent": ("Benigna / silenciosa", "Benign / silent"),
    "Not assessed": ("Sin evaluar", "Not assessed"),
    "all": ("Todas", "All"),
    # predictor calls
    "call_damaging": ("dañino", "damaging"),
    "call_benign": ("benigno", "benign"),
    "call_neutral": ("intermedio", "intermediate"),
}

# Prediction labels produced by predict.classify()
PREDICTION_ES: dict[str, str] = {
    "Likely damaging": "Probablemente dañina",
    "Likely damaging (single predictor)": "Probablemente dañina (un solo predictor)",
    "Likely benign": "Probablemente benigna",
    "Likely benign (single predictor)": "Probablemente benigna (un solo predictor)",
    "Uncertain (conflicting predictors)": "Incierta (predictores en conflicto)",
    "Uncertain (intermediate scores)": "Incierta (puntuaciones intermedias)",
    "Likely benign (common in population)": "Probablemente benigna (frecuente en la población)",
    "Likely damaging (limited evidence: SIFT/PolyPhen only)": (
        "Probablemente dañina (evidencia limitada: solo SIFT/PolyPhen)"
    ),
    "Likely benign (limited evidence: SIFT/PolyPhen only)": (
        "Probablemente benigna (evidencia limitada: solo SIFT/PolyPhen)"
    ),
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

PHRASES_ES: list[tuple[str, _Replacement]] = [
    (r"NMD escape likely \(single-exon gene\)", "Probable escape del NMD (gen de un solo exón)"),
    (r"NMD escape likely \(last exon\)", "Probable escape del NMD (último exón)"),
    (
        r"NMD uncertain \(penultimate exon; depends on distance to last junction\)",
        "NMD incierto (penúltimo exón; depende de la distancia a la última unión)",
    ),
    (r"NMD likely \(premature stop upstream of last exon\)", "NMD probable (parada prematura antes del último exón)"),
    (
        r"NMD uncertain \(the change spans an exon-intron boundary\)",
        "NMD incierto (el cambio abarca una unión exón-intrón)",
    ),
    (
        r"truncated at residue (\d+): (\d+) of (\d+) residues lost \((\d+)%\)",
        r"truncada en el residuo \1: se pierden \2 de \3 residuos (\4 %)",
    ),
    (
        r"(\d+) damaging / (\d+) benign predictors \(weighted (\d+)% damaging\)",
        r"\1 predictores dañinos / \2 benignos (ponderado: \3 % dañino)",
    ),
    (r"(\d+) damaging / (\d+) benign predictors", r"\1 predictores dañinos / \2 benignos"),
    (r"\bmodel score (\d\.\d+)", r"puntaje del modelo \1"),
    (
        r"common in gnomAD: grpmax (\S+) >= (\S+) \(protein-level call: (.+?)\); frequent variants can still cause "
        r"recessive or founder disease",
        lambda m: (
            f"frecuente en gnomAD: grpmax {m.group(1)} ≥ {m.group(2)} (predicción a nivel de proteína: "
            f"{PREDICTION_ES.get(m.group(3), m.group(3))}); una variante frecuente aún puede causar una enfermedad "
            "recesiva o de efecto fundador"
        ),
    ),
    (r"^common: BA1 \(grpmax (\S+)\)$", r"común: BA1 (grpmax \1)"),
    (
        r"^common: BA1 \(grpmax (\S+); allele numbers not checked in gnomAD\)$",
        r"común: BA1 (grpmax \1; números de alelos no comprobados en gnomAD)",
    ),
    (
        r"^frequent: BS1 candidate \(grpmax (\S+)\); BA1 not met: the gene's expert panel sets BA1 above (\S+)$",
        r"frecuente: posible BS1 (grpmax \1); BA1 no cumplido: el panel experto del gen fija BA1 por encima de \2",
    ),
    (
        r"^frequent, but on the ClinGen BA1 exception list \(grpmax (\S+)\)$",
        r"frecuente, pero está en la lista de excepciones a BA1 de ClinGen (grpmax \1)",
    ),
    (r"^frequent: BS1 candidate \(grpmax (\S+)\)$", r"frecuente: posible BS1 (grpmax \1)"),
    (r"^rare \(grpmax (\S+)\)$", r"rara (grpmax \1)"),
    (r"^not found in gnomAD$", "no encontrada en gnomAD"),
    (r"^absent from (gnomAD v[\d.]+)$", r"ausente de \1"),
    (
        r"^no gnomAD frequency in Ensembl's data \(not checked in gnomAD\)$",
        "sin frecuencia de gnomAD en los datos de Ensembl (no consultado en gnomAD)",
    ),
    (r"(\d+%) of exomes", r"\1 de los exomas"),
    (
        r"^PM2_Supporting \(absent from (gnomAD v[\d.]+); >= 20 reads in (.+)\)$",
        r"PM2_Supporting (ausente de \1; ≥ 20 lecturas en \2)",
    ),
    (
        r"^not applicable \(absent from (gnomAD v[\d.]+), but few reads there: >= 20 reads in (.+)\)$",
        r"no aplicable (ausente de \1, pero con pocas lecturas: ≥ 20 lecturas en \2)",
    ),
    (
        r"^depends on the disease \(present in gnomAD(, failing gnomAD quality filters)?; compare with its maximum "
        r"credible frequency\)$",
        lambda m: (
            "depende de la enfermedad (presente en gnomAD"
            + (", sin pasar sus filtros de calidad" if m.group(1) else "")
            + "; compárese con su frecuencia máxima creíble)"
        ),
    ),
    (r"^not met \(frequent in gnomAD\)$", "no se cumple (frecuente en gnomAD)"),
    (
        r"^none \(no pathogenic missense variant at this residue in ClinVar\)$",
        "ninguna (ClinVar no tiene variantes missense patogénicas en este residuo)",
    ),
    (
        r"^not applied \(SpliceAI (\S+): the variant may act on splicing\)$",
        r"no aplicada (SpliceAI \1: la variante podría actuar sobre el splicing)",
    ),
    (r"^not checked \(ClinVar turned off\)$", "no consultado (ClinVar desactivado)"),
    (r"^not checked \(offline\)$", "no consultado (sin conexión)"),
    (r"^not checked \(ClinVar did not respond\)$", "no consultado (ClinVar no respondió)"),
    (r"same change by (c\.\S+?):", r"mismo cambio por \1:"),
    (r": Pathogenic, (\d)★", r": patogénica, \1★"),
    (r": Likely pathogenic, (\d)★", r": probablemente patogénica, \1★"),
    (r"; and (\d+) more(?=[);])", r"; y \1 más"),
    (r"^none \((.+), indeterminate\)$", r"ninguna (\1, indeterminado)"),
    # ClinGen expert panel notes (vcep.py): one note per specification, its parts separated by "; ".
    # The patterns never cross a "; " or a " | ", so every note and every part of a cell is translated.
    (
        r"allows (PP3|BP4|PS1|PM5|PM2) up to ([A-Za-z+ ]+?): the tool's (\S+) \((\d) points\) is stronger",
        r"admite \1 hasta \2: el \3 de la herramienta (\4 puntos) es más fuerte",
    ),
    (
        r"uses (PP3|BP4|PS1|PM5|PM2) only at ([^:;|]+): the tool's (\S+) is weaker",
        r"usa \1 solo con fuerza \2: el \3 de la herramienta es más débil",
    ),
    (
        r"has no (\d)-point (PP3|BP4|PS1|PM5|PM2) \(([^)]+)\): the tool's (\S+) is not one of them",
        r"no tiene \2 de \1 puntos (\3): el \4 de la herramienta no es ninguna de ellas",
    ),
    (r"does not use (PP3|BP4) for missense variants \(“", r"no usa \1 en variantes missense (“"),
    (r"does not use (PP3|BP4|PS1|PM5|PM2)(?=; | \(http|$)", r"no usa \1"),
    (
        r"it sets its own predictor and thresholds for PP3/BP4 \(the tool's PP3/BP4 is its own calibration\)",
        r"fija su propio predictor y sus umbrales para PP3/BP4 (el PP3/BP4 de la herramienta es su propia calibración)",
    ),
    (r"its rules give neither PP3 nor BP4 \(REVEL", r"sus reglas no dan ni PP3 ni BP4 (REVEL"),
    (r"its rules give ", r"sus reglas dan "),
    (r"no PP3/BP4 for missense variants", r"sin PP3/BP4 en variantes missense"),
    (r"PP3/BP4 not evaluated \(no REVEL score\)", r"PP3/BP4 no evaluado (sin puntuación REVEL)"),
    (r"PM2 not assessed \(no gnomAD data\)", "PM2 no evaluado (sin datos de gnomAD)"),
    (
        r"PM2 not met \(absent from gnomAD, but the panel requires 20x coverage\)",
        "PM2 no cumplido (ausente de gnomAD, pero el panel exige una cobertura de 20x)",
    ),
    (r"\bPM2 not met\b", r"PM2 no cumplido"),
    (
        r"\bPM2 met \(absent from gnomAD, but few reads there: check\)",
        r"PM2 cumplido (ausente de gnomAD, pero con pocas lecturas: compruébelo)",
    ),
    (r"\b(PM2|BA1) met\b", r"\1 cumplido"),
    (
        r"BA1 only (>=?) ([\d.e+-]+%): grpmax ([\d.e+-]+%) does not meet it",
        r"BA1 solo \1 \2: el grpmax \3 no lo cumple",
    ),
    (
        r"grpmax ([\d.e+-]+%) is at or above its BA1 threshold \((>=?) ([\d.e+-]+%)\): check its other BA1 "
        r"conditions",
        r"el grpmax \1 alcanza o supera su umbral de BA1 (\2 \3): compruebe sus otras condiciones de BA1",
    ),
    (r" \(it applies it to the filtering allele frequency\)", r" (lo aplica a la frecuencia alélica de filtrado)"),
    (r"its own BA1 rule: “", r"su propia regla de BA1: “"),
    (r"(^|; |: )on (PP3|BP4|PS1|PM5|PM2): “", r"\1sobre \2: “"),
    (
        r"(^|; |: )on (PP3|BP4|PS1|PM5|PM2): (\d+) more conditions in the specification",
        r"\1sobre \2: \3 condiciones más en la especificación",
    ),
    (
        r"has no (PS1|PM5) below (\w+): the tool lowered (?:PS1|PM5) because the ClinVar reference is likely "
        r"pathogenic or has one review star, and the panel judges the reference under its own specification",
        r"no tiene \1 por debajo de \2: la herramienta rebajó \1 porque la referencia de ClinVar es probablemente "
        r"patogénica o tiene una estrella de revisión, y el panel juzga la referencia según su propia especificación",
    ),
    (r"uses (PP3|BP4) only at ([^:;|]+?) \(its own thresholds\)", r"usa \1 solo con fuerza \2 (sus propios umbrales)"),
    (
        r"grpmax allele frequency ([\d.e+-]+%) is above its BA1 threshold \((>=?) ([\d.e+-]+%)\), which it applies "
        r"to the filtering allele frequency \(not computed here\): check the filtering allele frequency and its other "
        r"BA1 conditions before BA1",
        r"la frecuencia alélica grpmax \1 supera su umbral de BA1 (\2 \3), que aplica a la frecuencia alélica de "
        r"filtrado (no calculada aquí): compruebe la frecuencia de filtrado y sus otras condiciones de BA1 antes de "
        r"aplicar BA1",
    ),
    (r" \(its other BA1 conditions: “", r" (sus otras condiciones de BA1: “"),
    (
        r"on ClinGen's BA1 exception list \(Ghosh et al\. 2018\): BA1 does not apply, whatever the panel's "
        r"threshold",
        "en la lista de excepciones de BA1 de ClinGen (Ghosh et al. 2018): BA1 no se aplica, sea cual sea el umbral "
        "del panel",
    ),
    (
        r"on the panel's BA1 exclusion list \(well-known pathogenic variants\): BA1 does not apply",
        "en la lista de exclusión de BA1 del panel (variantes patogénicas conocidas): BA1 no se aplica",
    ),
    # inheritance in a specification's name (notes) or description (… autosomal dominant (https://…))
    (r"\b(?<=[(,] )autosomal dominant(?= ?[;)]| or| \(http)", "autosómica dominante"),
    (r"\b(?<=[(,] )autosomal recessive(?=[;)]| \(http)", "autosómica recesiva"),
    (r"(?<=or )autosomal recessive(?=[;)]| \(http)", "autosómica recesiva"),
    (r"(?<=[(,] )X-linked(?=[;)]| or| \(http)", "ligada al X"),
    (r"(?<=or )X-linked(?=[;)]| \(http)", "ligada al X"),
    (r"(?<=[(,] )semidominant(?=[;)]| or| \(http)", "semidominante"),
    (
        r"; no SpliceAI score: check that neither variant acts on splicing\)",
        r"; sin puntuación de SpliceAI: compruebe que ninguna de las dos variantes actúe sobre el splicing)",
    ),
    (
        r"^none \((.+); (BP4\S*) withheld: SpliceAI (\S+) does not rule out a splicing effect\)$",
        r"ninguna (\1; \2 retenido: SpliceAI \3 no descarta un efecto sobre el splicing)",
    ),
    (
        r"^none \((.+); (BP4\S*) withheld: no SpliceAI score to rule out a splicing effect\)$",
        r"ninguna (\1; \2 retenido: no hay puntuación de SpliceAI para descartar un efecto sobre el splicing)",
    ),
    (
        r"^none \(only single-star likely pathogenic variants at this residue: ",
        r"ninguna (en este residuo solo hay variantes probablemente patogénicas de una estrella: ",
    ),
    (
        r"^frequent: BS1 candidate \(grpmax (\S+)\); BA1 not met: (\S+) at most in the gnomAD groups with >= 2,000 "
        r"alleles$",
        r"frecuente: posible BS1 (grpmax \1); no cumple BA1: como máximo \2 en los grupos de gnomAD con ≥ 2.000 alelos",
    ),
    (
        r"PP3 \((\d) points\) and (PS1\S*|PM5\S*) may partly reflect the same residue-level information \(Pejaver et "
        r"al\. 2022\): weigh their combination with care; if PM1 is added, PP3 \+ PM1 should not exceed 4 points",
        r"PP3 (\1 puntos) y \2 pueden reflejar en parte la misma información del residuo (Pejaver et al. 2022): "
        r"valore su combinación con cuidado; si se añade PM1, PP3 + PM1 no debe superar 4 puntos",
    ),
    (
        r"frameshift: downstream sequence changes, not reconstructed",
        "cambio del marco de lectura: la secuencia posterior cambia y no se reconstruye",
    ),
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
    (
        r"cysteine lost \(possible disulfide/metal site\)",
        "se pierde una cisteína (posible puente disulfuro o sitio metálico)",
    ),
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
    (
        r"UniProt canonical isoform differs from transcript; features mapped by local sequence match",
        "La isoforma canónica de UniProt difiere del transcrito; las anotaciones se asignaron por coincidencia local",
    ),
    (
        r"UniProt canonical numbering differs from this transcript; UniProt/AlphaFold features skipped",
        "La numeración de UniProt difiere de este transcrito; se omitieron UniProt/AlphaFold",
    ),
    (r"^Gene\(s\) not found in Ensembl: ", "Gen(es) no encontrados en Ensembl: "),
    (r"unknown gene or transcript symbol", "gen o transcrito desconocido"),
    (
        r"^(.+?) did not respond for (\d+) item\(s\) \((.*)\)\. Affected rows are marked in 'notes'\. Run the analysis "
        r"again to retry: everything already obtained is cached\.$",
        r"\1 no respondió para \2 elemento(s) (\3). Las filas afectadas están marcadas en «notas». Vuelve a ejecutar "
        r"el análisis para reintentar: lo ya obtenido está en caché.",
    ),
    (
        r"^(.+?) did not respond \((.*)\)\. Affected rows are marked in 'notes'\. Run the analysis again to retry: "
        r"everything already obtained is cached\.$",
        r"\1 no respondió (\2). Las filas afectadas están marcadas en «notas». Vuelve a ejecutar el análisis para "
        r"reintentar: lo ya obtenido está en caché.",
    ),
    (
        r"predictors unavailable \(service error; run again to retry\)",
        "predictores no disponibles (fallo del servicio; vuelve a ejecutar para reintentar)",
    ),
    (r"ClinVar unavailable \(service error\)", "ClinVar no disponible (fallo del servicio)"),
    (r"UniProt unavailable \(service error\)", "UniProt no disponible (fallo del servicio)"),
    (r"AlphaFold unavailable \(service error\)", "AlphaFold no disponible (fallo del servicio)"),
    (r"protein sequence unavailable \(service error\)", "secuencia de proteína no disponible (fallo del servicio)"),
    (
        r"Gene coordinates unavailable: the gene filter is applied after annotation \(slower for large files\)\.",
        "Coordenadas de genes no disponibles: el filtro de genes se aplica después de anotar "
        "(más lento en archivos grandes).",
    ),
    (r"Genome build not verified \(offline\); assuming (\S+)\.", r"Genoma no verificado (sin conexión); se asume \1."),
    (
        r"No internet connection: using local VEP only; ClinVar, UniProt and AlphaFold are skipped\.",
        "Sin conexión a internet: solo se usa VEP local; se omiten ClinVar, UniProt y AlphaFold.",
    ),
    (
        r"Ensembl cannot be reached: only answers already in the cache are used",
        "No se puede acceder a Ensembl: solo se usan las respuestas que ya están en la caché",
    ),
    (
        r"Cannot reach the Ensembl REST API \(rest\.ensembl\.org\)\. Check the internet connection, or use the local "
        r"VEP mode for VCF/genomic input\.",
        "No hay conexión con la API de Ensembl (rest.ensembl.org). Revisa la conexión a internet o usa el modo "
        "VEP local para VCF o coordenadas genómicas.",
    ),
    (
        r"local VEP needs genomic coordinates \(VCF or chrom:pos ref>alt\); HGVS and rsIDs need the online mode",
        "VEP local necesita coordenadas genómicas (VCF o chrom:pos ref>alt); HGVS y rsIDs requieren el modo en línea",
    ),
    (r"could not run local VEP", "no se pudo ejecutar VEP local"),
    (r"local VEP exited with code", "VEP local terminó con el código"),
    (
        r"VEP could not map this variant \(check the assembly and the alleles\)",
        "VEP no pudo ubicar esta variante (revisa el genoma y los alelos)",
    ),
    (
        r"VEP could not map this variant \(check notation/assembly\)",
        "VEP no pudo ubicar esta variante (revisa la notación y el genoma)",
    ),
    (
        r"reference allele does not match the genome/transcript",
        "el alelo de referencia no coincide con el genoma/transcrito",
    ),
    (
        r"needs a gene, e\.g\. 'BRAF V600E' \(or set the default gene\)",
        "falta el gen, p. ej. «BRAF V600E» (o define el gen por defecto)",
    ),
    (r"unrecognised variant format", "formato de variante no reconocido"),
    (
        r"REF alleles match (\S+) \((\d+%)\) better than (\S+) \((\d+%)\); using (\S+)\.",
        r"Los alelos REF coinciden mejor con \1 (\2) que con \3 (\4); se usa \5.",
    ),
    (
        r"Only (\d+%) of checked REF alleles match (\S+)(.*?)\. Check the genome build and coordinates\.",
        r"Solo el \1 de los alelos REF revisados coincide con \2\3. Revisa el genoma y las coordenadas.",
    ),
    (
        r"Multi-sample VCF \((\d+) samples\): carriers listed per variant\. Choose one sample to get its zygosity\.",
        r"VCF multimuestra (\1 muestras): se listan los portadores por variante. "
        r"Elige una muestra para ver su cigosidad.",
    ),
    (r"could not parse", "no se pudo interpretar"),
    (
        r"the amino acid at this position is not the one given \(check the position and the isoform\)",
        "el aminoácido de esa posición no es el indicado (revisa la posición y la isoforma)",
    ),
    (
        r"position outside the protein \(check the position and the isoform\)",
        "posición fuera de la proteína (revisa la posición y la isoforma)",
    ),
    (
        r"a frameshift written as a protein change cannot be traced to one DNA change: give the c\. change",
        "un desplazamiento del marco escrito como cambio de proteína no lleva a un único cambio de ADN: "
        "indica el cambio c.",
    ),
    (
        r"several DNA changes give this protein change and Ensembl did not resolve them: give the c\. change",
        "varios cambios de ADN producen este cambio de proteína y Ensembl no los resolvió: indica el cambio c.",
    ),
    (
        r"an insertion written as a protein change cannot be traced to one DNA change: give the c\. change",
        "una inserción escrita como cambio de proteína no lleva a un único cambio de ADN: indica el cambio c.",
    ),
    (
        r"Ensembl places this protein change outside the coding sequence of the gene's main transcript "
        r"\(check the position and the isoform\)",
        "Ensembl sitúa este cambio de proteína fuera de la secuencia codificante del transcrito principal del gen "
        "(revisa la posición y la isoforma)",
    ),
    (
        r"(\d+) nucleotide changes give this protein change; each is reported",
        r"\1 cambios de nucleótido producen este cambio de proteína; se informa cada uno",
    ),
    (
        r"protein notation: other DNA changes can give the same (\S+); this is the one Ensembl maps \(([^)]*)\)\. "
        r"Give the c\. or genomic change to analyse a specific one",
        r"notación de proteína: otros cambios de ADN pueden dar el mismo \1; este es el que asigna Ensembl (\2). "
        r"Indica el cambio c. o genómico para analizar uno concreto",
    ),
    (
        r"protein change given as (\S+) exists only on another isoform; on the reported transcript this DNA "
        r"change is (\S+): check the input",
        r"el cambio de proteína indicado, \1, solo existe en otra isoforma; en el transcrito informado este cambio "
        r"de ADN es \2: revisa la entrada",
    ),
    (
        r"protein change given as (\S+) \(another isoform's numbering\); on the reported transcript it is (\S+)",
        r"cambio de proteína indicado como \1 (numeración de otra isoforma); en el transcrito informado es \2",
    ),
    (
        r"not checked \(no MANE Select transcript, on which ClinVar names its variants\)",
        "no comprobado (sin transcrito MANE Select, en el que ClinVar nombra sus variantes)",
    ),
    (r"^Oncogenic$", "Oncogénica"),
    (r"^Likely oncogenic$", "Probablemente oncogénica"),
    (r"^VUS$", "Significado incierto (VUS)"),
    (r"^Benign$", "Benigna"),
    (
        r"hotspot (\S+): (\d+) tumours, (\d+) with (\S+); cancerhotspots\.org",
        r"hotspot \1: \2 tumores, \3 con \4; cancerhotspots.org",
    ),
    (
        r"in-frame indel hotspot (\d+-\d+): (\d+) tumours, (\d+) with (\S+)",
        r"hotspot de indels en fase \1: \2 tumores, \3 con \4",
    ),
    (r"same change by (\S+): ", r"mismo cambio por \1: "),
    (r"the MissensExplorer model calls it damaging", "el modelo de MissensExplorer la llama dañina"),
    (r"the MissensExplorer model calls it benign", "el modelo de MissensExplorer la llama benigna"),
    (r"calibrated computational evidence reaches (PP3|BP4)", r"la evidencia computacional calibrada alcanza \1"),
    (
        r"in an activating gene, a damage prediction does not tell activation from inactivation",
        "en un gen activador, una predicción de daño no distingue activación de inactivación",
    ),
    (r"null variant in a loss-of-function gene", "variante nula en un gen de pérdida de función"),
    (
        r"truncating variant expected to undergo nonsense-mediated decay",
        "variante truncante que se espera degradada por NMD",
    ),
    (
        r"the truncation removes (\d+)% of the protein \(ClinGen LoF guidance: > 10% strong\)",
        r"el truncamiento elimina el \1 % de la proteína (guía LoF de ClinGen: > 10 %, fuerte)",
    ),
    (r"the share of the protein removed is unknown", "se desconoce qué parte de la proteína se elimina"),
    (
        r"start codon lost: another in-frame start codon may be used \(ClinGen LoF guidance\)",
        "codón de inicio perdido: puede usarse otro codón de inicio en fase (guía LoF de ClinGen)",
    ),
    (
        r"canonical splice site; whether skipping the exon keeps the frame is unknown",
        "sitio de empalme canónico; se desconoce si saltar el exón mantiene la fase",
    ),
    (
        r"canonical splice site; skipping exon (\d+) would shift the reading frame",
        r"sitio de empalme canónico; saltar el exón \1 desplazaría la fase de lectura",
    ),
    (
        r"canonical splice site; skipping exon (\d+) would remove (\d+)% of the coding sequence in frame",
        r"sitio de empalme canónico; saltar el exón \1 eliminaría en fase el \2 % de la secuencia codificante",
    ),
    (r"OVS1 not applied \(not a loss-of-function gene: ", "OVS1 no aplicado (no es un gen de pérdida de función: "),
    (
        r"CIViC oncogenic evidence from clinical observations \(not functional studies\) not used",
        "evidencia oncogénica de CIViC de observaciones clínicas (no estudios funcionales) no usada",
    ),
    (r"in-frame change in a cancer gene \(", "cambio en fase en un gen del cáncer ("),
    (r"stop lost in a loss-of-function gene", "codón de parada perdido en un gen de pérdida de función"),
    (
        r"\b(activating|loss-of-function|ambiguous) \(UniProt: (tumour suppressor and proto-oncogene|tumour suppressor"
        r"|proto-oncogene)\)",
        lambda m: (
            {"activating": "activador", "loss-of-function": "de pérdida de función", "ambiguous": "ambiguo"}[m[1]]
            + " (UniProt: "
            + {
                "tumour suppressor": "supresor tumoral",
                "proto-oncogene": "protooncogén",
                "tumour suppressor and proto-oncogene": "supresor tumoral y protooncogén",
            }[m[2]]
            + ")"
        ),
    ),
    (
        r"\b(activating|loss-of-function|ambiguous) \(IntOGen: (\d+)% loss-of-function in (\d+) samples of cohorts "
        r"with a role, (\d+) cohorts?\)",
        lambda m: (
            {"activating": "activador", "loss-of-function": "de pérdida de función", "ambiguous": "ambiguo"}[m[1]]
            + f" (IntOGen: {m[2]} % de pérdida de función en {m[3]} muestras de cohortes con papel, {m[4]} cohortes)"
        ),
    ),
    (
        r"\b(activating|loss-of-function|ambiguous) \(IntOGen: (\d+) of (\d+) mutated samples in cohorts of "
        r"ambiguous role\)",
        lambda m: (
            {"activating": "activador", "loss-of-function": "de pérdida de función", "ambiguous": "ambiguo"}[m[1]]
            + f" (IntOGen: {m[2]} de {m[3]} muestras mutadas en cohortes de papel ambiguo)"
        ),
    ),
    (
        r"\b(activating|loss-of-function|ambiguous) \(IntOGen: (\d+)% loss-of-function in (\d+) samples, (\d+) "
        r"cohorts?; UniProt: (tumour suppressor|proto-oncogene): they disagree\)",
        lambda m: (
            {"activating": "activador", "loss-of-function": "de pérdida de función", "ambiguous": "ambiguo"}[m[1]]
            + f" (IntOGen: {m[2]} % de pérdida de función en {m[3]} muestras, {m[4]} cohortes; UniProt: "
            + {"tumour suppressor": "supresor tumoral", "proto-oncogene": "protooncogén"}[m[5]]
            + ": no coinciden)"
        ),
    ),
    (
        r"UniProt keywords unavailable: gene roles from IntOGen only",
        "palabras clave de UniProt no disponibles: papeles de los genes solo de IntOGen",
    ),
    (r"the standard's rule for OS1 would give (-?\d+) points", r"la regla de la norma para OS1 daría \1 puntos"),
    (
        r"(\S+) is not a known cancer gene \(UniProt, IntOGen\)",
        r"\1 no es un gen del cáncer conocido (UniProt, IntOGen)",
    ),
    (
        r"not assessed automatically: functional studies beyond CIViC, OM1 functional domain, OP2, tumour type",
        "no evaluados automáticamente: estudios funcionales fuera de CIViC, OM1 dominio funcional, OP2, tipo de tumor",
    ),
    (r"OS2/SBS2 not checked \(CIViC unavailable\)", "OS2/SBS2 no comprobados (CIViC no disponible)"),
    (
        r"CIViC functional evidence rated below (\d) stars not used",
        r"evidencia funcional de CIViC con menos de \1 estrellas no usada",
    ),
    (
        r"OS2/SBS2 not applied \(CIViC functional evidence conflicts\)",
        "OS2/SBS2 no aplicados (la evidencia funcional de CIViC es contradictoria)",
    ),
    (
        r"((?:OS3|OM3|OP3|OS2|OM4(?:_Supporting)?)(?:, (?:OS3|OM3|OP3|OS2|OM4(?:_Supporting)?))*) not counted "
        r"\((OS1\S*) uses this amino-acid change\)",
        lambda m: f"{m[1]} no contado{'s' if ',' in m[1] else ''} ({m[2]} usa este cambio de aminoácido)",
    ),
    (
        r"(OS1\S*) not counted \((.+?); ([A-Za-z0-9_ +]+) for this change (is stronger|counts as much and is the "
        r"standard's)\)",
        lambda m: (
            f"{m[1]} no contado ({m[2]}; {m[3]} para este cambio "
            + ("es más fuerte" if m[4] == "is stronger" else "cuenta lo mismo y es el de la norma")
            + ")"
        ),
    ),
    (
        r"OP4 not checked \(absence from gnomAD not confirmed\)",
        "OP4 no comprobado (ausencia en gnomAD no confirmada)",
    ),
    (
        r"OP4 not applied \(absent from gnomAD, but few reads there\)",
        "OP4 no aplicado (ausente de gnomAD, pero con pocas lecturas)",
    ),
    (
        r"\b(activating|loss-of-function|ambiguous) \(IntOGen: (\d+) mutated samples in (\d+) cohorts? with a "
        r"role: too few to decide\)",
        lambda m: (
            {"activating": "activador", "loss-of-function": "de pérdida de función", "ambiguous": "ambiguo"}[m[1]]
            + f" (IntOGen: {m[2]} muestras mutadas en {m[3]} cohorte{'s' if m[3] != '1' else ''} con papel: "
            + "demasiado pocas para decidir)"
        ),
    ),
    (r"\bDominant Negative\b", "dominante negativa"),
    (r"\bGain of Function\b", "ganancia de función"),
    (r"\bLoss of Function\b", "pérdida de función"),
    (r"\bNeomorphic\b", "neomórfica"),
    (r"\bUnaltered Function\b", "función inalterada"),
    (r"\boncogenicity supports\b", "apoya la oncogenicidad"),
    (r"\boncogenicity does not support\b", "no apoya la oncogenicidad"),
    (r"(\d) stars, EID", r"\1 estrellas, EID"),
    (r"(\d+) CIViC items \(level ", r"\1 evidencias en CIViC (nivel "),
    (
        r"; best: level ([A-E]) (predictive|diagnostic|prognostic)",
        lambda m: (
            f"; la mejor: nivel {m[1]}, "
            + {"predictive": "predictiva", "diagnostic": "diagnóstica", "prognostic": "pronóstica"}[m[2]]
        ),
    ),
    (r"Sensitivity/Response", "sensibilidad/respuesta"),
    (r"Reduced Sensitivity", "sensibilidad reducida"),
    (r"\bResistance\b", "resistencia"),
    (r"Poor Outcome", "mal pronóstico"),
    (r"Better Outcome", "mejor pronóstico"),
    (r" in ([^()]+) \(EID(\d+)\)$", r" en \1 (EID\2)"),
    (r"OP4 \(absent from (gnomAD[^)]*)\)", r"OP4 (ausente de \1)"),
    (r"gnomAD frequency (\S+) < 1e-5", r"frecuencia en gnomAD \1 < 1e-5"),
    (r"gnomAD grpmax (\S+) > (\d+)%", r"grpmax de gnomAD \1 > \2 %"),
    (
        r"OS1/OM4 not checked \(ClinVar oncogenicity search unavailable\)",
        "OS1/OM4 no comprobados (búsqueda de oncogenicidad en ClinVar no disponible)",
    ),
    (
        r"OM4 not applied \(this change is milder than the oncogenic ones at the residue\)",
        "OM4 no aplicado (este cambio es más leve que los oncogénicos del residuo)",
    ),
    (r"\bnot a known cancer gene\)", "no es un gen del cáncer conocido)"),
    (
        r"OM2 not applied \(not a known cancer gene: UniProt, IntOGen\)",
        "OM2 no aplicado (no es un gen del cáncer conocido: UniProt, IntOGen)",
    ),
    (
        r"(OM4\S*) not used \((p\.\S+): (.+?); the hotspot criterion OM3 is stronger\)",
        r"\1 no usado (\2: \3; el criterio de hotspot OM3 es más fuerte)",
    ),
    (
        r"not assessed \(cancerhotspots\.org or IntOGen unavailable\)",
        "no evaluada (cancerhotspots.org o IntOGen no disponibles)",
    ),
]

AA_CLASS_ES: dict[str, str] = {
    "hydrophobic": "hidrofóbico",
    "aromatic": "aromático",
    "polar": "polar",
    "positive": "positivo",
    "negative": "negativo",
    "glycine": "glicina",
    "proline": "prolina",
}

AA_NAMES_ES: dict[str, str] = {
    "A": "Alanina",
    "R": "Arginina",
    "N": "Asparagina",
    "D": "Aspartato",
    "C": "Cisteína",
    "Q": "Glutamina",
    "E": "Glutamato",
    "G": "Glicina",
    "H": "Histidina",
    "I": "Isoleucina",
    "L": "Leucina",
    "K": "Lisina",
    "M": "Metionina",
    "F": "Fenilalanina",
    "P": "Prolina",
    "S": "Serina",
    "T": "Treonina",
    "W": "Triptófano",
    "Y": "Tirosina",
    "V": "Valina",
    "*": "Stop",
}
_PHRASES = [(re.compile(a), b) for a, b in PHRASES_ES]


STAGES_ES: dict[str, str] = {
    "Reading input": "Leyendo la entrada",
    "Checking genome build": "Comprobando el genoma",
    "Annotating variants (VEP)": "Anotando variantes (VEP)",
    "Annotating variants (local VEP)": "Anotando variantes (VEP local)",
    "Retrying without predictors": "Reintentando sin predictores",
    "Annotating + predictors (VEP)": "Anotando y puntuando (VEP)",
    "Annotating variants (VEP, fast pass)": "Anotando variantes (VEP, pasada rápida)",
    "Scoring coding variants (predictors)": "Puntuando variantes codificantes",
    "ClinVar classification": "Clasificación ClinVar",
    "Oncogenicity (hotspots, IntOGen, CIViC)": "Oncogenicidad (hotspots, IntOGen, CIViC)",
    "Protein context (UniProt, AlphaFold)": "Contexto proteico (UniProt, AlphaFold)",
    "Done": "Listo",
}


CALLS_ES: dict[str, str] = {
    "likely_pathogenic": "probablemente patogénica",
    "likely_benign": "probablemente benigna",
    "ambiguous": "ambigua",
    "pathogenic": "patogénica",
    "benign": "benigna",
    "uncertain": "incierta",
    "deleterious": "deletérea",
    "deleterious_low_confidence": "deletérea (baja confianza)",
    "tolerated": "tolerada",
    "tolerated_low_confidence": "tolerada (baja confianza)",
    "probably_damaging": "probablemente dañina",
    "possibly_damaging": "posiblemente dañina",
    "unknown": "desconocida",
    "likely splice-altering": "probable alteración del splicing",
    "possible": "posible",
}


def tr_call(call: str | None) -> str:
    """Translate a predictor's call (e.g. ``likely_pathogenic``, ``deleterious``); ``""`` when there is none."""
    if not call:
        return ""
    if LANG == "es":
        return CALLS_ES.get(call.lower(), call.replace("_", " "))
    return call.replace("_", " ")


def tr_stage(stage: str) -> str:
    """Translate the name of a pipeline stage shown under the progress bar."""
    return STAGES_ES.get(stage, stage) if LANG == "es" else stage


def set_language(lang: str) -> None:
    """Switch the UI language: Spanish for any code starting with ``es`` (e.g. ``es_CO``), else English."""
    global LANG
    LANG = "es" if lang.lower().startswith("es") else "en"


def t(key: str, **kw: object) -> str:
    """Return the UI string ``key`` in the current language, formatted with ``kw``.

    Unknown keys are returned unchanged (formatted with ``kw`` when given).
    """
    pair = S.get(key)
    if pair is None:
        return key.format(**kw) if kw else key
    s = pair[0] if LANG == "es" else pair[1]
    return s.format(**kw) if kw else s


def tr_prediction(label: str | None) -> str:
    """Translate a prediction label from ``predict.classify()``; ``""`` when there is none."""
    if not label:
        return ""
    return PREDICTION_ES.get(label, label) if LANG == "es" else label


def aa_name(code: str) -> str:
    """Return the name, in the current language, of a one-letter amino-acid code (``*`` is the stop codon)."""
    from ..aa import AA_NAMES

    return (AA_NAMES_ES if LANG == "es" else AA_NAMES).get(code, code)


def aa_class(name: str | None) -> str:
    """Translate an amino-acid class (``hydrophobic``, ``polar``, ...); ``""`` when there is none."""
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

"use client";

import { ChangeEvent, FormEvent, useEffect, useMemo, useState } from "react";
import { AnalysisResult, buildAnalysisResult, isAnalysisResult, Variant, VariantPrediction } from "../lib/analysis";

type Result = VariantPrediction;
type View = "explore" | "effects" | "prioritize" | "vcf" | "genes" | "tissues" | "statistics" | "evidence" | "network" | "ask" | "reports" | "reproducibility" | "region";
type Job = { job_id: string; access_token?: string; created_at?: string; status: string; requested: number; completed?: number; results?: Result[]; failures?: { variant: Variant; error: string }[]; error?: string };
type HistoryEntry = { analysis_id: string; created_at: string; analysis_type: string; status: string; manifest: Record<string, unknown>; results: unknown };
type EvidenceScoreRow = { variant: string; avi?: number; gwas_score?: number; eqtl_p?: number; gnomad_af?: number };

const apiBase = process.env.NEXT_PUBLIC_API_URL ?? "http://localhost:8000";
const batchLimit = Number(process.env.NEXT_PUBLIC_MAX_BATCH_SIZE ?? 5000);
const historyStorageKey = "aga-local-analysis-history-v1";
const sections: { id: View; label: string }[] = [
  { id: "explore", label: "Explore" }, { id: "effects", label: "Effects" }, { id: "prioritize", label: "Prioritize" },
  { id: "vcf", label: "VCF analysis" }, { id: "genes", label: "Variant → gene" }, { id: "tissues", label: "Tissues" },
  { id: "statistics", label: "Statistics" }, { id: "evidence", label: "Evidence" }, { id: "network", label: "Network" },
  { id: "ask", label: "Ask AGA" }, { id: "reports", label: "Reports" }, { id: "reproducibility", label: "History" }, { id: "region", label: "Region" },
];

function parseVariant(input: string): Variant | null {
  const match = input.trim().match(/^(chr)?([0-9]+|X|Y|M|MT):([0-9]+)\s+([ACGT])>([ACGT])$/i);
  if (!match || Number(match[3]) < 1 || match[4].toUpperCase() === match[5].toUpperCase()) return null;
  return { chromosome: `chr${match[2].toUpperCase()}`, position: Number(match[3]), reference: match[4].toUpperCase(), alternate: match[5].toUpperCase(), genome_assembly: "GRCh38" };
}
function variantKey(variant?: Variant) { return variant ? `${variant.chromosome}:${variant.position} ${variant.reference}>${variant.alternate}` : ""; }
function parseEvidenceScores(text: string): EvidenceScoreRow[] {
  return text.split(/\r?\n/).map((line) => line.trim()).filter((line) => line && !line.startsWith("#"))
    .filter((line) => !/^variant[\t,]/i.test(line)).map((line) => {
      const [variant, gwas_score, eqtl_p, gnomad_af] = line.split(/[\t,]/).map((part) => part.trim());
      const parsed: EvidenceScoreRow = { variant };
      for (const [name, raw] of [["gwas_score", gwas_score], ["eqtl_p", eqtl_p], ["gnomad_af", gnomad_af]] as const) {
        if (raw !== "" && Number.isFinite(Number(raw))) parsed[name] = Number(raw);
      }
      return parsed;
    }).filter((row) => row.variant);
}

async function api<T>(path: string, init?: RequestInit): Promise<T> {
  const response = await fetch(`${apiBase}${path}`, init);
  const data = await response.json();
  if (!response.ok) throw new Error(data.detail ?? `Request failed (${response.status}).`);
  return data as T;
}

function download(name: string, data: unknown, type = "application/json") {
  const blob = new Blob([type === "text/html" ? String(data) : JSON.stringify(data, null, 2)], { type });
  const url = URL.createObjectURL(blob); const anchor = document.createElement("a"); anchor.href = url; anchor.download = name; anchor.click(); URL.revokeObjectURL(url);
}
function downloadHtml(name: string, title: string, payload: unknown) {
  const escaped = JSON.stringify(payload, null, 2).replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;");
  download(name, `<!doctype html><meta charset="utf-8"><title>${title}</title><style>body{font:14px system-ui;max-width:1000px;margin:40px auto;padding:0 20px;color:#183128}pre{white-space:pre-wrap;background:#f4f8f5;padding:20px;border-radius:8px;line-height:1.55}h1{font-weight:600}</style><h1>${title}</h1><p>Research predictions are not clinical interpretations.</p><pre>${escaped}</pre>`, "text/html");
}

function ResultCard({ result }: { result: Result }) {
  if (result.avi_score === undefined) return <div className="empty-state">No AlphaGenome prediction is attached to this record.</div>;
  const importance = result.feature_importance ?? [];
  const max = Math.max(0.0001, ...importance.map((item) => Math.abs(item.value)));
  return <div className="result-panel">
    <div className="result-top"><span className="eyebrow">{result.provider ?? "AlphaGenome Atlas"} · GRCh38</span><span className="result-variant">{result.variant && `${result.variant.chromosome}:${result.variant.position} ${result.variant.reference}>${result.variant.alternate}`}</span></div>
    <div className="score-row"><div><span className="score-label">AlphaGenome Variant Impact (AVI)</span><strong>{result.avi_score.toFixed(4)}</strong></div>{result.avi_quantile != null && <div><span className="score-label">Atlas quantile</span><strong>{(result.avi_quantile * 100).toFixed(1)}<small>%</small></strong></div>}</div>
    {importance.length > 0 && <div className="features"><span className="score-label">Feature attribution · raw value</span>{importance.slice(0, 30).map((item) => <div className="attribution-row" key={item.feature}><div className="attribution-label"><span>{item.feature.replace(/_/g, " ")}</span><b>{item.value.toFixed(4)}</b></div><div className="attribution-track"><i style={{ width: `${Math.max(1, Math.abs(item.value) / max * 100)}%` }} /></div></div>)}</div>}
    <p className="result-note">Prediction for research prioritization. It is not a clinical interpretation.</p>
    <details className="provenance"><summary>Analysis provenance</summary><dl><dt>Genome build</dt><dd>{result.variant?.genome_assembly ?? "GRCh38"}</dd><dt>Provider</dt><dd>{result.provider ?? "AlphaGenome Atlas"}</dd><dt>Queried</dt><dd>{new Date(result.created_at).toLocaleString()}</dd><dt>SDK version</dt><dd>{result.provider_sdk_version ?? "Not available"}</dd><dt>Model version</dt><dd>{result.provider_model_version ?? "Not exposed by the Atlas response"}</dd><dt>Analysis ID</dt><dd>{result.analysis_id}</dd></dl></details>
  </div>;
}

export default function Home() {
  const [view, setView] = useState<View>("explore");
  const [locusWidth, setLocusWidth] = useState(2000);
  const [query, setQuery] = useState(""); const [result, setResult] = useState<Result | null>(null);
  const [message, setMessage] = useState(""); const [busy, setBusy] = useState(false);
  const [vcfVariants, setVcfVariants] = useState<Variant[]>([]); const [consent, setConsent] = useState(false); const [job, setJob] = useState<Job | null>(null);
  const [geneQuery, setGeneQuery] = useState(""); const [rsid, setRsid] = useState(""); const [evidenceData, setEvidenceData] = useState<unknown>(null); const [evidenceQuery, setEvidenceQuery] = useState(""); const [evidenceRetrievedAt, setEvidenceRetrievedAt] = useState<string | null>(null);
  const [scoresA, setScoresA] = useState(""); const [scoresB, setScoresB] = useState(""); const [stats, setStats] = useState<Record<string, unknown> | null>(null); const [statsInput, setStatsInput] = useState<unknown>(null);
  const [genesInput, setGenesInput] = useState(""); const [backgroundInput, setBackgroundInput] = useState(""); const [geneSetsInput, setGeneSetsInput] = useState(""); const [enrichment, setEnrichment] = useState<unknown>(null);
  const [enrichmentInput, setEnrichmentInput] = useState<unknown>(null);
  const [evidenceScoreTsv, setEvidenceScoreTsv] = useState(""); const [priorityWeights, setPriorityWeights] = useState({ avi: 1, gwas_score: 1, eqtl_p: 1, gnomad_af: 1 });
  const [history, setHistory] = useState<HistoryEntry[]>([]); const [question, setQuestion] = useState(""); const [answer, setAnswer] = useState("");
  const [prioritized, setPrioritized] = useState<Result[]>([]);

  useEffect(() => {
    try { setHistory(JSON.parse(window.localStorage.getItem(historyStorageKey) ?? "[]") as HistoryEntry[]); }
    catch { setHistory([]); }
  }, []);
  function saveLocalHistory(entry: HistoryEntry) {
    setHistory((current) => {
      const updated = [entry, ...current.filter((item) => item.analysis_id !== entry.analysis_id)].slice(0, 50);
      try { window.localStorage.setItem(historyStorageKey, JSON.stringify(updated)); } catch { setMessage("Browser storage is full; download your report to keep a copy."); }
      return updated;
    });
  }
  useEffect(() => {
    if (!job || !["queued", "running"].includes(job.status)) return;
    const timer = window.setInterval(() => api<Job>(`/api/v1/jobs/${job.job_id}`, { headers: { Authorization: `Bearer ${job.access_token ?? ""}` } }).then((next) => {
      setJob((current) => current ? { ...next, access_token: current.access_token } : next);
      if (next.status === "complete" || next.status === "failed") {
        setPrioritized(next.results ?? []);
        const analysis = buildAnalysisResult({ analysis_id: next.job_id, created_at: next.created_at ?? new Date().toISOString(), kind: "vcf_batch", predictions: next.results ?? [], requested_variant_count: next.requested, batch: { status: next.status, completed: next.completed ?? 0, failures: next.failures ?? [] } });
        saveLocalHistory({ analysis_id: next.job_id, created_at: next.created_at ?? new Date().toISOString(), analysis_type: "vcf_batch", status: next.status, manifest: { genome_assembly: "GRCh38", variant_count: next.requested }, results: analysis });
      }
    }).catch((e) => setMessage(e.message)), 2500);
    return () => window.clearInterval(timer);
  }, [job]);

  async function submitVariant(event: FormEvent) {
    event.preventDefault(); const variant = parseVariant(query);
    if (!variant) { setMessage("Enter a GRCh38 SNV like chr22:36201698 A>C."); return; }
    setResult(null); setEvidenceData(null); setEvidenceQuery(""); setEvidenceRetrievedAt(null); setStats(null); setStatsInput(null); setEnrichment(null); setEnrichmentInput(null);
    setBusy(true); setMessage("");
    try { const data = await api<Result>("/api/v1/variants/analyze", { method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify(variant) }); setResult(data); const analysis = buildAnalysisResult({ analysis_id: data.analysis_id, created_at: data.created_at, kind: "single_variant", predictions: [data] }); saveLocalHistory({ analysis_id: data.analysis_id, created_at: data.created_at, analysis_type: "single_variant", status: "complete", manifest: { variant, genome_assembly: "GRCh38", provider: data.provider }, results: analysis }); setMessage("Atlas prediction saved in this browser’s local history."); }
    catch (e) { setMessage(e instanceof Error ? e.message : "Could not reach the analysis API."); } finally { setBusy(false); }
  }

  async function readVcf(event: ChangeEvent<HTMLInputElement>) {
    const file = event.target.files?.[0]; if (!file) return;
    const text = await file.text(); const parsed: Variant[] = []; const rejected: string[] = [];
    for (const line of text.split(/\r?\n/)) {
      if (!line || line.startsWith("#")) continue; const fields = line.split("\t");
      if (fields.length < 5 || fields[3].length !== 1 || fields[4].length !== 1 || !"ACGT".includes(fields[3].toUpperCase()) || !"ACGT".includes(fields[4].toUpperCase())) { if (rejected.length < 4) rejected.push(line); continue; }
      const value = parseVariant(`${fields[0]}:${fields[1]} ${fields[3]}>${fields[4]}`); if (value) parsed.push(value); else if (rejected.length < 4) rejected.push(line);
    }
    setVcfVariants(parsed); setJob(null); setConsent(false); setMessage(`${parsed.length.toLocaleString()} SNVs loaded locally. Confirm the file uses GRCh38 before analysis.${rejected.length ? ` ${rejected.length} example non-SNV/invalid rows were skipped.` : ""}`);
  }

  async function launchBatch() {
    if (!consent) { setMessage("Confirm the data-sharing notice before sending variants for analysis."); return; }
    if (vcfVariants.length > batchLimit) { setMessage(`This deployment accepts up to ${batchLimit.toLocaleString()} SNVs per batch.`); return; }
    setBusy(true); setMessage("");
    try { const data = await api<Job>("/api/v1/variants/analyze-batch", { method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify({ variants: vcfVariants }) }); setJob(data); setMessage("Batch submitted. Check progress below."); }
    catch (e) { setMessage(e instanceof Error ? e.message : "Could not submit batch."); } finally { setBusy(false); }
  }

  async function lookup(path: string) { setBusy(true); setEvidenceData(null); setEvidenceQuery(path); setEvidenceRetrievedAt(null); setMessage(""); try { setEvidenceData(await api(path)); setEvidenceRetrievedAt(new Date().toISOString()); } catch (e) { setMessage(e instanceof Error ? e.message : "Evidence lookup failed."); } finally { setBusy(false); } }
  async function compare(event: FormEvent) {
    event.preventDefault(); setStats(null); setStatsInput(null); setMessage("");
    const numbers = (value: string) => value.split(/[\s,]+/).filter(Boolean).map(Number);
    const caseScores = numbers(scoresA); const controlScores = numbers(scoresB);
    if (caseScores.length < 2 || controlScores.length < 2 || [...caseScores, ...controlScores].some((value) => !Number.isFinite(value))) {
      setMessage("Enter at least two finite numeric scores in each group, separated by commas or spaces."); return;
    }
    const request = { case_scores: caseScores, control_scores: controlScores, iterations: 10000, seed: 42 };
    try { setStats(await api("/api/v1/statistics/compare", { method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify(request) })); setStatsInput(request); }
    catch (e) { setMessage(e instanceof Error ? e.message : "Statistical analysis failed."); }
  }
  async function runEnrichment() {
    setEnrichment(null); setEnrichmentInput(null); setMessage("");
    const genes = genesInput.split(/[\s,]+/).filter(Boolean);
    const background = backgroundInput.split(/[\s,]+/).filter(Boolean);
    if (!genes.length || !background.length) { setMessage("Enter selected genes and the background universe used for your analysis."); return; }
    try {
      const geneSets: unknown = JSON.parse(geneSetsInput);
      if (!geneSets || typeof geneSets !== "object" || Array.isArray(geneSets) || Object.keys(geneSets).length === 0 || Object.values(geneSets).some((members) => !Array.isArray(members))) {
        setMessage("Enter gene sets as a JSON object whose values are arrays of gene symbols."); return;
      }
      const request = { genes, background, gene_sets: geneSets };
      setEnrichment(await api("/api/v1/statistics/enrichment", { method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify(request) })); setEnrichmentInput(request);
    } catch (e) { setMessage(e instanceof Error ? e.message : "Enrichment failed; check the gene-set JSON."); }
  }

  const analysisResult = useMemo<AnalysisResult | null>(() => {
    if (!result) return null;
    const evidence: AnalysisResult["evidence"] = evidenceData && typeof evidenceData === "object"
      ? (() => {
          const data = evidenceData as Record<string, unknown>;
          return [{
            provider: typeof data.source === "string" ? data.source : "External evidence provider",
            retrieved_at: evidenceRetrievedAt ?? new Date().toISOString(),
            source_url: typeof data.source_url === "string" ? data.source_url : undefined,
            query: evidenceQuery,
            data: evidenceData,
          }];
        })()
      : [];
    const statistics: AnalysisResult["statistics"] = [
      ...(stats ? [{ kind: "score_group_comparison" as const, created_at: String(stats.created_at ?? result.created_at), input: statsInput, data: stats }] : []),
      ...(enrichment ? [{ kind: "gene_set_enrichment" as const, created_at: new Date().toISOString(), input: enrichmentInput, data: enrichment }] : []),
    ];
    return buildAnalysisResult({
      analysis_id: result.analysis_id,
      created_at: result.created_at,
      kind: "single_variant",
      predictions: [result],
      evidence,
      statistics,
    });
  }, [result, evidenceData, evidenceQuery, evidenceRetrievedAt, stats, statsInput, enrichment, enrichmentInput]);

  useEffect(() => {
    if (!analysisResult || (!analysisResult.evidence.length && !analysisResult.statistics.length)) return;
    saveLocalHistory({
      analysis_id: analysisResult.analysis_id,
      created_at: analysisResult.created_at,
      analysis_type: "single_variant",
      status: "complete",
      manifest: { ...analysisResult.input, provider: analysisResult.alphagenome.provider },
      results: analysisResult,
    });
  // Persist attachments when they are added to the active analysis.
  }, [analysisResult]);
  const batchAnalysisResult = useMemo<AnalysisResult | null>(() => {
    if (!job?.results?.length) return null;
    return buildAnalysisResult({
      analysis_id: job.job_id,
      created_at: job.created_at ?? new Date().toISOString(),
      kind: "vcf_batch",
      predictions: job.results,
      requested_variant_count: job.requested,
      batch: { status: job.status, completed: job.completed ?? 0, failures: job.failures ?? [] },
    });
  }, [job]);
  const effects = useMemo(() => analysisResult?.annotations.effects ?? [], [analysisResult]);
  const tissueSummary = useMemo(() => {
    const groups = new Map<string, { n: number; sum: number; top: number }>();
    for (const effect of effects) if (effect.tissue) { const group = groups.get(effect.tissue) ?? { n: 0, sum: 0, top: 0 }; group.n++; group.sum += effect.raw_score; group.top = Math.max(group.top, Math.abs(effect.raw_score)); groups.set(effect.tissue, group); }
    return [...groups.entries()].map(([tissue, values]) => ({ tissue, ...values, mean: values.sum / values.n })).sort((a, b) => b.top - a.top).slice(0, 40);
  }, [effects]);
  const candidateGenes = useMemo(() => analysisResult?.annotations.genes ?? [], [analysisResult]);
  const networkEdges = useMemo(() => effects.filter((effect) => effect.gene && effect.feature).slice(0, 60), [effects]);
  const nearbyFeatures = useMemo(() => {
    const record = evidenceData && typeof evidenceData === "object" ? evidenceData as Record<string, unknown> : {};
    const features = Array.isArray(record.features) ? record.features as Record<string, unknown>[] : [];
    const position = result?.variant?.position;
    if (!position) return [];
    return features.map((feature) => {
      const start = Number(feature.start); const end = Number(feature.end);
      const distance = Number.isFinite(start) && Number.isFinite(end) ? position < start ? start - position : position > end ? position - end : 0 : Infinity;
      return { feature, distance };
    }).sort((a, b) => a.distance - b.distance).slice(0, 100);
  }, [evidenceData, result]);
  const modalityGroups = useMemo(() => {
    const rows = prioritized.length ? prioritized : result ? [result] : [];
    const counts = new Map<string, number>();
    for (const row of rows) {
      const totals = new Map<string, number>();
      for (const effect of row.effects ?? []) totals.set(effect.modality, (totals.get(effect.modality) ?? 0) + Math.abs(effect.raw_score));
      const dominant = [...totals.entries()].sort((a, b) => b[1] - a[1])[0]?.[0] ?? "No scorer effect returned";
      counts.set(dominant, (counts.get(dominant) ?? 0) + 1);
    }
    return [...counts.entries()].sort((a, b) => b[1] - a[1]);
  }, [prioritized, result]);
  const combinedRanking = useMemo(() => {
    const base = prioritized.length ? prioritized : result ? [result] : [];
    const evidenceRows = new Map(parseEvidenceScores(evidenceScoreTsv).map((row) => [row.variant, row]));
    const components = ["avi", "gwas_score", "eqtl_p", "gnomad_af"] as const;
    const directions = { avi: 1, gwas_score: 1, eqtl_p: -1, gnomad_af: -1 };
    if (base.length < 2) return [];
    const values = base.map((item) => {
      const ext = evidenceRows.get(variantKey(item.variant)) ?? {};
      return { item, row: { ...ext, ...(item.avi_score === undefined ? {} : { avi: item.avi_score }) } as EvidenceScoreRow };
    });
    const percentileMaps = new Map<string, Map<string, number>>();
    for (const component of components) {
      const present = values.filter((value) => Number.isFinite(value.row[component]));
      if (present.length < 2) { percentileMaps.set(component, new Map()); continue; }
      const sorted = present.map((value) => value.row[component] as number).sort((a, b) => a - b);
      percentileMaps.set(component, new Map(present.map(({ item, row }) => {
        const value = row[component] as number; const lo = sorted.findIndex((x) => x === value); const hi = sorted.lastIndexOf(value);
        const percentile = sorted.length < 2 ? 1 : (lo + hi) / (2 * (sorted.length - 1));
        return [item.analysis_id, directions[component] > 0 ? percentile : 1 - percentile];
      })));
    }
    return values.map(({ item }) => {
      const available = components.filter((component) => percentileMaps.get(component)?.has(item.analysis_id) && priorityWeights[component] > 0);
      const weightTotal = available.reduce((sum, component) => sum + priorityWeights[component], 0);
      const score = weightTotal ? available.reduce((sum, component) => sum + (percentileMaps.get(component)?.get(item.analysis_id) ?? 0) * priorityWeights[component], 0) / weightTotal : null;
      return { item, score, components: Object.fromEntries(available.map((component) => [component, percentileMaps.get(component)?.get(item.analysis_id)])), external_scores: evidenceRows.get(variantKey(item.variant)) ?? {} };
    }).filter((row) => row.score !== null).sort((a, b) => (b.score ?? 0) - (a.score ?? 0));
  }, [prioritized, result, evidenceScoreTsv, priorityWeights]);

  function answerQuestion() {
    const q = question.toLowerCase();
    if (!result) { setAnswer("Run a variant analysis or VCF batch first. Ask AGA answers from loaded results and does not generate biological claims independently."); return; }
    if (/tissue|cell/.test(q)) setAnswer(tissueSummary.length ? `Top tissues/cell contexts in the AlphaGenome result, ranked by largest absolute raw score:\n${tissueSummary.slice(0, 8).map((x) => `${x.tissue}: max |score| ${x.top.toPrecision(4)} across ${x.n} records`).join("\n")}` : "The loaded AlphaGenome result has no tissue labels to summarize.");
    else if (/gene/.test(q)) setAnswer(candidateGenes.length ? `Genes present in the loaded AlphaGenome effects (not a causal-gene determination):\n${candidateGenes.slice(0, 25).join(", ")}` : "No gene labels were returned in the loaded effect records.");
    else if (/impact|avi|score|variant/.test(q)) setAnswer(`Loaded result: AVI ${result.avi_score?.toPrecision(4) ?? "unavailable"}${result.avi_quantile == null ? "" : `; Atlas quantile ${(result.avi_quantile * 100).toFixed(1)}%`}. The score prioritizes predicted impact and does not establish disease relevance.`);
    else setAnswer("I can summarize the loaded AVI score, genes, or tissue/cell contexts. I only use analysis results already loaded in this app.");
  }

  function restoreHistoryEntry(entry: HistoryEntry) {
    const saved = entry.results;
    if (entry.analysis_type === "single_variant") {
      const prediction = isAnalysisResult(saved) ? saved.alphagenome.predictions[0] : saved as Result;
      if (prediction) {
        setResult(prediction);
        setPrioritized([]);
      }
      if (isAnalysisResult(saved)) {
        const latestEvidence = saved.evidence.at(-1);
        setEvidenceData(latestEvidence?.data ?? null);
        setEvidenceQuery(latestEvidence?.query ?? "");
        setEvidenceRetrievedAt(latestEvidence?.retrieved_at ?? null);
        setStats((saved.statistics.find((item) => item.kind === "score_group_comparison")?.data as Record<string, unknown> | undefined) ?? null);
        setStatsInput(saved.statistics.find((item) => item.kind === "score_group_comparison")?.input ?? null);
        setEnrichment(saved.statistics.find((item) => item.kind === "gene_set_enrichment")?.data ?? null);
        setEnrichmentInput(saved.statistics.find((item) => item.kind === "gene_set_enrichment")?.input ?? null);
      } else {
        setEvidenceData(null); setEvidenceQuery(""); setEvidenceRetrievedAt(null); setStats(null); setStatsInput(null); setEnrichment(null); setEnrichmentInput(null);
      }
    } else if (entry.analysis_type === "vcf_batch") {
      if (isAnalysisResult(saved)) {
        const predictions = saved.alphagenome.predictions;
        setPrioritized(predictions);
        setJob({ job_id: saved.analysis_id, created_at: saved.created_at, status: saved.batch?.status ?? "complete", requested: saved.input.requested_variant_count ?? saved.input.variants.length, completed: saved.batch?.completed ?? predictions.length, results: predictions, failures: saved.batch?.failures as Job["failures"] });
      } else {
        const previousJob = saved as Job;
        setJob(previousJob);
        setPrioritized(previousJob.results ?? []);
      }
    }
    setMessage(`Loaded ${entry.analysis_id} from this browser.`);
  }

  const title = sections.find((item) => item.id === view)?.label ?? "Explore";
  return <main className="shell">
    <header className="topbar"><a className="brand" href="#top" onClick={() => setView("explore")}><span className="brand-mark">A</span><span>AGA</span></a><nav className="main-nav" aria-label="Analysis workflows">{sections.map((section) => <button key={section.id} className={view === section.id ? "active" : ""} onClick={() => { setView(section.id); setMessage(""); }}>{section.label}</button>)}</nav><span className="status"><i /> Research workspace · GRCh38</span></header>
    <section className="hero compact" id="top"><div className="eyebrow"><span /> GENOMIC SIGNAL, MADE EXPLORABLE</div><h1>From genomic prediction<br/><em>to biological insight.</em></h1><p className="intro">Explore predictions, organize evidence, and preserve the analysis trail.</p></section>
    <section className="workspace" aria-live="polite"><div className="workspace-heading"><div><span className="eyebrow">WORKSPACE / {title.toUpperCase()}</span><h2>{title}</h2></div><span className="preview-tag">PREDICTIONS · EVIDENCE · REPRODUCIBILITY</span></div>

      {view === "explore" && <><p className="workspace-copy">Enter a GRCh38 SNV to query AlphaGenome Atlas. Searches are sent to the configured AlphaGenome service and saved to local analysis history.</p><form className="search" onSubmit={submitVariant}><span className="search-icon">⌕</span><input aria-label="Variant" placeholder="chr22:36201698 A>C" value={query} onChange={(e) => setQuery(e.target.value)}/><button disabled={busy}>{busy ? "Querying…" : "Analyze variant"} ↗</button></form>{result && <ResultCard result={result}/>}<p className={`feedback ${message ? "visible" : ""}`} role="status">{message || "Single nucleotide variants on GRCh38 are supported."}</p></>}
      {view === "effects" && <><p className="workspace-copy">AlphaGenome Atlas output separated by its returned scorer, feature, gene, tissue, raw score, and quantile where available.</p>{result ? <><ResultCard result={result}/><div className="result-table-wrap"><table className="result-table"><thead><tr><th>SCORER</th><th>FEATURE</th><th>GENE</th><th>TISSUE / CELL</th><th>RAW SCORE</th><th>QUANTILE</th></tr></thead><tbody>{effects.slice(0, 200).map((e, i) => <tr key={i}><td>{e.modality}</td><td>{e.feature ?? "—"}</td><td>{e.gene ?? "—"}</td><td>{e.tissue ?? "—"}</td><td>{e.raw_score.toPrecision(4)}</td><td>{e.quantile_score?.toPrecision(4) ?? "—"}</td></tr>)}</tbody></table></div></> : <div className="empty-state">Analyze a variant first.</div>}</>}
      {view === "prioritize" && <><p className="workspace-copy">AVI ranking is available by default. For a combined research ranking, paste source scores keyed by the exact variant string. Values are converted to within-dataset percentiles; AVI/GWAS are higher-is-stronger, while eQTL p-value and gnomAD allele frequency are lower-is-rarer/stronger. The displayed weights are user-controlled, and missing components are omitted per variant.</p>{(prioritized.length ? prioritized : result ? [result] : []).length ? <><div className="result-table-wrap"><table className="result-table"><thead><tr><th>RANK</th><th>VARIANT</th><th>AVI</th><th>ATLAS QUANTILE</th><th>ANALYSIS ID</th></tr></thead><tbody>{(prioritized.length ? [...prioritized] : [result!]).sort((a,b)=>(b.avi_score??-Infinity)-(a.avi_score??-Infinity)).map((r,i)=><tr key={r.analysis_id}><td>{i+1}</td><td>{variantKey(r.variant)}</td><td>{r.avi_score?.toPrecision(4) ?? "—"}</td><td>{r.avi_quantile == null ? "—" : `${(r.avi_quantile*100).toFixed(1)}%`}</td><td>{r.analysis_id}</td></tr>)}</tbody></table></div><h3 className="subheading">Combined evidence ranking</h3><p className="result-note">Paste rows as tab-separated values in this order: variant, GWAS score (higher is stronger), GTEx q/p value (lower is stronger), gnomAD AF (lower is rarer). AVI comes from the loaded AlphaGenome result.</p><textarea className="ranking-input" aria-label="External evidence scores" placeholder={`variant\tgwas_score\teqtl_p\tgnomad_af\nchr22:36201698 A>C\t0.8\t0.001\t0.0002`} value={evidenceScoreTsv} onChange={(e)=>setEvidenceScoreTsv(e.target.value)}/><div className="weight-grid">{([["avi","AVI weight"],["gwas_score","GWAS weight"],["eqtl_p","GTEx weight"],["gnomad_af","gnomAD weight"]] as const).map(([key,label])=><label key={key}>{label}<input type="number" min="0" max="10" step="0.1" value={priorityWeights[key]} onChange={(e)=>setPriorityWeights((old)=>({...old,[key]:Number(e.target.value)}))}/></label>)}</div>{combinedRanking.length ? <><div className="result-table-wrap"><table className="result-table"><thead><tr><th>RANK</th><th>VARIANT</th><th>COMPOSITE PERCENTILE</th><th>COMPONENT PERCENTILES</th></tr></thead><tbody>{combinedRanking.map((row,index)=><tr key={row.item.analysis_id}><td>{index+1}</td><td>{variantKey(row.item.variant)}</td><td>{(100*(row.score ?? 0)).toFixed(1)}%</td><td>{Object.entries(row.components).map(([key,value])=>`${key}: ${(100*Number(value)).toFixed(0)}%`).join(" · ")}</td></tr>)}</tbody></table></div><button className="primary-action" onClick={()=>download("aga-combined-prioritization.json",{method:"weighted mean of within-dataset percentiles",directions:{avi:"higher",gwas_score:"higher",eqtl_p:"lower",gnomad_af:"lower"},weights:priorityWeights,results:combinedRanking})}>Download ranking + method</button></> : <div className="empty-state">No combined ranking yet. Load variant results and provide any available external scores.</div>}<h3 className="subheading">Dominant scorer groups</h3><div className="cluster-grid">{modalityGroups.map(([name,count])=><div key={name}><strong>{count}</strong><span>{name}</span></div>)}</div></> : <div className="empty-state">Run a single variant or VCF analysis to create a ranked result list.</div>}</>}
      {view === "vcf" && <><p className="workspace-copy">VCF parsing is local. This prototype accepts uncompressed GRCh38 SNVs only (one base REF and ALT, A/C/G/T). Indels and other rows are skipped.</p><label className="upload-box"><span className="upload-icon">↑</span><b>Choose a VCF file</b><small>File contents are parsed in this browser; nothing is uploaded until you start analysis.</small><input type="file" accept=".vcf,text/plain" onChange={readVcf}/></label>{vcfVariants.length > 0 && <><div className="vcf-summary"><strong>{vcfVariants.length.toLocaleString()}</strong><span>eligible SNVs · limit {batchLimit.toLocaleString()} per batch</span><small>Parsed genotypes are discarded; only normalized coordinates are sent.</small></div><label className="consent"><input type="checkbox" checked={consent} onChange={(e)=>setConsent(e.target.checked)}/><span>I confirm this file uses GRCh38. I understand that variant coordinates will be sent to AlphaGenome Atlas and saved in this browser’s local history for this batch.</span></label><button className="primary-action" onClick={launchBatch} disabled={busy || !consent || !!job && ["queued","running"].includes(job.status)}>Start Atlas batch</button></>}{job && <div className="scope-note"><b>Batch {job.status}</b><span>{job.completed ?? 0} / {job.requested} variants processed · job {job.job_id}</span>{job.error && <span>{job.error}</span>}</div>}{job?.results?.length ? <><div className="dataset-stats"><div><strong>{job.results.length}</strong><span>SUCCEEDED</span></div><div><strong>{job.failures?.length ?? 0}</strong><span>FAILED</span></div></div><div className="heatmap"><div className="heatmap-row heatmap-head"><span>VARIANT</span><span>RELATIVE AVI MAGNITUDE</span><b>AVI</b></div>{[...job.results].sort((a,b)=>(b.avi_score??0)-(a.avi_score??0)).slice(0,100).map((r)=><div className="heatmap-row" key={r.analysis_id}><span>{r.variant && `${r.variant.chromosome}:${r.variant.position} ${r.variant.reference}>${r.variant.alternate}`}</span><i style={{width:`${Math.max(2,Math.abs(r.avi_score??0)/Math.max(1e-10,...job.results!.map((item)=>Math.abs(item.avi_score??0)))*100)}%`}}/><b>{r.avi_score?.toPrecision(3)}</b></div>)}</div><p className="result-note">Bar lengths are scaled to the largest absolute AVI in this batch; the table and reports preserve raw AVI values.</p><button className="primary-action" onClick={()=>download(`aga-${job.job_id}.json`, job)}>Download batch report</button></> : null}</>}
      {view === "genes" && <><p className="workspace-copy">Candidate labels from AlphaGenome are shown beside Ensembl genes and regulatory features near the selected variant. Proximity and model labels are evidence for prioritizing a gene, not proof of causality.</p><div className="ask-form"><input placeholder="Gene symbol, e.g. TP53" value={geneQuery} onChange={(e)=>setGeneQuery(e.target.value)}/><button onClick={()=>lookup(`/api/v1/evidence/gene/${encodeURIComponent(geneQuery)}`)} disabled={busy || !geneQuery}>Look up gene</button></div>{result?.variant && <button className="primary-action" onClick={()=>{const v=result.variant!;lookup(`/api/v1/evidence/region?chromosome=${encodeURIComponent(v.chromosome)}&start=${Math.max(1,v.position-100000)}&end=${v.position+100000}`);}}>Annotate 100 kb around loaded variant</button>}{candidateGenes.length > 0 && <p className="scope-note"><b>Genes in AlphaGenome output</b><span>{candidateGenes.join(", ")}</span></p>}{nearbyFeatures.length > 0 && <div className="result-table-wrap"><table className="result-table"><thead><tr><th>ENSEMBL FEATURE</th><th>TYPE</th><th>START–END</th><th>DISTANCE TO VARIANT</th></tr></thead><tbody>{nearbyFeatures.map(({feature,distance},index)=><tr key={`${String(feature.id)}-${index}`}><td>{String(feature.symbol ?? feature.gene_id ?? feature.id ?? "feature")}</td><td>{String(feature.feature_type ?? feature.biotype ?? "—")}</td><td>{String(feature.start ?? "—")}–{String(feature.end ?? "—")}</td><td>{Number.isFinite(distance) ? `${distance.toLocaleString()} bp` : "—"}</td></tr>)}</tbody></table></div>}{evidenceData && <pre className="answer-panel">{JSON.stringify(evidenceData,null,2)}</pre>}</>}
      {view === "tissues" && <><p className="workspace-copy">Tissue and cell labels are summarized from the returned AlphaGenome effect records. Only raw output scores are aggregated here.</p>{tissueSummary.length ? <div className="result-table-wrap"><table className="result-table"><thead><tr><th>TISSUE / CELL</th><th>EFFECT ROWS</th><th>MEAN RAW SCORE</th><th>MAX ABS RAW SCORE</th></tr></thead><tbody>{tissueSummary.map((t)=><tr key={t.tissue}><td>{t.tissue}</td><td>{t.n}</td><td>{t.mean.toPrecision(4)}</td><td>{t.top.toPrecision(4)}</td></tr>)}</tbody></table></div> : <div className="empty-state">Run a variant query; tissue summaries appear when labels are returned by Atlas.</div>}</>}
      {view === "statistics" && <><p className="workspace-copy">Compare score groups with a permutation test and bootstrap confidence interval, or test supplied gene sets against a stated background. Enter your own observed scores and gene lists; no sample data is preloaded.</p><form onSubmit={compare}><div className="form-grid"><label>Group A scores<textarea placeholder="At least 2 numbers, separated by commas or spaces" value={scoresA} onChange={(e)=>setScoresA(e.target.value)}/></label><label>Group B scores<textarea placeholder="At least 2 numbers, separated by commas or spaces" value={scoresB} onChange={(e)=>setScoresB(e.target.value)}/></label></div><button className="primary-action">Run comparison · 10,000 iterations</button></form>{stats && <pre className="answer-panel">{JSON.stringify(stats,null,2)}</pre>}<div className="form-grid"><label>Selected genes<textarea placeholder="Gene symbols from your selected set" value={genesInput} onChange={(e)=>setGenesInput(e.target.value)}/></label><label>Background universe<textarea placeholder="Gene symbols in the analysis background" value={backgroundInput} onChange={(e)=>setBackgroundInput(e.target.value)}/></label></div><label className="form-grid"><span>Gene sets as JSON: name → gene list<textarea placeholder={'{"My set":["GENE1","GENE2"]}'} value={geneSetsInput} onChange={(e)=>setGeneSetsInput(e.target.value)}/></span></label><button className="primary-action" onClick={runEnrichment}>Run enrichment</button>{enrichment && <pre className="answer-panel">{JSON.stringify(enrichment,null,2)}</pre>}</>}
      {view === "evidence" && <><p className="workspace-copy">Read-only lookups connect Ensembl gene/region/variation, GWAS Catalog, ClinVar, GTEx eQTLs, and gnomAD. Results include source provenance and run only when selected.</p><div className="ask-form"><input placeholder="rsID, e.g. rs12345" value={rsid} onChange={(e)=>setRsid(e.target.value)}/><button onClick={()=>lookup(`/api/v1/evidence/variant/${encodeURIComponent(rsid)}`)} disabled={busy || !rsid}>Ensembl</button><button onClick={()=>lookup(`/api/v1/evidence/gwas/${encodeURIComponent(rsid)}`)} disabled={busy || !rsid}>GWAS</button><button onClick={()=>lookup(`/api/v1/evidence/clinvar/${encodeURIComponent(rsid)}`)} disabled={busy || !rsid}>ClinVar</button><button onClick={()=>lookup(`/api/v1/evidence/gtex/${encodeURIComponent(rsid)}`)} disabled={busy || !rsid}>GTEx eQTL</button></div>{result?.variant && <button className="primary-action" onClick={()=>{const v=result.variant!;lookup(`/api/v1/evidence/gnomad?chromosome=${encodeURIComponent(v.chromosome)}&position=${v.position}&reference=${v.reference}&alternate=${v.alternate}`);}}>gnomAD frequencies for loaded variant</button>}{result?.variant && <button className="primary-action" onClick={()=>{const v=result.variant!;const start=Math.max(1,v.position-500);const end=v.position+500;lookup(`/api/v1/evidence/encode/region?chromosome=${encodeURIComponent(v.chromosome)}&start=${start}&end=${end}`);}}>ENCODE region lookup</button>}{evidenceData && <pre className="answer-panel">{JSON.stringify(evidenceData,null,2)}</pre>}{evidenceData && <p className="scope-note"><b>Lookup provenance</b><span>{typeof (evidenceData as Record<string, unknown>).source === "string" ? String((evidenceData as Record<string, unknown>).source) : "External evidence provider"} · retrieved {evidenceRetrievedAt ? new Date(evidenceRetrievedAt).toLocaleString() : "time unavailable"}</span><span>Request: {evidenceQuery}</span>{typeof (evidenceData as Record<string, unknown>).source_url === "string" && <span>Source URL: {String((evidenceData as Record<string, unknown>).source_url)}</span>}</p>}<p className="result-note">ENCODE region search requires a configured ENCODE-DCC genomic-data-service base URL in backend/.env. gnomAD uses coordinates from the loaded result. ClinVar assertions may conflict and must not be used as clinical advice.</p></>}
      {view === "network" && <><p className="workspace-copy">A lightweight graph view of feature-to-gene labels actually present in the latest AlphaGenome result. Edges show returned co-occurrence and do not establish regulation.</p>{networkEdges.length ? <div className="network-list">{networkEdges.map((e,i)=><div key={i}><span>{e.feature}</span><b>→</b><strong>{e.gene}</strong><small>{e.modality} · {e.raw_score.toPrecision(3)}</small></div>)}</div> : <div className="empty-state">No feature-to-gene labels are loaded for graphing.</div>}</>}
      {view === "ask" && <><p className="workspace-copy">Ask AGA summarizes loaded results using deterministic rules; it does not make external model calls or invent biological explanations.</p><div className="question-chips">{["Summarize the variant impact", "Which genes are listed?", "Which tissues have the strongest scores?"].map((q)=><button key={q} onClick={()=>setQuestion(q)}>{q}</button>)}</div><form className="ask-form" onSubmit={(e)=>{e.preventDefault();answerQuestion();}}><input placeholder="Ask about the loaded result…" value={question} onChange={(e)=>setQuestion(e.target.value)}/><button>Ask</button></form>{answer && <div className="answer-panel">{answer}</div>}</>}
      {view === "reports" && <><p className="workspace-copy">Download portable JSON records or a readable HTML report. Analysis reports use a shared structure for inputs, AlphaGenome output, annotations, evidence, statistics, and provenance.</p><div className="action-row"><button className="primary-action" disabled={!analysisResult} onClick={()=>analysisResult && downloadHtml(`aga-${analysisResult.analysis_id}.html`, "AlphaGenome Analytics research report", analysisResult)}>Download variant HTML report</button><button className="primary-action" disabled={!analysisResult} onClick={()=>analysisResult && download(`aga-${analysisResult.analysis_id}.json`, analysisResult)}>Download variant JSON</button><button className="primary-action" disabled={!batchAnalysisResult} onClick={()=>batchAnalysisResult && download(`aga-${batchAnalysisResult.analysis_id}.json`, batchAnalysisResult)}>Download batch JSON</button><button className="primary-action" disabled={!stats} onClick={()=>stats && download(`aga-statistics-${String(stats.analysis_id ?? "result")}.json`, stats)}>Download statistics JSON</button><button className="primary-action" disabled={!enrichment} onClick={()=>enrichment && download("aga-gene-enrichment.json", enrichment)}>Download enrichment JSON</button></div>{analysisResult && <div className="scope-note"><b>Analysis provenance</b><span>{analysisResult.alphagenome.provider} · GRCh38 · {analysisResult.created_at} · {analysisResult.alphagenome.sdk_version ?? "SDK version unavailable"} · Model version not exposed by Atlas</span>{analysisResult.evidence.map((entry,index)=><span key={`${entry.provider}-${index}`}>{entry.provider} · retrieved {new Date(entry.retrieved_at).toLocaleString()} · {entry.source_url ?? "source URL not returned"}</span>)}</div>}</>}
      {view === "reproducibility" && <><p className="workspace-copy">History is saved only in this browser. The public server does not keep a shared analysis history; clear this browser’s site data to remove these records.</p>{history.length ? <div className="history-list">{history.map((entry)=><div className="history-entry" key={entry.analysis_id}><button onClick={()=>restoreHistoryEntry(entry)}><b>{entry.analysis_type.replaceAll("_"," ")}</b><span>{entry.created_at} · {entry.status}</span><small>{entry.analysis_id}</small></button><details><summary>View input manifest</summary><pre className="answer-panel">{JSON.stringify(entry.manifest,null,2)}</pre></details></div>)}</div> : <div className="empty-state">No analysis records saved in this browser yet.</div>}</>}
      {view === "region" && <><p className="workspace-copy">Query Ensembl gene and regulatory features for a GRCh38 interval (maximum width 5 Mb). Coordinates are sent to Ensembl when you run the lookup.</p>{result?.variant && <><div className="locus-toolbar"><span className="locus-coordinate">{result.variant.chromosome}:{Math.max(1,result.variant.position-Math.floor(locusWidth/2)).toLocaleString()}–{(result.variant.position+Math.floor(locusWidth/2)).toLocaleString()}</span><div><button onClick={()=>setLocusWidth((n)=>Math.min(1000000,n*2))}>−</button><button onClick={()=>setLocusWidth((n)=>Math.max(200,Math.floor(n/2)))}>+</button></div></div><div className="locus-view"><div className="coordinate-ticks"><span>{Math.max(1,result.variant.position-Math.floor(locusWidth/2)).toLocaleString()}</span><span>{result.variant.position.toLocaleString()}</span><span>{(result.variant.position+Math.floor(locusWidth/2)).toLocaleString()}</span></div><div className="locus-track"><div className="track-line"/><div className="variant-marker" style={{left:"50%"}}><i/><span>{result.variant.reference}→{result.variant.alternate}</span></div></div><div className="track-label">SELECTED VARIANT<span>{result.variant.position.toLocaleString()}</span></div></div></>}<div className="form-grid"><label>Chromosome<input id="region-chrom" placeholder={result?.variant?.chromosome ?? "chr7"}/></label><label>Start<input id="region-start" type="number" placeholder="1000000"/></label><label>End<input id="region-end" type="number" placeholder="1010000"/></label></div><button className="primary-action" onClick={()=>{const c=(document.getElementById("region-chrom") as HTMLInputElement).value;const s=(document.getElementById("region-start") as HTMLInputElement).value;const e=(document.getElementById("region-end") as HTMLInputElement).value;lookup(`/api/v1/evidence/region?chromosome=${encodeURIComponent(c)}&start=${s}&end=${e}`);}} disabled={busy}>Load region annotations</button>{evidenceData && <pre className="answer-panel">{JSON.stringify(evidenceData,null,2)}</pre>}</>}
      {view === "region" && <button className="primary-action" onClick={()=>{const c=(document.getElementById("region-chrom") as HTMLInputElement).value;const s=(document.getElementById("region-start") as HTMLInputElement).value;const e=(document.getElementById("region-end") as HTMLInputElement).value;lookup(`/api/v1/evidence/encode/region?chromosome=${encodeURIComponent(c)}&start=${s}&end=${e}`);}} disabled={busy}>ENCODE region search · configured service</button>}
      {message && view !== "explore" && <p className="feedback visible" role="status">{message}</p>}
    </section>
    <footer><span>ALPHAGENOME ANALYTICS</span><span>Research predictions · Not for clinical decision-making</span></footer>
  </main>;
}

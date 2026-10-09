export type Variant = {
  chromosome: string;
  position: number;
  reference: string;
  alternate: string;
  genome_assembly: "GRCh38";
};

export type Effect = {
  modality: string;
  feature?: string | null;
  gene?: string | null;
  tissue?: string | null;
  raw_score: number;
  quantile_score?: number | null;
};

export type FeatureAttribution = { feature: string; value: number };

export type VariantPrediction = {
  analysis_id: string;
  created_at: string;
  status?: "complete";
  provider?: string;
  provider_sdk_version?: string | null;
  provider_model_version?: string | null;
  requested_scorers?: string[];
  variant?: Variant;
  reference_validation?: {
    status: "verified";
    assembly: "GRCh38";
    source: string;
    source_url: string;
    requested_reference: string;
    observed_reference: string;
    checked_at: string;
  };
  avi_score?: number;
  avi_quantile?: number | null;
  feature_importance?: FeatureAttribution[];
  effects?: Effect[];
  effects_total?: number;
  effects_total_by_scorer?: Record<string, number>;
  effects_truncated?: boolean;
  provider_response_available?: boolean;
  provider_response_endpoint?: string | null;
};

export type EvidenceAttachment = {
  provider: string;
  retrieved_at: string;
  source_url?: string;
  query?: string;
  data: unknown;
};

export type AnalysisStatistic = {
  kind: "score_group_comparison" | "gene_set_enrichment" | "tissue_comparison";
  created_at: string;
  input: unknown;
  data: unknown;
};

export type AnalysisResult = {
  schema_version: "1.0";
  analysis_id: string;
  created_at: string;
  input: {
    kind: "single_variant" | "vcf_batch";
    genome_assembly: "GRCh38";
    variants: Variant[];
    requested_variant_count?: number;
  };
  alphagenome: {
    provider: string;
    sdk_version?: string | null;
    model_version?: string | null;
    predictions: VariantPrediction[];
  };
  annotations: {
    genes: string[];
    tissues: string[];
    effects: Effect[];
    feature_importance: FeatureAttribution[];
  };
  evidence: EvidenceAttachment[];
  statistics: AnalysisStatistic[];
  batch?: {
    status: string;
    completed: number;
    failures: unknown[];
  };
  provenance: {
    generated_at: string;
    sources: Array<{ name: string; version?: string | null; url?: string }>;
  };
};

type BuildAnalysisInput = {
  analysis_id: string;
  created_at: string;
  kind: AnalysisResult["input"]["kind"];
  predictions: VariantPrediction[];
  requested_variant_count?: number;
  batch?: AnalysisResult["batch"];
  evidence?: EvidenceAttachment[];
  statistics?: AnalysisStatistic[];
};

export function buildAnalysisResult(input: BuildAnalysisInput): AnalysisResult {
  const effects = input.predictions.flatMap((prediction) => prediction.effects ?? []);
  const featureImportance = input.predictions.flatMap((prediction) => prediction.feature_importance ?? []);
  const providers = [...new Set(input.predictions.map((prediction) => prediction.provider).filter((value): value is string => Boolean(value)))];
  const evidence = input.evidence ?? [];

  return {
    schema_version: "1.0",
    analysis_id: input.analysis_id,
    created_at: input.created_at,
    input: {
      kind: input.kind,
      genome_assembly: "GRCh38",
      variants: input.predictions.flatMap((prediction) => prediction.variant ? [prediction.variant] : []),
      ...(input.requested_variant_count === undefined ? {} : { requested_variant_count: input.requested_variant_count }),
    },
    alphagenome: {
      provider: providers.join(", ") || "AlphaGenome Atlas",
      sdk_version: input.predictions.find((prediction) => prediction.provider_sdk_version)?.provider_sdk_version,
      model_version: input.predictions.find((prediction) => prediction.provider_model_version)?.provider_model_version,
      predictions: input.predictions,
    },
    annotations: {
      genes: [...new Set(effects.map((effect) => effect.gene).filter((gene): gene is string => Boolean(gene)))],
      tissues: [...new Set(effects.map((effect) => effect.tissue).filter((tissue): tissue is string => Boolean(tissue)))],
      effects,
      feature_importance: featureImportance,
    },
    evidence,
    statistics: input.statistics ?? [],
    ...(input.batch ? { batch: input.batch } : {}),
    provenance: {
      generated_at: new Date().toISOString(),
      sources: [
      ...providers.map((name) => ({
          name,
          version: input.predictions.find((prediction) => prediction.provider === name)?.provider_sdk_version,
          url: "https://deepmind.google.com/science/alphagenome/",
      })),
        ...[...new Map(input.predictions.flatMap((prediction) => {
          const check = prediction.reference_validation;
          return check ? [[check.source, { name: check.source, url: check.source_url, version: check.assembly }] as const] : [];
        })).values()],
        ...evidence.map((item) => ({ name: item.provider, url: item.source_url })),
      ],
    },
  };
}

export function isAnalysisResult(value: unknown): value is AnalysisResult {
  if (!value || typeof value !== "object") return false;
  const record = value as Partial<AnalysisResult>;
  return record.schema_version === "1.0"
    && typeof record.analysis_id === "string"
    && Boolean(record.alphagenome && Array.isArray(record.alphagenome.predictions));
}

# AlphaGenome Analytics

An early research tool scaffold with a Next.js interface and a FastAPI service boundary for variant analysis.

See [`docs/ROADMAP.md`](docs/ROADMAP.md) for the staged product plan.

## Project layout

- `frontend/` — Next.js + TypeScript user interface
- `backend/` — FastAPI API and analysis service boundary

## Run locally

### Backend

```sh
cd backend
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
uvicorn app.main:app --reload
```

The API is available at `http://localhost:8000`; interactive docs are at `/docs`.

### Frontend

```sh
cd frontend
npm install
npm run dev
```

The frontend expects the API at `http://localhost:8000` by default. Set `NEXT_PUBLIC_API_URL` to the deployed API URL for production; it is included in the Render Blueprint.

## AlphaGenome API key

Save `ALPHAGENOME_API_KEY=your_key` in `backend/.env`. The backend loads this file locally; `.env` is excluded from Git. The single-variant endpoint queries AlphaGenome Atlas for the AVI score and feature attributions. Use is subject to AlphaGenome's terms; predictions are for research and theoretical modelling, not clinical decision-making.

## Available workflows

The app includes single-variant prediction, effect review, AVI-based ranking, local VCF parsing and batch jobs, Ensembl/GWAS evidence lookups, tissue summaries, candidate gene labels, score statistics, gene-set enrichment, a feature-to-gene view, deterministic Ask AGA summaries, HTML/JSON exports, and local analysis history. See [`docs/ROADMAP.md`](docs/ROADMAP.md) for supported data sources and limitations.

The AlphaGenome key is used only by the backend. Variant and region lookups send queries to the selected providers. VCF contents are parsed locally in the browser; variant coordinates are sent to AlphaGenome Atlas only after the user confirms. In public deployment mode, server-side analysis persistence is disabled and each visitor's history is kept in that visitor's browser. Local development keeps SQLite persistence enabled by default.

The tool does not calculate clinical interpretations. It connects to Ensembl, GWAS Catalog, ClinVar, GTEx eQTL, gnomAD, and SCREEN's GraphQL API for ENCODE Registry cCRE annotations. SCREEN region search requires a `SCREEN_API_KEY` configured on the backend; save it in `backend/.env` locally and as a secret backend environment variable in Render. Keys expire after 90 days. cCRE spans and source signal annotations are kept separate from AlphaGenome predictions. No live genomic provider requests were made while implementing these workflows.

## Deploy to Render

The root [`render.yaml`](render.yaml) defines a public Next.js web service and a FastAPI web service. Set `ALPHAGENOME_API_KEY` and `SCREEN_API_KEY` as backend secrets in Render. Do not put either key in `render.yaml`, a frontend environment variable, or Git. Render assigns the service URLs from their names; if you change either service name or add a custom domain, update `NEXT_PUBLIC_API_URL` and the backend's `CORS_ORIGINS` to match.

Public mode disables server-side analysis history, protects batch polling with a per-job bearer token, and limits each batch to 100 SNVs. Browser history is local to that browser. Batch jobs run in memory and are lost when the API service restarts or scales to another instance. Render's free web services may sleep when idle, so the first request after inactivity can take longer. For a sustained community service, use a paid always-on API instance and a durable queue before raising batch limits.

The public API spends the configured AlphaGenome account's query capacity. Before sharing the URL widely, set appropriate usage limits for the AlphaGenome account and review AlphaGenome's terms. Analysis outputs remain research predictions and are not for clinical decision-making.

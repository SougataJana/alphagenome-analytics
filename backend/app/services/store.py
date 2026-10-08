"""Local persistence for reproducible analysis records and result payloads."""

import json
import os
import sqlite3
import threading
from pathlib import Path
from typing import Any

DB_PATH = Path(__file__).resolve().parents[2] / "data" / "aga.sqlite3"
_lock = threading.Lock()


def _persistence_enabled() -> bool:
    return os.getenv("ANALYSIS_PERSISTENCE", "true").lower() in {"1", "true", "yes"}


def _connect() -> sqlite3.Connection:
    DB_PATH.parent.mkdir(parents=True, exist_ok=True)
    connection = sqlite3.connect(DB_PATH, timeout=20)
    connection.row_factory = sqlite3.Row
    connection.execute(
        """CREATE TABLE IF NOT EXISTS analyses (
            analysis_id TEXT PRIMARY KEY,
            created_at TEXT NOT NULL,
            analysis_type TEXT NOT NULL,
            status TEXT NOT NULL,
            manifest_json TEXT NOT NULL,
            results_json TEXT NOT NULL
        )"""
    )
    return connection


def save_analysis(
    analysis_id: str,
    created_at: str,
    analysis_type: str,
    status: str,
    manifest: dict[str, Any],
    results: dict[str, Any],
) -> None:
    if not _persistence_enabled():
        return
    with _lock, _connect() as connection:
        connection.execute(
            """INSERT OR REPLACE INTO analyses
            (analysis_id, created_at, analysis_type, status, manifest_json, results_json)
            VALUES (?, ?, ?, ?, ?, ?)""",
            (
                analysis_id,
                created_at,
                analysis_type,
                status,
                json.dumps(manifest, allow_nan=False),
                json.dumps(results, allow_nan=False),
            ),
        )


def get_analysis(analysis_id: str) -> dict[str, Any] | None:
    if not _persistence_enabled():
        return None
    with _lock, _connect() as connection:
        row = connection.execute(
            "SELECT * FROM analyses WHERE analysis_id = ?", (analysis_id,)
        ).fetchone()
    if row is None:
        return None
    return {
        "analysis_id": row["analysis_id"],
        "created_at": row["created_at"],
        "analysis_type": row["analysis_type"],
        "status": row["status"],
        "manifest": json.loads(row["manifest_json"]),
        "results": json.loads(row["results_json"]),
    }


def list_analyses(limit: int = 50) -> list[dict[str, Any]]:
    if not _persistence_enabled():
        return []
    with _lock, _connect() as connection:
        rows = connection.execute(
            "SELECT analysis_id, created_at, analysis_type, status, manifest_json "
            "FROM analyses ORDER BY created_at DESC LIMIT ?",
            (limit,),
        ).fetchall()
    return [
        {
            "analysis_id": row["analysis_id"],
            "created_at": row["created_at"],
            "analysis_type": row["analysis_type"],
            "status": row["status"],
            "manifest": json.loads(row["manifest_json"]),
        }
        for row in rows
    ]

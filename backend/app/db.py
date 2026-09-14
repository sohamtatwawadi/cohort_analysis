"""DuckDB connection handling for the analytics store.

The store is derived and rebuildable. Every number the UI shows is computed by
SQL here, never in the browser and never by a model (spec D01).
"""
from __future__ import annotations

import threading
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence

import duckdb

from . import config

_lock = threading.RLock()
_conn: Optional[duckdb.DuckDBPyConnection] = None


def _schema_sql() -> str:
    """Lab Mode schema plus the Research Mode tables (Part II).

    One store, one connection: Part II §1 says the two modes share the cohort
    engine, governance, audit and export. Splitting the databases would mean
    duplicating all of that.
    """
    here = Path(__file__).parent
    return "\n".join([
        (here / "schema.sql").read_text(),
        (here / "research_schema.sql").read_text(),
    ])


def connect(path: Optional[Path] = None, fresh: bool = False) -> duckdb.DuckDBPyConnection:
    """Open (or reopen) the analytics store and ensure the schema exists."""
    global _conn
    with _lock:
        target = Path(path or config.DB_PATH)
        target.parent.mkdir(parents=True, exist_ok=True)
        if fresh and target.exists():
            close()
            target.unlink()
        if _conn is None:
            _conn = duckdb.connect(str(target))
            _conn.execute(_schema_sql())
        return _conn


def get() -> duckdb.DuckDBPyConnection:
    return _conn if _conn is not None else connect()


def close() -> None:
    global _conn
    with _lock:
        if _conn is not None:
            _conn.close()
            _conn = None


# --------------------------------------------------------------------- query --
def rows(sql: str, params: Sequence[Any] = ()) -> List[Dict[str, Any]]:
    """Run a query and return a list of dicts. All analytics go through here."""
    with _lock:
        cur = get().execute(sql, list(params))
        cols = [d[0] for d in cur.description]
        return [dict(zip(cols, r)) for r in cur.fetchall()]


def row(sql: str, params: Sequence[Any] = ()) -> Optional[Dict[str, Any]]:
    r = rows(sql, params)
    return r[0] if r else None


def scalar(sql: str, params: Sequence[Any] = (), default: Any = None) -> Any:
    with _lock:
        r = get().execute(sql, list(params)).fetchone()
    return default if r is None or r[0] is None else r[0]


def execute(sql: str, params: Sequence[Any] = ()) -> None:
    with _lock:
        get().execute(sql, list(params))


def executemany(sql: str, seq: Sequence[Sequence[Any]]) -> None:
    if not seq:
        return
    with _lock:
        get().executemany(sql, [list(s) for s in seq])


def insert_rows(table: str, columns: Sequence[str], data: Sequence[Sequence[Any]]) -> None:
    if not data:
        return
    placeholders = ", ".join(["?"] * len(columns))
    sql = "INSERT INTO {} ({}) VALUES ({})".format(table, ", ".join(columns), placeholders)
    executemany(sql, data)


def table_counts() -> Dict[str, int]:
    tables = [
        "subject", "family", "sample", "run", "run_scope", "run_fingerprint",
        "finding", "interpretation", "gene_disease", "test_code", "cohort_def",
    ]
    out: Dict[str, int] = {}
    for t in tables:
        try:
            out[t] = int(scalar("SELECT COUNT(*) FROM {}".format(t), default=0))
        except duckdb.Error:
            out[t] = 0
    return out


def meta_get(key: str, default: Optional[str] = None) -> Optional[str]:
    v = scalar("SELECT value FROM store_meta WHERE key = ?", [key])
    return default if v is None else str(v)


def meta_set(key: str, value: str) -> None:
    execute("DELETE FROM store_meta WHERE key = ?", [key])
    execute("INSERT INTO store_meta (key, value) VALUES (?, ?)", [key, value])

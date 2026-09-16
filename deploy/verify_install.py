"""Post-install smoke check: can the app be imported, and are its routes there?

Run before enabling the service, so a missing dependency surfaces as one clear
message instead of a unit that crash-loops in the background.

    PYTHONPATH=/opt/cohort /opt/cohort/.venv/bin/python deploy/verify_install.py

Counting len(app.routes) was the obvious check and the wrong one: FastAPI keeps
an included router as a single object rather than flattening its paths into the
parent, so a complete install reported "routes: 10" and looked broken. The
number depends on FastAPI's internals; the presence of specific endpoints does
not. So this walks into the routers and asserts on paths.
"""
from __future__ import annotations

import sys

# Endpoints that must exist. One per router, so a router failing to register is
# caught rather than inferred from a count.
REQUIRED = [
    "/health",                       # the app itself
    "/api/cohort/resolve",           # the lab router
    "/api/research/meta",            # the research router
]


def walk(router, prefix: str = ""):
    """Every concrete path under `router`, following included routers.

    FastAPI has changed how inclusion is represented more than once, so this
    looks for a nested router under either of the attribute names used, rather
    than assuming one shape.
    """
    for route in getattr(router, "routes", []):
        path = getattr(route, "path", None)
        if isinstance(path, str):
            yield prefix + path

        inner = getattr(route, "original_router", None)
        if inner is None:
            candidate = getattr(route, "app", None)
            if hasattr(candidate, "routes"):
                inner = candidate
        if inner is None or inner is router:
            continue

        ctx = getattr(route, "include_context", None)
        yield from walk(inner, prefix + (getattr(ctx, "prefix", "") or ""))


def main() -> int:
    try:
        from backend.app.main import app
    except Exception as exc:                            # noqa: BLE001
        print("  could not import the app: {}: {}".format(type(exc).__name__, exc))
        raise

    found = set(walk(app))
    print("  endpoints: {}".format(len(found)))

    missing = [r for r in REQUIRED if r not in found]
    if missing:
        print("  MISSING: {}".format(", ".join(missing)))
        print("  found instead: {}".format(
            ", ".join(sorted(found)[:12]) or "nothing"))
        return 1

    print("  required routes present")
    return 0


if __name__ == "__main__":
    sys.exit(main())

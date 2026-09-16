# Runtime image for the cohort platform.
#
# Slim rather than alpine: numpy and scipy publish manylinux wheels but not musl
# ones, so alpine would compile them from source — a long build and a fragile
# one. Slim installs the prebuilt wheels in seconds.

FROM python:3.11-slim

# Dependencies first, so a code change does not re-resolve numpy and scipy.
WORKDIR /app
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY backend/ ./backend/
COPY frontend/ ./frontend/

# The analytics store and any uploaded dataset live here. Mount a volume over
# it — without one, a container restart silently discards every uploaded
# dataset, and the app re-seeds synthetic data as though nothing was lost.
RUN mkdir -p /data
ENV COHORT_DB=/data/germline.duckdb \
    EXPORT_DIR=/data/exports \
    PYTHONUNBUFFERED=1 \
    PORT=8000
VOLUME ["/data"]

# Run unprivileged. The process needs to write only /data.
RUN useradd --create-home --uid 10001 cohort && chown -R cohort:cohort /data
USER cohort

EXPOSE 8000

# Container health is "can it answer", not "is the process alive" — the app
# seeds its store on first boot and is genuinely not ready until that finishes.
HEALTHCHECK --interval=30s --timeout=5s --start-period=180s --retries=3 \
  CMD python -c "import urllib.request,os,sys; \
sys.exit(0 if urllib.request.urlopen('http://127.0.0.1:'+os.environ.get('PORT','8000')+'/health',timeout=4).status==200 else 1)"

# One worker, deliberately. DuckDB takes a single writer lock on the database
# file, so a second worker process cannot open it — it would crash on startup
# with "Conflicting lock". Concurrency within the process is handled by the
# cohort lock and the job thread pool.
CMD ["sh", "-c", "exec uvicorn backend.app.main:app --host 0.0.0.0 --port ${PORT:-8000} --workers 1"]

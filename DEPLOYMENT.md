# Deploying on AWS

Measured figures, not estimates — every number here came from running the thing.

## Sizing

| What | Measured |
|---|---|
| Idle, after the store has seeded | **210 MB** |
| Loading the 150k-variant research dataset | **459 MB** |
| Peak during a GWAS on that dataset | **~1.0 GB** |
| DuckDB store + research payloads on disk | ~130 MB |
| GWAS wall time (multi-core dev machine) | ~4 min |

| Instance | Verdict |
|---|---|
| `t3.micro` (1 GB) | Lab Mode only. A research analysis will be OOM-killed. |
| **`t3.small` (2 GB)** | **Minimum for Research Mode.** Tight but workable. |
| **`t3.medium` (2 vCPU / 4 GB)** | **Recommended.** Headroom for the GWAS, and the second vCPU matters — the scan is BLAS-heavy. |
| `t3.large` (8 GB) | Only if you load datasets larger than the demo. |

Memory scales with `variants × samples`: a 150,000 × 3,000 matrix is 429 MB as
int8 and the analyses need working copies on top. Use burstable `t3`/`t4g`
rather than a fixed-small shape — analyses are spiky, idle is cheap.

> `t4g` (ARM Graviton) works and is cheaper. numpy, scipy and duckdb all ship
> arm64 wheels.

## Disk

Give it **20 GB** and keep the data on a volume that survives the instance:
`data/` holds genotype and clinical records, so losing it loses uploaded
datasets. Nothing under `data/` is in git — by design — so there is no restore
path if it is on ephemeral storage.

## Security group

Inbound **443 only** (and 22 from your own IP). Do **not** open 8000.

The application has **no authentication of any kind.** Anyone who reaches it can
resolve cohorts, upload datasets and export data. That is fine behind a login
and not fine on an open port, so nginx terminates TLS and carries basic auth
until real auth exists — see `deploy/nginx.conf`. If this ever holds real
patient data, a free-standing basic-auth box is not an adequate control and the
conversation is about a VPC, a BAA and audited access instead.

---

## Option A — Docker

```bash
sudo dnf install -y docker && sudo systemctl enable --now docker   # Amazon Linux 2023
sudo usermod -aG docker $USER && newgrp docker

git clone https://github.com/sohamtatwawadi/cohort_analysis.git && cd cohort_analysis
docker compose up -d --build
docker compose logs -f          # first boot seeds the store; wait for "startup complete"
```

The compose file publishes to `127.0.0.1:8000` only, mounts a named volume at
`/data`, and caps memory at 4 GB. One worker, deliberately: DuckDB takes a
single writer lock on the database file, so a second worker cannot open it and
exits with `Conflicting lock`.

> The Dockerfile has not been built and run here — no Docker on the dev machine.
> The systemd path below was verified end to end. If the image misbehaves, that
> is the more trusted route.

## Option B — systemd (verified)

```bash
sudo dnf install -y python3.11 python3.11-pip nginx git
sudo useradd --system --home /opt/cohort cohort
sudo mkdir -p /opt/cohort /var/lib/cohort && sudo chown -R cohort:cohort /opt/cohort /var/lib/cohort

sudo -u cohort git clone https://github.com/sohamtatwawadi/cohort_analysis.git /opt/cohort
cd /opt/cohort
sudo -u cohort python3.11 -m venv .venv
sudo -u cohort .venv/bin/pip install -r requirements.txt

sudo cp deploy/cohort.service /etc/systemd/system/
sudo systemctl daemon-reload && sudo systemctl enable --now cohort
sudo journalctl -u cohort -f
```

Then TLS and the proxy:

```bash
sudo cp deploy/nginx.conf /etc/nginx/conf.d/cohort.conf   # edit server_name first
sudo htpasswd -c /etc/nginx/.htpasswd demo
sudo dnf install -y certbot python3-certbot-nginx && sudo certbot --nginx
```

## First boot

The store seeds itself on an empty database — about a minute here, longer on a
small instance. `/health` returns row counts once it is ready, which is what the
container healthcheck waits on:

```bash
curl -s localhost:8000/health
```

## The research demo dataset

Not in the repo (it is data), so generate it on the instance. **Stop the app
first** — DuckDB allows one writer:

```bash
sudo systemctl stop cohort
sudo -u cohort COHORT_DB=/var/lib/cohort/germline.duckdb \
  /opt/cohort/.venv/bin/python -m backend.tools.make_research_fixture \
  --samples 3000 --variants 150000 --seed 20260915 \
  --name "Cardiomyopathy case-control · demo cohort"
sudo systemctl start cohort
```

Takes several minutes and needs ~1.5 GB free. On a 2 GB instance, drop to
`--variants 60000` — still `genome_wide`, so nothing locks.

## Environment variables

All have working defaults; none are required.

| Variable | Default | Set it when |
|---|---|---|
| `PORT` | `8000` | The platform assigns one |
| `COHORT_DB` | `data/germline.duckdb` | **Always** — point at persistent storage |
| `EXPORT_DIR` | `data/exports` | Alongside the database |
| `VARIMAT_DIR` / `CLINICAL_DIR` | `data/varimat`, `data/clinical` | Lab ingest reads elsewhere |
| `CORS_ORIGINS` | *unset — CORS off* | A different origin serves the UI |
| `COHORT_USER` | `analyst@impactomics` | Attribution in the audit log |
| `TENANT_ID` | `11` | Multi-tenant separation |
| `KB_SNAPSHOT_ID` | `kb-2026-07-28` | Pinning a knowledge-base release |

## Upgrading

```bash
cd /opt/cohort && sudo -u cohort git pull
sudo -u cohort .venv/bin/pip install -r requirements.txt
sudo systemctl restart cohort
```

The analytics store is derived and rebuildable; uploaded research datasets are
not. A restart is safe — jobs interrupted by it are reconciled to `failed` on
startup rather than being left to poll forever.

## Backups

```bash
sudo systemctl stop cohort
sudo tar czf /tmp/cohort-$(date +%F).tgz -C /var/lib cohort
sudo systemctl start cohort
```

Copy the DuckDB file only while stopped — a live copy can catch a partial WAL.

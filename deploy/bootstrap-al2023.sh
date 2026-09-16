#!/usr/bin/env bash
# One-shot server setup for Amazon Linux 2023.
#
# Idempotent: safe to re-run after a failure or to pick up a code change. Every
# step checks before acting, so a second run repairs rather than duplicates.
#
#   curl -fsSL https://raw.githubusercontent.com/sohamtatwawadi/cohort_analysis/main/deploy/bootstrap-al2023.sh | sudo bash
#
# or, from a clone:  sudo bash deploy/bootstrap-al2023.sh

set -euo pipefail

REPO="${REPO:-https://github.com/sohamtatwawadi/cohort_analysis.git}"
APP_DIR="${APP_DIR:-/opt/cohort}"
DATA_DIR="${DATA_DIR:-/var/lib/cohort}"
APP_USER="${APP_USER:-cohort}"
PY="${PY:-python3.11}"

say() { printf '\n\033[1;36m==> %s\033[0m\n' "$*"; }
die() { printf '\n\033[1;31mFAILED: %s\033[0m\n' "$*" >&2; exit 1; }

[[ $EUID -eq 0 ]] || die "run with sudo"

# ---------------------------------------------------------------- packages --
say "Installing packages"
dnf -q -y install git nginx httpd-tools "$PY" "${PY}-pip" tar gzip \
  || die "package install failed — check 'dnf search python3.11'"
"$PY" --version || die "$PY not on PATH after install"

# ------------------------------------------------------------------- user --
say "Creating service user and directories"
id -u "$APP_USER" >/dev/null 2>&1 || useradd --system --home "$APP_DIR" --shell /sbin/nologin "$APP_USER"
mkdir -p "$APP_DIR" "$DATA_DIR"
chown -R "$APP_USER:$APP_USER" "$APP_DIR" "$DATA_DIR"
# The database and any uploaded genotype data live here; nothing else needs to
# read them.
chmod 750 "$DATA_DIR"

# ------------------------------------------------------------------- code --
if [[ -d "$APP_DIR/.git" ]]; then
  say "Updating existing checkout"
  sudo -u "$APP_USER" git -C "$APP_DIR" pull --ff-only
else
  say "Cloning $REPO"
  # Clone into a temp dir and move, because useradd --home already created
  # $APP_DIR and git refuses to clone into a non-empty directory.
  rm -rf /tmp/cohort-clone
  sudo -u "$APP_USER" git clone --depth 1 "$REPO" /tmp/cohort-clone
  shopt -s dotglob
  mv /tmp/cohort-clone/* "$APP_DIR"/
  shopt -u dotglob
  rmdir /tmp/cohort-clone
  chown -R "$APP_USER:$APP_USER" "$APP_DIR"
fi

# ---------------------------------------------------------------- python --
say "Installing Python dependencies"
[[ -x "$APP_DIR/.venv/bin/python" ]] || sudo -u "$APP_USER" "$PY" -m venv "$APP_DIR/.venv"
sudo -u "$APP_USER" "$APP_DIR/.venv/bin/pip" install -q --upgrade pip
sudo -u "$APP_USER" "$APP_DIR/.venv/bin/pip" install -q -r "$APP_DIR/requirements.txt" \
  || die "pip install failed"

# Prove the dependency set is complete before wiring up a service that would
# otherwise crash-loop on an ImportError. The check lives in its own file so it
# can be tested rather than buried in a heredoc.
say "Verifying the app imports"
sudo -u "$APP_USER" env PYTHONPATH="$APP_DIR" \
  "$APP_DIR/.venv/bin/python" "$APP_DIR/deploy/verify_install.py" \
  || die "the app did not import cleanly — see above"

# --------------------------------------------------------------- service --
say "Installing the systemd unit"
install -m 644 "$APP_DIR/deploy/cohort.service" /etc/systemd/system/cohort.service
systemctl daemon-reload
systemctl enable cohort >/dev/null
systemctl restart cohort

# ----------------------------------------------------------------- nginx --
say "Configuring nginx"
install -m 644 "$APP_DIR/deploy/nginx-http.conf" /etc/nginx/conf.d/cohort.conf
# AL2023 ships a default server on :80 that would otherwise win the
# default_server slot and serve the nginx welcome page instead of the app.
if grep -q "default_server" /etc/nginx/nginx.conf 2>/dev/null; then
  sed -i 's/listen       80 default_server;/listen       80;/' /etc/nginx/nginx.conf || true
  sed -i 's/listen       \[::\]:80 default_server;/listen       [::]:80;/' /etc/nginx/nginx.conf || true
fi

if [[ ! -f /etc/nginx/.htpasswd ]]; then
  # openssl, not `tr </dev/urandom | head -c 16`. head closes the pipe after 16
  # bytes, tr dies with SIGPIPE, and `set -o pipefail` turns that into exit 141
  # — which `set -e` then acts on. The script died here, silently, because this
  # line had no `|| die` to report it.
  PASS="$(openssl rand -hex 8)"
  htpasswd -bc /etc/nginx/.htpasswd demo "$PASS" >/dev/null 2>&1
  chown root:nginx /etc/nginx/.htpasswd && chmod 640 /etc/nginx/.htpasswd
  printf '%s' "$PASS" > /root/cohort-demo-password
  chmod 600 /root/cohort-demo-password
  CREATED_PASS=1
fi

# SELinux is permissive on a stock AL2023, but if it has been set to enforcing
# this boolean is what stops nginx connecting to the app on loopback.
if command -v getenforce >/dev/null && [[ "$(getenforce)" == "Enforcing" ]]; then
  setsebool -P httpd_can_network_connect 1 || true
fi

nginx -t || die "nginx config is invalid"
systemctl enable nginx >/dev/null
systemctl restart nginx

# ----------------------------------------------------------------- verify --
say "Waiting for the app to answer (first boot seeds the store)"
for i in $(seq 1 120); do
  if curl -fsS --max-time 3 http://127.0.0.1:8000/health >/dev/null 2>&1; then
    printf '  ready after %ds\n' $((i * 3)); break
  fi
  [[ $i -eq 120 ]] && die "no response after 6 min — check: journalctl -u cohort -n 50"
  sleep 3
done

IP="$(curl -fsS --max-time 3 http://169.254.169.254/latest/meta-data/public-ipv4 2>/dev/null || echo '<instance-ip>')"

printf '\n\033[1;32mDone.\033[0m\n\n'
printf '  URL       http://%s/\n' "$IP"
printf '  Username  demo\n'
if [[ "${CREATED_PASS:-}" == "1" ]]; then
  printf '  Password  %s\n' "$(cat /root/cohort-demo-password)"
  printf '            (also saved to /root/cohort-demo-password)\n'
else
  printf '  Password  unchanged — see /root/cohort-demo-password\n'
fi
printf '\n  Security group must allow inbound 80 from your IP. Do not open 8000.\n'
printf '  Logs:  journalctl -u cohort -f\n\n'
printf '  Research Mode has no dataset yet. To add the demo cohort:\n'
printf '    sudo systemctl stop cohort\n'
printf '    sudo -u %s env COHORT_DB=%s/germline.duckdb PYTHONPATH=%s \\\n' \
  "$APP_USER" "$DATA_DIR" "$APP_DIR"
printf '      %s/.venv/bin/python -m backend.tools.make_research_fixture \\\n' "$APP_DIR"
printf '      --samples 3000 --variants 150000 --seed 20260915 \\\n'
printf '      --name "Cardiomyopathy case-control · demo cohort"\n'
printf '    sudo systemctl start cohort\n\n'

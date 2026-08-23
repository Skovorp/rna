#!/usr/bin/env bash
set -euo pipefail

ROOT="/opt/rna-atlas"
VENV="/opt/rna-venv"
cd "$ROOT"

git fetch origin main
before="$(git rev-parse HEAD)"
after="$(git rev-parse origin/main)"
if [[ "$before" == "$after" ]]; then
  exit 0
fi

git reset --hard origin/main

if ! git diff --quiet "$before" "$after" -- app/requirements.txt; then
  "$VENV/bin/pip" install -r app/requirements.txt
fi

# Rebuild the derived-dataset pickle before restarting so the first visitor
# never pays the parse.
(cd app && "$VENV/bin/python" -c \
  "from expression_explorer.data import load_datasets; load_datasets('$ROOT/expression')") || true

cp deploy/vps/rna-atlas.service deploy/vps/rna-atlas-update.service deploy/vps/rna-atlas-update.timer /etc/systemd/system/
systemctl daemon-reload
systemctl restart rna-atlas.service

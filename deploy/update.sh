#!/usr/bin/env bash
set -euo pipefail

ROOT="/home/pi-rus/Downloads/feral-remote/aedes-rna-atlas"
cd "$ROOT"

git fetch origin main
before="$(git rev-parse HEAD)"
after="$(git rev-parse origin/main)"
if [[ "$before" == "$after" ]]; then
  exit 0
fi

git pull --ff-only origin main

source /home/pi-rus/miniforge3/etc/profile.d/conda.sh
conda activate aedes-rna-atlas
python -m pip install -r app/requirements.txt
PYTHONPATH="$ROOT/app" python -m pytest -q app/tests
systemctl --user restart aedes-rna-atlas.service


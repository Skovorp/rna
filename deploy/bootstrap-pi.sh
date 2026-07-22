#!/usr/bin/env bash
set -euo pipefail

ROOT="/home/pi-rus/Downloads/feral-remote/aedes-rna-atlas"
source /home/pi-rus/miniforge3/etc/profile.d/conda.sh

if [[ ! -x /home/pi-rus/miniforge3/envs/aedes-rna-atlas/bin/python ]]; then
  conda create -n aedes-rna-atlas python=3.11 -y
fi

conda activate aedes-rna-atlas
python -m pip install -r "$ROOT/app/requirements.txt"
PYTHONPATH="$ROOT/app" python -m pytest -q "$ROOT/app/tests"

cp "$ROOT/deploy/aedes-rna-atlas.service" ~/.config/systemd/user/
cp "$ROOT/deploy/aedes-rna-atlas-update.service" ~/.config/systemd/user/
cp "$ROOT/deploy/aedes-rna-atlas-update.timer" ~/.config/systemd/user/
systemctl --user daemon-reload
systemctl --user enable --now aedes-rna-atlas.service aedes-rna-atlas-update.timer


#!/usr/bin/env bash
set -euo pipefail

APP_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RNA_DIR="$(cd "$APP_DIR/.." && pwd)"
VENV="$RNA_DIR/.venv"

if command -v python3.10 >/dev/null 2>&1; then
  PYTHON="$(command -v python3.10)"
else
  PYTHON="$(command -v python3)"
fi

if [[ ! -x "$VENV/bin/python" ]]; then
  "$PYTHON" -m venv "$VENV"
fi

"$VENV/bin/python" -m pip install --upgrade pip
"$VENV/bin/python" -m pip install -r "$APP_DIR/requirements.txt"

echo "Environment ready: $VENV"


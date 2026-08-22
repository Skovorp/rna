#!/usr/bin/env bash
set -euo pipefail

APP_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RNA_DIR="$(cd "$APP_DIR/.." && pwd)"

if [[ ! -x "$RNA_DIR/.venv/bin/streamlit" ]]; then
  echo "Missing local environment. Run: $APP_DIR/setup.sh" >&2
  exit 1
fi

cd "$APP_DIR"

# Pre-build the derived-dataset pickle so the first visitor after a data update
# skips the multi-second parse. Never let a warm-up failure block the server.
FERAL_EXPRESSION_DIR="$RNA_DIR/expression" "$RNA_DIR/.venv/bin/python" -c \
  "import os; from expression_explorer.data import load_datasets; load_datasets(os.environ['FERAL_EXPRESSION_DIR'])" \
  || true

exec "$RNA_DIR/.venv/bin/streamlit" run app.py "$@"


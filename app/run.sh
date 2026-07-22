#!/usr/bin/env bash
set -euo pipefail

APP_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RNA_DIR="$(cd "$APP_DIR/.." && pwd)"

if [[ ! -x "$RNA_DIR/.venv/bin/streamlit" ]]; then
  echo "Missing local environment. Run: $APP_DIR/setup.sh" >&2
  exit 1
fi

cd "$APP_DIR"
exec "$RNA_DIR/.venv/bin/streamlit" run app.py "$@"


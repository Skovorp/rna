#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=config.env
source "$SCRIPT_DIR/config.env"

SOFTWARE_DIR="$RNA_ROOT/software"
DOWNLOAD_DIR="$SOFTWARE_DIR/downloads"
STATE_DIR="$RNA_ROOT/state"
LOG_DIR="$RNA_ROOT/logs"

mkdir -p \
  "$DOWNLOAD_DIR" \
  "$STATE_DIR" \
  "$LOG_DIR" \
  "$RNA_ROOT/reference" \
  "$RNA_ROOT/raw" \
  "$RNA_ROOT/quant" \
  "$RNA_ROOT/code_snapshots"

apt-get update
DEBIAN_FRONTEND=noninteractive apt-get install -y \
  ca-certificates \
  curl \
  git \
  jq \
  parallel \
  pigz \
  procps \
  rsync \
  time \
  tmux \
  xz-utils

SALMON_ARCHIVE="salmon-cli-x86_64-unknown-linux-gnu.tar.xz"
SALMON_DOWNLOAD_DIR="$DOWNLOAD_DIR/salmon-$SALMON_VERSION"
SALMON_INSTALL_DIR="$SOFTWARE_DIR/salmon-$SALMON_VERSION"
mkdir -p "$SALMON_DOWNLOAD_DIR" "$SALMON_INSTALL_DIR"

if [[ ! -f "$SALMON_DOWNLOAD_DIR/$SALMON_ARCHIVE" ]]; then
  curl \
    --fail \
    --location \
    --retry 8 \
    --retry-all-errors \
    --output "$SALMON_DOWNLOAD_DIR/$SALMON_ARCHIVE" \
    "https://github.com/COMBINE-lab/salmon/releases/download/v$SALMON_VERSION/$SALMON_ARCHIVE"
fi

if [[ ! -f "$SALMON_DOWNLOAD_DIR/$SALMON_ARCHIVE.sha256" ]]; then
  curl \
    --fail \
    --location \
    --retry 8 \
    --retry-all-errors \
    --output "$SALMON_DOWNLOAD_DIR/$SALMON_ARCHIVE.sha256" \
    "https://github.com/COMBINE-lab/salmon/releases/download/v$SALMON_VERSION/$SALMON_ARCHIVE.sha256"
fi

(
  cd "$SALMON_DOWNLOAD_DIR"
  sha256sum --check "$SALMON_ARCHIVE.sha256"
)

if ! find "$SALMON_INSTALL_DIR" -type f -name salmon -perm -u+x -print -quit | grep -q .; then
  tar -xJf "$SALMON_DOWNLOAD_DIR/$SALMON_ARCHIVE" -C "$SALMON_INSTALL_DIR"
fi
SALMON_BIN="$(find "$SALMON_INSTALL_DIR" -type f -name salmon -perm -u+x -print -quit)"
if [[ ! -e /usr/local/bin/salmon ]]; then
  ln -s "$SALMON_BIN" /usr/local/bin/salmon
fi

GFFREAD_ARCHIVE="gffread-$GFFREAD_VERSION.Linux_x86_64.tar.gz"
GFFREAD_DOWNLOAD_DIR="$DOWNLOAD_DIR/gffread-$GFFREAD_VERSION"
GFFREAD_INSTALL_DIR="$SOFTWARE_DIR/gffread-$GFFREAD_VERSION"
mkdir -p "$GFFREAD_DOWNLOAD_DIR" "$GFFREAD_INSTALL_DIR"

if [[ ! -f "$GFFREAD_DOWNLOAD_DIR/$GFFREAD_ARCHIVE" ]]; then
  curl \
    --fail \
    --location \
    --retry 8 \
    --retry-all-errors \
    --output "$GFFREAD_DOWNLOAD_DIR/$GFFREAD_ARCHIVE" \
    "https://github.com/gpertea/gffread/releases/download/v$GFFREAD_VERSION/$GFFREAD_ARCHIVE"
fi

if ! find "$GFFREAD_INSTALL_DIR" -type f -name gffread -perm -u+x -print -quit | grep -q .; then
  tar -xzf "$GFFREAD_DOWNLOAD_DIR/$GFFREAD_ARCHIVE" -C "$GFFREAD_INSTALL_DIR"
fi
GFFREAD_BIN="$(find "$GFFREAD_INSTALL_DIR" -type f -name gffread -perm -u+x -print -quit)"
if [[ ! -e /usr/local/bin/gffread ]]; then
  ln -s "$GFFREAD_BIN" /usr/local/bin/gffread
fi

{
  date --iso-8601=seconds
  salmon --version
  gffread --version
  printf 'logical_cpus=%s\n' "$(nproc)"
  free -h
  df -h "$RNA_ROOT"
} | tee "$STATE_DIR/bootstrap.complete.txt"

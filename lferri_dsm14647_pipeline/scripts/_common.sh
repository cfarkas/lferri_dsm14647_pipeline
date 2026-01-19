#!/usr/bin/env bash
set -euo pipefail

# Repo root (scripts assumed to live in scripts/)
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
THREADS="${THREADS:-80}"

log() { echo "[${0##*/}] $*"; }

die() { echo "[${0##*/}] ERROR: $*" >&2; exit 1; }

require() {
  command -v "$1" >/dev/null 2>&1 || die "Missing dependency '$1' in PATH";
}

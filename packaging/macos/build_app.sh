#!/usr/bin/env bash
# Build MissensExplorer.app and its .dmg on macOS (see packaging/build.py, which builds every system).
#   packaging/macos/build_app.sh [--smoke]
set -euo pipefail
cd "$(dirname "$0")/../.."
PY="${PYTHON:-.venv/bin/python}"
"$PY" -m pip install -q ".[gui]" -r packaging/requirements-installer.txt
"$PY" packaging/build.py --installer "$@"

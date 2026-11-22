#!/usr/bin/env bash
# Build a standalone MissensExplorer.app (macOS) with PyInstaller.
#
#   packaging/macos/build_app.sh             -> build/macos/dist/MissensExplorer.app
#   packaging/macos/build_app.sh --install   -> also copies it to ~/Applications
set -euo pipefail
cd "$(dirname "$0")/../.."
PY="${PYTHON:-.venv/bin/python}"
"$PY" -m pip install -q -e ".[gui]" "pyinstaller>=6.10"
VERSION=$("$PY" -c "import missensexplorer; print(missensexplorer.__version__)")

"$PY" -m PyInstaller --noconfirm --clean --windowed \
  --name MissensExplorer \
  --icon "$PWD/packaging/macos/MissensExplorer.icns" \
  --osx-bundle-identifier org.missensexplorer.desktop \
  --distpath build/macos/dist --workpath build/macos/work --specpath build/macos \
  --hidden-import PySide6.QtWebEngineWidgets --hidden-import PySide6.QtWebEngineCore \
  --exclude-module streamlit --exclude-module altair --exclude-module pyarrow --exclude-module tornado \
  --exclude-module matplotlib --exclude-module IPython --exclude-module pytest --exclude-module tkinter \
  packaging/macos/launcher.py

APP=build/macos/dist/MissensExplorer.app
plutil -replace CFBundleShortVersionString -string "$VERSION" "$APP/Contents/Info.plist"
plutil -replace CFBundleVersion -string "$VERSION" "$APP/Contents/Info.plist"
plutil -replace NSHighResolutionCapable -bool true "$APP/Contents/Info.plist"
plutil -replace LSMinimumSystemVersion -string "12.0" "$APP/Contents/Info.plist"
codesign --force --deep --sign - "$APP" >/dev/null 2>&1 || true   # ad-hoc signature for local use
echo "Built $APP ($(du -sh "$APP" | cut -f1))"

if [[ "${1:-}" == "--install" ]]; then
  mkdir -p "$HOME/Applications"
  rm -rf "$HOME/Applications/MissensExplorer.app"
  cp -R "$APP" "$HOME/Applications/"
  echo "Installed to ~/Applications/MissensExplorer.app"
fi

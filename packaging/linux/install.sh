#!/bin/sh
# Install MissensExplorer for the current user (no root): ~/.local/share/missensexplorer,
# with a menu entry and an icon. Run it from the unpacked folder:  ./install.sh
set -e
HERE=$(cd "$(dirname "$0")" && pwd)
DEST="$HOME/.local/share/missensexplorer"
rm -rf "$DEST"
mkdir -p "$DEST" "$HOME/.local/share/applications" "$HOME/.local/share/icons/hicolor/256x256/apps" "$HOME/.local/bin"
cp -R "$HERE/MissensExplorer/." "$DEST/"
cp "$HERE/missensexplorer.png" "$HOME/.local/share/icons/hicolor/256x256/apps/missensexplorer.png"
sed "s|/opt/missensexplorer/MissensExplorer|$DEST/MissensExplorer|" "$HERE/missensexplorer.desktop" \
  > "$HOME/.local/share/applications/missensexplorer.desktop"
ln -sf "$DEST/MissensExplorer" "$HOME/.local/bin/missensexplorer-gui"
command -v update-desktop-database >/dev/null 2>&1 && update-desktop-database "$HOME/.local/share/applications" || true
echo "MissensExplorer is installed: open it from the applications menu, or run missensexplorer-gui."
echo "To uninstall: rm -rf \"$DEST\" ~/.local/share/applications/missensexplorer.desktop ~/.local/bin/missensexplorer-gui"

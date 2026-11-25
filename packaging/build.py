"""Build the MissensExplorer desktop app and its double-click installer for this operating system.

    python -m pip install ".[gui]" -r packaging/requirements-installer.txt
    python packaging/build.py --installer [--smoke]

What comes out, in build/installers/:

- macOS:   MissensExplorer-<version>-macOS-<AppleSilicon|Intel>.dmg  (drag the app to Applications)
- Windows: MissensExplorer-<version>-Windows-Setup.exe             (Inno Setup; no admin rights needed)
- Linux:   missensexplorer_<version>_amd64.deb                      (Ubuntu/Debian: double-click or apt)
           MissensExplorer-<version>-Linux-x86_64.tar.gz            (other distributions: ./install.sh)

Everything the app needs is inside: Python, Qt (with the 3D viewer) and the libraries. Each
installer is built on its own system (PyInstaller does not cross-compile); the GitHub workflow
.github/workflows/installers.yml builds all of them.

--smoke runs the packaged app's self-test (packaging/launcher.py): it predicts BRAF V600E online,
draws the 3D structure and saves a screenshot. It uses a temporary cache and settings, never the
user's. --test-install also installs the installer as a user would (the .dmg copied out, the
Setup.exe run silently into a temporary folder, the .deb with apt and the .tar.gz with install.sh),
runs the self-test of the installed app and uninstalls it. On Linux it uses sudo (meant for CI).
"""

from __future__ import annotations

import argparse
import os
import platform
import shutil
import subprocess
import sys
import tarfile
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
PKG = ROOT / "packaging"
NAME = "MissensExplorer"
BUNDLE_ID = "org.missensexplorer.desktop"
MACOS_MIN = "13.0"  # Qt 6.11 (packaging/requirements-installer.txt)
EXCLUDES = ["streamlit", "altair", "pyarrow", "tornado", "matplotlib", "IPython", "pytest", "tkinter", "sklearn"]
# Linux libraries Qt's xcb platform and Qt WebEngine load from the system (Ubuntu/Debian names)
DEB_DEPENDS = [
    "libxcb-cursor0", "libxkbcommon-x11-0", "libxcb-icccm4", "libxcb-image0", "libxcb-keysyms1",
    "libxcb-randr0", "libxcb-render-util0", "libxcb-shape0", "libxcb-xinerama0", "libxcb-xkb1",
    "libegl1", "libgl1", "libfontconfig1", "libdbus-1-3", "libnss3", "libxcomposite1", "libxdamage1",
    "libxrandr2", "libxtst6", "libxkbfile1", "libasound2 | libasound2t64",
]  # fmt: skip


def system() -> str:
    """'macos', 'windows' or 'linux'."""
    return {"Darwin": "macos", "Windows": "windows"}.get(platform.system(), "linux")


def version() -> str:
    """The installed package's version (the build uses what is installed)."""
    import missensexplorer

    return missensexplorer.__version__


def run(cmd: list[str], **kw: object) -> None:
    print("+", " ".join(str(c) for c in cmd), flush=True)
    subprocess.run([str(c) for c in cmd], check=True, **kw)  # type: ignore[call-overload]


def pyinstaller(osname: str, work: Path) -> Path:
    """Freeze the app (one folder; a .app bundle on macOS); return the app's path."""
    args = [
        sys.executable, "-m", "PyInstaller", "--noconfirm", "--clean", "--windowed", "--name", NAME,
        "--distpath", work / "dist", "--workpath", work / "work", "--specpath", work,
        "--hidden-import", "PySide6.QtWebEngineWidgets", "--hidden-import", "PySide6.QtWebEngineCore",
        "--hidden-import", "PySide6.QtSvg",
    ]  # fmt: skip
    for m in EXCLUDES:
        args += ["--exclude-module", m]
    if osname == "macos":
        args += ["--icon", PKG / "macos" / "MissensExplorer.icns", "--osx-bundle-identifier", BUNDLE_ID]
    elif osname == "windows":
        args += ["--icon", windows_icon(work)]
    run([*args, PKG / "launcher.py"])
    return work / "dist" / (f"{NAME}.app" if osname == "macos" else NAME)


def windows_icon(work: Path) -> Path:
    """A multi-size .ico from the PNG icons."""
    from PIL import Image

    work.mkdir(parents=True, exist_ok=True)
    ico = work / "MissensExplorer.ico"
    Image.open(PKG / "icons" / "missensexplorer-256.png").save(
        ico, sizes=[(16, 16), (24, 24), (32, 32), (48, 48), (64, 64), (128, 128), (256, 256)]
    )
    return ico


# --------------------------------------------------------------------------- installers


def macos_dmg(app: Path, ver: str, out: Path) -> Path:
    """Set the bundle's version, sign it ad hoc, and pack it with an Applications link in a .dmg."""
    plist = app / "Contents" / "Info.plist"
    for key, kind, value in (
        ("CFBundleShortVersionString", "-string", ver),
        ("CFBundleVersion", "-string", ver),
        ("NSHighResolutionCapable", "-bool", "true"),
        ("LSMinimumSystemVersion", "-string", MACOS_MIN),
        ("CFBundleDisplayName", "-string", NAME),
    ):
        run(["plutil", "-replace", key, kind, value, plist])
    run(["codesign", "--force", "--deep", "--sign", "-", app])  # ad hoc: no Apple Developer ID
    arch = "AppleSilicon" if platform.machine() == "arm64" else "Intel"
    dmg = out / f"{NAME}-{ver}-macOS-{arch}.dmg"
    with tempfile.TemporaryDirectory() as tmp:
        stage = Path(tmp) / NAME
        stage.mkdir()
        run(["ditto", app, stage / app.name])
        (stage / "Applications").symlink_to("/Applications")
        shutil.copy(PKG / "INSTALL.txt", stage / "Read me first.txt")
        dmg.unlink(missing_ok=True)
        run(["hdiutil", "create", "-volname", f"{NAME} {ver}", "-srcfolder", stage, "-ov", "-format", "UDZO", dmg])
    return dmg


def windows_setup(folder: Path, ver: str, out: Path) -> Path:
    """Inno Setup installer: per-user by default (no admin rights), Start menu and desktop shortcuts."""
    iscc = os.environ.get("ISCC") or shutil.which("iscc") or r"C:\Program Files (x86)\Inno Setup 6\ISCC.exe"
    work = folder.parent.parent
    readme = work / "INSTALL.txt"  # Inno Setup reads a text file as UTF-8 only with a byte-order mark
    readme.write_text((PKG / "INSTALL.txt").read_text(encoding="utf-8"), encoding="utf-8-sig")
    run([
        iscc, f"/DAppVersion={ver}", f"/DSourceDir={folder}", f"/DOutputDir={out}",
        f"/DIconFile={windows_icon(work)}", f"/DReadme={readme}",
        PKG / "windows" / "MissensExplorer.iss",
    ])  # fmt: skip
    return out / f"{NAME}-{ver}-Windows-Setup.exe"


def linux_packages(folder: Path, ver: str, out: Path) -> list[Path]:
    """A .deb (Ubuntu/Debian) and a .tar.gz with an install script (other distributions)."""
    made = []
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp) / "deb"
        shutil.copytree(folder, root / "opt" / "missensexplorer", symlinks=True)
        (root / "usr" / "bin").mkdir(parents=True)
        (root / "usr" / "bin" / "missensexplorer-gui").symlink_to("/opt/missensexplorer/MissensExplorer")
        apps = root / "usr" / "share" / "applications"
        apps.mkdir(parents=True)
        shutil.copy(PKG / "linux" / "missensexplorer.desktop", apps)
        for size in (48, 64, 128, 256, 512):
            d = root / "usr" / "share" / "icons" / "hicolor" / f"{size}x{size}" / "apps"
            d.mkdir(parents=True)
            shutil.copy(PKG / "icons" / f"missensexplorer-{size}.png", d / "missensexplorer.png")
        (root / "DEBIAN").mkdir()
        size_kb = sum(f.stat().st_size for f in root.rglob("*") if f.is_file() and not f.is_symlink()) // 1024
        (root / "DEBIAN" / "control").write_text(
            "Package: missensexplorer\n"
            f"Version: {ver}\n"
            "Section: science\n"
            "Priority: optional\n"
            "Architecture: amd64\n"
            f"Installed-Size: {size_kb}\n"
            f"Depends: {', '.join(DEB_DEPENDS)}\n"
            "Maintainer: MissensExplorer developers <https://github.com/juvenalyosa/missensexplorer>\n"
            "Homepage: https://github.com/juvenalyosa/missensexplorer\n"
            "Description: Predict how DNA variants change the encoded human protein\n"
            " Desktop app: type variants (BRAF V600E, HGVS, rsIDs, coordinates) or open VCF/CSV/Excel\n"
            " files, and get the protein change, a missense-effect prediction, ACMG evidence,\n"
            " ClinVar, gnomAD, UniProt and the AlphaFold structure. Needs an internet connection.\n"
        )
        deb = out / f"missensexplorer_{ver}_amd64.deb"
        run(["dpkg-deb", "--build", "--root-owner-group", root, deb])
        made.append(deb)
    tgz = out / f"{NAME}-{ver}-Linux-x86_64.tar.gz"
    with tarfile.open(tgz, "w:gz") as tar:
        tar.add(folder, arcname=f"{NAME}-{ver}/{NAME}")
        tar.add(PKG / "linux" / "install.sh", arcname=f"{NAME}-{ver}/install.sh")
        tar.add(PKG / "linux" / "missensexplorer.desktop", arcname=f"{NAME}-{ver}/missensexplorer.desktop")
        tar.add(PKG / "icons" / "missensexplorer-256.png", arcname=f"{NAME}-{ver}/missensexplorer.png")
        tar.add(PKG / "INSTALL.txt", arcname=f"{NAME}-{ver}/README.txt")
    made.append(tgz)
    return made


# --------------------------------------------------------------------------- self-test


def test_installs(osname: str, made: list[Path], allow_no_3d: bool) -> int:
    """Install each installer as a user would, self-test the installed app, uninstall; 0 if all pass."""
    worst = 0
    for f in made:
        tmp = Path(tempfile.mkdtemp(prefix="mx-install-"))
        print(f"== install test: {f.name}", flush=True)
        if f.suffix == ".dmg":
            mnt = tmp / "mnt"
            run(["hdiutil", "attach", "-nobrowse", "-readonly", "-mountpoint", mnt, f])
            try:
                run(["ditto", mnt / f"{NAME}.app", tmp / "Applications" / f"{NAME}.app"])
            finally:
                run(["hdiutil", "detach", mnt])
            code = smoke(tmp / "Applications" / f"{NAME}.app", osname)
        elif f.suffix == ".exe":
            target = tmp / "MissensExplorer"
            run([f, "/VERYSILENT", "/SUPPRESSMSGBOXES", "/NORESTART", f"/DIR={target}"])  # per-user by default
            code = smoke(target, osname)
            run([target / "unins000.exe", "/VERYSILENT", "/SUPPRESSMSGBOXES", "/NORESTART"])
        elif f.suffix == ".deb":
            run(["sudo", "apt-get", "install", "-y", f.resolve()])  # a path (with a slash): apt installs the file
            code = smoke(Path("/opt/missensexplorer"), osname)
            run(["sudo", "apt-get", "remove", "-y", "missensexplorer"])
        elif f.name.endswith(".tar.gz"):
            home = tmp / "home"
            home.mkdir()
            with tarfile.open(f) as tar:
                tar.extractall(tmp, filter="tar")
            folder = next(p for p in tmp.iterdir() if p.name.startswith(NAME))
            run(["sh", folder / "install.sh"], env={**os.environ, "HOME": str(home)})
            code = smoke(home / ".local" / "share" / "missensexplorer", osname)
        else:
            continue
        if code == 3 and allow_no_3d:
            code = 0
        print(f"== install test {f.name}: {'ok' if code == 0 else f'FAILED ({code})'}", flush=True)
        worst = worst or code
    return worst


def smoke(app: Path, osname: str) -> int:
    """Run the packaged app's self-test with a temporary cache; return its exit code."""
    exe = app / "Contents" / "MacOS" / NAME if osname == "macos" else app / (f"{NAME}.exe" if osname == "windows" else NAME)
    tmp = Path(tempfile.mkdtemp(prefix="mx-smoke-"))
    shot = tmp / "smoke.png"
    env = {**os.environ, "MISSENSEXPLORER_SMOKE": str(shot), "MISSENSEXPLORER_CACHE": str(tmp / "cache")}
    print("+ self-test:", exe, flush=True)
    code = subprocess.run([str(exe)], env=env, timeout=900).returncode
    log = Path(str(shot) + ".log")
    print(log.read_text() if log.exists() else "(no log)", flush=True)
    print(f"self-test exit code {code}; screenshot {shot}", flush=True)
    return code


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--installer", action="store_true", help="also build the installer")
    ap.add_argument("--smoke", action="store_true", help="run the packaged app's self-test")
    ap.add_argument("--test-install", action="store_true", help="install each installer and self-test it")
    ap.add_argument(
        "--allow-no-3d",
        action="store_true",
        help="accept a self-test whose 3D page loaded but drew nothing (no WebGL on CI machines)",
    )
    args = ap.parse_args(argv)
    osname, ver = system(), version()
    work = ROOT / "build" / osname
    out = ROOT / "build" / "installers"
    out.mkdir(parents=True, exist_ok=True)
    app = pyinstaller(osname, work)
    made: list[Path] = []
    if args.installer:
        if osname == "macos":
            made.append(macos_dmg(app, ver, out))
        elif osname == "windows":
            made.append(windows_setup(app, ver, out))
        else:
            made += linux_packages(app, ver, out)
    for f in made:
        print(f"built {f.name} ({f.stat().st_size / 1e6:.0f} MB)")
    code = 0
    if args.smoke:
        code = smoke(app, osname)
        if code == 3 and args.allow_no_3d:
            print("the 3D page and 3Dmol loaded but drew nothing (no WebGL on this machine); everything else passed")
            code = 0
    if args.test_install and made:
        code = code or test_installs(osname, made, args.allow_no_3d)
    return code


if __name__ == "__main__":
    raise SystemExit(main())

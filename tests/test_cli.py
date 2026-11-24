"""The command-line entry points that need no network: ``cache`` and the launch of the web app."""

import sys

import pytest

from missensexplorer import cli
from missensexplorer.cache import Cache


@pytest.fixture
def cache_dir(tmp_path, monkeypatch):
    """A scratch cache directory (never the user's), not created yet."""
    d = tmp_path / "cache"
    monkeypatch.setenv("MISSENSEXPLORER_CACHE", str(d))
    return d


def test_cache_help_prints_help_and_leaves_the_cache_alone(cache_dir, capsys):
    with pytest.raises(SystemExit) as e:
        cli.main(["cache", "--help"])
    assert e.value.code == 0
    out = capsys.readouterr().out
    assert "usage: missensexplorer cache" in out
    assert "clear" in out
    assert not cache_dir.exists()


@pytest.mark.parametrize("argv", [["cache"], ["cache", "path"]])
def test_cache_prints_its_path(cache_dir, capsys, argv):
    assert cli.main(argv) == 0
    assert capsys.readouterr().out.strip() == str(cache_dir / "cache.sqlite")


def test_cache_rejects_an_unknown_action(cache_dir, capsys):
    with pytest.raises(SystemExit) as e:
        cli.main(["cache", "purge"])
    assert e.value.code == 2
    assert "invalid choice" in capsys.readouterr().err
    assert not cache_dir.exists()


def test_cache_clear_empties_the_cache(cache_dir, capsys):
    with Cache() as c:
        c.set("k", {"v": 1})
    assert cli.main(["cache", "clear"]) == 0
    assert "Cache cleared" in capsys.readouterr().out
    with Cache() as c:
        assert c.get("k") is None


def test_web_passes_the_streamlit_settings_as_flags(monkeypatch):
    pytest.importorskip("streamlit")
    from streamlit.web import cli as stcli

    seen: list[str] = []
    monkeypatch.setattr(sys, "argv", ["missensexplorer"])  # restored after the test
    monkeypatch.setattr(stcli, "main", lambda: seen.extend(sys.argv) or 0)
    assert cli.main(["web", "--port", "8765"]) == 0
    flags = dict(zip(seen[3::2], seen[4::2], strict=True))
    assert seen[:2] == ["streamlit", "run"]
    assert seen[2].endswith("webapp.py")
    assert flags["--server.port"] == "8765"
    assert flags["--server.maxUploadSize"] == "1024"
    assert flags["--client.toolbarMode"] == "minimal"


def test_web_leaves_out_settings_the_installed_streamlit_lacks(monkeypatch):
    pytest.importorskip("streamlit")
    from streamlit import config

    real = config.get_option

    def older_streamlit(name):
        if name.startswith(("theme.light.", "theme.dark.")):
            raise RuntimeError(f'Config key "{name}" not defined.')
        return real(name)

    monkeypatch.setattr(config, "get_option", older_streamlit)
    settings = cli._streamlit_settings()
    assert "server.maxUploadSize" in settings
    assert not any(k.startswith("theme.") for k in settings)

import builtins
import types
import pytest

# Import the module under test
import galfitools.shell.commands_mge as cli  # adjust to real module path


def test_mainMGE_calls_mge2gal_with_parsed_args(monkeypatch):
    called = {}

    def fake_mge2gal(
        GalfitFile,
        Ds9regFile,
        center,
        psf,
        twist,
        gauss,
        freeser,
        freesky,
        numgauss,
        xypos=None,
        ellip=None,
        posang=None,
    ):
        called.update(locals())  # capture arguments

    # Patch external calls
    monkeypatch.setattr(cli, "mge2gal", fake_mge2gal)
    monkeypatch.setattr(cli, "printWelcome", lambda: None)

    argv = [
        "galfit.01",
        "galaxy.reg",
        "--twist",
        "--center",
        "--psf",
        "1.5",
        "--gauss",
        "--freeser",
        "--freesky",
        "--numgauss",
        "7",
        "--xypos",
        "100",
        "120",
        "--ellip",
        "0.3",
        "--posang",
        "45",
    ]
    rc = cli.mainMGE(argv)
    assert rc == 0

    # Assert the parsed values were forwarded correctly
    assert called["GalfitFile"] == "galfit.01"
    assert called["Ds9regFile"] == "galaxy.reg"
    assert called["center"] is True
    assert called["psf"] == 1.5
    assert called["twist"] is True
    assert called["gauss"] is True
    assert called["freeser"] is True
    assert called["freesky"] is True
    assert called["numgauss"] == 7
    assert called["xypos"] == [100, 120]
    assert called["ellip"] == 0.3
    assert called["posang"] == 45


def test_mainMGE_missing_required_args_exits(monkeypatch):
    monkeypatch.setattr(cli, "printWelcome", lambda: None)
    with pytest.raises(SystemExit) as excinfo:
        cli.mainMGE([])  # no args
    assert excinfo.value.code != 0  # argparse exits with error


def test_mainSbProf_calls_sbProf(monkeypatch):
    captured = {}

    def fake_sbProf(ns):
        # store a few attributes to verify parsing
        captured["Image"] = ns.Image
        captured["Ds9Region"] = ns.Ds9Region
        captured["mgzpt"] = ns.mgzpt
        captured["plate"] = ns.plate
        captured["center"] = ns.center
        captured["logx"] = ns.logx

    monkeypatch.setattr(cli, "sbProf", fake_sbProf)
    monkeypatch.setattr(cli, "printWelcome", lambda: None)

    argv = [
        "img.fits",
        "galaxy.reg",
        "--mgzpt",
        "26.0",
        "--plate",
        "0.396",
        "--center",
        "--logx",
    ]
    rc = cli.mainSbProf(argv)
    assert rc == 0
    assert captured == {
        "Image": "img.fits",
        "Ds9Region": "galaxy.reg",
        "mgzpt": 26.0,
        "plate": 0.396,
        "center": True,
        "logx": True,
    }

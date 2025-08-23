import pytest
import galfitools.shell.commands_galout as cli


def test_mainPhotDs9(monkeypatch, capsys):
    monkeypatch.setattr(cli, "photDs9", lambda i, r, m, zp, sk: (22.22, 300.0))
    monkeypatch.setattr(cli, "printWelcome", lambda: None)
    rc = cli.mainPhotDs9(["img.fits", "reg.reg", "--zeropoint", "26", "--sky", "1.5"])
    assert rc == 0
    out = capsys.readouterr().out
    assert "exposition time" in out and "22.22" in out


def test_mainGetBreak(monkeypatch, capsys):
    monkeypatch.setattr(cli, "getBreak", lambda *a, **k: (12.5, 3, 45.0))

    class Head:
        scale = 0.4

    monkeypatch.setattr(
        cli, "Galfit", lambda f: type("G", (), {"ReadHead": lambda self: Head()})()
    )
    monkeypatch.setattr(cli, "printWelcome", lambda: None)
    rc = cli.mainGetBreak(["galfit.01", "--numinitial", "2", "--quick"])
    assert rc == 0
    out = capsys.readouterr().out
    assert "break radius is 12.50" in out


def test_mainKappa2(monkeypatch, capsys):
    monkeypatch.setattr(cli, "getKappa2", lambda *a: (7.0, 4, 30.0))

    class Head:
        scale = 0.5

    monkeypatch.setattr(
        cli, "Galfit", lambda f: type("G", (), {"ReadHead": lambda self: Head()})()
    )
    monkeypatch.setattr(cli, "printWelcome", lambda: None)
    rc = cli.mainKappa2(["galfit.01"])
    assert rc == 0
    out = capsys.readouterr().out
    assert 'kappa radius is 7.00 pixels or 3.50 "' in out


def test_mainGetMeRad(monkeypatch, capsys):
    monkeypatch.setattr(cli, "getMeRad", lambda *a: (20.0, 18.5, 17.0, 5, 10.0))

    class Head:
        scale = 1.0

    monkeypatch.setattr(
        cli, "Galfit", lambda f: type("G", (), {"ReadHead": lambda self: Head()})()
    )
    monkeypatch.setattr(cli, "printWelcome", lambda: None)
    rc = cli.mainGetMeRad(["g.01", "--rad", "3"])
    assert rc == 0
    out = capsys.readouterr().out
    assert "Surface brightness at this radius" in out and "17.00" in out


def test_mainShowCube(monkeypatch):
    seen = {}

    def fake_display(cubeimage, outimage, dpi, br, co, cmap, scale, noplot):
        seen.update(locals())

    monkeypatch.setattr(cli, "displayCube", fake_display)
    monkeypatch.setattr(cli, "printWelcome", lambda: None)
    rc = cli.mainShowCube(
        ["cube.fits", "--outimage", "x.png", "--dotsinch", "120", "--noplot"]
    )
    assert rc == 0
    assert (
        seen["cubeimage"] == "cube.fits"
        and seen["outimage"] == "x.png"
        and seen["noplot"] is True
    )


def test_mainFitlog2CSV(monkeypatch, capsys):
    called = {}
    monkeypatch.setattr(
        cli, "log2csv", lambda n, o, path=None: called.update(n=n, o=o, path=path)
    )
    monkeypatch.setattr(cli, "printWelcome", lambda: None)
    rc = cli.mainFitlog2CSV(["-o", "out.csv", "-n", "3", "-p", "/tmp"])
    assert rc == 0 and called == {"n": 3, "o": "out.csv", "path": "/tmp"}

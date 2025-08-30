import pytest
import re
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


def test_mainGetBreak2(monkeypatch, capsys):
    # Patch the correct function used by mainGetBreak2
    monkeypatch.setattr(cli, "getBreak2", lambda *a, **k: (12.5, 3, 45.0))

    # Stub Galfit only for the post-processing (arcsec conversion)
    class Head:
        scale = 0.4

    monkeypatch.setattr(
        cli, "Galfit", lambda f: type("G", (), {"ReadHead": lambda self: Head()})()
    )

    monkeypatch.setattr(cli, "printWelcome", lambda: None)

    rc = cli.mainGetBreak2(["galfit.01"])
    assert rc == 0

    out = capsys.readouterr().out
    assert "The break radius is 12.50" in out


def test_mainFWHM(monkeypatch, capsys):
    # Stub getFWHM so we don’t need real galfit files
    monkeypatch.setattr(cli, "getFWHM", lambda *a, **k: (8.5, 3, 42.0))

    # Fake Galfit.ReadHead().scale
    class Head:
        scale = 0.5

    monkeypatch.setattr(
        cli, "Galfit", lambda f: type("G", (), {"ReadHead": lambda self: Head()})()
    )

    # Silence welcome banner
    monkeypatch.setattr(cli, "printWelcome", lambda: None)

    # Call with minimal arguments
    rc = cli.mainFWHM(["galfit.01", "-d", "5", "-n", "2", "-pa", "37.5"])
    assert rc == 0

    out = capsys.readouterr().out
    # Validate key parts of output
    assert "number of model components:  3" in out
    assert "position angle: 42.00 degrees" in out
    assert "FWHM is 8.50 pixels or 4.25" in out  # 8.5 * 0.5 = 4.25


def test_mainKappa(monkeypatch, capsys):
    monkeypatch.setattr(cli, "getKappa", lambda *a: (7.0, 4, 30.0))

    class Head:
        scale = 0.5

    monkeypatch.setattr(
        cli, "Galfit", lambda f: type("G", (), {"ReadHead": lambda self: Head()})()
    )
    monkeypatch.setattr(cli, "printWelcome", lambda: None)

    rc = cli.mainKappa(["galfit.01"])
    assert rc == 0

    out = capsys.readouterr().out
    norm = " ".join(out.split()).lower()  # collapse spaces, lower-case
    assert 'kappa radius is 7.00 pixels or 3.50 "' in norm


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
    norm = " ".join(out.split()).lower()
    assert 'kappa radius is 7.00 pixels or 3.50 "' in norm


def test_mainGetReComp_default_half(monkeypatch, capsys):
    """
    fracrad defaults to 0.5 → the 'Mean Surface Brightness at effective radius' line should print.
    """
    # Stub heavy computation
    # Returns: EffRad, totmag, meanme, me, N, theta
    monkeypatch.setattr(
        cli, "getReComp", lambda *a, **k: (12.0, 19.25, 21.10, 20.50, 3, 37.0)
    )

    # Fake Galfit.ReadHead().scale
    class Head:
        scale = 0.4  # arcsec / pixel

    monkeypatch.setattr(
        cli, "Galfit", lambda f: type("G", (), {"ReadHead": lambda self: Head()})()
    )

    # Silence banner
    monkeypatch.setattr(cli, "printWelcome", lambda: None)

    # Minimal CLI args (use defaults: dis=10, fracrad=0.5, numcomp=1, angle=None)
    rc = cli.mainGetReComp(["galfit.01"])
    assert rc == 0

    out = capsys.readouterr().out
    # Key checks
    assert "number of model components:  3" in out
    assert "position angle: 37.00 degrees" in out
    assert "Total Magnitude of the galaxy: 19.25" in out
    assert "Surface brightness at radius of 50% of light" in out
    assert "(μr): 20.50 mag/''" in out
    # Because fracrad == 0.5, this line should be present
    assert "Mean Surface Brightness at effective radius (<μ>e): 21.10 mag/''" in out
    # 12.0 px × 0.4 = 4.80 arcsec
    assert 'The radius at 50% of light is 12.00 pixels or 4.80 "' in out


def test_mainGetReComp_nonhalf_fraction(monkeypatch, capsys):
    """
    fracrad != 0.5 → the 'Mean Surface Brightness at effective radius' line should NOT print.
    """
    monkeypatch.setattr(
        cli, "getReComp", lambda *a, **k: (8.0, 18.00, 22.00, 21.50, 2, 10.0)
    )

    class Head:
        scale = 0.5

    monkeypatch.setattr(
        cli, "Galfit", lambda f: type("G", (), {"ReadHead": lambda self: Head()})()
    )
    monkeypatch.setattr(cli, "printWelcome", lambda: None)

    # Explicitly set a different fraction (e.g., 0.3 → 30%)
    rc = cli.mainGetReComp(
        ["galfit.01", "-fr", "0.3", "-d", "7", "-n", "2", "-pa", "10"]
    )
    assert rc == 0

    out = capsys.readouterr().out
    assert "number of model components:  2" in out
    assert "position angle: 10.00 degrees" in out
    assert "Total Magnitude of the galaxy: 18.00" in out
    assert "Surface brightness at radius of 30% of light" in out
    assert "(μr): 21.50 mag/''" in out
    # Because fracrad != 0.5, this line should be absent
    assert "Mean Surface Brightness at effective radius (<μ>e)" not in out
    # 8.0 px × 0.5 = 4.00 arcsec
    assert 'The radius at 30% of light is 8.00 pixels or 4.00 "' in out


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


def test_maingetSlope_default(monkeypatch, capsys):
    """Uses all defaults: dis=10, numcomp=1, angle=None, slope=0.5, no plot/ranx."""
    # Stub heavy computation
    # Returns: rgam, N, theta
    monkeypatch.setattr(cli, "getSlope", lambda *a, **k: (15.0, 3, 25.0))

    # Fake Galfit.ReadHead().scale
    class Head:
        scale = 0.4  # arcsec/pixel

    monkeypatch.setattr(
        cli, "Galfit", lambda f: type("G", (), {"ReadHead": lambda self: Head()})()
    )

    # Silence banner
    monkeypatch.setattr(cli, "printWelcome", lambda: None)

    rc = cli.maingetSlope(["galfit.01"])
    assert rc == 0

    out = capsys.readouterr().out
    assert "number of model components:  3" in out
    assert "position angle: 25.00 degrees" in out
    # 15.0 px * 0.4 = 6.0 arcsec
    assert 'The radius with slope 0.50 is 15.00 pixels, or 6.00 "' in out


def test_maingetSlope_with_options(monkeypatch, capsys):
    """Explicit options: dis, numcomp, angle, slope, plot flag, and ranx."""
    monkeypatch.setattr(cli, "getSlope", lambda *a, **k: (8.0, 2, 10.0))

    class Head:
        scale = 0.5

    monkeypatch.setattr(
        cli, "Galfit", lambda f: type("G", (), {"ReadHead": lambda self: Head()})()
    )
    monkeypatch.setattr(cli, "printWelcome", lambda: None)

    rc = cli.maingetSlope(
        [
            "galfit.01",
            "-d",
            "7",
            "-n",
            "2",
            "-a",
            "35",
            "-s",
            "0.75",
            "--plot",
            "-rx",
            "5",
            "50",
        ]
    )
    assert rc == 0

    out = capsys.readouterr().out
    assert "number of model components:  2" in out
    assert "position angle: 10.00 degrees" in out
    # 8.0 px * 0.5 = 4.0 arcsec
    assert 'The radius with slope 0.75 is 8.00 pixels, or 4.00 "' in out


def test_maingetN_default(monkeypatch, capsys):
    """Defaults: dis=10, numcomp=1, radfrac=0.2, angle=None, no plot, const=0."""
    # getN returns: sersic, meanser, stdser, totmag, N, theta
    monkeypatch.setattr(cli, "getN", lambda *a, **k: (3.20, 3.10, 0.25, 18.75, 4, 22.0))
    monkeypatch.setattr(cli, "printWelcome", lambda: None)

    rc = cli.maingetN(["galfit.01"])
    assert rc == 0

    out = capsys.readouterr().out
    # Be tolerant to spacing/newlines
    norm = " ".join(out.split())
    assert "number of model components: 4" in norm
    assert "position angle: 22.00 degrees" in norm
    assert "Total Magnitude of the galaxy: 18.75" in norm
    assert "Mean Surface Brightness at effective radius: 3.20" in norm
    assert "Sersic index mean: 3.10 Standard deviation: 0.25" in norm


def test_maingetN_with_options(monkeypatch, capsys):
    """Explicit options exercise flags: -d, -n, -pa, -rf, -p, -c."""
    monkeypatch.setattr(cli, "getN", lambda *a, **k: (2.50, 2.45, 0.10, 19.00, 2, 10.0))
    monkeypatch.setattr(cli, "printWelcome", lambda: None)

    rc = cli.maingetN(
        [
            "galfit.01",
            "-d",
            "7",
            "-n",
            "2",
            "-pa",
            "10",
            "-rf",
            "0.3",
            "--plot",
            "-c",
            "0.15",
        ]
    )
    assert rc == 0

    out = capsys.readouterr().out
    norm = " ".join(out.split())
    assert "number of model components: 2" in norm
    assert "position angle: 10.00 degrees" in norm
    assert "Total Magnitude of the galaxy: 19.00" in norm
    # Spot-check both index lines
    assert "Mean Surface Brightness at effective radius: 2.50" in norm
    assert "Sersic index mean: 2.45 Standard deviation: 0.10" in norm


def test_maingetBT(monkeypatch, capsys):
    """Simple path: bulge/total + total magnitude printed."""
    monkeypatch.setattr(cli, "getBT", lambda *a, **k: (0.37, 17.80, 3))
    monkeypatch.setattr(cli, "printWelcome", lambda: None)

    rc = cli.maingetBT(["galfit.01", "-d", "6", "-n", "3"])
    assert rc == 0

    out = capsys.readouterr().out
    norm = " ".join(out.split())
    assert "number of model components: 3" in norm
    assert "Total Magnitude of the galaxy: 17.80" in norm
    assert "Bulge to total luminosity ratio: 0.37" in norm


def test_mainMissingLight_with_options(monkeypatch, capsys):
    # Stub: returns (magmiss, N1, N2)
    monkeypatch.setattr(cli, "getMissLight", lambda *a, **k: (0.85, 5, 7))
    monkeypatch.setattr(cli, "printWelcome", lambda: None)

    rc = cli.mainMissingLight(
        [
            "coreless.01",
            "core.01",
            "15",  # rad (float)
            "-d",
            "6",
            "-n",
            "2",
        ]
    )
    assert rc == 0

    out = capsys.readouterr().out
    # Tolerant checks (spacing/newlines may vary)
    assert re.search(r"number of model components coreless model:\s+5", out)
    assert re.search(r"number of model components core model:\s+7", out)
    assert re.search(r"the missing light is\s+0\.85\s+mag", out, re.IGNORECASE)


def test_mainMissingLight_defaults(monkeypatch, capsys):
    # Defaults: dis=10, numcomp=1
    monkeypatch.setattr(cli, "getMissLight", lambda *a, **k: (1.23, 2, 3))
    monkeypatch.setattr(cli, "printWelcome", lambda: None)

    rc = cli.mainMissingLight(["g1.01", "g2.01", "20"])
    assert rc == 0

    out = capsys.readouterr().out
    assert re.search(r"coreless model:\s+2", out)
    assert re.search(r"core model:\s+3", out)
    assert re.search(r"missing light is\s+1\.23\s+mag", out, re.IGNORECASE)


def test_maingetCOW_defaults(monkeypatch, capsys):
    """Defaults: dis=10, plotfile=cow.png, fracrad=0.95, numcomp=1, dpi=100."""
    # Stub: returns (totmag, N, theta)
    monkeypatch.setattr(cli, "getCOW", lambda *a, **k: (18.75, 3, 25.0))
    monkeypatch.setattr(cli, "printWelcome", lambda: None)

    rc = cli.maingetCOW(["galfit.01"])
    assert rc == 0

    out = capsys.readouterr().out
    norm = " ".join(out.split())
    assert "number of model components: 3" in norm
    assert "position angle: 25.00 degrees" in norm
    assert "Total Magnitude of the galaxy: 18.75" in norm
    assert "plot file: cow.png" in norm


def test_maingetPeak_defaults(monkeypatch, capsys):
    """No optional flags: center=False, mask=None."""
    # Stub getPeak → returns X, Y, q, PA
    monkeypatch.setattr(cli, "getPeak", lambda *a, **k: (123.4, 56.7, 0.82, 37.5))
    # Silence banner
    monkeypatch.setattr(cli, "printWelcome", lambda: None)

    rc = cli.maingetPeak(["image.fits", "ellipse.reg"])
    assert rc == 0

    out = capsys.readouterr().out
    # Tolerant checks (whitespace/newlines can vary)
    assert re.search(r"peak is at \(x,\s*y\)\s*=\s*123\.4\s+56\.7", out)
    assert re.search(r"axis ratio q\s*=\s*0\.82", out)
    assert re.search(r"position angle\s*=\s*37\.50", out)
    # assert "position angle = 37.50 degrees" in out


def test_maingetPeak_with_flags(monkeypatch, capsys):
    """Use --center and a mask file to exercise optional args."""
    monkeypatch.setattr(cli, "getPeak", lambda *a, **k: (200.0, 150.0, 0.60, 12.0))
    monkeypatch.setattr(cli, "printWelcome", lambda: None)

    rc = cli.maingetPeak(
        [
            "image.fits",
            "ellipse.reg",
            "--center",
            "-m",
            "mask.fits",
        ]
    )
    assert rc == 0

    out = capsys.readouterr().out
    assert re.search(r"peak is at \(x,\s*y\)\s*=\s*200\.0\s+150\.0", out)
    assert re.search(r"axis ratio q\s*=\s*0\.60", out)
    assert re.search(r"position angle\s*=\s*12\.00", out)
    # assert "position angle = 12.00 degrees" in out


def test_mainGetBulgeRad(monkeypatch, capsys):
    """Two files; ensure both component counts and radius conversion print."""
    # getBulgeRad returns: rbulge, N1, N2, theta
    monkeypatch.setattr(cli, "getBulgeRad", lambda *a, **k: (12.5, 4, 6, 33.0))

    # Fake Galfit.ReadHead().scale for arcsec conversion
    class Head:
        scale = 0.4  # arcsec/pixel

    monkeypatch.setattr(
        cli, "Galfit", lambda f: type("G", (), {"ReadHead": lambda self: Head()})()
    )

    monkeypatch.setattr(cli, "printWelcome", lambda: None)

    rc = cli.mainGetBulgeRad(
        [
            "coreless.01",
            "core.01",
            "-d",
            "8",
            "-n",
            "2",
            "-pa",
            "33",
            "--plot",
            "-rx",
            "5",
            "100",
        ]
    )
    assert rc == 0

    out = capsys.readouterr().out
    # Use regex to be resilient to spacing/newlines
    assert re.search(r"number of model components for the bulge:\s+4", out)
    assert re.search(r"number of model components for the rest of the galaxy:\s+6", out)
    assert re.search(r"position angle:\s+33\.00 degrees", out)
    # 12.5 px * 0.4 = 5.00 arcsec
    assert re.search(r"The bulge radius is 12\.50 pixels or 5\.00 \"", out)


def test_mainGetBarSize(monkeypatch, capsys):
    # stub heavy functions
    monkeypatch.setattr(cli, "getBarSize", lambda *a: (20.0, 1, 17))

    class Head:
        scale = 1.0

    monkeypatch.setattr(
        cli, "Galfit", lambda f: type("G", (), {"ReadHead": lambda self: Head()})()
    )
    monkeypatch.setattr(cli, "printWelcome", lambda: None)

    # minimal valid args (no ranx, no plot)
    rc = cli.mainGetBarSize(["g.01", "-d", "5", "-n", "1", "-o", "bar.reg"])
    assert rc == 0
    out = capsys.readouterr().out
    assert "bar size is 20.00" in out

    # with plot flag and ranx range
    rc = cli.mainGetBarSize(
        ["g.01", "-d", "5", "-n", "1", "-o", "bar.reg", "--plot", "-rx", "10", "100"]
    )
    assert rc == 0


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

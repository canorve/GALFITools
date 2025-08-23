import pytest

# adjust to your real import path
import galfitools.shell.commands_sky as cli


def test_mainGalfitSky_calls_galfitSky(monkeypatch):
    called = {}

    def fake_galfitSky(image, mask, mgzpt, scale, x, y, initsky):
        called.update(
            image=image, mask=mask, mgzpt=mgzpt, scale=scale, x=x, y=y, initsky=initsky
        )

    monkeypatch.setattr(cli, "galfitSky", fake_galfitSky)
    monkeypatch.setattr(cli, "printWelcome", lambda: None)

    argv = [
        "img.fits",
        "mask.fits",
        "--mgzpt",
        "26.0",
        "--scale",
        "0.4",
        "--xpos",
        "10",
        "--ypos",
        "12",
        "--initsky",
        "5.5",
    ]
    rc = cli.mainGalfitSky(argv)
    assert rc == 0
    assert called == dict(
        image="img.fits",
        mask="mask.fits",
        mgzpt=26.0,
        scale=0.4,
        x=10.0,
        y=12.0,
        initsky=5.5,
    )


def test_mainGalfitSky_missing_args(monkeypatch):
    monkeypatch.setattr(cli, "printWelcome", lambda: None)
    with pytest.raises(SystemExit):
        cli.mainGalfitSky([])


def test_mainSkyDs9_calls_SkyDs9_and_prints(monkeypatch, capsys):
    def fake_SkyDs9(ImageFile, RegFile, mask, outliers):
        assert ImageFile == "img.fits"
        assert RegFile == "regions.reg"
        assert mask == "mask.fits"
        assert outliers is True
        return 12.3456, 0.789

    monkeypatch.setattr(cli, "SkyDs9", fake_SkyDs9)
    monkeypatch.setattr(cli, "printWelcome", lambda: None)

    rc = cli.mainSkyDs9(
        ["img.fits", "regions.reg", "--mask", "mask.fits", "--outliers"]
    )
    assert rc == 0
    out = capsys.readouterr().out
    assert "mean sky: 12.346" in out
    assert "std sky: 0.789" in out


def test_mainSkyDs9_requires_positionals(monkeypatch):
    monkeypatch.setattr(cli, "printWelcome", lambda: None)
    with pytest.raises(SystemExit):
        cli.mainSkyDs9([])


def test_mainSkyRing_calls_SkyRing_and_prints(monkeypatch, capsys):
    def fake_SkyRing(Image, mask, Ds9regFile, width, center, outliers):
        assert Image == "img.fits"
        assert mask is None
        assert Ds9regFile == "ell.reg"
        assert width == 25
        assert center is True
        assert outliers is False
        return 10.0, 2.0, 9.5, 30.25

    monkeypatch.setattr(cli, "SkyRing", fake_SkyRing)
    monkeypatch.setattr(cli, "printWelcome", lambda: None)

    rc = cli.mainSkyRing(["img.fits", "ell.reg", "--width", "25", "--center"])
    assert rc == 0
    out = capsys.readouterr().out
    assert "Major axis of ellipse is used as initial radius." in out
    assert "mean =  10.00" in out
    assert "std=2.00" in out
    assert "median = 9.50" in out
    assert "radius 30.25" in out


def test_mainSkyRing_requires_positionals(monkeypatch):
    monkeypatch.setattr(cli, "printWelcome", lambda: None)
    with pytest.raises(SystemExit):
        cli.mainSkyRing([])

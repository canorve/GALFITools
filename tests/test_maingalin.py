import pytest

import galfitools.shell.commands_galin as cli
import numpy as np
import re


def test_mainGetStar(monkeypatch, capsys):
    called = {}

    def fake_getStar(image, regfile, imsize, center, sky, imout, sigma, sigout):
        called.update(locals())

    monkeypatch.setattr(cli, "getStar", fake_getStar)
    monkeypatch.setattr(cli, "printWelcome", lambda: None)

    argv = [
        "img.fits",
        "star.reg",
        "64",
        "--center",
        "--sky",
        "10.0",
        "--out",
        "cut.fits",
        "--sigma",
        "sig.fits",
    ]
    rc = cli.mainGetStar(argv)
    assert rc == 0
    assert called["image"] == "img.fits"
    assert called["regfile"] == "star.reg"
    assert called["imsize"] == 64
    assert called["center"] is True
    assert called["sky"] == 10.0
    assert called["imout"] == "cut.fits"
    assert called["sigma"] == "sig.fits"
    out = capsys.readouterr().out
    assert "Object fits file: cut.fits" in out


def test_mainInitGal(monkeypatch):
    captured = {}

    def fake_InitGal(GalfitFile, number, p3, p4, p5, p6, p7, p8, p9, p10, numcomp):
        captured.update(file=GalfitFile, number=number, p3=p3, numcomp=numcomp)

    monkeypatch.setattr(cli, "InitGal", fake_InitGal)
    monkeypatch.setattr(cli, "printWelcome", lambda: None)

    argv = ["galfit.01", "-n", "3", "--param3", "10", "20", "--numcomp", "2"]
    rc = cli.mainInitGal(argv)
    assert rc == 0
    assert captured["file"] == "galfit.01"
    assert captured["number"] == 3
    assert captured["p3"] == [10.0, 20.0]
    assert captured["numcomp"] == 2


def test_mainMaskDs9(monkeypatch):
    called = {}

    def fake_maskDs9(
        MaskFile, RegFile, fill, image, border, borValue, skymean=None, skystd=None
    ):
        called.update(locals())

    monkeypatch.setattr(cli, "maskDs9", fake_maskDs9)
    monkeypatch.setattr(cli, "printWelcome", lambda: None)

    argv = [
        "mask.fits",
        "regions.reg",
        "--fill",
        "5",
        "--image",
        "img.fits",
        "--border",
        "--borValue",
        "2.0",
        "--skymean",
        "1.1",
        "--skystd",
        "0.9",
    ]
    rc = cli.mainMaskDs9(argv)
    assert rc == 0
    assert called["MaskFile"] == "mask.fits"
    assert called["RegFile"] == "regions.reg"
    assert called["fill"] == 5
    assert called["image"] == "img.fits"
    assert called["border"] is True
    assert called["borValue"] == 2.0
    assert called["skymean"] == 1.1
    assert called["skystd"] == 0.9


def test_mainMaskSky(monkeypatch, capsys):
    got = {}

    def fake_skyRem(image, mask, sky_mean, sky_sig, nsig, borValue, bor_flag):
        got.update(locals())

    def fake_maskDs9(*args, **kwargs):
        got["maskDs9_called"] = True

    monkeypatch.setattr(cli, "skyRem", fake_skyRem)
    monkeypatch.setattr(cli, "maskDs9", fake_maskDs9)
    monkeypatch.setattr(cli, "printWelcome", lambda: None)

    argv = [
        "img.fits",
        "mask.fits",
        "--skymean",
        "0.5",
        "--skysigma",
        "0.2",
        "--numbersig",
        "3",
        "--borValue",
        "1.0",
        "--border",
        "--region",
        "remove.reg",
    ]
    rc = cli.mainMaskSky(argv)
    assert rc == 0
    assert got["image"] == "img.fits"
    assert got["mask"] == "mask.fits"
    assert got["sky_mean"] == 0.5
    assert got["sky_sig"] == 0.2
    assert got["nsig"] == 3.0
    assert got["borValue"] == 1.0
    assert got["bor_flag"] is True
    assert got.get("maskDs9_called") is True


def test_mainxy2fits(monkeypatch):
    class FakeXY:
        def MakeFits(self, ImageFile, AsciiFile, Value):
            self.args = (ImageFile, AsciiFile, Value)

    fx = FakeXY()
    monkeypatch.setattr(cli, "xy2fits", lambda: fx)
    monkeypatch.setattr(cli, "printWelcome", lambda: None)

    rc = cli.mainxy2fits(["img.fits", "pts.txt", "--val", "7"])
    assert rc == 0
    assert fx.args == ("img.fits", "pts.txt", 7)


def test_maincheckFile(monkeypatch, capsys):
    class Head:
        inputimageflag = False
        inputimage = "A"
        sigimageflag = False
        sigimage = "B"
        psfimageflag = True
        convxflag = True
        convyflag = True
        maskimageflag = True
        constraintsflag = True
        xsizeflag = True
        ysizeflag = True

    galax = np.array([1, 1, 2])
    mag = np.array([20.0, 21.0, 22.0])

    def fake_checkFile(path, dis):
        return Head(), galax, mag, 7

    monkeypatch.setattr(cli, "checkFile", fake_checkFile)
    monkeypatch.setattr(cli, "printWelcome", lambda: None)

    rc = cli.maincheckFile(["galfit.01"])
    assert rc == 0
    out = capsys.readouterr().out
    assert "File A not found" in out
    assert "GALFIT will create a sigma file" in out
    assert "Total number of model components:" in out
    assert "Total number of free parameters: 7" in out


def test_mainGetBoxSizeDs9(monkeypatch, capsys):
    monkeypatch.setattr(cli, "getBoxSizeDs9", lambda f: (1, 2, 3, 4))
    monkeypatch.setattr(cli, "printWelcome", lambda: None)
    rc = cli.mainGetBoxSizeDs9(["box.reg"])
    assert rc == 0
    out = capsys.readouterr().out
    assert "1 2 3 4" in out


@pytest.mark.parametrize(
    "argv, expected",
    [
        # 1) Minimal: only required args → all defaults
        (
            ["img.fits", "ell.reg"],
            # image, regfile, center, mask, zeropoint, sky, noprint, bulgetot, bards9
            ("img.fits", "ell.reg", False, None, 25.0, 0.0, False, None, None),
        ),
        # 2) Typical flags: custom zeropoint/sky + center + mask + noprint
        (
            [
                "img.fits",
                "ell.reg",
                "--zeropoint",
                "26",
                "--sky",
                "1.0",
                "--center",
                "-m",
                "mask.fits",
                "--noprint",
            ],
            ("img.fits", "ell.reg", True, "mask.fits", 26.0, 1.0, True, None, None),
        ),
        # 3) Bulge/Disk with a bar region: requires --bulgetot and --bards9
        (
            ["img.fits", "ell.reg", "--bulgetot", "0.35", "--bards9", "bar.reg"],
            ("img.fits", "ell.reg", False, None, 25.0, 0.0, False, 0.35, "bar.reg"),
        ),
        # 4) Everything at once (different values)
        (
            [
                "img.fits",
                "ell.reg",
                "-zp",
                "24.5",
                "-sk",
                "0.2",
                "-c",
                "-n",
                "-m",
                "m.fits",
                "-bt",
                "0.6",
                "-b",
                "bar2.reg",
            ],
            ("img.fits", "ell.reg", True, "m.fits", 24.5, 0.2, True, 0.6, "bar2.reg"),
        ),
    ],
)
def test_maingetSersic(monkeypatch, argv, expected):
    seen = {}

    def fake_getSersic(
        image, regfile, center, mask, zeropoint, sky, noprint, bulgetot, bards9
    ):
        seen["args"] = (
            image,
            regfile,
            center,
            mask,
            zeropoint,
            sky,
            noprint,
            bulgetot,
            bards9,
        )

    monkeypatch.setattr(cli, "getSersic", fake_getSersic)
    monkeypatch.setattr(cli, "printWelcome", lambda: None)

    rc = cli.maingetSersic(argv)
    assert rc == 0
    assert seen["args"] == expected


def test_main_imarith(monkeypatch):
    seen = {}

    def fake_imarith(ImageFile, output, image2, add, mul, div, sub):
        seen.update(locals())

    monkeypatch.setattr(cli, "imarith", fake_imarith)
    monkeypatch.setattr(cli, "printWelcome", lambda: None)

    rc = cli.main_imarith(["img.fits", "--output", "out.fits", "--add", "2.5"])
    assert rc == 0
    assert seen["ImageFile"] == "img.fits"
    assert seen["output"] == "out.fits"
    assert seen["add"] == 2.5
    assert seen["mul"] is None and seen["div"] is None and seen["sub"] is None


def test_mainMakePSF(monkeypatch):
    seen = {}

    def fake_makePSF(
        GalfitFile, image, Ds9regFile, center, out, sigma, twist, numgauss
    ):
        seen.update(locals())

    monkeypatch.setattr(cli, "makePSF", fake_makePSF)
    monkeypatch.setattr(cli, "printWelcome", lambda: None)

    rc = cli.mainMakePSF(
        [
            "img.fits",
            "galfit.01",
            "star.reg",
            "--center",
            "--out",
            "psf.fits",
            "--twist",
            "--numgauss",
            "5",
        ]
    )
    assert rc == 0
    assert seen["image"] == "img.fits"
    assert seen["GalfitFile"] == "galfit.01"
    assert seen["Ds9regFile"] == "star.reg"
    assert seen["center"] is True
    assert seen["out"] == "psf.fits"
    assert seen["twist"] is True
    assert seen["numgauss"] == 5


def test_mainMakeMask_defaults(monkeypatch, capsys):
    """Call with only required args; defaults for maskout, satds9, scale."""
    called = {}

    def fake_makeMask(sexfile, imagefile, maskout, scale, satds9):
        called.update(
            {
                "sexfile": sexfile,
                "imagefile": imagefile,
                "maskout": maskout,
                "scale": scale,
                "satds9": satds9,
            }
        )

    monkeypatch.setattr(cli, "makeMask", fake_makeMask)
    monkeypatch.setattr(cli, "printWelcome", lambda: None)

    rc = cli.mainMakeMask(["catalog.sex", "image.fits"])
    assert rc == 0

    out = capsys.readouterr().out
    assert "Done. Mask image created" in out
    # Verify defaults passed correctly
    assert called["maskout"] == "masksex.fits"
    assert called["satds9"] == "ds9sat.reg"
    assert called["scale"] == 1


def test_mainMakeMask_with_options(monkeypatch, capsys):
    """Call with explicit overrides for optional args."""
    called = {}
    monkeypatch.setattr(
        cli,
        "makeMask",
        lambda sexfile, imagefile, maskout, scale, satds9: called.update(
            dict(
                sexfile=sexfile,
                imagefile=imagefile,
                maskout=maskout,
                scale=scale,
                satds9=satds9,
            )
        ),
    )
    monkeypatch.setattr(cli, "printWelcome", lambda: None)

    rc = cli.mainMakeMask(
        [
            "cat.sex",
            "img.fits",
            "-o",
            "custom_mask.fits",
            "-sf",
            "custom_sat.reg",
            "-s",
            "2.5",
        ]
    )
    assert rc == 0

    out = capsys.readouterr().out
    assert re.search(r"Done\. Mask image created", out)

    # Ensure options were passed through
    assert called["maskout"] == "custom_mask.fits"
    assert called["satds9"] == "custom_sat.reg"
    assert called["scale"] == 2.5

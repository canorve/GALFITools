import pytest

import galfitools.shell.commands_sex as cli
import re

import shutil
from pathlib import Path


def test_mainMakeMask_defaults(monkeypatch, capsys):
    """Call with only required args; defaults for maskout, satds9, scale."""
    region_dir = Path("kron_regions")
    called = {}

    def fake_makeMask(sexfile, imagefile, maskout, scale, satds9, regdir, regid):
        called.update(
            {
                "sexfile": sexfile,
                "imagefile": imagefile,
                "maskout": maskout,
                "scale": scale,
                "satds9": satds9,
                "regdir": regdir,
                "regid": regid,
            }
        )

    monkeypatch.setattr(cli, "makeMask", fake_makeMask)
    monkeypatch.setattr(cli, "printWelcome", lambda: None)

    try:

        rc = cli.mainMakeMask(["catalog.sex", "image.fits"])

        assert rc == 0

        out = capsys.readouterr().out
        assert "Done. Mask image created" in out
        # Verify defaults passed correctly
        assert called["maskout"] == "masksex.fits"
        assert called["satds9"] == "ds9sat.reg"
        assert called["scale"] == 1

    finally:
        if region_dir.exists():
            shutil.rmtree(region_dir)


def test_mainMakeMask_with_options(monkeypatch, capsys):
    """Call with explicit overrides for optional args."""
    region_dir = Path("kron_regions")
    called = {}
    monkeypatch.setattr(
        cli,
        "makeMask",
        lambda sexfile, imagefile, maskout, scale, satds9, region_dir, region_id: called.update(
            dict(
                sexfile=sexfile,
                imagefile=imagefile,
                maskout=maskout,
                scale=scale,
                satds9=satds9,
                region_dir=region_dir,
                region_id=region_id,
            )
        ),
    )
    monkeypatch.setattr(cli, "printWelcome", lambda: None)

    try:
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
                "--region_dir",
                "kron_regions",
                "--region_id",
                "1",
            ]
        )
        assert rc == 0

        out = capsys.readouterr().out
        assert re.search(r"Done\. Mask image created", out)

        # Ensure options were passed through
        assert called["maskout"] == "custom_mask.fits"
        assert called["satds9"] == "custom_sat.reg"
        assert called["scale"] == 2.5
        assert called["region_dir"] == "kron_regions"
        assert called["region_id"] == 1

    finally:
        if region_dir.exists():
            shutil.rmtree(region_dir)

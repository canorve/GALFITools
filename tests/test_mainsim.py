import pytest

# adjust import to your actual module path
import galfitools.cli_make_sim as cli


def test_mainMakeSim_calls_makeSim_with_args(monkeypatch):
    called = {}

    def fake_makeSim(image, gain, skymean, skystd, newimage):
        called.update(
            image=image, gain=gain, skymean=skymean, skystd=skystd, newimage=newimage
        )

    monkeypatch.setattr(cli, "makeSim", fake_makeSim)
    monkeypatch.setattr(cli, "printWelcome", lambda: None)

    argv = [
        "model.fits",
        "sim.fits",
        "--sky",
        "12.3",
        "--std",
        "0.9",
        "--gain",
        "2.0",
    ]
    rc = cli.mainMakeSim(argv)
    assert rc == 0

    assert called == {
        "image": "model.fits",
        "gain": 2.0,
        "skymean": 12.3,
        "skystd": 0.9,
        "newimage": "sim.fits",
    }


def test_mainMakeSim_requires_positional_args(monkeypatch):
    monkeypatch.setattr(cli, "printWelcome", lambda: None)
    # Missing required positionals should trigger argparse SystemExit
    with pytest.raises(SystemExit):
        cli.mainMakeSim([])

    # Provide only one positional -> still error
    with pytest.raises(SystemExit):
        cli.mainMakeSim(["model.fits"])

import galfitools.shell.prt as prt


def test_printWelcome(capsys):
    # Call the function
    prt.printWelcome()

    # Capture output
    captured = capsys.readouterr()

    # Check that key strings appear in the output
    assert "GALFITools: a library for GALFIT" in captured.out
    assert "Version:" in captured.out
    assert "https://github.com/canorve/GALFITools" in captured.out

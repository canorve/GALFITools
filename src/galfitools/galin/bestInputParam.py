import argparse
from pathlib import Path

# Import your real function here
from galfitools.galout.getChiNu import getChiNu


def read_galfile_list(list_file):
    """
    Read a text file containing one GALFIT file per line.

    Empty lines and lines starting with '#' are ignored.
    """
    galfiles = []

    with open(list_file, "r", encoding="utf-8") as file:
        for line in file:
            galfile = line.strip()

            if not galfile:
                continue

            if galfile.startswith("#"):
                continue

            galfiles.append(galfile)

    return galfiles


def run_get_chinu_for_list(list_file, numcomp, fracrad, regfile, delete):
    """
    Apply getChiNu to each GALFIT file in the input list.

    Returns
    -------
    results : list of dict
        List with the results for each GALFIT file.
    best_result : dict
        Result with the lowest chi_nu value.
    """
    galfiles = read_galfile_list(list_file)

    results = []

    for galfile in galfiles:
        chinu, aic, bic, numparfree = getChiNu(
            galfile,
            numcomp,
            fracrad,
            regfile,
            delete,
        )

        result = {
            "galfile": galfile,
            "chinu": chinu,
            "aic": aic,
            "bic": bic,
            "numparfree": numparfree,
        }

        results.append(result)

    if not results:
        raise ValueError("The input list does not contain valid GALFIT files.")

    best_result = min(results, key=lambda item: item["chinu"])

    return results, best_result


def parse_args():
    """
    Parse command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description="Find the GALFIT file with the lowest chi_nu value."
    )

    parser.add_argument(
        "list_file",
        help="Text file containing one GALFIT file per line.",
    )

    parser.add_argument(
        "-n",
        "--numcomp",
        type=int,
        default=1,
        required=True,
        help="Component Number",
    )

    parser.add_argument(
        "-r",
        "--regfile",
        required=True,
        help="DS9 ellipse region file.",
    )

    parser.add_argument(
        "-d",
        "--delete",
        action="store_true",
        help="deletes the sigma image used to compute chinu",
    )

    return parser.parse_args()


def mainBestInputParam():
    """
    Main program.
    """
    args = parse_args()

    fracrad = 0.99
    results, best_result = run_get_chinu_for_list(
        args.list_file,
        args.numcomp,
        fracrad,
        args.regfile,
        args.delete,
    )

    chinu_list = [item["chinu"] for item in results]

    print("Chi_nu values:")
    print(chinu_list)

    print()
    print("Best GALFIT file:")
    print(f"galfile = {best_result['galfile']}")
    print(f"chi_nu  = {best_result['chinu']}")
    print(f"AIC     = {best_result['aic']}")
    print(f"BIC     = {best_result['bic']}")
    print(f"Nfree   = {best_result['numparfree']}")


if __name__ == "__main__":
    mainBestInputParam()

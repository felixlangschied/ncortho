"""
ncOrtho - Targeted ortholog search for miRNAs
Copyright (C) 2021 Felix Langschied

ncOrtho is a free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

ncOrtho is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with ncOrtho.  If not, see <http://www.gnu.org/licenses/>.
"""

# Create a core set of orthologs
# Find the corresponding syntenic regions in reference and core species
# Search for core orthologs by reciprocal BLAST search
# Create Stockholm structural alignment

import argparse
import logging
import multiprocessing as mp
import os
import sys
from typing import Dict, List, Optional, Tuple

from pyfiglet import Figlet
from importlib.metadata import version

try:
    from createcm import create_cm, create_phmm
    from core_reblast import blastsearch
    from locate_position import categorize_mirna_position
    from synteny import analyze_synteny
    import coreset_utils as cu
except ModuleNotFoundError:
    from ncOrtho.coreset.createcm import create_cm, create_phmm
    from ncOrtho.coreset.core_reblast import blastsearch
    from ncOrtho.coreset.locate_position import categorize_mirna_position
    from ncOrtho.coreset.synteny import analyze_synteny
    import ncOrtho.coreset.coreset_utils as cu


# ---------------------------------------------------------------------------
# Logging setup
# ---------------------------------------------------------------------------
logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Custom argparse helpers
# ---------------------------------------------------------------------------
VALID_ID_TYPES = ("ID", "Name", "GeneID", "gene_id", "CDS")


def _positive_int(value: str) -> int:
    """Argparse type: strictly positive integer."""
    try:
        ivalue = int(value)
    except ValueError:
        raise argparse.ArgumentTypeError(f"'{value}' is not an integer")
    if ivalue < 1:
        raise argparse.ArgumentTypeError(f"'{value}' must be a positive integer (>= 1)")
    return ivalue


def _yes_no(value: str) -> bool:
    """Argparse type: convert 'yes'/'no' strings to bool."""
    lower = value.strip().lower()
    if lower == "yes":
        return True
    if lower == "no":
        return False
    raise argparse.ArgumentTypeError(
        f"Invalid value '{value}'. Expected 'yes' or 'no'."
    )


def _valid_id_type(value: str) -> str:
    """Argparse type: validate --idtype against allowed values."""
    if value not in VALID_ID_TYPES:
        raise argparse.ArgumentTypeError(
            f"Invalid ID type '{value}'. Must be one of: {', '.join(VALID_ID_TYPES)}"
        )
    return value


def _existing_file(value: str) -> str:
    """Argparse type: validate that a file exists."""
    path = os.path.realpath(value)
    if not os.path.isfile(path):
        raise argparse.ArgumentTypeError(f"File not found: {path}")
    return path


# ---------------------------------------------------------------------------
# Argument parser
# ---------------------------------------------------------------------------
def build_parser() -> argparse.ArgumentParser:
    """Build and return the argument parser for the coreset pipeline."""
    parser = argparse.ArgumentParser(
        prog="ncOrtho-coreset",
        description="Build covariance models of reference miRNAs from a core set of orthologs.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser._action_groups.pop()
    required = parser.add_argument_group("Required arguments")
    optional = parser.add_argument_group("Optional arguments")

    # --- Required -----------------------------------------------------------
    required.add_argument(
        "-p", "--parameters",
        metavar="<path>",
        type=_existing_file,
        required=True,
        help="Path to the parameters file in YAML format.",
    )
    required.add_argument(
        "-n", "--ncrna",
        metavar="<path>",
        type=_existing_file,
        required=True,
        help="Path to tab-separated file of reference miRNA information.",
    )
    required.add_argument(
        "-o", "--output",
        metavar="<path>",
        type=str,
        required=True,
        help="Path for the output folder.",
    )

    # --- Optional -----------------------------------------------------------
    max_cpu = mp.cpu_count()
    optional.add_argument(
        "--threads",
        metavar="N",
        type=_positive_int,
        default=max_cpu,
        help=f"Number of CPU cores to use (default: all available = {max_cpu}).",
    )
    optional.add_argument(
        "--mgi",
        metavar="N",
        type=_positive_int,
        default=3,
        help="Maximum number of gene insertions in the core species (default: 3).",
    )
    optional.add_argument(
        "--dust",
        metavar="yes/no",
        type=_yes_no,
        default=False,
        help='Use BLASTn dust filter; skips repeat regions when enabled (default: no).',
    )
    optional.add_argument(
        "--create_model",
        metavar="yes/no",
        type=_yes_no,
        default=True,
        help='Set to "no" to only create the alignment without building a model (default: yes).',
    )
    optional.add_argument(
        "--phmm",
        metavar="yes/no",
        type=_yes_no,
        default=False,
        help='Create a profile HMM instead of a covariance model (default: no).',
    )
    optional.add_argument(
        "--rcoffee",
        metavar="yes/no",
        type=_yes_no,
        default=True,
        help='Use r_coffee for alignment; set to "no" for default t_coffee (default: yes).',
    )
    optional.add_argument(
        "--redo",
        metavar="yes/no",
        type=_yes_no,
        default=True,
        help='Rebuild models even if they already exist at the output destination (default: yes).',
    )
    optional.add_argument(
        "--idtype",
        metavar="TYPE",
        type=_valid_id_type,
        default="ID",
        choices=VALID_ID_TYPES,
        help=(
            "ID field in the reference GFF to match against the pairwise-orthologs file. "
            f"Options: {', '.join(VALID_ID_TYPES)} (default: ID)."
        ),
    )
    optional.add_argument(
        "--max_anchor_dist",
        metavar="N",
        type=_positive_int,
        default=3,
        help=(
            "Number of additional genes left and right of the reference miRNA "
            "to consider as syntenic anchors (default: 3)."
        ),
    )
    optional.add_argument(
        "--verbose",
        action="store_true",
        default=False,
        help="Print additional diagnostic information.",
    )

    return parser


# ---------------------------------------------------------------------------
# Validation helpers
# ---------------------------------------------------------------------------
def validate_args(args: argparse.Namespace) -> None:
    """Run cross-field validation that argparse types alone cannot cover."""
    available_cpu = mp.cpu_count()
    if args.threads > available_cpu:
        raise SystemExit(
            f"Error: Requested {args.threads} threads, but only {available_cpu} are available."
        )


# ---------------------------------------------------------------------------
# Pipeline steps (thin wrappers for clarity)
# ---------------------------------------------------------------------------
def setup_output_dirs(output: str) -> str:
    """Create the output and tmp directories. Return the tmp path."""
    os.makedirs(output, exist_ok=True)
    tmpout = os.path.join(output, "tmp")
    os.makedirs(tmpout, exist_ok=True)
    return tmpout


def load_pairwise_orthologs(
    core_dict: dict, ref_annotation: str, idtype: str
) -> dict:
    """Load pairwise orthologs, dispatching on idtype."""
    logger.info("Reading pairwise orthologs")
    if idtype == "CDS":
        return cu.pairwise_orthologs_from_cds(core_dict, ref_annotation)
    return cu.read_pairwise_orthologs(core_dict)


def locate_mirna_positions(
    mirnas: list,
    ref_dict: dict,
    ortho_dict: dict,
    max_anchor_dist: int,
    verbose: bool,
) -> Dict[str, Tuple]:
    """Determine the genomic position category and core orthologs for each miRNA."""
    logger.info("Determining the position of each miRNA and its neighboring gene(s)")
    mirna_positions: Dict[str, Tuple] = {}
    for mirna in mirnas:
        mirid, chromo, start, end, strand = cu.mirna_position(mirna)
        syntype, core_orthos = categorize_mirna_position(
            mirid, chromo, start, end, strand,
            ref_dict, ortho_dict, max_anchor_dist, verbose,
        )
        if not syntype:
            logger.warning("Could not localize %s", mirid)
            continue
        mirna_positions[mirid] = (syntype, core_orthos)
    return mirna_positions


def write_synteny_fastas(
    syntenyregion_per_mirna: dict, tmpout: str
) -> None:
    """Write per-miRNA FASTA files for syntenic regions."""
    for mirid, fastalist in syntenyregion_per_mirna.items():
        if not fastalist:
            logger.warning(
                "No syntenic region found in any core species for %s. "
                "Verify that annotation IDs match the orthologs file.",
                mirid,
            )
            continue
        miroutdir = os.path.join(tmpout, mirid)
        os.makedirs(miroutdir, exist_ok=True)
        fasta_path = os.path.join(miroutdir, f"synteny_regions_{mirid}.fa")
        with open(fasta_path, "w") as fh:
            for line in fastalist:
                fh.write(line)


def build_models(
    mirnas: list,
    ref_genome: str,
    tmpout: str,
    output: str,
    cpu: int,
    dust: bool,
    verbose: bool,
    rcoffee: bool,
    create_model: bool,
    phmm: bool,
    redo: bool,
) -> None:
    """Run reciprocal BLAST, alignment, and model construction for each miRNA."""
    logger.info("Collecting core orthologs and training models")

    # Map bool back to the string expected by downstream functions that
    # have not yet been refactored.
    dust_str = "yes" if dust else "no"
    rcoffee_str = "yes" if rcoffee else "no"

    for mirna in mirnas:
        mirid = mirna[0]
        sto_path = blastsearch(mirna, ref_genome, tmpout, cpu, dust_str, verbose, rcoffee_str)

        if not create_model or sto_path is None:
            continue

        if phmm:
            logger.info("Constructing pHMM for %s", mirid)
            create_phmm(sto_path, output, mirid, cpu)
        else:
            model_out = os.path.join(output, f"{mirid}.cm")
            if not redo and os.path.isfile(model_out):
                logger.info(
                    "Model for %s already exists at %s — skipping (--redo is off).",
                    mirid, model_out,
                )
                continue
            logger.info("Constructing covariance model for %s", mirid)
            create_cm(sto_path, output, mirid, cpu)


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------
def main() -> None:
    """Entry point for the ncOrtho coreset pipeline."""
    # Print banner
    custom_fig = Figlet(font="stop")
    print(custom_fig.renderText("ncOrtho")[:-3], flush=True)
    v = version("ncOrtho")
    print(f"Version: {v}\n", flush=True)

    parser = build_parser()

    # Show help when called without arguments.
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()

    # Configure logging
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(levelname)s: %(message)s",
    )

    validate_args(args)

    # Resolve paths
    output = os.path.realpath(args.output)
    tmpout = setup_output_dirs(output)

    # Load parameters
    core_dict, ref_paths, all_paths = cu.parse_yaml(args.parameters)
    missing = [p for p in all_paths if not os.path.isfile(p)]
    if missing:
        raise FileNotFoundError(
            "The following required files are missing:\n  " + "\n  ".join(missing)
        )

    # Pipeline
    ortho_dict = load_pairwise_orthologs(core_dict, ref_paths["annotation"], args.idtype)

    logger.info("Reading miRNA data")
    mirnas = cu.read_mirnas(args.ncrna)
    if not mirnas:
        raise SystemExit("Error: No miRNAs loaded from the input file.")

    logger.info("Reading reference annotation")
    ref_dict = cu.parse_annotation(ref_paths["annotation"], args.idtype)

    mirna_positions = locate_mirna_positions(
        mirnas, ref_dict, ortho_dict, args.max_anchor_dist, args.verbose,
    )

    logger.info("Identifying syntenic regions in core species")
    syntenyregion_per_mirna = analyze_synteny(
        core_dict, mirna_positions, tmpout, args.idtype, args.mgi, args.verbose,
    )
    if not syntenyregion_per_mirna:
        raise SystemExit("Error: No regions of conserved synteny found in any core species.")

    write_synteny_fastas(syntenyregion_per_mirna, tmpout)

    build_models(
        mirnas=mirnas,
        ref_genome=ref_paths["genome"],
        tmpout=tmpout,
        output=output,
        cpu=args.threads,
        dust=args.dust,
        verbose=args.verbose,
        rcoffee=args.rcoffee,
        create_model=args.create_model,
        phmm=args.phmm,
        redo=args.redo,
    )

    logger.info("Construction of core set finished")


if __name__ == "__main__":
    main()
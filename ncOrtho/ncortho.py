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

# Find orthologs of reference miRNAs in the genome of a query species.
# Pipeline: CM/pHMM search -> candidate extraction -> reverse BLAST verification

import argparse
import logging
import multiprocessing as mp
import os
import sys
from time import time
from typing import Dict, List, Optional, Tuple

from tqdm import tqdm
from pyfiglet import Figlet
from importlib.metadata import version

try:
    from blastparser import ReBlastParser
    from genparser import GenomeParser
    from cmsearch import model_search
    from utils import check_blastdb, make_blastndb, write_output, str2bool
    from reblast import perform_reblast
    from rna_object import rna_maker
except ModuleNotFoundError:
    from ncOrtho.blastparser import ReBlastParser
    from ncOrtho.genparser import GenomeParser
    from ncOrtho.cmsearch import model_search
    from ncOrtho.utils import check_blastdb, make_blastndb, write_output, str2bool
    from ncOrtho.reblast import perform_reblast
    from ncOrtho.rna_object import rna_maker


logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Custom argparse types
# ---------------------------------------------------------------------------
def _positive_int(value: str) -> int:
    """Argparse type: strictly positive integer."""
    try:
        ivalue = int(value)
    except ValueError:
        raise argparse.ArgumentTypeError(f"'{value}' is not an integer")
    if ivalue < 1:
        raise argparse.ArgumentTypeError(f"'{value}' must be a positive integer (>= 1)")
    return ivalue


def _ratio_float(value: str) -> float:
    """Argparse type: float in range [0.0, 1.0] (inclusive)."""
    try:
        fvalue = float(value)
    except ValueError:
        raise argparse.ArgumentTypeError(f"'{value}' is not a number")
    if not 0.0 <= fvalue <= 1.0:
        raise argparse.ArgumentTypeError(f"'{value}' must be between 0.0 and 1.0")
    return fvalue


def _non_negative_float(value: str) -> float:
    """Argparse type: float >= 0."""
    try:
        fvalue = float(value)
    except ValueError:
        raise argparse.ArgumentTypeError(f"'{value}' is not a number")
    if fvalue < 0:
        raise argparse.ArgumentTypeError(f"'{value}' must be >= 0")
    return fvalue


def _existing_file(value: str) -> str:
    """Argparse type: validate that a file exists."""
    path = os.path.realpath(value)
    if not os.path.isfile(path):
        raise argparse.ArgumentTypeError(f"File not found: {path}")
    return path


def _existing_dir(value: str) -> str:
    """Argparse type: validate that a directory exists."""
    path = os.path.realpath(value)
    if not os.path.isdir(path):
        raise argparse.ArgumentTypeError(f"Directory not found: {path}")
    return path


def _optional_positive_int(value: str) -> Optional[int]:
    """Argparse type: positive integer or None (from 'None'/'none'/empty)."""
    if value is None or value.strip().lower() == "none":
        return None
    return _positive_int(value)


def _dust_choice(value: str) -> str:
    """Argparse type: validate dust filter value."""
    cleaned = value.strip().lower()
    if cleaned not in ("yes", "no"):
        raise argparse.ArgumentTypeError(
            f"Invalid value '{value}' for --dust. Expected 'yes' or 'no'."
        )
    return cleaned


# ---------------------------------------------------------------------------
# Argument parser
# ---------------------------------------------------------------------------
def build_parser() -> argparse.ArgumentParser:
    """Build and return the argument parser."""
    parser = argparse.ArgumentParser(
        prog="ncOrtho",
        description="Find orthologs of reference miRNAs in the genome of a query species.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser._action_groups.pop()
    required = parser.add_argument_group("Required arguments")
    optional = parser.add_argument_group("Optional arguments")

    # --- Required -----------------------------------------------------------
    required.add_argument(
        "-m", "--models",
        metavar="<path>",
        type=_existing_dir,
        required=True,
        help="Path to directory containing covariance models (.cm).",
    )
    required.add_argument(
        "-n", "--ncrna",
        metavar="<path>",
        type=_existing_file,
        required=True,
        help="Path to tab-separated file with reference miRNA information.",
    )
    required.add_argument(
        "-o", "--output",
        metavar="<path>",
        type=str,
        required=True,
        help="Path to the output directory.",
    )
    required.add_argument(
        "-q", "--query",
        metavar="<.fa>",
        type=_existing_file,
        required=True,
        help="Path to query genome in FASTA format.",
    )
    required.add_argument(
        "-r", "--reference",
        metavar="<.fa>",
        type=_existing_file,
        required=True,
        help="Path to reference genome in FASTA format.",
    )

    # --- Optional -----------------------------------------------------------
    max_cpu = mp.cpu_count()
    optional.add_argument(
        "--queryname",
        metavar="STR",
        type=str,
        default="",
        help="Name for the output directory (recommended). "
             "Defaults to the query filename without extension.",
    )
    optional.add_argument(
        "--cpu",
        metavar="N",
        type=_positive_int,
        default=3,
        help=f"Number of CPU cores to use (default: 3).",
    )
    optional.add_argument(
        "--cm_cutoff",
        metavar="FLOAT",
        type=_ratio_float,
        default=0.5,
        help="CMsearch bitscore cutoff as ratio of the CM score against the "
             "reference species (default: 0.5).",
    )
    optional.add_argument(
        "--minlength",
        metavar="FLOAT",
        type=_ratio_float,
        default=0.7,
        help="Minimum length of a CMsearch hit as ratio of the reference "
             "pre-miRNA length (default: 0.7).",
    )
    optional.add_argument(
        "--heuristic",
        type=str2bool,
        metavar="True/False",
        default=True,
        help="Perform a preliminary BLAST search to narrow candidate regions "
             "for CMsearch. Greatly improves speed (default: True).",
    )
    optional.add_argument(
        "--heur_blast_evalue",
        type=_non_negative_float,
        metavar="FLOAT",
        default=0.5,
        help="E-value filter for the heuristic BLASTn search "
             "(default: 0.5; set to 10 to effectively disable).",
    )
    optional.add_argument(
        "--heur_blast_length",
        type=_ratio_float,
        metavar="FLOAT",
        default=0.5,
        help="Length cutoff for heuristic BLASTn as ratio of reference "
             "pre-miRNA length (default: 0.5; set to 0 to disable).",
    )
    optional.add_argument(
        "--sensitive_heuristic",
        type=str2bool,
        metavar="True/False",
        default=False,
        help="If no candidate region is found via BLASTn, search with CM in "
             "the full query genome (default: False).",
    )
    optional.add_argument(
        "--cleanup",
        type=str2bool,
        metavar="True/False",
        default=True,
        help="Remove temporary files after completion (default: True).",
    )
    optional.add_argument(
        "--phmm",
        type=str2bool,
        metavar="True/False",
        default=False,
        help="Use pHMM instead of CM for ortholog search (default: False).",
    )
    optional.add_argument(
        "--refblast",
        type=str,
        metavar="<path>",
        default="",
        help="Path to BLASTdb of the reference species.",
    )
    optional.add_argument(
        "--queryblast",
        type=str,
        metavar="<path>",
        default="",
        help="Path to BLASTdb of the query species.",
    )
    optional.add_argument(
        "--maxcmhits",
        type=_optional_positive_int,
        metavar="N",
        default=None,
        help="Maximum number of CMsearch hits to examine. Useful when "
             "reference miRNA is in a repeat region (default: no limit).",
    )
    optional.add_argument(
        "--dust",
        metavar="yes/no",
        type=_dust_choice,
        default="no",
        help="Use BLASTn dust filter during re-BLAST. Reduces runtime for "
             "miRNAs in repeat regions, but will skip those miRNAs (default: no).",
    )
    optional.add_argument(
        "--checkCoorthologsRef",
        type=str2bool,
        metavar="True/False",
        default=False,
        help="If the re-BLAST best hit is not the original reference miRNA, "
             "check whether it is a likely co-ortholog (default: False).",
    )

    return parser


# ---------------------------------------------------------------------------
# Validation and setup helpers
# ---------------------------------------------------------------------------
def validate_args(args: argparse.Namespace) -> None:
    """Cross-field validation that argparse types alone cannot cover."""
    available_cpu = mp.cpu_count()
    if args.cpu > available_cpu:
        raise SystemExit(
            f"Error: Requested {args.cpu} threads, but only {available_cpu} are available."
        )


def resolve_query_name(args: argparse.Namespace) -> str:
    """Derive a query name from --queryname or the query filename."""
    if args.queryname:
        return args.queryname
    basename = os.path.basename(args.query)
    # Strip final extension: "species.genome.fa" -> "species.genome"
    name, _ = os.path.splitext(basename)
    return name


def setup_output(output_root: str, qname: str) -> str:
    """Create output directory tree. Return final output path."""
    output = os.path.realpath(output_root)
    if os.path.basename(output) != qname:
        output = os.path.join(output, qname)
    os.makedirs(output, exist_ok=True)
    return output


def setup_data_dir(output: str, qname: str, query: str) -> str:
    """Create data directory with a symlink to the query genome. Return symlink path."""
    q_data = os.path.join(output, "data")
    os.makedirs(q_data, exist_ok=True)
    qlink = os.path.join(q_data, f"{qname}.fa")
    if not os.path.islink(qlink):
        os.symlink(query, qlink)
    return qlink


def resolve_blastdb(
    label: str,
    user_path: str,
    fallback_path: str,
    data_dir_path: str,
    source_genome: str,
    create_if_missing: bool = True,
) -> str:
    """
    Resolve a BLASTdb path with fallback logic:
      1. Use user-supplied path (validate it exists).
      2. Check if a db exists at fallback_path (e.g. alongside the genome).
      3. Create one at data_dir_path from source_genome.
    """
    if user_path:
        if not check_blastdb(user_path):
            raise SystemExit(f"Error: {label} BLASTdb not found at: {user_path}")
        return user_path

    if check_blastdb(fallback_path):
        return fallback_path

    if create_if_missing:
        if not check_blastdb(data_dir_path):
            logger.info("Creating %s BLASTdb", label)
            make_blastndb(source_genome, data_dir_path)
        return data_dir_path

    return fallback_path


# ---------------------------------------------------------------------------
# Logging setup
# ---------------------------------------------------------------------------
def setup_logging(output: str, qname: str, verbose: bool = False) -> None:
    """Configure file and console logging."""
    logout = os.path.join(output, f"{qname}_extended.log")
    # Remove stale log from previous run
    if os.path.isfile(logout):
        os.remove(logout)

    ncortho_logger = logging.getLogger("ncortho")
    ncortho_logger.setLevel(logging.DEBUG)

    fh = logging.FileHandler(logout)
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(logging.Formatter("%(levelname)s: %(message)s"))
    ncortho_logger.addHandler(fh)


# ---------------------------------------------------------------------------
# Core ortholog search
# ---------------------------------------------------------------------------
def search_orthologs(
    mirna_dict: dict,
    cm_cutoff: float,
    cpu: int,
    msl: float,
    models: str,
    qlink: str,
    qblast: str,
    output: str,
    cleanup: bool,
    heuristic: Tuple,
    phmm: bool,
    refblast: str,
    dust: str,
    check_coorth_ref: bool,
    max_hits: Optional[int],
) -> Tuple[List[str], List[str], int]:
    """
    Run the CM/pHMM search + re-BLAST verification loop for all miRNAs.

    Returns:
        mirna_orthologs: FASTA lines for verified orthologs.
        shortlog: Per-miRNA summary lines.
        success_count: Number of miRNAs with at least one verified ortholog.
    """
    ncortho_logger = logging.getLogger("ncortho")
    qname = os.path.basename(output)
    shortlog: List[str] = ["# miRNA\tSeconds\tStatus\n"]
    mirna_orthologs: List[str] = []
    success_count = 0

    for mirna in tqdm(mirna_dict, file=sys.stdout):
        start_time = time()
        restricted = False
        sys.stdout.flush()
        mirna_data = mirna_dict[mirna]

        # miRNA objects for which no CM was found are empty
        if not mirna_data:
            elapsed = f"{time() - start_time:.3f}"
            shortlog.append(f"{mirna}\t{elapsed}\tNo CM\n")
            continue

        ncortho_logger.info("\n### %s", mirna)

        # Create per-miRNA output folder when not cleaning up
        if not heuristic[0] or not cleanup:
            outdir = os.path.join(output, mirna)
        else:
            outdir = output
        os.makedirs(outdir, exist_ok=True)

        # CM/pHMM search
        cm_results, exitstatus = model_search(
            mirna_data, cm_cutoff, cpu, msl, models, qlink, qblast,
            outdir, cleanup, heuristic, phmm,
        )

        if not cm_results:
            ncortho_logger.info("# %s for %s", exitstatus, mirna)
            elapsed = f"{time() - start_time:.3f}"
            shortlog.append(f"{mirna}\t{elapsed}\t{exitstatus}\n")
            continue

        # Optionally restrict to top N hits
        if max_hits and len(cm_results) > max_hits:
            ncortho_logger.warning(
                "Maximum CMsearch hits reached. Restricting to best %d hits.", max_hits
            )
            cm_results = dict(list(cm_results.items())[:max_hits])
            restricted = True

        # Extract candidate sequences
        gp = GenomeParser(qlink, cm_results.values())
        candidates = gp.extract_sequences()
        ncortho_logger.info(
            "Covariance model search successful, found %d ortholog candidate(s).",
            len(candidates),
        )

        # Reverse BLAST verification
        reblast_hits: Dict[str, str] = {}
        for candidate, sequence in candidates.items():
            # Skip near-entirely masked sequences
            if sequence.count("N") >= 0.9 * len(sequence):
                continue
            reblasthit = perform_reblast(
                sequence, refblast, cpu, outdir, candidate,
                mirna_data, dust, msl, cleanup, check_coorth_ref,
            )
            if reblasthit:
                reblast_hits[candidate] = reblasthit

        if not reblast_hits:
            elapsed = f"{time() - start_time:.3f}"
            ncortho_logger.info("None of the candidates for %s could be verified.", mirna)
            if restricted:
                shortlog.append(
                    f"{mirna}\t{elapsed}\tNo re-BLAST after restricting to {max_hits} CMsearch hits\n"
                )
            else:
                shortlog.append(f"{mirna}\t{elapsed}\tNo re-BLAST\n")
            continue

        # Evaluate co-orthologs if multiple hits
        nr_accepted = len(reblast_hits)
        if nr_accepted == 1:
            ncortho_logger.info("ncOrtho found 1 verified ortholog.")
            out_dict = reblast_hits
        else:
            ncortho_logger.info(
                "ncOrtho found %d potential co-orthologs. Evaluating distances.",
                nr_accepted,
            )
            rbp = ReBlastParser(mirna_data, reblast_hits)
            out_dict = rbp.verify_coorthologs(outdir)
            ncortho_logger.info(
                "ncOrtho found %d verified %s.",
                len(out_dict),
                "ortholog" if len(out_dict) == 1 else "co-orthologs",
            )

        # Collect output FASTA entries
        for hit, seq in out_dict.items():
            cmres = list(cm_results[hit])
            cmres.insert(0, qname)
            header = "|".join(str(entry) for entry in cmres)
            mirna_orthologs.append(f">{header}\n{seq}\n")

        elapsed = f"{time() - start_time:.3f}"
        shortlog.append(f"{mirna}\t{elapsed}\tSUCCESS\n")
        success_count += 1

    return mirna_orthologs, shortlog, success_count


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------
def main() -> None:
    """Entry point for the ncOrtho ortholog search pipeline."""
    # Print banner
    custom_fig = Figlet(font="stop")
    print(custom_fig.renderText("ncOrtho")[:-3], flush=True)
    v = version("ncOrtho")
    print(f"Version: {v}\n", flush=True)

    parser = build_parser()

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()
    validate_args(args)

    # Resolve names and paths
    qname = resolve_query_name(args)
    output = setup_output(args.output, qname)
    qlink = setup_data_dir(output, qname, args.query)

    # Resolve BLASTdbs
    q_data = os.path.join(output, "data")
    refname = os.path.basename(args.reference)

    refblast = resolve_blastdb(
        label="Reference",
        user_path=args.refblast,
        fallback_path=args.reference,
        data_dir_path=os.path.join(q_data, refname),
        source_genome=args.reference,
    )

    qblast = ""
    if args.heuristic:
        qblast = resolve_blastdb(
            label="Query",
            user_path=args.queryblast,
            fallback_path=qlink,
            data_dir_path=qlink,
            source_genome=args.query,
        )

    # Build heuristic tuple expected by downstream functions
    heuristic = (
        args.heuristic,
        args.heur_blast_evalue,
        args.heur_blast_length,
        args.sensitive_heuristic,
    )

    # Setup logging
    setup_logging(output, qname)

    print(f"# Starting ncOrtho run for {qname}\n", flush=True)

    # Create miRNA objects
    mirna_dict = rna_maker(args.ncrna, args.models, args.phmm, args.cpu)

    # Run ortholog search
    print("\n# Ortholog search", flush=True)
    mirna_orthologs, shortlog, success_count = search_orthologs(
        mirna_dict=mirna_dict,
        cm_cutoff=args.cm_cutoff,
        cpu=args.cpu,
        msl=args.minlength,
        models=args.models,
        qlink=qlink,
        qblast=qblast,
        output=output,
        cleanup=args.cleanup,
        heuristic=heuristic,
        phmm=args.phmm,
        refblast=refblast,
        dust=args.dust,
        check_coorth_ref=args.checkCoorthologsRef,
        max_hits=args.maxcmhits,
    )

    # Write results
    write_output(mirna_orthologs, os.path.join(output, f"{qname}_orthologs.fa"))
    write_output(shortlog, os.path.join(output, f"{qname}_summary.log"))

    print(
        f"\n### ncOrtho found orthologs for {success_count} out of "
        f"{len(mirna_dict)} ncRNAs",
        flush=True,
    )


if __name__ == "__main__":
    main()
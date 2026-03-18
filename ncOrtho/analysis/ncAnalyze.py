"""
ncOrtho - Targeted ortholog search for miRNAs
Copyright (C) 2021 Felix Langschied

ncOrtho is free software: you can redistribute it and/or modify
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

# Analyse ncOrtho results:
#   - Write PhyloProfile input
#   - Build per-miRNA alignments (MUSCLE 5)
#   - Construct a supermatrix and infer a species tree (IQ-TREE)

import argparse
import glob
import inspect
import logging
import os
import re
import shutil
import subprocess as sp
import sys
import tempfile
from collections import Counter
from typing import Dict, List, Optional, Set, Tuple

from pyfiglet import Figlet
from importlib.metadata import version


logger = logging.getLogger("ncortho")


# ---------------------------------------------------------------------------
# Custom argparse types
# ---------------------------------------------------------------------------
def _ratio_float(value: str) -> float:
    """Argparse type: float in [0.0, 1.0]."""
    try:
        fval = float(value)
    except ValueError:
        raise argparse.ArgumentTypeError(f"'{value}' is not a number")
    if not 0.0 <= fval <= 1.0:
        raise argparse.ArgumentTypeError(
            f"'{value}' must be between 0.0 and 1.0"
        )
    return fval


def _existing_dir(value: str) -> str:
    """Argparse type: directory must exist."""
    path = os.path.realpath(value)
    if not os.path.isdir(path):
        raise argparse.ArgumentTypeError(f"Directory not found: {path}")
    return path


def _existing_file(value: str) -> str:
    """Argparse type: file must exist."""
    path = os.path.realpath(value)
    if not os.path.isfile(path):
        raise argparse.ArgumentTypeError(f"File not found: {path}")
    return path


# ---------------------------------------------------------------------------
# MUSCLE version check
# ---------------------------------------------------------------------------
def _check_muscle_version() -> str:
    """
    Verify that MUSCLE >= 5.0 is available. Returns the version string.

    Raises SystemExit with an informative message if MUSCLE is missing
    or is version 3.x (commonly shipped with Biopython).
    """
    try:
        result = sp.run(
            ["muscle", "-version"],
            capture_output=True, text=True, check=False,
        )
        # MUSCLE 5 prints e.g. "muscle 5.1.linux64 ..."
        # MUSCLE 3.8 prints e.g. "MUSCLE v3.8.31 by Robert C. Edgar"
        output = (result.stdout + result.stderr).strip()
    except FileNotFoundError:
        raise SystemExit(
            "Error: 'muscle' not found on PATH.\n"
            "ncOrtho requires MUSCLE 5 (https://github.com/rcedgar/muscle).\n"
            "The MUSCLE 3.8 bundled with some Biopython installations is "
            "not compatible."
        )

    # Try to extract a major version number
    match = re.search(r"v?(\d+)\.\d+", output)
    if match:
        major = int(match.group(1))
        if major < 5:
            raise SystemExit(
                f"Error: MUSCLE {output.splitlines()[0]} detected, but "
                f"ncOrtho requires MUSCLE >= 5.0.\n"
                f"The '-align' flag used by ncOrtho was introduced in "
                f"MUSCLE 5. MUSCLE 3.8 (often shipped with Biopython) "
                f"uses different command-line syntax and is not compatible.\n"
                f"Please install MUSCLE 5: "
                f"https://github.com/rcedgar/muscle"
            )
        return output.splitlines()[0]

    # Could not parse — warn but continue
    logger.warning(
        "Could not determine MUSCLE version from: %s. "
        "Proceeding, but MUSCLE 5+ is required.",
        output,
    )
    return output


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _arglistcheck(arg: str) -> List[str]:
    """
    Parse a user-supplied argument that is either a comma-separated
    string or a path to a newline-separated file.  Returns a list of
    stripped, non-empty strings, or an empty list.
    """
    if not arg:
        return []
    if os.path.isfile(arg):
        with open(arg) as fh:
            return [el.strip() for el in fh.read().split("\n") if el.strip()]
    return [el.strip() for el in arg.split(",") if el.strip()]


def _parse_matseq(path: str) -> Dict[str, str]:
    """
    Read a miRNA info file and extract the 7-nt seed (positions 2–8)
    from the mature sequence column.

    Returns ``{mirid: seed_sequence}``.
    """
    m2s: Dict[str, str] = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 2:
                continue
            mature = fields[-1].replace("U", "T")
            mirid = fields[1].split("_")[0]
            if len(mature) < 8:
                raise ValueError(
                    f"Mature sequence too short for {mirid} in ncrna_file "
                    f"(need >= 8 nt, got {len(mature)})"
                )
            m2s[mirid] = mature[1:8]
    return m2s


def _make_supermatrix(out: str) -> str:
    """
    Concatenate per-miRNA alignments into a supermatrix, de-gap, and
    return the path to the processed alignment.
    """
    align_out = os.path.join(out, "alignments")
    tree_out = os.path.join(out, "supermatrix")
    curr_dir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
    os.makedirs(tree_out, exist_ok=True)

    print("# Creating supermatrix alignment", flush=True)
    concat_script = os.path.join(curr_dir, "concat_alignments_dmp.pl")
    cmd = ["perl", concat_script, "-in", align_out, "-out", "supermatrix.aln"]
    res = sp.run(cmd, capture_output=True, text=True, check=False)
    if res.returncode != 0:
        raise RuntimeError(f"Supermatrix concatenation failed:\n{res.stderr}")

    print("# De-gapping alignment", flush=True)
    degap_script = os.path.join(curr_dir, "degapper.pl")
    supermatrix_aln = os.path.join(out, "supermatrix.aln")
    cmd = ["perl", degap_script, "-in", supermatrix_aln]
    sp.run(cmd, capture_output=True, text=True, check=False)

    # Move outputs into the supermatrix directory
    for fname in ("supermatrix.aln", "subalignment_positions.txt", "supermatrix.aln.proc"):
        src = os.path.join(out, fname)
        if os.path.isfile(src):
            shutil.move(src, os.path.join(tree_out, fname))

    return os.path.join(tree_out, "supermatrix.aln.proc")


# ---------------------------------------------------------------------------
# Argument parser
# ---------------------------------------------------------------------------
def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="ncOrtho-analyze",
        description=(
            "Analyse ncOrtho results: PhyloProfile input, supermatrix "
            "alignment, and species tree."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser._action_groups.pop()
    required = parser.add_argument_group("Required arguments")
    optional = parser.add_argument_group("Optional arguments")

    required.add_argument(
        "-r", "--results",
        metavar="<path>",
        type=_existing_dir,
        required=True,
        help="Path to ncOrtho output directory to analyse.",
    )
    required.add_argument(
        "-o", "--output",
        metavar="<path>",
        type=str,
        required=True,
        help="Path to analysis output directory.",
    )
    required.add_argument(
        "-m", "--mapping",
        metavar="<path>",
        type=_existing_file,
        required=True,
        help=(
            "Tab-separated mapping: NCBI taxonomy ID <TAB> species/directory name "
            '(e.g. "9606\\tHomo_sapiens").'
        ),
    )

    optional.add_argument(
        "--skip",
        metavar="STR",
        type=str,
        default="",
        help="Comma-separated list or path to file of species to skip.",
    )
    optional.add_argument(
        "--include",
        metavar="STR",
        type=str,
        default="",
        help=(
            "Comma-separated list or path to file of species to include "
            "(all others are skipped)."
        ),
    )
    optional.add_argument(
        "--auto_skip",
        metavar="FLOAT",
        type=_ratio_float,
        default=0.0,
        help=(
            "Exclude species with fewer than this fraction of reference "
            "miRNAs detected (e.g. 0.5 = at least 50%%). Default: off."
        ),
    )
    optional.add_argument(
        "--iqtree",
        metavar="STR",
        type=str,
        default="",
        help=(
            "Custom iqtree command string. Default: "
            '"iqtree -s <supermatrix> -bb 1000 -alrt 1000 -nt AUTO -redo '
            '-pre <tree_out>/species_tree".'
        ),
    )
    optional.add_argument(
        "--mirnas",
        metavar="STR",
        type=str,
        default="",
        help=(
            "Comma-separated list or path to file of miRNAs to include "
            "in the reconstruction."
        ),
    )
    optional.add_argument(
        "--ncrna_file",
        metavar="<path>",
        type=str,
        default="",
        help=(
            "Path to reference miRNA info file (as used in ncCreate / "
            "ncSearch). Enables seed-conservation annotation in "
            "PhyloProfile output."
        ),
    )

    return parser


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> None:
    # Banner
    custom_fig = Figlet(font="stop")
    print(custom_fig.renderText("ncOrtho")[:-3], flush=True)
    v = version("ncOrtho")
    print(f"Version: {v}\n", flush=True)

    parser = build_parser()
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()

    # Check MUSCLE version early
    muscle_version = _check_muscle_version()
    print(f"# Using MUSCLE: {muscle_version}", flush=True)

    # Parse list arguments
    spec_to_skip: List[str] = _arglistcheck(args.skip)
    spec_include: List[str] = _arglistcheck(args.include)

    overlap = set(spec_to_skip) & set(spec_include)
    if overlap:
        raise SystemExit(
            f"Error: Species present in both --skip and --include: "
            f"{', '.join(sorted(overlap))}"
        )

    mirlist: List[str] = _arglistcheck(args.mirnas)
    mirid2seed: Dict[str, str] = {}
    if args.ncrna_file:
        if not os.path.isfile(args.ncrna_file):
            raise SystemExit(f"Error: --ncrna_file not found: {args.ncrna_file}")
        mirid2seed = _parse_matseq(args.ncrna_file)

    # Paths
    res_dir = args.results
    outdir = os.path.realpath(args.output)
    os.makedirs(outdir, exist_ok=True)

    # -----------------------------------------------------------------------
    # Read results
    # -----------------------------------------------------------------------
    ortholog_files = glob.glob(os.path.join(res_dir, "*", "*_orthologs.fa"))
    if not ortholog_files:
        ortholog_files = glob.glob(os.path.join(res_dir, "*.fa"))
    if not ortholog_files:
        raise SystemExit(
            f"Error: No ortholog result files found in {res_dir}"
        )

    # Read mapping
    name_2_taxid: Dict[str, str] = {}
    with open(args.mapping) as fh:
        for lineno, line in enumerate(fh, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                logger.warning(
                    "Skipping malformed line %d in mapping file.", lineno
                )
                continue
            taxid, name = parts[0], parts[1]
            name_2_taxid[name] = taxid

    # -----------------------------------------------------------------------
    # Build ortholog dict and write PhyloProfile input
    # -----------------------------------------------------------------------
    ortho_dict: Dict[str, Dict[str, str]] = {}
    spec_list: List[str] = []
    pp_path = os.path.join(outdir, "PhyloProfile.long")

    with open(pp_path, "w") as pp:
        header = "geneID\tncbiID\torthoID"
        if mirid2seed:
            header += "\tseedConservation"
        pp.write(header + "\n")

        for filepath in ortholog_files:
            current_spec = ""
            current_mirna = ""
            current_seq_parts: List[str] = []

            def _flush_record():
                """Write the current record to PhyloProfile and ortho_dict."""
                if not current_spec or not current_mirna:
                    return
                seq = "".join(current_seq_parts)
                if not seq:
                    return
                group_str = current_mirna.split("_")[0]
                taxid = name_2_taxid.get(current_spec)
                if taxid is None:
                    logger.warning(
                        "Species '%s' not in mapping file; skipping.",
                        current_spec,
                    )
                    return
                taxstr = f"ncbi{taxid}"
                line_parts = [group_str, taxstr, current_mirna]
                if mirid2seed:
                    seed = mirid2seed.get(group_str, "")
                    seedcon = 1 if seed and seed in seq else 0
                    line_parts.append(str(seedcon))
                pp.write("\t".join(line_parts) + "\n")

                if group_str not in ortho_dict:
                    ortho_dict[group_str] = {}
                if current_spec not in ortho_dict[group_str]:
                    ortho_dict[group_str][current_spec] = seq
                    spec_list.append(current_spec)

            with open(filepath) as fh:
                for line in fh:
                    line = line.strip()
                    if not line:
                        continue
                    if line.startswith(">"):
                        _flush_record()
                        parts = line.lstrip(">").split("|")
                        current_spec = parts[0] if len(parts) > 0 else ""
                        current_mirna = parts[1] if len(parts) > 1 else ""
                        current_seq_parts = []
                    else:
                        current_seq_parts.append(line)
                _flush_record()

    print("# Finished writing PhyloProfile input", flush=True)

    # Optional miRNA filter
    if mirlist:
        ortho_dict = {
            mid: seqs for mid, seqs in ortho_dict.items() if mid in mirlist
        }

    # -----------------------------------------------------------------------
    # Auto-skip species with few orthologs
    # -----------------------------------------------------------------------
    if args.auto_skip > 0:
        spec_counter = Counter(spec_list)
        max_count = spec_counter.most_common(1)[0][1]
        cutoff = max_count * args.auto_skip
        print(
            f"# Only species with >= {round(cutoff)} identified orthologs "
            "are used for tree calculation",
            flush=True,
        )
        for species, count in spec_counter.items():
            if count < cutoff and species not in spec_to_skip:
                spec_to_skip.append(species)
        if spec_to_skip:
            print("# Skipping:", flush=True)
            for sp_name in spec_to_skip:
                print(f"  {sp_name}", flush=True)

    if not spec_include:
        spec_include = list(set(spec_list))

    # -----------------------------------------------------------------------
    # Per-miRNA alignments (MUSCLE 5)
    # -----------------------------------------------------------------------
    print("# Starting alignments", flush=True)
    align_out = os.path.join(outdir, "alignments")
    os.makedirs(align_out, exist_ok=True)

    skip_set: Set[str] = set(spec_to_skip)
    include_set: Set[str] = set(spec_include)

    for mirna, spec_seqs in ortho_dict.items():
        aln_path = os.path.join(align_out, f"{mirna}.aln")

        # Filter species
        filtered = {
            sp_name: seq
            for sp_name, seq in spec_seqs.items()
            if sp_name not in skip_set or sp_name in include_set
        }

        if len(filtered) < 2:
            # Cannot align fewer than 2 sequences; write as-is
            with open(aln_path, "w") as fh:
                for sp_name, seq in filtered.items():
                    fh.write(f">{sp_name}\n{seq}\n")
            continue

        with tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False) as tmp:
            tmp_path = tmp.name
            for sp_name, seq in filtered.items():
                tmp.write(f">{sp_name}\n{seq}\n")

        try:
            muscle_cmd = ["muscle", "-align", tmp_path, "-output", aln_path]
            res = sp.run(muscle_cmd, capture_output=True, text=True, check=False)
            if res.returncode != 0:
                raise RuntimeError(
                    f"MUSCLE alignment failed for {mirna}:\n{res.stderr}"
                )
        finally:
            os.remove(tmp_path)

    sys.stdout.flush()

    # -----------------------------------------------------------------------
    # Supermatrix and species tree
    # -----------------------------------------------------------------------
    superm_path = _make_supermatrix(outdir)

    print("# Starting tree calculation", flush=True)
    tree_out = os.path.join(outdir, "species_tree")
    os.makedirs(tree_out, exist_ok=True)

    if args.iqtree:
        tree_cmd = args.iqtree
        # Still use shell=True for custom commands since the user controls the string
        res = sp.run(tree_cmd, shell=True, capture_output=True, text=True, check=False)
    else:
        tree_prefix = os.path.join(tree_out, "species_tree")
        tree_cmd_list = [
            "iqtree",
            "-s", superm_path,
            "-bb", "1000",
            "-alrt", "1000",
            "-nt", "AUTO",
            "-redo",
            "-pre", tree_prefix,
        ]
        res = sp.run(tree_cmd_list, capture_output=True, text=True, check=False)

    if res.stdout:
        print(res.stdout, flush=True)
    if res.returncode != 0:
        logger.error("IQ-TREE failed:\n%s", res.stderr)

    # -----------------------------------------------------------------------
    # Write PhyloProfile-compatible tree
    # -----------------------------------------------------------------------
    treepath = os.path.join(tree_out, "species_tree.contree")
    if os.path.isfile(treepath):
        with open(treepath) as fh:
            treestr = fh.read().strip()
        for name, taxid in name_2_taxid.items():
            treestr = treestr.replace(name, f"ncbi{taxid}")
        pptree = os.path.join(outdir, "PhyloProfile_tree.contree")
        with open(pptree, "w") as fh:
            fh.write(treestr)
        print("# PhyloProfile tree written", flush=True)
    else:
        logger.warning(
            "IQ-TREE consensus tree not found at %s; "
            "PhyloProfile tree not written.",
            treepath,
        )


if __name__ == "__main__":
    main()
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

# Utility functions for the coreset pipeline:
#   - Annotation parsing (GTF, GFF3, TSV)
#   - YAML parameter loading
#   - Pairwise-ortholog I/O
#   - miRNA data I/O

import logging
import os
import re
import sys
from typing import Any, Dict, List, Optional, Set, Tuple

import yaml


logger = logging.getLogger("ncortho")


# ---------------------------------------------------------------------------
# Annotation parsers
# ---------------------------------------------------------------------------
# All three parsers return a shared "chr_dict" structure:
#
#   chr_dict[gene_id]  = (contig, position_index)
#   chr_dict[contig]   = {position_index: (gene_id, start, end, strand)}
#
# where *position_index* is the 1-based ordinal of the gene on that contig
# (i.e. the first protein-coding gene on contig "chr1" has index 1, the
# second has index 2, etc.).


def gtf_parser(gtf: str) -> Dict[str, Any]:
    """Parse an Ensembl GTF file for protein-coding gene coordinates."""
    chr_dict: Dict[str, Any] = {}
    chromo = ""
    i = 0 

    with open(gtf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 3 or fields[2] != "gene":
                continue
            if "gene_biotype" not in line:
                continue
            biotype = line.split('gene_biotype')[1].split('"')[1]
            if biotype != "protein_coding":
                continue

            linedata = line.strip().split("\t")
            contig = linedata[0]
            geneid = linedata[-1].split('"')[1]
            start = int(linedata[3])
            end = int(linedata[4])
            strand = linedata[6]

            if contig != chromo:
                i = 1
                chromo = contig
            else:
                i += 1

            chr_dict[geneid] = (contig, i)
            if contig not in chr_dict:
                chr_dict[contig] = {}
            chr_dict[contig][i] = (geneid, start, end, strand)

    return chr_dict


def gff_parser(gff: str, id_type: str) -> Dict[str, Any]:
    """
    Parse a RefSeq GFF3 file for protein-coding gene coordinates.

    Parameters
    ----------
    gff : str
        Path to the GFF3 file.
    id_type : str
        Which attribute field to use as the gene identifier.
        One of ``"GeneID"``, ``"ID"``, ``"Name"``, ``"CDS"``.
    """
    chr_dict: Dict[str, Any] = {}
    chromo = ""
    i = 0

    with open(gff) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            if "gene_biotype=protein_coding" not in line:
                continue
            linedata = line.strip().split("\t")
            if linedata[2] != "gene":
                continue

            attributes = linedata[-1]
            if id_type == "GeneID":
                geneid = re.split(r"[;,]", attributes.split(f"{id_type}:")[1])[0]
            elif id_type in ("ID", "Name", "CDS"):
                geneid = attributes.split(f"{id_type}=")[1].split(";")[0]
                if id_type in ("ID", "CDS"):
                    geneid = "-".join(geneid.split("-")[1:])
            else:
                raise ValueError(f'Unknown identifier type "{id_type}"')

            contig = linedata[0]
            start = int(linedata[3])
            end = int(linedata[4])
            strand = linedata[6]

            if contig != chromo:
                i = 1
                chromo = contig
            else:
                i += 1

            chr_dict[geneid] = (contig, i)
            if contig not in chr_dict:
                chr_dict[contig] = {}
            chr_dict[contig][i] = (geneid, start, end, strand)

    return chr_dict


def table_parser(path: str) -> Dict[str, Any]:
    """Parse a simple TSV annotation (contig, start, end, strand, gene_id)."""
    chr_dict: Dict[str, Any] = {}
    chromo = ""
    i = 0

    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            contig, start, end, strand, geneid = line.strip().split("\t")

            if contig != chromo:
                i = 1
                chromo = contig
            else:
                i += 1

            chr_dict[geneid] = (contig, i)
            if contig not in chr_dict:
                chr_dict[contig] = {}
            chr_dict[contig][i] = (geneid, int(start), int(end), strand)

    return chr_dict


def parse_annotation(path: str, idtype: str) -> Dict[str, Any]:
    """
    Dispatch to the correct annotation parser based on file extension.

    Supported extensions: ``.gtf``, ``.gff``, ``.gff3``, ``.tsv``, ``.txt``.
    """
    ext = os.path.splitext(path)[1].lstrip(".").lower()
    if ext == "gtf":
        return gtf_parser(path)
    if ext in ("gff3", "gff"):
        return gff_parser(path, idtype)
    if ext in ("tsv", "txt"):
        return table_parser(path)

    raise ValueError(
        f'Unsupported annotation extension ".{ext}". '
        "Expected .gtf, .gff, .gff3, .tsv, or .txt."
    )


# ---------------------------------------------------------------------------
# YAML parameter parsing
# ---------------------------------------------------------------------------
def parse_yaml(path: str) -> Tuple[Dict, Dict, List[str]]:
    """
    Parse the ncOrtho YAML parameter file.

    Returns
    -------
    core_dict : dict
        ``{species_name: {"orthologs": path, "genome": path, "annotation": path}}``.
    ref_dict : dict
        ``{"annotation": path, "genome": path}`` for the reference.
    all_paths : list[str]
        Every file path mentioned (for existence checks).
    """
    core_dict: Dict[str, Any] = {}
    all_paths: List[str] = []
    ref_dict: Dict[str, str] = {}

    with open(path) as fh:
        for entry in yaml.safe_load_all(fh):
            in_type = entry.pop("type")
            if in_type == "reference":
                ref_dict["annotation"] = entry["annotation"]
                ref_dict["genome"] = entry["genome"]
            elif in_type == "core":
                species = entry.pop("name")
                core_dict[species] = entry
            else:
                logger.warning("Unknown entry type '%s' in YAML; skipping.", in_type)
                continue
            all_paths.extend([entry.get("annotation", ""), entry.get("genome", "")])

    # Filter out empty strings from paths
    all_paths = [p for p in all_paths if p]

    return core_dict, ref_dict, all_paths


# ---------------------------------------------------------------------------
# Pairwise ortholog I/O
# ---------------------------------------------------------------------------
def read_pairwise_orthologs(pathdict: Dict[str, Dict]) -> Dict[str, Dict[str, List[str]]]:
    """
    Read pairwise-ortholog files for each core species.

    Returns
    -------
    dict
        ``{taxon: {ref_gene: [ortho1, ortho2, ...], ...}, ...}``.
    """
    od: Dict[str, Dict[str, List[str]]] = {}
    for taxon, pd in pathdict.items():
        taxd: Dict[str, List[str]] = {}
        orthofile = pd["orthologs"]
        with open(orthofile) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                parts = line.strip().split()
                if len(parts) < 2:
                    continue
                ref, core = parts[:2]
                taxd.setdefault(ref, []).append(core)
        od[taxon] = taxd
    return od


def gene_from_cds(gff: str) -> Dict[str, str]:
    """
    Map CDS accessions to gene identifiers from a GFF3 file.

    Returns ``{cds_accession: gene_id}``.
    """
    cds2gene: Dict[str, str] = {}
    cds_mrna_map: Dict[str, str] = {}

    with open(gff) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            data = line.strip().split("\t")
            if len(data) < 9:
                continue
            if data[2] != "CDS":
                continue

            cdsid = (
                data[-1].split(";")[0]
                .split(".")[0]
                .replace("ID=cds-", "")
                .replace("ID=id-", "")
            )
            parent = data[-1].split("Parent=")[1].split(";")[0].split(".")[0]

            if parent.startswith("gene-") or parent.startswith("id-"):
                cds2gene[cdsid] = parent.replace("gene-", "").replace("id-", "")
            else:
                cds_mrna_map[cdsid] = parent
                cds_mrna_map[parent] = cdsid

    # Second pass: resolve CDS→mRNA→gene links
    with open(gff) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            data = line.strip().split("\t")
            if len(data) < 9:
                continue
            if data[2] != "mRNA":
                continue

            mrnaid = data[-1].split("ID=")[1].split(".")[0]
            parent = data[-1].split("Parent=")[1].split(".")[0]

            if mrnaid not in cds_mrna_map:
                continue
            cdsid = cds_mrna_map[mrnaid]
            if parent.startswith("gene-"):
                cds2gene[cdsid] = parent.split(";")[0].replace("gene-", "")

    return cds2gene


def pairwise_orthologs_from_cds(
    cd: Dict[str, Dict], refanno: str
) -> Dict[str, Dict[str, List[str]]]:
    """
    Build pairwise-ortholog mappings when ortholog files use CDS
    accessions rather than gene IDs.

    Returns ``{core_species: {ref_gene: [ortho_gene1, ...], ...}, ...}``.
    """
    ref_cds2gene = gene_from_cds(refanno)
    od: Dict[str, Dict[str, List[str]]] = {}

    for cspec, spec_data in cd.items():
        cds2gene = gene_from_cds(spec_data["annotation"])
        od[cspec] = {}
        orthopath = spec_data["orthologs"]

        with open(orthopath) as fh:
            for line in fh:
                parts = line.strip().split()
                if len(parts) < 2:
                    continue
                cds_refprot, cds_orthoprot = parts[:2]

                refprot = ref_cds2gene.get(cds_refprot)
                if refprot is None:
                    # CDS may have been removed in a later annotation version
                    continue

                orthoprot = cds2gene.get(cds_orthoprot)
                if orthoprot is None:
                    logger.warning(
                        "CDS '%s' not found in %s annotation; skipping.",
                        cds_orthoprot, cspec,
                    )
                    continue

                od[cspec].setdefault(refprot, []).append(orthoprot)

    return od


# ---------------------------------------------------------------------------
# miRNA I/O
# ---------------------------------------------------------------------------
def read_mirnas(path: str) -> List[List[str]]:
    """
    Read a tab-separated miRNA information file.

    Lines starting with ``#`` are ignored. Returns a list of rows, each
    row being a list of tab-separated fields.
    """
    mirnas: List[List[str]] = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            mirnas.append(line.split("\t"))

    if not mirnas:
        logger.warning("No miRNAs loaded from %s", path)

    return mirnas


def mirna_position(mirlist: List[str]) -> Tuple[str, str, int, int, str]:
    """
    Extract miRNA positional data from a row of the miRNA info file.

    Returns ``(mirid, chromosome, start, end, strand)`` with
    coordinates as integers and the ``chr`` prefix stripped.
    """
    mirid, rawchrom, start, end, strand = mirlist[:5]
    chromo = rawchrom.replace("chr", "")
    return mirid, chromo, int(start), int(end), strand
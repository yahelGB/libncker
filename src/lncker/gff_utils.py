from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Dict, Iterator, Set, Tuple

from lncker.io_utils import open_text

if TYPE_CHECKING:
    from lncker.validate import DropLog


@dataclass(frozen=True)
class GffFeature:
    contig: str
    ftype: str
    start: int
    end: int
    strand: str
    attrs: Dict[str, str]


def parse_attrs(s: str) -> Dict[str, str]:
    out: Dict[str, str] = {}
    for item in s.split(";"):
        if "=" not in item:
            continue
        k, v = item.split("=", 1)
        out[k] = v
    return out


def norm_rna_id(raw: str) -> str:
    # "rna-XR_..." -> "XR_..."
    return raw[4:] if raw.startswith("rna-") else raw


def iter_gff(gff: Path, drops: "DropLog | None" = None) -> Iterator[GffFeature]:
    """Yield every well-formed feature.

    A record that cannot be parsed is reported to the drop log and skipped --
    never silently ignored, and never allowed to raise. Pass ``drops`` to
    collect the counts; otherwise they go to the process-wide log.
    """
    from lncker.validate import default_drop_log, parse_gff_record

    log = drops if drops is not None else default_drop_log()
    with open_text(gff) as f:
        for lineno, line in enumerate(f, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            rec = parse_gff_record(line, lineno, log)
            if rec is None:
                continue
            contig, _src, ftype, start, end, strand, _phase, attrs_s = rec
            yield GffFeature(
                contig=contig,
                ftype=ftype,
                start=start,
                end=end,
                strand=strand,
                attrs=parse_attrs(attrs_s),
            )


def extract_lncrna_ids_from_gff(gff: Path, level: str = "transcript") -> Tuple[Set[str], Dict[str, int]]:
    """
    NCBI/RefSeq GFF3 typical:
      - gene: ID=gene-...; gene_biotype=lncRNA   (may be missing)
      - lnc_RNA: ID=rna-XR_/NR_; Parent=gene-...

    level:
      - transcript: returns XR_/NR_ set
      - gene: returns gene-... set

    stats:
      - genes_lncRNA
      - transcripts_lncRNA
      - transcripts_prefix_XR / NR / others...
    """
    if level not in {"transcript", "gene"}:
        raise ValueError("level must be 'transcript' or 'gene'")

    genes_by_biotype: Set[str] = set()
    genes_by_parent: Set[str] = set()
    transcripts: Set[str] = set()
    prefix_counts: Dict[str, int] = {}

    # pass 1: collect gene IDs that explicitly say gene_biotype=lncRNA
    for feat in iter_gff(gff):
        if feat.ftype != "gene":
            continue
        if feat.attrs.get("gene_biotype", "") != "lncRNA":
            continue
        gid = feat.attrs.get("ID", "")
        if gid:
            genes_by_biotype.add(gid)

    # pass 2: collect lnc_RNA transcripts and their parents (robust even if gene_biotype missing)
    for feat in iter_gff(gff):
        if feat.ftype != "lnc_RNA":
            continue
        parent = feat.attrs.get("Parent", "")
        if parent:
            genes_by_parent.add(parent)

        tid = norm_rna_id(feat.attrs.get("ID", ""))
        if tid:
            transcripts.add(tid)
            pref = tid.split("_", 1)[0] if "_" in tid else tid[:2]
            prefix_counts[pref] = prefix_counts.get(pref, 0) + 1

    genes_lnc = genes_by_biotype | genes_by_parent

    stats: Dict[str, int] = {
        "genes_lncRNA": len(genes_lnc),
        "genes_lncRNA_by_biotype": len(genes_by_biotype),
        "genes_lncRNA_by_parent": len(genes_by_parent),
        "transcripts_lncRNA": len(transcripts),
    }
    for k, v in sorted(prefix_counts.items()):
        stats[f"transcripts_prefix_{k}"] = v

    if level == "gene":
        return genes_lnc, stats
    return transcripts, stats


def build_transcript_to_gene_map(gff: Path) -> Dict[str, str]:
    """
    Map transcript accession IDs to Parent gene IDs.

    Supports NCBI style where mRNA/lnc_RNA have:
      ID=rna-XM_/NM_/XR_/NR_ ; Parent=gene-...

    Also maps if ID already without rna- prefix.
    """
    m: Dict[str, str] = {}
    for feat in iter_gff(gff):
        if feat.ftype not in {"mRNA", "lnc_RNA"}:
            continue
        tid = norm_rna_id(feat.attrs.get("ID", ""))
        parent = feat.attrs.get("Parent", "")
        if tid and parent:
            m[tid] = parent
    return m


def build_lnc_gene_to_transcripts(gff: Path) -> Dict[str, Set[str]]:
    """
    gene-... -> set(XR_/NR_/...) inferred from lnc_RNA features directly.
    Robust even when gene features do not have gene_biotype=lncRNA.
    """
    out: Dict[str, Set[str]] = {}
    for feat in iter_gff(gff):
        if feat.ftype != "lnc_RNA":
            continue
        parent = feat.attrs.get("Parent", "")
        if not parent:
            continue
        tid = norm_rna_id(feat.attrs.get("ID", ""))
        if tid:
            out.setdefault(parent, set()).add(tid)
    return out


def write_extract_ids_summary(
    path: Path,
    gff: Path,
    level: str,
    stats: Dict[str, int],
    prov: "object | None" = None,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        if prov is not None:
            prov.write_header(f)
        f.write("key\tvalue\n")
        f.write(f"gff\t{gff}\n")
        f.write(f"level\t{level}\n")
        for k in sorted(stats.keys()):
            f.write(f"{k}\t{stats[k]}\n")

from __future__ import annotations

import csv
import json
import time
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

_NCBI_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
_DEFAULT_DELAY = 0.4  # seconds between NCBI requests (safe without an API key)
_FIELDNAMES = [
    "mRNA_ID",
    "product",
    "lncRNA_IDs",
    "mRNA_positions",
    "ncbi_summary",
    "go_biological_process",
    "go_molecular_function",
    "go_cellular_component",
]

_GO_HEADINGS: dict[str, str] = {
    "Biological Process": "biological_process",
    "Molecular Function": "molecular_function",
    "Cellular Component": "cellular_component",
}


@dataclass(frozen=True)
class _MRNARecord:
    mRNA_ID: str
    product: str
    lncRNA_IDs: list[str]
    mRNA_positions: list[str]


# ---------------------------------------------------------------------------
# NCBI helpers (stdlib urllib + xml.etree only)
# ---------------------------------------------------------------------------

def _http_get(url: str) -> bytes:
    with urllib.request.urlopen(url, timeout=30) as resp:  # noqa: S310
        return resp.read()


def _fetch_gene_info(
    accession: str,
    api_key: str,
    delay: float,
) -> tuple[str, dict[str, list[str]]]:
    """XM_ accession → (ncbi_summary, go_terms).

    go_terms keys: biological_process, molecular_function, cellular_component.
    Each value is a list of GO term names (e.g. ["metabolic process", "proteolysis"]).
    Returns ("", empty dict) on any error.
    """
    empty_go: dict[str, list[str]] = {
        "biological_process": [],
        "molecular_function": [],
        "cellular_component": [],
    }

    def _p(d: dict[str, str]) -> str:
        if api_key:
            d["api_key"] = api_key
        return urllib.parse.urlencode(d)

    try:
        # 1) esearch accession directly in gene DB — more reliable than
        #    nucleotide→gene elink for predicted RefSeq records (XM_/XR_)
        acc_bare = accession.split(".")[0]  # strip version suffix
        data = json.loads(_http_get(
            f"{_NCBI_BASE}/esearch.fcgi?"
            + _p({"db": "gene", "term": f"{acc_bare}[accn]", "retmode": "json", "retmax": "1"})
        ))
        gene_ids = data.get("esearchresult", {}).get("idlist", [])
        if not gene_ids:
            return "", empty_go
        time.sleep(delay)

        # 2) efetch XML: gene UID → full record (summary + GO terms)
        xml_bytes = _http_get(
            f"{_NCBI_BASE}/efetch.fcgi?"
            + _p({"db": "gene", "id": str(gene_ids[0]), "retmode": "xml"})
        )
        return _parse_gene_xml(xml_bytes)

    except (urllib.error.URLError, KeyError, IndexError, json.JSONDecodeError, OSError):
        return "", empty_go


def _parse_gene_xml(xml_bytes: bytes) -> tuple[str, dict[str, list[str]]]:
    """Parse NCBI Gene XML and return (summary, go_terms)."""
    go_terms: dict[str, list[str]] = {
        "biological_process": [],
        "molecular_function": [],
        "cellular_component": [],
    }
    summary = ""

    try:
        root = ET.fromstring(xml_bytes)
    except ET.ParseError:
        return summary, go_terms

    # Gene summary text
    el = root.find(".//Entrezgene_summary")
    if el is not None and el.text:
        summary = el.text.strip()

    # GO terms: find Gene-commentary elements whose heading matches a GO category
    # Structure: Gene-commentary_heading → Gene-commentary_comment →
    #            Gene-commentary → Gene-commentary_source →
    #            Other-source → Other-source_anchor (term name)
    for commentary in root.iter("Gene-commentary"):
        heading_el = commentary.find("Gene-commentary_heading")
        if heading_el is None or heading_el.text not in _GO_HEADINGS:
            continue
        go_key = _GO_HEADINGS[heading_el.text]
        for anchor in commentary.iter("Other-source_anchor"):
            if anchor.text and anchor.text.strip():
                term = anchor.text.strip()
                if term not in go_terms[go_key]:
                    go_terms[go_key].append(term)

    return summary, go_terms


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def _load_de_records_for_tissue(path: Path, positions: set[str]) -> list[_MRNARecord]:
    """Read one cis output file, keep DE rows at the given positions, return unique mRNAs."""
    lnc_map: dict[str, set[str]] = defaultdict(set)
    pos_map: dict[str, set[str]] = defaultdict(set)
    product_map: dict[str, str] = {}

    with path.open(encoding="utf-8") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {col: i for i, col in enumerate(header)}
        required = {"mRNA_ID", "lncRNA_ID", "mRNA_position", "DE_status", "product"}
        missing = required - idx.keys()
        if missing:
            raise ValueError(f"{path.name}: missing columns {missing}")
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            if cols[idx["DE_status"]] != "DE":
                continue
            if cols[idx["mRNA_position"]] not in positions:
                continue
            mRNA_ID = cols[idx["mRNA_ID"]]
            lnc_map[mRNA_ID].add(cols[idx["lncRNA_ID"]])
            pos_map[mRNA_ID].add(cols[idx["mRNA_position"]])
            if mRNA_ID not in product_map:
                product_map[mRNA_ID] = cols[idx["product"]]

    return [
        _MRNARecord(
            mRNA_ID=mid,
            product=product_map[mid],
            lncRNA_IDs=sorted(lnc_map[mid]),
            mRNA_positions=sorted(pos_map[mid]),
        )
        for mid in sorted(product_map)
    ]


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

_DEFAULT_POSITIONS = {"1_U", "1_D"}


def run_annotate(
    cis_dir: Path,
    outdir: Path,
    ncbi_api_key: str = "",
    delay: float = _DEFAULT_DELAY,
    positions: set[str] = _DEFAULT_POSITIONS,
    use_emapper: bool = False,
    gff: Path | None = None,
    emapper_db: Path | None = None,
    tax_scope: int | None = None,
    emapper_cpu: int = 4,
) -> None:
    cis_files = sorted(cis_dir.glob("*_cis_regulation_module_output.txt"))
    if not cis_files:
        raise ValueError(f"No *_cis_regulation_module_output.txt files found in {cis_dir}")

    outdir.mkdir(parents=True, exist_ok=True)
    print(f"Filtering to positions: {', '.join(sorted(positions))}")

    for cis_path in cis_files:
        tissue = cis_path.name.replace("_cis_regulation_module_output.txt", "")
        records = _load_de_records_for_tissue(cis_path, positions=positions)
        out = outdir / f"{tissue}_de_mrna_annotations.tsv"

        print(f"\n[{tissue}] {len(records)} unique DE mRNAs → {out.name}")

        # Optional: run eggNOG-mapper once per tissue for all mRNAs (batch)
        emapper_go: dict[str, dict[str, list[str]]] = {}
        if use_emapper:
            from libncker.emapper import run_emapper_annotation
            print(f"  Running eggNOG-mapper for {len(records)} mRNAs…")
            emapper_go = run_emapper_annotation(
                mrna_ids=[r.mRNA_ID for r in records],
                gff_path=gff,
                api_key=ncbi_api_key,
                delay=delay,
                emapper_db=emapper_db,
                tax_scope_override=tax_scope,
                cpu=emapper_cpu,
                outdir=outdir,
                tissue=tissue,
            )

        with out.open("w", encoding="utf-8", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=_FIELDNAMES, delimiter="\t")
            writer.writeheader()

            for i, rec in enumerate(records, 1):
                print(f"  [{i}/{len(records)}] {rec.mRNA_ID}  {rec.product[:55]}")

                summary, ncbi_go = _fetch_gene_info(rec.mRNA_ID, api_key=ncbi_api_key, delay=delay)
                time.sleep(delay)

                # Prefer eggNOG-mapper GO terms when available; fall back to NCBI
                go = emapper_go.get(rec.mRNA_ID, ncbi_go)

                writer.writerow(
                    {
                        "mRNA_ID": rec.mRNA_ID,
                        "product": rec.product,
                        "lncRNA_IDs": ";".join(rec.lncRNA_IDs),
                        "mRNA_positions": ";".join(rec.mRNA_positions),
                        "ncbi_summary": summary,
                        "go_biological_process": ";".join(go["biological_process"]),
                        "go_molecular_function": ";".join(go["molecular_function"]),
                        "go_cellular_component": ";".join(go["cellular_component"]),
                    }
                )
                fh.flush()

    print("\nDone.")

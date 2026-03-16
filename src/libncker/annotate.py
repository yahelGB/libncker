from __future__ import annotations

import csv
import json
import time
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass
from pathlib import Path

_NCBI_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
_DEFAULT_DELAY = 0.4  # seconds between NCBI requests (safe without an API key)
_FIELDNAMES = ["mRNA_ID", "product", "lncRNA_IDs", "mRNA_positions", "ncbi_summary"]


@dataclass(frozen=True)
class _MRNARecord:
    mRNA_ID: str
    product: str
    lncRNA_IDs: list[str]
    mRNA_positions: list[str]


def _fetch_gene_summary(accession: str, api_key: str, delay: float) -> str:
    """Convert an XM_ accession → nuccore UID → Gene UID → gene summary text."""
    params: dict[str, str] = {
        "db": "nucleotide",
        "term": accession,
        "retmode": "json",
        "retmax": "1",
    }
    if api_key:
        params["api_key"] = api_key

    try:
        url = f"{_NCBI_BASE}/esearch.fcgi?" + urllib.parse.urlencode(params)
        with urllib.request.urlopen(url, timeout=15) as resp:  # noqa: S310
            data = json.loads(resp.read().decode("utf-8"))
        nuccore_ids = data.get("esearchresult", {}).get("idlist", [])
        if not nuccore_ids:
            return ""
        time.sleep(delay)

        link_params: dict[str, str] = {
            "dbfrom": "nucleotide",
            "db": "gene",
            "id": nuccore_ids[0],
            "retmode": "json",
        }
        if api_key:
            link_params["api_key"] = api_key
        url = f"{_NCBI_BASE}/elink.fcgi?" + urllib.parse.urlencode(link_params)
        with urllib.request.urlopen(url, timeout=15) as resp:  # noqa: S310
            link_data = json.loads(resp.read().decode("utf-8"))

        try:
            gene_ids: list[str] = link_data["linksets"][0]["linksetdbs"][0]["links"]
        except (KeyError, IndexError):
            return ""
        if not gene_ids:
            return ""
        time.sleep(delay)

        sum_params: dict[str, str] = {"db": "gene", "id": str(gene_ids[0]), "retmode": "json"}
        if api_key:
            sum_params["api_key"] = api_key
        url = f"{_NCBI_BASE}/esummary.fcgi?" + urllib.parse.urlencode(sum_params)
        with urllib.request.urlopen(url, timeout=15) as resp:  # noqa: S310
            sum_data = json.loads(resp.read().decode("utf-8"))

        gene_result = sum_data.get("result", {}).get(str(gene_ids[0]), {})
        return gene_result.get("summary", "")

    except (urllib.error.URLError, KeyError, IndexError, json.JSONDecodeError, OSError):
        return ""


def _load_de_records_for_tissue(path: Path) -> list[_MRNARecord]:
    """Read one cis output file, keep DE rows, return unique mRNAs for that tissue."""
    from collections import defaultdict

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


def run_annotate(
    cis_dir: Path,
    outdir: Path,
    ncbi_api_key: str = "",
    delay: float = _DEFAULT_DELAY,
) -> None:
    cis_files = sorted(cis_dir.glob("*_cis_regulation_module_output.txt"))
    if not cis_files:
        raise ValueError(f"No *_cis_regulation_module_output.txt files found in {cis_dir}")

    outdir.mkdir(parents=True, exist_ok=True)

    for cis_path in cis_files:
        tissue = cis_path.name.replace("_cis_regulation_module_output.txt", "")
        records = _load_de_records_for_tissue(cis_path)
        out = outdir / f"{tissue}_de_mrna_annotations.tsv"

        print(f"\n[{tissue}] {len(records)} unique DE mRNAs → {out.name}")

        with out.open("w", encoding="utf-8", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=_FIELDNAMES, delimiter="\t")
            writer.writeheader()

            for i, rec in enumerate(records, 1):
                print(f"  [{i}/{len(records)}] {rec.mRNA_ID}  {rec.product[:55]}")

                ncbi_summary = _fetch_gene_summary(rec.mRNA_ID, api_key=ncbi_api_key, delay=delay)
                time.sleep(delay)

                writer.writerow(
                    {
                        "mRNA_ID": rec.mRNA_ID,
                        "product": rec.product,
                        "lncRNA_IDs": ";".join(rec.lncRNA_IDs),
                        "mRNA_positions": ";".join(rec.mRNA_positions),
                        "ncbi_summary": ncbi_summary,
                    }
                )
                fh.flush()

    print("\nDone.")

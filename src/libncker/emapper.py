"""eggNOG-mapper web API integration (stdlib only).

Workflow:
  1. Detect the organism's taxonomic scope from the GFF assembly accession
  2. Fetch protein FASTAs for XM_ mRNA IDs via NCBI (elink nucleotide→protein)
  3. Submit FASTAs to eggNOG-mapper web API with the detected tax_scope
  4. Poll until done, fetch annotations
  5. Look up GO term names and aspects via QuickGO API
  6. Return per-mRNA GO annotation dicts for all three GO categories
"""
from __future__ import annotations

import json
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

_NCBI_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
_EMAPPER_BASE = "http://eggnog-mapper.embl.de"
_QUICKGO_BASE = "https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms"

# eggNOG taxonomic scope IDs that the web server recognises.
# Keys are NCBI taxon IDs of lineage nodes; values are the same ID
# (used as tax_scope in the API call).
_EGGNOC_SCOPE_IDS: list[int] = [
    6656,    # Arthropoda  ← closest for Crustacea / L. vannamei
    33208,   # Metazoa
    2759,    # Eukaryota
    1,       # root
]

# Hard-coded taxon→scope shortcuts for common non-model organisms.
# Saves extra NCBI taxonomy walk calls.
_TAXON_SCOPE: dict[int, int] = {
    6689: 6656,   # Litopenaeus vannamei → Arthropoda
    6683: 6656,   # Penaeus monodon      → Arthropoda
}


# ---------------------------------------------------------------------------
# Low-level HTTP helpers
# ---------------------------------------------------------------------------

def _http_get(url: str, timeout: int = 60) -> bytes:
    with urllib.request.urlopen(url, timeout=timeout) as resp:  # noqa: S310
        return resp.read()


def _http_post(url: str, payload: bytes, content_type: str = "application/json") -> bytes:
    req = urllib.request.Request(
        url,
        data=payload,
        headers={"Content-Type": content_type, "Accept": "application/json"},
    )
    with urllib.request.urlopen(req, timeout=120) as resp:  # noqa: S310
        return resp.read()


def _qs(params: dict, api_key: str = "") -> str:
    if api_key:
        params = {**params, "api_key": api_key}
    return urllib.parse.urlencode(params)


# ---------------------------------------------------------------------------
# Step 1 – Detect taxonomic scope from GFF
# ---------------------------------------------------------------------------

def detect_tax_scope(gff_path: Path, api_key: str = "", delay: float = 0.4) -> int:
    """Return an eggNOG-compatible taxon ID for the organism in *gff_path*.

    Reads the GFF header comments for a ``genome-build-accession`` line,
    queries NCBI Assembly to get the species taxon ID, and then walks up
    the NCBI lineage until it finds a node in *_EGGNOC_SCOPE_IDS*.

    Falls back to Arthropoda (6656) or Metazoa (33208) on any error.
    """
    accession = _parse_assembly_accession(gff_path)
    if not accession:
        print("  [emapper] Could not detect assembly accession; defaulting tax_scope=Metazoa")
        return 33208

    try:
        taxon_id = _accession_to_taxon(accession, api_key, delay)
    except Exception as exc:
        print(f"  [emapper] NCBI taxonomy lookup failed ({exc}); defaulting tax_scope=Metazoa")
        return 33208

    scope = _taxon_to_eggnoc_scope(taxon_id, api_key, delay)
    print(f"  [emapper] Assembly {accession} → taxon {taxon_id} → eggNOG scope {scope}")
    return scope


def _parse_assembly_accession(gff_path: Path) -> str:
    """Extract GCF_*/GCA_* accession from GFF header comments or filename."""
    try:
        with gff_path.open(encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if not line.startswith("#"):
                    break
                # e.g. "#!genome-build-accession NCBI_Assembly:GCF_003789085.1"
                if "genome-build-accession" in line:
                    for token in line.split():
                        if ":" in token:
                            token = token.split(":")[-1]
                        if token.startswith(("GCF_", "GCA_")):
                            return token
    except OSError:
        pass
    # Fallback: parse from filename
    for token in gff_path.name.replace("-", "_").split("_"):
        if token.startswith(("GCF", "GCA")) and len(token) > 3:
            # Reconstruct e.g. "GCF_003789085"
            pass
    # Try a broader scan of the filename parts
    name = gff_path.name
    for i, char in enumerate(name):
        if name[i:i+4] in ("GCF_", "GCA_"):
            end = name.find("_", i + 4)
            end2 = name.find(".", i + 4)
            end = min(e for e in (end, end2, len(name)) if e > i)
            return name[i:end]
    return ""


def _accession_to_taxon(accession: str, api_key: str, delay: float) -> int:
    """Query NCBI Assembly esummary to get the species taxon ID."""
    data = json.loads(_http_get(
        f"{_NCBI_BASE}/esearch.fcgi?"
        + _qs({"db": "assembly", "term": accession, "retmode": "json", "retmax": "1"}, api_key)
    ))
    ids = data.get("esearchresult", {}).get("idlist", [])
    if not ids:
        return 33208
    time.sleep(delay)

    summary = json.loads(_http_get(
        f"{_NCBI_BASE}/esummary.fcgi?"
        + _qs({"db": "assembly", "id": ids[0], "retmode": "json"}, api_key)
    ))
    uid_data = summary.get("result", {}).get(ids[0], {})
    return int(uid_data.get("taxid", 33208))


def _taxon_to_eggnoc_scope(taxon_id: int, api_key: str, delay: float) -> int:
    """Map a NCBI taxon ID to the closest eggNOG-compatible scope ID."""
    if taxon_id in _TAXON_SCOPE:
        return _TAXON_SCOPE[taxon_id]
    if taxon_id in _EGGNOC_SCOPE_IDS:
        return taxon_id

    # Walk up the NCBI taxonomy lineage
    try:
        import xml.etree.ElementTree as ET
        xml_bytes = _http_get(
            f"{_NCBI_BASE}/efetch.fcgi?"
            + _qs({"db": "taxonomy", "id": str(taxon_id), "retmode": "xml"}, api_key)
        )
        root = ET.fromstring(xml_bytes)
        scope_set = set(_EGGNOC_SCOPE_IDS)
        for taxon_el in root.iter("Taxon"):
            tid_el = taxon_el.find("TaxId")
            if tid_el is not None and tid_el.text:
                tid = int(tid_el.text)
                if tid in scope_set:
                    return tid
    except Exception:
        pass

    return 33208  # Metazoa fallback


# ---------------------------------------------------------------------------
# Step 2 – Fetch protein FASTAs for XM_ mRNAs via NCBI elink
# ---------------------------------------------------------------------------

def fetch_protein_fastas(
    mrna_ids: list[str],
    api_key: str = "",
    delay: float = 0.4,
) -> dict[str, str]:
    """Return ``{mrna_id: protein_fasta}`` for each XM_/XR_ accession.

    Uses NCBI elink (nucleotide → protein) to find the corresponding
    protein, then efetch to retrieve the FASTA sequence.
    Skips accessions for which no protein link is found.
    The FASTA header of each entry is replaced with the original *mrna_id*
    so eggNOG-mapper returns results keyed by mRNA ID.
    """
    print(f"  [emapper] Fetching protein FASTAs for {len(mrna_ids)} mRNAs…")

    # 1) esearch every accession in nucleotide DB to get numeric UIDs
    nucl_uid_map: dict[str, str] = {}  # mrna_id → nucleotide UID
    for mrna_id in mrna_ids:
        bare = mrna_id.split(".")[0]
        try:
            data = json.loads(_http_get(
                f"{_NCBI_BASE}/esearch.fcgi?"
                + _qs({"db": "nucleotide", "term": bare, "retmode": "json", "retmax": "1"}, api_key)
            ))
            ids = data.get("esearchresult", {}).get("idlist", [])
            if ids:
                nucl_uid_map[mrna_id] = ids[0]
        except Exception:
            pass
        time.sleep(delay)

    if not nucl_uid_map:
        return {}

    # 2) elink nucleotide → protein (batch, up to 50 at a time)
    prot_uid_map: dict[str, str] = {}  # mrna_id → protein UID
    nucl_items = list(nucl_uid_map.items())
    for i in range(0, len(nucl_items), 50):
        batch = nucl_items[i : i + 50]
        ids_str = ",".join(uid for _, uid in batch)
        try:
            link_data = json.loads(_http_get(
                f"{_NCBI_BASE}/elink.fcgi?"
                + _qs(
                    {"dbfrom": "nucleotide", "db": "protein", "id": ids_str,
                     "retmode": "json", "cmd": "neighbor"},
                    api_key,
                )
            ))
            for ls in link_data.get("linksets", []):
                input_uid = str(ls.get("ids", [None])[0])
                # find which mrna_id this belongs to
                mrna_for_uid = next(
                    (m for m, u in nucl_uid_map.items() if u == input_uid), None
                )
                if mrna_for_uid is None:
                    continue
                for lsdb in ls.get("linksetdbs", []):
                    prot_links = lsdb.get("links", [])
                    if prot_links:
                        prot_uid_map[mrna_for_uid] = str(prot_links[0])
                        break
        except Exception:
            pass
        time.sleep(delay)

    if not prot_uid_map:
        return {}

    # 3) efetch protein FASTAs in batches
    prot_items = list(prot_uid_map.items())
    uid_to_fasta: dict[str, str] = {}

    for i in range(0, len(prot_items), 50):
        batch = prot_items[i : i + 50]
        uid_list = [uid for _, uid in batch]
        try:
            fasta_text = _http_get(
                f"{_NCBI_BASE}/efetch.fcgi?"
                + _qs(
                    {"db": "protein", "id": ",".join(uid_list),
                     "rettype": "fasta", "retmode": "text"},
                    api_key,
                )
            ).decode("utf-8", errors="replace")

            # Split into individual FASTA records
            records: list[tuple[str, list[str]]] = []
            header = ""
            seq_lines: list[str] = []
            for line in fasta_text.splitlines():
                if line.startswith(">"):
                    if header:
                        records.append((header, seq_lines))
                    header = line
                    seq_lines = []
                elif line.strip():
                    seq_lines.append(line.strip())
            if header:
                records.append((header, seq_lines))

            # Map protein UID → fasta (header encodes the UID as gi|NNN| or accession)
            for j, (mrna_id, uid) in enumerate(batch):
                if j < len(records):
                    hdr, seq = records[j]
                    uid_to_fasta[uid] = hdr + "\n" + "\n".join(seq)
        except Exception:
            pass
        time.sleep(delay)

    # 4) Rename FASTA headers to mRNA IDs for easy post-lookup
    result: dict[str, str] = {}
    for mrna_id, puid in prot_uid_map.items():
        fasta = uid_to_fasta.get(puid, "")
        if fasta:
            lines = fasta.splitlines()
            lines[0] = f">{mrna_id}"
            result[mrna_id] = "\n".join(lines)

    print(f"  [emapper] Retrieved protein sequences for {len(result)}/{len(mrna_ids)} mRNAs")
    return result


# ---------------------------------------------------------------------------
# Step 3–4 – Submit job to eggNOG-mapper and poll
# ---------------------------------------------------------------------------

def submit_emapper_job(fasta_str: str, tax_scope: int = 33208) -> str | None:
    """Submit FASTA sequences to the eggNOG-mapper web API.

    Returns the job ID string, or None on failure.
    """
    payload = json.dumps(
        {
            "seqs": fasta_str,
            "tax_scope": str(tax_scope),
            "ortho_method": "diamond",
            "target_orthologs": "all",
            "go_evidence": "experimental",
        }
    ).encode("utf-8")

    try:
        resp = json.loads(_http_post(f"{_EMAPPER_BASE}/api/job/", payload))
        job_id = resp.get("jobid") or (resp.get("data") or {}).get("jobid")
        if job_id:
            print(f"  [emapper] Job submitted: {job_id}")
        return job_id
    except Exception as exc:
        print(f"  [emapper] Job submission failed: {exc}")
        return None


def poll_emapper_job(
    job_id: str,
    max_wait: int = 7200,
    poll_interval: int = 30,
) -> bool:
    """Poll the eggNOG-mapper API until the job finishes or *max_wait* seconds elapse.

    Returns True on success, False on failure/timeout.
    """
    deadline = time.time() + max_wait
    print(f"  [emapper] Waiting for job {job_id} (max {max_wait // 60} min)…")
    while time.time() < deadline:
        try:
            status = json.loads(_http_get(f"{_EMAPPER_BASE}/api/job/{job_id}/"))
            state = str(status.get("status", "")).upper()
            elapsed = int(time.time() % 60)
            if state in ("FINISHED", "DONE", "SUCCESS", "COMPLETED"):
                print(f"  [emapper] Job {job_id} finished.")
                return True
            if state in ("ERROR", "FAILED", "FAILURE"):
                print(f"  [emapper] Job {job_id} failed: {status.get('message', '')}")
                return False
            print(f"  [emapper] Status: {state}…")
        except Exception as exc:
            print(f"  [emapper] Poll error: {exc}")
        time.sleep(poll_interval)
    print(f"  [emapper] Timed out waiting for job {job_id}")
    return False


# ---------------------------------------------------------------------------
# Step 5 – Fetch annotations and look up GO term names
# ---------------------------------------------------------------------------

def fetch_emapper_annotations(job_id: str) -> dict[str, list[str]]:
    """Fetch the annotations TSV from eggNOG-mapper and return GO IDs per query.

    Returns ``{mrna_id: ["GO:0001234", …]}``.
    """
    try:
        resp = json.loads(_http_get(f"{_EMAPPER_BASE}/api/job/{job_id}/?resultType=annotations"))
        tsv_data: str = resp.get("data") or resp.get("annotations") or ""
    except Exception as exc:
        print(f"  [emapper] Could not fetch annotations: {exc}")
        return {}

    result: dict[str, list[str]] = {}
    for line in tsv_data.splitlines():
        if line.startswith("#") or not line.strip():
            continue
        cols = line.split("\t")
        if len(cols) < 7:
            continue
        query = cols[0]
        go_col = cols[9] if len(cols) > 9 else cols[-1]  # GOs column (index 9)
        if not go_col or go_col in ("-", "NA"):
            continue
        go_ids = [g.strip() for g in go_col.split(",") if g.strip().startswith("GO:")]
        if go_ids:
            result[query] = go_ids

    return result


def lookup_go_terms(
    go_ids: list[str],
) -> dict[str, tuple[str, str]]:
    """Look up GO term names and aspects for a list of GO IDs via QuickGO.

    Returns ``{go_id: (name, aspect)}`` where aspect is one of
    ``"biological_process"``, ``"molecular_function"``, ``"cellular_component"``.
    """
    if not go_ids:
        return {}

    name_aspect: dict[str, tuple[str, str]] = {}
    batch_size = 100
    for i in range(0, len(go_ids), batch_size):
        batch = go_ids[i : i + batch_size]
        ids_str = ",".join(batch)
        url = f"{_QUICKGO_BASE}/{urllib.parse.quote(ids_str)}"
        try:
            resp = json.loads(
                urllib.request.urlopen(  # noqa: S310
                    urllib.request.Request(url, headers={"Accept": "application/json"}),
                    timeout=30,
                ).read()
            )
            for term in resp.get("results", []):
                go_id = term.get("id", "")
                name = term.get("name", "")
                aspect = term.get("aspect", "").lower().replace(" ", "_")
                if go_id and name:
                    name_aspect[go_id] = (name, aspect)
        except Exception:
            # Fall back to storing raw IDs with unknown aspect
            pass
        time.sleep(0.1)

    return name_aspect


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

def run_emapper_annotation(
    mrna_ids: list[str],
    gff_path: Path | None = None,
    api_key: str = "",
    delay: float = 0.4,
) -> dict[str, dict[str, list[str]]]:
    """Full eggNOG-mapper annotation pipeline for a list of mRNA accessions.

    Returns a dict keyed by mRNA ID, each value being::

        {
            "biological_process": ["term name", …],
            "molecular_function": ["term name", …],
            "cellular_component": ["term name", …],
        }

    Returns an empty dict on any unrecoverable error.
    """
    empty: dict[str, dict[str, list[str]]] = {}

    # 1) Tax scope
    if gff_path is not None:
        tax_scope = detect_tax_scope(gff_path, api_key=api_key, delay=delay)
    else:
        print("  [emapper] No GFF provided; defaulting tax_scope=Arthropoda")
        tax_scope = 6656

    # 2) Protein FASTAs
    fasta_map = fetch_protein_fastas(mrna_ids, api_key=api_key, delay=delay)
    if not fasta_map:
        print("  [emapper] No protein sequences retrieved; skipping eggNOG-mapper.")
        return empty

    fasta_str = "\n".join(fasta_map.values())

    # 3) Submit job
    job_id = submit_emapper_job(fasta_str, tax_scope=tax_scope)
    if not job_id:
        return empty

    # 4) Poll
    if not poll_emapper_job(job_id):
        return empty

    # 5) Fetch annotations (GO IDs)
    go_id_map = fetch_emapper_annotations(job_id)
    if not go_id_map:
        print("  [emapper] No GO annotations returned.")
        return empty

    # 6) Look up GO term names and aspects
    all_go_ids: list[str] = sorted({gid for ids in go_id_map.values() for gid in ids})
    print(f"  [emapper] Looking up {len(all_go_ids)} GO terms via QuickGO…")
    term_info = lookup_go_terms(all_go_ids)

    # 7) Build per-mRNA annotation dicts
    result: dict[str, dict[str, list[str]]] = {}
    for mrna_id, go_ids in go_id_map.items():
        bp: list[str] = []
        mf: list[str] = []
        cc: list[str] = []
        for gid in go_ids:
            name, aspect = term_info.get(gid, (gid, ""))
            if aspect == "biological_process":
                bp.append(name)
            elif aspect == "molecular_function":
                mf.append(name)
            elif aspect == "cellular_component":
                cc.append(name)
            # else: unknown aspect — skip
        result[mrna_id] = {
            "biological_process": bp,
            "molecular_function": mf,
            "cellular_component": cc,
        }

    return result

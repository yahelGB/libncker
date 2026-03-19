"""eggNOG-mapper integration (stdlib only) — web API or local installation.

Web mode (default):
  1. Detect the organism's taxonomic scope from the GFF assembly accession
  2. Fetch protein FASTAs for XM_ mRNA IDs via NCBI (fasta_cds_aa)
  3. Submit FASTAs to eggNOG-mapper web API with the detected tax_scope
  4. Poll until done, fetch annotations
  5. Look up GO term names/aspects via QuickGO API
  6. Return per-mRNA GO annotation dicts for all three GO categories

Local mode (--emapper-db):
  Steps 1–2 are the same; instead of the web API, emapper.py is run locally
  using the provided database directory.
"""
from __future__ import annotations

import json
import shutil
import subprocess
import tempfile
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

# Hard-coded GCF/GCA assembly accession → NCBI taxon ID.
# Avoids network calls for known assemblies (useful on HPC with no internet).
_ASSEMBLY_TAXON: dict[str, int] = {
    "GCF_003789085": 6689,   # L. vannamei ASM378908v1
    "GCF_003789085.1": 6689,
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
        print("  [emapper] Could not detect assembly accession; defaulting tax_scope=Arthropoda")
        return 6656

    # 1) Check hard-coded map first — no network needed
    taxon_id = _ASSEMBLY_TAXON.get(accession) or _ASSEMBLY_TAXON.get(accession.split(".")[0])
    if taxon_id:
        scope = _TAXON_SCOPE.get(taxon_id, 33208)
        print(f"  [emapper] Assembly {accession} → taxon {taxon_id} → eggNOG scope {scope} (cached)")
        return scope

    # 2) Fall back to NCBI lookup
    try:
        taxon_id = _accession_to_taxon(accession, api_key, delay)
    except Exception as exc:
        print(f"  [emapper] NCBI taxonomy lookup failed ({exc}); defaulting tax_scope=Arthropoda")
        return 6656

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
# Step 2 – Fetch protein FASTAs for XM_ mRNAs via fasta_cds_aa
# ---------------------------------------------------------------------------

def fetch_protein_fastas(
    mrna_ids: list[str],
    api_key: str = "",
    delay: float = 0.4,
) -> dict[str, str]:
    """Return ``{mrna_id: protein_fasta}`` for each XM_/XR_ accession.

    Uses NCBI efetch with ``rettype=fasta_cds_aa`` directly on the mRNA
    accession — no elink needed. This is reliable for predicted RefSeq
    records (XM_) that often lack nucleotide→protein elink entries.
    The FASTA header is replaced with the mRNA ID so eggNOG-mapper
    returns results keyed by mRNA ID.
    """
    print(f"  [emapper] Fetching protein FASTAs for {len(mrna_ids)} mRNAs…")

    result: dict[str, str] = {}

    # Fetch in batches of 20 (fasta_cds_aa can be large per record)
    for i in range(0, len(mrna_ids), 20):
        batch = mrna_ids[i : i + 20]
        ids_str = ",".join(batch)
        try:
            fasta_text = _http_get(
                f"{_NCBI_BASE}/efetch.fcgi?"
                + _qs(
                    {"db": "nucleotide", "id": ids_str,
                     "rettype": "fasta_cds_aa", "retmode": "text"},
                    api_key,
                )
            ).decode("utf-8", errors="replace")

            # Split multi-FASTA into individual records
            current_header = ""
            current_seq: list[str] = []

            def _save(hdr: str, seq: list[str]) -> None:
                if not hdr or not seq:
                    return
                # Try to match header back to one of the batch accessions
                matched = next(
                    (m for m in batch if m.split(".")[0] in hdr or m in hdr), None
                )
                if matched:
                    result[matched] = f">{matched}\n" + "\n".join(seq)

            for line in fasta_text.splitlines():
                if line.startswith(">"):
                    _save(current_header, current_seq)
                    current_header = line
                    current_seq = []
                elif line.strip():
                    current_seq.append(line.strip())
            _save(current_header, current_seq)

        except Exception as exc:
            print(f"  [emapper] FASTA fetch error for batch {i}: {exc}")
        time.sleep(delay)

    print(f"  [emapper] Retrieved protein sequences for {len(result)}/{len(mrna_ids)} mRNAs")
    return result


# ---------------------------------------------------------------------------
# Step 3–4 – Submit job to eggNOG-mapper and poll
# ---------------------------------------------------------------------------

def submit_emapper_job(fasta_str: str, tax_scope: int = 33208) -> str | None:
    """Submit FASTA sequences to the eggNOG-mapper web API.

    Returns the job ID string, or None on failure.
    Retries up to 3 times with a 10-second pause between attempts.
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

    for attempt in range(1, 4):
        try:
            resp = json.loads(_http_post(f"{_EMAPPER_BASE}/api/job/", payload))
            job_id = resp.get("jobid") or (resp.get("data") or {}).get("jobid")
            if job_id:
                print(f"  [emapper] Job submitted: {job_id}")
                return job_id
        except Exception as exc:
            print(f"  [emapper] Submission attempt {attempt}/3 failed: {exc}")
            if attempt < 3:
                time.sleep(10)

    print(f"  [emapper] All submission attempts failed — eggNOG-mapper server may be busy.")
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


# Mapping from numeric eggNOG scope ID → name accepted by emapper.py --tax_scope
_SCOPE_ID_TO_NAME: dict[int, str] = {
    1: "root",
    2: "Bacteria",
    2157: "Archaea",
    2759: "Eukaryota",
    33208: "Metazoa",
    6656: "Arthropoda",
    50557: "Insecta",
    7742: "Vertebrata",
    33090: "Viridiplantae",
    4751: "Fungi",
}


# ---------------------------------------------------------------------------
# Local eggNOG-mapper runner
# ---------------------------------------------------------------------------

def _run_emapper_local(
    fasta_str: str,
    tax_scope: int,
    emapper_db: Path,
    cpu: int = 4,
    outdir: Path | None = None,
) -> dict[str, list[str]]:
    """Run emapper.py locally and return {mrna_id: [GO_id, …]}.

    Writes a temp FASTA, runs ``emapper.py``, parses the ``.annotations`` file.
    Requires ``emapper.py`` to be on PATH (e.g. conda activate eggnog-mapper).
    The annotations file is also copied to *outdir* for inspection.
    """
    emapper_cmd = shutil.which("emapper.py")
    if not emapper_cmd:
        print("  [emapper] emapper.py not found on PATH; skipping local run.")
        return {}

    # emapper.py --tax_scope expects a name, not a numeric ID
    scope_name = _SCOPE_ID_TO_NAME.get(tax_scope, "auto")

    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)
        fasta_file = tmp / "proteins.fasta"
        fasta_file.write_text(fasta_str)

        cmd = [
            emapper_cmd,
            "-i", str(fasta_file),
            "--itype", "proteins",
            "-o", "emapper_out",
            "--output_dir", str(tmp),
            "--data_dir", str(emapper_db),
            "--tax_scope", scope_name,
            "--cpu", str(cpu),
            "--override",
        ]
        print(f"  [emapper] Running locally (tax_scope={scope_name}, cpu={cpu})…")
        try:
            subprocess.run(cmd, check=True)
        except subprocess.CalledProcessError as exc:
            print(f"  [emapper] Local run failed (exit {exc.returncode}).")
            return {}

        ann_file = tmp / "emapper_out.emapper.annotations"
        if not ann_file.exists():
            print("  [emapper] Annotations file not found after local run.")
            return {}

        tsv = ann_file.read_text()

        # Copy annotations to outdir for inspection
        if outdir is not None:
            dest = outdir / "emapper_out.emapper.annotations"
            dest.write_text(tsv)
            print(f"  [emapper] Annotations saved to {dest}")

        return _parse_annotations_tsv(tsv)


def _parse_annotations_tsv(tsv: str) -> dict[str, list[str]]:
    """Parse eggNOG-mapper annotations TSV → {query: [GO_id, …]}."""
    result: dict[str, list[str]] = {}
    for line in tsv.splitlines():
        if line.startswith("#") or not line.strip():
            continue
        cols = line.split("\t")
        if len(cols) < 10:
            continue
        query = cols[0]
        go_col = cols[9]
        if not go_col or go_col in ("-", "NA"):
            continue
        go_ids = [g.strip() for g in go_col.split(",") if g.strip().startswith("GO:")]
        if go_ids:
            result[query] = go_ids
    return result


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

def run_emapper_annotation(
    mrna_ids: list[str],
    gff_path: Path | None = None,
    api_key: str = "",
    delay: float = 0.4,
    emapper_db: Path | None = None,
    tax_scope_override: int | None = None,
    cpu: int = 4,
    outdir: Path | None = None,
) -> dict[str, dict[str, list[str]]]:
    """Full eggNOG-mapper annotation pipeline for a list of mRNA accessions.

    Returns a dict keyed by mRNA ID, each value being::

        {
            "biological_process": ["term name", …],
            "molecular_function": ["term name", …],
            "cellular_component": ["term name", …],
        }

    When *emapper_db* is set, runs ``emapper.py`` locally instead of the
    web API — useful on HPC nodes with restricted internet or when the
    database is already downloaded.

    Returns an empty dict on any unrecoverable error.
    """
    empty: dict[str, dict[str, list[str]]] = {}

    # 1) Tax scope
    if tax_scope_override is not None:
        tax_scope = tax_scope_override
        print(f"  [emapper] Using explicit tax_scope={tax_scope}")
    elif gff_path is not None:
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

    # 3a) Local mode
    if emapper_db is not None:
        go_id_map = _run_emapper_local(fasta_str, tax_scope, emapper_db, cpu=cpu, outdir=outdir)
        if not go_id_map:
            return empty
    else:
        # 3b) Web API mode
        job_id = submit_emapper_job(fasta_str, tax_scope=tax_scope)
        if not job_id:
            return empty
        if not poll_emapper_job(job_id):
            return empty
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

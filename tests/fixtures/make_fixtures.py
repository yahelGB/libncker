#!/usr/bin/env python3
"""Regenerate the synthetic regression fixtures in this directory.

    python3 tests/fixtures/make_fixtures.py

The fixture is a hand-sized NCBI-RefSeq-shaped GFF3 built to exercise the cases
that break real annotations:

  NW_TEST01  multi-isoform gene (2 mRNAs, one Parent) and an overlapping
             antisense lncRNA on the '-' strand.
  NW_TEST02  a single-gene scaffold (contributes no intergenic distance).
  NW_TEST03  features flush against both scaffold edges (starts at base 1,
             ends at the last base).

Expected geometry, asserted by the tests:
  XR_001.1 (5000-6000) -> 1_U = XM_002.1 at 2.0 kb, 1_D = XM_003.1 at 2.0 kb
  XR_002.1 (8500-9500) overlaps XM_003.1 (8000-9000) by 501 bp -> -0.501 kb
  XR_003.1 (1-300) has no upstream neighbour; 1_D = XM_005.1 at 0.3 kb
"""

from pathlib import Path

HERE = Path(__file__).parent

SCAFFOLDS = [("NW_TEST01.1", 20000), ("NW_TEST02.1", 5000), ("NW_TEST03.1", 1000)]


def attrs(**kw):
    return ";".join(f"{k}={v}" for k, v in kw.items())


def row(seq, src, ftype, start, end, strand, phase, a):
    return "\t".join([seq, src, ftype, str(start), str(end), ".", strand, phase, a])


def build_gff() -> str:
    L = ["##gff-version 3", "#!gff-spec-version 1.21", "#!processor NCBI annotwriter",
         "#!genome-build TEST1", "#!genome-build-accession NCBI_Assembly:GCF_000000000.0",
         "#!annotation-source lncker synthetic test fixture"]
    for name, length in SCAFFOLDS:
        L.append(f"##sequence-region {name} 1 {length}")
    for name, length in SCAFFOLDS:
        L.append(row(name, "RefSeq", "region", 1, length, "+", ".",
                     attrs(ID=f"{name}:1..{length}", Dbxref="taxon:9999", gbkey="Src", mol_type="genomic DNA")))

    s = "NW_TEST01.1"
    # multi-isoform protein-coding gene: two mRNAs sharing one Parent
    L.append(row(s, "Gnomon", "gene", 1000, 3000, "+", ".",
                 attrs(ID="gene-LOC001", Name="LOC001", gbkey="Gene", gene="LOC001", gene_biotype="protein_coding")))
    L.append(row(s, "Gnomon", "mRNA", 1000, 3000, "+", ".",
                 attrs(ID="rna-XM_001.1", Parent="gene-LOC001", gbkey="mRNA", gene="LOC001",
                       product="test protein isoform X1", transcript_id="XM_001.1")))
    L.append(row(s, "Gnomon", "exon", 1000, 1500, "+", ".", attrs(ID="exon-XM_001.1-1", Parent="rna-XM_001.1", gbkey="mRNA")))
    L.append(row(s, "Gnomon", "exon", 2000, 3000, "+", ".", attrs(ID="exon-XM_001.1-2", Parent="rna-XM_001.1", gbkey="mRNA")))
    L.append(row(s, "Gnomon", "CDS", 1000, 1500, "+", "0", attrs(ID="cds-XP_001.1", Parent="rna-XM_001.1", gbkey="CDS")))
    L.append(row(s, "Gnomon", "CDS", 2000, 2900, "+", "0", attrs(ID="cds-XP_001.1-2", Parent="rna-XM_001.1", gbkey="CDS")))
    L.append(row(s, "Gnomon", "mRNA", 1200, 3000, "+", ".",
                 attrs(ID="rna-XM_002.1", Parent="gene-LOC001", gbkey="mRNA", gene="LOC001",
                       product="test protein isoform X2", transcript_id="XM_002.1")))
    L.append(row(s, "Gnomon", "exon", 1200, 1500, "+", ".", attrs(ID="exon-XM_002.1-1", Parent="rna-XM_002.1", gbkey="mRNA")))
    L.append(row(s, "Gnomon", "exon", 2000, 3000, "+", ".", attrs(ID="exon-XM_002.1-2", Parent="rna-XM_002.1", gbkey="mRNA")))
    L.append(row(s, "Gnomon", "CDS", 1200, 1500, "+", "0", attrs(ID="cds-XP_002.1", Parent="rna-XM_002.1", gbkey="CDS")))

    # lncRNA gene, isolated
    L.append(row(s, "Gnomon", "gene", 5000, 6000, "+", ".",
                 attrs(ID="gene-LOC002", Name="LOC002", gbkey="Gene", gene="LOC002", gene_biotype="lncRNA")))
    L.append(row(s, "Gnomon", "lnc_RNA", 5000, 6000, "+", ".",
                 attrs(ID="rna-XR_001.1", Parent="gene-LOC002", gbkey="ncRNA",
                       product="uncharacterized LOC002", transcript_id="XR_001.1")))
    L.append(row(s, "Gnomon", "exon", 5000, 5400, "+", ".", attrs(ID="exon-XR_001.1-1", Parent="rna-XR_001.1", gbkey="ncRNA")))
    L.append(row(s, "Gnomon", "exon", 5600, 6000, "+", ".", attrs(ID="exon-XR_001.1-2", Parent="rna-XR_001.1", gbkey="ncRNA")))

    # protein-coding gene overlapped antisense by a lncRNA
    L.append(row(s, "Gnomon", "gene", 8000, 9000, "+", ".",
                 attrs(ID="gene-LOC003", Name="LOC003", gbkey="Gene", gene="LOC003", gene_biotype="protein_coding")))
    L.append(row(s, "Gnomon", "mRNA", 8000, 9000, "+", ".",
                 attrs(ID="rna-XM_003.1", Parent="gene-LOC003", gbkey="mRNA", gene="LOC003",
                       product="overlapped test protein", transcript_id="XM_003.1")))
    L.append(row(s, "Gnomon", "exon", 8000, 9000, "+", ".", attrs(ID="exon-XM_003.1-1", Parent="rna-XM_003.1", gbkey="mRNA")))
    L.append(row(s, "Gnomon", "CDS", 8000, 8900, "+", "0", attrs(ID="cds-XP_003.1", Parent="rna-XM_003.1", gbkey="CDS")))
    L.append(row(s, "Gnomon", "gene", 8500, 9500, "-", ".",
                 attrs(ID="gene-LOC004", Name="LOC004", gbkey="Gene", gene="LOC004", gene_biotype="lncRNA")))
    L.append(row(s, "Gnomon", "lnc_RNA", 8500, 9500, "-", ".",
                 attrs(ID="rna-XR_002.1", Parent="gene-LOC004", gbkey="ncRNA",
                       product="antisense to LOC003", transcript_id="XR_002.1")))
    L.append(row(s, "Gnomon", "exon", 8500, 9500, "-", ".", attrs(ID="exon-XR_002.1-1", Parent="rna-XR_002.1", gbkey="ncRNA")))

    # single-gene scaffold
    s = "NW_TEST02.1"
    L.append(row(s, "Gnomon", "gene", 1000, 2000, "+", ".",
                 attrs(ID="gene-LOC005", Name="LOC005", gbkey="Gene", gene="LOC005", gene_biotype="protein_coding")))
    L.append(row(s, "Gnomon", "mRNA", 1000, 2000, "+", ".",
                 attrs(ID="rna-XM_004.1", Parent="gene-LOC005", gbkey="mRNA", gene="LOC005",
                       product="lonely test protein", transcript_id="XM_004.1")))
    L.append(row(s, "Gnomon", "exon", 1000, 2000, "+", ".", attrs(ID="exon-XM_004.1-1", Parent="rna-XM_004.1", gbkey="mRNA")))
    L.append(row(s, "Gnomon", "CDS", 1000, 1900, "+", "0", attrs(ID="cds-XP_004.1", Parent="rna-XM_004.1", gbkey="CDS")))

    # scaffold-edge features: starts at base 1, ends at the final base
    s = "NW_TEST03.1"
    L.append(row(s, "Gnomon", "gene", 1, 300, "+", ".",
                 attrs(ID="gene-LOC006", Name="LOC006", gbkey="Gene", gene="LOC006", gene_biotype="lncRNA")))
    L.append(row(s, "Gnomon", "lnc_RNA", 1, 300, "+", ".",
                 attrs(ID="rna-XR_003.1", Parent="gene-LOC006", gbkey="ncRNA",
                       product="edge lncRNA", transcript_id="XR_003.1")))
    L.append(row(s, "Gnomon", "exon", 1, 300, "+", ".", attrs(ID="exon-XR_003.1-1", Parent="rna-XR_003.1", gbkey="ncRNA")))
    L.append(row(s, "Gnomon", "gene", 600, 1000, "+", ".",
                 attrs(ID="gene-LOC007", Name="LOC007", gbkey="Gene", gene="LOC007", gene_biotype="protein_coding")))
    L.append(row(s, "Gnomon", "mRNA", 600, 1000, "+", ".",
                 attrs(ID="rna-XM_005.1", Parent="gene-LOC007", gbkey="mRNA", gene="LOC007",
                       product="edge test protein", transcript_id="XM_005.1")))
    L.append(row(s, "Gnomon", "exon", 600, 1000, "+", ".", attrs(ID="exon-XM_005.1-1", Parent="rna-XM_005.1", gbkey="mRNA")))
    L.append(row(s, "Gnomon", "CDS", 600, 900, "+", "0", attrs(ID="cds-XP_005.1", Parent="rna-XM_005.1", gbkey="CDS")))
    return "\n".join(L) + "\n"


# Three tissues -> three pairwise comparisons -> strict requires N-1 = 2 hits.
INTERSECTS = {
    "AvsB": [
        ("XR_001.1", "Up_A_Down_B"),   # exclusive to A (also up in AvsC)
        ("XR_002.1", "Up_A_Down_B"),   # only one hit -> "up consistent but incomplete"
        ("XR_003.1", "Down_A_Up_B"),   # up in B here, up in C later -> mixed, never exclusive
        ("XM_002.1", "Up_A_Down_B"),   # cis_score 1 for A
        ("XM_003.1", "Up_A_Down_B"),   # cis_score 2 for A (also up in AvsC)
        ("XM_005.1", "Down_A_Up_B"),
    ],
    "AvsC": [
        ("XR_001.1", "Up_A_Down_C"),
        ("XR_003.1", "Down_A_Up_C"),
        ("XM_003.1", "Up_A_Down_C"),
        ("XM_005.1", "Down_A_Up_C"),
    ],
    "BvsC": [
        ("XM_004.1", "Up_B_Down_C"),
        ("XM_005.1", "Down_B_Up_C"),
    ],
}

INTERSECT_COLUMNS = ["ID", "baseMean", "log2FoldChange", "padj", "Expression"]


def build_intersect(rows):
    out = ["\t".join(INTERSECT_COLUMNS)]
    for i, (rid, expr) in enumerate(rows):
        lfc = 2.5 if expr.startswith("Up_") else -2.5
        out.append("\t".join([rid, f"{100 + i * 10}.0", f"{lfc}", "0.001", expr]))
    return "\n".join(out) + "\n"


def main() -> None:
    (HERE / "mini.gff").write_text(build_gff(), encoding="utf-8")
    for tag, rows in INTERSECTS.items():
        (HERE / f"{tag}_intersect.txt").write_text(build_intersect(rows), encoding="utf-8")

    # Contract-violation fixtures.
    (HERE / "no_version.gff").write_text(
        "NW_X.1\tRefSeq\tgene\t1\t100\t.\t+\t.\tID=gene-1\n", encoding="utf-8"
    )
    (HERE / "ensembl_style.gff").write_text(
        "##gff-version 3\n"
        "1\tensembl\tgene\t1\t1000\t.\t+\t.\tID=gene:ENSG01;biotype=protein_coding\n"
        "1\tensembl\ttranscript\t1\t1000\t.\t+\t.\tID=transcript:ENST01;Parent=gene:ENSG01\n"
        "1\tensembl\texon\t1\t1000\t.\t+\t.\tID=exon:ENSE01;Parent=transcript:ENST01\n",
        encoding="utf-8",
    )
    (HERE / "malformed.gff").write_text(
        "##gff-version 3\n"
        "NW_M.1\tRefSeq\tregion\t1\t5000\t.\t+\t.\tID=NW_M.1:1..5000\n"
        "NW_M.1\tGnomon\tgene\t100\t200\t.\t+\t.\tID=gene-OK;gene_biotype=protein_coding\n"
        "NW_M.1\tGnomon\tmRNA\t100\t200\t.\t+\t.\tID=rna-XM_OK.1;Parent=gene-OK;product=ok\n"
        "NW_M.1\tGnomon\texon\t100\t200\t.\t+\t.\tID=exon-1;Parent=rna-XM_OK.1\n"
        "NW_M.1\tGnomon\tCDS\t100\t190\t.\t+\t0\tID=cds-1;Parent=rna-XM_OK.1\n"
        "NW_M.1\tGnomon\tlnc_RNA\t300\t400\t.\t+\t.\tID=rna-XR_OK.1;Parent=gene-OK\n"
        "NW_M.1\tGnomon\tgene\tNOTANUMBER\t500\t.\t+\t.\tID=gene-BAD1\n"
        "NW_M.1\tGnomon\tgene\t900\t800\t.\t+\t.\tID=gene-BAD2\n"
        "NW_M.1\tGnomon\tgene\t0\t50\t.\t+\t.\tID=gene-BAD3\n"
        "NW_M.1\tGnomon\tgene\t100\n",
        encoding="utf-8",
    )
    # Same shape as malformed.gff but every Parent dangles.
    (HERE / "dangling_parent.gff").write_text(
        "##gff-version 3\n"
        "NW_D.1\tGnomon\tgene\t100\t200\t.\t+\t.\tID=gene-OK;gene_biotype=protein_coding\n"
        "NW_D.1\tGnomon\tmRNA\t100\t200\t.\t+\t.\tID=rna-XM_OK.1;Parent=gene-MISSING\n"
        "NW_D.1\tGnomon\texon\t100\t200\t.\t+\t.\tID=exon-1;Parent=rna-XM_OK.1\n"
        "NW_D.1\tGnomon\tCDS\t100\t190\t.\t+\t0\tID=cds-1;Parent=rna-XM_OK.1\n"
        "NW_D.1\tGnomon\tlnc_RNA\t300\t400\t.\t+\t.\tID=rna-XR_OK.1;Parent=gene-OK\n",
        encoding="utf-8",
    )
    print(f"fixtures written to {HERE}")


if __name__ == "__main__":
    main()

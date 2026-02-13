from Bio import Entrez
import time
import csv
import re

# =========================
# 必填：改成你自己的邮箱
# =========================
Entrez.email = "your_email@example.com"

# =========================
# 可选：NCBI API Key（强烈建议）
# 有key可以更稳定、更不容易被限速
# =========================
# Entrez.api_key = "YOUR_NCBI_API_KEY"

# 让NCBI日志里更友好（可选）
Entrez.tool = "TRIM_fetcher"

# Biopython 自带的重试参数（再叠加我们自己的重试）
Entrez.max_tries = 6
Entrez.sleep_between_tries = 5


def esearch_gene_ids(term, retmax=10000):
    """Search Gene IDs by term."""
    handle = Entrez.esearch(db="gene", term=term, retmax=retmax)
    rec = Entrez.read(handle)
    handle.close()
    return rec["IdList"]


def esummary_gene(gene_id):
    """Get gene symbol and description from Gene ID."""
    h = Entrez.esummary(db="gene", id=str(gene_id))
    r = Entrez.read(h)
    h.close()
    doc = r["DocumentSummarySet"]["DocumentSummary"][0]
    symbol = str(doc.get("NomenclatureSymbol") or doc.get("Name") or "")
    desc = str(doc.get("Description") or "")
    return symbol, desc


def elink_gene_to_refseq_rna_uids(gene_id, max_retry=6):
    """
    Link Gene -> nuccore RefSeq RNA UIDs (gene_nuccore_refseqrna).
    Add retry + exponential backoff to handle NCBI transient errors.
    """
    gene_id = str(gene_id)
    for attempt in range(1, max_retry + 1):
        try:
            handle = Entrez.elink(
                dbfrom="gene",
                db="nuccore",
                id=gene_id,
                linkname="gene_nuccore_refseqrna",
            )
            r = Entrez.read(handle)
            handle.close()

            uids = []
            for linkset in r:
                for lsdb in linkset.get("LinkSetDb", []):
                    for l in lsdb.get("Link", []):
                        uids.append(l["Id"])
            return uids

        except Exception as e:
            wait = min(2 ** (attempt - 1), 30)  # 1,2,4,8,16,30...
            print(f"[WARN] elink failed for GeneID={gene_id} (try {attempt}/{max_retry}): {e}")
            time.sleep(wait)

    # 多次失败：返回空，避免整个程序崩掉
    return []


def esummary_nuccore(uids, max_retry=6):
    """
    Get (accession.version, title) from nuccore UIDs via esummary.
    Also add retry to handle transient errors.
    """
    if not uids:
        return []

    uid_str = ",".join(map(str, uids))
    for attempt in range(1, max_retry + 1):
        try:
            h = Entrez.esummary(db="nuccore", id=uid_str)
            r = Entrez.read(h)
            h.close()

            out = []
            for doc in r:
                acc = str(doc.get("Caption") or "")  # accession.version e.g. XM_015278376.4
                title = str(doc.get("Title") or "")
                if acc:
                    out.append((acc, title))
            return out

        except Exception as e:
            wait = min(2 ** (attempt - 1), 30)
            print(f"[WARN] esummary(nuccore) failed (try {attempt}/{max_retry}): {e}")
            time.sleep(wait)

    return []


def pick_best_accession(acc_title_list):
    """
    Pick the best RefSeq mRNA accession:
    1) Prefer isoform X1 / transcript variant X1
    2) If none, prefer NM_ (curated) over XM_ (model)
    """
    if not acc_title_list:
        return ("", "")

    # Prefer X1
    x1 = [
        x for x in acc_title_list
        if re.search(r"(isoform|transcript variant)\s+X1\b", x[1], re.IGNORECASE)
    ]
    if x1:
        x1_nm = [x for x in x1 if x[0].startswith("NM_")]
        return x1_nm[0] if x1_nm else x1[0]

    # Otherwise NM_ then XM_
    nm = [x for x in acc_title_list if x[0].startswith("NM_")]
    if nm:
        return nm[0]
    xm = [x for x in acc_title_list if x[0].startswith("XM_")]
    if xm:
        return xm[0]

    return acc_title_list[0]


def main():
    # 更稳的检索方式：合并两种查询
    term1 = 'TRIM*[Gene] AND "Gallus gallus"[Organism]'
    term2 = '"tripartite motif"[All Fields] AND "Gallus gallus"[Organism]'

    gene_ids = set(esearch_gene_ids(term1))
    gene_ids |= set(esearch_gene_ids(term2))
    gene_ids = list(gene_ids)

    print("Gene hits (raw) =", len(gene_ids))

    rows = []
    kept = 0
    processed = 0

    for gid in gene_ids:
        processed += 1

        try:
            symbol, desc = esummary_gene(gid)
        except Exception as e:
            print(f"[WARN] esummary(gene) failed for GeneID={gid}: {e}")
            continue

        # 只保留 TRIM 开头（避免 term2 引入非 TRIM 基因）
        if not symbol.upper().startswith("TRIM"):
            continue

        kept += 1

        rna_uids = elink_gene_to_refseq_rna_uids(gid)
        best_acc, best_title = "", ""

        if rna_uids:
            acc_title = esummary_nuccore(rna_uids)
            if acc_title:
                best_acc, best_title = pick_best_accession(acc_title)

        rows.append([gid, symbol, desc, best_acc, best_title])

        # 限速：建议更慢一点更稳
        time.sleep(0.7)

        if processed % 20 == 0:
            print(f"Processed {processed}/{len(gene_ids)} | TRIM kept {kept} | rows {len(rows)}")

    out_file = "gallus_gallus_TRIM_refseq_mrna.tsv"
    with open(out_file, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["GeneID", "Symbol", "Description", "Picked_mRNA_Accession", "Picked_Record_Title"])
        w.writerows(rows)

    print(f"Done. Saved: {out_file}")
    print(f"Total TRIM genes kept: {len(rows)}")


if __name__ == "__main__":
    main()

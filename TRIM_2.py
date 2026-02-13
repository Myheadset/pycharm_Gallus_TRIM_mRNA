from Bio import Entrez
import time
import csv

# =========================
# 必填：改成你自己的邮箱
# =========================
Entrez.email = "your_email@example.com"

# =========================
# 可选：NCBI API Key（强烈建议，更稳定/更不容易限速）
# =========================
# Entrez.api_key = "YOUR_NCBI_API_KEY"

Entrez.tool = "TRIM_all_protein_fetcher"
Entrez.max_tries = 6
Entrez.sleep_between_tries = 5


def esearch_gene_ids(term, retmax=10000):
    h = Entrez.esearch(db="gene", term=term, retmax=retmax)
    r = Entrez.read(h)
    h.close()
    return r["IdList"]


def esummary_gene(gene_id):
    h = Entrez.esummary(db="gene", id=str(gene_id))
    r = Entrez.read(h)
    h.close()
    doc = r["DocumentSummarySet"]["DocumentSummary"][0]
    symbol = str(doc.get("NomenclatureSymbol") or doc.get("Name") or "")
    desc = str(doc.get("Description") or "")
    return symbol, desc


def elink_gene_to_protein_uids(gene_id, linkname, max_retry=6):
    gene_id = str(gene_id)
    for attempt in range(1, max_retry + 1):
        try:
            h = Entrez.elink(
                dbfrom="gene",
                db="protein",
                id=gene_id,
                linkname=linkname
            )
            r = Entrez.read(h)
            h.close()

            uids = []
            for linkset in r:
                for lsdb in linkset.get("LinkSetDb", []):
                    for l in lsdb.get("Link", []):
                        uids.append(l["Id"])
            return uids

        except Exception as e:
            wait = min(2 ** (attempt - 1), 30)
            print(f"[WARN] elink({linkname}) failed for GeneID={gene_id} try {attempt}/{max_retry}: {e}")
            time.sleep(wait)

    return []


def esummary_protein(uids, batch_size=200, max_retry=6):
    """
    Batch esummary for protein UIDs.
    Returns list of dicts: {accession, title}
    """
    out = []
    if not uids:
        return out

    for i in range(0, len(uids), batch_size):
        chunk = uids[i:i + batch_size]
        uid_str = ",".join(map(str, chunk))

        for attempt in range(1, max_retry + 1):
            try:
                h = Entrez.esummary(db="protein", id=uid_str)
                r = Entrez.read(h)
                h.close()

                for doc in r:
                    acc = str(doc.get("Caption") or "")  # NP_... or XP_...
                    title = str(doc.get("Title") or "")
                    if acc:
                        out.append({"accession": acc, "title": title})
                break

            except Exception as e:
                wait = min(2 ** (attempt - 1), 30)
                print(f"[WARN] esummary(protein) failed try {attempt}/{max_retry}: {e}")
                time.sleep(wait)

        # 给服务器喘口气
        time.sleep(0.3)

    return out


def main():
    # 两种检索合并更稳：TRIM* + tripartite motif
    term1 = 'TRIM*[Gene] AND "Gallus gallus"[Organism]'
    term2 = '"tripartite motif"[All Fields] AND "Gallus gallus"[Organism]'

    gene_ids = set(esearch_gene_ids(term1))
    gene_ids |= set(esearch_gene_ids(term2))
    gene_ids = list(gene_ids)

    print("Gene hits (raw) =", len(gene_ids))

    rows = []
    processed = 0
    kept = 0
    protein_count = 0

    for gid in gene_ids:
        processed += 1

        try:
            symbol, desc = esummary_gene(gid)
        except Exception as e:
            print(f"[WARN] esummary(gene) failed for GeneID={gid}: {e}")
            continue

        # 只保留 TRIM 开头（严格符合 TRIM 家族）
        if not symbol.upper().startswith("TRIM"):
            continue

        kept += 1

        # 优先 RefSeq protein，若为空再 fallback 到全部 protein
        uids = elink_gene_to_protein_uids(gid, linkname="gene_protein_refseq")
        if not uids:
            uids = elink_gene_to_protein_uids(gid, linkname="gene_protein")

        prot_infos = esummary_protein(uids) if uids else []

        # 每个蛋白一行：不丢任何 isoform
        if prot_infos:
            for p in prot_infos:
                rows.append([gid, symbol, desc, p["accession"], p["title"]])
            protein_count += len(prot_infos)
        else:
            # 即使没有蛋白也记录一行（方便你后续检查）
            rows.append([gid, symbol, desc, "", ""])

        # 限速：更稳
        time.sleep(0.7)

        if processed % 20 == 0:
            print(f"Processed {processed}/{len(gene_ids)} | TRIM genes kept {kept} | proteins {protein_count}")

    out_file = "gallus_gallus_TRIM_all_proteins.tsv"
    with open(out_file, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["GeneID", "Symbol", "Description", "Protein_Accession", "Protein_Title"])
        w.writerows(rows)

    print(f"Done. Saved: {out_file}")
    print(f"Total TRIM genes kept: {kept}")
    print(f"Total protein records (NP_/XP_ etc.): {protein_count}")


if __name__ == "__main__":
    main()

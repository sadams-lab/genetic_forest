import sys

from pysam import VariantFile

def join_p(p):
    return "\t".join(p)

vcf_data = VariantFile(sys.argv[1])
with open(sys.argv[2], "r") as pheno_in:
    pheno_in.readline()
    pheno_lines = {l.strip().split("\t")[0]: [l.strip().split("\t")[1]] for l in pheno_in.readlines()}

variants = []
index = 0
for record in vcf_data.fetch():
    variants.append((str(index), record.id))
    for i, p in pheno_lines.items():
        rec = sum(record.samples[i]["GT"])
        p.append(str(rec) if rec <= 1 else "1")
    index += 1

with open(sys.argv[1] + ".vars", "w") as out_vars:
    for i in variants:
        out_vars.write(join_p(i) + "\n")


with open(sys.argv[1] + ".tsv", "w") as out_data:
    for i, p in pheno_lines.items():
        out_data.write(f"{i}\t{join_p(p)}\n")

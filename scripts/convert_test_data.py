import sys

def join_p(p):
    return "\t".join(p)

def joinvar(v):
    if int(v) == 0:
        return "1"
    elif int(v) == 2:
        return "0"
    else:
        return v

with open(sys.argv[2], "r") as phenos:
    phenos = [p.strip().split("\t")[1] for p in phenos.readlines()]

with open(sys.argv[1], "r") as in_data:
    variants = in_data.readline().strip().split("\t")[6:]
    lines = in_data.readlines()
    with open(sys.argv[1] + ".tsv", "w") as out_data:
        i = 1
        for l in lines:
            nl = [str(i), phenos[i]] + [joinvar(v) for v in l.strip().split("\t")[6:]]
            out_data.write("\t".join(nl) + "\n")
            i += 1

with open(sys.argv[1] + ".vars", "w") as out_vars:
    i = 0
    for v in variants:
        out_vars.write(f"{str(i)}\t{v}\n")
        i += 1


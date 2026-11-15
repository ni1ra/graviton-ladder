"""ODD-SECOND self-check 2: primary sup enclosures (log-derivative bracket, Arb) against the
independent branch-and-bound cell cover (mpmath iv).  Pure decimal comparison, stdlib only.
Also plants a wrong sup (Sodd with sinh(2 pi p) -> sinh(pi p) would give a different number;
here simply the primary Codd is swapped in for Sodd) to show the comparison can fail.
Writes logs/03_compare.log (UTF-8)."""
import json, os, sys
from decimal import Decimal as Dec, getcontext

getcontext().prec = 60
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
P = json.load(open(os.path.join(ROOT, "out", "primary.json"), encoding="utf-8"))["quantities"]
C = json.load(open(os.path.join(ROOT, "out", "sup_cover.json"), encoding="utf-8"))
lines = []
ok = True


def cmp(name, prim, cov):
    plo, phi = Dec(prim["lo"]), Dec(prim["hi"])
    clo, chi = Dec(cov["lo"]), Dec(cov["hi"])
    inside = clo <= plo and phi <= chi          # tight primary must sit inside the wide cover
    rel = (chi - clo) / clo
    lines.append(f"{name:5s} primary [{prim['lo']}, {prim['hi']}]")
    lines.append(f"      cover   [{cov['lo']}, {cov['hi']}]  cover rel width {rel:.2E}  -> "
                 f"{'AGREE (primary inside cover)' if inside else 'MISMATCH'}")
    return inside


for k in ("Sodd", "Codd", "s_3"):
    ok &= cmp(k, P[k], C[k])
lines.append(f"all sup checks agree: {ok}")
plant = not cmp("PLANT", P["Codd"], C["Sodd"])
lines.append(f"planted sup mismatch caught: {plant}")
txt = "\n".join(lines)
print(txt)
open(os.path.join(ROOT, "logs", "03_compare.log"), "w", encoding="utf-8").write(txt + "\n")
sys.exit(0 if (ok and plant) else 1)

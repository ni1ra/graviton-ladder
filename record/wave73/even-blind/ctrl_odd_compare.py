"""Control (i): this lane's copy of the odd pipeline against hinge-blind's stored raw JSON."""
import glob, json, os
HERE = os.path.dirname(os.path.abspath(__file__))
HB = os.path.join(os.path.dirname(HERE), "hinge-blind", "results")
for n, l in ((24, 8), (32, 4), (48, 16)):
    mine = json.load(open(glob.glob(os.path.join(HERE, "results", "raw", "E_n%d_N%d_l%d_p*_odd.json" % (n, 4 * n, l)))[0]))
    ref = json.load(open(glob.glob(os.path.join(HB, "raw", "D2_n%d_N%d_l%d_p*_mid-site.json" % (n, 4 * n, l)))[0]))
    for c in mine["cases"]:
        r = [x for x in ref["cases"] if x["set"] == "mid" and x["g"] == c["g"]][0]
        k = 0
        while k < min(len(c["Q"]), len(r["Q"])) and c["Q"][k] == r["Q"][k]:
            k += 1
        print("n=%d l=%d g=%s  mine=%s  stored=%s  identical leading chars=%d of %d  prec mine/stored=%d/%d" % (
            n, l, c["g"], c["Q"], r["Q"], k, len(r["Q"]), mine["prec"], ref["prec"]))
    p = glob.glob(os.path.join(HB, "pt2", "D2_n%d_N%d_l%d_p*.json" % (n, 4 * n, l)))
    if p:
        q = json.load(open(p[0]))["Q2_mid"]
        print("      Q2_mid mine=%s stored=%s  equal=%s" % (mine["q2"]["odd_mid"], q, mine["q2"]["odd_mid"] == q))

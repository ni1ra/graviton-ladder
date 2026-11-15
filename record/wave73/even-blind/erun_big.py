"""Launcher for the large rows, in stages so that a coherent lattice set exists early and peak
memory stays bounded (about 1.5 GB per n = 192 worker).   python erun_big.py <stage> <workers>"""
import os, sys, time
import erun

STAGES = {
    "A": erun.grid((112, 128, 144)),
    "B192": erun.grid((192,)),
    "B160": erun.grid((160,)),
    "B176": erun.grid((176,)),
    "local192": erun.PLANS["local192"],
}

if __name__ == "__main__":
    import multiprocessing as mp
    jobs = STAGES[sys.argv[1]]
    workers = int(sys.argv[2])
    print("stage %s: %d jobs, %d workers, pid %d" % (sys.argv[1], len(jobs), workers, os.getpid()), flush=True)
    with mp.Pool(workers, maxtasksperchild=1) as pool:
        for line in pool.imap_unordered(erun._job, jobs):
            print(time.strftime("%H:%M:%S"), line, flush=True)
    print("stage %s done" % sys.argv[1], flush=True)

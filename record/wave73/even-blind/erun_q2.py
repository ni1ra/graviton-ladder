"""Addendum-1 item A1: q2-only rows at n = 224 and n = 256.   python erun_q2.py <stage> <workers>"""
import os, sys, time
import erun

STAGES = {
    "224": [(224, 4, l, "q2only") for l in (16, 24, 32, 12, 20, 8, 48, 40, 10)],
    "256": [(256, 4, l, "q2only") for l in (16, 24, 32, 20, 48, 40)],
    "256a": [(256, 4, l, "q2only") for l in (16, 32, 48)],
    "256b": [(256, 4, l, "q2only") for l in (24, 20, 40)],
}

if __name__ == "__main__":
    import multiprocessing as mp
    jobs = STAGES[sys.argv[1]]
    workers = int(sys.argv[2])
    print("stage q2only %s: %d jobs, %d workers, pid %d" % (sys.argv[1], len(jobs), workers, os.getpid()), flush=True)
    with mp.Pool(workers, maxtasksperchild=1) as pool:
        for line in pool.imap_unordered(erun._job, jobs):
            print(time.strftime("%H:%M:%S"), line, flush=True)
    print("stage q2only %s done" % sys.argv[1], flush=True)

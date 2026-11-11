#!/usr/bin/env python3
"""
Diagnostic, not physics (10 Nov 2026): why do 16 workers on the author's machine take about as long as 2 in the sandbox?
Runs one fixed, identical task (a fixed number of Lagrangian evaluations of gravscalar_operator_derive.py on one background)
many times under pools of 1, 2, 4, 8, 16, ... workers and prints, for each pool size:
    startup      seconds from creating the pool to all workers ready (import + building the background, once per worker)
    wall         seconds for the tasks, after startup
    task mean / max   seconds measured INSIDE the workers
    speedup      (tasks x single-worker task time) / wall      -- the ideal is the number of workers
If 'task mean' grows with the pool size, the cores slow each other down (efficiency cores, memory, thermal limits, hyper-threads).
If 'task mean' stays flat but 'speedup' is poor, the time is going to the pool itself (startup, pickling, scheduling).
If one pool size shows task max >> task mean, a few workers sit on slow cores.

Run from graviton_loop (about 2-4 minutes):
    python3 pool_timing.py
Nothing is written to disk; paste the printed table back.
"""
import os, sys, time, random, platform
import sympy as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gravscalar_operator_derive as D

NPROBES = 24          # Lagrangian evaluations per task: about 2 s on one sandbox core
_W = {}


def _init(Rd, v, W, cub):
    t = time.time()
    _W['bg'] = D.Background(Rd, v, W, cub); _W['probes'] = D.probe_fields()[:NPROBES]
    _W['init_s'] = time.time() - t


def _ready(_):
    return os.getpid(), _W['init_s']


def _task(_):
    t = time.time()
    bg = _W['bg']
    base = D.comp_field(D.COMPS[1], D.XR[0] * D.XR[2])
    for kind, comp, key, (ph_h, ph_p) in _W['probes']:
        h = [[base[0][a][b] + D.DL * ph_h[a][b] for b in D.idx] for a in D.idx]
        bg.L2(h, base[1] + D.DL * ph_p, 1)
    return os.getpid(), time.time() - t


def main():
    import multiprocessing as mp
    ncpu = os.cpu_count()
    print(f"# pool_timing: {platform.system()} {platform.machine()}, python {platform.python_version()}, os.cpu_count() = {ncpu}, start method {mp.get_start_method()}")
    print(D.version_line())
    rng = random.Random(1)
    R = D.riemann_sample(rng); v, W, cub = D.scalar_sample(rng)
    # single process, no pool: the reference task time
    _init(R, v, W, cub)
    ref = [_task(0)[1] for _ in range(3)]
    t1 = sum(ref) / len(ref)
    print(f"no pool: building the background {_W['init_s']:.2f} s; one task {t1:.2f} s (three runs: {', '.join(f'{x:.2f}' for x in ref)})")
    print(f"{'workers':>8} {'tasks':>6} {'startup':>8} {'wall':>8} {'task mean':>10} {'task max':>9} {'speedup':>8} {'ideal':>6} {'pids used':>10}")
    sizes = [j for j in (1, 2, 4, 8, 12, 16, 22, 32) if j <= ncpu]
    for jobs in sizes:
        ntask = 3 * jobs
        t0 = time.time()
        with mp.Pool(jobs, initializer=_init, initargs=(R, v, W, cub)) as pool:
            ready = pool.map(_ready, range(4 * jobs), chunksize=1)          # touches every worker (almost surely); returns when they are up
            ts = time.time() - t0
            t0 = time.time()
            res = pool.map(_task, range(ntask), chunksize=1)
            wall = time.time() - t0
        times = [x[1] for x in res]
        print(f"{jobs:>8} {ntask:>6} {ts:>8.1f} {wall:>8.1f} {sum(times) / len(times):>10.2f} {max(times):>9.2f} {ntask * t1 / wall:>8.1f} {jobs:>6} {len({x[0] for x in res}):>10}")
    print("END pool_timing")


if __name__ == '__main__':
    main()

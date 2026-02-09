# Lecture 11 Exercises — Parallel Performance Lab & Wrap-up (1-week workload)

Objective: run a complete scaling study (strong + weak), present results clearly, and reflect on what limits your speedup. Keep your lab log in `codes/lab11/lab-notebook.md`.

## Exercise 1 — Setup and correctness

1. Create `codes/lab11/<name.surname>/` in your sandbox and copy the starter materials from `codes/lab11/python/`.
2. Run serial first and check correctness:
   - `python bench_pi.py --samples 2000000 --processes 1`
3. Record your environment in the lab log (CPU/OS/Python).

## Exercise 2 — Strong scaling

1. Fix the total work (choose a large enough `--samples` so runs last at least a few tenths of a second).
2. Run a sweep and save a CSV:
   - `python run_scaling.py strong --samples 20000000 --max-procs 16 --trials 5 --out strong.csv`
3. Analyze and paste the output table into your lab log:
   - `python analyze_scaling.py strong.csv`

## Exercise 3 — Weak scaling

1. Keep work per process constant and scale total samples with `p`.
2. Run a sweep and save a CSV:
   - `python run_scaling.py weak --samples-per-proc 2000000 --max-procs 16 --trials 5 --out weak.csv`
3. Analyze and paste the output table into your lab log:
   - `python analyze_scaling.py weak.csv`

## Exercise 4 — Interpretation (short write-up)

Write a short interpretation including:

- where efficiency drops and why you think it drops
- whether weak scaling is near-ideal; if not, likely causes
- one “do next” optimization idea you would try (e.g., reduce overhead, tune chunk sizes, reduce communication)

## Exercise 5 — MPI as a parallel launcher for bash scripts

Objective: use MPI to run a bash script **N times in parallel**, and use the rank to split work.

Mini reminder:

- threads: shared memory, locks
- MPI: separate processes (separate memory), message passing, scales across nodes

1. From `codes/lab11/mpi/`, run the “hello from rank” example:
   - `mpiexec -n 4 bash hello_rank.sh`
2. Create a small set of dummy inputs:
   - `mkdir -p inputs && for i in {0..15}; do echo $i > inputs/in_$i.txt; done`
3. Run the parallel splitter:
   - `mpiexec -n 4 bash run_parallel.sh inputs/in_*.txt`
4. In your lab notebook, explain:
   - what “rank” and “world size” mean
   - which environment variables your MPI provides (e.g., `OMPI_COMM_WORLD_RANK` or `PMI_RANK`)
   - how the script assigns tasks to ranks

## Exercise 6 — Wrap-up checklist

Confirm you can:

- explain the difference between strong and weak scaling
- compute and interpret $S(p)$ and $E(p)$
- describe serial vs parallel parts in your program
- run reproducible experiments (commands + parameters + data)

## Submission checklist

- `codes/lab11/<name.surname>/lab-notebook.md`
- `codes/lab11/<name.surname>/strong.csv` and `codes/lab11/<name.surname>/weak.csv`
- your scripts (copied and/or modified from the starter)

# Lecture 11 — Parallel Performance Lab & Wrap-up (3h)

Objectives

- Run practical parallel experiments (e.g., Monte Carlo π)
- Measure runtimes, compute speedup and efficiency
- Prepare final mini-project deliverables and next steps

Content summary

Hands-on lab: serial and parallel implementations (Python multiprocessing, optional C++/OpenMP). Collect timings for 1,2,4,... cores, tabulate speedup, and analyze results vs theoretical expectations. Course summary and integration of tools.

## MPI launcher for parallel bash scripts

In many HPC environments, MPI is available even when you are not writing MPI code. You can use it to **launch multiple copies of a bash script in parallel**:

- MPI is different from multithreading:
  - threads share memory and synchronize with locks
  - MPI uses *separate processes* (separate memory) and coordinates with message passing; it also runs across nodes
- `mpiexec -n N bash script.sh ...` starts **N processes**
- each process has a unique **rank** in `0..N-1`
- rank and world size are exposed via environment variables (implementation-dependent):
  - OpenMPI: `OMPI_COMM_WORLD_RANK`, `OMPI_COMM_WORLD_SIZE`
  - MPICH/PMI: `PMI_RANK`, `PMI_SIZE`

In Lab 11 you will find example scripts under `codes/lab11/mpi/` showing how to:

- print rank/size (`hello_rank.sh`)
- split a list of input files across ranks (`run_parallel.sh`)

<!-- FOOTER START -->
<iframe src="/slideshow/slides11.html" width="100%" height="800px" style="border: none;"></iframe>

---

```{admonition} 🎬 View Slides
:class: tip

**[Open slides in full screen](/slideshow/slides11.html)** for the best viewing experience.
```
<!-- FOOTER END -->

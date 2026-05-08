# EDQC

This repository provides the implementation, experimental scripts, and demo results for EDQC.

Our paper titled **“Cohesive Group Discovery in Interaction Graphs under Explicit Density Constraints”** has been accepted by SIGIR 2026.

## Overview

EDQC is an algorithm for discovering cohesive groups in interaction graphs under explicit density constraints. This repository contains the source code and scripts needed to reproduce the demo experiments included with the repository.

For storage reasons, this repository includes only three demo dataset instances for a lightweight reproducibility check. The full datasets used in the paper are publicly available from the SNAP and the Network Data Repository benchmarks.

## Compilation

Compile the EDQC executable with:

```bash
g++ edqc.cpp -O3 -o edqc
```

## Quick Start

To reproduce the demo results, run:

```bash
bash run_edqc.sh
```

The script executes EDQC on the three included demo instances and reports the corresponding results.

## Datasets

The repository includes three demo dataset instances.

The full experimental evaluation in the paper uses datasets from the SNAP and the Network Data Repository benchmarks. To reproduce the complete experimental results, please download the corresponding datasets and place them in the expected input format (.dimacs) and directory structure used by the scripts.

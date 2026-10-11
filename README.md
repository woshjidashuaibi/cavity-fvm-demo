# Cavity FVM Demo
2D lid-driven cavity solver based on Finite Volume Method (FVM), written in pure C.

## Overview
This is a simple incompressible Navier-Stokes solver for lid-driven cavity flow.
- Spatial discretization: Finite Volume Method
- Pressure-velocity coupling: SIMPLE algorithm
- Collocated grid with Rhie-Chow interpolation
- Support steady-state simulation

## Features
✅ 2D lid-driven cavity classic benchmark
✅ Pure C implementation, no external CFD library dependency
✅ Simple iterative solver
✅ Output velocity and pressure field for post-processing
✅ Run script for Windows (run.bat) / Linux (run.sh)

## Compile & Run
### Windows
```bash
gcc cavity_fvm.c -o cavity_fvm.exe
cavity_fvm.exe

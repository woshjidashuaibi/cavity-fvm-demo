\# cavity-fvm-demo

A minimal single-file C implementation of \*\*Lid-Driven Cavity\*\* using Finite Volume Method(FVM).

Single source file, no external library dependency, compile \& run directly with GCC.



\## ✨ Features

\- Lid-driven cavity flow (classic CFD benchmark)

\- Pure C, only depends on standard library + math library

\- Single `cavity\_fvm.c`, easy to read and modify

\- Output flow field data for post-processing (VTK/dat)



\## 📋 Environment

\- GCC (>=7.5)

\- Linux / WSL / MinGW(Windows)



\## 🚀 Quick Start

```bash

\# clone repository

git clone https://github.com/你的用户名/cavity-fvm-demo.git

cd cavity-fvm-demo



\# compile

gcc cavity\_fvm.c -o cavity\_fvm -lm



\# run

./cavity\_fvm




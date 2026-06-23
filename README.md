# Moment Generator

Moment Generator is a Qt + KaTeX desktop tool for deriving exact finite-sample bias formulas of sample central moments via integer partitions.

Given the moment order `k`, the program enumerates every integer partition of `k` without parts equal to `1`, computes the exact coefficient of each population central-moment product, and renders the formula

```tex
\mathbb{E}[M_k]
=
\sum_{\lambda \vdash k,\;1\notin\lambda}
c_\lambda(n)\prod_q v_q^{m_q}.
```

The first version focuses on formula generation and inspection:

- integer partitions without `1`;
- exact symbolic coefficient formulas in `n`;
- optional numeric coefficient evaluation for a chosen sample size;
- KaTeX preview through Qt WebEngine;
- copyable LaTeX output for paper writing.

## Build

This project is set up for the local Qt installation:

```powershell
D:\Qt\6.9.3\msvc2022_64\bin\qt-cmake.bat -S . -B build -G Ninja
D:\Qt\Tools\CMake_64\bin\cmake.exe --build build
```

If MSVC is not already available in the shell, run the commands from a Visual Studio Developer PowerShell or open the folder in Qt Creator with the `Qt 6.9.3 msvc2022_64` kit.

## Notes

The formula renderer currently loads KaTeX from jsDelivr. A later version can vendor KaTeX assets under `resources/katex/` for offline rendering.

# CompHENS

[![Build Status](https://github.com/avinashresearch1/CompHENS.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/avinashresearch1/CompHENS.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Build Status](https://travis-ci.com/avinashresearch1/CompHENS.jl.svg?branch=main)](https://travis-ci.com/avinashresearch1/CompHENS.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/avinashresearch1/CompHENS.jl?svg=true)](https://ci.appveyor.com/project/avinashresearch1/CompHENS-jl)
[![Coverage](https://codecov.io/gh/avinashresearch1/CompHENS.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/avinashresearch1/CompHENS.jl)
[![Coverage](https://coveralls.io/repos/github/avinashresearch1/CompHENS.jl/badge.svg?branch=main)](https://coveralls.io/github/avinashresearch1/CompHENS.jl?branch=main)
[![DOI](https://zenodo.org/badge/511286580.svg)](https://zenodo.org/badge/latestdoi/511286580)

This software provides a Julia-based toolkit for synthesis of Heat Exchanger Networks (HENs) using a mathematical programming framework. Currently, a sequential algorithm is implemented whereby an LP is formulated to determine the minimum utility consumption, MILP to determine the minimum number of units and the stream matches. Finally, an NLP is formulated to generate the network and calculate the HEN area. 

If you use this toolkit, please cite:
Avinash Subramanian, Flemming Holtorf, Rahul Anantharaman, Truls Gundersen. (2023). CompHENS: Computational Tools for Heat Exchanger Network Synthesis (Version v0.1.0) [Computer software]
DOI: 10.5281/zenodo.7545869

:warning: **Warning** <br>
This package is currently under development. The user interface may change substantially prior to the first stable release.

## No code usage:
1. Download the 2 interface files from: https://github.com/avinashresearch1/CompHENS.jl/tree/main/NoCode_Interface and put them in a folder of your choice (same files from email). 
2. Type in the stream data in the `InputData.xlsx` file. Note that all streams must have a sensible temperature difference (say use a 1 C temperature difference for condensing steam).
3. Download and Install Julia from [Julia website](https://julialang.org/downloads/). The No-Code interface only requires the Julia REPL.
4. From the Julia REPL, access the Package Manager by typing: `]`. Install `CompHENS` and `Pluto`
![image](https://user-images.githubusercontent.com/90404321/217259675-2c48f58c-bd7a-4a86-9d76-1da82989c559.png)
4. Once everything is installed, exit the package manager by typing backspace. Type `using Pluto; Pluto.run()`. This will launch the browser.
5. Navigate to your directory to the `Interface.jl file`. Pluto will launch. 
6. The slider can be used to move the composite curves. The curves update automatically with changing the `DT_min`.

**Note:** It may be necessary to re-run the Pluto notebook in case of errors: Click out of any cell, and `Ctrl+A`followed by `Shift+Enter`.


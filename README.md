# CompHENS

[![Build Status](https://github.com/avinashresearch1/CompHENS.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/avinashresearch1/CompHENS.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Build Status](https://travis-ci.com/avinashresearch1/CompHENS.jl.svg?branch=main)](https://travis-ci.com/avinashresearch1/CompHENS.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/avinashresearch1/CompHENS.jl?svg=true)](https://ci.appveyor.com/project/avinashresearch1/CompHENS-jl)
[![Coverage](https://codecov.io/gh/avinashresearch1/CompHENS.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/avinashresearch1/CompHENS.jl)
[![Coverage](https://coveralls.io/repos/github/avinashresearch1/CompHENS.jl/badge.svg?branch=main)](https://coveralls.io/github/avinashresearch1/CompHENS.jl?branch=main)

This software provides a Julia-based toolkit for synthesis of Heat Exchanger Networks (HENs) using a mathematical programming framework. Currently, a sequential algorithm is implemented whereby an LP is formulated to determine the minimum utility consumption, MILP to determine the minimum number of units and the stream matches. Finally, an NLP is formulated to generate the network and calculate the HEN area. 

If you use this toolkit, please cite:


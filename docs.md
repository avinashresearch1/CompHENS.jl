## Subproblem 4: Network Generation

* Plasmo.jl will be used to model the optimization problem.
* Each stream has a superstructure that is a function of the number of matches for that stream as determined in the match generation subproblem 3. Need to use named nodes for
  Source -> Major splitter -> Minor Mixers -> HXs -> Minor Splitters -> Major Mixer -> Sink. The number of Minor Mixers, HXs and Minor Splitters depend on the number of stream matches. Thus,

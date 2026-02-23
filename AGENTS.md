# Procedure for all cases: 
1. Research:
   - Read this folder and all files in it in depth, understand how it works deeply, what it does and all its specificities. when that’s done, write a detailed report of your learnings and findings in research.md
2. Planning:
   - For the prompt given, write a detailed plan.md document outlining how to implement this. include code snippets. read source files before suggesting changes, base the plan on the actual codebase.
3. Implementation:
   - implement it all. when you’re done with a task or phase, mark it as completed in the plan document. do not stop until all tasks and phases are completed.
   - For each function you add or modify. Add detailed but concise documentation on what you have done. A developer and reviewer should be able to read the code and understand what was removed, modified and added and why.


# How to run Julia Code (relevant for Julia only)
- For Julia work, use a single DaemonMode server per Julia project to avoid repeated precompile costs.
- Start (or ensure) the daemon from the target Julia project directory (the one with `Project.toml`) before running any Julia code: `bash scripts/julia-daemon-start`.
- Run files via the daemon: `bash scripts/julia-daemon-run path/to/file.jl [args...]`.
- Run ad-hoc expressions via the daemon: `bash scripts/julia-daemon-expr 'expr'`.
- Avoid calling `julia` directly unless explicitly requested.
- Stop the daemon when asked: `bash scripts/julia-daemon-stop`.
- If working in a subproject, set `workdir` to that project's root or set `JULIA_PROJECT_DIR` to the project path.
- Port is stored in `.agent/julia-daemon.port` (default 3000); override with `JULIA_DAEMON_PORT`.
- **IMPORTANT**: All Julia commands require `all` permissions because the `juliaup` launcher needs write access to `~/.julia/juliaup/` (lockfile) which is outside the workspace sandbox.
- You must never directly run `julia` commands. Always use the Julia daemon server and client. 

You should never need to create a new Julia depot anywhere. Inform the user if you are stuck. 
You should never need to install another Julia version. Always ask for permission if you think you do.

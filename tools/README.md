# Process tools

Helpers for the research workflow described in `HANDOFF.md` (§2, §4). None of them is used by the experiments or the tests.

| tool | what it does |
|---|---|
| `merge_parts.py` | merges parts of an experiment's raw from private run directories into the repository's raw (then regenerate the md with `--from-raw`, twice, and `cmp`) |
| `mutants.py` | single-point mutation testing of one experiment script against its tests (`pytest -x`) and its `--from-raw` md, each slot a private copy |
| `mutants/` | mutant specs (`(id, what, old, new, k)` tuples): round 84's own list for `exp_width_mismatch.py` and the round-84 referee's |
| `numcheck.py` | the numbers a draft's added lines quote that the results md does not print (drafts quote the md verbatim) |
| `extract_report.py` | extracts a referee subagent's report (from its marker line on) out of the agent's JSONL output into `review/REVIEW_roundNN_<id>.md` |
| `rewrap.py` | rewraps a markdown file's prose paragraphs at width 140 (prefer rewrapping only the changed blocks) |
| `requirements-pinned.txt` | the package versions the current raws were produced with |

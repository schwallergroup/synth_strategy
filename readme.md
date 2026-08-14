# SynthStrategy

Distills synthetic strategies from reaction routes into executable Python functions, then uses
that library to annotate, cluster, and retrieve routes by strategy.

## Install

```bash
git clone https://github.com/schwallergroup/synth_strategy
cd synth_strategy
conda create -n synth_strategy python=3.10 && conda activate synth_strategy
conda install -c conda-forge rdkit
pip install -e .
export OPENROUTER_API_KEY="..."   # only needed for function generation/filtering
```

## API

Apply the function library to a route:

```python
from synth_strategy.api import SynthStrategyAPI

api = SynthStrategyAPI()
annotated = api.annotate_strategies(routes=[route], functions_dir="data/strategy_function_library")
matched = annotated[0]["passing_functions"]  # function name -> detected evidence
```

`route` is an AiZynthFinder-style tree of `mol`/`reaction` nodes:

```json
{"type": "mol", "smiles": "CC(=O)O", "metadata": {"target_smiles": "CC(=O)O"},
 "children": [{"type": "reaction", "metadata": {"mapped_reaction_smiles": "CCO>>CC(=O)O"}, "children": [...]}]}
```

For convergent steps (a reaction node with more than one child branch), you do not need to order
`children` yourself: `annotate_strategies` sorts each node's children by subtree depth (deepest
branch first) before running any function, so the main synthetic line is always index 0 regardless
of input order.

Same thing from the CLI:

```bash
synth-strategy annotate --input-file routes.json --functions-dir data/strategy_function_library --output results/ [--cluster] [--visualize]
synth-strategy cluster --input-file results/annotated_routes.json --functions-dir data/strategy_function_library --output results/
synth-strategy retrieve --query "late-stage Suzuki coupling with an intact piperidine core" --metadata-db data/function_metadata_database.json --route-db data/routes/ --top-k 10
synth-strategy visualize --type clustering --input results/clustering_results.json --annotated_dir data/annotated_routes/
```

`retrieve` rewrites your query into a structured JSON query (`operator`/`queries`/`filters`, with
`AND`/`OR`/`negate`) via an LLM, then runs semantic + exact-match filtering over the function
library. Pass `--query-file query.json` to supply that structure directly instead of natural
language. New route databases just need to be annotated first (`synth-strategy annotate`) before
`--route-db` can point at them.

Add `--config config.yaml` to any command to set defaults instead of repeating flags; see
`--help` on each command for the full flag list.

## Scripts

Generating the function library from scratch — each script's output feeds the next:

- `scripts/01_run_function_generation.sh` — an LLM generates candidate strategy functions from a
  batch of routes. Edit the config block at the top (input routes, model, sample count) before
  running.
- `scripts/02_extract_functions_from_json.py` — unpacks that output into individual `.py` files.
- `scripts/03_run_filtering_pipeline.py` — 4-stage refinement: Gemini-Flash filter → Gemini-Pro
  filter → metadata extraction → return-value normalization. `--source-code-dir` (the output of
  step 2) is required; `--run-stages 1 2 3 4` controls which stages run.

For large-scale annotation/clustering on HPC: `submit_annotation_jobs.sh` and
`submit_clustering_job.sh` (SLURM — edit the `USER CONFIGURATION` block at the top first). Without
SLURM, `run_local.sh` runs the same two steps in parallel on one machine.

## Data

```bash
python scripts/download_uspto_data.py
```

Downloads ~2.5GB of pre-annotated USPTO routes to `data/uspto_st/` (or grab it manually from
[https://doi.org/10.6084/m9.figshare.30146374](https://doi.org/10.6084/m9.figshare.30146374)).

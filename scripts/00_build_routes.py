"""
Build multi-step synthesis routes from single-step reactions, restricted to within-patent links.

This is a standalone re-implementation of the route-construction procedure used to create the
USPTO-Train/Val route sets (originally run with SynTrees, github.com/schwallergroup/SynTrees,
commits 9b60d38/cafbaa8). It follows the within-patent reaction-network approach of Mo et al.
(Chem. Sci. 2021) used for the PaRoutes benchmark (Genheden & Bjerrum, Digital Discovery 2022).

Procedure
  1. Reactions are grouped by patent. Optionally, reactions whose hash is in an exclusion list
     (e.g. the PaRoutes n1/n5 test reactions) are removed first.
  2. Within each reaction, only atom-mapped species (i.e. species contributing atoms to the product)
     count as reactants; unmapped reagents/solvents are ignored. Atom maps are then removed and every
     molecule is written as canonical, stereo-preserving RDKit SMILES. Reactions where a reactant is
     identical to the product are dropped.
  3. Two reactions are linked only if the product of one is string-identical to a reactant of the
     other, and both come from the same patent.
  4. Molecules that are produced but never consumed within a patent are route targets. A depth-first
     search expands each molecule through the patent reaction that produces it; if several
     reactions produce the same molecule, the one giving the deepest sub-route is kept (one reaction
     per molecule, as in the PaRoutes route format). Any route that revisits a molecule along a path
     is a loop and is discarded.
  5. Routes are cleaned: duplicates are removed (and, optionally, routes shorter than
     --min-reactions); within each patent the
     PaRoutes non-overlap rule is applied (drop a route if its target is another route's target or
     intermediate, or its leaves are another route's intermediates; one route per target); routes
     whose largest leaf has >= 40 heavy atoms or > 2x the target's heavy atoms are removed; finally
     the non-overlap rule is applied across all patents.

Input formats (auto-detected, or set --patent-col/--rxn-col/--hash-col):
  * AiZynthTrain/rxnutils-processed USPTO (tab-separated): columns `id` (patent;;idx),
    `rsmi_processed`, `PseudoHash`.
  * Raw USPTO-STEREO (Schwaller et al. 2019, https://ibm.ent.box.com/v/ReactionSeq2SeqDataset,
    US_patents_1976-Sep2016_1product_reactions_*.csv): columns `PatentNumber`, `OriginalReaction`.

Output: JSON Lines, one route per line, in the AiZynthFinder/PaRoutes tree format used by
synth_strategy (reaction nodes carry metadata.mapped_reaction_smiles, so the routes can be passed
straight to SynthStrategyAPI.annotate_strategies).

Examples
  # rebuild the routes for a handful of patents (e.g. for a case study)
  python 00_build_routes.py --input US_patents_1976-Sep2016_1product_reactions_*.csv --output routes.jsonl --patents US08415349B2 US09150547B2

  # full run with the PaRoutes test reactions excluded and a patent-level train/val split
  python 00_build_routes.py --input selected_reactions_all.csv --output-dir routes/ \
      --exclude-hashes paroutes_test_hashes.txt --val-fraction 0.05
"""
import argparse
import hashlib
import json
import os
import random
import re
import sys
from collections import defaultdict
from typing import Dict, List, Optional, Set, Tuple

import pandas as pd
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")
sys.setrecursionlimit(10000)

ATOM_MAP_RE = re.compile(r"\[[^\[\]]+:\d+\]")
MAX_LEAF_HEAVY_ATOMS = 40
MAX_LEAF_TO_ROOT_RATIO = 2.0


# --------------------------------------------------------------------------- molecules / reactions


def canonical_unmapped(smiles: str) -> Optional[str]:
    """Canonical, stereo-preserving SMILES with atom-map numbers removed (None if unparsable)."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(0)
    return Chem.MolToSmiles(mol)


def heavy_atoms(smiles: str) -> int:
    mol = Chem.MolFromSmiles(smiles)
    return mol.GetNumHeavyAtoms() if mol else 0


def inchikey(smiles: str) -> str:
    mol = Chem.MolFromSmiles(smiles)
    return Chem.MolToInchiKey(mol) if mol else smiles


def parse_reaction(rxn_smiles: str) -> Optional[Tuple[List[str], str]]:
    """Return (reactants, product) for an atom-mapped reaction, keeping only mapped reactants."""
    rxn_smiles = rxn_smiles.split(" |")[0].strip()  # drop CXSMILES extensions, e.g. " |f:2.3.4|"
    parts = rxn_smiles.split(">")
    if len(parts) != 3:
        return None
    reactant_part, _, product_part = parts
    reactants = []
    for smi in reactant_part.split("."):
        if smi and ATOM_MAP_RE.search(smi):
            canon = canonical_unmapped(smi)
            if canon:
                reactants.append(canon)
    products = [canonical_unmapped(smi) for smi in product_part.split(".") if smi]
    products = [p for p in products if p]
    if not reactants or not products:
        return None
    product = products[0]
    if product in reactants:
        return None
    return sorted(set(reactants), key=len, reverse=True), product


# --------------------------------------------------------------------------- route extraction


class PatentNetwork:
    def __init__(self, reactions: List[dict]):
        self.producers: Dict[str, List[dict]] = defaultdict(list)
        consumed: Set[str] = set()
        produced: Set[str] = set()
        for rxn in reactions:
            parsed = parse_reaction(rxn["rxn_smiles"])
            if parsed is None:
                continue
            reactants, product = parsed
            self.producers[product].append({**rxn, "reactants": reactants})
            consumed.update(reactants)
            produced.add(product)
        self.targets = sorted(produced - consumed)
        self._depth_cache: Dict[str, int] = {}

    def _depth(self, mol: str, path: frozenset) -> int:
        """Depth (reactions) of the deepest loop-free expansion of `mol`; -1 marks a loop."""
        if mol in path:
            return -1
        if mol not in self.producers:
            return 0
        best = 0
        for rxn in self.producers[mol]:
            sub = [self._depth(r, path | {mol}) for r in rxn["reactants"]]
            if any(d < 0 for d in sub):
                continue
            best = max(best, 1 + max(sub))
        return best

    def build(self, mol: str, path: frozenset = frozenset()) -> Optional[dict]:
        node = {"smiles": mol, "type": "mol", "in_stock": False, "children": []}
        if mol in path:
            return None  # loop
        candidates = []
        for rxn in self.producers.get(mol, []):
            sub = [self._depth(r, path | {mol}) for r in rxn["reactants"]]
            if any(d < 0 for d in sub):
                continue
            candidates.append((1 + max(sub), rxn))
        if not candidates:
            return node  # leaf
        _, rxn = max(candidates, key=lambda c: c[0])  # stable: first of equally deep reactions
        children = [self.build(r, path | {mol}) for r in rxn["reactants"]]
        if any(c is None for c in children):
            return None
        node["children"].append(
            {
                "type": "reaction",
                "smiles": "",
                "metadata": {
                    "mapped_reaction_smiles": rxn["rxn_smiles"],
                    "smiles": rxn["rxn_smiles"],
                    "reaction_hash": rxn.get("rxn_hash"),
                    "retro_template": rxn.get("retro_template"),
                },
                "children": children,
            }
        )
        return node


def route_stats(rt: dict, route_id: str, patent_id: str) -> dict:
    leaves, intermediates, reactions = set(), set(), 0

    def walk(node, is_root):
        nonlocal reactions
        if not node["children"]:
            leaves.add(node["smiles"])
            return 0
        if not is_root:
            intermediates.add(node["smiles"])
        reactions += 1
        return 1 + max(walk(c, False) for c in node["children"][0]["children"])

    llr = walk(rt, True)
    return {
        "rt": rt,
        "id": route_id,
        "patent_id": patent_id,
        "root": inchikey(rt["smiles"]),
        "leaves": sorted({inchikey(s) for s in leaves}),
        "intermediates": sorted({inchikey(s) for s in intermediates}),
        "nreactions": reactions,
        "nleaves": len(leaves),
        "llr": llr,
        "root_size": heavy_atoms(rt["smiles"]),
        "biggest_leave_size": max(heavy_atoms(s) for s in leaves),
    }


# --------------------------------------------------------------------------- cleaning


def non_overlapping(routes: List[dict]) -> List[dict]:
    """PaRoutes non-overlap rule (greedy, in input order), keeping one route per target."""
    taken, kept = set(), []
    for i, route in enumerate(routes):
        if i in taken:
            continue
        taken.add(i)
        overlap = False
        for j in range(i + 1, len(routes)):
            if j in taken:
                continue
            other = routes[j]
            if (
                route["root"] == other["root"]
                or route["root"] in other["intermediates"]
                or other["root"] in route["intermediates"]
                or set(route["leaves"]) & set(other["intermediates"])
                or set(other["leaves"]) & set(route["intermediates"])
            ):
                overlap = True
                taken.add(j)
                break
        if not overlap:
            kept.append(route)
    return list({r["root"]: r for r in kept}.values())


def global_non_overlapping(routes: List[dict]) -> List[dict]:
    """Across-patent pass: drop a route whose target or leaves collide with other routes."""
    inter_count, root_count = defaultdict(int), defaultdict(int)
    for r in routes:
        root_count[r["root"]] += 1
        for m in r["intermediates"]:
            inter_count[m] += 1
    kept = []
    for r in routes:
        own_inter = set(r["intermediates"])
        if root_count[r["root"]] > 1 or inter_count[r["root"]] > (r["root"] in own_inter):
            continue
        if any(inter_count[leaf] > (leaf in own_inter) for leaf in r["leaves"]):
            continue
        kept.append(r)
    return kept


def size_ok(route: dict) -> bool:
    return (
        route["biggest_leave_size"] < MAX_LEAF_HEAVY_ATOMS
        and route["biggest_leave_size"] <= MAX_LEAF_TO_ROOT_RATIO * route["root_size"]
    )


def dedup_key(rt: dict) -> str:
    return hashlib.md5(json.dumps(rt, sort_keys=True).encode()).hexdigest()


# --------------------------------------------------------------------------- driver


def load_reactions(args) -> pd.DataFrame:
    frames = []
    for path in args.input:
        opener = __import__("gzip").open if path.endswith(".gz") else open
        with opener(path, "rt") as f:
            n_comment, header = 0, f.readline()
            while header.startswith("#"):  # USPTO-STEREO files start with two '#' provenance lines
                n_comment, header = n_comment + 1, f.readline()
        sep = "\t" if header.count("\t") > header.count(",") else ","
        frames.append(pd.read_csv(path, sep=sep, skiprows=n_comment, on_bad_lines="skip"))
    df = pd.concat(frames, ignore_index=True)
    cols = set(df.columns)
    patent_col = args.patent_col or ("id" if "id" in cols else "PatentNumber")
    rxn_col = args.rxn_col or next(c for c in ("rsmi_processed", "OriginalReaction", "rsmi") if c in cols)
    hash_col = args.hash_col or ("PseudoHash" if "PseudoHash" in cols else None)
    tmpl_col = args.template_col or ("retro_smarts" if "retro_smarts" in cols else None)
    out = pd.DataFrame(
        {
            "patent_id": df[patent_col].astype(str).str.split(";").str[0],
            "rxn_smiles": df[rxn_col].astype(str),
            "rxn_hash": df[hash_col] if hash_col else None,
            "retro_template": df[tmpl_col] if tmpl_col else None,
        }
    )
    if args.exclude_hashes:
        if hash_col is None:
            sys.exit("--exclude-hashes needs a reaction-hash column (--hash-col)")
        with open(args.exclude_hashes) as f:
            excluded = {line.strip() for line in f if line.strip()}
        n_before = len(out)
        out = out[~out["rxn_hash"].isin(excluded)]
        print(f"Excluded {n_before - len(out)} reactions by hash", flush=True)
    if args.patents:
        out = out[out["patent_id"].isin(set(args.patents))]
    print(f"Loaded {len(out)} reactions from {out['patent_id'].nunique()} patents", flush=True)
    return out


def build_routes_for_patents(grouped: Dict[str, List[dict]], global_pass: bool, min_reactions: int = 1) -> List[dict]:
    routes, seen = [], set()
    n_loops = 0
    for patent_id, reactions in grouped.items():
        network = PatentNetwork(reactions)
        patent_routes = []
        for idx, target in enumerate(network.targets):
            rt = network.build(target)
            if rt is None:
                n_loops += 1
                continue
            key = dedup_key(rt)
            if not rt["children"] or key in seen:
                continue
            seen.add(key)
            rt["patent_id"] = patent_id
            stats = route_stats(rt, f"{patent_id}@{idx}", patent_id)
            if stats["nreactions"] >= min_reactions:
                patent_routes.append(stats)
        patent_routes = [r for r in non_overlapping(patent_routes) if size_ok(r)]
        routes.extend(patent_routes)
    print(f"Routes after within-patent cleaning: {len(routes)} (discarded {n_loops} looped routes)", flush=True)
    if global_pass:
        routes = global_non_overlapping(routes)
        print(f"Routes after cross-patent non-overlap: {len(routes)}", flush=True)
    return routes


def write_jsonl(routes: List[dict], path: str, with_stats: bool) -> None:
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w") as f:
        for r in routes:
            f.write(json.dumps(r if with_stats else r["rt"]) + "\n")
    print(f"Wrote {len(routes)} routes to {path}", flush=True)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", required=True, nargs="+", help="reaction table(s) (CSV/TSV); several files are concatenated")
    p.add_argument("--output", help="output JSONL (single set)")
    p.add_argument("--output-dir", help="write train.jsonl / val.jsonl with a patent-level split")
    p.add_argument("--val-fraction", type=float, default=0.05)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--patents", nargs="+", help="only build routes for these patent IDs")
    p.add_argument("--exclude-hashes", help="file with one reaction hash per line to exclude")
    p.add_argument("--min-reactions", type=int, default=1, help="drop routes with fewer reactions (1 = keep single-step routes, as in the released USPTO-Train/Val sets)")
    p.add_argument("--no-global-overlap", action="store_true", help="skip the cross-patent non-overlap pass")
    p.add_argument("--with-stats", action="store_true", help="write route statistics alongside each tree")
    p.add_argument("--patent-col"), p.add_argument("--rxn-col"), p.add_argument("--hash-col"), p.add_argument("--template-col")
    args = p.parse_args()
    if not (args.output or args.output_dir):
        p.error("give --output or --output-dir")

    df = load_reactions(args)
    grouped = {pid: g.to_dict("records") for pid, g in df.groupby("patent_id", sort=True)}
    global_pass = not args.no_global_overlap

    if args.output:
        write_jsonl(build_routes_for_patents(grouped, global_pass, args.min_reactions), args.output, args.with_stats)
        return

    patent_ids = sorted(grouped)
    random.Random(args.seed).shuffle(patent_ids)
    n_val = int(round(args.val_fraction * len(patent_ids)))
    splits = {"val": patent_ids[:n_val], "train": patent_ids[n_val:]}
    for name, ids in splits.items():
        routes = build_routes_for_patents({pid: grouped[pid] for pid in ids}, global_pass, args.min_reactions)
        write_jsonl(routes, os.path.join(args.output_dir, f"{name}.jsonl"), args.with_stats)


if __name__ == "__main__":
    main()

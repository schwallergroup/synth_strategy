from typing import Tuple, Dict, List
import copy
from rdkit.Chem import AllChem, rdFMCS
from collections import deque
import rdkit.Chem as Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdChemReactions
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
import rdkit.Chem.rdFMCS
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem, rdMolDescriptors
from rdkit.Chem import AllChem, Descriptors, Lipinski
from rdkit.Chem import rdmolops
import re
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import AllChem, Descriptors
import traceback
import rdkit
from collections import Counter
from steerable_retro.utils.check import Check
from steerable_retro.utils import fuzzy_dict, check

from pathlib import Path
root_data = Path(__file__).parent.parent

fg_args = {
    "file_path": f"{root_data}/patterns/functional_groups.json",
    "value_field": "pattern",
    "key_field": "name",
}
reaction_class_args = {
    "file_path": f"{root_data}/patterns/smirks.json",
    "value_field": "smirks",
    "key_field": "name",
}
ring_smiles_args = {
    "file_path": f"{root_data}/patterns/chemical_rings_smiles.json",
    "value_field": "smiles",
    "key_field": "name",
}
functional_groups = fuzzy_dict.FuzzyDict.from_json(**fg_args)
reaction_classes = fuzzy_dict.FuzzyDict.from_json(**reaction_class_args)
ring_smiles = fuzzy_dict.FuzzyDict.from_json(**ring_smiles_args)

checker = check.Check(
    fg_dict=functional_groups, reaction_dict=reaction_classes, ring_dict=ring_smiles
)


HETEROCYCLE_TYPES = [
    "thiazole",
    "oxazole",
    "imidazole",
    "pyrazole",
    "isoxazole",
    "isothiazole",
    "triazole",
    "tetrazole",
    "pyrrole",
    "furan",
    "thiophene",
    "pyridine",
    "pyrimidine",
    "pyrazine",
    "pyridazine",
    "indole",
    "benzothiazole",
    "benzoxazole",
    "benzimidazole",
]

RELATED_HETEROCYCLES = {
    "azoles": [
        "thiazole",
        "oxazole",
        "imidazole",
        "pyrazole",
        "isoxazole",
        "isothiazole",
        "triazole",
        "tetrazole",
    ],
    "benzazoles": ["benzothiazole", "benzoxazole", "benzimidazole"],
    "pyridines": ["pyridine", "pyrimidine", "pyrazine", "pyridazine"],
    "five_membered": ["pyrrole", "furan", "thiophene", "indole"],
}

def main(route) -> Tuple[bool, Dict]:
    """
    Detects syntheses that involve the formation of multiple heterocyclic rings from a predefined list.
    """
    findings_template = {
        "atomic_checks": {
            "named_reactions": [],
            "ring_systems": [],
            "functional_groups": []
        },
        "structural_constraints": []
    }
    findings_json = copy.deepcopy(findings_template)

    heterocycle_formations = []

    def dfs_traverse(reaction, depth=0):
        if reaction["type"] == "reaction":
            if "rsmi" in reaction.get("metadata", {}):
                rsmi = reaction["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0]
                product_smiles = rsmi.split(">")[-1]

                try:
                    for heterocycle in HETEROCYCLE_TYPES:
                        if checker.check_ring(heterocycle, product_smiles):
                            reactants_list = reactants_smiles.split(".")
                            reactant_has_ring = any(
                                checker.check_ring(heterocycle, reactant)
                                for reactant in reactants_list
                            )

                            if not reactant_has_ring:
                                heterocycle_formations.append(
                                    {"type": heterocycle, "depth": depth, "rsmi": rsmi}
                                )
                                # Record atomic check for ring system formation
                                if heterocycle not in findings_json["atomic_checks"]["ring_systems"]:
                                    findings_json["atomic_checks"]["ring_systems"].append(heterocycle)
                                # Record atomic check for ring_formation (named_reaction)
                                if "ring_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                                    findings_json["atomic_checks"]["named_reactions"].append("ring_formation")
                except Exception as e:
                    print(f"Error processing reaction SMILES: {e}")

        for child in reaction.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if reaction["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    dfs_traverse(route)

    heterocycle_families = set()
    for formation in heterocycle_formations:
        for family, members in RELATED_HETEROCYCLES.items():
            if formation["type"] in members:
                heterocycle_families.add(family)
                break

    result = len(heterocycle_formations) >= 2 or len(heterocycle_families) >= 2

    # Record structural constraints
    if len(heterocycle_formations) >= 2:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "heterocycle_formation",
                "operator": ">=",
                "value": 2,
                "options": {
                    "condition": "A"
                }
            }
        })
    if len(heterocycle_families) >= 2:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "heterocycle_family_formation",
                "operator": ">=",
                "value": 2,
                "options": {
                    "condition": "B"
                }
            }
        })

    return result, findings_json

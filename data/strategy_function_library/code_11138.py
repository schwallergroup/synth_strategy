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
from synth_strategy.utils.check import Check
from synth_strategy.utils import fuzzy_dict, check

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


# Refactoring for Enumeration: Isolate the list of reactions.
HETEROCYCLE_FORMATION_REACTIONS = [
    "Formation of NOS Heterocycles",
    "Paal-Knorr pyrrole synthesis",
    "benzimidazole_derivatives_carboxylic-acid/ester",
    "benzimidazole_derivatives_aldehyde",
    "benzothiazole",
    "benzoxazole_arom-aldehyde",
    "benzoxazole_carboxylic-acid",
    "thiazole",
    "tetrazole_terminal",
    "tetrazole_connect_regioisomere_1",
    "tetrazole_connect_regioisomere_2",
    "1,2,4-triazole_acetohydrazide",
    "1,2,4-triazole_carboxylic-acid/ester",
    "pyrazole",
    "Fischer indole",
    "benzofuran",
    "benzothiophene",
    "indole",
    "oxadiazole",
    "imidazole",
]

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects a synthesis route that forms a complex heterocycle
    in the middle of the synthesis rather than using it as a starting material.
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

    heterocycle_formation_midroute = False

    def dfs_traverse(node, depth=0):
        nonlocal heterocycle_formation_midroute, findings_json

        if node["type"] == "reaction" and depth > 0:
            rsmi = node["metadata"].get("rsmi", "")
            if not rsmi:
                return

            # Check for heterocycle formation reaction
            for reaction_type in HETEROCYCLE_FORMATION_REACTIONS:
                if checker.check_reaction(reaction_type, rsmi):
                    heterocycle_formation_midroute = True
                    findings_json["atomic_checks"]["named_reactions"].append(reaction_type)
                    # This structural constraint is met if any heterocycle formation reaction is found at depth > 0
                    findings_json["structural_constraints"].append({
                        "type": "positional",
                        "details": {
                            "target": "any_heterocycle_formation_reaction",
                            "position": "not_last_stage"
                        }
                    })
                    return

        # Determine the next depth based on the current node's type
        next_depth = depth
        if node["type"] != "reaction": # If current node is chemical, depth increases
            next_depth = depth + 1

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, next_depth)

    # Start traversal from the root
    dfs_traverse(route)

    return heterocycle_formation_midroute, findings_json

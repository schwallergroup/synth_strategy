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

HETEROCYCLES_OF_INTEREST = [
    "furan",
    "pyrrole",
    "pyridine",
    "pyrazole",
    "imidazole",
    "oxazole",
    "thiazole",
    "pyrimidine",
    "pyrazine",
    "pyridazine",
    "triazole",
    "tetrazole",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if specific heterocycles found in the target molecule are also present in all intermediate precursors. The heterocycles checked are defined in the HETEROCYCLES_OF_INTEREST list: furan, pyrrole, pyridine, pyrazole, imidazole, oxazole, thiazole, pyrimidine, pyrazine, pyridazine, triazole, and tetrazole.
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
    
    result = True

    def get_heterocycles(mol_smiles):
        heterocycles = []
        for ring_name in HETEROCYCLES_OF_INTEREST:
            if checker.check_ring(ring_name, mol_smiles):
                heterocycles.append(ring_name)
                # Record detected heterocycles in findings_json
                if ring_name not in findings_json["atomic_checks"]["ring_systems"]:
                    findings_json["atomic_checks"]["ring_systems"].append(ring_name)
        return heterocycles

    # First identify target heterocycles
    target_heterocycles = []
    if route["type"] == "mol":
        target_heterocycles = get_heterocycles(route["smiles"])
        print(f"Target molecule contains heterocycles: {target_heterocycles}")

    # If no heterocycles in target, strategy doesn't apply
    if not target_heterocycles:
        print("No heterocycles in target molecule")
        return True, findings_json # Return True as the condition for the strategy is not met

    def dfs(node, depth=0):
        nonlocal result
        if node["type"] == "mol" and not node.get("in_stock", False):
            current_heterocycles = get_heterocycles(node["smiles"])

            # Check if all target heterocycles are present
            for cycle in target_heterocycles:
                if cycle not in current_heterocycles:
                    print(f"Heterocycle {cycle} not preserved at depth {depth}")
                    result = False # Set the main result to False
                    # Add the structural constraint for negation of ring_formation
                    if {"type": "negation", "details": {"target": "ring_formation"}} not in findings_json["structural_constraints"]:
                        findings_json["structural_constraints"].append({"type": "negation", "details": {"target": "ring_formation"}})
                    return False

        # Continue DFS traversal
        all_preserved = True
        for child in node.get("children", []):
            # Determine the new depth based on the current node's type
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            
            if not dfs(child, new_depth):
                all_preserved = False

        return all_preserved

    result = dfs(route)
    return result, findings_json

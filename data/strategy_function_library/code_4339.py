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

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if specific aromatic heterocycles (thiazole, pyridine) are present
    in the final product and preserved throughout the entire synthesis.
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

    # Track heterocycles in final product
    final_product_heterocycles = {"thiazole": False, "pyridine": False}

    # Check if the final product contains any of the target heterocycles
    if route["type"] == "mol":
        final_product_smiles = route["smiles"]
        if checker.check_ring("thiazole", final_product_smiles):
            final_product_heterocycles["thiazole"] = True
            findings_json["atomic_checks"]["ring_systems"].append("thiazole")
        if checker.check_ring("pyridine", final_product_smiles):
            final_product_heterocycles["pyridine"] = True
            findings_json["atomic_checks"]["ring_systems"].append("pyridine")

    # If no heterocycles in final product, strategy doesn't apply
    if not any(final_product_heterocycles.values()):
        return False, findings_json

    # Track if heterocycles are modified
    heterocycles_modified = False

    def dfs_traverse(node):
        nonlocal heterocycles_modified, findings_json

        if node["type"] == "mol" and "children" not in node:
            # This is a leaf node (starting material)
            return

        if node["type"] == "mol" and "children" in node:
            mol_smiles = node["smiles"]
            product_has_thiazole = checker.check_ring("thiazole", mol_smiles)
            product_has_pyridine = checker.check_ring("pyridine", mol_smiles)

            # Check each reaction that produced this molecule
            for child in node.get("children", []):
                if child["type"] == "reaction":
                    rsmi = child["metadata"]["rsmi"]
                    reactants_smiles = rsmi.split(">")[0].split(".")

                    # Check if heterocycles in product are preserved in reactants
                    reactant_has_thiazole = any(
                        checker.check_ring("thiazole", smi) for smi in reactants_smiles
                    )
                    reactant_has_pyridine = any(
                        checker.check_ring("pyridine", smi) for smi in reactants_smiles
                    )

                    # Check if any heterocycle was lost (in retrosynthetic direction)
                    if (
                        (product_has_thiazole and not reactant_has_thiazole)
                        or (product_has_pyridine and not reactant_has_pyridine)
                    ):
                        heterocycles_modified = True

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Strategy is found if final product has heterocycles and they were preserved
    has_heterocycles = any(final_product_heterocycles.values())
    preserved = not heterocycles_modified

    # Add structural constraints based on findings
    if has_heterocycles:
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "thiazole",
                    "pyridine"
                ],
                "note": "At least one of the target rings must be present in the final product for the strategy to apply."
            }
        })
    
    if preserved:
        findings_json["structural_constraints"].append({
            "type": "negation",
            "details": {
                "target": "ring_formation",
                "scope": [
                    "thiazole",
                    "pyridine"
                ]
            }
        })

    return has_heterocycles and preserved, findings_json
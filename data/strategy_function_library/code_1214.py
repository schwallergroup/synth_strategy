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


PYRIDINE_FUSED_SYSTEMS = ["quinoline", "isoquinoline"]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if the synthetic route involves building complexity on a
    pyridine-containing scaffold by forming specific fused ring systems.
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

    # Track if we found a reaction that builds a heterocycle on a pyridine scaffold
    heterocycle_elaboration_found = False

    def dfs_traverse(node, depth=0):
        nonlocal heterocycle_elaboration_found, findings_json

        # Check reaction nodes
        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants_part = rsmi.split(">")[0]
            product_part = rsmi.split(">")[-1]

            # Check if reactants contain pyridine
            reactant_has_pyridine = False
            for reactant in reactants_part.split("."):
                if checker.check_ring("pyridine", reactant):
                    reactant_has_pyridine = True
                    if "pyridine" not in findings_json["atomic_checks"]["ring_systems"]:
                        findings_json["atomic_checks"]["ring_systems"].append("pyridine")
                    break

            # Check if product has a fused system from the list
            product_has_fused_system = False
            for system in PYRIDINE_FUSED_SYSTEMS:
                if checker.check_ring(system, product_part):
                    product_has_fused_system = True
                    if system not in findings_json["atomic_checks"]["ring_systems"]:
                        findings_json["atomic_checks"]["ring_systems"].append(system)

            # If the reaction starts with pyridine and results in a product with a fused system
            if reactant_has_pyridine and product_has_fused_system:
                heterocycle_elaboration_found = True
                # Add the structural constraint
                findings_json["structural_constraints"].append({
                    "type": "sequence",
                    "details": {
                        "before": [
                            "pyridine"
                        ],
                        "after": [
                            "quinoline",
                            "isoquinoline"
                        ],
                        "scope": "reaction_step"
                    }
                })

        # Continue traversing
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    # Start traversal from the root
    dfs_traverse(route)

    return heterocycle_elaboration_found, findings_json

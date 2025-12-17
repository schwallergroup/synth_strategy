from typing import Tuple, Dict, List
import copy
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
    Detects a synthetic strategy where a nitrile group is maintained through
    multiple steps before being transformed in a late stage.
    """
    findings_template = {
        "atomic_checks": {"named_reactions": [], "ring_systems": [], "functional_groups": []},
        "structural_constraints": []
    }
    findings_json = copy.deepcopy(findings_template)

    # Track nitrile presence through synthesis
    nitrile_steps = []
    final_nitrile_transformation = False
    total_steps = 0

    def dfs_traverse(node, depth=0):
        nonlocal nitrile_steps, final_nitrile_transformation, total_steps, findings_json

        if node["type"] == "reaction":
            total_steps += 1
            if "mapped_reaction_smiles" in node.get("metadata", {}):
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_str = rsmi.split(">")[0]
                reactants = reactants_str.split(".")
                product_str = rsmi.split(">")[-1]

                # Check if nitrile is present in any reactant
                reactant_has_nitrile = False
                for reactant in reactants:
                    if checker.check_fg("Nitrile", reactant):
                        reactant_has_nitrile = True
                        if "Nitrile" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Nitrile")
                        break

                # Check if nitrile is present in product
                product_has_nitrile = checker.check_fg("Nitrile", product_str)
                if product_has_nitrile:
                    if "Nitrile" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Nitrile")
                    nitrile_steps.append(depth)
                    print(f"Detected nitrile in product at depth {depth}")

                # Check if nitrile is transformed in this step (present in reactant but not in product)
                if reactant_has_nitrile and not product_has_nitrile:
                    # Late stage corresponds to low depth in retrosynthetic traversal
                    if depth <= 1:  # Consider depth 0 or 1 as late-stage
                        final_nitrile_transformation = True
                        print(f"Detected final nitrile transformation at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    # Start traversal from the root
    dfs_traverse(route)

    # Strategy criteria: nitrile present in at least 2 steps and transformed in final step
    strategy_present = len(nitrile_steps) >= 2 and final_nitrile_transformation

    if len(nitrile_steps) >= 2:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "reaction_with_product_functional_group",
                "on_group": "Nitrile",
                "operator": ">=",
                "value": 2
            }
        })

    if final_nitrile_transformation:
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "functional_group_transformation",
                "on_group": "Nitrile",
                "position": "late_stage"
            }
        })

    if strategy_present:
        print(f"Nitrile-directed synthesis detected: maintained for {len(nitrile_steps)} steps")
    else:
        print(
            f"Nitrile-directed synthesis NOT detected: nitrile steps = {len(nitrile_steps)}, final transformation = {final_nitrile_transformation}"
        )

    return strategy_present, findings_json

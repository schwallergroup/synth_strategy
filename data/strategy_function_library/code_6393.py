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


HETEROCYCLES_FOR_CF3_INSTALLATION = [
    "pyridine",
    "pyrimidine",
    "pyrazine",
    "pyridazine",
    "triazine",
    "furan",
    "thiophene",
    "pyrrole",
    "imidazole",
    "oxazole",
    "thiazole",
    "indole",
    "benzimidazole",
    "quinoline",
    "isoquinoline",
]

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects a strategy involving early installation of a CF3 group
    on an aromatic heterocycle.
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

    has_cf3_installation = False
    early_step_threshold = 2  # Consider depth >= 2 as early steps

    def dfs_traverse(node, depth=0):
        nonlocal has_cf3_installation, findings_json

        if node["type"] == "reaction" and depth >= early_step_threshold:
            # Extract reactants and product
            rsmi = node["metadata"].get("rsmi", "")
            if not rsmi:
                return

            product_smiles = rsmi.split(">")[-1]

            try:
                # Check for CF3 group in product
                has_cf3 = checker.check_fg("Trifluoro group", product_smiles)
                if has_cf3:
                    if "Trifluoro group" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Trifluoro group")

                # Check if product contains a heterocycle
                has_heterocycle = False
                detected_heterocycles = []
                for ring in HETEROCYCLES_FOR_CF3_INSTALLATION:
                    if checker.check_ring(ring, product_smiles):
                        has_heterocycle = True
                        if ring not in findings_json["atomic_checks"]["ring_systems"]:
                            findings_json["atomic_checks"]["ring_systems"].append(ring)
                        detected_heterocycles.append(ring)

                # Check for fluorination reactions
                is_fluorination = checker.check_reaction("Fluorination", rsmi)
                if is_fluorination:
                    if "Fluorination" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("Fluorination")

                if has_cf3 and has_heterocycle and is_fluorination:
                    has_cf3_installation = True
                    # Record structural constraint: co-occurrence
                    findings_json["structural_constraints"].append({
                        "type": "co-occurrence",
                        "details": {
                            "description": "A single reaction step must be a Fluorination reaction whose product contains both a Trifluoro group and one of the specified heterocycles.",
                            "targets": [
                                "Fluorination",
                                "Trifluoro group",
                                "one_of_specified_heterocycles"
                            ]
                        }
                    })
                    # Record structural constraint: positional
                    findings_json["structural_constraints"].append({
                        "type": "positional",
                        "details": {
                            "description": "The reaction installing the CF3 group must occur at least two steps before the final product.",
                            "target": "Fluorination",
                            "position": "depth >= 2"
                        }
                    })

            except Exception as e:
                print(f"Error processing reaction at depth {depth}: {e}")

        # Process children
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    return has_cf3_installation, findings_json
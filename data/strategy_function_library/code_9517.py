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
    This function detects if cyano and chloro functional groups are present and
    maintained throughout the synthesis.
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

    # Initialize tracking
    molecules_with_cyano = []
    molecules_with_chloro = []
    starting_materials = []
    final_product = None

    def dfs_traverse(node, depth=0):
        nonlocal final_product, molecules_with_cyano, molecules_with_chloro, starting_materials, findings_json

        if node["type"] == "mol":
            mol_smiles = node["smiles"]

            # Store final product (depth 0)
            if depth == 0:
                final_product = mol_smiles

            # Store starting materials
            if node.get("in_stock", False):
                starting_materials.append((mol_smiles, depth))

            # Check for cyano groups
            if checker.check_fg("Nitrile", mol_smiles):
                molecules_with_cyano.append((mol_smiles, depth))
                if "Nitrile" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Nitrile")

            # Check for chloro groups
            if checker.check_fg(["Alkyl chloride", "Aryl chloride"], mol_smiles):
                molecules_with_chloro.append((mol_smiles, depth))
                if "Alkyl chloride" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Alkyl chloride")
                if "Aryl chloride" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Aryl chloride")

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # Depth increases only if current node is not 'reaction'
                new_depth = depth + 1
            dfs_traverse(child, new_depth)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if both functional groups are present in the final product
    final_has_cyano = any(final_product == mol for mol, _ in molecules_with_cyano)
    final_has_chloro = any(final_product == mol for mol, _ in molecules_with_chloro)

    if final_has_cyano:
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "Nitrile",
                "position": "final_product"
            }
        })
    if final_has_chloro:
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": [
                    "Alkyl chloride",
                    "Aryl chloride"
                ],
                "logic": "OR",
                "position": "final_product"
            }
        })

    # Check if both functional groups are present in at least one starting material
    starting_has_cyano = any(
        mol in [sm for sm, _ in starting_materials] for mol, _ in molecules_with_cyano
    )
    starting_has_chloro = any(
        mol in [sm for sm, _ in starting_materials] for mol, _ in molecules_with_chloro
    )

    # Calculate coverage through the synthesis
    all_depths = set(d for _, d in molecules_with_cyano).union(
        set(d for _, d in molecules_with_chloro)
    )
    max_depth = max(all_depths) if all_depths else 0

    cyano_depths = set(d for _, d in molecules_with_cyano)
    chloro_depths = set(d for _, d in molecules_with_chloro)

    cyano_coverage = len(cyano_depths) / (max_depth + 1) if max_depth > 0 else 0
    chloro_coverage = len(chloro_depths) / (max_depth + 1) if max_depth > 0 else 0

    print(f"Molecules with cyano: {len(molecules_with_cyano)}")
    print(f"Molecules with chloro: {len(molecules_with_chloro)}")
    print(f"Starting materials: {len(starting_materials)}")
    print(f"Max depth: {max_depth}")
    print(f"Cyano coverage: {cyano_coverage:.2f}")
    print(f"Chloro coverage: {chloro_coverage:.2f}")

    # Determine if both groups are maintained throughout synthesis
    # Criteria:
    # 1. Present in final product
    # 2. Present in at least one starting material OR have good coverage (>50%)
    cyano_maintained = final_has_cyano and (starting_has_cyano or cyano_coverage >= 0.5)
    chloro_maintained = final_has_chloro and (starting_has_chloro or chloro_coverage >= 0.5)

    if cyano_maintained:
        findings_json["structural_constraints"].append({
            "type": "conditional_origin",
            "details": {
                "target": "Nitrile",
                "logic": "OR",
                "conditions": [
                    {
                        "type": "positional",
                        "position": "starting_material"
                    },
                    {
                        "type": "coverage",
                        "operator": ">=",
                        "value": 0.5
                    }
                ]
            }
        })
    if chloro_maintained:
        findings_json["structural_constraints"].append({
            "type": "conditional_origin",
            "details": {
                "target": [
                    "Alkyl chloride",
                    "Aryl chloride"
                ],
                "target_logic": "OR",
                "logic": "OR",
                "conditions": [
                    {
                        "type": "positional",
                        "position": "starting_material"
                    },
                    {
                        "type": "coverage",
                        "operator": ">=",
                        "value": 0.5
                    }
                ]
            }
        })

    print(f"Cyano maintained: {cyano_maintained}")
    print(f"Chloro maintained: {chloro_maintained}")

    result = cyano_maintained and chloro_maintained

    # Both groups must be maintained
    if result:
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "description": "The maintenance conditions for both cyano and chloro groups must be met simultaneously.",
                "targets": [
                    "Nitrile",
                    [
                        "Alkyl chloride",
                        "Aryl chloride"
                    ]
                ]
            }
        })

    return result, findings_json

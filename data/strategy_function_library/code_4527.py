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


def main(route) -> Tuple[bool, Dict]:
    """
    This function detects if the synthesis includes at least one ring formation step.
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

    found_ring_formation = False

    def dfs_traverse(node, depth=0):
        nonlocal found_ring_formation, findings_json

        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]

            # Check by counting rings
            try:
                reactants_part = rsmi.split(">")[0]
                product_part = rsmi.split(">")[-1]

                reactants = [Chem.MolFromSmiles(r) for r in reactants_part.split(".")]
                product = Chem.MolFromSmiles(product_part)

                if all(r is not None for r in reactants) and product is not None:
                    # Count rings in reactants and product
                    reactant_rings = sum(r.GetRingInfo().NumRings() for r in reactants)
                    product_rings = product.GetRingInfo().NumRings()

                    if product_rings > reactant_rings:
                        found_ring_formation = True
                        findings_json["atomic_checks"]["named_reactions"].append("ring_formation")
                        # Add the structural constraint if not already present
                        constraint_obj = {
                            "type": "count",
                            "details": {
                                "target": "ring_formation",
                                "operator": ">=",
                                "value": 1
                            }
                        }
                        if constraint_obj not in findings_json["structural_constraints"]:
                            findings_json["structural_constraints"].append(constraint_obj)
                        print(
                            f"Detected ring formation by ring count at {node.get('metadata', {}).get('ID', '')}"
                        )
            except Exception as e:
                print(f"Error processing reaction: {e}")

        for child in node.get("children", []):
            new_depth = depth
            if node["type"] != "reaction": # From chemical to reaction
                new_depth = depth + 1
            # Else (from reaction to chemical), new_depth remains 'depth'
            dfs_traverse(child, new_depth)

    dfs_traverse(route)
    return found_ring_formation, findings_json

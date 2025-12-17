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
    Detects if a cyano group is preserved throughout the synthesis.

    This function checks if a cyano group present in the final product
    is preserved (not created or modified) throughout the synthesis route.
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

    # Track if we've found a cyano group that's preserved
    final_product_has_cyano = False
    cyano_created_in_synthesis = False

    def dfs_traverse(node, depth=0):
        nonlocal final_product_has_cyano, cyano_created_in_synthesis, findings_json

        # For molecule nodes
        if node["type"] == "mol":
            has_cyano = checker.check_fg("Nitrile", node["smiles"])

            # If this is the final product (depth 0)
            if depth == 0:
                if has_cyano:
                    final_product_has_cyano = True
                    findings_json["atomic_checks"]["functional_groups"].append("Nitrile")
                    findings_json["structural_constraints"].append({
                        "type": "positional",
                        "details": {
                            "target": "Nitrile",
                            "position": "last_stage"
                        }
                    })
                else:
                    # No need to continue if final product has no cyano
                    return

        # For reaction nodes when final product has cyano
        elif node["type"] == "reaction" and final_product_has_cyano:
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if product has cyano
            product_has_cyano = checker.check_fg("Nitrile", product)

            # Check if any reactant has cyano
            reactant_has_cyano = any(checker.check_fg("Nitrile", r) for r in reactants)

            # If product has cyano but no reactant does, cyano was created
            if product_has_cyano and not reactant_has_cyano:
                cyano_created_in_synthesis = True
                # Assuming 'Nitrile_formation' is the name for this event in the strategy
                findings_json["atomic_checks"]["named_reactions"].append("Nitrile_formation")
                findings_json["structural_constraints"].append({
                    "type": "negation",
                    "details": {
                        "target": "Nitrile_formation"
                    }
                })

        # Process children
        for child in node.get("children", []):
            # Determine the new depth based on the current node's type
            new_depth = depth
            if node["type"] != "reaction": # If current node is 'mol' (chemical)
                new_depth = depth + 1
            # If current node is 'reaction', new_depth remains 'depth'

            dfs_traverse(child, new_depth)

    # Start traversal from the root
    dfs_traverse(route)

    # Cyano is preserved if it's in the final product and wasn't created during synthesis
    preserved = final_product_has_cyano and not cyano_created_in_synthesis

    return preserved, findings_json

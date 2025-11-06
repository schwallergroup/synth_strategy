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
    This function detects a synthetic strategy where a lactam core structure
    is preserved throughout the synthesis while peripheral modifications are made.
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

    # Track if lactam core is preserved
    lactam_core_preserved = True
    lactam_core_present = False

    # Helper function to identify lactam structures
    def is_lactam(smiles):
        # A lactam is a cyclic amide. The checker library is the robust
        # and correct tool for identifying this functional group.
        return checker.check_fg("lactam", smiles)

    def dfs_traverse(node, depth=0):
        nonlocal lactam_core_preserved, lactam_core_present, findings_json

        if node["type"] == "mol":
            mol_smiles = node["smiles"]
            has_lactam = is_lactam(mol_smiles)

            # Check if the final product (depth=0) has a lactam core
            if depth == 0 and has_lactam:
                lactam_core_present = True
                findings_json["atomic_checks"]["functional_groups"].append("lactam")
                findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "lactam", "position": "last_stage"}})
                print(f"Found lactam core in final product: {mol_smiles}")
            elif has_lactam:
                findings_json["atomic_checks"]["functional_groups"].append("lactam")
                print(f"Found lactam core in intermediate: {mol_smiles}")

        elif node["type"] == "reaction":
            # For reaction nodes, check if the lactam core is preserved
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if product has lactam
                product_has_lactam = is_lactam(product)

                # If product has lactam, at least one reactant should have it for preservation
                if product_has_lactam:
                    reactant_has_lactam = any(is_lactam(r) for r in reactants)
                    if not reactant_has_lactam:
                        # No reactant has lactam core - core is not preserved
                        lactam_core_preserved = False
                        # This implies a lactam formation, which is a negation constraint
                        findings_json["structural_constraints"].append({"type": "negation", "details": {"target": "lactam_formation"}})
                        print(f"Lactam core not preserved in reaction: {rsmi}")
                    else:
                        print(f"Lactam core preserved in reaction: {rsmi}")
            except Exception as e:
                print(f"Error analyzing reaction: {e}")

        # Traverse children
        for child in node.get("children", []):
            # New depth calculation logic
            if node["type"] == "reaction":
                # If current node is a reaction, depth remains the same for children
                dfs_traverse(child, depth)
            else:
                # If current node is not a reaction (e.g., 'mol'), depth increases
                dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)

    # Strategy is present if lactam core is both present and preserved
    strategy_present = lactam_core_present and lactam_core_preserved

    if strategy_present:
        print("Detected lactam core preservation strategy")
    else:
        if not lactam_core_present:
            print("No lactam core detected in the final product")
        elif not lactam_core_preserved:
            print("Lactam core was detected but not preserved throughout the synthesis")

    return strategy_present, findings_json

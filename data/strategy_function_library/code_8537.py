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
    This function detects a strategy involving transformation of a nitro group
    during heterocycle formation.
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

    nitro_in_heterocycle_formation = False

    def dfs_traverse(node, depth=0):
        nonlocal nitro_in_heterocycle_formation, findings_json

        if node["type"] == "reaction":
            try:
                # Extract reactants and product
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check for nitro group in reactants
                has_nitro_in_reactants = False
                for r in reactants_smiles:
                    if checker.check_fg("Nitro group", r):
                        has_nitro_in_reactants = True
                        if "Nitro group" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Nitro group")
                        break

                if has_nitro_in_reactants:
                    # Check if nitro group is absent in product (transformed)
                    nitro_transformed = not checker.check_fg("Nitro group", product_smiles)

                    if nitro_transformed:
                        # Check for heterocycle formation reaction
                        is_heterocycle_formation = checker.check_reaction(
                            "Formation of NOS Heterocycles", rsmi
                        )

                        if is_heterocycle_formation:
                            if "Formation of NOS Heterocycles" not in findings_json["atomic_checks"]["named_reactions"]:
                                findings_json["atomic_checks"]["named_reactions"].append("Formation of NOS Heterocycles")
                            nitro_in_heterocycle_formation = True
            except Exception:
                # Silently ignore errors in reaction processing
                pass

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, increase depth
                new_depth = depth + 1
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    if nitro_in_heterocycle_formation:
        # Add the structural constraint if the main condition is met
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "Formation of NOS Heterocycles",
                    "Nitro group consumption"
                ],
                "scope": "reaction_step"
            }
        })

    return nitro_in_heterocycle_formation, findings_json

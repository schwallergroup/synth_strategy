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
    Detects if a trifluoromethoxy group is present in the final product and is carried through the synthesis from a precursor, rather than being formed *de novo*. It verifies that:
    1. The final target molecule contains a trifluoromethoxy group.
    2. For any reaction step that yields a product with a trifluoromethoxy group, at least one of the reactants must also contain that group.
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

    # Track if the trifluoromethoxy group is maintained
    maintained = True

    def dfs_traverse(node, depth=0):
        nonlocal maintained, findings_json

        # For molecule nodes
        if node["type"] == "mol":
            mol_smiles = node["smiles"]
            has_trifluoromethoxy = checker.check_fg("Trifluoromethoxy group", mol_smiles)

            # If this is the target molecule (depth 0), it must have a trifluoromethoxy group
            if depth == 0:
                if not has_trifluoromethoxy:
                    maintained = False
                else:
                    # If target has trifluoromethoxy, record it as an atomic check
                    if "Trifluoromethoxy group" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Trifluoromethoxy group")
                    # Record the positional constraint
                    findings_json["structural_constraints"].append({
                        "type": "positional",
                        "details": {
                            "target": "Trifluoromethoxy group",
                            "position": "last_stage"
                        }
                    })

        # For reaction nodes
        elif node["type"] == "reaction":
            # Get the reaction SMILES
            if "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check if product has trifluoromethoxy
                product_has_trifluoromethoxy = checker.check_fg("Trifluoromethoxy group", product_smiles)

                # If product has trifluoromethoxy, at least one reactant should have it
                if product_has_trifluoromethoxy:
                    if "Trifluoromethoxy group" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Trifluoromethoxy group")

                    reactant_has_trifluoromethoxy = any(
                        checker.check_fg("Trifluoromethoxy group", reactant)
                        for reactant in reactants_smiles
                    )

                    if not reactant_has_trifluoromethoxy:
                        maintained = False
                    else:
                        # If the group is maintained, it means the negation constraint is met for this step
                        findings_json["structural_constraints"].append({
                            "type": "negation",
                            "details": {
                                "target": "functional_group_formation",
                                "details": {
                                    "functional_group": "Trifluoromethoxy group"
                                }
                            }
                        })

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else: # node["type"] == "mol"
                dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    return maintained, findings_json

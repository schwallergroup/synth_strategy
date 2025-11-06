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
    This function detects the strategy of installing a methylthio group early in the synthesis
    that later serves as a leaving group in a nucleophilic substitution.
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

    # Track if we find methylthio installation and its use as leaving group
    methylthio_installed = False
    methylthio_as_leaving_group = False
    installation_depth = -1
    leaving_group_depth = -1

    def dfs_traverse(node, depth=0):
        nonlocal methylthio_installed, methylthio_as_leaving_group, installation_depth, leaving_group_depth, findings_json

        if node["type"] == "reaction":
            try:
                # Extract reactants and products
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for methylthio installation
                if not methylthio_installed:
                    # Look for reaction where methylthio group is added
                    # Check if the group is absent in reactants but present in the product
                    if not any(
                        checker.check_fg("Monosulfide", r) for r in reactants
                    ) and checker.check_fg("Monosulfide", product):
                        # Verify this is an installation by checking if the methylthio is new in the product
                        methylthio_installed = True
                        installation_depth = depth
                        findings_json["atomic_checks"]["functional_groups"].append("Monosulfide")
                        findings_json["atomic_checks"]["named_reactions"].append("Monosulfide_formation")
                        print(f"Methylthio group installed at depth {depth}, reaction: {rsmi}")

                # Check for methylthio group being used as leaving group in a nucleophilic substitution
                if any(
                    checker.check_fg("Monosulfide", r) for r in reactants
                ) and not checker.check_fg("Monosulfide", product):
                    # The methylthio group was present in reactant but not in product
                    # This indicates it was used as a leaving group
                    methylthio_as_leaving_group = True
                    leaving_group_depth = depth
                    findings_json["atomic_checks"]["functional_groups"].append("Monosulfide")
                    findings_json["atomic_checks"]["named_reactions"].append("Monosulfide_destruction")
                    print(
                        f"Methylthio group used as leaving group at depth {depth}, reaction: {rsmi}"
                    )
            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    # Strategy is present if methylthio was installed early and used as leaving group later
    # In retrosynthetic direction, early steps have higher depth
    strategy_present = (
        methylthio_installed
        and methylthio_as_leaving_group
        and installation_depth > leaving_group_depth
    )

    if methylthio_installed and methylthio_as_leaving_group:
        # Add co-occurrence constraint if both events happened
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "Monosulfide_formation",
                    "Monosulfide_destruction"
                ]
            }
        })

    if methylthio_installed and methylthio_as_leaving_group and installation_depth > leaving_group_depth:
        # Add sequence constraint if the depth condition is met
        findings_json["structural_constraints"].append({
            "type": "sequence",
            "details": {
                "before": "Monosulfide_formation",
                "after": "Monosulfide_destruction"
            }
        })

    print(f"Methylthio leaving group strategy detected: {strategy_present}")
    if strategy_present:
        print(
            f"  Installed at depth {installation_depth}, used as leaving group at depth {leaving_group_depth}"
        )

    return strategy_present, findings_json

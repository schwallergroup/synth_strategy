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
    This function detects sulfonamide cleavage to form a primary amine.
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

    sulfonamide_cleavage_found = False

    def dfs_traverse(node, depth=0):
        nonlocal sulfonamide_cleavage_found, findings_json

        if node["type"] == "reaction":
            if "mapped_reaction_smiles" in node.get("metadata", {}):
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for forward direction: sulfonamide cleavage
                has_sulfonamide_reactant = False
                for r in reactants:
                    if r and checker.check_fg("Sulfonamide", r):
                        has_sulfonamide_reactant = True
                        if "Sulfonamide" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Sulfonamide")
                        break

                has_primary_amine_product = False
                if product and checker.check_fg("Primary amine", product):
                    has_primary_amine_product = True
                    if "Primary amine" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Primary amine")

                # Check for relevant reaction types
                is_cleavage_reaction = checker.check_reaction(
                    "Hydrolysis or Hydrogenolysis of amides/imides/carbamates", rsmi
                )
                if is_cleavage_reaction:
                    if "Hydrolysis or Hydrogenolysis of amides/imides/carbamates" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("Hydrolysis or Hydrogenolysis of amides/imides/carbamates")

                # Forward direction: sulfonamide in reactant, primary amine in product
                if (
                    has_sulfonamide_reactant
                    and has_primary_amine_product
                    and is_cleavage_reaction
                ):
                    sulfonamide_cleavage_found = True
                    # Add the structural constraint if all conditions are met
                    constraint_obj = {
                        "type": "co-occurrence",
                        "details": {
                            "targets": [
                                "Hydrolysis or Hydrogenolysis of amides/imides/carbamates",
                                "Sulfonamide",
                                "Primary amine"
                            ],
                            "context": "A single reaction step must be a 'Hydrolysis or Hydrogenolysis of amides/imides/carbamates' where a 'Sulfonamide' is present in a reactant and a 'Primary amine' is present in the product."
                        }
                    }
                    if constraint_obj not in findings_json["structural_constraints"]:
                        findings_json["structural_constraints"].append(constraint_obj)

        # Determine the new depth for recursive calls
        new_depth = depth
        if node["type"] != "reaction": # This means it's a chemical node
            new_depth = depth + 1

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, new_depth)

    # Start traversal from root
    dfs_traverse(route)
    return sulfonamide_cleavage_found, findings_json

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
    """This function detects a late-stage ketone reduction strategy, occurring within the last three steps of a synthesis. It identifies reactions where a ketone is reduced to either a secondary alcohol or a methylene group by checking for specific named reactions (e.g., Wolff-Kishner, Clemmensen) or by observing the disappearance of a ketone and the appearance of a secondary alcohol."""
    findings_template = {
        "atomic_checks": {
            "named_reactions": [],
            "ring_systems": [],
            "functional_groups": []
        },
        "structural_constraints": []
    }
    findings_json = copy.deepcopy(findings_template)

    found_ketone_reduction = False

    def dfs_traverse(node, depth=0):
        nonlocal found_ketone_reduction, findings_json

        if node["type"] == "reaction" and depth <= 3:
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_part = rsmi.split(">")[0]
                product = rsmi.split(">")[-1]
                reactants = reactants_part.split(".")

                ketone_found_in_reactant = False
                for reactant_smiles in reactants:
                    if checker.check_fg("Ketone", reactant_smiles):
                        ketone_found_in_reactant = True
                        if "Ketone" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Ketone")
                        break

                if ketone_found_in_reactant:
                    # Check for specific named reduction reactions
                    named_reactions_to_check = [
                        "Reduction of ketone to secondary alcohol",
                        "Reduction of aldehydes and ketones to alcohols",
                        "Wolff-Kishner reduction",
                        "Clemmensen reduction"
                    ]
                    named_reaction_found = False
                    for r_name in named_reactions_to_check:
                        if checker.check_reaction(r_name, rsmi):
                            found_ketone_reduction = True
                            if r_name not in findings_json["atomic_checks"]["named_reactions"]:
                                findings_json["atomic_checks"]["named_reactions"].append(r_name)
                            named_reaction_found = True
                            break

                    # Generic check for reduction to an alcohol
                    if not named_reaction_found and not checker.check_fg("Ketone", product) and checker.check_fg("Secondary alcohol", product):
                        found_ketone_reduction = True
                        if "Secondary alcohol" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Secondary alcohol")

                    if found_ketone_reduction:
                        # Add the positional constraint if a ketone reduction is found within depth <= 3
                        positional_constraint = {
                            "type": "positional",
                            "details": {
                                "target": "ketone_reduction",
                                "position": "depth <= 3"
                            }
                        }
                        if positional_constraint not in findings_json["structural_constraints"]:
                            findings_json["structural_constraints"].append(positional_constraint)

        for child in node.get("children", []):
            new_depth = depth
            if node["type"] != "reaction":  # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    dfs_traverse(route)
    return found_ketone_reduction, findings_json

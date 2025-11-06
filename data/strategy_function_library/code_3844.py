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
    """Detects the late-stage introduction of an aminoalkoxy group via ether synthesis. This is identified by checking the final synthetic step (depth=1) for the formation of a new ether bond, where one reactant contains both a tertiary amine and an alcohol or halide, and the tertiary amine is present in the product."""
    findings_template = {
        "atomic_checks": {
            "named_reactions": [],
            "ring_systems": [],
            "functional_groups": []
        },
        "structural_constraints": []
    }
    findings_json = copy.deepcopy(findings_template)

    final_step_is_alkylation = False

    def dfs_traverse(node, depth=0):
        nonlocal final_step_is_alkylation, findings_json

        # Check for alkylation in the final step
        if node["type"] == "reaction" and depth == 1:  # Final step
            if "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                print(f"Examining reaction at depth {depth}: {rsmi}")

                # Check for aminoalkoxy reagent (like dimethylaminoethanol)
                has_aminoalkoxy_reagent = False
                aminoalkoxy_reactant = None

                for reactant in reactants:
                    # Check for tertiary amine with alcohol or halide (alkylating agent)
                    fg_checks = []
                    if checker.check_fg("Tertiary amine", reactant):
                        fg_checks.append("Tertiary amine")
                    
                    alcohol_or_halide_found = False
                    if checker.check_fg("Primary alcohol", reactant):
                        fg_checks.append("Primary alcohol")
                        alcohol_or_halide_found = True
                    if checker.check_fg("Secondary alcohol", reactant):
                        fg_checks.append("Secondary alcohol")
                        alcohol_or_halide_found = True
                    if checker.check_fg("Primary halide", reactant):
                        fg_checks.append("Primary halide")
                        alcohol_or_halide_found = True
                    if checker.check_fg("Secondary halide", reactant):
                        fg_checks.append("Secondary halide")
                        alcohol_or_halide_found = True

                    if "Tertiary amine" in fg_checks and alcohol_or_halide_found:
                        has_aminoalkoxy_reagent = True
                        aminoalkoxy_reactant = reactant
                        for fg in fg_checks:
                            if fg not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append(fg)
                        print(f"Found potential aminoalkyl reagent: {reactant}")
                        break

                # Verify the amine group is transferred to the product and an ether is formed
                if has_aminoalkoxy_reagent and checker.check_fg("Tertiary amine", product):
                    if "Tertiary amine" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Tertiary amine")

                    # Check if an ether was formed in the product
                    if checker.check_fg("Ether", product):
                        if "Ether" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Ether")

                        # Check if the ether is newly formed (not present in all reactants)
                        ether_in_all_reactants = True
                        for r in reactants:
                            if not checker.check_fg("Ether", r):
                                ether_in_all_reactants = False
                                break

                        if not ether_in_all_reactants:
                            print(
                                f"Confirmed ether formation with aminoalkyl group at depth {depth}"
                            )
                            final_step_is_alkylation = True
                            # Record structural constraints if all conditions are met
                            findings_json["structural_constraints"].append({
                                "type": "positional",
                                "details": {
                                    "target": "aminoalkoxy_ether_synthesis",
                                    "position": "last_stage"
                                }
                            })
                            findings_json["structural_constraints"].append({
                                "type": "co-occurrence",
                                "details": {
                                    "targets": [
                                        "Ether",
                                        "Tertiary amine"
                                    ]
                                }
                            })
                            findings_json["structural_constraints"].append({
                                "type": "negation",
                                "details": {
                                    "target": "Ether group is present in all reactants"
                                }
                            })

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

    print(f"Late-stage alkylation strategy detected: {final_step_is_alkylation}")
    return final_step_is_alkylation, findings_json

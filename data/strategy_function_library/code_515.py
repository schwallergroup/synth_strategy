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


AROMATIC_NITRATION_REACTIONS = [
    "Aromatic nitration with HNO3",
    "Aromatic nitration with NO3 salt",
    "Aromatic nitration with NO2 salt",
    "Aromatic nitration with alkyl NO2",
]

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects a synthetic strategy where an aromatic nitro group is introduced
    and subsequently reduced to an amine in a later step of the synthesis. The specific
    aromatic nitration reactions checked are defined in the AROMATIC_NITRATION_REACTIONS list.
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

    # Track if we found nitration and reduction steps
    nitration_found = False
    nitro_reduction_found = False

    # Track the depth of these reactions
    nitration_depth = -1
    reduction_depth = -1

    print("Starting analysis of synthesis route")

    def dfs_traverse(node, depth=0):
        nonlocal nitration_found, nitro_reduction_found, nitration_depth, reduction_depth, findings_json

        if node["type"] == "reaction":
            try:
                # Extract reactants and products
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                print(f"Analyzing reaction at depth {depth}: {rsmi}")

                # Check for nitration reaction
                is_nitration = False
                for nitration_type in AROMATIC_NITRATION_REACTIONS:
                    if checker.check_reaction(nitration_type, rsmi):
                        is_nitration = True
                        findings_json["atomic_checks"]["named_reactions"].append(nitration_type)
                        print(f"Found nitration reaction: {nitration_type}")
                        break

                if is_nitration:
                    nitration_found = True
                    nitration_depth = depth
                    print(f"Confirmed nitration reaction at depth {depth}")

                # Check for nitro reduction (nitro to amine)
                is_nitro_reduction = checker.check_reaction(
                    "Reduction of nitro groups to amines", rsmi
                )
                if is_nitro_reduction:
                    findings_json["atomic_checks"]["named_reactions"].append("Reduction of nitro groups to amines")

                # If not a standard reduction, check for functional group changes
                if not is_nitro_reduction:
                    reactants_have_nitro = False
                    for r in reactants_smiles:
                        if checker.check_fg("Nitro group", r):
                            reactants_have_nitro = True
                            findings_json["atomic_checks"]["functional_groups"].append("Nitro group")
                            break

                    product_has_nitro = checker.check_fg("Nitro group", product_smiles)
                    if product_has_nitro:
                        findings_json["atomic_checks"]["functional_groups"].append("Nitro group")

                    product_has_amine = False
                    if checker.check_fg("Primary amine", product_smiles):
                        product_has_amine = True
                        findings_json["atomic_checks"]["functional_groups"].append("Primary amine")
                    if checker.check_fg("Aniline", product_smiles):
                        product_has_amine = True
                        findings_json["atomic_checks"]["functional_groups"].append("Aniline")

                    if reactants_have_nitro and not product_has_nitro and product_has_amine:
                        is_nitro_reduction = True
                        print("Detected nitro reduction based on functional group changes")

                if is_nitro_reduction:
                    nitro_reduction_found = True
                    reduction_depth = depth
                    print(f"Confirmed nitro reduction at depth {depth}")

            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Determine the next depth based on the current node's type
        next_depth = depth
        if node["type"] != "reaction": # This means it's a 'chemical' node
            next_depth = depth + 1

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, next_depth)

    # Start traversal from the root
    dfs_traverse(route)

    print(
        f"Analysis complete: nitration_found={nitration_found} at depth {nitration_depth}, "
        f"nitro_reduction_found={nitro_reduction_found} at depth {reduction_depth}"
    )

    result = False
    # Check if we found both nitration and reduction in the correct order
    # Remember: lower depth is later in the synthesis (closer to final product)
    if nitration_found and nitro_reduction_found and reduction_depth < nitration_depth:
        print("Late-stage nitro reduction strategy detected")
        result = True
        # Add the structural constraint if the condition is met
        findings_json["structural_constraints"].append({
            "type": "sequence",
            "details": {
                "before": "aromatic_nitration",
                "after": "Reduction of nitro groups to amines"
            }
        })

    return result, findings_json

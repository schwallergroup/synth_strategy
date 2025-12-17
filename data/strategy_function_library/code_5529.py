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


AMIDE_FORMATION_REACTION_TYPES = [
    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_OS",
    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
    "Acyl chloride with ammonia to amide",
    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
    "Acyl chloride with primary amine to imide",
    "Acyl chloride with secondary amine to amide",
    "Carboxylic acid with primary amine to amide",
    "Ester with ammonia to amide",
    "Ester with primary amine to amide",
    "Ester with secondary amine to amide",
    "Schotten-Baumann to ester",
    "Schotten-Baumann_amide",
    "Acylation of primary amines",
    "Acylation of secondary amines",
    "Nitrile to amide",
    "Carboxylic acid to amide conversion",
]

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects if at least two amide bond formations occur anywhere in the synthesis.
    An amide formation is identified either by matching a known reaction name from the predefined
    list (AMIDE_FORMATION_REACTION_TYPES) or by observing the net formation of a primary,
    secondary, or tertiary amide functional group in a reaction step.
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

    # Track amide formations with their depths
    amide_formation_depths = []

    def dfs_traverse(node, depth=0):
        nonlocal amide_formation_depths, findings_json
        # Process reaction nodes
        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]

            # Check if this is an amide formation reaction
            is_amide_formation = False

            # Check for specific amide formation reaction types
            for reaction_type in AMIDE_FORMATION_REACTION_TYPES:
                if checker.check_reaction(reaction_type, rsmi):
                    is_amide_formation = True
                    findings_json["atomic_checks"]["named_reactions"].append(reaction_type)
                    print(f"Found amide formation reaction: {reaction_type} at depth {depth}")
                    break

            # If not found by reaction type, check for amide formation by functional group change
            if not is_amide_formation:
                try:
                    reactants_smiles = rsmi.split(">")[0].split(".")
                    product_smiles = rsmi.split(">")[-1]

                    # Check if product has an amide group
                    product_has_amide = False
                    for fg_name in ["Primary amide", "Secondary amide", "Tertiary amide"]:
                        if checker.check_fg(fg_name, product_smiles):
                            product_has_amide = True
                            findings_json["atomic_checks"]["functional_groups"].append(fg_name)
                            break

                    # Check if reactants have amide groups
                    reactants_have_amide = False
                    for reactant_smiles in reactants_smiles:
                        for fg_name in ["Primary amide", "Secondary amide", "Tertiary amide"]:
                            if checker.check_fg(fg_name, reactant_smiles):
                                reactants_have_amide = True
                                # Do not add to findings_json here, as it's about *formation*, not presence in reactants
                                break
                        if reactants_have_amide: # Optimization: if one reactant has it, no need to check others
                            break

                    # If product has amide but reactants don't, it's an amide formation
                    if product_has_amide and not reactants_have_amide:
                        is_amide_formation = True
                        print(f"Found amide formation by FG change at depth {depth}")
                except Exception as e:
                    print(f"Error analyzing amide formation: {e}")

            if is_amide_formation:
                amide_formation_depths.append(depth)

        # Continue traversing
        for child in node.get("children", []):
            # Determine the new depth based on the current node's type
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    # Sort depths to check for sequential formations
    amide_formation_depths.sort()
    print(
        f"Found {len(amide_formation_depths)} amide formations at depths: {amide_formation_depths}"
    )

    result = len(amide_formation_depths) >= 2

    # Add structural constraint if the condition is met
    if result:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "amide_formation",
                "operator": ">=",
                "value": 2
            }
        })

    return result, findings_json

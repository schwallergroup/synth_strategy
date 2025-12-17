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


AZIDE_FORMATION_REACTIONS = [
    "Formation of Azides from halogens",
    "Formation of Azides from boronic acids",
    "Amine to azide",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Checks for the late-stage introduction of an azide group. Azide formation is identified by specific named reactions (see `AZIDE_FORMATION_REACTIONS` list) or by the de novo appearance of the azide functional group in a reaction product. The introduction is considered 'late-stage' if it occurs within the final half of the total synthetic steps.
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

    azide_introduction_depth = None
    max_depth = 0
    result = False

    def dfs_traverse(node, current_depth=0):
        nonlocal azide_introduction_depth, max_depth, findings_json

        # Track maximum depth
        max_depth = max(max_depth, current_depth)

        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for azide introduction
            has_azide_reactants = False
            for reactant in reactants:
                if checker.check_fg("Azide", reactant):
                    has_azide_reactants = True
                    if "Azide" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Azide")
                    break

            has_azide_product = checker.check_fg("Azide", product)
            if has_azide_product and "Azide" not in findings_json["atomic_checks"]["functional_groups"]:
                findings_json["atomic_checks"]["functional_groups"].append("Azide")

            # Check if this is an azide introduction reaction
            is_azide_introduction = False
            for reaction_name in AZIDE_FORMATION_REACTIONS:
                if checker.check_reaction(reaction_name, rsmi):
                    is_azide_introduction = True
                    if reaction_name not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append(reaction_name)
                    break

            # Alternative check: if reaction type check fails, verify azide appears in product but not reactants
            if not is_azide_introduction and not has_azide_reactants and has_azide_product:
                is_azide_introduction = True

            if is_azide_introduction:
                # Track the latest (lowest depth number) azide introduction
                if azide_introduction_depth is None or current_depth < azide_introduction_depth:
                    azide_introduction_depth = current_depth

        # Continue traversing
        for child in node.get("children", []):
            # New logic for depth calculation:
            # Depth increases only when going from chemical to reaction.
            # Depth remains the same when going from reaction to chemical.
            if node["type"] == "reaction":
                dfs_traverse(child, current_depth)
            else: # node["type"] == "chemical"
                dfs_traverse(child, current_depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if azide was introduced in the second half of the synthesis (late-stage)
    # A smaller depth value means later in the synthesis.
    if azide_introduction_depth is not None and azide_introduction_depth <= max_depth / 2:
        result = True
        # Add the structural constraint if the condition is met
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "azide_introduction",
                "position": "late_stage_half"
            }
        })

    return result, findings_json

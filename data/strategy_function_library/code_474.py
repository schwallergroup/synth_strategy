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


SUZUKI_REACTION_TYPES = [
    "Suzuki coupling with boronic acids",
    "Suzuki coupling with boronic acids OTf",
    "Suzuki coupling with boronic esters",
    "Suzuki coupling with boronic esters OTf",
    "Suzuki coupling with sulfonic esters",
    "Suzuki",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a late-stage (depth <= 3) Suzuki coupling that forms a biaryl system
    by joining a fragment containing a pyrazole ring with another fragment
    containing a trifluoromethoxy group. The specific reaction names are defined
    in the SUZUKI_REACTION_TYPES list.
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

    suzuki_reaction_node = None
    suzuki_depth = float("inf")
    result = False

    def dfs_traverse(node, depth=0):
        nonlocal suzuki_reaction_node, suzuki_depth, findings_json

        if node.get("type") == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            is_suzuki = False
            for r_type in SUZUKI_REACTION_TYPES:
                if checker.check_reaction(r_type, rsmi):
                    is_suzuki = True
                    if r_type not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append(r_type)
                    break

            if is_suzuki:
                # Find the Suzuki coupling that occurs at the latest stage (minimum depth)
                if depth < suzuki_depth:
                    suzuki_reaction_node = node
                    suzuki_depth = depth

        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node.get("type") != "reaction": # If current node is chemical, increase depth
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    dfs_traverse(route)

    # If no Suzuki coupling was found, or if it's not a late-stage reaction, exit.
    if not suzuki_reaction_node or suzuki_depth > 3:
        return False, findings_json

    # Add positional constraint if met
    if suzuki_depth <= 3:
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "Suzuki coupling",
                "constraint": "must occur at depth <= 3"
            }
        })

    # Verify that the identified late-stage Suzuki joins the correct fragments.
    rsmi = suzuki_reaction_node["metadata"]["mapped_reaction_smiles"]
    reactants_part, _, product_part = rsmi.partition('>>')
    reactants = reactants_part.split(".")

    reactant_with_pyrazole = None
    reactant_with_trifluoromethoxy = None

    for r in reactants:
        if checker.check_ring("pyrazole", r):
            reactant_with_pyrazole = r
            if "pyrazole" not in findings_json["atomic_checks"]["ring_systems"]:
                findings_json["atomic_checks"]["ring_systems"].append("pyrazole")
        if checker.check_fg("Trifluoromethoxy group", r):
            reactant_with_trifluoromethoxy = r
            if "Trifluoromethoxy group" not in findings_json["atomic_checks"]["functional_groups"]:
                findings_json["atomic_checks"]["functional_groups"].append("Trifluoromethoxy group")

    # Check if the product has both groups, confirming they were joined.
    product_has_both_pyrazole = checker.check_ring("pyrazole", product_part)
    if product_has_both_pyrazole and "pyrazole" not in findings_json["atomic_checks"]["ring_systems"]:
        findings_json["atomic_checks"]["ring_systems"].append("pyrazole")

    product_has_both_trifluoromethoxy = checker.check_fg(
        "Trifluoromethoxy group", product_part
    )
    if product_has_both_trifluoromethoxy and "Trifluoromethoxy group" not in findings_json["atomic_checks"]["functional_groups"]:
        findings_json["atomic_checks"]["functional_groups"].append("Trifluoromethoxy group")

    product_has_both = product_has_both_pyrazole and product_has_both_trifluoromethoxy

    # The strategy is present if two different reactants, one with pyrazole and one
    # with trifluoromethoxy, are joined to form a product containing both.
    result = (
        reactant_with_pyrazole is not None
        and reactant_with_trifluoromethoxy is not None
        and reactant_with_pyrazole != reactant_with_trifluoromethoxy
        and product_has_both
    )

    if result:
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "scope": "single_reaction",
                "reaction": "Suzuki coupling",
                "constraint": "A reactant must contain 'pyrazole' and a different reactant must contain 'Trifluoromethoxy group', and the product must contain both."
            }
        })

    return result, findings_json

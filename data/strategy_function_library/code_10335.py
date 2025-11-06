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
    Detects a synthetic strategy where a trifluoromethyl-containing fragment is incorporated in the final synthetic step (late stage) via acylation of a hydrazide or acylhydrazide, forming a trifluoromethyl-acylhydrazide.
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

    # Initialize tracking variables
    has_trifluoromethyl_final = False
    has_acylhydrazide_final = False
    has_acylation_reaction = False

    def dfs_traverse(node, depth=0):
        nonlocal has_trifluoromethyl_final, has_acylhydrazide_final, has_acylation_reaction, findings_json

        if node["type"] == "mol" and depth == 0:  # Final product
            # Check for trifluoromethyl group and acylhydrazide in final product
            if checker.check_fg("Trifluoro group", node["smiles"]):
                has_trifluoromethyl_final = True
                findings_json["atomic_checks"]["functional_groups"].append("Trifluoro group")
                findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "Trifluoro group", "position": "last_stage"}})
                print(f"Found trifluoromethyl group in final product: {node['smiles']}")

            if checker.check_fg("Acylhydrazine", node["smiles"]):
                has_acylhydrazide_final = True
                findings_json["atomic_checks"]["functional_groups"].append("Acylhydrazine")
                findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "Acylhydrazine", "position": "last_stage"}})
                print(f"Found acylhydrazide in final product: {node['smiles']}")

        elif node["type"] == "reaction":
            if "mapped_reaction_smiles" in node["metadata"]:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for acylation reaction introducing CF3 to a hydrazide in the final step
                if depth == 1 and checker.check_fg("Trifluoro group", product) and checker.check_fg(
                    "Acylhydrazine", product
                ):
                    has_cf3_reactant = False
                    has_hydrazide_reactant = False

                    for reactant in reactants:
                        if checker.check_fg("Trifluoro group", reactant):
                            has_cf3_reactant = True
                            if "Trifluoro group" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Trifluoro group")
                        if checker.check_fg("Acylhydrazine", reactant):
                            has_hydrazide_reactant = True
                            if "Acylhydrazine" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Acylhydrazine")
                        elif checker.check_fg("Hydrazine", reactant):
                            has_hydrazide_reactant = True
                            if "Hydrazine" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Hydrazine")

                    # Check if one reactant has CF3 and another has hydrazide/hydrazine
                    if has_cf3_reactant and has_hydrazide_reactant:
                        has_acylation_reaction = True
                        if "acylation" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("acylation")
                        findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "acylation", "position": "last_stage"}})
                        print(f"Found acylation reaction introducing CF3 group: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            # New depth calculation logic
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    # Check if the strategy is present
    strategy_present = (
        has_trifluoromethyl_final
        and has_acylhydrazide_final
        and has_acylation_reaction
    )

    if strategy_present:
        # Add the co-occurrence constraint if all conditions are met
        findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"targets": ["Trifluoro group", "Acylhydrazine", "acylation"]}})

    print(
        f"Strategy components found: CF3 in final: {has_trifluoromethyl_final}, "
        + f"Acylhydrazide in final: {has_acylhydrazide_final}, "
        + f"Acylation reaction: {has_acylation_reaction}"
    )

    return strategy_present, findings_json

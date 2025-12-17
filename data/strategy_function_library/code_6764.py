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


def main(route) -> Tuple[bool, Dict]:
    """
    Detects if the synthesis involves late-stage formation of a hydrazide from an ester and hydrazine, typically occurring in the final two steps (depth <= 1).
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

    late_stage_hydrazide = False

    def dfs_traverse(node, depth=0):
        nonlocal late_stage_hydrazide, findings_json

        if node["type"] == "reaction" and depth <= 1:
            try:
                if "mapped_reaction_smiles" in node.get("metadata", {}):
                    rsmi = node["metadata"]["mapped_reaction_smiles"]
                    reactants = rsmi.split(">")[0].split(".")
                    product = rsmi.split(">")[-1]

                    # Check if any reactant contains an ester group
                    has_ester = False
                    for r in reactants:
                        if r and checker.check_fg("Ester", r):
                            has_ester = True
                            if "Ester" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Ester")
                            break

                    # Check if any reactant contains hydrazine
                    has_hydrazine = False
                    for r in reactants:
                        if r and checker.check_fg("Hydrazine", r):
                            has_hydrazine = True
                            if "Hydrazine" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Hydrazine")
                            break

                    # Check if product contains hydrazide
                    has_hydrazide = checker.check_fg("Acylhydrazine", product)
                    if has_hydrazide:
                        if "Acylhydrazine" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Acylhydrazine")

                    # Verify that this is a hydrazide formation based on functional group transformation
                    if has_ester and has_hydrazine and has_hydrazide:
                        late_stage_hydrazide = True
                        if "hydrazide_formation_from_ester" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("hydrazide_formation_from_ester")

            except Exception:
                # Silently ignore errors in reaction analysis
                pass

        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction":  # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same

            dfs_traverse(child, new_depth)

    dfs_traverse(route)

    if late_stage_hydrazide:
        # Add the structural constraint if the condition is met
        findings_json["structural_constraints"].append(
            {
                "type": "positional",
                "details": {
                    "target": "hydrazide_formation_from_ester",
                    "position": "last_two_stages"
                }
            }
        )

    return late_stage_hydrazide, findings_json

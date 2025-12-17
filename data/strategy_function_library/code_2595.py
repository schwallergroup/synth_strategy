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


HYDRAZIDE_DERIVATIVES = ["Acylhydrazine", "Hydrazone amide"]
ACYL_SOURCES = ["Acyl halide", "Ester", "Carboxylic acid", "Anhydride"]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects the late-stage formation of specific hydrazide derivatives (Acylhydrazine, Hydrazone amide) in the final or penultimate step. The check verifies that one of these groups is formed in the product, was not present in the reactants, and that the reactants include a hydrazine and a suitable acyl source (Acyl halide, Ester, Carboxylic acid, or Anhydride).
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

    has_late_stage_hydrazide = False

    def dfs_traverse(node, depth=0):
        nonlocal has_late_stage_hydrazide, findings_json

        if has_late_stage_hydrazide: # Optimization: stop searching if already found
            return

        if node["type"] == "reaction":
            # Check if this is a late-stage step (depth 1 or 2 in some conventions, here 0 or 1)
            if depth <= 1:
                try:
                    if "rsmi" not in node.get("metadata", {}):
                        return

                    rsmi = node["metadata"]["mapped_reaction_smiles"]
                    reactants_part = rsmi.split(">")[0]
                    reactants = reactants_part.split(".")
                    product = rsmi.split(">")[-1]

                    # Check if product contains a target hydrazide derivative
                    product_has_hydrazide = False
                    for fg in HYDRAZIDE_DERIVATIVES:
                        if checker.check_fg(fg, product):
                            product_has_hydrazide = True
                            if fg not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append(fg)

                    # Check if any reactant contains a target hydrazide derivative
                    reactants_have_hydrazide = False
                    for reactant in reactants:
                        for fg in HYDRAZIDE_DERIVATIVES:
                            if checker.check_fg(fg, reactant):
                                reactants_have_hydrazide = True
                                # This is a negation check, so we don't add it to findings_json if found in reactants
                                break
                        if reactants_have_hydrazide: break

                    # Check for hydrazine in reactants
                    has_hydrazine = False
                    for reactant in reactants:
                        if checker.check_fg("Hydrazine", reactant):
                            has_hydrazine = True
                            if "Hydrazine" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Hydrazine")
                            break

                    # Check for a suitable acyl source in reactants
                    has_acyl_source = False
                    for reactant in reactants:
                        for fg in ACYL_SOURCES:
                            if checker.check_fg(fg, reactant):
                                has_acyl_source = True
                                if fg not in findings_json["atomic_checks"]["functional_groups"]:
                                    findings_json["atomic_checks"]["functional_groups"].append(fg)
                                break
                        if has_acyl_source: break

                    # Verify hydrazide formation: product has a target group, reactants don't,
                    # and appropriate reagents (hydrazine + acyl source) are present.
                    if (
                        product_has_hydrazide
                        and not reactants_have_hydrazide
                        and has_hydrazine
                        and has_acyl_source
                    ):
                        has_late_stage_hydrazide = True
                        # Add named reaction if detected
                        if "hydrazide_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("hydrazide_formation")

                        # Add structural constraints
                        # Positional constraint
                        findings_json["structural_constraints"].append({
                            "type": "positional",
                            "details": {
                                "target": "hydrazide_formation",
                                "position": "last_or_penultimate_stage"
                            }
                        })
                        # Co-occurrence constraint
                        findings_json["structural_constraints"].append({
                            "type": "co-occurrence",
                            "details": {
                                "targets": [
                                    "Hydrazine",
                                    [
                                        "Acyl halide",
                                        "Ester",
                                        "Carboxylic acid",
                                        "Anhydride"
                                    ]
                                ],
                                "scope": "reactants"
                            }
                        })
                        # Negation constraint
                        findings_json["structural_constraints"].append({
                            "type": "negation",
                            "details": {
                                "target": [
                                    "Acylhydrazine",
                                    "Hydrazone amide"
                                ],
                                "scope": "reactants"
                            }
                        })

                except Exception:
                    # Silently ignore errors in reaction processing
                    pass

        # Traverse children
        for child in node.get("children", []):
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    return has_late_stage_hydrazide, findings_json

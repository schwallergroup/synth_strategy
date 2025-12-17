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
    Detects if an amide formation occurs in the late stage of the synthesis.
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

    late_amide_formation = False

    def dfs_traverse(node, depth=0):
        nonlocal late_amide_formation, findings_json

        if node["type"] == "reaction" and depth <= 3:  # Late stage (depth ≤ 3)
            print(f"Examining reaction at depth {depth}")
            # Check if this is an amide formation
            if "mapped_reaction_smiles" in node["metadata"]:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                print(f"Reaction SMILES: {rsmi}")

                # Check for reactants and products
                try:
                    reactants = rsmi.split(">")[0].split(".")
                    product = rsmi.split(">")[-1]

                    print(f"Reactants: {reactants}")
                    print(f"Product: {product}")

                    # Check for amine in reactants
                    amine_reactants = []
                    for r in reactants:
                        if checker.check_fg("Primary amine", r):
                            amine_reactants.append(r)
                            if "Primary amine" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Primary amine")
                            print(f"Found amine in reactant: {r}")
                        if checker.check_fg("Secondary amine", r):
                            amine_reactants.append(r)
                            if "Secondary amine" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Secondary amine")
                            print(f"Found amine in reactant: {r}")

                    # Check for acyl source in reactants
                    acyl_reactants = []
                    for r in reactants:
                        if checker.check_fg("Acyl halide", r):
                            acyl_reactants.append(r)
                            if "Acyl halide" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Acyl halide")
                            print(f"Found acyl source in reactant: {r}")
                        if checker.check_fg("Carboxylic acid", r):
                            acyl_reactants.append(r)
                            if "Carboxylic acid" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Carboxylic acid")
                            print(f"Found acyl source in reactant: {r}")
                        if checker.check_fg("Ester", r):
                            acyl_reactants.append(r)
                            if "Ester" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Ester")
                            print(f"Found acyl source in reactant: {r}")
                        if checker.check_fg("Anhydride", r):
                            acyl_reactants.append(r)
                            if "Anhydride" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Anhydride")
                            print(f"Found acyl source in reactant: {r}")

                    # Check for amide in product
                    amide_types = ["Primary amide", "Secondary amide", "Tertiary amide"]
                    product_amides = []
                    for amide_type in amide_types:
                        if checker.check_fg(amide_type, product):
                            product_amides.append(amide_type)
                            if amide_type not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append(amide_type)
                            print(f"Found {amide_type} in product")

                    # Count amides in reactants
                    reactant_amides = []
                    for r in reactants:
                        for amide_type in amide_types:
                            if checker.check_fg(amide_type, r):
                                reactant_amides.append((r, amide_type))
                                # Do not add to findings_json if it's a reactant amide, only product/newly formed
                                print(f"Found {amide_type} in reactant: {r}")

                    # Verify that a new amide is formed
                    if (
                        len(amine_reactants) > 0
                        and len(acyl_reactants) > 0
                        and len(product_amides) > 0
                        and len(reactant_amides) < len(product_amides)
                    ):
                        print("Found amide formation based on functional group analysis")
                        late_amide_formation = True
                        if "amide_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("amide_formation")
                        
                        # Add structural constraint if late stage and amide formation detected
                        if depth <= 3:
                            constraint_obj = {
                                "type": "positional",
                                "details": {
                                    "target": "amide_formation",
                                    "position": "late_stage"
                                }
                            }
                            if constraint_obj not in findings_json["structural_constraints"]:
                                findings_json["structural_constraints"].append(constraint_obj)

                except Exception as e:
                    print(f"Error analyzing reaction: {e}")

        # Traverse children
        for child in node.get("children", []):
            # New depth calculation logic
            if node["type"] == "reaction":
                # If current node is a reaction, depth remains the same for child (chemical)
                dfs_traverse(child, depth)
            else:
                # If current node is a chemical, depth increases for child (reaction)
                dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(f"Late-stage amide formation: {late_amide_formation}")
    return late_amide_formation, findings_json

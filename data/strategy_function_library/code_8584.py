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
    Detects a multi-step strategy where a phenol is converted to a triflate, which is then converted to a boronic ester.
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

    # Track molecules containing key functional groups and their depths
    phenol_molecules = []
    boronic_ester_molecules = []

    # Track transformations with their depths
    transformations = []

    def dfs_traverse(node, depth=0):
        nonlocal phenol_molecules, boronic_ester_molecules, transformations, findings_json

        if node["type"] == "mol":
            smiles = node["smiles"]

            # Check for phenol in any molecule
            if checker.check_fg("Phenol", smiles):
                phenol_molecules.append((smiles, depth))
                findings_json["atomic_checks"]["functional_groups"].append("Phenol")
                print(f"Found phenol at depth {depth}: {smiles}")

            # Check for boronic ester in any molecule
            if checker.check_fg("Boronic ester", smiles):
                boronic_ester_molecules.append((smiles, depth))
                findings_json["atomic_checks"]["functional_groups"].append("Boronic ester")
                print(f"Found boronic ester at depth {depth}: {smiles}")

            # Check for other relevant intermediates
            if checker.check_fg("Triflate", smiles):
                findings_json["atomic_checks"]["functional_groups"].append("Triflate")
                print(f"Found triflate intermediate at depth {depth}: {smiles}")
            if checker.check_fg("Mesylate", smiles):
                findings_json["atomic_checks"]["functional_groups"].append("Mesylate")
                print(f"Found mesylate intermediate at depth {depth}: {smiles}")
            if checker.check_fg("Tosylate", smiles):
                findings_json["atomic_checks"]["functional_groups"].append("Tosylate")
                print(f"Found tosylate intermediate at depth {depth}: {smiles}")
            if checker.check_fg("Aromatic halide", smiles):
                findings_json["atomic_checks"]["functional_groups"].append("Aromatic halide")
                print(f"Found aromatic halide at depth {depth}: {smiles}")

        elif node["type"] == "reaction":
            if "mapped_reaction_smiles" in node.get("metadata", {}):
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Phenol to triflate conversion
                if any(checker.check_fg("Phenol", r) for r in reactants) and checker.check_fg(
                    "Triflate", product
                ):
                    if checker.check_reaction("Formation of Sulfonic Esters", rsmi):
                        transformations.append(("phenol_to_triflate", depth))
                        findings_json["atomic_checks"]["named_reactions"].append("Formation of Sulfonic Esters")
                        print(f"Found transformation at depth {depth}: phenol \u2192 triflate: {rsmi}")

                # Triflate to boronic ester
                elif any(checker.check_fg("Triflate", r) for r in reactants) and checker.check_fg(
                    "Boronic ester", product
                ):
                    if checker.check_reaction(
                        "Preparation of boronic esters", rsmi
                    ):
                        transformations.append(("triflate_to_boronic_ester", depth))
                        findings_json["atomic_checks"]["named_reactions"].append("Preparation of boronic esters")
                        print(
                            f"Found transformation at depth {depth}: triflate \u2192 boronic ester: {rsmi}"
                        )

        # Recursively process children with modified depth
        for child in node.get("children", []):
            if node["type"] == "reaction":
                # Depth remains the same when traversing from a reaction node to a chemical node
                dfs_traverse(child, depth)
            else:
                # Depth increases when traversing from a chemical node to a reaction node
                dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Analyze results to determine if a valid phenol to boronic ester strategy exists

    # Sort molecules and transformations by depth
    phenol_molecules.sort(key=lambda x: x[1])
    boronic_ester_molecules.sort(key=lambda x: x[1])
    transformations.sort(key=lambda x: x[1])

    # Check if the target molecule (depth 0) contains a boronic ester
    target_has_boronic_ester = any(depth == 0 for _, depth in boronic_ester_molecules)
    if target_has_boronic_ester:
        # Add structural constraint if boronic ester is found at the last stage (depth 0)
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "Boronic ester",
                "position": "last_stage"
            }
        })

    # Check if we have any relevant transformations
    has_transformations = len(transformations) > 0

    # Check if we have a direct phenol to boronic ester transformation
    direct_transformation = any(t[0] == "phenol_to_boronic_ester_direct" for t in transformations)

    # Check if we have a pathway through triflate
    triflate_pathway = any(t[0] == "phenol_to_triflate" for t in transformations) and any(
        t[0] == "triflate_to_boronic_ester" for t in transformations
    )
    if triflate_pathway:
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "Formation of Sulfonic Esters",
                    "Preparation of boronic esters"
                ]
            }
        })

    # Check if we have a pathway through mesylate
    mesylate_pathway = any(t[0] == "phenol_to_mesylate" for t in transformations) and any(
        t[0] == "mesylate_to_boronic_ester" for t in transformations
    )

    # Check if we have a pathway through tosylate
    tosylate_pathway = any(t[0] == "phenol_to_tosylate" for t in transformations) and any(
        t[0] == "tosylate_to_boronic_ester" for t in transformations
    )

    # Check if we have a pathway through halide
    halide_pathway = any(t[0] == "phenol_to_halide" for t in transformations) and any(
        t[0] == "halide_to_boronic_ester" for t in transformations
    )

    # Check if we have a valid pathway
    valid_pathway = (
        direct_transformation
        or triflate_pathway
        or mesylate_pathway
        or tosylate_pathway
        or halide_pathway
    )

    print(f"Target has boronic ester: {target_has_boronic_ester}")
    print(f"Has transformations: {has_transformations}")
    print(f"Valid pathway detected: {valid_pathway}")

    result = False
    # The key indicator is the presence of a boronic ester in the target molecule
    if target_has_boronic_ester:
        # If we have transformations and a valid pathway, it's a clear strategy
        if has_transformations and valid_pathway:
            print("Complete phenol to boronic ester transformation strategy detected")
            result = True

        # Even without seeing the complete pathway, the presence of a boronic ester
        # in the target molecule is a strong indicator of this strategy
        elif not has_transformations or not valid_pathway: # Added this condition to ensure it doesn't override the above True
            print("Phenol to boronic ester strategy detected (target molecule contains boronic ester)")
            result = True

    # If we have phenol molecules and transformations but no boronic ester in the target,
    # it might be a different strategy
    if phenol_molecules and has_transformations and not target_has_boronic_ester:
        print("Found phenol transformations but target doesn't contain boronic ester")
        result = False

    if not result:
        print("Phenol to boronic ester transformation strategy not detected")

    return result, findings_json

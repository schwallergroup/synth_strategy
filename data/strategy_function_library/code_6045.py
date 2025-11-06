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
    This function detects a strategy involving thionation (O=C to S=C)
    followed by ring opening.
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

    # Track reaction sequence
    reaction_sequence = []
    result = False

    def dfs_traverse(node, depth=0):
        nonlocal result, findings_json
        if node["type"] == "reaction":
            try:
                # Extract reactants and product
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check if this is a thionation reaction
                is_thionation = False
                for r in reactants_smiles:
                    if (
                        (
                            checker.check_fg("Ketone", r)
                            and checker.check_fg("Thiocarbonyl", product_smiles)
                        )
                        or (
                            checker.check_fg("Aldehyde", r)
                            and checker.check_fg("Thiocarbonyl", product_smiles)
                        )
                        or (
                            checker.check_fg("Primary amide", r)
                            and checker.check_fg("Thioamide", product_smiles)
                        )
                        or (
                            checker.check_fg("Secondary amide", r)
                            and checker.check_fg("Thioamide", product_smiles)
                        )
                    ):
                        is_thionation = True
                        # Record functional groups involved in thionation
                        if checker.check_fg("Ketone", r): findings_json["atomic_checks"]["functional_groups"].append("Ketone")
                        if checker.check_fg("Aldehyde", r): findings_json["atomic_checks"]["functional_groups"].append("Aldehyde")
                        if checker.check_fg("Primary amide", r): findings_json["atomic_checks"]["functional_groups"].append("Primary amide")
                        if checker.check_fg("Secondary amide", r): findings_json["atomic_checks"]["functional_groups"].append("Secondary amide")
                        if checker.check_fg("Thiocarbonyl", product_smiles): findings_json["atomic_checks"]["functional_groups"].append("Thiocarbonyl")
                        if checker.check_fg("Thioamide", product_smiles): findings_json["atomic_checks"]["functional_groups"].append("Thioamide")
                        break

                if is_thionation:
                    reaction_sequence.append(("thionation", depth))
                    findings_json["atomic_checks"]["named_reactions"].append("thionation")
                    print(f"Thionation detected at depth {depth}")

                # Check for ring opening
                # In forward direction: ring → open chain
                # So in retrosynthesis: open chain in product, ring in reactants
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles]
                product_mol = Chem.MolFromSmiles(product_smiles)

                if all(r is not None for r in reactant_mols) and product_mol is not None:
                    # Count rings in reactants and product
                    reactant_ring_count = sum(
                        [mol.GetRingInfo().NumRings() for mol in reactant_mols]
                    )
                    product_ring_count = product_mol.GetRingInfo().NumRings()

                    # Check if any rings were opened (more rings in reactants than in product)
                    if product_ring_count < reactant_ring_count:
                        reaction_sequence.append(("ring_opening", depth))
                        findings_json["atomic_checks"]["named_reactions"].append("ring_opening")
                        print(f"Ring opening detected at depth {depth}")

                        # Check if both thionation and ring opening occur in the same reaction
                        if is_thionation:
                            reaction_sequence.append(("thionation_and_ring_opening", depth))
                            findings_json["atomic_checks"]["named_reactions"].append("thionation_and_ring_opening")
                            print(f"Thionation and ring opening in same reaction at depth {depth}")

            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Process children
        for child in node.get("children", []):
            # Determine the new depth based on the current node's type
            new_depth = depth if node["type"] == "reaction" else depth + 1
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    # Check if both thionation and ring opening are present in the route
    has_thionation = any(reaction[0] == "thionation" for reaction in reaction_sequence)
    has_ring_opening = any(reaction[0] == "ring_opening" for reaction in reaction_sequence)
    has_combined = any(
        reaction[0] == "thionation_and_ring_opening" for reaction in reaction_sequence
    )

    # If we found a reaction that does both thionation and ring opening, return True
    if has_combined:
        print("Found reaction that performs both thionation and ring opening")
        result = True
        # Add the combined structural constraint if it's explicitly defined in the strategy
        # For this specific strategy, 'thionation_and_ring_opening' is a named reaction, not a structural constraint.
        # The structural constraints are 'co-occurrence' and 'sequence'.

    if has_thionation and has_ring_opening:
        print("Both thionation and ring opening detected in the route")
        # Add co-occurrence constraint
        findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"targets": ["thionation", "ring_opening"]}})

        # In forward synthesis: thionation → ring opening
        # In retrosynthesis: we see ring opening at lower depth, thionation at higher depth
        for i in range(len(reaction_sequence) - 1):
            for j in range(i + 1, len(reaction_sequence)):
                # Check if thionation is at a higher depth (earlier in synthesis) than ring opening
                if (
                    reaction_sequence[i][0] == "thionation"
                    and reaction_sequence[j][0] == "ring_opening"
                    and reaction_sequence[i][1] > reaction_sequence[j][1]
                ):
                    print(
                        f"Found thionation at depth {reaction_sequence[i][1]} and ring opening at depth {reaction_sequence[j][1]} in retrosynthesis"
                    )
                    print(
                        f"This corresponds to thionation followed by ring opening in forward synthesis"
                    )
                    result = True
                    # Add sequence constraint
                    findings_json["structural_constraints"].append({"type": "sequence", "details": {"before": "thionation", "after": "ring_opening"}})
                    break # Found the sequence, no need to check further pairs
            if result: # If sequence found, break outer loop too
                break

    return result, findings_json

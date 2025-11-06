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
    This function detects nitrile to thioamide functional group conversion.
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

    conversion_detected = False

    def dfs_traverse(node, depth, max_depth):
        nonlocal conversion_detected, findings_json

        if conversion_detected:
            return

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                print(f"Analyzing reaction: {rsmi}")

                # First check: nitrile to thioamide conversion
                # Check for nitrile in reactants
                nitrile_reactant_idx = None
                for idx, r_smiles in enumerate(reactants_smiles):
                    try:
                        if checker.check_fg("Nitrile", r_smiles):
                            nitrile_reactant_idx = idx
                            if "Nitrile" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Nitrile")
                            print(f"Found nitrile in reactant {idx}: {r_smiles}")
                            break
                    except Exception as e:
                        print(f"Error checking nitrile in reactant: {e}")
                        continue

                # If nitrile found in reactants, check for thioamide in product
                if nitrile_reactant_idx is not None:
                    try:
                        # Check if product has thioamide
                        has_thioamide_in_product = checker.check_fg("Thioamide", product_smiles)

                        print(f"Product has thioamide: {has_thioamide_in_product}")

                        # Verify other reactants don't already have thioamide
                        other_reactants = [
                            r for i, r in enumerate(reactants_smiles) if i != nitrile_reactant_idx
                        ]
                        has_thioamide_in_other_reactants = any(
                            checker.check_fg("Thioamide", r_smiles) for r_smiles in other_reactants
                        )

                        print(f"Other reactants have thioamide: {has_thioamide_in_other_reactants}")

                        if has_thioamide_in_product:
                            if "Thioamide" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Thioamide")

                        if (
                            has_thioamide_in_product
                            and not has_thioamide_in_other_reactants
                        ):
                            print(
                                "Functional group pattern matches nitrile to thioamide conversion"
                            )
                            findings_json["structural_constraints"].append({
                                "type": "co-occurrence",
                                "details": {
                                    "targets": [
                                        "Nitrile",
                                        "Thioamide"
                                    ],
                                    "description": "A reactant must contain a Nitrile group and the corresponding product must contain a Thioamide group within the same reaction step."
                                }
                            })
                            findings_json["structural_constraints"].append({
                                "type": "negation",
                                "details": {
                                    "target": "Thioamide",
                                    "description": "Other reactants in the same step must not contain a Thioamide group."
                                }
                            })

                            # Get the nitrile reactant and product molecules
                            nitrile_reactant = reactants_smiles[nitrile_reactant_idx]
                            nitrile_mol = Chem.MolFromSmiles(nitrile_reactant)
                            product_mol = Chem.MolFromSmiles(product_smiles)

                            if nitrile_mol and product_mol:
                                # Get atom indices for nitrile in reactant and thioamide in product
                                nitrile_indices = checker.get_fg_atom_indices(
                                    "Nitrile", nitrile_reactant
                                )
                                thioamide_indices = checker.get_fg_atom_indices(
                                    "Thioamide", product_smiles
                                )

                                print(f"Nitrile indices: {nitrile_indices}")
                                print(f"Thioamide indices: {thioamide_indices}")

                                if nitrile_indices and thioamide_indices:
                                    # Extract atom mapping numbers from the nitrile carbon
                                    nitrile_carbon_maps = []
                                    for group_indices in nitrile_indices:
                                        for atom_idx in group_indices:
                                            atom = nitrile_mol.GetAtomWithIdx(atom_idx)
                                            if atom.GetSymbol() == "C" and atom.HasProp(
                                                "molAtomMapNumber"
                                            ):
                                                map_num = atom.GetProp("molAtomMapNumber")
                                                nitrile_carbon_maps.append(map_num)
                                                print(
                                                    f"Found nitrile carbon with map number: {map_num}"
                                                )

                                    # Extract atom mapping numbers from the thioamide carbon
                                    thioamide_carbon_maps = []
                                    for group_indices in thioamide_indices:
                                        for atom_idx in group_indices:
                                            atom = product_mol.GetAtomWithIdx(atom_idx)
                                            if atom.GetSymbol() == "C" and atom.HasProp(
                                                "molAtomMapNumber"
                                            ):
                                                map_num = atom.GetProp("molAtomMapNumber")
                                                thioamide_carbon_maps.append(map_num)
                                                print(
                                                    f"Found thioamide carbon with map number: {map_num}"
                                                )

                                    # Check if any nitrile carbon maps to a thioamide carbon
                                    for map_num in nitrile_carbon_maps:
                                        if map_num in thioamide_carbon_maps:
                                            print(
                                                f"Confirmed nitrile carbon (map {map_num}) converted to thioamide carbon"
                                            )
                                            conversion_detected = True
                                            return

                    except Exception as e:
                        print(f"Error checking thioamide in product: {e}")

        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                dfs_traverse(child, depth, max_depth) # Depth remains the same
            else:
                dfs_traverse(child, depth + 1, max_depth) # Depth increases

    dfs_traverse(route, 0, float('inf')) # Initial call with depth 0 and infinite max_depth
    return conversion_detected, findings_json
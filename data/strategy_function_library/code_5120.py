from typing import Tuple, Dict, List
import copy
from rdkit import Chem

# Refactoring for Enumeration: Isolate and refine lists
REDUCTION_REACTION_NAMES = [
    "Reduction of aldehydes and ketones to alcohols",
    "Reduction of ketone to secondary alcohol",
    "Reduction of carboxylic acid to primary alcohol",
    "Reduction of ester to primary alcohol",
]

CARBONYL_GROUPS = [
    "Ketone",
    "Aldehyde",
    "Ester",
    "Carboxylic acid",
]

ALCOHOL_GROUPS = [
    "Primary alcohol",
    "Secondary alcohol",
]

# Assuming checker is initialized globally as in the original context
# from synth_strategy.utils.check import Check
# from synth_strategy.utils import fuzzy_dict, check
# from pathlib import Path
root_data = Path(__file__).parent.parent
# fg_args = {"file_path": f"{root_data}/patterns/functional_groups.json", "value_field": "pattern", "key_field": "name"}
# reaction_class_args = {"file_path": f"{root_data}/patterns/smirks.json", "value_field": "smirks", "key_field": "name"}
# ring_smiles_args = {"file_path": f"{root_data}/patterns/chemical_rings_smiles.json", "value_field": "smiles", "key_field": "name"}
# functional_groups = fuzzy_dict.FuzzyDict.from_json(**fg_args)
# reaction_classes = fuzzy_dict.FuzzyDict.from_json(**reaction_class_args)
# ring_smiles = fuzzy_dict.FuzzyDict.from_json(**ring_smiles_args)
# checker = check.Check(fg_dict=functional_groups, reaction_dict=reaction_classes, ring_dict=ring_smiles)

# Placeholder for checker if not run in original environment
class MockChecker:
    def check_reaction(self, name, rsmi): return False
    def check_fg(self, fg_name, smi): return False
    def get_fg_atom_indices(self, fg_name, smi): return []
checker = MockChecker()

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects if the synthetic route involves reduction of multiple carbonyl groups
    (ketone/aldehyde/ester/carboxylic acid) to hydroxyl groups in the final step of the synthesis.
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

    # Track if we found a late-stage carbonyl reduction
    found_late_stage_reduction = False

    def dfs_traverse(node, depth=0):
        nonlocal found_late_stage_reduction, findings_json

        if node["type"] == "reaction" and depth <= 2:  # Final or near-final reaction
            # Extract reactants and product from reaction SMILES
            if "mapped_reaction_smiles" in node.get("metadata", {}):
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                print(f"Analyzing late-stage reaction at depth {depth}: {rsmi}")
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Parse molecules
                try:
                    reactant_mols = [Chem.MolFromSmiles(smi) for smi in reactants_smiles if smi]
                    product_mol = Chem.MolFromSmiles(product_smiles)

                    if not all(reactant_mols) or not product_mol:
                        print("Failed to parse some molecules")
                        return

                    # Check if this is a known reduction reaction using the refined list
                    is_known_reduction = False
                    for name in REDUCTION_REACTION_NAMES:
                        if checker.check_reaction(name, rsmi):
                            is_known_reduction = True
                            if name not in findings_json["atomic_checks"]["named_reactions"]:
                                findings_json["atomic_checks"]["named_reactions"].append(name)
                            break

                    print(f"Is this a known carbonyl reduction reaction? {is_known_reduction}")

                    # Count carbonyl groups in reactants using the refined list
                    total_carbonyls_reactants = 0
                    carbonyl_atoms_in_reactants = {}  # Track atom indices with atom mapping

                    for reactant_idx, reactant_smi in enumerate(reactants_smiles):
                        for carbonyl_type in CARBONYL_GROUPS:
                            if checker.check_fg(carbonyl_type, reactant_smi):
                                if carbonyl_type not in findings_json["atomic_checks"]["functional_groups"]:
                                    findings_json["atomic_checks"]["functional_groups"].append(carbonyl_type)
                                carbonyl_indices = checker.get_fg_atom_indices(
                                    carbonyl_type, reactant_smi
                                )
                                total_carbonyls_reactants += len(carbonyl_indices)
                                print(
                                    f"Found {len(carbonyl_indices)} {carbonyl_type}s in reactant {reactant_smi}"
                                )

                                # Store the atom indices with their atom mapping
                                reactant_mol = Chem.MolFromSmiles(reactant_smi)
                                if reactant_mol:
                                    for atom_indices in carbonyl_indices:
                                        for atom_idx in atom_indices:
                                            if atom_idx < reactant_mol.GetNumAtoms():
                                                atom = reactant_mol.GetAtomWithIdx(atom_idx)
                                                if atom.HasProp("molAtomMapNumber"):
                                                    map_num = atom.GetProp("molAtomMapNumber")
                                                    carbonyl_atoms_in_reactants[map_num] = (
                                                        reactant_idx,
                                                        carbonyl_type,
                                                        atom_idx,
                                                    )

                    print(f"Total carbonyl groups in reactants: {total_carbonyls_reactants}")

                    # Count hydroxyl groups in product using the refined list
                    hydroxyl_count_product = 0
                    hydroxyl_atoms_in_product = {}  # Track atom indices with atom mapping

                    for alcohol_type in ALCOHOL_GROUPS:
                        if checker.check_fg(alcohol_type, product_smiles):
                            if alcohol_type not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append(alcohol_type)
                            alcohol_indices = checker.get_fg_atom_indices(
                                alcohol_type, product_smiles
                            )
                            hydroxyl_count_product += len(alcohol_indices)
                            print(f"Found {len(alcohol_indices)} {alcohol_type}s in product")

                            # Store the atom indices with their atom mapping
                            for atom_indices in alcohol_indices:
                                for atom_idx in atom_indices:
                                    if atom_idx < product_mol.GetNumAtoms():
                                        atom = product_mol.GetAtomWithIdx(atom_idx)
                                        if atom.HasProp("molAtomMapNumber"):
                                            map_num = atom.GetProp("molAtomMapNumber")
                                            hydroxyl_atoms_in_product[map_num] = (
                                                alcohol_type,
                                                atom_idx,
                                            )

                    print(f"Total hydroxyl groups in product: {hydroxyl_count_product}")

                    # Count hydroxyl groups in reactants using the refined list
                    hydroxyl_count_reactants = 0
                    hydroxyl_atoms_in_reactants = {}

                    for reactant_idx, reactant_smi in enumerate(reactants_smiles):
                        for alcohol_type in ALCOHOL_GROUPS:
                            if checker.check_fg(alcohol_type, reactant_smi):
                                if alcohol_type not in findings_json["atomic_checks"]["functional_groups"]:
                                    findings_json["atomic_checks"]["functional_groups"].append(alcohol_type)
                                alcohol_indices = checker.get_fg_atom_indices(
                                    alcohol_type, reactant_smi
                                )
                                hydroxyl_count_reactants += len(alcohol_indices)
                                print(
                                    f"Found {len(alcohol_indices)} {alcohol_type}s in reactant {reactant_smi}"
                                )

                                # Store the atom indices with their atom mapping
                                reactant_mol = Chem.MolFromSmiles(reactant_smi)
                                if reactant_mol:
                                    for atom_indices in alcohol_indices:
                                        for atom_idx in atom_indices:
                                            if atom_idx < reactant_mol.GetNumAtoms():
                                                atom = reactant_mol.GetAtomWithIdx(atom_idx)
                                                if atom.HasProp("molAtomMapNumber"):
                                                    map_num = atom.GetProp("molAtomMapNumber")
                                                    hydroxyl_atoms_in_reactants[map_num] = (
                                                        reactant_idx,
                                                        alcohol_type,
                                                        atom_idx,
                                                    )

                    print(f"Total hydroxyl groups in reactants: {hydroxyl_count_reactants}")

                    # Calculate the number of new hydroxyl groups formed
                    new_hydroxyls = hydroxyl_count_product - hydroxyl_count_reactants
                    print(f"New hydroxyl groups formed: {new_hydroxyls}")

                    # Check for carbonyl reduction by examining atom mapping
                    carbonyl_to_hydroxyl_conversions = 0

                    # Look for carbonyl carbons in reactants that become hydroxyl-bearing carbons in product
                    for map_num in carbonyl_atoms_in_reactants:
                        if map_num in hydroxyl_atoms_in_product:
                            print(f"Found carbonyl to hydroxyl conversion for atom map {map_num}")
                            carbonyl_to_hydroxyl_conversions += 1

                    print(
                        f"Confirmed carbonyl to hydroxyl conversions: {carbonyl_to_hydroxyl_conversions}"
                    )

                    # If there are confirmed carbonyl reductions or if it's a known reduction reaction with new hydroxyls
                    if carbonyl_to_hydroxyl_conversions >= 1 or (
                        is_known_reduction and new_hydroxyls >= 1
                    ):
                        # Add positional constraint if met
                        if {"type": "positional", "details": {"target": "carbonyl_reduction", "position": "near_final_stage (depth <= 2)"}} not in findings_json["structural_constraints"]:
                            findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "carbonyl_reduction", "position": "near_final_stage (depth <= 2)"}})

                        # Check if multiple carbonyls were reduced or if this is a significant reduction
                        if carbonyl_to_hydroxyl_conversions >= 2:
                            print(
                                f"Found late-stage carbonyl reduction: {carbonyl_to_hydroxyl_conversions} carbonyls reduced at depth {depth}"
                            )
                            found_late_stage_reduction = True
                            if {"type": "count", "details": {"target": "carbonyl_to_hydroxyl_conversion_in_step", "operator": ">=", "value": 2}} not in findings_json["structural_constraints"]:
                                findings_json["structural_constraints"].append({"type": "count", "details": {"target": "carbonyl_to_hydroxyl_conversion_in_step", "operator": ">=", "value": 2}})
                        
                        if carbonyl_to_hydroxyl_conversions >= 1 and depth <= 1:
                            print(
                                f"Found late-stage carbonyl reduction: {carbonyl_to_hydroxyl_conversions} carbonyls reduced at depth {depth}"
                            )
                            found_late_stage_reduction = True
                            if {"type": "co-occurrence", "details": {"targets": ["single_carbonyl_reduction", "very_late_stage_reaction (depth <= 1)"]}} not in findings_json["structural_constraints"]:
                                findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"targets": ["single_carbonyl_reduction", "very_late_stage_reaction (depth <= 1)"]}})

                        if is_known_reduction and new_hydroxyls >= 1:
                            if {"type": "co-occurrence", "details": {"targets": ["known_reduction_reaction", "hydroxyl_group_formation"]}} not in findings_json["structural_constraints"]:
                                findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"targets": ["known_reduction_reaction", "hydroxyl_group_formation"]}})

                        if is_known_reduction and (new_hydroxyls >= 2 or (new_hydroxyls >= 1 and depth <= 1)):
                            print(
                                f"Found late-stage carbonyl reduction via known reaction: {new_hydroxyls} new hydroxyls at depth {depth}"
                            )
                            found_late_stage_reduction = True
                            if new_hydroxyls >= 2:
                                if {"type": "count", "details": {"target": "carbonyl_to_hydroxyl_conversion_in_step", "operator": ">=", "value": 2}} not in findings_json["structural_constraints"]:
                                    findings_json["structural_constraints"].append({"type": "count", "details": {"target": "carbonyl_to_hydroxyl_conversion_in_step", "operator": ">=", "value": 2}})
                            if new_hydroxyls >= 1 and depth <= 1:
                                if {"type": "co-occurrence", "details": {"targets": ["single_carbonyl_reduction", "very_late_stage_reaction (depth <= 1)"]}} not in findings_json["structural_constraints"]:
                                    findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"targets": ["single_carbonyl_reduction", "very_late_stage_reaction (depth <= 1)"]}})

                except Exception as e:
                    print(f"Error in processing reaction: {e}")

        # Process children
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                # If current node is a reaction, depth remains the same for its children (chemicals)
                dfs_traverse(child, depth)
            else:
                # If current node is a chemical, depth increases for its children (reactions)
                dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)
    print(f"Final result: {found_late_stage_reduction}")

    return found_late_stage_reduction, findings_json
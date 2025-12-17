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
    Detects a synthetic strategy where a nitrile is converted to a trifluoromethyl ketone
    or trifluoromethyl alcohol while maintaining a halogenated aromatic scaffold.
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

    # Track key features
    has_nitrile_intermediate = False
    has_aldehyde_intermediate = False
    has_alcohol_intermediate = False
    has_final_trifluoromethyl_compound = False
    has_brominated_aromatic = False

    # Track reaction sequence
    reaction_sequence = []

    def dfs_traverse(node, depth=0):
        nonlocal has_nitrile_intermediate, has_aldehyde_intermediate
        nonlocal has_alcohol_intermediate, has_final_trifluoromethyl_compound
        nonlocal has_brominated_aromatic, reaction_sequence, findings_json

        if node["type"] == "mol":
            # Check molecule features
            mol_smiles = node["smiles"]

            # Check for brominated aromatic
            if checker.check_fg("Aromatic halide", mol_smiles):
                has_brominated_aromatic = True
                if "Aromatic halide" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Aromatic halide")
                print(f"Found brominated aromatic at depth {depth}: {mol_smiles}")

            # Check for nitrile
            if checker.check_fg("Nitrile", mol_smiles):
                has_nitrile_intermediate = True
                if "Nitrile" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Nitrile")
                print(f"Found nitrile intermediate at depth {depth}: {mol_smiles}")

            # Check for aldehyde
            if checker.check_fg("Aldehyde", mol_smiles):
                has_aldehyde_intermediate = True
                if "Aldehyde" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Aldehyde")
                print(f"Found aldehyde intermediate at depth {depth}: {mol_smiles}")
            if checker.check_fg("Formaldehyde", mol_smiles):
                has_aldehyde_intermediate = True
                if "Formaldehyde" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Formaldehyde")
                print(f"Found aldehyde intermediate at depth {depth}: {mol_smiles}")

            # Check for alcohol
            if checker.check_fg("Primary alcohol", mol_smiles):
                has_alcohol_intermediate = True
                if "Primary alcohol" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Primary alcohol")
                print(f"Found alcohol intermediate at depth {depth}: {mol_smiles}")
            if checker.check_fg("Secondary alcohol", mol_smiles):
                has_alcohol_intermediate = True
                if "Secondary alcohol" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Secondary alcohol")
                print(f"Found alcohol intermediate at depth {depth}: {mol_smiles}")
            if checker.check_fg("Tertiary alcohol", mol_smiles):
                has_alcohol_intermediate = True
                if "Tertiary alcohol" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Tertiary alcohol")
                print(f"Found alcohol intermediate at depth {depth}: {mol_smiles}")
            if checker.check_fg("Aromatic alcohol", mol_smiles):
                has_alcohol_intermediate = True
                if "Aromatic alcohol" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Aromatic alcohol")
                print(f"Found alcohol intermediate at depth {depth}: {mol_smiles}")

            # Check for trifluoromethyl compounds at depth 0 (final product)
            if depth == 0 and checker.check_fg("Trifluoro group", mol_smiles):
                if "Trifluoro group" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Trifluoro group")

                if checker.check_fg("Ketone", mol_smiles):
                    has_final_trifluoromethyl_compound = True
                    if "Ketone" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Ketone")
                    print(f"Found final trifluoromethyl compound (Ketone) at depth {depth}: {mol_smiles}")
                elif (checker.check_fg("Primary alcohol", mol_smiles) or
                      checker.check_fg("Secondary alcohol", mol_smiles) or
                      checker.check_fg("Tertiary alcohol", mol_smiles) or
                      checker.check_fg("Aromatic alcohol", mol_smiles)):
                    has_final_trifluoromethyl_compound = True
                    # FGs already added above if found
                    print(f"Found final trifluoromethyl compound (Alcohol) at depth {depth}: {mol_smiles}")

        elif node["type"] == "reaction":
            # Analyze reaction
            if "mapped_reaction_smiles" in node.get("metadata", {}):
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0]
                product = rsmi.split(">")[-1]

                # Check for nitrile to aldehyde conversion
                if checker.check_fg("Nitrile", reactants) and (
                    checker.check_fg("Aldehyde", product)
                    or checker.check_fg("Formaldehyde", product)
                ):
                    reaction_sequence.append("nitrile_to_aldehyde")
                    if "nitrile_to_aldehyde" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("nitrile_to_aldehyde")
                    print(f"Detected nitrile to aldehyde reaction: {rsmi}")

                # Check for aldehyde to alcohol conversion
                elif (
                    checker.check_fg("Aldehyde", reactants)
                    or checker.check_fg("Formaldehyde", reactants)
                ) and (
                    checker.check_fg("Primary alcohol", product)
                    or checker.check_fg("Secondary alcohol", product)
                    or checker.check_fg("Tertiary alcohol", product)
                    or checker.check_fg("Aromatic alcohol", product)
                ):
                    # Check if the product also has a trifluoro group
                    if checker.check_fg("Trifluoro group", product):
                        reaction_sequence.append("aldehyde_to_trifluoromethyl_alcohol")
                        if "aldehyde_to_trifluoromethyl_alcohol" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("aldehyde_to_trifluoromethyl_alcohol")
                        print(f"Detected aldehyde to trifluoromethyl alcohol reaction: {rsmi}")
                    else:
                        reaction_sequence.append("aldehyde_to_alcohol")
                        if "aldehyde_to_alcohol" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("aldehyde_to_alcohol")
                        print(f"Detected aldehyde to alcohol reaction: {rsmi}")

                # Check for alcohol to trifluoromethyl ketone conversion
                elif (
                    (
                        checker.check_fg("Primary alcohol", reactants)
                        or checker.check_fg("Secondary alcohol", reactants)
                        or checker.check_fg("Tertiary alcohol", reactants)
                        or checker.check_fg("Aromatic alcohol", reactants)
                    )
                    and checker.check_fg("Ketone", product)
                    and checker.check_fg("Trifluoro group", product)
                ):
                    reaction_sequence.append("alcohol_to_trifluoromethyl_ketone")
                    if "alcohol_to_trifluoromethyl_ketone" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("alcohol_to_trifluoromethyl_ketone")
                    print(f"Detected alcohol to trifluoromethyl ketone reaction: {rsmi}")

                # Check for direct nitrile to trifluoromethyl ketone conversion
                elif (
                    checker.check_fg("Nitrile", reactants)
                    and checker.check_fg("Ketone", product)
                    and checker.check_fg("Trifluoro group", product)
                ):
                    reaction_sequence.append("nitrile_to_trifluoromethyl_ketone")
                    if "nitrile_to_trifluoromethyl_ketone" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("nitrile_to_trifluoromethyl_ketone")
                    print(f"Detected direct nitrile to trifluoromethyl ketone reaction: {rsmi}")

                # Check for direct nitrile to trifluoromethyl alcohol conversion
                elif (
                    checker.check_fg("Nitrile", reactants)
                    and (
                        checker.check_fg("Primary alcohol", product)
                        or checker.check_fg("Secondary alcohol", product)
                        or checker.check_fg("Tertiary alcohol", product)
                        or checker.check_fg("Aromatic alcohol", product)
                    )
                    and checker.check_fg("Trifluoro group", product)
                ):
                    reaction_sequence.append("nitrile_to_trifluoromethyl_alcohol")
                    if "nitrile_to_trifluoromethyl_alcohol" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("nitrile_to_trifluoromethyl_alcohol")
                    print(f"Detected direct nitrile to trifluoromethyl alcohol reaction: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else: # Assuming 'mol' type or any other type that should increase depth
                dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if the strategy is present - be more flexible with the sequence
    has_valid_sequence = False

    # Check for the full pathway: nitrile -> aldehyde -> trifluoromethyl alcohol
    if (
        "nitrile_to_aldehyde" in reaction_sequence
        and "aldehyde_to_trifluoromethyl_alcohol" in reaction_sequence
    ):
        has_valid_sequence = True
        # Add structural constraint for this pathway
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "description": "The route must contain one of several valid reaction pathways converting a nitrile to a trifluoromethyl ketone or alcohol.",
                "pathway_logic": {
                    "operator": "AND",
                    "operands": [
                        "nitrile_to_aldehyde",
                        "aldehyde_to_trifluoromethyl_alcohol"
                    ]
                }
            }
        })

    # Check for the full pathway: nitrile -> aldehyde -> alcohol -> trifluoromethyl ketone
    elif (
        "nitrile_to_aldehyde" in reaction_sequence
        and "aldehyde_to_alcohol" in reaction_sequence
        and "alcohol_to_trifluoromethyl_ketone" in reaction_sequence
    ):
        has_valid_sequence = True
        # Add structural constraint for this pathway
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "description": "The route must contain one of several valid reaction pathways converting a nitrile to a trifluoromethyl ketone or alcohol.",
                "pathway_logic": {
                    "operator": "AND",
                    "operands": [
                        "nitrile_to_aldehyde",
                        "aldehyde_to_alcohol",
                        "alcohol_to_trifluoromethyl_ketone"
                    ]
                }
            }
        })

    # Check for direct nitrile to trifluoromethyl ketone or alcohol
    elif (
        "nitrile_to_trifluoromethyl_ketone" in reaction_sequence
        or "nitrile_to_trifluoromethyl_alcohol" in reaction_sequence
    ):
        has_valid_sequence = True
        # Add structural constraint for this pathway
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "description": "The route must contain one of several valid reaction pathways converting a nitrile to a trifluoromethyl ketone or alcohol.",
                "pathway_logic": {
                    "operator": "OR",
                    "operands": [
                        "nitrile_to_trifluoromethyl_ketone",
                        "nitrile_to_trifluoromethyl_alcohol"
                    ]
                }
            }
        })

    strategy_present = (
        has_nitrile_intermediate
        and has_final_trifluoromethyl_compound
        and has_brominated_aromatic
        and has_valid_sequence
    )

    # Add structural constraints if conditions are met
    if has_nitrile_intermediate and has_brominated_aromatic:
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "description": "The synthesis route must involve molecules containing Nitrile and Aromatic halide functional groups.",
                "targets": [
                    "Nitrile",
                    "Aromatic halide"
                ]
            }
        })

    if has_final_trifluoromethyl_compound:
        # This constraint is more complex as it depends on which specific FGs were found in the final product
        # We'll add a simplified version based on the flag, assuming the FGs were already added to atomic_checks
        # when has_final_trifluoromethyl_compound was set.
        # The full logic for 'target_logic' would require re-evaluating the FGs of the final product here.
        # For simplicity, we'll just add the positional constraint if the flag is true.
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "description": "The final product must contain a Trifluoro group and either a Ketone or an Alcohol group.",
                "position": "last_stage",
                "target_logic": {
                    "operator": "AND",
                    "operands": [
                        "Trifluoro group",
                        {
                            "operator": "OR",
                            "operands": [
                                "Ketone",
                                "Primary alcohol",
                                "Secondary alcohol",
                                "Tertiary alcohol",
                                "Aromatic alcohol"
                            ]
                        }
                    ]
                }
            }
        })

    print(f"Trifluoromethyl compound from nitrile strategy detected: {strategy_present}")
    print(
        f"Key intermediates: Nitrile: {has_nitrile_intermediate}, Aldehyde: {has_aldehyde_intermediate}, Alcohol: {has_alcohol_intermediate}"
    )
    print(f"Final trifluoromethyl compound: {has_final_trifluoromethyl_compound}")
    print(f"Brominated aromatic maintained: {has_brominated_aromatic}")
    print(f"Reaction sequence: {reaction_sequence}")

    return strategy_present, findings_json

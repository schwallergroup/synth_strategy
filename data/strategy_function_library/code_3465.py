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


AMIDE_FORMATION_REACTIONS = [
    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
    "Acyl chloride with secondary amine to amide",
    "Carboxylic acid with primary amine to amide",
    "Ester with primary amine to amide",
    "Ester with secondary amine to amide",
    "Schotten-Baumann_amide",
]

AMINE_DEPROTECTION_REACTIONS = [
    "Boc amine deprotection",
    "N-glutarimide deprotection",
    "Phthalimide deprotection",
    "Tert-butyl deprotection of amine",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a specific synthetic sequence where an amine on a piperidine scaffold is deprotected and then subsequently acylated to form an amide. The strategy is flagged only if a trifluoromethyl group is present on one of the reactants in the amide formation step. The function uses specific lists of known reaction types for both the deprotection (e.g., 'Boc amine deprotection') and amide formation (e.g., 'Acyl chloride with secondary amine to amide') steps.
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

    # Track if we found the key reactions
    amide_formation_found = False
    amine_deprotection_found = False
    piperidine_scaffold_present = False
    trifluoromethyl_present = False

    # Track the sequence of reactions
    amide_formation_depth = -1
    deprotection_depth = -1

    # Track the final product
    final_product_smiles = None

    def dfs_traverse(node, depth=0):
        nonlocal amide_formation_found, amine_deprotection_found, piperidine_scaffold_present
        nonlocal trifluoromethyl_present, amide_formation_depth, deprotection_depth, final_product_smiles
        nonlocal findings_json

        if node["type"] == "mol":
            mol_smiles = node["smiles"]

            # Check if molecule contains piperidine scaffold
            if checker.check_ring("piperidine", mol_smiles):
                piperidine_scaffold_present = True
                if "piperidine" not in findings_json["atomic_checks"]["ring_systems"]:
                    findings_json["atomic_checks"]["ring_systems"].append("piperidine")
                print(f"Found piperidine scaffold in molecule: {mol_smiles}")

            # Check if molecule contains trifluoromethyl group
            if checker.check_fg("Trifluoro group", mol_smiles):
                trifluoromethyl_present = True
                if "Trifluoro group" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Trifluoro group")
                print(f"Found trifluoromethyl group in molecule: {mol_smiles}")

            # Store final product SMILES (depth 0)
            if depth == 0:
                final_product_smiles = mol_smiles
                print(f"Final product: {mol_smiles}")

        elif node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                print(f"\nAnalyzing reaction at depth {depth}: {rsmi}")

                # Check for amide formation reactions
                is_amide_formation = False
                for r_name in AMIDE_FORMATION_REACTIONS:
                    if checker.check_reaction(r_name, rsmi):
                        is_amide_formation = True
                        if r_name not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(r_name)
                        break

                if is_amide_formation:
                    print(f"Detected standard amide formation reaction")
                    # Check if the product contains both piperidine and an amide
                    product_has_piperidine = checker.check_ring("piperidine", product_smiles)
                    product_has_amide = False
                    if checker.check_fg("Secondary amide", product_smiles):
                        product_has_amide = True
                        if "Secondary amide" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Secondary amide")
                    if checker.check_fg("Primary amide", product_smiles):
                        product_has_amide = True
                        if "Primary amide" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Primary amide")

                    if product_has_piperidine and product_has_amide:
                        if "piperidine" not in findings_json["atomic_checks"]["ring_systems"]:
                            findings_json["atomic_checks"]["ring_systems"].append("piperidine")

                        # Check if any reactant has a trifluoromethyl group
                        has_trifluoromethyl_reactant = False
                        for r in reactants_smiles:
                            if checker.check_fg("Trifluoro group", r):
                                has_trifluoromethyl_reactant = True
                                if "Trifluoro group" not in findings_json["atomic_checks"]["functional_groups"]:
                                    findings_json["atomic_checks"]["functional_groups"].append("Trifluoro group")
                                break

                        # Check if any reactant has a piperidine
                        has_piperidine_reactant = False
                        for r in reactants_smiles:
                            if checker.check_ring("piperidine", r):
                                has_piperidine_reactant = True
                                if "piperidine" not in findings_json["atomic_checks"]["ring_systems"]:
                                    findings_json["atomic_checks"]["ring_systems"].append("piperidine")
                                break

                        if has_piperidine_reactant and has_trifluoromethyl_reactant:
                            print(f"Found relevant amide formation at depth {depth}: {rsmi}")
                            amide_formation_found = True
                            amide_formation_depth = depth
                            # Structural constraint: co-occurrence of amide_formation and reactant_with_Trifluoro_group
                            findings_json["structural_constraints"].append({
                                "type": "co-occurrence",
                                "details": {
                                    "scope": "reaction_step",
                                    "targets": [
                                        "amide_formation",
                                        "reactant_with_Trifluoro_group"
                                    ]
                                }
                            })

                # Check for amine deprotection reactions
                is_deprotection = False
                for r_name in AMINE_DEPROTECTION_REACTIONS:
                    if checker.check_reaction(r_name, rsmi):
                        is_deprotection = True
                        if r_name not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(r_name)
                        break

                if is_deprotection:
                    print(f"Detected standard deprotection/reductive amination reaction")
                    # Check if the product has a piperidine with a primary or secondary amine
                    product_has_piperidine = checker.check_ring("piperidine", product_smiles)
                    product_has_amine = False
                    if checker.check_fg("Primary amine", product_smiles):
                        product_has_amine = True
                        if "Primary amine" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Primary amine")
                    if checker.check_fg("Secondary amine", product_smiles):
                        product_has_amine = True
                        if "Secondary amine" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Secondary amine")

                    if product_has_piperidine and product_has_amine:
                        if "piperidine" not in findings_json["atomic_checks"]["ring_systems"]:
                            findings_json["atomic_checks"]["ring_systems"].append("piperidine")
                        print(f"Found relevant amine deprotection at depth {depth}: {rsmi}")
                        amine_deprotection_found = True
                        deprotection_depth = depth

            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                # If current node is a reaction, depth remains the same for its children (reactants)
                dfs_traverse(child, depth)
            else:
                # If current node is a chemical (mol), depth increases for its children (reactions)
                dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)

    # Final check: verify that the final product has both piperidine and trifluoromethyl
    if final_product_smiles:
        if checker.check_ring("piperidine", final_product_smiles):
            if "piperidine" not in findings_json["atomic_checks"]["ring_systems"]:
                findings_json["atomic_checks"]["ring_systems"].append("piperidine")
        if checker.check_fg("Trifluoro group", final_product_smiles):
            if "Trifluoro group" not in findings_json["atomic_checks"]["functional_groups"]:
                findings_json["atomic_checks"]["functional_groups"].append("Trifluoro group")

        if checker.check_ring("piperidine", final_product_smiles) and checker.check_fg(
            "Trifluoro group", final_product_smiles
        ):
            print(
                "Confirmed: Final product contains both piperidine scaffold and trifluoromethyl group"
            )
        else:
            print("Final product does not contain both required structural features")
            piperidine_scaffold_present = checker.check_ring("piperidine", final_product_smiles)
            trifluoromethyl_present = checker.check_fg("Trifluoro group", final_product_smiles)

    # Check if the sequence is correct (amide formation should be at a lower depth than deprotection)
    correct_sequence = (
        amide_formation_depth != -1
        and deprotection_depth != -1
        and amide_formation_depth < deprotection_depth
    )

    if correct_sequence:
        print(
            f"Correct sequence: Deprotection (depth {deprotection_depth}) followed by amide formation (depth {amide_formation_depth})"
        )
        # Structural constraint: sequence
        findings_json["structural_constraints"].append({
            "type": "sequence",
            "details": {
                "before": "amine_deprotection",
                "after": "amide_formation"
            }
        })
    else:
        print(
            f"Incorrect sequence: Amide formation depth = {amide_formation_depth}, Deprotection depth = {deprotection_depth}"
        )

    # Check if all required elements of the strategy are present
    strategy_present = (
        amide_formation_found
        and amine_deprotection_found
        and piperidine_scaffold_present
        and trifluoromethyl_present
        and correct_sequence
    )

    if strategy_present:
        print(
            "Detected strategy: Late-stage amide formation after deprotection on piperidine scaffold"
        )
        # Structural constraint: co-occurrence of amide_formation, amine_deprotection, piperidine
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "amide_formation",
                    "amine_deprotection",
                    "piperidine"
                ]
            }
        })
    else:
        print(
            f"Strategy not fully detected: amide={amide_formation_found}, deprotection={amine_deprotection_found}, piperidine={piperidine_scaffold_present}, CF3={trifluoromethyl_present}, sequence={correct_sequence}"
        )

    return strategy_present, findings_json

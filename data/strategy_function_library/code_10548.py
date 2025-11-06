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
    This function detects a strategy involving late-stage deprotection of
    t-butyl protected groups (Boc for amine, t-butyl for carboxylic acid).
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

    # Track if we found the relevant deprotection reactions
    boc_deprotection_found = False
    tbutyl_ester_deprotection_found = False

    # Track the depth of the deprotections to ensure they're late-stage
    boc_deprotection_depth = -1
    tbutyl_ester_deprotection_depth = -1

    # Track maximum depth to determine what "late-stage" means
    max_depth = 0

    result = False # Initialize the main boolean result

    def dfs_traverse(node, depth=0):
        nonlocal boc_deprotection_found, tbutyl_ester_deprotection_found
        nonlocal boc_deprotection_depth, tbutyl_ester_deprotection_depth
        nonlocal max_depth, findings_json

        # Update max depth
        max_depth = max(max_depth, depth)

        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for Boc deprotection
                if checker.check_reaction("Boc amine deprotection", rsmi):
                    print(f"Found Boc deprotection reaction at depth {depth}")
                    boc_deprotection_found = True
                    boc_deprotection_depth = depth
                    if "Boc amine deprotection" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("Boc amine deprotection")
                elif any(checker.check_fg("Boc", r) for r in reactants) and not checker.check_fg(
                    "Boc", product
                ):
                    # Check if reactant has Boc but product doesn't
                    # Also verify product has primary or secondary amine
                    if checker.check_fg("Primary amine", product):
                        print(f"Found Boc deprotection (by FG analysis) at depth {depth}")
                        boc_deprotection_found = True
                        boc_deprotection_depth = depth
                        if "Boc" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Boc")
                        if "Primary amine" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Primary amine")
                    if checker.check_fg(
                        "Secondary amine", product
                    ):
                        print(f"Found Boc deprotection (by FG analysis) at depth {depth}")
                        boc_deprotection_found = True
                        boc_deprotection_depth = depth
                        if "Boc" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Boc")
                        if "Secondary amine" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Secondary amine")

                # Alternative check for t-butyl ester deprotection by functional group analysis
                if not tbutyl_ester_deprotection_found:
                    tbutyl_pattern = Chem.MolFromSmarts("CC(C)(C)OC(=O)")

                    # Check if any reactant has t-butyl ester
                    reactant_has_tbutyl = False
                    for r in reactants:
                        try:
                            mol = Chem.MolFromSmiles(r)
                            if (
                                mol
                                and mol.HasSubstructMatch(tbutyl_pattern)
                                and checker.check_fg("Ester", r)
                            ):
                                reactant_has_tbutyl = True
                                if "Ester" not in findings_json["atomic_checks"]["functional_groups"]:
                                    findings_json["atomic_checks"]["functional_groups"].append("Ester")
                                break
                        except Exception as e:
                            print(f"Error checking reactant {r}: {e}")
                            continue

                    # Check if product has carboxylic acid but no t-butyl ester
                    try:
                        product_mol = Chem.MolFromSmiles(product)
                        product_has_acid = checker.check_fg("Carboxylic acid", product)
                        if product_has_acid and "Carboxylic acid" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Carboxylic acid")

                        product_has_tbutyl_ester = False

                        if product_mol:
                            product_has_tbutyl_ester = product_mol.HasSubstructMatch(
                                tbutyl_pattern
                            ) and checker.check_fg("Ester", product)

                        if (
                            reactant_tbutyl_ester_found
                            and product_has_acid
                            and not product_has_tbutyl_ester
                        ):
                            print(
                                f"Found t-butyl ester deprotection (by FG analysis) at depth {depth}"
                            )
                            tbutyl_ester_deprotection_found = True
                            tbutyl_ester_deprotection_depth = depth
                            if "t-butyl ester deprotection" not in findings_json["atomic_checks"]["named_reactions"]:
                                findings_json["atomic_checks"]["named_reactions"].append("t-butyl ester deprotection")
                    except Exception as e:
                        print(f"Error checking product {product}: {e}")
                        pass

            except Exception as e:
                print(f"Error processing reaction at depth {depth}: {e}")

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                # Depth remains the same when traversing from a reaction node to a chemical node
                dfs_traverse(child, depth)
            else:
                # Depth increases when traversing from a chemical node to a reaction node
                dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    print(f"Max depth: {max_depth}")
    print(f"Boc deprotection found: {boc_deprotection_found} at depth {boc_deprotection_depth}")
    print(
        f"t-butyl ester deprotection found: {tbutyl_ester_deprotection_found} at depth {tbutyl_ester_deprotection_depth}"
    )

    # Define what "late-stage" means - typically first third of the synthesis
    # Use a more lenient threshold to capture more potential cases
    late_stage_threshold = max(1, max_depth // 3 + 1)  # Adding 1 to include more reactions
    print(f"Late stage threshold: {late_stage_threshold}")

    # Check if both deprotections were found and they're in the late stage (low depth)
    if boc_deprotection_found and tbutyl_ester_deprotection_found:
        # Check if both deprotections occur in the late stage of the synthesis
        if (
            boc_deprotection_depth <= late_stage_threshold
            and tbutyl_ester_deprotection_depth <= late_stage_threshold
        ):
            print(
                "Strategy detected: Late-stage deprotection of orthogonally protected functional groups"
            )
            result = True
            if {"type": "positional", "details": {"target": "Boc amine deprotection", "position": "late_stage_first_third"}} not in findings_json["structural_constraints"]:
                findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "Boc amine deprotection", "position": "late_stage_first_third"}})
            if {"type": "positional", "details": {"target": "t-butyl ester deprotection", "position": "late_stage_first_third"}} not in findings_json["structural_constraints"]:
                findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "t-butyl ester deprotection", "position": "late_stage_first_third"}})

    # If only one deprotection was found in late stage, it might still be a partial strategy
    if boc_deprotection_found and boc_deprotection_depth <= late_stage_threshold:
        print("Partial strategy detected: Late-stage Boc deprotection")
        result = True
        if {"type": "positional", "details": {"target": "Boc amine deprotection", "position": "late_stage_first_third"}} not in findings_json["structural_constraints"]:
            findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "Boc amine deprotection", "position": "late_stage_first_third"}})

    if tbutyl_ester_deprotection_found and tbutyl_ester_deprotection_depth <= late_stage_threshold:
        print("Partial strategy detected: Late-stage t-butyl ester deprotection")
        result = True
        if {"type": "positional", "details": {"target": "t-butyl ester deprotection", "position": "late_stage_first_third"}} not in findings_json["structural_constraints"]:
            findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "t-butyl ester deprotection", "position": "late_stage_first_third"}})

    return result, findings_json

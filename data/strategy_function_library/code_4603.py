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
    Detects a protection-deprotection sequence with late-stage carbamate formation.
    Specifically looks for:
    1. Amine deprotection (from sulfonamide, Boc, etc.)
    2. Late-stage carbamate formation as one of the final steps
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

    # Initialize tracking variables
    deprotected_amine_smiles = None
    deprotection_depth = float("inf")
    has_late_carbamate_formation = False
    carbamate_reactant_smiles = None
    carbamate_formation_depth = float("inf")
    result = False # Initialize the main boolean result

    def dfs_traverse(node, depth=0):
        nonlocal deprotected_amine_smiles, has_late_carbamate_formation, carbamate_reactant_smiles
        nonlocal deprotection_depth, carbamate_formation_depth, findings_json

        if node["type"] == "reaction":
            try:
                # Extract reactants and product
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check for carbamate formation in the late stages (depth 0 or 1)
                if depth <= 1:
                    is_carbamate_formation = False

                    if checker.check_reaction(
                        "Boc amine protection", rsmi
                    ):
                        is_carbamate_formation = True
                        if "Boc amine protection" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("Boc amine protection")

                    if checker.check_reaction("Carbamic ester", rsmi):
                        is_carbamate_formation = True
                        if "Carbamic ester" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("Carbamic ester")

                    # This is a more general check for carbamate formation
                    if checker.check_fg("Carbamic ester", product_smiles) and not any(
                        checker.check_fg("Carbamic ester", r) for r in reactants_smiles
                    ):
                        is_carbamate_formation = True
                        if "Carbamic ester" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Carbamic ester")
                        if "carbamate_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("carbamate_formation")

                    if is_carbamate_formation:
                        has_late_carbamate_formation = True
                        carbamate_formation_depth = depth
                        # Add positional constraint if met
                        if {"type": "positional", "details": {"target": "carbamate_formation", "position": "last_n_stages", "value": 2}} not in findings_json["structural_constraints"]:
                            findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "carbamate_formation", "position": "last_n_stages", "value": 2}})

                        for reactant in reactants_smiles:
                            if checker.check_fg("Primary amine", reactant):
                                carbamate_reactant_smiles = reactant
                                if "Primary amine" not in findings_json["atomic_checks"]["functional_groups"]:
                                    findings_json["atomic_checks"]["functional_groups"].append("Primary amine")
                                break
                            if checker.check_fg("Secondary amine", reactant):
                                carbamate_reactant_smiles = reactant
                                if "Secondary amine" not in findings_json["atomic_checks"]["functional_groups"]:
                                    findings_json["atomic_checks"]["functional_groups"].append("Secondary amine")
                                break

                # Check for amine deprotection at any depth
                is_deprotection = False

                if checker.check_reaction("Boc amine deprotection", rsmi):
                    is_deprotection = True
                    if "Boc amine deprotection" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("Boc amine deprotection")

                if checker.check_reaction("Sulfonamide deprotection", rsmi):
                    is_deprotection = True
                    if "Sulfonamide deprotection" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("Sulfonamide deprotection")

                # General deprotection check
                if (
                    (any(checker.check_fg("Sulfonamide", reactant) for reactant in reactants_smiles) or
                     any(checker.check_fg("Boc", reactant) for reactant in reactants_smiles))
                    and (
                        checker.check_fg("Primary amine", product_smiles) or
                        checker.check_fg("Secondary amine", product_smiles)
                    )
                    and not any(checker.check_fg("Sulfonamide", product_smiles) for _ in [0]) # Ensure not present in product
                    and not any(checker.check_fg("Boc", product_smiles) for _ in [0]) # Ensure not present in product
                ):
                    is_deprotection = True
                    if any(checker.check_fg("Sulfonamide", reactant) for reactant in reactants_smiles):
                        if "Sulfonamide" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Sulfonamide")
                    if any(checker.check_fg("Boc", reactant) for reactant in reactants_smiles):
                        if "Boc" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Boc")
                    if checker.check_fg("Primary amine", product_smiles):
                        if "Primary amine" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Primary amine")
                    if checker.check_fg("Secondary amine", product_smiles):
                        if "Secondary amine" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Secondary amine")
                    if "amine_deprotection" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("amine_deprotection")

                if is_deprotection:
                    deprotected_amine_smiles = product_smiles
                    deprotection_depth = depth

            except Exception as e:
                pass

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

        return False

    # Start traversal from root
    dfs_traverse(route)

    # Check if the deprotected amine is used in the carbamate formation
    if deprotected_amine_smiles and carbamate_reactant_smiles:
        # Canonicalize SMILES for comparison
        try:
            deprotected_mol = Chem.MolFromSmiles(deprotected_amine_smiles)
            carbamate_reactant_mol = Chem.MolFromSmiles(carbamate_reactant_smiles)

            if deprotected_mol and carbamate_reactant_mol:
                deprotected_canonical = Chem.MolToSmiles(deprotected_mol)
                carbamate_reactant_canonical = Chem.MolToSmiles(carbamate_reactant_mol)

                # Check if the deprotected amine is used in carbamate formation
                # and that deprotection happens before carbamate formation
                if (
                    deprotected_canonical == carbamate_reactant_canonical
                    and deprotection_depth > carbamate_formation_depth
                ):
                    result = True
                    # Add sequence constraint if met
                    if {"type": "sequence", "details": {"before": "amine_deprotection", "after": "carbamate_formation", "constraint": "product_is_reactant"}} not in findings_json["structural_constraints"]:
                        findings_json["structural_constraints"].append({"type": "sequence", "details": {"before": "amine_deprotection", "after": "carbamate_formation", "constraint": "product_is_reactant"}})

        except Exception as e:
            pass

    return result, findings_json

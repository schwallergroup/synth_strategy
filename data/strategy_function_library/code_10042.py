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
    Detects if the synthetic route involves opening of a cyclic anhydride
    to form amide bonds with an amine.
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

    has_anhydride_opening = False

    def dfs_traverse(node, depth=0):
        nonlocal has_anhydride_opening, findings_json

        if node["type"] == "reaction":
            try:
                rsmi = node.get("metadata", {}).get("rsmi", "")
                if not rsmi:
                    return

                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                has_anhydride = False
                has_amine = False
                anhydride_reactant = None
                amine_reactant = None

                for reactant_smiles in reactants_smiles:
                    # Check for anhydride
                    if checker.check_fg("Anhydride", reactant_smiles):
                        # Verify it's a cyclic anhydride by checking for rings
                        reactant_mol = Chem.MolFromSmiles(reactant_smiles)
                        if reactant_mol and reactant_mol.GetRingInfo().NumRings() > 0:
                            has_anhydride = True
                            anhydride_reactant = reactant_smiles
                            if "Anhydride" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Anhydride")

                    # Check for primary or secondary amine
                    if checker.check_fg("Primary amine", reactant_smiles):
                        has_amine = True
                        amine_reactant = reactant_smiles
                        if "Primary amine" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Primary amine")
                    elif checker.check_fg("Secondary amine", reactant_smiles):
                        has_amine = True
                        amine_reactant = reactant_smiles
                        if "Secondary amine" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Secondary amine")

                # Check if product has amide bonds (products of anhydride opening)
                if has_anhydride and has_amine:
                    # Check for amide formation in product
                    has_primary_amide = checker.check_fg("Primary amide", product_smiles)
                    has_secondary_amide = checker.check_fg("Secondary amide", product_smiles)
                    has_tertiary_amide = checker.check_fg("Tertiary amide", product_smiles)
                    has_amide = has_primary_amide or has_secondary_amide or has_tertiary_amide

                    if has_primary_amide and "Primary amide" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Primary amide")
                    if has_secondary_amide and "Secondary amide" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Secondary amide")
                    if has_tertiary_amide and "Tertiary amide" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Tertiary amide")

                    # Check if the product has a carboxylic acid (another product of anhydride opening)
                    has_carboxylic_acid = checker.check_fg("Carboxylic acid", product_smiles)
                    if has_carboxylic_acid and "Carboxylic acid" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Carboxylic acid")

                    if has_amide and has_carboxylic_acid:
                        has_anhydride_opening = True
                        # Add the structural constraint if all conditions are met
                        if {
                            "type": "co-occurrence",
                            "details": {
                                "targets": [
                                    "Anhydride",
                                    "Primary amine",
                                    "Secondary amine",
                                    "Primary amide",
                                    "Secondary amide",
                                    "Tertiary amide",
                                    "Carboxylic acid"
                                ],
                                "notes": "The strategy requires a single reaction step where reactants contain a cyclic anhydride and an amine (primary or secondary), and the products contain an amide and a carboxylic acid, indicating a ring-opening amidation."
                            }
                        } not in findings_json["structural_constraints"]:
                            findings_json["structural_constraints"].append({
                                "type": "co-occurrence",
                                "details": {
                                    "targets": [
                                        "Anhydride",
                                        "Primary amine",
                                        "Secondary amine",
                                        "Primary amide",
                                        "Secondary amide",
                                        "Tertiary amide",
                                        "Carboxylic acid"
                                    ],
                                    "notes": "The strategy requires a single reaction step where reactants contain a cyclic anhydride and an amine (primary or secondary), and the products contain an amide and a carboxylic acid, indicating a ring-opening amidation."
                                }
                            })
            except Exception as e:
                pass

        # Continue traversing
        if has_anhydride_opening:
            return
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else:
                dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    return has_anhydride_opening, findings_json
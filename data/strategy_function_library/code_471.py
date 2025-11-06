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


PROTECTION_REACTIONS = [
    "Boc amine protection",
    "Boc amine protection explicit",
    "Boc amine protection with Boc anhydride",
    "Boc amine protection (ethyl Boc)",
    "Boc amine protection of secondary amine",
    "Boc amine protection of primary amine",
    "Alcohol protection with silyl ethers",
    "Protection of carboxylic acid",
    "Aldehyde or ketone acetalization",
    "Diol acetalization",
]

DEPROTECTION_REACTIONS = [
    "Boc amine deprotection",
    "Boc amine deprotection of guanidine",
    "Boc amine deprotection to NH-NH2",
    "Alcohol deprotection from silyl ethers",
    "Alcohol deprotection from silyl ethers (double)",
    "Alcohol deprotection from silyl ethers (diol)",
    "Cleavage of methoxy ethers to alcohols",
    "Cleavage of alkoxy ethers to alcohols",
    "Ether cleavage to primary alcohol",
    "COOH ethyl deprotection",
    "Hydroxyl benzyl deprotection",
    "Carboxyl benzyl deprotection",
    "Tert-butyl deprotection of amine",
    "N-glutarimide deprotection",
    "Phthalimide deprotection",
    "TMS deprotection from alkyne",
    "Acetal hydrolysis to diol",
    "Acetal hydrolysis to aldehyde",
    "Ketal hydrolysis to ketone",
    "Deprotection of carboxylic acid",
    "Ester saponification (methyl deprotection)",
    "Ester saponification (alkyl deprotection)",
    "Hydrogenolysis of amides/imides/carbamates",
    "Hydrolysis of amides/imides/carbamates",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if a route involves both protection and deprotection steps. It first checks for named reactions defined in the module-level PROTECTION_REACTIONS and DEPROTECTION_REACTIONS lists. It then uses functional group analysis as a fallback to identify common protection schemes (carbamate, acetal/ketal formation) and deprotection schemes (acetal/ketal hydrolysis, ester-to-acid hydrolysis, amide/carbamate-to-amine hydrolysis).
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

    protection_found = False
    deprotection_found = False

    def dfs_traverse(node, depth=0):
        nonlocal protection_found, deprotection_found, findings_json

        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for protection reactions
                for rxn_name in PROTECTION_REACTIONS:
                    if checker.check_reaction(rxn_name, rsmi):
                        print(f"Protection reaction detected: {rxn_name}")
                        protection_found = True
                        if rxn_name not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(rxn_name)
                        break

                # Check for deprotection reactions
                for rxn_name in DEPROTECTION_REACTIONS:
                    if checker.check_reaction(rxn_name, rsmi):
                        print(f"Deprotection reaction detected: {rxn_name}")
                        deprotection_found = True
                        if rxn_name not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(rxn_name)
                        break

                # Check for Carbamic ester protection (Cbz is a type of carbamic ester)
                has_carbamic_ester_in_product = checker.check_fg("Carbamic ester", product)

                if has_carbamic_ester_in_product:
                    has_carbamic_ester_in_reactants = any(
                        checker.check_fg("Carbamic ester", r) for r in reactants
                    )

                    if not has_carbamic_ester_in_reactants:
                        print("Carbamic ester protection detected through FG analysis")
                        protection_found = True
                        if "Carbamic ester" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Carbamic ester")

                # Check for acetal/ketal protection
                has_acetal_in_product = checker.check_fg("Acetal/Ketal", product)
                if has_acetal_in_product:
                    has_acetal_in_reactants = any(
                        checker.check_fg("Acetal/Ketal", r) for r in reactants
                    )
                    has_aldehyde_or_ketone_in_reactants = any(
                        checker.check_fg("Aldehyde", r) or checker.check_fg("Ketone", r)
                        for r in reactants
                    )

                    if not has_acetal_in_reactants and has_aldehyde_or_ketone_in_reactants:
                        print("Acetal/ketal protection detected through FG analysis")
                        protection_found = True
                        if "Acetal/Ketal" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Acetal/Ketal")
                        if "Aldehyde" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Aldehyde")
                        if "Ketone" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Ketone")

                # Check for acetal/ketal deprotection
                product_has_carbonyl = checker.check_fg(
                    "Aldehyde", product
                ) or checker.check_fg("Ketone", product)
                if product_has_carbonyl:
                    reactants_have_acetal = any(
                        checker.check_fg("Acetal/Ketal", r) for r in reactants
                    )
                    if reactants_have_acetal and not any(
                        checker.check_fg("Aldehyde", r) or checker.check_fg("Ketone", r)
                        for r in reactants
                    ):
                        print("Acetal/ketal deprotection detected through FG analysis")
                        deprotection_found = True
                        if "Acetal/Ketal" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Acetal/Ketal")
                        if "Aldehyde" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Aldehyde")
                        if "Ketone" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Ketone")

                # Check for carboxylic acid deprotection from esters
                product_has_carboxylic_acid = checker.check_fg("Carboxylic acid", product)
                if product_has_carboxylic_acid:
                    reactants_have_ester = any(checker.check_fg("Ester", r) for r in reactants)
                    if reactants_have_ester and not any(
                        checker.check_fg("Carboxylic acid", r) for r in reactants
                    ):
                        print(
                            "Ester deprotection to carboxylic acid detected through FG analysis"
                        )
                        deprotection_found = True
                        if "Carboxylic acid" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Carboxylic acid")
                        if "Ester" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Ester")

                # Check for amine deprotection from carbamates or amides
                product_has_amine = (
                    checker.check_fg("Primary amine", product)
                    or checker.check_fg("Secondary amine", product)
                    or checker.check_fg("Tertiary amine", product)
                )
                if product_has_amine:
                    reactants_have_protected_amine = any(
                        checker.check_fg("Carbamic ester", r)
                        or checker.check_fg("Primary amide", r)
                        or checker.check_fg("Secondary amide", r)
                        or checker.check_fg("Tertiary amide", r)
                        or checker.check_fg("Boc", r)
                        for r in reactants
                    )
                    if reactants_have_protected_amine and not any(
                        checker.check_fg("Primary amine", r)
                        or checker.check_fg("Secondary amine", r)
                        or checker.check_fg("Tertiary amine", r)
                        for r in reactants
                    ):
                        print("Amine deprotection detected through FG analysis")
                        deprotection_found = True
                        if "Primary amine" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Primary amine")
                        if "Secondary amine" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Secondary amine")
                        if "Tertiary amine" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Tertiary amine")
                        if "Carbamic ester" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Carbamic ester")
                        if "Primary amide" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Primary amide")
                        if "Secondary amide" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Secondary amide")
                        if "Tertiary amide" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Tertiary amide")
                        if "Boc" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Boc")

            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Determine the next depth based on the current node's type
        next_depth = depth
        if node["type"] != "reaction":  # If current node is 'chemical'
            next_depth = depth + 1

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, next_depth)

    dfs_traverse(route)

    print(f"Protection found: {protection_found}, Deprotection found: {deprotection_found}")
    # Return True if both protection and deprotection are found
    result = protection_found and deprotection_found

    if result:
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "protection_found",
                    "deprotection_found"
                ],
                "description": "The route must contain at least one protection event and at least one deprotection event. These events can be identified either by name (e.g., 'Boc amine protection') or by functional group transformation (e.g., formation of a carbamic ester from an amine)."
            }
        })

    return result, findings_json

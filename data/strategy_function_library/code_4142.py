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


BOC_PROTECTION_REACTIONS = [
    "Boc amine protection",
    "Boc amine protection explicit",
    "Boc amine protection with Boc anhydride",
    "Boc amine protection (ethyl Boc)",
    "Boc amine protection of secondary amine",
    "Boc amine protection of primary amine",
]

BOC_DEPROTECTION_REACTIONS = [
    "Boc amine deprotection",
    "Boc amine deprotection of guanidine",
    "Boc amine deprotection to NH-NH2",
    "Tert-butyl deprotection of amine",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a full Boc-protection and deprotection sequence for amines within a synthetic route. 
    It identifies protection/deprotection steps by checking against a defined list of named reactions 
    or by observing the formation/removal of the Boc group on an amine-containing molecule.
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

    # Track molecules that get protected and deprotected
    protected_molecules = set()
    deprotected_molecules = set()

    # Track all molecules with Boc groups to catch protection that might be missed
    molecules_with_boc = set()

    def dfs_traverse(node, depth=0):
        nonlocal protected_molecules, deprotected_molecules, molecules_with_boc, findings_json

        if node["type"] == "mol" and "smiles" in node:
            # Track molecules with Boc groups
            if checker.check_fg("Boc", node["smiles"]):
                molecules_with_boc.add(node["smiles"])
                if "Boc" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Boc")

        elif node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]

            try:
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for amine protection reactions
                for r_name in BOC_PROTECTION_REACTIONS:
                    if checker.check_reaction(r_name, rsmi):
                        protected_molecules.add(product)
                        if r_name not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(r_name)
                        break # Found a match, no need to check other protection reactions

                # Manual check for protection
                elif checker.check_fg("Boc", product) and not any(
                    checker.check_fg("Boc", r) for r in reactants
                ):
                    if "Boc" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Boc")

                    has_amine = False
                    for r in reactants:
                        if checker.check_fg("Primary amine", r):
                            if "Primary amine" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Primary amine")
                            has_amine = True
                            break
                        if checker.check_fg("Secondary amine", r):
                            if "Secondary amine" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Secondary amine")
                            has_amine = True
                            break

                    if has_amine:
                        protected_molecules.add(product)

                # Check for amine deprotection reactions
                for r_name in BOC_DEPROTECTION_REACTIONS:
                    if checker.check_reaction(r_name, rsmi):
                        deprotected_molecules.add(product)
                        if r_name not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(r_name)
                        break # Found a match, no need to check other deprotection reactions

                # Manual check for deprotection
                elif any(checker.check_fg("Boc", r) for r in reactants) and not checker.check_fg(
                    "Boc", product
                ):
                    if "Boc" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Boc")

                    has_amine = False
                    if checker.check_fg("Primary amine", product):
                        if "Primary amine" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Primary amine")
                        has_amine = True
                    if checker.check_fg("Secondary amine", product):
                        if "Secondary amine" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Secondary amine")
                        has_amine = True

                    if has_amine:
                        deprotected_molecules.add(product)

            except Exception:
                # Silently ignore errors in reaction processing
                pass

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                # Depth remains the same when traversing from reaction to chemical
                dfs_traverse(child, depth)
            else:
                # Depth increases when traversing from chemical to reaction
                dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if we found both protection and deprotection
    has_protection = len(protected_molecules) > 0
    has_deprotection = len(deprotected_molecules) > 0

    # If we didn't find explicit protection but found molecules with Boc and deprotection,
    # we can infer that protection happened somewhere
    if not has_protection and len(molecules_with_boc) > 0 and has_deprotection:
        has_protection = True

    result = has_protection and has_deprotection

    if result:
        # Add the structural constraint if both protection and deprotection are found
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "Boc_protection_of_amine",
                    "Boc_deprotection_of_amine"
                ],
                "description": "The route must contain at least one Boc protection event and at least one Boc deprotection event. An event can be a named reaction or inferred from the appearance/disappearance of a Boc group on an amine."
            }
        })

    return result, findings_json

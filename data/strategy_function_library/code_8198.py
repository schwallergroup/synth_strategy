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
    This function detects a silicon-based protection-deprotection sequence
    for alcohols (specifically TBDPS protection).
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

    # Track protected molecules and their protection/deprotection status
    protected_molecules = {}

    def dfs_traverse(node, depth=0):
        nonlocal findings_json
        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for protection: R-OH + silyl-Cl → R-O-silyl
            is_protection_reaction = checker.check_reaction(
                "Alcohol protection with silyl ethers", rsmi
            )

            if is_protection_reaction:
                findings_json["atomic_checks"]["named_reactions"].append("Alcohol protection with silyl ethers")
                # The product is the silyl-protected molecule.
                # We use its SMILES as the key to track it.
                protected_molecules[product] = {"protected": True, "depth": depth}

            # Check for deprotection: R-O-silyl → R-OH
            is_deprotection_reaction_1 = checker.check_reaction("Alcohol deprotection from silyl ethers", rsmi)
            is_deprotection_reaction_2 = checker.check_reaction("Alcohol deprotection from silyl ethers (double)", rsmi)
            is_deprotection_reaction_3 = checker.check_reaction("Alcohol deprotection from silyl ethers (diol)", rsmi)

            if is_deprotection_reaction_1:
                findings_json["atomic_checks"]["named_reactions"].append("Alcohol deprotection from silyl ethers")
            if is_deprotection_reaction_2:
                findings_json["atomic_checks"]["named_reactions"].append("Alcohol deprotection from silyl ethers (double)")
            if is_deprotection_reaction_3:
                findings_json["atomic_checks"]["named_reactions"].append("Alcohol deprotection from silyl ethers (diol)")

            if is_deprotection_reaction_1 or is_deprotection_reaction_2 or is_deprotection_reaction_3:
                # Find which reactant has the silyl group
                for r in reactants:
                    # Ensure we are deprotecting a molecule we have previously tracked as protected.
                    if r in protected_molecules and (
                        checker.check_fg("Silyl protective group", r)
                        or checker.check_fg("TMS ether protective group", r)
                    ):
                        if checker.check_fg("Silyl protective group", r):
                            findings_json["atomic_checks"]["functional_groups"].append("Silyl protective group")
                        if checker.check_fg("TMS ether protective group", r):
                            findings_json["atomic_checks"]["functional_groups"].append("TMS ether protective group")

                        protected_molecules[r]["deprotected"] = True
                        protected_molecules[r]["deprotection_depth"] = depth
                        break

        # Traverse children
        for child in node.get("children", []):
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    # Check if we have at least one molecule that was both protected and deprotected
    protection_deprotection_found = False
    for mol, status in protected_molecules.items():
        if status.get("protected") and status.get("deprotected"):
            protection_deprotection_found = True
            break

    if protection_deprotection_found:
        # Add the structural constraint if the condition is met
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "description": "The synthesis must contain both a protection reaction and a subsequent deprotection reaction for the same silyl-protected alcohol intermediate. The deprotection can be one of several related reaction types.",
                "targets": [
                    {
                        "event": "Alcohol protection with silyl ethers",
                        "required": True
                    },
                    {
                        "any_of": [
                            "Alcohol deprotection from silyl ethers",
                            "Alcohol deprotection from silyl ethers (double)",
                            "Alcohol deprotection from silyl ethers (diol)"
                        ],
                        "required": True
                    }
                ]
            }
        })

    return protection_deprotection_found, findings_json

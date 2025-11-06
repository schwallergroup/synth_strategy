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
    Detects if the synthesis involves a piperazine scaffold with a late-stage Boc protection/deprotection strategy.
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

    # Track if we found the pattern
    found_pattern = False

    # Track the depth at which we find key transformations
    boc_protection_present = False
    deprotection_depth = None
    piperazine_present = False

    # Track molecules for analysis
    protected_molecules = []
    deprotected_molecules = []

    def dfs_traverse(node, depth=0):
        nonlocal found_pattern, boc_protection_present, deprotection_depth, piperazine_present
        nonlocal protected_molecules, deprotected_molecules, findings_json

        if node["type"] == "mol":
            mol_smiles = node["smiles"]

            # Check for piperazine scaffold
            if checker.check_ring("piperazine", mol_smiles):
                piperazine_present = True
                if "piperazine" not in findings_json["atomic_checks"]["ring_systems"]:
                    findings_json["atomic_checks"]["ring_systems"].append("piperazine")
                print(f"Found piperazine scaffold in molecule: {mol_smiles}")

                # Check for Boc protection
                if checker.check_fg("Boc", mol_smiles):
                    boc_protection_present = True
                    protected_molecules.append(mol_smiles)
                    if "Boc" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Boc")
                    print(f"Found Boc protection on piperazine in molecule: {mol_smiles}")

        elif node["type"] == "reaction":
            # Check for deprotection reactions
            if "mapped_reaction_smiles" in node["metadata"]:
                rxn_smiles = node["metadata"]["mapped_reaction_smiles"]

                # Extract reactants and product
                reactants = rxn_smiles.split(">")[0].split(".")
                product = rxn_smiles.split(">")[-1]

                # Check if piperazine is involved in the reaction
                reactant_has_piperazine = any(
                    checker.check_ring("piperazine", r) for r in reactants
                )
                product_has_piperazine = checker.check_ring("piperazine", product)

                if reactant_has_piperazine and product_has_piperazine:
                    # Check for Boc deprotection
                    boc_deprotection_reactions = [
                        "Boc amine deprotection",
                        "Boc amine deprotection of guanidine",
                        "Boc amine deprotection to NH-NH2"
                    ]
                    for r_name in boc_deprotection_reactions:
                        if checker.check_reaction(r_name, rxn_smiles):
                            deprotection_depth = depth
                            deprotected_molecules.append(product)
                            if r_name not in findings_json["atomic_checks"]["named_reactions"]:
                                findings_json["atomic_checks"]["named_reactions"].append(r_name)
                            print(f"Found Boc deprotection of piperazine at depth {depth}")
                            break # Found one, no need to check others

        # Traverse children
        for child in node.get("children", []):
            # New logic: depth increases only from chemical to reaction node
            new_depth = depth
            if node["type"] != "reaction": # If current node is 'mol' (chemical)
                new_depth = depth + 1
            
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    # Check if we found the pattern
    if (
        piperazine_present
        and boc_protection_present
        and deprotection_depth is not None
    ):
        # Deprotection should happen at a low depth (late in synthesis)
        if deprotection_depth <= 3:  # Allow slightly deeper reactions as "late-stage"
            found_pattern = True
            # Add structural constraints if the pattern is found
            findings_json["structural_constraints"].append({
                "type": "co-occurrence",
                "details": {
                    "targets": [
                        "piperazine",
                        "Boc",
                        [
                            "Boc amine deprotection",
                            "Boc amine deprotection of guanidine",
                            "Boc amine deprotection to NH-NH2"
                        ]
                    ]
                }
            })
            findings_json["structural_constraints"].append({
                "type": "positional",
                "details": {
                    "target": [
                        "Boc amine deprotection",
                        "Boc amine deprotection of guanidine",
                        "Boc amine deprotection to NH-NH2"
                    ],
                    "position": "late_stage",
                    "condition": {
                        "operator": "<=",
                        "value": 3
                    }
                }
            })
            print("Found piperazine scaffold with protection/deprotection strategy")
            print(
                f"Protection type: {'Boc' if boc_protection_present else ''}"
            )
            print(f"Deprotection depth: {deprotection_depth}")

    return found_pattern, findings_json

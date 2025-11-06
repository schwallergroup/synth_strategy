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


# Refactoring for Enumeration: Isolate the lists of chemical entities.
N_HETEROCYCLES_OF_INTEREST = [
    "piperazine", "piperidine", "pyrrolidine", "morpholine",
    "pyridine", "pyrrole", "imidazole", "triazole", "tetrazole"
]

BOC_PROTECTION_REACTIONS = [
    "Boc amine protection", "Boc amine protection explicit",
    "Boc amine protection with Boc anhydride", "Boc amine protection (ethyl Boc)",
    "Boc amine protection of secondary amine", "Boc amine protection of primary amine"
]

BOC_DEPROTECTION_REACTIONS = [
    "Boc amine deprotection", "Boc amine deprotection of guanidine",
    "Boc amine deprotection to NH-NH2", "Tert-butyl deprotection of amine"
]

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects if the synthetic route uses Boc protection/deprotection strategy,
    particularly for nitrogen-containing heterocycles like piperazine.
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

    has_boc_protected = False
    has_deprotection = False
    has_n_heterocycle = False
    deprotection_depth = None

    def dfs_traverse(node, depth=0):
        nonlocal has_boc_protected, has_deprotection, has_n_heterocycle, deprotection_depth, findings_json

        # Check for nitrogen heterocycles in molecule nodes
        if node["type"] == "mol":
            mol_smiles = node["smiles"]
            try:
                # Check for nitrogen heterocycles
                for ring in N_HETEROCYCLES_OF_INTEREST:
                    if checker.check_ring(ring, mol_smiles):
                        has_n_heterocycle = True
                        if ring not in findings_json["atomic_checks"]["ring_systems"]:
                            findings_json["atomic_checks"]["ring_systems"].append(ring)

            except Exception as e:
                print(f"Error checking molecule {mol_smiles}: {e}")

        # Check for reactions
        elif node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]

                # Check for Boc protection reaction
                for rxn in BOC_PROTECTION_REACTIONS:
                    if checker.check_reaction(rxn, rsmi):
                        has_boc_protected = True
                        if rxn not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(rxn)
                        print(f"Found Boc protection reaction: {rsmi}")
                        break # Only need to find one

                # Check for Boc deprotection reaction
                for rxn in BOC_DEPROTECTION_REACTIONS:
                    if checker.check_reaction(rxn, rsmi):
                        has_deprotection = True
                        deprotection_depth = depth
                        if rxn not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(rxn)
                        print(f"Found Boc deprotection reaction at depth {depth}: {rsmi}")
                        break # Only need to find one

            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Traverse children
        for child in node.get("children", []):
            # New logic: depth increases only when traversing from chemical to reaction node
            # Depth remains the same when traversing from reaction to chemical node
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else: # node["type"] == "mol"
                dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if it's a late-stage deprotection (depth 0 or 1)
    late_stage_deprotection = (
        has_deprotection and deprotection_depth is not None and deprotection_depth <= 1
    )

    print(f"Boc protection detected: {has_boc_protected}")
    print(f"N-heterocycle detected: {has_n_heterocycle}")
    print(f"Boc deprotection detected: {has_deprotection}")
    print(f"Deprotection depth: {deprotection_depth}")
    print(f"Late-stage deprotection: {late_stage_deprotection}")

    # Determine the final result
    result = has_boc_protected and has_deprotection and late_stage_deprotection and has_n_heterocycle

    # Add structural constraints to findings_json if conditions are met
    if has_boc_protected and has_deprotection and has_n_heterocycle:
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "Boc protection reaction",
                    "Boc deprotection reaction",
                    "N-heterocycle"
                ]
            }
        })
    if has_deprotection and late_stage_deprotection:
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "Boc deprotection reaction",
                "position": "depth <= 1"
            }
        })

    # Return true if all conditions are met
    return result, findings_json

from typing import Tuple, Dict, List
import copy
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
    This function detects a strategy involving both ring formation and ring opening steps.
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

    ring_formation_detected = False
    ring_opening_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal ring_formation_detected, ring_opening_detected, findings_json

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check using ring count as a general method
                try:
                    reactant_mols = [Chem.MolFromSmiles(smi) for smi in reactants_smiles if smi]
                    product_mol = Chem.MolFromSmiles(product_smiles) if product_smiles else None

                    if all(reactant_mols) and product_mol:
                        reactant_ring_count = sum(
                            mol.GetRingInfo().NumRings() for mol in reactant_mols
                        )
                        product_ring_count = product_mol.GetRingInfo().NumRings()

                        if product_ring_count > reactant_ring_count:
                            ring_formation_detected = True
                            if "ring_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                                findings_json["atomic_checks"]["named_reactions"].append("ring_formation")

                        if reactant_ring_count > product_ring_count:
                            ring_opening_detected = True
                            if "ring_destruction" not in findings_json["atomic_checks"]["named_reactions"]:
                                findings_json["atomic_checks"]["named_reactions"].append("ring_destruction")
                except Exception as e:
                    # Handle any errors in RDKit operations
                    pass

                # Check for specific ring-forming reactions
                if checker.check_reaction("Diels-Alder", rsmi):
                    ring_formation_detected = True
                    if "Diels-Alder" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("Diels-Alder")
                if checker.check_reaction("Paal-Knorr pyrrole", rsmi):
                    ring_formation_detected = True
                    if "Paal-Knorr pyrrole" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("Paal-Knorr pyrrole")

                # Check for specific ring-opening reactions
                if checker.check_reaction("Ring opening of epoxide with amine", rsmi):
                    ring_opening_detected = True
                    if "Ring opening of epoxide with amine" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("Ring opening of epoxide with amine")

        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    dfs_traverse(route)

    result = ring_formation_detected and ring_opening_detected

    # Add structural constraint if both conditions are met
    if result:
        # This corresponds to the structural constraint in the metadata JSON
        # {"type": "co-occurrence", "details": {"targets": ["ring_formation", "ring_destruction"]}}
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "ring_formation",
                    "ring_destruction"
                ]
            }
        })

    # Return True if both ring formation and opening are detected
    return result, findings_json

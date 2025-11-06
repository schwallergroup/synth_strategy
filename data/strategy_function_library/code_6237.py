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
    Detects if key functional groups (CF3 and ester) are maintained throughout the synthesis.
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

    has_cf3_and_ester_in_final = False
    all_main_intermediates_have_cf3_and_ester = True

    def is_main_reactant(product_smiles, reactant_smiles):
        """Determine if a reactant is the main contributor to the product structure"""
        # Use MCS to find the maximum common substructure
        from rdkit import Chem
        from rdkit.Chem import rdFMCS

        product_mol = Chem.MolFromSmiles(product_smiles)
        reactant_mol = Chem.MolFromSmiles(reactant_smiles)

        if product_mol is None or reactant_mol is None:
            return False

        # Find the maximum common substructure
        mcs_result = rdFMCS.FindMCS(
            [product_mol, reactant_mol],
            atomCompare=rdFMCS.AtomCompare.CompareElements,
            bondCompare=rdFMCS.BondCompare.CompareOrder,
            completeRingsOnly=True,
            ringMatchesRingOnly=True,
        )

        if mcs_result.numAtoms == 0:
            return False

        # Calculate the fraction of product atoms that are in the MCS
        reactant_atom_count = reactant_mol.GetNumHeavyAtoms()
        product_atom_count = product_mol.GetNumHeavyAtoms()

        # If the reactant contributes significantly to the product structure
        # and has a substantial size relative to the product
        if reactant_atom_count > 0 and product_atom_count > 0:
            fraction = mcs_result.numAtoms / reactant_atom_count
            size_ratio = reactant_atom_count / product_atom_count
            return fraction > 0.5 and size_ratio > 0.5

        return False

    def dfs_traverse(node, depth=0, in_main_path=True):
        nonlocal has_cf3_and_ester_in_final, all_main_intermediates_have_cf3_and_ester, findings_json

        if node["type"] == "mol" and "smiles" in node:
            smiles = node["smiles"]

            # Check for CF3 and ester groups using the checker functions
            has_cf3 = checker.check_fg("Trifluoromethyl group", smiles)
            has_ester = checker.check_fg("Ester", smiles)

            if has_cf3 and "Trifluoromethyl group" not in findings_json["atomic_checks"]["functional_groups"]:
                findings_json["atomic_checks"]["functional_groups"].append("Trifluoromethyl group")
            if has_ester and "Ester" not in findings_json["atomic_checks"]["functional_groups"]:
                findings_json["atomic_checks"]["functional_groups"].append("Ester")

            if depth == 0:  # Final product
                has_cf3_and_ester_in_final = has_cf3 and has_ester
            elif in_main_path and not (has_cf3 and has_ester):
                all_main_intermediates_have_cf3_and_ester = False

        # Determine the next depth based on the current node's type
        next_depth = depth
        if node["type"] != "reaction": # If current node is 'mol', depth increases
            next_depth = depth + 1

        # For reaction nodes, determine which child is in the main synthetic path
        if (
            node["type"] == "reaction"
            and "children" in node
            and "metadata" in node
            and "mapped_reaction_smiles" in node["metadata"]
        ):
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            product = rsmi.split(">")[-1]

            # Determine main reactants for the next traversal
            main_reactant_indices = []
            for i, child in enumerate(node.get("children", [])):
                if child["type"] == "mol" and "smiles" in child:
                    if is_main_reactant(product, child["smiles"]):
                        main_reactant_indices.append(i)

            # If we couldn't determine main reactants, consider all as main
            if not main_reactant_indices:
                main_reactant_indices = list(range(len(node.get("children", []))))

            # Traverse children
            for i, child in enumerate(node.get("children", [])):
                child_in_main_path = i in main_reactant_indices
                # Depth remains the same when going from reaction to chemical
                dfs_traverse(child, depth, child_in_main_path)
        else:
            # Traverse children
            for child in node.get("children", []):
                # Depth increases when going from chemical to reaction
                dfs_traverse(child, depth + 1, in_main_path)

    dfs_traverse(route)

    result = has_cf3_and_ester_in_final and all_main_intermediates_have_cf3_and_ester

    if has_cf3_and_ester_in_final:
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "Trifluoromethyl group",
                    "Ester"
                ],
                "position": "last_stage"
            }
        })
    
    if all_main_intermediates_have_cf3_and_ester:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "main_path_intermediate_lacking_cooccurring_groups",
                "operator": "==",
                "value": 0,
                "cooccurring_groups": [
                    "Trifluoromethyl group",
                    "Ester"
                ]
            }
        })

    return result, findings_json

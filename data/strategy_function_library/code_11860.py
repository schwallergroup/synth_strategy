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


# This list is moved outside the function as per the 'Refactoring for Enumeration' rule.
N_FUNCTIONALIZATION_REACTIONS = [
    "N-methylation",
    "Methylation with MeI_primary",
    "Methylation with MeI_secondary",
    "Methylation with MeI_tertiary",
    "Eschweiler-Clarke Primary Amine Methylation",
    "Eschweiler-Clarke Secondary Amine Methylation",
    "Reductive methylation of primary amine with formaldehyde",
    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
    "Acylation of primary amines",
    "Acylation of secondary amines",
    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
    "Boc amine protection",
    "Boc amine protection explicit",
    "Boc amine protection with Boc anhydride",
    "Boc amine protection (ethyl Boc)",
    "Boc amine protection of secondary amine",
    "Boc amine protection of primary amine",
    "Alkylation of amines",
    "N-alkylation of primary amines with alkyl halides",
    "N-alkylation of secondary amines with alkyl halides",
    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
    "Sulfonamide synthesis (Schotten-Baumann) primary amine",
    "Sulfonamide synthesis (Schotten-Baumann) secondary amine",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if the same nitrogen atom is sequentially modified in a route. It identifies N-functionalization steps by checking against a defined list of reaction types (e.g., N-methylation, Boc protection, N-acylation) and confirms the sequential modification by tracking the atom-mapping number of the nitrogen atom across multiple reaction steps.
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
    
    result = False

    # Dictionary to track N-functionalization by atom mapping
    n_functionalization_by_atom = {}

    def dfs_traverse(node, depth=0):
        nonlocal result
        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants_str = rsmi.split(">")[0]
            product_str = rsmi.split(">")[-1]

            # Check if this is an N-functionalization reaction
            is_n_functionalization = False
            reaction_type = ""

            # Check against known N-functionalization reaction types
            for rxn_type in N_FUNCTIONALIZATION_REACTIONS:
                if checker.check_reaction(rxn_type, rsmi):
                    is_n_functionalization = True
                    reaction_type = rxn_type
                    if rxn_type not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append(rxn_type)
                    print(f"Found {rxn_type} at depth {depth}: {rsmi}")
                    break

            # If it's an N-functionalization reaction, track the N atoms being modified
            if is_n_functionalization:
                # Parse reactants and product to find mapped N atoms
                reactants_mol = Chem.MolFromSmiles(reactants_str)
                product_mol = Chem.MolFromSmiles(product_str)

                if reactants_mol and product_mol:
                    # Find N atoms in reactants
                    for atom in reactants_mol.GetAtoms():
                        if atom.GetSymbol() == "N" and atom.GetAtomMapNum() > 0:
                            map_num = atom.GetAtomMapNum()

                            # Check if this N atom is also in the product
                            for p_atom in product_mol.GetAtoms():
                                if p_atom.GetSymbol() == "N" and p_atom.GetAtomMapNum() == map_num:
                                    # This N atom is involved in the reaction
                                    if map_num not in n_functionalization_by_atom:
                                        n_functionalization_by_atom[map_num] = []

                                    n_functionalization_by_atom[map_num].append(
                                        (depth, reaction_type)
                                    )
                                    print(
                                        f"Tracked N atom with map number {map_num} at depth {depth}"
                                    )
                                    break

        # Continue traversing
        for child in node.get("children", []):
            new_depth = depth
            if node["type"] != "reaction": # Only increase depth if current node is not 'reaction'
                new_depth = depth + 1
            dfs_traverse(child, new_depth)

    dfs_traverse(route)

    # Check if any N atom has been functionalized sequentially
    for map_num, functionalizations in n_functionalization_by_atom.items():
        if len(functionalizations) >= 2:
            # Sort by depth
            functionalizations.sort()

            # Check if depths are consecutive or close
            for i in range(len(functionalizations) - 1):
                depth_diff = abs(functionalizations[i][0] - functionalizations[i + 1][0])
                if depth_diff <= 2:  # Allow for intermediate steps
                    print(
                        f"Found sequential N-functionalization on atom {map_num}: {functionalizations[i][1]} at depth {functionalizations[i][0]} followed by {functionalizations[i+1][1]} at depth {functionalizations[i+1][0]}"
                    )
                    result = True
                    # Add the structural constraint to findings_json
                    structural_constraint_obj = {
                        "type": "sequence",
                        "details": {
                            "description": "Finds a sequence of at least two N-functionalization reactions on the same nitrogen atom, where consecutive reactions are separated by at most one intermediate step (depth difference <= 2).",
                            "event_class": "N-functionalization",
                            "anchor_entity": "nitrogen_atom",
                            "min_events": 2,
                            "max_depth_difference": 2
                        }
                    }
                    if structural_constraint_obj not in findings_json["structural_constraints"]:
                        findings_json["structural_constraints"].append(structural_constraint_obj)
                    # Since we found one, we can break and return, or continue to find all instances.
                    # The original function returns True on first finding, so we'll stick to that.
                    return result, findings_json

    return result, findings_json
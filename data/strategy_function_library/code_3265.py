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


N_FUNCTIONALIZATION_REACTIONS = [
    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
    "Acylation of primary amines",
    "Acylation of secondary amines",
    "Acylation of secondary amines with anhydrides",
    "Sulfonamide synthesis (Schotten-Baumann) primary amine",
    "Sulfonamide synthesis (Schotten-Baumann) secondary amine",
    "N-alkylation of primary amines with alkyl halides",
    "N-alkylation of secondary amines with alkyl halides",
    "Urea synthesis via isocyanate and primary amine",
    "Urea synthesis via isocyanate and secondary amine",
    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
    "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
    "Methylation with MeI_primary",
    "Methylation with MeI_secondary",
    "Eschweiler-Clarke Primary Amine Methylation",
    "Eschweiler-Clarke Secondary Amine Methylation",
    "N-methylation",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if the synthesis involves multiple sequential nitrogen functionalization steps by checking for specific reaction types. It identifies sequences where the same nitrogen atom is modified in multiple steps. The targeted reactions include N-acylation, N-alkylation, N-arylation, sulfonamide formation, and urea synthesis, as defined in the `N_FUNCTIONALIZATION_REACTIONS` list.
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

    # Dictionary to track functionalization of specific nitrogen atoms
    n_functionalizations = {}

    def dfs_traverse(node, depth=0):
        nonlocal n_functionalizations, findings_json
        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            product = rsmi.split(">")[-1]

            # Check if this is an N-functionalization reaction
            is_n_functionalization = False
            for rxn_type in N_FUNCTIONALIZATION_REACTIONS:
                if checker.check_reaction(rxn_type, rsmi):
                    is_n_functionalization = True
                    if rxn_type not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append(rxn_type)
                    break

            if is_n_functionalization:
                # Try to identify the nitrogen atoms being functionalized using atom mapping
                try:
                    # Extract atom-mapped nitrogen atoms from the product
                    product_mol = Chem.MolFromSmiles(product)
                    if product_mol:
                        for atom in product_mol.GetAtoms():
                            if (
                                atom.GetAtomicNum() == 7 and atom.HasProp("molAtomMapNumber")
                            ):
                                map_num = atom.GetProp("molAtomMapNumber")
                                if map_num not in n_functionalizations:
                                    n_functionalizations[map_num] = 0
                                n_functionalizations[map_num] += 1
                except Exception as e:
                    print(f"Error processing atom mapping: {e}")
                    # If we can't track specific atoms, just count the reaction
                    if "generic" not in n_functionalizations:
                        n_functionalizations["generic"] = 0
                    n_functionalizations["generic"] += 1

        # Traverse children
        for child in node.get("children", []):
            # New depth calculation logic
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    # Count how many nitrogen atoms have been functionalized multiple times
    sequential_count = 0
    for map_num, count in n_functionalizations.items():
        if count >= 2:
            sequential_count += 1

    # If we couldn't track specific atoms but found multiple functionalization steps
    if "generic" in n_functionalizations and n_functionalizations["generic"] >= 2:
        sequential_count = 1  # At least one sequence detected

    print(f"N-functionalization tracking: {n_functionalizations}")
    print(f"Sequential N-functionalization count: {sequential_count}")

    # Determine the final boolean result
    result = sequential_count > 0

    # If the structural constraint is met, add it to findings_json
    if result:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "N-functionalization on a single nitrogen atom",
                "operator": ">=",
                "value": 2
            }
        })

    # Return True if we found at least one nitrogen with multiple functionalizations
    # or if we detected multiple generic N-functionalization steps
    return result, findings_json

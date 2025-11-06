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
    "Acylation of primary amines",
    "Acylation of secondary amines",
    "Sulfonamide synthesis (Schotten-Baumann) primary amine",
    "Sulfonamide synthesis (Schotten-Baumann) secondary amine",
    "Urea synthesis via isocyanate and primary amine",
    "Urea synthesis via isocyanate and secondary amine",
    "N-alkylation of primary amines with alkyl halides",
    "N-alkylation of secondary amines with alkyl halides",
    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
    "Schotten-Baumann_amide",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a linear synthesis strategy (steps are primarily bimolecular) involving two or more sequential N-functionalization reactions. The specific reactions checked are defined in the `N_FUNCTIONALIZATION_REACTIONS` list and include various acylations, sulfonations, alkylations, and arylations of nitrogen.
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

    n_functionalizations = 0
    max_reactants_per_step = 0

    def dfs_traverse(node, depth=0):
        nonlocal n_functionalizations, max_reactants_per_step, findings_json

        if node["type"] == "reaction":
            if "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")

                # Count reactants to check for linear vs convergent synthesis
                num_reactants = len([r for r in reactants if r and r.strip()])
                max_reactants_per_step = max(max_reactants_per_step, num_reactants)

                # Check for N-functionalization reactions
                for rxn_type in N_FUNCTIONALIZATION_REACTIONS:
                    if checker.check_reaction(rxn_type, rsmi):
                        n_functionalizations += 1
                        findings_json["atomic_checks"]["named_reactions"].append(rxn_type)
                        break

        # Determine the new depth for the recursive call
        new_depth = depth
        if node["type"] != "reaction": # If current node is chemical, depth increases
            new_depth = depth + 1

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    # Strategy is present if we have multiple N-functionalizations and mostly binary reactions (linear)
    result = n_functionalizations >= 2 and max_reactants_per_step <= 2

    # Record structural constraints if met
    if n_functionalizations >= 2:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "N-functionalization_reaction",
                "operator": ">=",
                "value": 2
            }
        })
    if max_reactants_per_step <= 2:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "max_reactants_per_step",
                "operator": "<=",
                "value": 2
            }
        })

    return result, findings_json

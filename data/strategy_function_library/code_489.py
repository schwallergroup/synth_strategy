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
from synth_strategy.utils.check import Check
from synth_strategy.utils import fuzzy_dict, check

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


PIPERAZINE_FUNCTIONALIZATIONS = [
    "N-alkylation of primary amines with alkyl halides",
    "N-alkylation of secondary amines with alkyl halides",
    "Acylation of primary amines",
    "Acylation of secondary amines",
    "Reductive amination with aldehyde",
    "Reductive amination with ketone",
    "Sulfonamide synthesis (Schotten-Baumann) primary amine",
    "Sulfonamide synthesis (Schotten-Baumann) secondary amine",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a synthetic strategy where a piperazine scaffold is present in the final product and is functionalized in at least two separate synthetic steps. The functionalization is identified by checking if the piperazine core is present in both a reactant and the product of a reaction that matches one of the following types: N-alkylation, acylation, reductive amination, or sulfonamide synthesis, as defined in the `PIPERAZINE_FUNCTIONALIZATIONS` list.
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

    piperazine_in_final = False
    functionalization_reactions = []

    def dfs_traverse(node, depth=0):
        nonlocal piperazine_in_final, functionalization_reactions, findings_json

        if node["type"] == "mol":
            mol_smiles = node["smiles"]
            if depth == 0:
                if checker.check_ring("piperazine", mol_smiles):
                    piperazine_in_final = True
                    if "piperazine" not in findings_json["atomic_checks"]["ring_systems"]:
                        findings_json["atomic_checks"]["ring_systems"].append("piperazine")

        elif node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            try:
                reactants_part = rsmi.split(">")[0]
                product_part = rsmi.split(">")[-1]
                reactants = reactants_part.split(".")

                if checker.check_ring("piperazine", product_part):
                    if "piperazine" not in findings_json["atomic_checks"]["ring_systems"]:
                        findings_json["atomic_checks"]["ring_systems"].append("piperazine")

                    piperazine_in_reactants = False
                    for reactant in reactants:
                        if checker.check_ring("piperazine", reactant):
                            piperazine_in_reactants = True
                            if "piperazine" not in findings_json["atomic_checks"]["ring_systems"]:
                                findings_json["atomic_checks"]["ring_systems"].append("piperazine")
                            break

                    if piperazine_in_reactants:
                        for rxn_type in PIPERAZINE_FUNCTIONALIZATIONS:
                            if checker.check_reaction(rxn_type, rsmi):
                                functionalization_reactions.append((depth, rxn_type, rsmi))
                                if rxn_type not in findings_json["atomic_checks"]["named_reactions"]:
                                    findings_json["atomic_checks"]["named_reactions"].append(rxn_type)
                                break
            except Exception:
                pass

        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    dfs_traverse(route)

    unique_depths = len(set(depth for depth, _, _ in functionalization_reactions))

    result = piperazine_in_final and unique_depths >= 2

    if piperazine_in_final:
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "piperazine",
                "position": "last_stage"
            }
        })

    if unique_depths >= 2:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "piperazine_functionalization_step",
                "operator": ">=",
                "value": 2
            }
        })

    return result, findings_json

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


HETEROCYCLES_OF_INTEREST = [
    "pyrrole",
    "pyridine",
    "pyrazole",
    "imidazole",
    "oxazole",
    "thiazole",
    "pyrimidine",
    "pyrazine",
    "triazole",
    "tetrazole",
    "indole",
    "quinoline",
    "isoquinoline",
    "benzimidazole",
    "benzoxazole",
    "benzothiazole",
    "furan",
    "thiophene",
    "oxadiazole",
    "thiadiazole",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects the late-stage (final two steps) formation of a heterocycle from a predefined list. The function identifies reactions where a specific heterocycle (e.g., pyrrole, indole, furan) is present in the product but absent from all reactants.
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

    # Track heterocycle formation
    heterocycle_formation = [False]

    def dfs(node, depth=0):
        nonlocal heterocycle_formation, findings_json
        if node["type"] == "reaction" and depth <= 2:  # Late-stage (near root)
            try:
                rxn_smiles = node["metadata"]["mapped_reaction_smiles"]
                product = rxn_smiles.split(">")[-1]
                reactants = rxn_smiles.split(">")[0].split(".")

                # Check if product contains a heterocycle that wasn't in the reactants
                heterocycle_in_product = False
                for ring_name in HETEROCYCLES_OF_INTEREST:
                    if checker.check_ring(ring_name, product):
                        heterocycle_in_product = True
                        # Check if this heterocycle was formed in this reaction
                        heterocycle_in_reactants = False
                        for reactant in reactants:
                            if checker.check_ring(ring_name, reactant):
                                heterocycle_in_reactants = True
                                break

                        if not heterocycle_in_reactants:
                            heterocycle_formation[0] = True
                            if ring_name not in findings_json["atomic_checks"]["ring_systems"]:
                                findings_json["atomic_checks"]["ring_systems"].append(ring_name)
                            # Add the structural constraint if the condition is met
                            if {"type": "positional", "details": {"target": "ring_formation", "position": "within_first_3_stages"}} not in findings_json["structural_constraints"]:
                                findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "ring_formation", "position": "within_first_3_stages"}})
                            break
            except Exception as e:
                print(f"Error in heterocycle check: {e}")

        # Continue traversal
        for child in node.get("children", []):
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs(child, new_depth)

    dfs(route)
    return heterocycle_formation[0], findings_json

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


HETEROCYCLES_OF_INTEREST = [
    "pyridine",
    "pyrimidine",
    "pyrazine",
    "pyridazine",
    "triazine",
    "furan",
    "thiophene",
    "pyrrole",
    "imidazole",
    "oxazole",
    "thiazole",
    "indole",
    "benzimidazole",
    "benzoxazole",
    "benzothiazole",
    "quinoline",
    "isoquinoline",
    "quinazoline",
    "piperidine",
    "piperazine",
    "morpholine",
    "thiomorpholine",
    "pyrrolidine",
    "oxazolidine",
    "thiazolidine",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects the de novo formation of a specific heterocycle from the predefined `HETEROCYCLES_OF_INTEREST` list.
    This strategic event must occur in the middle stages of the synthesis, defined as being between 25% and 75% of the total synthesis depth.
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

    found_heterocycle_formation = False
    max_depth = 0

    def get_max_depth(node, current_depth=0):
        nonlocal max_depth
        max_depth = max(max_depth, current_depth)

        for child in node.get("children", []):
            if node["type"] == "reaction":
                get_max_depth(child, current_depth)
            else:
                get_max_depth(child, current_depth + 1)

    get_max_depth(route)
    print(f"Maximum depth of synthesis route: {max_depth}")

    def is_heterocycle_formation(reaction_node):
        nonlocal findings_json
        try:
            rsmi = reaction_node.get("metadata", {}).get("rsmi", "")
            if not rsmi:
                return False

            reactants_part = rsmi.split(">")[0]
            product_part = rsmi.split(">")[-1]

            reactants = reactants_part.split(".")
            product = product_part

            for ring in HETEROCYCLES_OF_INTEREST:
                if checker.check_ring(ring, product):
                    reactant_has_ring = False
                    for reactant in reactants:
                        if checker.check_ring(ring, reactant):
                            reactant_has_ring = True
                            break

                    if not reactant_has_ring:
                        print(f"Found formation of {ring} heterocycle")
                        if "ring_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("ring_formation")
                        if ring not in findings_json["atomic_checks"]["ring_systems"]:
                            findings_json["atomic_checks"]["ring_systems"].append(ring)
                        return True

            return False
        except Exception as e:
            print(f"Error in is_heterocycle_formation: {e}")
            return False

    def dfs_traverse(node, depth=0):
        nonlocal found_heterocycle_formation, findings_json

        if node["type"] == "reaction":
            lower_bound = max(1, int(max_depth * 0.25))
            upper_bound = int(max_depth * 0.75)

            if lower_bound <= depth <= upper_bound and is_heterocycle_formation(node):
                print(f"Found heterocycle formation at depth {depth} (middle stage)")
                found_heterocycle_formation = True
                # Add the structural constraint if the condition is met
                if {"type": "positional", "details": {"target": "ring_formation", "position": "middle_stages", "min_percent_depth": 25, "max_percent_depth": 75}} not in findings_json["structural_constraints"]:
                    findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "ring_formation", "position": "middle_stages", "min_percent_depth": 25, "max_percent_depth": 75}})

        for child in node.get("children", []):
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else:
                dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    print(f"Heterocycle formation strategy detected: {found_heterocycle_formation}")
    return found_heterocycle_formation, findings_json

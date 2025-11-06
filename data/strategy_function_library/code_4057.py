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


HETEROCYCLE_RINGS_OF_INTEREST = [
    "furan", "pyrrole", "thiophene", "pyridine", "pyrazole", "imidazole",
    "oxazole", "thiazole", "pyrimidine", "pyrazine", "pyridazine", "triazole",
    "tetrazole", "indole", "benzimidazole", "benzoxazole", "benzothiazole",
    "quinoline", "isoquinoline", "purine", "piperidine", "piperazine",
    "morpholine", "thiomorpholine", "oxathiolane", "dioxathiolane",
    "thiazolidine", "oxazolidine", "isoxazole", "isothiazole", "oxadiazole",
    "thiadiazole",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Checks for the formation of specific heterocyclic rings via a multi-component reaction (3+ reactants).
    The list of heterocycles is defined in HETEROCYCLE_RINGS_OF_INTEREST.
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

    strategy_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal strategy_detected, findings_json

        if strategy_detected:
            return  # Early exit if strategy already detected

        if node["type"] == "reaction":
            if "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if we have 3 or more reactants (multi-component)
                if len(reactants) >= 3:
                    product_heterocycles = set()
                    for ring in HETEROCYCLE_RINGS_OF_INTEREST:
                        if checker.check_ring(ring, product):
                            product_heterocycles.add(ring)
                            if ring not in findings_json["atomic_checks"]["ring_systems"]:
                                findings_json["atomic_checks"]["ring_systems"].append(ring)

                    if product_heterocycles:
                        reactant_heterocycles = set()
                        for reactant in reactants:
                            for ring in HETEROCYCLE_RINGS_OF_INTEREST:
                                if checker.check_ring(ring, reactant):
                                    reactant_heterocycles.add(ring)

                        # Check if any new heterocycles are formed
                        new_heterocycles = product_heterocycles - reactant_heterocycles
                        if new_heterocycles:
                            strategy_detected = True
                            if "ring_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                                findings_json["atomic_checks"]["named_reactions"].append("ring_formation")
                            findings_json["structural_constraints"].append({
                                "type": "co-occurrence",
                                "details": {
                                    "targets": [
                                        "ring_formation",
                                        "multi-component reaction (>=3 reactants)"
                                    ],
                                    "notes": "The two events must occur in the same reaction step."
                                }
                            })
                            return

        # Continue traversing the tree
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    # Start traversal from the root
    dfs_traverse(route)

    return strategy_detected, findings_json

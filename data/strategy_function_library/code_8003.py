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


# Refactored list of heteroaromatic rings
HETEROAROMATIC_RINGS_OF_INTEREST = [
    "pyridine", "pyrrole", "furan", "thiophene", "imidazole", "oxazole",
    "thiazole", "pyrimidine", "pyrazine", "pyridazine", "indole",
    "benzimidazole", "benzoxazole", "benzothiazole", "quinoline",
    "isoquinoline", "purine", "triazole", "tetrazole", "isoxazole",
    "isothiazole", "oxadiazole", "thiadiazole",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Checks for a convergent reaction step where at least two distinct reactant molecules are coupled, and each of these reactants contains a heteroaromatic ring from a predefined list (e.g., pyridine, pyrrole, furan).
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

    # Track if we found the pattern
    found_convergent_step = False
    heteroaromatic_fragments = 0

    def dfs_traverse(node, depth=0):
        nonlocal found_convergent_step, heteroaromatic_fragments, findings_json

        if node["type"] == "reaction":
            if "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check for convergent step (multiple complex fragments)
                if len(reactants_smiles) >= 2:
                    # Count heteroaromatic fragments
                    heteroaromatic_count = 0
                    for reactant in reactants_smiles:
                        # Check if molecule has heteroaromatic rings
                        has_heteroaromatic_ring = False
                        for ring_name in HETEROAROMATIC_RINGS_OF_INTEREST:
                            if checker.check_ring(ring_name, reactant):
                                has_heteroaromatic_ring = True
                                if ring_name not in findings_json["atomic_checks"]["ring_systems"]:
                                    findings_json["atomic_checks"]["ring_systems"].append(ring_name)
                                break

                        if has_heteroaromatic_ring:
                            heteroaromatic_count += 1

                    if heteroaromatic_count >= 2:
                        found_convergent_step = True
                        heteroaromatic_fragments = heteroaromatic_count
                        # Add structural constraint if met
                        findings_json["structural_constraints"].append({
                            "type": "count",
                            "details": {
                                "target": "heteroaromatic_reactants_in_single_reaction",
                                "operator": ">=",
                                "value": 2
                            }
                        })

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                # If current node is a reaction, depth remains the same for children (chemicals)
                dfs_traverse(child, depth)
            else:
                # If current node is not a reaction (e.g., chemical), depth increases
                dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Return True if we found a convergent step with at least 2 heteroaromatic fragments
    result = found_convergent_step and heteroaromatic_fragments >= 2
    return result, findings_json

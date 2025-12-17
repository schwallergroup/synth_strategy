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
    "imidazole",
    "pyrrole",
    "pyrazole",
    "oxazole",
    "thiazole",
    "triazole",
    "tetrazole",
    "furan",
    "thiophene",
    "pyridine",
    "pyrimidine",
    "pyrazine",
    "pyridazine",
    "benzimidazole",
    "benzoxazole",
    "benzothiazole",
    "indole",
    "quinoline",
    "isoquinoline",
    "isoxazole",
    "isothiazole",
    "oxadiazole",
    "thiadiazole",
    "indazole",
    "benzotriazole",
]

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects if the synthesis includes a mid-stage heterocycle formation.
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

    has_heterocycle_formation = False
    total_depth = 0

    # First pass to determine total depth
    def get_max_depth(node, current_depth=0):
        nonlocal total_depth
        total_depth = max(total_depth, current_depth)

        for child in node.get("children", []):
            # The depth passed to the recursive call should be `current_depth` if the current node's type is 'reaction',
            # and `current_depth + 1` if the current node's type is not 'reaction' (e.g., 'chemical').
            new_depth = current_depth if node["type"] == "reaction" else current_depth + 1
            get_max_depth(child, new_depth)

    # Second pass to check for heterocycle formation
    def dfs_traverse(node, depth=0):
        nonlocal has_heterocycle_formation, findings_json

        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            # We consider "mid-stage" to be in the middle third of the synthesis
            is_mid_stage = depth >= total_depth // 3 and depth <= (2 * total_depth) // 3

            if is_mid_stage:
                try:
                    # Extract reactants and product
                    rsmi = node["metadata"]["mapped_reaction_smiles"]
                    reactants_smiles = rsmi.split(">")[0].split(".")
                    product_smiles = rsmi.split(">")[-1]

                    # Check for heterocycle formation
                    for heterocycle in HETEROCYCLES_OF_INTEREST:
                        # Check if heterocycle is in product
                        if checker.check_ring(heterocycle, product_smiles):
                            # Check if heterocycle wasn't present in any reactant
                            heterocycle_in_reactants = False
                            for reactant_smiles in reactants_smiles:
                                if checker.check_ring(heterocycle, reactant_smiles):
                                    heterocycle_in_reactants = True
                                    break

                            if not heterocycle_in_reactants:
                                has_heterocycle_formation = True
                                findings_json["atomic_checks"]["ring_systems"].append(heterocycle)
                                # Add the structural constraint if it's a mid-stage formation
                                findings_json["structural_constraints"].append({
                                    "type": "positional",
                                    "details": {
                                        "target": "ring_formation",
                                        "position": "middle_third"
                                    }
                                })
                                # Since we found a heterocycle formation and recorded it, we can break from this loop
                                # and potentially from the outer loop if only one such finding is needed.
                                break
                except Exception:
                    pass

        # Process children
        for child in node.get("children", []):
            # The depth passed to the recursive call should be `depth` if the current node's type is 'reaction',
            # and `depth + 1` if the current node's type is not 'reaction' (e.g., 'chemical').
            new_depth = depth if node["type"] == "reaction" else depth + 1
            dfs_traverse(child, new_depth)

    # Get total depth first
    get_max_depth(route)

    # Start traversal
    dfs_traverse(route)

    return has_heterocycle_formation, findings_json

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
    Detects a strategy using a cyclic anhydride as a building block in a convergent synthesis.
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

    # Track anhydride reactions and convergent points
    anhydride_reactions = []
    convergent_nodes = []
    max_depth = 0

    def count_branches(node):
        """Count the number of branches in the synthesis tree from this node"""
        if not node.get("children", []):
            return 1

        branch_count = 0
        for child in node.get("children", []):
            branch_count += count_branches(child)

        return branch_count

    def dfs_traverse(node, depth=0, path=None):
        nonlocal max_depth, findings_json
        if path is None:
            path = []

        max_depth = max(max_depth, depth)

        if node["type"] == "reaction":
            try:
                # Extract reactants and product
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check if this is a convergent step (multiple reactants)
                if len(reactants_smiles) > 1:
                    # Check for cyclic anhydride in reactants
                    for reactant_smiles in reactants_smiles:
                        if checker.check_fg("Anhydride", reactant_smiles):
                            # Verify it's cyclic by checking for rings
                            mol = Chem.MolFromSmiles(reactant_smiles)
                            if mol and mol.GetRingInfo().NumRings() > 0:
                                # Store information about this reaction
                                anhydride_reactions.append(
                                    {"depth": depth, "reactant": reactant_smiles, "reaction": rsmi}
                                )
                                if "Anhydride" not in findings_json["atomic_checks"]["functional_groups"]:
                                    findings_json["atomic_checks"]["functional_groups"].append("Anhydride")

                    # Count branches from this node
                    branches = count_branches(node)
                    if branches > 1:
                        convergent_nodes.append(
                            {"depth": depth, "branches": branches, "reaction": rsmi}
                        )
            except Exception as e:
                pass

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth, path + [node])

    # Start traversal
    dfs_traverse(route)

    # Strategy is present if:
    # 1. We found at least one cyclic anhydride in a reaction
    # 2. The synthesis is convergent (has multiple branches)
    # 3. The anhydride is used in early to mid stages

    has_cyclic_anhydride = len(anhydride_reactions) > 0
    is_convergent = len(convergent_nodes) > 0

    # Check if anhydrides are used in early to mid stages
    early_to_mid_stage = False
    for reaction in anhydride_reactions:
        # Consider early to mid stage if depth is in the first 2/3 of the synthesis
        if reaction["depth"] > max_depth / 3:
            early_to_mid_stage = True
            break

    strategy_present = has_cyclic_anhydride and is_convergent and early_to_mid_stage

    if has_cyclic_anhydride:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "cyclic_anhydride_reactant",
                "operator": ">",
                "value": 0
            }
        })
    if is_convergent:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "convergent_step",
                "operator": ">",
                "value": 0
            }
        })
    if early_to_mid_stage:
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "cyclic_anhydride_reactant",
                "position": "early_to_mid_stage"
            }
        })

    return strategy_present, findings_json

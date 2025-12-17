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


def main(route) -> Tuple[bool, Dict]:
    """
    This function detects if a morpholine group is introduced in the second half of the synthesis.
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

    morpholine_introduction_detected = False
    max_depth = 0

    # First pass to determine the maximum depth
    def get_max_depth(node, depth=0):
        nonlocal max_depth
        max_depth = max(max_depth, depth)
        for child in node.get("children", []):
            get_max_depth(child, depth + 1)

    get_max_depth(route)

    # Second pass to detect morpholine introduction
    def dfs_traverse(node, depth=0):
        nonlocal morpholine_introduction_detected, findings_json

        if node["type"] == "reaction" and depth < max_depth / 2:  # Second half of synthesis
            if "mapped_reaction_smiles" in node.get("metadata", {}):
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for morpholine in reactants
                morpholine_pattern = Chem.MolFromSmarts("[N]1CCO[C][C]1")
                morpholine_in_reactants = False
                for reactant in reactants:
                    if (
                        reactant
                        and Chem.MolFromSmiles(reactant)
                        and Chem.MolFromSmiles(reactant).HasSubstructMatch(morpholine_pattern)
                    ):
                        morpholine_in_reactants = True
                        break

                # Check for morpholine in product but not in all reactants
                morpholine_in_product = False
                if (
                    product
                    and Chem.MolFromSmiles(product)
                    and Chem.MolFromSmiles(product).HasSubstructMatch(morpholine_pattern)
                ):
                    morpholine_in_product = True

                if morpholine_in_product and not morpholine_in_reactants:
                    morpholine_introduction_detected = True
                    findings_json["atomic_checks"]["ring_systems"].append("morpholine")
                    findings_json["structural_constraints"].append({
                        "type": "positional",
                        "details": {
                            "target": "ring_formation",
                            "context": {
                                "ring_system": "morpholine"
                            },
                            "position": "second_half"
                        }
                    })

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else:  # chemical node
                dfs_traverse(child, depth + 1)

    # Call dfs_traverse on the root node
    dfs_traverse(route)

    return morpholine_introduction_detected, findings_json
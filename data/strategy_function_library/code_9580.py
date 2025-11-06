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
    Detects if the synthetic route involves introducing fluorinated aromatic rings
    in a late-stage diversification step.
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

    has_fluorinated_aromatic = False
    fluorinated_aromatic_depth = -1

    def has_fluorinated_aromatic_pattern(smiles):
        # Pattern for fluorinated aromatic
        pattern = Chem.MolFromSmarts("[#6]1:[#6]:[#6]:[#6](F):[#6]:[#6]:1")
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False
        return mol.HasSubstructMatch(pattern)

    def dfs_traverse(node, depth=0):
        nonlocal has_fluorinated_aromatic, fluorinated_aromatic_depth, findings_json

        if node["type"] == "reaction":
            if "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
                reaction_smiles = node["metadata"]["mapped_reaction_smiles"]
                reactants = reaction_smiles.split(">")[0]
                product = reaction_smiles.split(">")[-1]

                # Check if product has fluorinated aromatic but any reactant doesn't
                if has_fluorinated_aromatic_pattern(product):
                    reactant_list = reactants.split(".")
                    fluorinated_in_reactants = any(
                        has_fluorinated_aromatic_pattern(r) for r in reactant_list
                    )

                    # If we're introducing a fluorinated aromatic
                    if not fluorinated_in_reactants:
                        has_fluorinated_aromatic = True
                        fluorinated_aromatic_depth = depth
                        findings_json["atomic_checks"]["named_reactions"].append("introduction_of_fluorinated_aromatic_ring")
                        print(f"Found fluorinated aromatic introduction at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                # If current node is a reaction, depth remains the same for child (chemical)
                dfs_traverse(child, depth)
            else:
                # If current node is a chemical, depth increases for child (reaction)
                dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    result = False
    # The strategy is present if we found fluorinated aromatic introduction at depth 0 or 1
    if has_fluorinated_aromatic and fluorinated_aromatic_depth <= 1:
        print(
            f"Detected fluorinated aromatics diversification strategy at depth {fluorinated_aromatic_depth}"
        )
        result = True
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "introduction_of_fluorinated_aromatic_ring",
                "position": "final_two_stages"
            }
        })

    return result, findings_json
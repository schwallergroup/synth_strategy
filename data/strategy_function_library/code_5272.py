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
    Detects the construction of a specific fused heterocycle (defined by SMARTS `[#7]1[#6][#7][#6]2[#6][#6][#6][#7][#6]12`). The strategy is positive if the product contains this scaffold and none of the reactants do.
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

    imidazopyridine_pattern = Chem.MolFromSmarts("[#7]1[#6][#7][#6]2[#6][#6][#6][#7][#6]12")
    ring_formation_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal ring_formation_detected, findings_json

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Convert to RDKit molecules
            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles]
            product_mol = Chem.MolFromSmiles(product_smiles)

            # Check if product contains imidazopyridine
            if product_mol and product_mol.HasSubstructMatch(imidazopyridine_pattern):
                # Check if any reactant does NOT contain imidazopyridine
                reactants_without_scaffold = [
                    r
                    for r in reactant_mols
                    if r and not r.HasSubstructMatch(imidazopyridine_pattern)
                ]

                if len(reactants_without_scaffold) == len(reactant_mols):
                    print("Imidazopyridine scaffold construction detected")
                    ring_formation_detected = True
                    findings_json["atomic_checks"]["ring_systems"].append("[#7]1[#6][#7][#6]2[#6][#6][#6][#7][#6]12")
                    findings_json["atomic_checks"]["named_reactions"].append("ring_formation")

        # Continue traversing
        for child in node.get("children", []):
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    if ring_formation_detected:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "ring_formation",
                "operator": ">=",
                "value": 1
            }
        })

    return ring_formation_detected, findings_json

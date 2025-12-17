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


from rdkit import Chem

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if the synthesis route involves amide formation in the late stage (low depth).
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

    amide_formation_detected = False
    max_depth = 0

    def dfs_traverse(node, depth=0):
        nonlocal amide_formation_detected, max_depth, findings_json

        max_depth = max(max_depth, depth)

        if node["type"] == "reaction" and depth <= 1:  # Late stage (depth 0 or 1)
            if "mapped_reaction_smiles" in node.get("metadata", {}):
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for amide formation
                product_mol = Chem.MolFromSmiles(product)
                amide_pattern = Chem.MolFromSmarts("[#6]-[#6](=[#8])-[#7]")

                # Check if amide exists in product
                if product_mol and product_mol.HasSubstructMatch(amide_pattern):
                    findings_json["atomic_checks"]["functional_groups"].append("amide")
                    # Check if amide bond is formed in this reaction
                    amide_in_reactants = False
                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol and reactant_mol.HasSubstructMatch(amide_pattern):
                            amide_in_reactants = True
                            break

                    if not amide_in_reactants:
                        print(f"Late-stage amide formation detected at depth {depth}")
                        amide_formation_detected = True
                        findings_json["atomic_checks"]["named_reactions"].append("amide_formation")
                        findings_json["structural_constraints"].append({
                            "type": "positional",
                            "details": {
                                "target": "amide_formation",
                                "position": "late_stage"
                            }
                        })

        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is not a reaction (e.g., chemical)
                new_depth = depth + 1
            dfs_traverse(child, new_depth)

    dfs_traverse(route)
    print(f"Maximum depth in route: {max_depth}")
    return amide_formation_detected, findings_json
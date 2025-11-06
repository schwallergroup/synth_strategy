from typing import Tuple, Dict, List
import copy
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
    Detects if the synthesis route incorporates a sulfonamide group in the late stage.
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

    sulfonamide_incorporation = False
    max_depth_for_late_stage = 1

    def dfs_traverse(node, depth=0):
        nonlocal sulfonamide_incorporation, findings_json

        if node["type"] == "reaction" and depth <= max_depth_for_late_stage:
            rsmi = node.get("metadata", {}).get("rsmi", "")
            if rsmi:
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for sulfonamide pattern
                sulfonamide_pattern = Chem.MolFromSmarts("[#16](=[#8])(=[#8])[#7]")

                sulfonamide_in_reactants = False
                for reactant in reactants:
                    try:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol and mol.HasSubstructMatch(sulfonamide_pattern):
                            sulfonamide_in_reactants = True
                            findings_json["atomic_checks"]["functional_groups"].append("sulfonamide")
                            break
                    except:
                        continue

                sulfonamide_in_product = False
                try:
                    mol = Chem.MolFromSmiles(product)
                    if mol and mol.HasSubstructMatch(sulfonamide_pattern):
                        sulfonamide_in_product = True
                        findings_json["atomic_checks"]["functional_groups"].append("sulfonamide")
                except:
                    pass

                if not sulfonamide_in_reactants and sulfonamide_in_product:
                    sulfonamide_incorporation = True
                    findings_json["atomic_checks"]["named_reactions"].append("sulfonamide_formation")
                    findings_json["structural_constraints"].append({
                        "type": "positional",
                        "details": {
                            "target": "sulfonamide_formation",
                            "position": "late_stage"
                        }
                    })

        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else:
                dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return sulfonamide_incorporation, findings_json
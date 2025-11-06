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
    Detects if the synthesis involves a trifluoroethyl group as a key structural element
    throughout the synthesis.
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

    trifluoroethyl_present = False
    mol_count = 0

    def dfs_traverse(node, depth=0):
        nonlocal trifluoroethyl_present, mol_count, findings_json

        if node["type"] == "mol" and "smiles" in node:
            mol_count += 1
            try:
                mol = Chem.MolFromSmiles(node["smiles"])
                if mol:
                    trifluoroethyl_pattern = Chem.MolFromSmarts("[#6][#6]([F])([F])[F]")
                    if mol.HasSubstructMatch(trifluoroethyl_pattern):
                        trifluoroethyl_present = True
                        if "trifluoroethyl" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("trifluoroethyl")
            except:
                print("Error processing SMILES in trifluoroethyl_containing_synthesis")

        # Process children
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else:
                dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Only return true if trifluoroethyl is present in at least one molecule
    # and we've examined at least 3 molecules
    result = trifluoroethyl_present and mol_count >= 3
    if result:
        print("Found trifluoroethyl group as a key structural element")
        # Add structural constraints if the overall condition is met
        if trifluoroethyl_present:
            findings_json["structural_constraints"].append({
                "type": "count",
                "details": {
                    "target": "trifluoroethyl",
                    "operator": ">=",
                    "value": 1
                }
            })
        if mol_count >= 3:
            findings_json["structural_constraints"].append({
                "type": "count",
                "details": {
                    "target": "molecule",
                    "operator": ">=",
                    "value": 3
                }
            })

    return result, findings_json
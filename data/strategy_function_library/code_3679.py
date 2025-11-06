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
    Detects if the synthesis introduces a cyclopropyl group in a late stage.
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

    cyclopropyl_late_stage = False

    def dfs_traverse(node, depth=0):
        nonlocal cyclopropyl_late_stage, findings_json

        if node["type"] == "reaction" and depth <= 1:
            # Extract reaction SMILES
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for cyclopropyl pattern
                cyclopropyl_pattern = Chem.MolFromSmarts("[#6]1[#6][#6]1")

                # Check if cyclopropyl is in product but not in reactants
                product_mol = Chem.MolFromSmiles(product)
                if product_mol and product_mol.HasSubstructMatch(cyclopropyl_pattern):
                    findings_json["atomic_checks"]["ring_systems"].append("cyclopropane")
                    has_cyclopropyl_reactant = False

                    for reactant in reactants:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol:
                            if mol.HasSubstructMatch(cyclopropyl_pattern):
                                has_cyclopropyl_reactant = True
                                break # optimization: no need to check further

                    if not has_cyclopropyl_reactant:
                        cyclopropyl_late_stage = True
                        findings_json["structural_constraints"].append({
                            "type": "positional",
                            "details": {
                                "target": "cyclopropane_formation",
                                "position": "late_stage"
                            }
                        })

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else:
                dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)
    return cyclopropyl_late_stage, findings_json
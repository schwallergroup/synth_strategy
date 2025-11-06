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
    This function detects if the synthesis route incorporates a trifluoromethyl group
    in a late stage (depth 0 or 1).
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

    late_cf3_incorporation = False

    def dfs_traverse(node, depth=0):
        nonlocal late_cf3_incorporation, findings_json

        if node["type"] == "reaction" and depth <= 1:  # Late stage
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if product contains CF3 group
            product_mol = Chem.MolFromSmiles(product)
            cf3_pattern = Chem.MolFromSmarts("[C]([F])([F])[F]")

            if product_mol and product_mol.HasSubstructMatch(cf3_pattern):
                findings_json["atomic_checks"]["functional_groups"].append("trifluoromethyl")
                # Check if CF3 comes from one of the reactants
                cf3_in_reactants = False
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol and reactant_mol.HasSubstructMatch(cf3_pattern):
                        cf3_in_reactants = True
                        break

                if not cf3_in_reactants:
                    late_cf3_incorporation = True
                    # Add the structural constraint if the condition is met
                    findings_json["structural_constraints"].append({
                        "type": "positional",
                        "details": {
                            "target": "trifluoromethyl_group_formation",
                            "position": "late_stage"
                        }
                    })

        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases when going to reaction
                new_depth = depth + 1
            # If current node is reaction, depth remains same when going to chemical
            
            dfs_traverse(child, new_depth)

    dfs_traverse(route)
    return late_cf3_incorporation, findings_json

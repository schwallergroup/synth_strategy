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
    Detects the introduction of a trifluoromethyl (CF3) group onto a carbon atom
    in a late-stage reaction (depth <= 1).
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

    found_late_cf3_introduction = False

    def dfs_traverse(node, depth=0):
        nonlocal found_late_cf3_introduction, findings_json

        if node["type"] == "reaction" and depth <= 1:
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants_smiles = rsmi.split(">")[0]
            product_smiles = rsmi.split(">")[-1]

            reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles.split(".")]
            product = Chem.MolFromSmiles(product_smiles)

            if product is None or any(r is None for r in reactants):
                return

            # Check for CF3 group in product
            cf3_pattern = Chem.MolFromSmarts("[#6]-[C]([F])([F])[F]")
            product_has_cf3 = product.HasSubstructMatch(cf3_pattern)

            # Check if any reactant has CF3 group
            reactant_has_cf3 = any(r.HasSubstructMatch(cf3_pattern) for r in reactants)

            if product_has_cf3 and not reactant_has_cf3:
                found_late_cf3_introduction = True
                if "trifluoromethyl" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("trifluoromethyl")
                if "trifluoromethylation" not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append("trifluoromethylation")
                
                # Add the structural constraint if not already present
                constraint_obj = {
                    "type": "positional",
                    "details": {
                        "target": "trifluoromethylation",
                        "position": "late_stage (depth <= 1)"
                    }
                }
                if constraint_obj not in findings_json["structural_constraints"]:
                    findings_json["structural_constraints"].append(constraint_obj)

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else:
                dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return found_late_cf3_introduction, findings_json

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
    Detects synthesis routes where a trifluoromethylphenyl group is present in the product of every reaction step.
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

    trifluoromethyl_count = 0
    total_reactions = 0
    reactions_lacking_trifluoromethyl = 0

    def dfs_traverse(node, depth=0):
        nonlocal trifluoromethyl_count, total_reactions, reactions_lacking_trifluoromethyl, findings_json

        if node["type"] == "reaction":
            total_reactions += 1
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            product_smiles = rsmi.split(">")[-1]

            # Create RDKit molecule
            product = Chem.MolFromSmiles(product_smiles) if product_smiles else None

            if product:
                # Check for trifluoromethylphenyl group
                # The SMARTS pattern c-[#6]([#9])([#9])[#9] matches a trifluoromethyl group attached to an aromatic carbon.
                # This corresponds to 'trifluoromethyl_on_aromatic_ring' in the strategy.
                trifluoromethyl_pattern = Chem.MolFromSmarts("c-[#6]([#9])([#9])[#9]")
                if product.HasSubstructMatch(trifluoromethyl_pattern):
                    trifluoromethyl_count += 1
                    if "trifluoromethyl_on_aromatic_ring" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("trifluoromethyl_on_aromatic_ring")
                else:
                    reactions_lacking_trifluoromethyl += 1

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else: # chemical node
                dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if trifluoromethyl group is present in all reactions
    is_preserved = trifluoromethyl_count > 0 and trifluoromethyl_count == total_reactions

    # Populate structural constraints
    if total_reactions > 0:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "reaction",
                "operator": ">",
                "value": 0
            }
        })
    if reactions_lacking_trifluoromethyl == 0:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "reaction_lacking_trifluoromethyl_on_aromatic_ring_in_product",
                "operator": "==",
                "value": 0
            }
        })

    return is_preserved, findings_json

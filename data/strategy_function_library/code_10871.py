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
    Detects if the route contains and preserves a 2,4-difluorophenoxy group
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

    # Track if the motif is present in the final product
    final_product_has_motif = False
    # Track if any reaction modifies the motif
    motif_modified = False

    def is_24_difluorophenoxy(smiles):
        """Check if molecule contains a 2,4-difluorophenoxy group"""
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return False

        # Check for the specific 2,4-difluoro pattern with phenoxy
        pattern = Chem.MolFromSmarts("c1(F)cc(F)ccc1Oc")
        return mol.HasSubstructMatch(pattern)

    def dfs_traverse(node, depth=0):
        nonlocal final_product_has_motif, motif_modified, findings_json

        if node["type"] == "mol":
            # Check if this molecule has the 2,4-difluorophenoxy group
            if node["smiles"]:
                has_motif = is_24_difluorophenoxy(node["smiles"])

                # If this is the final product (depth 0), record if it has the motif
                if depth == 0:
                    final_product_has_motif = has_motif
                    if has_motif:
                        findings_json["atomic_checks"]["functional_groups"].append("2,4-difluorophenoxy")

        elif node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            # For reaction nodes, check if the motif is preserved
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if product has the motif
            product_has_motif = is_24_difluorophenoxy(product)

            # Check if any reactant has the motif
            reactant_has_motif = any(is_24_difluorophenoxy(r) for r in reactants)

            # If a reactant has the motif but the product doesn't, the motif is modified
            if reactant_has_motif and not product_has_motif:
                motif_modified = True

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                # If current node is reaction, depth remains same for child (chemical)
                dfs_traverse(child, depth)
            else:
                # If current node is chemical, depth increases for child (reaction)
                dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # The route contains and preserves the motif if it's in the final product
    # and is not modified in any reaction
    result = final_product_has_motif and not motif_modified

    if final_product_has_motif:
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "2,4-difluorophenoxy",
                "position": "last_stage"
            }
        })
    if not motif_modified:
        findings_json["structural_constraints"].append({
            "type": "negation",
            "details": {
                "target": "modification",
                "scope": "2,4-difluorophenoxy"
            }
        })

    return result, findings_json

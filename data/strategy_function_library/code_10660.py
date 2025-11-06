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
    Detects if a nitro group is introduced early and removed later,
    serving as a temporary functional group to enable other transformations.
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

    nitro_addition = False
    nitro_removal = False
    nitro_addition_depth = -1
    nitro_removal_depth = -1

    # SMARTS pattern for nitro group
    nitro_pattern = Chem.MolFromSmarts("[#8]-[#7+](=[#8])-[#6]")

    def dfs_traverse(node, depth=0):
        nonlocal nitro_addition, nitro_removal, nitro_addition_depth, nitro_removal_depth, findings_json

        if node["type"] == "reaction":
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles]
            product = Chem.MolFromSmiles(product_smiles)

            # Check for nitro addition
            if not nitro_addition and product is not None:
                if product.HasSubstructMatch(nitro_pattern):
                    has_nitro_in_reactants = False
                    for r in reactants:
                        if r is not None and r.HasSubstructMatch(nitro_pattern):
                            has_nitro_in_reactants = True
                            break

                    if not has_nitro_in_reactants:
                        nitro_addition = True
                        nitro_addition_depth = depth
                        findings_json["atomic_checks"]["named_reactions"].append("nitro_group_formation")
                        findings_json["atomic_checks"]["functional_groups"].append("nitro group")

            # Check for nitro removal
            if not nitro_removal and product is not None:
                has_nitro_in_reactants = False
                for r in reactants:
                    if r is not None and r.HasSubstructMatch(nitro_pattern):
                        has_nitro_in_reactants = True
                        break

                if has_nitro_in_reactants and not product.HasSubstructMatch(nitro_pattern):
                    nitro_removal = True
                    nitro_removal_depth = depth
                    findings_json["atomic_checks"]["named_reactions"].append("nitro_group_removal")
                    findings_json["atomic_checks"]["functional_groups"].append("nitro group")

        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                # If current node is a reaction, depth remains the same for children (chemicals)
                dfs_traverse(child, depth)
            else:
                # If current node is a chemical, depth increases for children (reactions)
                dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Return True if nitro is added early and removed later
    result = nitro_addition and nitro_removal and nitro_addition_depth > nitro_removal_depth

    if nitro_addition and nitro_removal:
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "nitro_group_formation",
                    "nitro_group_removal"
                ]
            }
        })
    if result:
        findings_json["structural_constraints"].append({
            "type": "sequence",
            "details": {
                "before": "nitro_group_formation",
                "after": "nitro_group_removal"
            }
        })

    return result, findings_json
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
    Detects if a late-stage reaction involves the reduction of an aromatic nitro group. This is checked by verifying that the first reactant molecule contains an aromatic nitro group ([c][N+](=O)[O-]), the product contains an aromatic amine ([c][NH2]), and that the number of aromatic nitro groups decreases during the reaction. Note: This function only checks the first reactant and is limited to aromatic systems.
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

    nitro_reduction_found = False
    max_depth_for_late_stage = 1

    def dfs_traverse(node, depth=0):
        nonlocal nitro_reduction_found, findings_json

        if node["type"] == "reaction" and depth <= max_depth_for_late_stage:
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if this is a nitro reduction reaction
            reactant_mol = Chem.MolFromSmiles(reactants[0])
            product_mol = Chem.MolFromSmiles(product)

            if reactant_mol and product_mol:
                nitro_pattern = Chem.MolFromSmarts("[c][N+](=O)[O-]")
                amine_pattern = Chem.MolFromSmarts("[c][NH2]")

                reactant_has_nitro = reactant_mol.HasSubstructMatch(nitro_pattern)
                product_has_amine = product_mol.HasSubstructMatch(amine_pattern)

                if reactant_has_nitro:
                    if "aromatic nitro group" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("aromatic nitro group")
                if product_has_amine:
                    if "aromatic amine" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("aromatic amine")

                if reactant_has_nitro and product_has_amine:
                    # Confirm nitro group is reduced to amine
                    nitro_matches = reactant_mol.GetSubstructMatches(nitro_pattern)
                    if not product_mol.HasSubstructMatch(nitro_pattern) or len(
                        product_mol.GetSubstructMatches(nitro_pattern)
                    ) < len(nitro_matches):
                        nitro_reduction_found = True
                        if "aromatic_nitro_reduction" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("aromatic_nitro_reduction")
                        print(f"Late-stage nitro reduction detected at depth {depth}")

        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                # If current node is a reaction, depth remains the same for children
                dfs_traverse(child, depth)
            else:
                # If current node is not a reaction (e.g., chemical), depth increases
                dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    if nitro_reduction_found:
        # Add the structural constraint if the nitro reduction was found
        findings_json["structural_constraints"].append(
            {
                "type": "positional",
                "details": {
                    "target": "aromatic_nitro_reduction",
                    "position": "late_stage"
                }
            }
        )

    return nitro_reduction_found, findings_json

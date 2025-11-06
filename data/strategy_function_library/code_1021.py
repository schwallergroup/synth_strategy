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
    This function detects if a trifluoromethyl group is introduced in the late stage of synthesis.
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

    trifluoromethyl_introduced = False

    def dfs_traverse(node, depth=0):
        nonlocal trifluoromethyl_introduced, findings_json

        if node["type"] == "reaction" and depth <= 1:  # Late stage (depth 0 or 1)
            # Extract reactants and product
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check if product has trifluoromethyl
            product_mol = Chem.MolFromSmiles(product_smiles)
            trifluoromethyl_pattern = Chem.MolFromSmarts("[C]([F])([F])[F]")

            product_has_cf3 = False
            if product_mol and product_mol.HasSubstructMatch(trifluoromethyl_pattern):
                product_has_cf3 = True
                if "trifluoromethyl" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("trifluoromethyl")

            if product_has_cf3:
                # Check if any reactant doesn't have trifluoromethyl
                reactant_without_cf3 = False
                for r_smiles in reactants_smiles:
                    r_mol = Chem.MolFromSmiles(r_smiles)
                    if r_mol and not r_mol.HasSubstructMatch(trifluoromethyl_pattern):
                        reactant_without_cf3 = True
                        break

                if reactant_without_cf3:
                    print(f"Late-stage trifluoromethyl introduction detected in reaction: {rsmi}")
                    trifluoromethyl_introduced = True
                    if "trifluoromethyl_introduction" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("trifluoromethyl_introduction")
                    
                    # Add structural constraint if not already present
                    constraint_found = False
                    for constraint in findings_json["structural_constraints"]:
                        if constraint.get("type") == "positional" and \
                           constraint.get("details", {}).get("target") == "trifluoromethyl_introduction" and \
                           constraint.get("details", {}).get("position") == "late_stage":
                            constraint_found = True
                            break
                    if not constraint_found:
                        findings_json["structural_constraints"].append({
                            "type": "positional",
                            "details": {
                                "target": "trifluoromethyl_introduction",
                                "position": "late_stage"
                            }
                        })

        # Traverse children with incremented depth
        for child in node.get("children", []):
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else:
                dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)
    return trifluoromethyl_introduced, findings_json

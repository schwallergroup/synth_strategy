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
    Detects a strategy involving the formation of a quaternary carbon center.
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

    found_quaternary_carbon_formation = False

    def dfs_traverse(node, depth=0):
        nonlocal found_quaternary_carbon_formation, findings_json

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Convert to RDKit molecules
            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles if r]
            product_mol = Chem.MolFromSmiles(product_smiles)

            if product_mol and all(r for r in reactant_mols):
                # Check for quaternary carbon in product
                quat_carbon_patt = Chem.MolFromSmarts("[#6D4](-[#6])(-[#6])(-[#6])(-[#6])")

                if product_mol.HasSubstructMatch(quat_carbon_patt):
                    findings_json["atomic_checks"]["functional_groups"].append("quaternary carbon")
                    # Check if quaternary carbon was formed in this step
                    reactants_have_quat_carbon = False
                    for r in reactant_mols:
                        if r.HasSubstructMatch(quat_carbon_patt):
                            reactants_have_quat_carbon = True
                            break

                    if not reactants_have_quat_carbon:
                        found_quaternary_carbon_formation = True
                        findings_json["atomic_checks"]["named_reactions"].append("quaternary_carbon_formation")
                        print(f"Found quaternary carbon formation at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            # New logic: depth increases only when traversing from chemical to reaction
            # Depth remains the same when traversing from reaction to chemical
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else:
                dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    if found_quaternary_carbon_formation:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "quaternary_carbon_formation",
                "operator": ">=",
                "value": 1
            }
        })

    return found_quaternary_carbon_formation, findings_json

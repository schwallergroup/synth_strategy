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
    Detects if the synthesis involves formation of an aryl-aryl amine linkage (biaryl amine).
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

    biaryl_amine_found = False

    def dfs_traverse(node, depth=0):
        nonlocal biaryl_amine_found, findings_json

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check for biaryl amine formation
                try:
                    # Create molecules from SMILES
                    product_mol = Chem.MolFromSmiles(product_smiles)

                    if product_mol:
                        # Look for aryl-NH-aryl pattern in product
                        biaryl_amine_pattern = Chem.MolFromSmarts("[c]-[NH]-[c]")
                        if product_mol.HasSubstructMatch(biaryl_amine_pattern):
                            # Check if this pattern exists in any reactant
                            pattern_in_reactants = False
                            for r_smi in reactants_smiles:
                                r_mol = Chem.MolFromSmiles(r_smi)
                                if r_mol and r_mol.HasSubstructMatch(biaryl_amine_pattern):
                                    pattern_in_reactants = True
                                    break

                            if not pattern_in_reactants:
                                biaryl_amine_found = True
                                findings_json["atomic_checks"]["functional_groups"].append("biaryl amine")
                                findings_json["atomic_checks"]["named_reactions"].append("biaryl_amine_formation")
                                findings_json["structural_constraints"].append({
                                    "type": "count",
                                    "details": {
                                        "target": "biaryl_amine_formation",
                                        "operator": ">=",
                                        "value": 1
                                    }
                                })
                except:
                    pass

        # Continue traversal
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else: # Assuming 'chemical' or other non-reaction type
                dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return biaryl_amine_found, findings_json
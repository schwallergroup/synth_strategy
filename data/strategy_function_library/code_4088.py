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
    Detects if the synthesis route involves vinyl sulfone formation as a key step.
    Specifically looks for C=C-SO2 formation in the final step.
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

    vinyl_sulfone_pattern = Chem.MolFromSmarts("[C]=[C][S](=[O])(=[O])[#6]")
    late_stage_formation = False

    def dfs_traverse(node, depth=0):
        nonlocal late_stage_formation, findings_json

        if node["type"] == "reaction" and depth <= 1:  # Late-stage reaction
            if "mapped_reaction_smiles" in node.get("metadata", {}):
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0]
                product = rsmi.split(">")[-1]

                reactant_mol = Chem.MolFromSmiles(reactants)
                product_mol = Chem.MolFromSmiles(product)

                if (
                    product_mol
                    and product_mol.HasSubstructMatch(vinyl_sulfone_pattern)
                    and (
                        not reactant_mol
                        or not reactant_mol.HasSubstructMatch(vinyl_sulfone_pattern)
                    )
                ):
                    late_stage_formation = True
                    findings_json["atomic_checks"]["functional_groups"].append("vinyl sulfone")
                    findings_json["atomic_checks"]["named_reactions"].append("vinyl_sulfone_formation")
                    findings_json["structural_constraints"].append({
                        "type": "positional",
                        "details": {
                            "target": "vinyl_sulfone_formation",
                            "position": "last_stage"
                        }
                    })

        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            
            dfs_traverse(child, new_depth)

    dfs_traverse(route)
    return late_stage_formation, findings_json
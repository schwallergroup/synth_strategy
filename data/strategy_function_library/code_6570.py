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
from steerable_retro.utils.check import Check
from steerable_retro.utils import fuzzy_dict, check

from pathlib import Path
root_data = Path(__file__).parent.parent

fg_args = {
    "file_path": f"{root_data}/patterns/functional_groups.json",
    "value_field": "pattern",
    "key_field": "name",
}
reaction_class_args = {
    "file_path": f"{root_data}/patterns/smirks.json",
    "value_field": "smirks",
    "key_field": "name",
}
ring_smiles_args = {
    "file_path": f"{root_data}/patterns/chemical_rings_smiles.json",
    "value_field": "smiles",
    "key_field": "name",
}
functional_groups = fuzzy_dict.FuzzyDict.from_json(**fg_args)
reaction_classes = fuzzy_dict.FuzzyDict.from_json(**reaction_class_args)
ring_smiles = fuzzy_dict.FuzzyDict.from_json(**ring_smiles_args)

checker = check.Check(
    fg_dict=functional_groups, reaction_dict=reaction_classes, ring_dict=ring_smiles
)


def main(route) -> Tuple[bool, Dict]:
    """
    Detects the presence of 2,4-dichlorophenyl thioether motif in early stages of synthesis
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

    found_motif = False
    early_stage_threshold = 3  # Define what "early stage" means in terms of depth

    def dfs_traverse(node, depth=0):
        nonlocal found_motif, findings_json

        if node["type"] == "mol":
            mol_smiles = node["smiles"]

            # Create RDKit molecule for more detailed analysis
            mol = Chem.MolFromSmiles(mol_smiles)
            if mol:
                # Look specifically for 2,4-dichlorophenyl connected to sulfur
                # Pattern for 2,4-dichlorophenyl thioether: chlorines at positions 2 and 4, sulfur connected to position 1
                pattern = Chem.MolFromSmarts("c1c(Cl)cc(Cl)cc1S")
                if mol.HasSubstructMatch(pattern):
                    findings_json["atomic_checks"]["functional_groups"].append("2,4-dichlorophenyl thioether")
                    # Check if we're in early stages of synthesis
                    if depth >= early_stage_threshold:
                        found_motif = True
                        findings_json["structural_constraints"].append({
                            "type": "positional",
                            "details": {
                                "target": "2,4-dichlorophenyl thioether",
                                "position": "early_stage"
                            }
                        })

        # Determine the new depth for the recursive call
        new_depth = depth
        if node["type"] != "reaction":  # If current node is chemical, depth increases
            new_depth = depth + 1

        # Traverse children (going deeper into the synthesis route)
        for child in node.get("children", []):
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    return found_motif, findings_json

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
    Checks if the final product of the synthesis is a molecule containing both a piperidine ring and at least one stereocenter.
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

    # Track stereocenters at each step and store SMILES
    stereocenters_by_depth = {}
    stereocenters_by_depth_smiles = {}

    def dfs_traverse(node, depth=0):
        if node["type"] == "mol" and "smiles" in node:
            smiles = node["smiles"]
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                # Properly count stereocenters using RDKit
                chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=False)
                stereo_count = len(chiral_centers)
                stereocenters_by_depth[depth] = stereo_count
                stereocenters_by_depth_smiles[depth] = smiles

        # Determine the next depth based on the current node's type
        next_depth = depth
        if node["type"] != "reaction":  # If current node is 'chemical' (or 'mol'), depth increases
            next_depth = depth + 1

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, next_depth)

    # Start traversal
    dfs_traverse(route)

    result = False

    if not stereocenters_by_depth:
        return result, findings_json

    # Get the depths in order
    depths = sorted(stereocenters_by_depth.keys())
    if not depths:
        return result, findings_json

    # In retrosynthesis, depth 0 is the final product (target molecule)
    final_product_stereocenters = stereocenters_by_depth.get(min(depths), 0)

    # Check if the final product has stereocenters and contains piperidine
    final_product_smiles = stereocenters_by_depth_smiles.get(min(depths), "")
    has_piperidine_scaffold = checker.check_ring("piperidine", final_product_smiles)

    if final_product_stereocenters > 0:
        findings_json["atomic_checks"]["functional_groups"].append("stereocenter")
        findings_json["structural_constraints"].append({"type": "count", "details": {"target": "stereocenter", "operator": ">", "value": 0}})

    if has_piperidine_scaffold:
        findings_json["atomic_checks"]["ring_systems"].append("piperidine")
        findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "piperidine", "position": "last_stage"}})

    if final_product_stereocenters > 0 and has_piperidine_scaffold:
        result = True
        findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"targets": ["piperidine", "stereocenter"], "location": "last_stage"}})

    return result, findings_json

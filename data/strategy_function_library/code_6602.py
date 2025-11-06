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
    Detects if the synthesis involves formation of C-S bonds, particularly in late stages.
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

    cs_bond_formation = False
    cs_bond_depth = None

    # List of reactions that typically form C-S bonds
    cs_bond_reactions = [
        "S-alkylation of thiols",
        "S-alkylation of thiols (ethyl)",
        "S-alkylation of thiols with alcohols",
        "S-alkylation of thiols with alcohols (ethyl)",
        "thia-Michael addition",
        "Newman-Kwart rearrangement",
        "thioether_nucl_sub",
    ]

    # List of functional groups containing C-S bonds
    cs_functional_groups = [
        "Thioamide",
        "Thiourea",
        "Monosulfide",
        "Disulfide",
        "Thiocarbonyl",
        "Carbo-thioester",
        "Thiocyanate",
        "Isothiocyanate",
    ]

    def dfs_traverse(node, depth=0):
        nonlocal cs_bond_formation, cs_bond_depth, findings_json

        if node["type"] == "reaction":
            if "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                # The original implementation using reaction name lists and functional group
                # checks is prone to both false positives (e.g., FG interconversion)
                # and false negatives (e.g., unlisted reactions).
                # A direct check for bond formation using atom-mapping is the most
                # robust and chemically accurate method.
                if checker.check_bond_formation('C-S', rsmi):
                    cs_bond_formation = True
                    findings_json["atomic_checks"]["named_reactions"].append("C-S_bond_formation")
                    if cs_bond_depth is None or depth < cs_bond_depth:
                        cs_bond_depth = depth

        # Determine the depth for the recursive call based on the current node's type
        new_depth = depth
        if node["type"] != "reaction": # This means it's a 'chemical' node or similar
            new_depth = depth + 1

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, new_depth)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if C-S bond formation occurs in late stages (depth <= 1)
    late_stage_cs_formation = cs_bond_formation and cs_bond_depth is not None and cs_bond_depth <= 1
    
    if late_stage_cs_formation:
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "C-S_bond_formation",
                "position": "late_stage",
                "definition": "depth <= 1"
            }
        })

    print(
        f"C-S bond formation: {cs_bond_formation}, depth: {cs_bond_depth}, late stage: {late_stage_cs_formation}"
    )
    return late_stage_cs_formation, findings_json

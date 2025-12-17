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
from synth_strategy.utils.check import Check
from synth_strategy.utils import fuzzy_dict, check

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


BOC_PROTECTION_REACTIONS = [
    "Boc amine protection",
    "Boc amine protection explicit",
    "Boc amine protection with Boc anhydride",
    "Boc amine protection (ethyl Boc)",
    "Boc amine protection of secondary amine",
    "Boc amine protection of primary amine",
]
BOC_DEPROTECTION_REACTIONS = [
    "Boc amine deprotection",
    "Boc amine deprotection of guanidine",
    "Boc amine deprotection to NH-NH2",
    "Tert-butyl deprotection of amine",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a 'long-lived' Boc protecting group strategy. The function identifies routes
    where a Boc group is introduced early (depth >= 3) and maintained until late stages
    (depth <= 2). It confirms this by checking for the presence of Boc-protected
    intermediates and identifying reactions from a defined list of Boc protection and
    deprotection reaction types.
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

    boc_info = {
        "protected_molecules": [],
        "protection_reactions": [],
        "deprotection_reactions": [],
        "depths": {},
    }

    def dfs_traverse(node, depth=0):
        if node["type"] == "mol":
            mol_smiles = node["smiles"]
            if checker.check_fg("Boc", mol_smiles):
                boc_info["protected_molecules"].append(mol_smiles)
                if depth not in boc_info["depths"]:
                    boc_info["depths"][depth] = []
                boc_info["depths"][depth].append(mol_smiles)
                if "Boc" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Boc")

        elif node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rxn_smiles = node["metadata"]["mapped_reaction_smiles"]

            for r in BOC_PROTECTION_REACTIONS:
                if checker.check_reaction(r, rxn_smiles):
                    boc_info["protection_reactions"].append((depth, rxn_smiles))
                    if r not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append(r)
                    break # Assuming only one match per reaction node

            for r in BOC_DEPROTECTION_REACTIONS:
                if checker.check_reaction(r, rxn_smiles):
                    boc_info["deprotection_reactions"].append((depth, rxn_smiles))
                    if r not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append(r)
                    break # Assuming only one match per reaction node

        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth if node["type"] == "reaction" else depth + 1
            dfs_traverse(child, new_depth)

    dfs_traverse(route)

    depths = list(boc_info["depths"].keys())

    result = False

    if not depths:
        return result, findings_json

    early_protection = max(depths) >= 3
    late_stage_presence = min(depths) <= 2
    multiple_depths = len(depths) >= 2

    if early_protection:
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "Boc",
                "position_type": "max_depth",
                "operator": ">=",
                "value": 3
            }
        })
    if late_stage_presence:
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "Boc",
                "position_type": "min_depth",
                "operator": "<=",
                "value": 2
            }
        })
    if multiple_depths:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "distinct_depths_with_Boc",
                "operator": ">=",
                "value": 2
            }
        })

    if early_protection and late_stage_presence and multiple_depths:
        result = True

    return result, findings_json
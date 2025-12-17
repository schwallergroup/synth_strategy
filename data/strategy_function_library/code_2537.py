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

def main(route) -> Tuple[bool, Dict]:
    """
    Detects the protection of a ketone as a ketal, or the deprotection of a ketal back to a ketone, within a synthetic route.
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

    ketone_protected = False

    def dfs_traverse(node, depth=0):
        nonlocal ketone_protected, findings_json

        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants_smiles = rsmi.split(">")[0].split(".")

            # Check for ketone protection (ketone -> ketal)
            if checker.check_reaction(
                "Aldehyde or ketone acetalization", rsmi
            ):
                if any(checker.check_fg("Ketone", r) for r in reactants_smiles):
                    ketone_protected = True
                    findings_json["atomic_checks"]["named_reactions"].append("Aldehyde or ketone acetalization")
                    findings_json["atomic_checks"]["functional_groups"].append("Ketone")
            elif checker.check_reaction("Diol acetalization", rsmi):
                if any(checker.check_fg("Ketone", r) for r in reactants_smiles):
                    ketone_protected = True
                    findings_json["atomic_checks"]["named_reactions"].append("Diol acetalization")
                    findings_json["atomic_checks"]["functional_groups"].append("Ketone")

            # Check for ketal deprotection (ketal -> ketone)
            elif checker.check_reaction(
                "Ketal hydrolysis to ketone", rsmi
            ):
                ketone_protected = True
                findings_json["atomic_checks"]["named_reactions"].append("Ketal hydrolysis to ketone")
            elif checker.check_reaction("Acetal hydrolysis to ketone", rsmi):
                ketone_protected = True
                findings_json["atomic_checks"]["named_reactions"].append("Acetal hydrolysis to ketone")

        # Traverse children
        for child in node.get("children", []):
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth += 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    if ketone_protected:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "ketone_protection_or_deprotection_event",
                "operator": ">=",
                "value": 1
            }
        })

    return ketone_protected, findings_json

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


# Refactoring for Enumeration: Isolate the list of reaction types.
LATE_STAGE_AMIDE_COUPLING_REACTIONS = [
    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
    "Acyl chloride with ammonia to amide",
    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
    "Acyl chloride with secondary amine to amide",
    "Carboxylic acid with primary amine to amide",
    "Ester with ammonia to amide",
    "Ester with primary amine to amide",
    "Ester with secondary amine to amide",
    "Schotten-Baumann_amide",
    "Acylation of primary amines",
    "Acylation of secondary amines",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if the synthesis route involves a late-stage amide coupling
    (formation of an amide bond in the final or penultimate step).
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

    found_late_amide = False
    max_depth = 1  # Consider only reactions at depth 1 as "late-stage"

    def dfs_traverse(node, current_depth=0):
        nonlocal found_late_amide, findings_json

        if node["type"] == "reaction":
            if "mapped_reaction_smiles" in node.get("metadata", {}):
                rsmi = node["metadata"]["mapped_reaction_smiles"]

                # Only check reactions at depth 1 (the final synthetic step)
                if current_depth <= max_depth:
                    # Check for specific amide coupling reactions
                    for reaction_type in LATE_STAGE_AMIDE_COUPLING_REACTIONS:
                        if checker.check_reaction(reaction_type, rsmi):
                            found_late_amide = True
                            findings_json["atomic_checks"]["named_reactions"].append(reaction_type)
                            # Add structural constraint if this is the first time it's met
                            if {"type": "positional", "details": {"target": "any_amide_coupling_reaction", "position": "depth <= 1"}} not in findings_json["structural_constraints"]:
                                findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "any_amide_coupling_reaction", "position": "depth <= 1"}})
                            break

        # Traverse children with incremented depth based on node type
        for child in node.get("children", []):
            if not found_late_amide:  # Stop traversal if we already found what we're looking for
                new_depth = current_depth
                if node["type"] != "reaction": # If current node is chemical, depth increases
                    new_depth += 1
                # If current node is reaction, depth remains the same
                dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    return found_late_amide, findings_json

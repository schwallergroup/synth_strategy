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


AMIDE_COUPLING_REACTIONS = [
    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_OS",
    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
    "Acyl chloride with secondary amine to amide",
    "Carboxylic acid with primary amine to amide",
    "Ester with primary amine to amide",
    "Ester with secondary amine to amide",
    "Schotten-Baumann_amide",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if the synthesis route uses a late-stage amide coupling strategy. This is identified by checking if a reaction in the final step (depth <= 1) matches any of the reaction types defined in the `AMIDE_COUPLING_REACTIONS` list.
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

    amide_coupling_found = False

    def dfs_traverse(node, depth=0):
        nonlocal amide_coupling_found, findings_json

        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            # Check if this is at depth 0 or 1 (late stage)
            if depth <= 1:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                # Check for various amide coupling reactions
                is_amide_coupling = False
                for reaction_name in AMIDE_COUPLING_REACTIONS:
                    if checker.check_reaction(reaction_name, rsmi):
                        is_amide_coupling = True
                        findings_json["atomic_checks"]["named_reactions"].append(reaction_name)
                        break

                if is_amide_coupling:
                    amide_coupling_found = True
                    # Add the structural constraint if found at late stage
                    findings_json["structural_constraints"].append({
                        "type": "positional",
                        "details": {
                            "target": AMIDE_COUPLING_REACTIONS,
                            "position": "late_stage (depth <= 1)"
                        }
                    })
                    return  # Early return once found

        # Only continue traversal if we haven't found an amide coupling yet
        if not amide_coupling_found:
            # Traverse children with incremented depth
            for child in node.get("children", []):
                # New logic for depth calculation
                new_depth = depth
                if node["type"] != "reaction": # If current node is chemical, depth increases
                    new_depth = depth + 1
                # If current node is reaction, depth remains the same
                
                dfs_traverse(child, new_depth)
                # Early return if amide coupling found in children
                if amide_coupling_found:
                    return

    # Start traversal from the root
    dfs_traverse(route)
    return amide_coupling_found, findings_json

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


AMIDE_COUPLING_REACTIONS = [
    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_OS",
    "Acyl chloride with ammonia to amide",
    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
    "Acyl chloride with primary amine to imide",
    "Acyl chloride with secondary amine to amide",
    "Carboxylic acid with primary amine to amide",
    "Ester with ammonia to amide",
    "Ester with primary amine to amide",
    "Ester with secondary amine to amide",
    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
    "Schotten-Baumann_amide",
    "Acylation of primary amines",
    "Acylation of secondary amines",
    "Carboxylic acid to amide conversion",
    "Nitrile and hydrogen peroxide to amide",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if the final synthetic step involves an amide formation reaction from a predefined list.
    This function checks if the reaction at depth=1 matches any of the specified named reactions
    or reaction types, such as 'Acylation of Nitrogen Nucleophiles by Carboxylic Acids' or 'Schotten-Baumann_amide'.
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

    amide_formation_at_final_step = False

    def dfs_traverse(reaction, depth, max_depth):
        nonlocal amide_formation_at_final_step, findings_json

        if amide_formation_at_final_step:
            return

        if reaction.get("type") == "reaction" and depth == 1:
            rsmi = reaction.get("metadata", {}).get("rsmi")
            if not rsmi:
                return

            for reaction_type in AMIDE_COUPLING_REACTIONS:
                if checker.check_reaction(reaction_type, rsmi):
                    amide_formation_at_final_step = True
                    findings_json["atomic_checks"]["named_reactions"].append(reaction_type)
                    # Add the structural constraint if the condition is met
                    findings_json["structural_constraints"].append({
                        "type": "positional",
                        "details": {
                            "target": AMIDE_COUPLING_REACTIONS,
                            "position": "last_stage"
                        }
                    })
                    return

        for child in reaction.get("children", []):
            # Determine the new depth based on the current node's type
            new_depth = depth
            if reaction.get("type") != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            
            dfs_traverse(child, new_depth, max_depth)

    # The analysis framework is expected to provide max_depth.
    # For this function, its value is not used, so a placeholder is acceptable.
    max_d = 0
    dfs_traverse(route, 0, max_d)
    return amide_formation_at_final_step, findings_json

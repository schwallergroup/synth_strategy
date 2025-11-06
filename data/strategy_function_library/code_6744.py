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
    "Carboxylic acid with primary amine to amide",
    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
    "Acyl chloride with secondary amine to amide",
    "Ester with primary amine to amide",
    "Ester with secondary amine to amide",
    "Schotten-Baumann_amide",
    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
    "Acyl chloride with ammonia to amide",
    "Acyl chloride with primary amine to imide",
    "Ester with ammonia to amide",
]

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects if the synthesis uses a late-stage amide coupling strategy.
    It checks if the final reaction (depth 0) or penultimate reaction (depth 1) involves formation of an amide bond.
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

    late_stage_amide_coupling = False

    def dfs_traverse(node, depth=0):
        nonlocal late_stage_amide_coupling, findings_json

        if node["type"] == "reaction" and depth <= 1:  # Check depths 0 and 1 for late-stage
            # Check if this is an amide coupling reaction
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]

                is_amide_coupling = False
                for rxn in AMIDE_COUPLING_REACTIONS:
                    if checker.check_reaction(rxn, rsmi):
                        is_amide_coupling = True
                        findings_json["atomic_checks"]["named_reactions"].append(rxn)
                        break

                if is_amide_coupling:
                    late_stage_amide_coupling = True
                    # Add structural constraint if an amide coupling reaction is found at depth <= 1
                    findings_json["structural_constraints"].append({
                        "type": "positional",
                        "details": {
                            "target": "amide_coupling_reaction",
                            "position": "depth <= 1"
                        }
                    })
            except Exception as e:
                pass

        # Traverse children
        for child in node.get("children", []):
            # New logic: depth increases only when going from chemical to reaction
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else: # node["type"] == "chemical"
                dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)
    return late_stage_amide_coupling, findings_json

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
    "Acylation of primary amines",
    "Acylation of secondary amines",
    "Schotten-Baumann_amide",
    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_OS",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Checks for the use of a late-stage amide coupling strategy. This is confirmed if the final reaction step matches one of the specified amide formation reactions in the AMIDE_COUPLING_REACTIONS list.
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

    final_amide_coupling = False

    def dfs_traverse(node, current_depth=0):
        nonlocal final_amide_coupling, findings_json

        if node["type"] == "reaction" and node.get("metadata", {}).get("rsmi"):
            depth = node.get("metadata", {}).get("depth")
            if depth is None:
                depth = current_depth
            else:
                try:
                    depth = int(depth)
                except:
                    depth = current_depth

            rsmi = node["metadata"]["mapped_reaction_smiles"]

            if depth <= 1:
                is_amide_coupling = False
                for rxn_type in AMIDE_COUPLING_REACTIONS:
                    if checker.check_reaction(rxn_type, rsmi):
                        is_amide_coupling = True
                        findings_json["atomic_checks"]["named_reactions"].append(rxn_type)
                        break

                if is_amide_coupling:
                    final_amide_coupling = True

        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                # If current node is a reaction, depth remains the same for children (chemicals)
                dfs_traverse(child, current_depth)
            else:
                # If current node is a chemical, depth increases for children (reactions)
                dfs_traverse(child, current_depth + 1)

    dfs_traverse(route)

    if final_amide_coupling:
        # Add the structural constraint if the condition is met
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "amide_coupling_reaction",
                "position": "last_stage"
            }
        })

    return final_amide_coupling, findings_json

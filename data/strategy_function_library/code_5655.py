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


CYCLOADDITION_REACTIONS_OF_INTEREST = [
    "Diels-Alder",
    "Diels-Alder (ON bond)",
    "Huisgen alkyne-azide 1,3 dipolar cycloaddition",
    "Huisgen 1,3 dipolar cycloaddition",
    "Huisgen alkene-azide 1,3 dipolar cycloaddition",
    "[3+2]-cycloaddition of hydrazone and alkyne",
    "[3+2]-cycloaddition of hydrazone and alkene",
    "[3+2]-cycloaddition of diazoalkane and alkyne",
    "[3+2]-cycloaddition of diazoalkane and alkene",
    "[3+2]-cycloaddition of diazoalkane and alpha-alkyne",
    "[3+2]-cycloaddition of diazoalkane and alpha-alkene",
    "Retro-Diels-Alder from oxazole",
    "Azide-nitrile click cycloaddition to tetrazole",
    "Azide-nitrile click cycloaddition to triazole",
    "{Huisgen_Cu-catalyzed_1,4-subst}",
    "{Huisgen_Ru-catalyzed_1,5_subst}",
    "{Huisgen_disubst-alkyne}",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if the synthesis uses a specific, named cycloaddition reaction from an enumerated list.
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

    has_cycloaddition = False

    def dfs_traverse(node, depth=0):
        nonlocal has_cycloaddition, findings_json

        if has_cycloaddition:
            return

        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]

            for reaction_type in CYCLOADDITION_REACTIONS_OF_INTEREST:
                if checker.check_reaction(reaction_type, rsmi):
                    has_cycloaddition = True
                    findings_json["atomic_checks"]["named_reactions"].append(reaction_type)
                    break

        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else: # Assuming 'chemical' or other types
                dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    return has_cycloaddition, findings_json

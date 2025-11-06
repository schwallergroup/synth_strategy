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


AZIDE_CLICK_CHEMISTRY_REACTIONS = [
    "Huisgen alkyne-azide 1,3 dipolar cycloaddition",
    "Huisgen 1,3 dipolar cycloaddition",
    "Huisgen alkene-azide 1,3 dipolar cycloaddition",
    "Huisgen_Cu-catalyzed_1,4-subst",
    "Huisgen_Ru-catalyzed_1,5_subst",
    "Huisgen_disubst-alkyne",
]

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects a synthetic strategy involving azide chemistry,
    specifically the conversion of an amine to an azide group, azide to amine reduction,
    or azide-based click chemistry reactions.
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

    azide_strategy_found = False

    def dfs_traverse(node, depth=0):
        nonlocal azide_strategy_found, findings_json

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["mapped_reaction_smiles"]

                # 1. Amine to azide conversion
                if checker.check_reaction("Amine to azide", rsmi):
                    azide_strategy_found = True
                    findings_json["atomic_checks"]["named_reactions"].append("Amine to azide")

                # 2. Formation of azides from halogens
                elif checker.check_reaction("Formation of Azides from halogens", rsmi):
                    azide_strategy_found = True
                    findings_json["atomic_checks"]["named_reactions"].append("Formation of Azides from halogens")

                # 3. Formation of azides from boronic acids
                elif checker.check_reaction("Formation of Azides from boronic acids", rsmi):
                    azide_strategy_found = True
                    findings_json["atomic_checks"]["named_reactions"].append("Formation of Azides from boronic acids")

                # 4. Azide reduction (Staudinger)
                elif checker.check_reaction("Azide to amine reduction (Staudinger)", rsmi):
                    azide_strategy_found = True
                    findings_json["atomic_checks"]["named_reactions"].append("Azide to amine reduction (Staudinger)")

                # 5. Huisgen cycloaddition (click chemistry)
                elif any(
                    checker.check_reaction(rxn, rsmi)
                    for rxn in AZIDE_CLICK_CHEMISTRY_REACTIONS
                ):
                    for rxn in AZIDE_CLICK_CHEMISTRY_REACTIONS:
                        if checker.check_reaction(rxn, rsmi):
                            azide_strategy_found = True
                            findings_json["atomic_checks"]["named_reactions"].append(rxn)
                            break # Only add the first matching reaction

        # Traverse children
        for child in node.get("children", []):
            # New logic: depth increases only when traversing from chemical to reaction node
            # Depth remains the same when traversing from reaction to chemical node
            new_depth = depth + 1 if node["type"] != "reaction" else depth
            dfs_traverse(child, new_depth)

    dfs_traverse(route)

    if azide_strategy_found:
        # Add the structural constraint if any azide chemistry reaction was found
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "azide_chemistry_reaction",
                "operator": ">=",
                "value": 1
            }
        })

    return azide_strategy_found, findings_json

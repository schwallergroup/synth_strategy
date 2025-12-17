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


HETEROCYCLE_FORMATION_REACTIONS = [
    "benzimidazole_derivatives_carboxylic-acid/ester",
    "benzimidazole_derivatives_aldehyde",
    "benzothiazole",
    "benzoxazole_arom-aldehyde",
    "benzoxazole_carboxylic-acid",
    "thiazole",
    "tetrazole_terminal",
    "tetrazole_connect_regioisomere_1",
    "tetrazole_connect_regioisomere_2",
    "1,2,4-triazole_acetohydrazide",
    "1,2,4-triazole_carboxylic-acid/ester",
    "pyrazole",
    "Paal-Knorr pyrrole",
    "triaryl-imidazole",
    "Fischer indole",
    "Friedlaender chinoline",
    "benzofuran",
    "benzothiophene",
    "indole",
    "oxadiazole",
    "Pictet-Spengler",
    "Huisgen_Cu-catalyzed_1,4-subst",
    "Huisgen_Ru-catalyzed_1,5_subst",
    "Huisgen_disubst-alkyne",
    "Huisgen alkyne-azide 1,3 dipolar cycloaddition",
    "Huisgen 1,3 dipolar cycloaddition",
    "Huisgen alkene-azide 1,3 dipolar cycloaddition",
    "Azide-nitrile click cycloaddition to tetrazole",
    "Azide-nitrile click cycloaddition to triazole",
    "Formation of NOS Heterocycles",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if the final synthetic step involves a named reaction that forms a heterocycle.
    The function checks against a predefined list of known heterocycle-forming reactions,
    such as the Fischer indole synthesis, Paal-Knorr synthesis, and various cycloadditions.
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

    has_late_stage_ring_formation = False

    def dfs_traverse(node, depth=0):
        nonlocal has_late_stage_ring_formation, findings_json

        if node["type"] == "reaction" and depth <= 1:  # Final or penultimate reaction in synthesis
            if "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
                rsmi = node["metadata"]["mapped_reaction_smiles"]

                # Check for specific heterocycle formation reactions
                for reaction in HETEROCYCLE_FORMATION_REACTIONS:
                    if checker.check_reaction(reaction, rsmi):
                        has_late_stage_ring_formation = True
                        findings_json["atomic_checks"]["named_reactions"].append(reaction)
                        # Add structural constraint if detected at the correct depth
                        findings_json["structural_constraints"].append({
                            "type": "positional",
                            "details": {
                                "target": "heterocycle_formation_reaction",
                                "position": "final_or_penultimate_stage"
                            }
                        })
                        return

        # Traverse children
        for child in node.get("children", []):
            if has_late_stage_ring_formation:
                return
            
            # Determine the new depth based on the current node's type
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            
            dfs_traverse(child, new_depth)

    # Start traversal from root
    dfs_traverse(route)

    return has_late_stage_ring_formation, findings_json

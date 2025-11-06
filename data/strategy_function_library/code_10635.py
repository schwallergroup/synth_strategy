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


def main(route) -> Tuple[bool, Dict]:
    """
    Detects if the synthesis route involves construction of a benzothiophene core.
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

    benzothiophene_found = False

    def dfs_traverse(reaction, depth):
        nonlocal benzothiophene_found, findings_json

        if reaction["type"] == "reaction" and reaction.get("metadata", {}).get("rsmi"):
            rsmi = reaction["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check if product contains benzothiophene
            if checker.check_ring("benzothiophene", product_smiles):
                findings_json["atomic_checks"]["ring_systems"].append("benzothiophene")

                # Check if any reactant contains benzothiophene
                reactant_has_benzothiophene = False
                for reactant_smiles in reactants_smiles:
                    if checker.check_ring("benzothiophene", reactant_smiles):
                        reactant_has_benzothiophene = True
                        findings_json["atomic_checks"]["ring_systems"].append("benzothiophene")
                        break

                # If benzothiophene is in product but not in any reactant, it was constructed
                if not reactant_has_benzothiophene:
                    benzothiophene_found = True
                    # Add the structural constraint when benzothiophene is found to be constructed
                    findings_json["structural_constraints"].append({
                        "type": "count",
                        "details": {
                            "target": "formation of benzothiophene",
                            "operator": ">=",
                            "value": 1
                        }
                    })

        # Continue traversal
        for child in reaction.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if reaction["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    dfs_traverse(route, 1)
    return benzothiophene_found, findings_json

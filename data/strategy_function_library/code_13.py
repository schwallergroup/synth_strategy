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


AMIDE_FGS = ["Primary amide", "Secondary amide", "Tertiary amide"]

def main(route) -> Tuple[bool, Dict]:
    """
    Identifies syntheses with at least two reaction steps where a substrate containing both an aromatic halide and an amide group is transformed while preserving both groups. This highlights a repeated chemoselective strategy. The specific amide types checked are Primary amide, Secondary amide, and Tertiary amide.
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

    multiple_transformations = 0

    def dfs_traverse(node, depth=0):
        nonlocal multiple_transformations, findings_json

        if node["type"] == "reaction":
            if "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                product_has_halogen = checker.check_fg("Aromatic halide", product)
                if product_has_halogen:
                    if "Aromatic halide" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Aromatic halide")

                product_has_amide = False
                for fg in AMIDE_FGS:
                    if checker.check_fg(fg, product):
                        product_has_amide = True
                        if fg not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append(fg)
                
                if product_has_halogen and product_has_amide:
                    for reactant in reactants:
                        reactant_has_halogen = checker.check_fg("Aromatic halide", reactant)
                        if reactant_has_halogen:
                            if "Aromatic halide" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Aromatic halide")

                        reactant_has_amide = False
                        for fg in AMIDE_FGS:
                            if checker.check_fg(fg, reactant):
                                reactant_has_amide = True
                                if fg not in findings_json["atomic_checks"]["functional_groups"]:
                                    findings_json["atomic_checks"]["functional_groups"].append(fg)

                        if reactant_has_halogen and reactant_has_amide:
                            multiple_transformations += 1
                            break

        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                # If current node is a reaction, depth remains the same for children (chemical nodes)
                dfs_traverse(child, depth)
            else:
                # If current node is a chemical, depth increases for children (reaction nodes)
                dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    result = multiple_transformations >= 2
    if result:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "reaction_step_where_reactant_and_product_both_contain_an_aromatic_halide_and_an_amide",
                "operator": ">=",
                "value": 2
            }
        })

    return result, findings_json

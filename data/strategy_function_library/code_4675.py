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


LATE_STAGE_FGI_REACTIONS = [
    "Reduction of carboxylic acid to primary alcohol",
    "Reduction of ester to primary alcohol",
    "Reduction of nitrile to amine",
    "Reduction of aldehydes and ketones to alcohols",
    "Reduction of primary amides to amines",
    "Reduction of secondary amides to amines",
    "Reduction of tertiary amides to amines",
    "Reduction of nitro groups to amines",
    "Oxidation of aldehydes to carboxylic acids",
    "Oxidation of alcohol to carboxylic acid",
    "Oxidation of nitrile to carboxylic acid",
    "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters",
    "Nitrile to amide",
    "Esterification of Carboxylic Acids",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Identifies routes that employ a convergent strategy combined with late-stage functionalization. A route is flagged if it contains at least one convergent step (>= 2 reactants) and at least one late-stage (final step) reaction matching a defined list of functional group interconversions, such as reductions, oxidations, and hydrolysis.
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

    has_convergent_step = False
    has_late_stage_functionalization = False

    def dfs_traverse(node, depth=0):
        nonlocal has_convergent_step, has_late_stage_functionalization, findings_json

        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants = rsmi.split(">")[0].split(".")

            # Check for convergent synthesis (multiple reactants)
            if len(reactants) >= 2:
                has_convergent_step = True
                findings_json["structural_constraints"].append({
                    "type": "count",
                    "details": {
                        "target": "reactants_in_a_step",
                        "operator": ">=",
                        "value": 2
                    }
                })
                # print(f"Detected convergent step at depth {depth}")

            # Check for late-stage functionalization (depth <= 1)
            if depth <= 1:
                try:
                    # Check for reaction types
                    for reaction_type in LATE_STAGE_FGI_REACTIONS:
                        if checker.check_reaction(reaction_type, rsmi):
                            # print(
                            #     f"Detected late-stage functionalization reaction: {reaction_type} at depth {depth}"
                            # )
                            has_late_stage_functionalization = True
                            findings_json["atomic_checks"]["named_reactions"].append(reaction_type)
                            findings_json["structural_constraints"].append({
                                "type": "positional",
                                "details": {
                                    "target": "any_LATE_STAGE_FGI_REACTIONS",
                                    "position": "last_or_penultimate_step"
                                }
                            })
                            break

                except Exception as e:
                    print(
                        f"Error processing reaction SMILES for late-stage functionalization detection: {e}"
                    )

        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    dfs_traverse(route)

    result = has_convergent_step and has_late_stage_functionalization

    if result:
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "convergent_step",
                    "late_stage_functionalization"
                ]
            }
        })

    return result, findings_json

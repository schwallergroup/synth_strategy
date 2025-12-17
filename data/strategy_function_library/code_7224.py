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
rng_smiles_args = {
    "file_path": f"{root_data}/patterns/chemical_rings_smiles.json",
    "value_field": "smiles",
    "key_field": "name",
}
functional_groups = fuzzy_dict.FuzzyDict.from_json(**fg_args)
reaction_classes = fuzzy_dict.FuzzyDict.from_json(**reaction_class_args)
ring_smiles = fuzzy_dict.FuzzyDict.from_json(**rng_smiles_args)

checker = check.Check(
    fg_dict=functional_groups, reaction_dict=reaction_classes, ring_dict=ring_smiles
)


FG_INTERCONVERSION_REACTION_TYPES = [
    "Oxidation of aldehydes to carboxylic acids",
    "Alcohol protection with silyl ethers",
    "Alcohol deprotection from silyl ethers",
    "Oxidation or Dehydrogenation of Alcohols to Aldehydes and Ketones",
    "Williamson Ether Synthesis",
    "Reduction of ester to primary alcohol",
    "Reduction of ketone to secondary alcohol",
    "Reduction of carboxylic acid to primary alcohol",
    "Oxidation of alcohol to carboxylic acid",
    "Esterification of Carboxylic Acids",
    "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters",
    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
    "Reduction of nitro groups to amines",
    "Reduction of nitrile to amide",
    "Reduction of nitrile to amine",
    "Boc amine protection",
    "Boc amine deprotection",
]

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects if the synthetic strategy focuses on functional group
    interconversions, identified from a specific list of reaction types,
    rather than carbon skeleton construction.
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

    fg_interconversions = 0
    carbon_skeleton_changes = 0
    total_reactions = 0

    def dfs_traverse(node, depth=0):
        nonlocal fg_interconversions, carbon_skeleton_changes, total_reactions, findings_json

        node["depth"] = depth # Assign depth to the current node

        if node["type"] == "reaction":
            total_reactions += 1

            try:
                # Extract reactants and product from reaction SMILES
                rsmi = node["metadata"]["mapped_reaction_smiles"]

                # Check if this is a known functional group interconversion reaction
                fg_interconversion_detected = False
                for rxn_type in FG_INTERCONVERSION_REACTION_TYPES:
                    if checker.check_reaction(rxn_type, rsmi):
                        fg_interconversions += 1
                        findings_json["atomic_checks"]["named_reactions"].append(rxn_type)
                        print(
                            f"Detected functional group interconversion reaction at depth {node.get('depth', 'unknown')}"
                        )
                        fg_interconversion_detected = True
                        break

                if not fg_interconversion_detected:
                    carbon_skeleton_changes += 1
                    print(
                        f"Detected carbon skeleton change at depth {node.get('depth', 'unknown')}"
                    )
            except Exception as e:
                print(f"Error processing reaction: {e}")
                # If we can't analyze the reaction, we'll count it as a carbon skeleton change
                carbon_skeleton_changes += 1

        # Determine the next depth based on the current node's type
        next_depth = depth
        if node["type"] != "reaction": # If current node is chemical, depth increases
            next_depth = depth + 1

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, next_depth)

    # Start traversal
    dfs_traverse(route)

    # Strategy is present if most reactions are functional group interconversions
    fg_focus = (fg_interconversions > carbon_skeleton_changes) and (total_reactions > 0)
    print(
        f"FG interconversions: {fg_interconversions}, Carbon skeleton changes: {carbon_skeleton_changes}, Total reactions: {total_reactions}"
    )

    # Add structural constraints to findings_json if met
    if fg_interconversions > carbon_skeleton_changes:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "count(functional_group_interconversion) - count(carbon_skeleton_change)",
                "operator": ">",
                "value": 0
            }
        })
    if total_reactions > 0:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "total_reactions",
                "operator": ">",
                "value": 0
            }
        })

    return fg_focus, findings_json

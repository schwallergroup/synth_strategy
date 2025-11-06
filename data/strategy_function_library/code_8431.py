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


# Refactoring for Enumeration: Isolate the lists of reaction types.
OXIDATION_REACTIONS = [
    "Oxidation of aldehydes to carboxylic acids",
    "Oxidation of ketone to carboxylic acid",
    "Oxidation of alcohol to carboxylic acid",
    "Oxidation of primary alcohols",
    "Oxidation or Dehydrogenation of Alcohols to Aldehydes and Ketones",
    "Oxidation of alkene to aldehyde",
    "Oxidation of alkene to carboxylic acid",
    "Oxidation of alcohol and aldehyde to ester",
    "Oxidative esterification of primary alcohols",
    "Quinone formation",
    "Aromatic hydroxylation",
    "Sulfanyl to sulfinyl_peroxide",
    "Sulfanyl to sulfinyl_H2O2",
]

REDUCTION_REACTIONS = [
    "Reduction of aldehydes and ketones to alcohols",
    "Reduction of ester to primary alcohol",
    "Reduction of ketone to secondary alcohol",
    "Reduction of carboxylic acid to primary alcohol",
    "Reduction of nitro groups to amines",
    "Reduction of nitrile to amine",
    "Reduction of primary amides to amines",
    "Reduction of secondary amides to amines",
    "Reduction of tertiary amides to amines",
    "Nef reaction (nitro to ketone)",
    "Azide to amine reduction (Staudinger)",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if a synthetic route contains both an oxidation and a reduction step.
    The identification is performed by checking each reaction against a predefined list of named oxidation reactions and a predefined list of named reduction reactions.
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

    oxidation_found = False
    reduction_found = False

    def dfs_traverse(node, depth=0):
        nonlocal oxidation_found, reduction_found, findings_json

        if node["type"] == "reaction":
            if "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
                rsmi = node["metadata"]["mapped_reaction_smiles"]

                # Check for known oxidation reactions
                for rxn_type in OXIDATION_REACTIONS:
                    if checker.check_reaction(rxn_type, rsmi):
                        oxidation_found = True
                        if rxn_type not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(rxn_type)
                        break

                # Check for known reduction reactions
                for rxn_type in REDUCTION_REACTIONS:
                    if checker.check_reaction(rxn_type, rsmi):
                        reduction_found = True
                        if rxn_type not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(rxn_type)
                        break

        # Traverse children
        for child in node.get("children", []):
            # Stop traversing if both have been found
            if oxidation_found and reduction_found:
                break
            
            # New depth calculation logic
            if node["type"] == "reaction":
                # Depth remains the same when traversing from a reaction node to a chemical node
                dfs_traverse(child, depth)
            else:
                # Depth increases when traversing from a chemical node to a reaction node
                dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    result = oxidation_found and reduction_found
    
    if result:
        # Add the structural constraint if both oxidation and reduction are found
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "oxidation_reaction",
                    "reduction_reaction"
                ],
                "description": "The route must contain at least one reaction from a predefined list of oxidations and at least one reaction from a predefined list of reductions. The specific reactions checked are listed in atomic_checks."
            }
        })

    return result, findings_json

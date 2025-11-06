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

def main(route) -> Tuple[bool, Dict]:
    """
    Detects the masked amine strategy using azide as an intermediate:
    alcohol → sulfonate → azide → amine → amide

    In retrosynthetic analysis, we'll see this in reverse:
    amide → amine → azide → sulfonate → alcohol
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

    # Track if we've found each step in the sequence
    found_alcohol_activation = False
    found_azide_formation = False
    found_azide_reduction = False
    found_amide_coupling = False

    def dfs_traverse(node, depth=0):
        nonlocal found_alcohol_activation, found_azide_formation, found_azide_reduction, found_amide_coupling, findings_json

        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]

                # Check for amide formation reactions (retro: amide cleavage)
                amide_reactions = [
                    "Carboxylic acid with primary amine to amide",
                    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
                    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                    "Ester with primary amine to amide",
                    "Acylation of primary amines",
                    "Acylation of secondary amines"
                ]
                for r_name in amide_reactions:
                    if checker.check_reaction(r_name, rsmi):
                        found_amide_coupling = True
                        if r_name not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(r_name)

                # Check for azide reduction reactions (retro: amine to azide)
                azide_reduction_reactions = [
                    "Azide to amine reduction (Staudinger)"
                ]
                for r_name in azide_reduction_reactions:
                    if checker.check_reaction(r_name, rsmi):
                        found_azide_reduction = True
                        if r_name not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(r_name)

                # Check for azide formation reactions (retro: azide to activated alcohol/halide)
                azide_formation_reactions = [
                    "Formation of Azides from halogens",
                    "Formation of Azides from boronic acids",
                    "Alcohol to azide",
                    "Amine to azide"
                ]
                for r_name in azide_formation_reactions:
                    if checker.check_reaction(r_name, rsmi):
                        found_azide_formation = True
                        if r_name not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(r_name)

                # Check for alcohol activation reactions (retro: deactivation to alcohol)
                alcohol_activation_reactions = [
                    "Formation of Sulfonic Esters",
                    "Alcohol to triflate conversion",
                    "Formation of Sulfonic Esters on TMS protected alcohol",
                    "Alcohol to chloride_sulfonyl chloride",
                    "Alcohol to chloride_SOCl2",
                    "Alcohol to chloride_CHCl3",
                    "Alcohol to chloride_CH2Cl2",
                    "Alcohol to chloride_PCl5_ortho",
                    "Alcohol to chloride_POCl3_ortho",
                    "Alcohol to chloride_POCl3_para",
                    "Alcohol to chloride_POCl3",
                    "Alcohol to chloride_HCl",
                    "Alcohol to chloride_Salt",
                    "Alcohol to chloride_Other",
                    "Appel reaction"
                ]
                for r_name in alcohol_activation_reactions:
                    if checker.check_reaction(r_name, rsmi):
                        found_alcohol_activation = True
                        if r_name not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(r_name)

            except Exception:
                pass

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else:  # node["type"] == "chemical"
                dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if we found the complete sequence
    result = (
        found_alcohol_activation
        and found_azide_formation
        and found_azide_reduction
        and found_amide_coupling
    )

    if result:
        # Add the structural constraint if all conditions are met
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "alcohol_activation",
                    "azide_formation",
                    "azide_reduction",
                    "amide_coupling"
                ],
                "description": "The route must contain at least one reaction from each of the four key stages: alcohol activation (e.g., sulfonate formation), azide formation, azide reduction (e.g., Staudinger), and amide coupling."
            }
        })

    return result, findings_json

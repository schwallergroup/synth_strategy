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
    This function detects if the synthesis includes a halogenation followed by dehalogenation sequence.
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

    found_halogenation = False
    found_dehalogenation = False
    halogenation_depth = float("inf")
    dehalogenation_depth = float("inf")

    # List of halogenation reactions
    halogenation_reactions = [
        "Alcohol to chloride_SOCl2",
        "Alcohol to chloride_PCl5_ortho",
        "Alcohol to chloride_POCl3_ortho",
        "Alcohol to chloride_POCl3_para",
        "Alcohol to chloride_POCl3",
        "Alcohol to chloride_HCl",
        "Alcohol to chloride_Salt",
        "Alcohol to chloride_Other",
        "Alcohol to chloride_sulfonyl chloride",
        "Alcohol to chloride_CHCl3",
        "Alcohol to chloride_CH2Cl2",
        "Primary amine to chloride",
        "Primary amine to fluoride",
        "Primary amine to bromide",
        "Primary amine to iodide",
        "Aromatic chlorination",
        "Chlorination",
        "Aromatic fluorination",
        "Fluorination",
        "Aromatic bromination",
        "Bromination",
        "Aromatic iodination",
        "Iodination",
        "Wohl-Ziegler bromination benzyl primary",
        "Wohl-Ziegler bromination benzyl secondary",
        "Wohl-Ziegler bromination benzyl tertiary",
        "Wohl-Ziegler bromination allyl primary",
        "Wohl-Ziegler bromination allyl secondary",
        "Wohl-Ziegler bromination allyl tertiary",
        "Wohl-Ziegler bromination carbonyl primary",
        "Wohl-Ziegler bromination carbonyl secondary",
        "Wohl-Ziegler bromination carbonyl tertiary",
        "Halodeboronation of boronic acids",
        "Halodeboronation of boronic esters",
    ]

    # List of dehalogenation reactions
    dehalogenation_reactions = ["Dehalogenation", "Aromatic dehalogenation"]

    # List of halogen-containing functional groups
    halogen_fgs = [
        "Primary halide",
        "Secondary halide",
        "Tertiary halide",
        "Aromatic halide",
        "Alkenyl halide",
        "Haloalkyne",
        "Acyl halide",
        "Magnesium halide",
        "Zinc halide",
        "Sulfonyl halide",
    ]

    def dfs_traverse(node, depth=0):
        nonlocal found_halogenation, found_dehalogenation, halogenation_depth, dehalogenation_depth, findings_json

        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]

            # Check for halogenation reactions
            for reaction_name in halogenation_reactions:
                if checker.check_reaction(reaction_name, rsmi):
                    found_halogenation = True
                    halogenation_depth = min(halogenation_depth, depth)
                    if reaction_name not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append(reaction_name)
                    print(f"Found halogenation ({reaction_name}) at depth {depth}")
                    break

            # Check for dehalogenation reactions
            for reaction_name in dehalogenation_reactions:
                if checker.check_reaction(reaction_name, rsmi):
                    found_dehalogenation = True
                    dehalogenation_depth = min(dehalogenation_depth, depth)
                    if reaction_name not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append(reaction_name)
                    print(f"Found dehalogenation ({reaction_name}) at depth {depth}")
                    break

            # If we haven't identified specific reaction types, check for patterns
            if not found_halogenation or not found_dehalogenation:
                try:
                    reactants = rsmi.split(">")[0].split(".")
                    product = rsmi.split(">")[-1]

                    # Check for halogenation patterns
                    if not found_halogenation:
                        # Check if product contains halide
                        halide_in_product = False
                        for halogen_fg in halogen_fgs:
                            if checker.check_fg(halogen_fg, product):
                                halide_in_product = True
                                if halogen_fg not in findings_json["atomic_checks"]["functional_groups"]:
                                    findings_json["atomic_checks"]["functional_groups"].append(halogen_fg)
                        
                        # Check if reactants don't contain halide but product does
                        no_halide_in_reactants = True
                        for r in reactants:
                            for halogen_fg in halogen_fgs:
                                if checker.check_fg(halogen_fg, r):
                                    no_halide_in_reactants = False
                                    if halogen_fg not in findings_json["atomic_checks"]["functional_groups"]:
                                        findings_json["atomic_checks"]["functional_groups"].append(halogen_fg)
                                    break
                            if not no_halide_in_reactants: # If halide found in any reactant, break outer loop
                                break

                        if no_halide_in_reactants and halide_in_product:
                            found_halogenation = True
                            halogenation_depth = min(halogenation_depth, depth)
                            print(f"Found halogenation pattern at depth {depth}")

                    # Check for dehalogenation patterns
                    if not found_dehalogenation:
                        # Check if reactants contain halide
                        halide_in_reactants = False
                        for r in reactants:
                            for halogen_fg in halogen_fgs:
                                if checker.check_fg(halogen_fg, r):
                                    halide_in_reactants = True
                                    if halogen_fg not in findings_json["atomic_checks"]["functional_groups"]:
                                        findings_json["atomic_checks"]["functional_groups"].append(halogen_fg)
                                    break
                            if halide_in_reactants: # If halide found in any reactant, break outer loop
                                break

                        # Check if product doesn't contain halide
                        no_halide_in_product = True
                        for halogen_fg in halogen_fgs:
                            if checker.check_fg(halogen_fg, product):
                                no_halide_in_product = False
                                if halogen_fg not in findings_json["atomic_checks"]["functional_groups"]:
                                    findings_json["atomic_checks"]["functional_groups"].append(halogen_fg)
                                break

                        if halide_in_reactants and no_halide_in_product:
                            found_dehalogenation = True
                            dehalogenation_depth = min(dehalogenation_depth, depth)
                            print(f"Found dehalogenation pattern at depth {depth}")
                except Exception as e:
                    print(f"Error processing reaction SMILES: {e}")

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                # If current node is a reaction, depth remains the same for children (chemical nodes)
                dfs_traverse(child, depth)
            else:
                # If current node is a chemical, depth increases for children (reaction nodes)
                dfs_traverse(child, depth + 1)

    # Call dfs_traverse on the root node
    dfs_traverse(route)

    # In forward synthesis, halogenation happens first, then dehalogenation
    # In retrosynthetic analysis, halogenation should appear at a higher depth than dehalogenation
    correct_sequence = (
        found_halogenation and found_dehalogenation and halogenation_depth > dehalogenation_depth
    )

    if found_halogenation:
        print(f"Found halogenation at depth {halogenation_depth}")
    if found_dehalogenation:
        print(f"Found dehalogenation at depth {dehalogenation_depth}")

    if correct_sequence:
        print(
            f"Found halogenation-dehalogenation sequence: halogenation at depth {halogenation_depth}, dehalogenation at depth {dehalogenation_depth}"
        )
        findings_json["structural_constraints"].append({
            "type": "sequence",
            "details": {
                "first_event": "halogenation",
                "second_event": "dehalogenation",
                "description": "A halogenation event must occur earlier in the synthesis (higher depth in retrosynthesis tree) than a dehalogenation event."
            }
        })
    else:
        if found_halogenation and found_dehalogenation:
            print(
                f"Found both reactions but in wrong order: halogenation at depth {halogenation_depth}, dehalogenation at depth {dehalogenation_depth}"
            )
        elif not found_halogenation:
            print("No halogenation reaction found")
        elif not found_dehalogenation:
            print("No dehalogenation reaction found")

    return correct_sequence, findings_json

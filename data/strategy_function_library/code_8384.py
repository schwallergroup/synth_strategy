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
    Detects a synthesis with at least four distinct, sequential functionalization steps chosen from a predefined list of reaction types (e.g., amide/ester formation, halogenation, S_NAr, etc.). The strategy is identified if these steps occur in a generally sequential order throughout the synthesis.
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

    # Track functionalization steps with their depths
    functionalization_steps = []
    result = False

    def dfs_traverse(node, depth=0):
        nonlocal functionalization_steps, findings_json
        # Process reaction nodes
        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for common functionalization reactions
            reaction_detected = False

            # Amide formation
            amide_reactions = [
                "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
                "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                "Carboxylic acid with primary amine to amide",
                "Ester with primary amine to amide",
                "{Schotten-Baumann_amide}"
            ]
            for r_name in amide_reactions:
                if checker.check_reaction(r_name, rsmi):
                    functionalization_steps.append(("amide_formation", depth))
                    findings_json["atomic_checks"]["named_reactions"].append(r_name)
                    reaction_detected = True
                    print(f"Detected amide formation at depth {depth}: {rsmi}")
                    break

            # Ester formation
            if not reaction_detected:
                ester_reactions = [
                    "Esterification of Carboxylic Acids",
                    "Schotten-Baumann to ester",
                    "Transesterification",
                    "Oxidative esterification of primary alcohols",
                    "Acetic anhydride and alcohol to ester"
                ]
                for r_name in ester_reactions:
                    if checker.check_reaction(r_name, rsmi):
                        functionalization_steps.append(("ester_formation", depth))
                        findings_json["atomic_checks"]["named_reactions"].append(r_name)
                        reaction_detected = True
                        print(f"Detected ester formation at depth {depth}: {rsmi}")
                        break

            # Phenol alkylation / Ether formation
            if not reaction_detected:
                ether_reactions = [
                    "Williamson Ether Synthesis",
                    "{Williamson ether}",
                    "Alcohol to ether"
                ]
                for r_name in ether_reactions:
                    if checker.check_reaction(r_name, rsmi):
                        # Check if one of the reactants contains a phenol
                        is_phenol_alkylation = False
                        for reactant in reactants:
                            if checker.check_fg("Phenol", reactant):
                                functionalization_steps.append(("phenol_alkylation", depth))
                                findings_json["atomic_checks"]["functional_groups"].append("Phenol")
                                findings_json["atomic_checks"]["named_reactions"].append(r_name)
                                reaction_detected = True
                                is_phenol_alkylation = True
                                print(f"Detected phenol alkylation at depth {depth}: {rsmi}")
                                break

                        # If not a phenol alkylation but still an ether formation
                        if not is_phenol_alkylation:
                            functionalization_steps.append(("ether_formation", depth))
                            findings_json["atomic_checks"]["named_reactions"].append(r_name)
                            reaction_detected = True
                            print(f"Detected ether formation at depth {depth}: {rsmi}")
                        break

            # Nucleophilic aromatic substitution
            if not reaction_detected:
                nuc_sub_reactions = [
                    "nucl_sub_aromatic_ortho_nitro",
                    "nucl_sub_aromatic_para_nitro",
                    "heteroaromatic_nuc_sub"
                ]
                for r_name in nuc_sub_reactions:
                    if checker.check_reaction(r_name, rsmi):
                        functionalization_steps.append(("nucleophilic_substitution", depth))
                        findings_json["atomic_checks"]["named_reactions"].append(r_name)
                        reaction_detected = True
                        print(f"Detected nucleophilic substitution at depth {depth}: {rsmi}")
                        break

            # Halogenation
            if not reaction_detected:
                halogenation_reactions = [
                    "Aromatic chlorination",
                    "Aromatic bromination",
                    "Aromatic fluorination",
                    "Aromatic iodination",
                    "Chlorination",
                    "Bromination",
                    "Fluorination",
                    "Iodination"
                ]
                for r_name in halogenation_reactions:
                    if checker.check_reaction(r_name, rsmi):
                        functionalization_steps.append(("halogenation", depth))
                        findings_json["atomic_checks"]["named_reactions"].append(r_name)
                        reaction_detected = True
                        print(f"Detected halogenation at depth {depth}: {rsmi}")
                        break

            # Oxidation reactions
            if not reaction_detected:
                oxidation_reactions = [
                    "Oxidation of aldehydes to carboxylic acids",
                    "Oxidation of alcohol to carboxylic acid",
                    "Oxidation of ketone to carboxylic acid",
                    "Oxidation or Dehydrogenation of Alcohols to Aldehydes and Ketones",
                    "Oxidation of alkene to carboxylic acid",
                    "Oxidation of nitrile to carboxylic acid"
                ]
                for r_name in oxidation_reactions:
                    if checker.check_reaction(r_name, rsmi):
                        functionalization_steps.append(("oxidation", depth))
                        findings_json["atomic_checks"]["named_reactions"].append(r_name)
                        reaction_detected = True
                        print(f"Detected oxidation reaction at depth {depth}: {rsmi}")
                        break

            # Reduction reactions
            if not reaction_detected:
                reduction_reactions = [
                    "Reduction of aldehydes and ketones to alcohols",
                    "Reduction of ester to primary alcohol",
                    "Reduction of carboxylic acid to primary alcohol",
                    "Reduction of nitrile to amine",
                    "Reduction of primary amides to amines",
                    "Reduction of nitro groups to amines"
                ]
                for r_name in reduction_reactions:
                    if checker.check_reaction(r_name, rsmi):
                        functionalization_steps.append(("reduction", depth))
                        findings_json["atomic_checks"]["named_reactions"].append(r_name)
                        reaction_detected = True
                        print(f"Detected reduction reaction at depth {depth}: {rsmi}")
                        break

            # Alkylation reactions
            if not reaction_detected:
                alkylation_reactions = [
                    "Alkylation of amines",
                    "N-alkylation of primary amines with alkyl halides",
                    "N-alkylation of secondary amines with alkyl halides",
                    "S-alkylation of thiols",
                    "O-alkylation of carboxylic acids with diazo compounds"
                ]
                for r_name in alkylation_reactions:
                    if checker.check_reaction(r_name, rsmi):
                        functionalization_steps.append(("alkylation", depth))
                        findings_json["atomic_checks"]["named_reactions"].append(r_name)
                        reaction_detected = True
                        print(f"Detected alkylation reaction at depth {depth}: {rsmi}")
                        break

            # Protection reactions
            if not reaction_detected:
                protection_reactions = [
                    "Alcohol protection with silyl ethers",
                    "Boc amine protection",
                    "Protection of carboxylic acid"
                ]
                for r_name in protection_reactions:
                    if checker.check_reaction(r_name, rsmi):
                        functionalization_steps.append(("protection", depth))
                        findings_json["atomic_checks"]["named_reactions"].append(r_name)
                        reaction_detected = True
                        print(f"Detected protection reaction at depth {depth}: {rsmi}")
                        break

            # Deprotection reactions
            if not reaction_detected:
                deprotection_reactions = [
                    "Alcohol deprotection from silyl ethers",
                    "Boc amine deprotection",
                    "Deprotection of carboxylic acid",
                    "Ester saponification (methyl deprotection)",
                    "Ester saponification (alkyl deprotection)"
                ]
                for r_name in deprotection_reactions:
                    if checker.check_reaction(r_name, rsmi):
                        functionalization_steps.append(("deprotection", depth))
                        findings_json["atomic_checks"]["named_reactions"].append(r_name)
                        reaction_detected = True
                        print(f"Detected deprotection reaction at depth {depth}: {rsmi}")
                        break

        # Traverse children (in retrosynthetic direction)
        for child in node.get("children", []):
            # New logic: depth increases only when traversing from chemical to reaction node
            # Depth remains the same when traversing from reaction to chemical node
            new_depth = depth + 1 if node["type"] != "reaction" else depth
            dfs_traverse(child, new_depth)

    # Start traversal from the root
    dfs_traverse(route)

    # Sort functionalization steps by depth to check if they're sequential
    functionalization_steps.sort(key=lambda x: x[1])

    # Check if we have at least 4 different functionalization steps
    unique_steps = set(step[0] for step in functionalization_steps)

    if len(unique_steps) >= 4:
        # Add to findings_json if this constraint is met
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "unique_functionalization_step",
                "operator": ">=",
                "value": 4
            }
        })
        print(f"Found at least 4 unique functionalization steps: {unique_steps}")

        # Check if the steps are generally sequential
        # We don't require strictly increasing depths, just a logical progression
        depths = [step[1] for step in functionalization_steps]

        # Check if depths are non-decreasing or have at most one decrease
        decreases = sum(1 for i in range(len(depths) - 1) if depths[i] > depths[i + 1])
        is_sequential = decreases <= 1

        if is_sequential:
            # Add to findings_json if this constraint is met
            findings_json["structural_constraints"].append({
                "type": "sequence",
                "details": {
                    "target": "functionalization_steps",
                    "condition": "The sequence of reaction depths must be non-decreasing, allowing for at most one inversion."
                }
            })
            print(f"Found sequential functionalization with steps: {unique_steps}")
            print(f"Steps in order: {functionalization_steps}")
            result = True
        else:
            print(f"Found {len(unique_steps)} functionalization steps, but they are not sequential")
    else:
        print(f"Found only {len(unique_steps)} different functionalization steps, need at least 4")

    return result, findings_json

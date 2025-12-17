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


# Refactored lists of reactions
EARLY_STAGE_HALOGENATIONS = [
    "Aromatic chlorination", "Chlorination", "Alkyl chlorides from alcohols",
    "Aromatic fluorination", "Aromatic bromination", "Aromatic iodination",
    "Fluorination", "Bromination", "Iodination"
]

MIDDLE_STAGE_CN_FORMATIONS = [
    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
    "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
    "Reductive amination with aldehyde", "Reductive amination with ketone",
    "Reductive amination with alcohol",
    "N-alkylation of primary amines with alkyl halides",
    "N-alkylation of secondary amines with alkyl halides",
    "Goldberg coupling", "Ullmann-Goldberg Substitution amine"
]

LATE_STAGE_CO_FORMATIONS = [
    "Esterification of Carboxylic Acids", "Williamson Ether Synthesis",
    "Oxidation or Dehydrogenation of Alcohols to Aldehydes and Ketones",
    "Oxidation of aldehydes to carboxylic acids",
    "Oxidation of alcohol to carboxylic acid",
    "Reduction of aldehydes and ketones to alcohols",
    "Reduction of ester to primary alcohol",
    "Reduction of ketone to secondary alcohol",
    "Schotten-Baumann to ester", "Acetic anhydride and alcohol to ester"
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a sequential strategy using specific reaction types at different synthesis stages:
    early-stage halogenations (from EARLY_STAGE_HALOGENATIONS), middle-stage C-N couplings
    (from MIDDLE_STAGE_CN_FORMATIONS), and late-stage C-O transformations (from LATE_STAGE_CO_FORMATIONS).
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

    # Track bond formations by depth
    bond_formations = {}
    max_depth = 0
    result = False

    # First pass to determine maximum depth
    def get_max_depth(node, current_depth=0):
        nonlocal max_depth
        max_depth = max(max_depth, current_depth)
        for child in node.get("children", []):
            get_max_depth(child, current_depth + 1)

    get_max_depth(route)
    print(f"Maximum depth in route: {max_depth}")

    def dfs_traverse(node, current_depth=0):
        nonlocal findings_json, bond_formations
        print(f"Traversing node type: {node['type']} at depth {current_depth}")

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["mapped_reaction_smiles"]

            # Calculate synthesis stage (early, middle, late)
            # In retrosynthesis: high depth = early stage, low depth = late stage
            synthesis_stage = None
            if max_depth >= 2:  # If we have enough depth for three stages
                if current_depth >= max_depth - 1:  # Early stage (deepest or second deepest)
                    synthesis_stage = "early"
                elif current_depth == 1:  # Late stage (closest to target)
                    synthesis_stage = "late"
                else:  # Middle stage
                    synthesis_stage = "middle"
            else:  # Not enough depth, use relative positioning
                if current_depth == max_depth:  # Deepest = early
                    synthesis_stage = "early"
                elif current_depth == 0:  # Root = late
                    synthesis_stage = "late"
                else:  # Anything in between = middle
                    synthesis_stage = "middle"

            print(
                f"Examining reaction at depth {current_depth}, stage: {synthesis_stage}, RSMI: {rsmi}"
            )

            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for specific bond formations based on synthesis stage
            if synthesis_stage == "early":
                # Check for halogenation reactions
                for rxn in EARLY_STAGE_HALOGENATIONS:
                    if checker.check_reaction(rxn, rsmi):
                        print(f"Found halogenation reaction '{rxn}' at depth {current_depth} (early stage)")
                        bond_formations["early"] = "C-X"
                        findings_json["atomic_checks"]["named_reactions"].append(rxn)
                        break

            elif synthesis_stage == "middle":
                # Check for C-N bond formation reactions
                for rxn in MIDDLE_STAGE_CN_FORMATIONS:
                    if checker.check_reaction(rxn, rsmi):
                        print(
                            f"Found C-N bond formation reaction '{rxn}' at depth {current_depth} (middle stage)"
                        )
                        bond_formations["middle"] = "C-N"
                        findings_json["atomic_checks"]["named_reactions"].append(rxn)
                        break

            elif synthesis_stage == "late":
                # Check for C-O bond formation reactions
                for rxn in LATE_STAGE_CO_FORMATIONS:
                    if checker.check_reaction(rxn, rsmi):
                        print(
                            f"Found C-O bond formation reaction '{rxn}' at depth {current_depth} (late stage)"
                        )
                        bond_formations["late"] = "C-O"
                        findings_json["atomic_checks"]["named_reactions"].append(rxn)
                        break

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation:
            # Depth increases only when traversing from a chemical node to a reaction node.
            # Depth remains the same when traversing from a reaction node to a chemical node.
            if node['type'] == 'reaction':
                # If current node is a reaction, depth remains the same for its children (which are chemicals)
                dfs_traverse(child, current_depth)
            else:
                # If current node is a chemical, depth increases for its children (which are reactions)
                dfs_traverse(child, current_depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(f"Bond formations found: {bond_formations}")

    # Check if we have the correct sequence of bond formations
    expected_sequence = {"late": "C-O", "middle": "C-N", "early": "C-X"}
    all_found = True
    for stage, bond_type in expected_sequence.items():
        if stage not in bond_formations or bond_formations[stage] != bond_type:
            all_found = False
            break

    if all_found:
        print("Found sequential C-X bond formation strategy (C-X → C-N → C-O)")
        result = True
        # Add the structural constraint if the sequence is found
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "early_stage_halogenation",
                    "middle_stage_C-N_coupling",
                    "late_stage_C-O_transformation"
                ],
                "description": "The synthesis must contain a halogenation in an early stage, a C-N coupling in a middle stage, and a C-O transformation in a late stage. Stages are determined by reaction depth in the synthesis tree."
            }
        })

    return result, findings_json

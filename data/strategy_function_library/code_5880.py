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

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a one-pot Curtius rearrangement, identified by a direct 'Carboxylic acid to carbamate conversion' reaction type. This function also attempts to find a multi-step sequence of intermediates (Carboxylic Acid -> Azide/Isocyanate -> Carbamate) across a full synthesis path.
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

    # Track if we've found the Curtius rearrangement pattern
    has_curtius = False
    
    # Track if the direct reaction constraint was met
    direct_reaction_constraint_met = False

    # Track if the multi-step sequence constraint was met
    multi_step_sequence_constraint_met = False

    # Track molecules in the synthesis path to detect multi-step sequences
    synthesis_paths = []

    def dfs_traverse(node, path=None, depth=0):
        nonlocal has_curtius, findings_json, direct_reaction_constraint_met

        if path is None:
            path = []

        # For molecule nodes, add to the current path
        if node["type"] == "mol":
            current_path = path + [node["smiles"]]

            # Add this path to our collection of synthesis paths
            if node.get("in_stock", False):
                synthesis_paths.append(current_path)
        else:
            current_path = path

        # For reaction nodes, check for Curtius rearrangement characteristics
        if node["type"] == "reaction":
            try:
                # Extract reactants and product
                rsmi = node["metadata"]["mapped_reaction_smiles"]

                # Check for specific reaction types that might be part of Curtius
                if checker.check_reaction("Carboxylic acid to carbamate conversion", rsmi):
                    print(f"Detected Curtius rearrangement reaction: {rsmi}")
                    has_curtius = True
                    if "Carboxylic acid to carbamate conversion" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("Carboxylic acid to carbamate conversion")
                    
                    # Mark the direct reaction constraint as met
                    direct_reaction_constraint_met = True

            except Exception as e:
                print(f"Error analyzing reaction: {e}")

        # Determine the depth for the next recursive call
        next_depth = depth
        if node["type"] != "reaction": # This means it's a chemical node (mol)
            next_depth = depth + 1

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, current_path, next_depth)

    # Start traversal
    dfs_traverse(route)

    # Analyze synthesis paths for multi-step Curtius patterns
    for path in synthesis_paths:
        # Check if the path contains a sequence of:
        # 1. A molecule with carboxylic acid
        # 2. A molecule with azide or isocyanate (intermediate)
        # 3. A molecule with carbamic ester (final product)

        has_carboxylic_acid_step = False
        has_intermediate_step = False
        has_carbamic_ester_step = False

        for mol_smiles in path:
            if checker.check_fg("Carboxylic acid", mol_smiles):
                has_carboxylic_acid_step = True
                if "Carboxylic acid" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Carboxylic acid")

            if checker.check_fg("Azide", mol_smiles):
                if "Azide" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Azide")
                if has_carboxylic_acid_step:
                    has_intermediate_step = True
            
            if checker.check_fg("Isocyanate", mol_smiles):
                if "Isocyanate" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Isocyanate")
                if has_carboxylic_acid_step:
                    has_intermediate_step = True

            if checker.check_fg("Carbamic ester", mol_smiles):
                if "Carbamic ester" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Carbamic ester")
                if has_carboxylic_acid_step and has_intermediate_step:
                    has_carbamic_ester_step = True

        if has_carboxylic_acid_step and has_intermediate_step and has_carbamic_ester_step:
            print(f"Detected Curtius rearrangement sequence in synthesis path")
            has_curtius = True
            multi_step_sequence_constraint_met = True

    # Record structural constraints if met
    if direct_reaction_constraint_met:
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "Carboxylic acid to carbamate conversion"
                ],
                "min_occurrences": 1,
                "description": "The route must contain at least one instance of the 'Carboxylic acid to carbamate conversion' reaction."
            }
        })
    
    if multi_step_sequence_constraint_met:
        findings_json["structural_constraints"].append({
            "type": "sequence",
            "details": {
                "ordered_events": [
                    "Carboxylic acid",
                    [
                        "Azide",
                        "Isocyanate"
                    ],
                    "Carbamic ester"
                ],
                "event_type": "functional_group",
                "description": "The route must contain a sequence of molecules where a Carboxylic acid appears, followed by an Azide or Isocyanate, followed by a Carbamic ester."
            }
        })

    return has_curtius, findings_json

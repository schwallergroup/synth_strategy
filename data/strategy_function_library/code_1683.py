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
    This function detects if the synthesis involves a specific sequence of
    heteroatom bond formations: C-O → C-B → C-C → C-N → C-N
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

    # Define patterns for each bond type
    co_pattern = "C-O formation"  # C-O bond formation
    cb_pattern = "C-B formation"  # Borylation
    cc_pattern = "C-C formation"  # C-C bond formation
    cn1_pattern = "C-N formation 1"  # First C-N bond formation
    cn2_pattern = "C-N formation 2"  # Second C-N bond formation

    # Track bond formations in order (from early to late)
    bond_formation_sequence = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                print(f"Analyzing reaction at depth {depth}: {rsmi}")

                # Extract reactants and product for manual analysis
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for C-O bond formation
                co_formation = False
                co_reactions = [
                    "Esterification of Carboxylic Acids",
                    "Williamson Ether Synthesis",
                    "Alcohol protection with silyl ethers",
                    "O-alkylation of carboxylic acids with diazo compounds",
                    "Oxidative esterification of primary alcohols",
                    "Acetic anhydride and alcohol to ester",
                    "Mitsunobu esterification",
                    "Mitsunobu aryl ether"
                ]
                for r_name in co_reactions:
                    if checker.check_reaction(r_name, rsmi):
                        co_formation = True
                        if r_name not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(r_name)

                if co_formation:
                    bond_formation_sequence.append((depth, co_pattern))
                    print(f"C-O bond formation detected at depth {depth}")

                # Check for C-B bond formation
                cb_formation = False
                cb_reactions = [
                    "Preparation of boronic acids",
                    "Preparation of boronic acids without boronic ether",
                    "Preparation of boronic acids from trifluoroborates",
                    "Preparation of boronic esters",
                    "Synthesis of boronic acids"
                ]
                for r_name in cb_reactions:
                    if checker.check_reaction(r_name, rsmi):
                        cb_formation = True
                        if r_name not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(r_name)

                if cb_formation:
                    bond_formation_sequence.append((depth, cb_pattern))
                    print(f"C-B bond formation detected at depth {depth}")

                # Check for C-C bond formation
                cc_formation = False
                cc_reactions = [
                    "Suzuki coupling with boronic acids",
                    "Suzuki coupling with boronic esters",
                    "Suzuki coupling with boronic acids OTf",
                    "Suzuki coupling with boronic esters OTf",
                    "Suzuki coupling with sulfonic esters",
                    "Negishi coupling",
                    "Heck terminal vinyl",
                    "Stille reaction_aryl",
                    "Grignard from aldehyde to alcohol",
                    "Grignard from ketone to alcohol",
                    "Wittig reaction with triphenylphosphorane",
                    "Aldol condensation",
                    "Suzuki"
                ]
                for r_name in cc_reactions:
                    if checker.check_reaction(r_name, rsmi):
                        cc_formation = True
                        if r_name not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(r_name)

                if cc_formation:
                    bond_formation_sequence.append((depth, cc_pattern))
                    print(f"C-C bond formation detected at depth {depth}")

                # Check for first C-N bond formation
                cn1_formation = False
                cn1_reactions = [
                    "Reductive amination with aldehyde",
                    "Reductive amination with ketone",
                    "Reductive amination with alcohol",
                    "N-alkylation of primary amines with alkyl halides",
                    "N-alkylation of secondary amines with alkyl halides",
                    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
                    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
                    "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
                    "Alkylation of amines",
                    "reductive amination"
                ]
                for r_name in cn1_reactions:
                    if checker.check_reaction(r_name, rsmi):
                        cn1_formation = True
                        if r_name not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(r_name)

                if cn1_formation:
                    bond_formation_sequence.append((depth, cn1_pattern))
                    print(f"First C-N bond formation detected at depth {depth}")

                # Check for second C-N bond formation (amide formation)
                cn2_formation = False
                cn2_reactions = [
                    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                    "Acyl chloride with secondary amine to amide",
                    "Carboxylic acid with primary amine to amide",
                    "Ester with primary amine to amide",
                    "Ester with secondary amine to amide",
                    "Acyl chloride with ammonia to amide",
                    "Acylation of primary amines",
                    "Acylation of secondary amines",
                    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
                    "Schotten-Baumann_amide"
                ]
                for r_name in cn2_reactions:
                    if checker.check_reaction(r_name, rsmi):
                        cn2_formation = True
                        if r_name not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(r_name)

                if cn2_formation:
                    bond_formation_sequence.append((depth, cn2_pattern))
                    print(f"Second C-N bond formation (amide) detected at depth {depth}")

        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            dfs_traverse(child, new_depth)

    dfs_traverse(route)

    # Sort by depth (high to low, as high depth is early in synthesis)
    bond_formation_sequence.sort(key=lambda x: x[0], reverse=True)
    print(f"Bond formation sequence (sorted by depth): {bond_formation_sequence}")

    # Extract just the patterns in sequence
    patterns = [p for _, p in bond_formation_sequence]
    print(f"Patterns in sequence: {patterns}")

    # Check if the sequence matches our expected pattern
    expected_sequence = [co_pattern, cb_pattern, cc_pattern, cn1_pattern, cn2_pattern]
    print(f"Expected sequence: {expected_sequence}")

    # Check if the patterns appear in the correct order (not necessarily consecutive)
    last_found_idx = -1
    sequence_valid = True

    # Track which patterns we've found
    found_expected_patterns = []

    for pattern in expected_sequence:
        if pattern in patterns:
            current_idx = patterns.index(pattern)
            found_expected_patterns.append(pattern)
            if current_idx <= last_found_idx:
                print(
                    f"Heteroatom bond formation sequence not in expected order: {pattern} found at position {current_idx}, after position {last_found_idx}"
                )
                sequence_valid = False
                break
            last_found_idx = current_idx

    # Check if we found at least 3 of the expected patterns in the correct order
    found_patterns = len(found_expected_patterns)
    print(f"Found {found_patterns} patterns from the expected sequence: {found_expected_patterns}")

    result = False
    if found_patterns >= 3 and sequence_valid:
        print(
            f"Sequential heteroatom bond formation detected with {found_patterns} steps in the expected order"
        )
        result = True
        # Add structural constraints to findings_json
        findings_json["structural_constraints"].append({
            "type": "sequence",
            "details": {
                "ordered_events": [
                    "C-O formation",
                    "C-B formation",
                    "C-C formation",
                    "C-N formation 1",
                    "C-N formation 2"
                ],
                "consecutive": False,
                "description": "Checks that specific bond formation events occur in a predefined order, starting from the earliest steps of the synthesis. The events are grouped by type (e.g., 'C-O formation' includes various esterification and etherification reactions)."
            }
        })
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "events_from_sequence",
                "operator": ">=",
                "value": 3,
                "description": "Requires that at least 3 of the 5 specified bond formation events are found in the route."
            }
        })

    if not result:
        print(f"Insufficient heteroatom bond formations in the expected sequence")
    return result, findings_json

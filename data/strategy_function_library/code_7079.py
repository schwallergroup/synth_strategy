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


TETRAZOLE_FORMATION_REACTIONS = [
    "Azide-nitrile click cycloaddition to tetrazole",
    "{tetrazole_connect_regioisomere_1}",
    "{tetrazole_connect_regioisomere_2}",
    "{tetrazole_terminal}",
]

TETRAZOLE_ALKYLATION_REACTIONS = [
    "{Mitsunobu_tetrazole_1}",
    "{Mitsunobu_tetrazole_2}",
    "{Mitsunobu_tetrazole_3}",
    "{Mitsunobu_tetrazole_4}",
]

def main(route) -> Tuple[bool, Dict]:
    """Detects a two-step strategy involving the formation of a tetrazole from a nitrile precursor, followed by a late-stage N-alkylation of the tetrazole ring. Tetrazole formation is confirmed using a specific list of nitrile-azide cycloaddition reaction checkers (TETRAZOLE_FORMATION_REACTIONS). The N-alkylation is confirmed using a specific list of tetrazole alkylation reaction checkers (TETRAZOLE_ALKYLATION_REACTIONS)."""
    findings_template = {
        "atomic_checks": {
            "named_reactions": [],
            "ring_systems": [],
            "functional_groups": []
        },
        "structural_constraints": []
    }
    findings_json = copy.deepcopy(findings_template)

    # Track if we found the key features
    found_tetrazole_formation = False
    found_nitrile_intermediate = False
    found_late_stage_alkylation = False
    max_depth = 0
    alkylation_depth = float("inf")

    def dfs_traverse(node, depth=0):
        nonlocal found_tetrazole_formation, found_nitrile_intermediate, found_late_stage_alkylation
        nonlocal max_depth, alkylation_depth, findings_json

        # Update max depth to determine early vs late stage
        max_depth = max(max_depth, depth)

        if node["type"] == "reaction":
            rsmi = node.get("metadata", {}).get("rsmi", "")
            if not rsmi:
                return

            reactants_part = rsmi.split(">")[0]
            product_part = rsmi.split(">")[-1]
            reactants = reactants_part.split(".")

            # Check for tetrazole formation from nitrile and azide using a specific list of reactions
            is_tetrazole_formation = False
            for r in TETRAZOLE_FORMATION_REACTIONS:
                if checker.check_reaction(r, rsmi):
                    is_tetrazole_formation = True
                    findings_json["atomic_checks"]["named_reactions"].append(r)
                    break

            if is_tetrazole_formation:
                found_tetrazole_formation = True
                # Confirm a nitrile was a reactant in this specific step
                nitrile_found_in_reactants = False
                for r_mol_smi in reactants:
                    if checker.check_fg("Nitrile", r_mol_smi):
                        nitrile_found_in_reactants = True
                        findings_json["atomic_checks"]["functional_groups"].append("Nitrile")
                        break
                if nitrile_found_in_reactants:
                    found_nitrile_intermediate = True

            # Check for N-alkylation of tetrazole
            has_tetrazole_reactant = False
            for r_mol_smi in reactants:
                if checker.check_ring("tetrazole", r_mol_smi):
                    has_tetrazole_reactant = True
                    findings_json["atomic_checks"]["ring_systems"].append("tetrazole")
                    break

            has_tetrazole_product = checker.check_ring("tetrazole", product_part)
            if has_tetrazole_product:
                if "tetrazole" not in findings_json["atomic_checks"]["ring_systems"]:
                    findings_json["atomic_checks"]["ring_systems"].append("tetrazole")

            if has_tetrazole_reactant and has_tetrazole_product:
                is_alkylation = False
                for r in TETRAZOLE_ALKYLATION_REACTIONS:
                    if checker.check_reaction(r, rsmi):
                        is_alkylation = True
                        findings_json["atomic_checks"]["named_reactions"].append(r)
                        break
                if is_alkylation:
                    found_late_stage_alkylation = True
                    alkylation_depth = min(alkylation_depth, depth)

        # Traverse children
        for child in node.get("children", []):
            new_depth = depth
            if node["type"] != "reaction": # Only increase depth if current node is not 'reaction'
                new_depth = depth + 1
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    # Determine if alkylation is late-stage (in first third of synthesis, where depth=1 is the last step)
    is_late_stage = alkylation_depth <= max_depth / 3

    strategy_detected = (
        found_tetrazole_formation
        and found_nitrile_intermediate
        and found_late_stage_alkylation
        and is_late_stage
    )

    # Populate structural constraints if detected
    if found_tetrazole_formation and found_nitrile_intermediate and found_late_stage_alkylation:
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "description": "The strategy requires both a tetrazole formation event from a nitrile precursor and a subsequent tetrazole alkylation event.",
                "targets": [
                    "tetrazole_formation_from_nitrile",
                    "tetrazole_alkylation"
                ]
            }
        })

    if found_late_stage_alkylation and is_late_stage:
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "description": "The tetrazole alkylation reaction must occur late in the synthesis, defined as within the first third of the forward synthesis steps (depth <= max_depth / 3).",
                "target": "tetrazole_alkylation",
                "position": "late_stage"
            }
        })

    return strategy_detected, findings_json

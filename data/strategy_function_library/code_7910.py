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
    This function detects a strategy involving heterocyclic ring formation (particularly quinolin-2-one)
    via intramolecular cyclization, followed by fragment coupling through N-alkylation and completed
    with late-stage SNAr reaction to incorporate a piperazine moiety.
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

    # Initialize tracking variables
    has_quinolinone_formation = False
    has_n_alkylation = False
    has_late_stage_snar_with_piperazine = False

    # Track the maximum depth for relative positioning
    max_depth = 0

    def dfs_traverse(node, current_depth=0):
        nonlocal has_quinolinone_formation, has_n_alkylation, has_late_stage_snar_with_piperazine, max_depth, findings_json

        # Update max depth
        max_depth = max(max_depth, current_depth)

        if node["type"] == "reaction":
            rsmi = node.get("metadata", {}).get("rsmi", "")

            if rsmi:
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                print(f"Analyzing reaction at depth {current_depth}: {rsmi}")

                # Check for quinolinone formation via cyclization (early stage)
                if current_depth >= 5:  # Early stage reaction based on stdout
                    # Check if product contains a lactam structure in a bicyclic system
                    if "c1c(=O)[nH]c2ccccc12" in product:
                        # Check if reactants don't have the quinolinone structure
                        reactant_has_quinolinone = any(
                            "c1c(=O)[nH]c2ccccc12" in r for r in reactants
                        )

                        if not reactant_has_quinolinone:
                            print(f"Detected quinolinone formation at depth {current_depth}")
                            has_quinolinone_formation = True
                            findings_json["atomic_checks"]["named_reactions"].append("quinolin-2-one_ring_formation")
                            findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "quinolin-2-one_ring_formation", "position": "early_stage", "min_depth": 5}})
                            findings_json["structural_constraints"].append({"type": "negation", "details": {"target": "quinolin-2-one", "scope": "reactants", "context": "Ensures de novo formation for the quinolin-2-one_ring_formation event"}})

                # Check for N-alkylation (middle stage)
                if 3 <= current_depth <= 7:  # Expanded middle stage range
                    # Check for N-alkylation reactions or benzyl halide + amine pattern
                    n_alkylation_reactions = [
                        "N-alkylation of primary amines with alkyl halides",
                        "N-alkylation of secondary amines with alkyl halides",
                        "Alkylation of amines"
                    ]
                    for r_name in n_alkylation_reactions:
                        if checker.check_reaction(r_name, rsmi):
                            print(f"Detected N-alkylation ({r_name}) at depth {current_depth}")
                            has_n_alkylation = True
                            findings_json["atomic_checks"]["named_reactions"].append(r_name)
                            findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "N-alkylation", "position": "middle_stage", "min_depth": 3, "max_depth": 7}})
                            break

                # Check for late-stage SNAr with piperazine (late stage)
                if current_depth <= 2:  # Late stage reaction based on stdout
                    # Check if piperazine is a reactant
                    piperazine_reactant_found = False
                    for r in reactants:
                        if checker.check_ring("piperazine", r):
                            piperazine_reactant_found = True
                            findings_json["atomic_checks"]["ring_systems"].append("piperazine")
                            break

                    if piperazine_reactant_found:
                        print(f"Detected piperazine as reactant at depth {current_depth}")

                        # Check for aromatic halide in reactants
                        has_aromatic_halide = False
                        for r in reactants:
                            if checker.check_fg("Aromatic halide", r):
                                has_aromatic_halide = True
                                findings_json["atomic_checks"]["functional_groups"].append("Aromatic halide")
                                break

                        if has_aromatic_halide:
                            findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"targets": ["piperazine", "Aromatic halide"], "scope": "reactants", "context": "Defines the late-stage_SNAr_with_piperazine event"}})
                            # Verify piperazine is in product
                            if checker.check_ring("piperazine", product):
                                print(
                                    f"Confirmed late-stage SNAr with piperazine at depth {current_depth}"
                                )
                                has_late_stage_snar_with_piperazine = True
                                findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "late-stage_SNAr_with_piperazine", "position": "late_stage", "max_depth": 2}})

        # Determine the next depth based on the current node's type
        next_depth = current_depth
        if node["type"] != "reaction": # This means it's a chemical node
            next_depth = current_depth + 1

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, next_depth)

    # Start traversal
    dfs_traverse(route)

    print(f"Strategy detection results:")
    print(f"- Quinolinone formation: {has_quinolinone_formation}")
    print(f"- N-alkylation: {has_n_alkylation}")
    print(f"- Late-stage SNAr with piperazine: {has_late_stage_snar_with_piperazine}")

    # Return True if all key elements of the strategy are present
    result = has_quinolinone_formation and has_n_alkylation and has_late_stage_snar_with_piperazine

    if result:
        findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"targets": ["quinolin-2-one_ring_formation", "N-alkylation", "late-stage_SNAr_with_piperazine"], "scope": "route"}})

    return result, findings_json

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


PYRAZOLE_FORMATION_REACTIONS = [
    "pyrazole",  # Generic name for pyrazole formation, e.g., Knorr synthesis
    "[3+2]-cycloaddition of hydrazone and alkyne",
    "[3+2]-cycloaddition of hydrazone and alkene",
    "[3+2]-cycloaddition of diazoalkane and alkyne",
    "[3+2]-cycloaddition of diazoalkane and alkene",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a strategy involving late-stage pyrazole ring formation via specific, known reaction types.
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

    has_late_pyrazole = False

    def dfs_traverse(node, depth=0):
        nonlocal has_late_pyrazole, findings_json

        if node["type"] == "reaction" and depth <= 3:  # Late stage (expanded to depth 3)
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check for pyrazole formation
                product_has_pyrazole = checker.check_ring("pyrazole", product_smiles)
                if product_has_pyrazole:
                    if "pyrazole" not in findings_json["atomic_checks"]["ring_systems"]:
                        findings_json["atomic_checks"]["ring_systems"].append("pyrazole")

                reactants_have_pyrazole = any(
                    checker.check_ring("pyrazole", r) for r in reactants_smiles
                )

                if product_has_pyrazole and not reactants_have_pyrazole:
                    # To avoid false positives, confirm the reaction is a known pyrazole formation type
                    is_known_formation_rxn = False
                    for name in PYRAZOLE_FORMATION_REACTIONS:
                        if checker.check_reaction(name, rsmi):
                            is_known_formation_rxn = True
                            if name not in findings_json["atomic_checks"]["named_reactions"]:
                                findings_json["atomic_checks"]["named_reactions"].append(name)
                            break

                    if is_known_formation_rxn:
                        print(
                            f"Found late-stage pyrazole formation at depth {depth} (reaction type match)"
                        )
                        print(f"Reaction SMILES: {rsmi}")
                        has_late_pyrazole = True
                        # Record structural constraints
                        if {"type": "co-occurrence", "details": {"scope": "single_reaction_step", "targets": ["ring_formation_of_pyrazole", "known_pyrazole_synthesis_reaction"], "description": "A pyrazole ring must be formed de novo (i.e., not present in reactants) and the reaction must match a known pyrazole synthesis type."}} not in findings_json["structural_constraints"]:
                            findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"scope": "single_reaction_step", "targets": ["ring_formation_of_pyrazole", "known_pyrazole_synthesis_reaction"], "description": "A pyrazole ring must be formed de novo (i.e., not present in reactants) and the reaction must match a known pyrazole synthesis type."}})
                        if {"type": "positional", "details": {"target": "pyrazole_ring_formation_by_known_synthesis", "position": "late_stage", "condition": "depth <= 3", "description": "The qualifying pyrazole formation reaction must occur within the first 4 steps from the final product."}} not in findings_json["structural_constraints"]:
                            findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "pyrazole_ring_formation_by_known_synthesis", "position": "late_stage", "condition": "depth <= 3", "description": "The qualifying pyrazole formation reaction must occur within the first 4 steps from the final product."}})

            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    print(f"Late-stage pyrazole formation: {has_late_pyrazole}")

    return has_late_pyrazole, findings_json

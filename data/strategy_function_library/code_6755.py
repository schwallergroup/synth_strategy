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
    This function detects a linear functional group interconversion sequence on a pyridine scaffold,
    featuring a halogen-to-carbonyl transformation pathway (Br → CN → COOH → COOMe)
    with the final step being esterification.
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

    # Track the transformation sequence and positions
    transformations = []

    def dfs_traverse(node, depth=0):
        nonlocal transformations, findings_json

        if node["type"] == "reaction":
            if "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Br to CN transformation (cyanation)
                br_to_cn_detected = False
                if (
                    any(
                        checker.check_fg("Aromatic halide", r) and checker.check_ring("pyridine", r)
                        for r in reactants
                    )
                    and checker.check_fg("Nitrile", product)
                    and checker.check_ring("pyridine", product)
                ):
                    transformations.append(("br_to_cn", depth))
                    br_to_cn_detected = True
                    if "br_to_cn" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("br_to_cn")
                    if "Aromatic halide" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Aromatic halide")
                    if "Nitrile" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Nitrile")
                    if "pyridine" not in findings_json["atomic_checks"]["ring_systems"]:
                        findings_json["atomic_checks"]["ring_systems"].append("pyridine")

                # CN to COOH transformation (hydrolysis)
                cn_to_cooh_detected = False
                if (
                    any(
                        checker.check_fg("Nitrile", r) and checker.check_ring("pyridine", r)
                        for r in reactants
                    )
                    and checker.check_fg("Carboxylic acid", product)
                    and checker.check_ring("pyridine", product)
                ):
                    transformations.append(("cn_to_cooh", depth))
                    cn_to_cooh_detected = True
                    if "cn_to_cooh" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("cn_to_cooh")
                    if "Nitrile" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Nitrile")
                    if "Carboxylic acid" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Carboxylic acid")
                    if "pyridine" not in findings_json["atomic_checks"]["ring_systems"]:
                        findings_json["atomic_checks"]["ring_systems"].append("pyridine")

                # COOH to COOMe transformation (esterification)
                cooh_to_coome_detected = False
                if (
                    any(
                        checker.check_fg("Carboxylic acid", r) and checker.check_ring("pyridine", r)
                        for r in reactants
                    )
                    and checker.check_fg("Ester", product)
                    and checker.check_ring("pyridine", product)
                ):
                    if checker.check_reaction("Esterification of Carboxylic Acids", rsmi):
                        transformations.append(("cooh_to_coome", depth))
                        cooh_to_coome_detected = True
                        if "cooh_to_coome" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("cooh_to_coome")
                        if "Esterification of Carboxylic Acids" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("Esterification of Carboxylic Acids")
                        if "Carboxylic acid" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Carboxylic acid")
                        if "Ester" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Ester")
                        if "pyridine" not in findings_json["atomic_checks"]["ring_systems"]:
                            findings_json["atomic_checks"]["ring_systems"].append("pyridine")

        # Traverse children
        for child in node.get("children", []):
            # New logic: depth increases only from chemical to reaction node
            # Depth remains same from reaction to chemical node
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, next is reaction, so depth increases
                new_depth = depth + 1
            # If current node is reaction, next is chemical, so depth remains the same
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    # Check if we have all the transformation reactions
    has_br_to_cn = any(t[0] == "br_to_cn" for t in transformations)
    has_cn_to_cooh = any(t[0] == "cn_to_cooh" for t in transformations)
    has_cooh_to_coome = any(t[0] == "cooh_to_coome" for t in transformations)

    # If all transformations are present, check their sequence
    correct_sequence = False
    if has_br_to_cn and has_cn_to_cooh and has_cooh_to_coome:
        # Record co-occurrence constraint
        if {"type": "co-occurrence", "details": {"targets": ["br_to_cn", "cn_to_cooh", "cooh_to_coome"]}} not in findings_json["structural_constraints"]:
            findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"targets": ["br_to_cn", "cn_to_cooh", "cooh_to_coome"]}})

        br_to_cn_depth = min(t[1] for t in transformations if t[0] == "br_to_cn")
        cn_to_cooh_depth = min(t[1] for t in transformations if t[0] == "cn_to_cooh")
        cooh_to_coome_depth = min(t[1] for t in transformations if t[0] == "cooh_to_coome")

        # In retrosynthetic order: cooh_to_coome -> cn_to_cooh -> br_to_cn
        # So depths should increase in that order (or be equal)
        if cooh_to_coome_depth <= cn_to_cooh_depth <= br_to_cn_depth:
            correct_sequence = True
            # Record sequence constraint
            if {"type": "sequence", "details": {"ordered_events": ["br_to_cn", "cn_to_cooh", "cooh_to_coome"], "direction": "forward_synthesis"}} not in findings_json["structural_constraints"]:
                findings_json["structural_constraints"].append({"type": "sequence", "details": {"ordered_events": ["br_to_cn", "cn_to_cooh", "cooh_to_coome"], "direction": "forward_synthesis"}})

    return correct_sequence, findings_json

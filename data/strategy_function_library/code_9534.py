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
    This function detects synthesis routes that involve the formation or modification of piperazine or piperidine rings.

    A relevant strategy involves:
    1. The presence of piperazine or late-stage piperidine scaffolds
    2. Reactions that form or modify these specific rings
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

    piperazine_present = False
    piperazine_reactions = False

    def dfs_traverse(node, depth=0):
        nonlocal piperazine_present, piperazine_reactions, findings_json

        if node["type"] == "mol" and node.get("smiles"):
            mol_smiles = node["smiles"]

            has_piperazine = checker.check_ring("piperazine", mol_smiles)
            if has_piperazine:
                piperazine_present = True
                if "piperazine" not in findings_json["atomic_checks"]["ring_systems"]:
                    findings_json["atomic_checks"]["ring_systems"].append("piperazine")

            has_piperidine = checker.check_ring("piperidine", mol_smiles)
            if has_piperidine and depth <= 2:  # Focus on late-stage molecules
                piperazine_present = True
                if "piperidine" not in findings_json["atomic_checks"]["ring_systems"]:
                    findings_json["atomic_checks"]["ring_systems"].append("piperidine")
                # Add structural constraint for late-stage piperidine
                if {"type": "positional", "details": {"target": "piperidine", "position": "depth <= 2"}} not in findings_json["structural_constraints"]:
                    findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "piperidine", "position": "depth <= 2"}})

        elif node["type"] == "reaction" and node.get("metadata", {}).get("rsmi"):
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            product_has_piperazine = checker.check_ring("piperazine", product)
            reactant_has_piperazine = any(checker.check_ring("piperazine", r) for r in reactants)

            if product_has_piperazine and not reactant_has_piperazine:
                piperazine_reactions = True
                if "ring_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append("ring_formation")
            elif product_has_piperazine and reactant_has_piperazine:
                piperazine_reactions = True
                if "ring_modification" not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append("ring_modification")

            product_has_piperidine = checker.check_ring("piperidine", product)
            reactant_has_piperidine = any(checker.check_ring("piperidine", r) for r in reactants)

            if (product_has_piperidine and not reactant_has_piperidine) or (
                product_has_piperidine and reactant_has_piperidine
            ):
                piperazine_reactions = True
                if product_has_piperidine and not reactant_has_piperidine:
                    if "ring_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("ring_formation")
                elif product_has_piperidine and reactant_has_piperidine:
                    if "ring_modification" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("ring_modification")

        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                # If current node is a reaction, depth remains the same for children (chemical nodes)
                dfs_traverse(child, depth)
            else:
                # If current node is a chemical, depth increases for children (reaction nodes)
                dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    result = piperazine_present and piperazine_reactions
    if result:
        # Add co-occurrence constraint if both conditions are met
        if {"type": "co-occurrence", "details": {"targets": ["presence_of_key_scaffold", "reaction_on_key_scaffold"]}} not in findings_json["structural_constraints"]:
            findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"targets": ["presence_of_key_scaffold", "reaction_on_key_scaffold"]}})

    return result, findings_json

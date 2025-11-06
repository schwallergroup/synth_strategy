from typing import Tuple, Dict, List
import copy
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
    Detects if the synthesis route involves a late-stage cyanation (introduction of a cyano group
    in the final or near-final step).
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

    cyano_added_in_late_stage = False

    def dfs_traverse(node, depth=0):
        nonlocal cyano_added_in_late_stage, findings_json

        if node["type"] == "reaction" and depth <= 1:
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_part = rsmi.split(">")[0]
                product_part = rsmi.split(">")[-1]

                cyano_in_product = checker.check_fg("Nitrile", product_part)

                if cyano_in_product:
                    if "Nitrile" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Nitrile")

                    reactants = reactants_part.split(".")

                    cyanide_sources = []
                    substrate_reactants = []

                    for reactant in reactants:
                        # Heuristic: Small molecules with nitrile are likely cyanide sources
                        if checker.check_fg("Nitrile", reactant) and len(reactant) <= 5:
                            cyanide_sources.append(reactant)
                        else:
                            substrate_reactants.append(reactant)

                    # Check if any substrate reactant doesn't have a nitrile group
                    for substrate in substrate_reactants:
                        if not checker.check_fg("Nitrile", substrate):
                            # This substrate gained a nitrile group. This is a cyanation.
                            if cyanide_sources or len(reactants) == 1:
                                cyano_added_in_late_stage = True
                                if "cyanation" not in findings_json["atomic_checks"]["named_reactions"]:
                                    findings_json["atomic_checks"]["named_reactions"].append("cyanation")
                                if {"type": "positional", "details": {"target": "cyanation", "position": "late_stage", "condition": "depth <= 1"}} not in findings_json["structural_constraints"]:
                                    findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "cyanation", "position": "late_stage", "condition": "depth <= 1"}})
                                break
            except Exception:
                pass

        # Traverse children
        for child in node.get("children", []):
            if cyano_added_in_late_stage:
                return
            
            # New depth calculation logic
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            
            dfs_traverse(child, new_depth)

    # Start traversal from root
    dfs_traverse(route)
    return cyano_added_in_late_stage, findings_json

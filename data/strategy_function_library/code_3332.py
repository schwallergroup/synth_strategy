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


CORE_SCAFFOLDS_OF_INTEREST = ["indole"]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if a synthesis maintains a specific core scaffold throughout the route. This strategy is confirmed if a scaffold from the CORE_SCAFFOLDS_OF_INTEREST list is present in the final product and is never destroyed in any synthetic step. The core is allowed to be formed during the synthesis.
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

    core_present_at_depths = set()
    core_preserved = True

    def dfs_traverse(node, depth=0):
        nonlocal core_preserved, findings_json

        if node["type"] == "mol" and "smiles" in node:
            mol_smiles = node["smiles"]
            has_core = False
            for core in CORE_SCAFFOLDS_OF_INTEREST:
                if checker.check_ring(core, mol_smiles):
                    has_core = True
                    if core not in findings_json["atomic_checks"]["ring_systems"]:
                        findings_json["atomic_checks"]["ring_systems"].append(core)
            if has_core:
                core_present_at_depths.add(depth)

        elif node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                product = rsmi.split(">")[-1]
                reactants = rsmi.split(">")[0].split(".")

                product_has_core = False
                for core in CORE_SCAFFOLDS_OF_INTEREST:
                    if checker.check_ring(core, product):
                        product_has_core = True
                        if core not in findings_json["atomic_checks"]["ring_systems"]:
                            findings_json["atomic_checks"]["ring_systems"].append(core)

                reactants_with_core = []
                reactant_has_core = False
                for r in reactants:
                    for core in CORE_SCAFFOLDS_OF_INTEREST:
                        if checker.check_ring(core, r):
                            reactants_with_core.append(r)
                            reactant_has_core = True
                            if core not in findings_json["atomic_checks"]["ring_systems"]:
                                findings_json["atomic_checks"]["ring_systems"].append(core)
                            break # Found core in this reactant, move to next reactant

                if not product_has_core and reactant_has_core:
                    core_preserved = False
                    # This implies ring destruction
                    if "ring_destruction" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("ring_destruction")
                    # Add the negation constraint if it's the first time this condition is met
                    if {"type": "negation", "details": {"target": "ring_destruction", "scope": "indole"}} not in findings_json["structural_constraints"]:
                        findings_json["structural_constraints"].append({"type": "negation", "details": {"target": "ring_destruction", "scope": "indole"}})

            except Exception:
                pass

        for child in node.get("children", []):
            # New depth calculation logic
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    dfs_traverse(route)

    result = len(core_present_at_depths) >= 2 and 0 in core_present_at_depths and core_preserved

    # Populate structural constraints based on final result
    if 0 in core_present_at_depths:
        if {"type": "positional", "details": {"target": "indole", "position": "last_stage"}} not in findings_json["structural_constraints"]:
            findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "indole", "position": "last_stage"}})

    if len(core_present_at_depths) >= 2:
        if {"type": "count", "details": {"target": "stages_with_indole", "operator": ">=", "value": 2}} not in findings_json["structural_constraints"]:
            findings_json["structural_constraints"].append({"type": "count", "details": {"target": "stages_with_indole", "operator": ">=", "value": 2}})

    return result, findings_json

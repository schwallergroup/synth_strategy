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
    Detects a specific convergent strategy involving thiazole and furan cores
    with a nitrile-containing fragment in the final step.
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

    # Track key features
    has_thiazole_core = False
    has_furan_core = False
    has_nitrile_fragment = False
    has_convergent_final = False

    # Track the final product
    final_product_smiles = None

    def dfs_traverse(node, depth=0):
        nonlocal has_thiazole_core, has_furan_core, has_nitrile_fragment, has_convergent_final, findings_json

        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            if depth == 0:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_str = rsmi.split(">")[0]
                product_str = rsmi.split(">")[-1]
                reactants_smiles = reactants_str.split(".")

                product_has_thiazole = checker.check_ring("thiazole", product_str)
                product_has_furan = checker.check_ring("furan", product_str)
                product_has_nitrile = checker.check_fg("Nitrile", product_str)

                if product_has_thiazole and product_has_furan and product_has_nitrile:
                    has_thiazole_core = True
                    has_furan_core = True
                    has_nitrile_fragment = True
                    if "thiazole" not in findings_json["atomic_checks"]["ring_systems"]:
                        findings_json["atomic_checks"]["ring_systems"].append("thiazole")
                    if "furan" not in findings_json["atomic_checks"]["ring_systems"]:
                        findings_json["atomic_checks"]["ring_systems"].append("furan")
                    if "Nitrile" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Nitrile")
                    
                    # Structural constraint: all_strategy_checks at last_stage
                    findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "all_strategy_checks", "position": "last_stage"}})
                    # Structural constraint: co-occurrence of thiazole, furan, Nitrile
                    findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"targets": ["thiazole", "furan", "Nitrile"]}})

                    if len(reactants_smiles) >= 2:
                        # Structural constraint: count reactants >= 2
                        findings_json["structural_constraints"].append({"type": "count", "details": {"target": "reactants", "operator": ">=", "value": 2}})

                        thiazole_reactants = 0
                        furan_reactants = 0
                        nitrile_reactants = 0

                        for r_smiles in reactants_smiles:
                            if checker.check_ring("thiazole", r_smiles):
                                thiazole_reactants += 1
                            if checker.check_ring("furan", r_smiles):
                                furan_reactants += 1
                            if checker.check_fg("Nitrile", r_smiles):
                                nitrile_reactants += 1

                        if (
                            (thiazole_reactants > 0 and furan_reactants > 0)
                            or (thiazole_reactants > 0 and nitrile_reactants > 0)
                            or (furan_reactants > 0 and nitrile_reactants > 0)
                        ):
                            has_convergent_final = True
                            # Structural constraint: distinct_key_fragments_in_reactants >= 2
                            findings_json["structural_constraints"].append({"type": "count", "details": {"target": "distinct_key_fragments_in_reactants", "operator": ">=", "value": 2}})

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

    # This block is now inert as final_product_smiles is never set.
    # The original code had this block, but final_product_smiles is always None.
    # Keeping it for consistency with original logic flow, though it won't trigger.
    if final_product_smiles:
        if checker.check_ring("thiazole", final_product_smiles) and checker.check_ring(
            "furan", final_product_smiles
        ):
            has_thiazole_core = True
            has_furan_core = True
            if "thiazole" not in findings_json["atomic_checks"]["ring_systems"]:
                findings_json["atomic_checks"]["ring_systems"].append("thiazole")
            if "furan" not in findings_json["atomic_checks"]["ring_systems"]:
                findings_json["atomic_checks"]["ring_systems"].append("furan")

        if checker.check_fg("Nitrile", final_product_smiles):
            has_nitrile_fragment = True
            if "Nitrile" not in findings_json["atomic_checks"]["functional_groups"]:
                findings_json["atomic_checks"]["functional_groups"].append("Nitrile")

        if (
            checker.check_ring("thiazole", final_product_smiles)
            and checker.check_ring("furan", final_product_smiles)
            and checker.check_fg("Nitrile", final_product_smiles)
        ):
            has_convergent_final = True
            # Structural constraint: all_strategy_checks at last_stage
            findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "all_strategy_checks", "position": "last_stage"}})
            # Structural constraint: co-occurrence of thiazole, furan, Nitrile
            findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"targets": ["thiazole", "furan", "Nitrile"]}})

    # The print statements are for debugging and not part of the core logic, removed for production code.
    # print(f"Has thiazole core: {has_thiazole_core}")
    # print(f"Has furan core: {has_furan_core}")
    # print(f"Has nitrile fragment: {has_nitrile_fragment}")
    # print(f"Has convergent final: {has_convergent_final}")

    result = has_thiazole_core and has_furan_core and has_nitrile_fragment and has_convergent_final
    return result, findings_json
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


BENZOFUSED_HETEROCYCLES_OF_INTEREST = [
    "benzimidazole",
    "quinazoline",
    "benzoxazole",
    "benzothiazole",
]

BENZOFUSED_HETEROCYCLE_FORMATION_REACTIONS = [
    "benzimidazole_derivatives_carboxylic-acid/ester",
    "benzimidazole_derivatives_aldehyde",
    "Niementowski_quinazoline",
    "benzoxazole",
    "benzoxazole_arom-aldehyde",
    "benzoxazole_carboxylic-acid",
    "benzothiazole",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a synthetic strategy focused on specific benzofused heterocyclic cores.
    The strategy is flagged if either of the following conditions is met:
    1. A reaction from the BENZOFUSED_HETEROCYCLE_FORMATION_REACTIONS list is used to construct a core.
    2. A core from the BENZOFUSED_HETEROCYCLES_OF_INTEREST list is present in the final product and is preserved in at least half of the synthetic intermediates.
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

    # Track if benzimidazole or related cores are present in molecules
    molecules_with_core = 0
    total_molecules = 0

    # Track if we're building the core in the synthesis
    building_core = False

    # Track if we have the core in the final product
    final_product_has_core = False

    def dfs_traverse(node, depth=0):
        nonlocal molecules_with_core, total_molecules, building_core, final_product_has_core, findings_json

        if node["type"] == "mol" and not node.get("in_stock", False):
            total_molecules += 1
            mol_smiles = node["smiles"]

            # Check for specific benzofused heterocyclic cores
            has_target_core = False
            for core in BENZOFUSED_HETEROCYCLES_OF_INTEREST:
                if checker.check_ring(core, mol_smiles):
                    has_target_core = True
                    if core not in findings_json["atomic_checks"]["ring_systems"]:
                        findings_json["atomic_checks"]["ring_systems"].append(core)

            if has_target_core:
                molecules_with_core += 1
                print(f"Molecule with target heterocyclic core found: {mol_smiles}")

                # If this is the root node (final product), mark it
                if depth == 0:
                    final_product_has_core = True
            else:
                print(f"Molecule without target heterocyclic core found: {mol_smiles}")

        # Check if this is a reaction that forms one of the target cores
        elif node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rxn_smiles = node["metadata"]["mapped_reaction_smiles"]

            # Check if this reaction is in the list of formation reactions
            for rxn in BENZOFUSED_HETEROCYCLE_FORMATION_REACTIONS:
                if checker.check_reaction(rxn, rxn_smiles):
                    building_core = True
                    if rxn not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append(rxn)
                    print(f"Found reaction forming target core: {rxn_smiles}")
                    break # Found one, no need to check others

        # If the node contains children, traverse each recursively
        for child in node.get("children", []):
            # Determine the new depth based on the current node's type
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            
            dfs_traverse(child, new_depth)

    # Start traversal from the root
    dfs_traverse(route)

    # Strategy is valid if:
    # 1. We're building the core during synthesis, OR
    # 2. The final product has the core and a significant portion of intermediates preserve it
    core_preservation_ratio = molecules_with_core / total_molecules if total_molecules > 0 else 0

    print(f"Core building: {building_core}")
    print(f"Final product has core: {final_product_has_core}")
    print(
        f"Core preservation ratio: {core_preservation_ratio} ({molecules_with_core}/{total_molecules})"
    )

    result = building_core or (final_product_has_core and core_preservation_ratio >= 0.5)

    # Record structural constraints if met
    if building_core:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "BENZOFUSED_HETEROCYCLE_FORMATION_REACTIONS",
                "operator": ">=",
                "value": 1
            }
        })
    
    if final_product_has_core:
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "BENZOFUSED_HETEROCYCLES_OF_INTEREST",
                "position": "last_stage"
            }
        })

    if core_preservation_ratio >= 0.5:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "core_preservation_ratio",
                "operator": ">=",
                "value": 0.5
            }
        })

    if final_product_has_core and core_preservation_ratio >= 0.5:
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "core_in_final_product",
                    "high_core_preservation_ratio"
                ]
            }
        })

    return result, findings_json

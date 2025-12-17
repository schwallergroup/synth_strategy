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


HETEROCYCLES_OF_INTEREST = [
    "pyrimidine",
    "pyrazole",
    "pyridine",
    "pyrrole",
    "imidazole",
    "oxazole",
    "thiazole",
    "triazole",
    "tetrazole",
    "furan",
    "thiophene",
    "isoxazole",
    "isothiazole",
    "oxadiazole",
    "thiadiazole",
]

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects the incorporation of multiple heterocyclic fragments
    in the synthesis.
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

    # Track which heterocycles have been found
    found_heterocycles = set()

    def dfs_traverse(node, depth=0):
        nonlocal found_heterocycles, findings_json
        if node["type"] == "reaction":
            if "mapped_reaction_smiles" in node.get("metadata", {}):
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                try:
                    # Extract reactants and product
                    reactants_part = rsmi.split(">")[0]
                    product_part = rsmi.split(">")[-1]

                    # Check reactants for heterocycles
                    for reactant_smiles in reactants_part.split("."):
                        for heterocycle in HETEROCYCLES_OF_INTEREST:
                            if heterocycle not in found_heterocycles and checker.check_ring(
                                heterocycle, reactant_smiles
                            ):
                                found_heterocycles.add(heterocycle)
                                findings_json["atomic_checks"]["ring_systems"].append(heterocycle)
                                print(
                                    f"Detected {heterocycle} heterocycle in reactant at depth {depth}"
                                )

                    # Check product for heterocycles that might have been formed
                    for heterocycle in HETEROCYCLES_OF_INTEREST:
                        if heterocycle not in found_heterocycles and checker.check_ring(
                            heterocycle, product_part
                        ):
                            # Verify this heterocycle wasn't in any of the reactants
                            reactant_has_heterocycle = any(
                                checker.check_ring(heterocycle, r)
                                for r in reactants_part.split(".")
                            )
                            if not reactant_has_heterocycle:
                                found_heterocycles.add(heterocycle)
                                findings_json["atomic_checks"]["ring_systems"].append(heterocycle)
                                print(
                                    f"Detected {heterocycle} heterocycle formation in product at depth {depth}"
                                )

                except Exception as e:
                    print(f"Error processing reaction SMILES at depth {depth}: {e}")

        # Process children (depth-first)
        for child in node.get("children", []):
            # New logic: depth increases only when traversing from chemical to reaction
            # Depth remains the same when traversing from reaction to chemical
            new_depth = depth + 1 if node["type"] != "reaction" else depth
            dfs_traverse(child, new_depth)

    # Start traversal from the root
    dfs_traverse(route)

    # Determine the final result
    result = len(found_heterocycles) >= 2

    # Add structural constraint if met
    if result:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "unique_heterocycle_incorporation",
                "operator": ">=",
                "value": 2
            }
        })

    # Return True if at least 2 different heterocycles are incorporated
    print(f"Total unique heterocycles found: {len(found_heterocycles)}")
    print(f"Heterocycles found: {found_heterocycles}")
    return result, findings_json

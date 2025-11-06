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


# List of heterocycles to check
HETEROCYCLES_OF_INTEREST = [
    "furan",
    "pyrrole",
    "thiophene",
    "pyridine",
    "pyrazine",
    "pyrimidine",
    "pyridazine",
    "oxazole",
    "thiazole",
    "imidazole",
    "isoxazole",
    "isothiazole",
    "triazole",
    "tetrazole",
    "oxadiazole",
    "thiadiazole",
    "piperidine",
    "piperazine",
    "morpholine",
    "thiomorpholine",
    "indole",
    "benzofuran",
    "benzothiophene",
    "benzimidazole",
    "benzoxazole",
    "benzothiazole",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a linear synthesis strategy that includes a late-stage formation of a specific heterocycle. The list of heterocycles checked is defined in the module-level constant HETEROCYCLES_OF_INTEREST.
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

    # Track if we found the key features
    late_stage_heterocycle = False
    linear_synthesis = True
    max_depth = 0

    def dfs_traverse(node, depth=0):
        nonlocal late_stage_heterocycle, linear_synthesis, max_depth, findings_json

        max_depth = max(max_depth, depth)

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node.get("metadata", {}).get("rsmi", "")
            if not rsmi:
                return

            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check if this is a linear step (only one main reactant)
            # Exclude atom-mapped reagents and use a more reasonable length threshold
            main_reactants = [r for r in reactants_smiles if not r.startswith("[") and len(r) > 5]
            if len(main_reactants) > 1:
                print(
                    f"Non-linear step found at depth {depth}: {len(main_reactants)} main reactants"
                )
                linear_synthesis = False
                # Add structural constraint for non-linearity
                if {"type": "negation", "details": {"target": "convergent_reaction", "description": "The synthesis must be linear, meaning no reaction step can have more than one main reactant."}} not in findings_json["structural_constraints"]:
                    findings_json["structural_constraints"].append({"type": "negation", "details": {"target": "convergent_reaction", "description": "The synthesis must be linear, meaning no reaction step can have more than one main reactant."}})

            # Check for heterocycle formation at late stage (depth 0, 1, or 2)
            if depth <= 2:
                # Check if any heterocycle exists in product but not in any reactants
                for heterocycle in HETEROCYCLES_OF_INTEREST:
                    product_has_heterocycle = checker.check_ring(heterocycle, product_smiles)

                    if product_has_heterocycle:
                        # Check if any reactants have this heterocycle
                        reactants_with_heterocycle = sum(
                            1 for r in reactants_smiles if checker.check_ring(heterocycle, r)
                        )

                        # If product has the heterocycle and no reactants do, it's formed in this step
                        if reactants_with_heterocycle == 0:
                            print(f"Found {heterocycle} formation at depth {depth}")
                            late_stage_heterocycle = True
                            # Add atomic check for ring system
                            if heterocycle not in findings_json["atomic_checks"]["ring_systems"]:
                                findings_json["atomic_checks"]["ring_systems"].append(heterocycle)
                            # Add structural constraint for late-stage formation
                            if {"type": "positional", "details": {"target": "ring_formation", "position": "late_stage", "description": "A heterocycle from the specified list must be formed in a late stage of the synthesis (depth <= 2)."}} not in findings_json["structural_constraints"]:
                                findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "ring_formation", "position": "late_stage", "description": "A heterocycle from the specified list must be formed in a late stage of the synthesis (depth <= 2)."}})
                            break # Break from heterocycle loop once one is found

        # Determine the next depth based on the current node's type
        next_depth = depth
        if node["type"] != "reaction": # This means it's a 'chemical' node
            next_depth = depth + 1

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, next_depth)

    # Start traversal
    dfs_traverse(route)

    print(f"Linear synthesis: {linear_synthesis}, Late-stage heterocycle: {late_stage_heterocycle}")

    # Return True if all key features are found
    result = linear_synthesis and late_stage_heterocycle
    return result, findings_json

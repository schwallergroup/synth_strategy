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

HETEROCYCLE_TYPES = [
    "furan",
    "pyran",
    "dioxane",
    "tetrahydrofuran",
    "tetrahydropyran",
    "oxirane",
    "oxetane",
    "oxolane",
    "oxane",
    "pyrrole",
    "pyridine",
    "pyrazole",
    "imidazole",
    "oxazole",
    "thiazole",
    "pyrimidine",
    "pyrazine",
    "pyridazine",
    "triazole",
    "tetrazole",
    "pyrrolidine",
    "piperidine",
    "piperazine",
    "morpholine",
    "thiomorpholine",
    "indole",
    "quinoline",
    "isoquinoline",
    "benzoxazole",
    "benzothiazole",
    "benzimidazole",
    "thiophene",
    "isoxazole",
    "isothiazole",
    "oxadiazole",
    "thiadiazole",
]

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects a linear synthetic strategy for constructing heterocycles
    with no convergent steps and a final cyclization.
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
    reaction_count = 0
    has_final_cyclization = False
    max_reactants_per_step = 0
    heterocycle_formed = False
    final_product_smiles = ""

    # Get the final product SMILES from the root node
    if route["type"] == "mol":
        final_product_smiles = route["smiles"]
        print(f"Final product: {final_product_smiles}")

    def dfs_traverse(node, depth=0):
        nonlocal reaction_count, has_final_cyclization, max_reactants_per_step, heterocycle_formed, findings_json

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node.get("metadata", {}).get("rsmi", "")

            if not rsmi:
                print(f"No reaction SMILES found at depth {depth}")
                return

            parts = rsmi.split(">")
            if len(parts) < 3:
                print(f"Invalid reaction SMILES format at depth {depth}: {rsmi}")
                return

            reactants_smiles = [r for r in parts[0].split(".") if r]
            product_smiles = parts[2]

            print(f"Processing reaction at depth {depth}: {rsmi}")

            # Count this reaction
            reaction_count += 1

            # Track maximum number of reactants in any step
            max_reactants_per_step = max(max_reactants_per_step, len(reactants_smiles))
            print(f"Reactants in this step: {len(reactants_smiles)}")

            try:
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles if r]
                product_mol = Chem.MolFromSmiles(product_smiles) if product_smiles else None

                if not product_mol or not all(reactant_mols):
                    print("Invalid molecules detected, skipping reaction")
                    return

                # Check for ring formation in final step (depth 1 in retrosynthetic tree)
                if depth == 1:
                    print(f"Analyzing final synthetic step at depth {depth}")

                    # Count rings in product and reactants
                    product_ring_count = len(Chem.GetSSSR(product_mol))
                    reactant_ring_counts = [len(Chem.GetSSSR(r)) for r in reactant_mols]
                    max_reactant_ring_count = max(reactant_ring_counts, default=0)

                    print(
                        f"Final step - Product rings: {product_ring_count}, Max reactant rings: {max_reactant_ring_count}"
                    )

                    # Check if a new ring was formed
                    if product_ring_count > max_reactant_ring_count:
                        print("New ring formed in final step")
                        # Record ring formation
                        if "ring_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("ring_formation")

                        # Check if the product contains a heterocycle from the predefined list
                        for heterocycle in HETEROCYCLE_TYPES:
                            if checker.check_ring(heterocycle, product_smiles):
                                # Record heterocycle in product
                                if heterocycle not in findings_json["atomic_checks"]["ring_systems"]:
                                    findings_json["atomic_checks"]["ring_systems"].append(heterocycle)

                                # Check if this heterocycle wasn't in any reactant
                                heterocycle_in_reactants = any(
                                    checker.check_ring(heterocycle, r) for r in reactants_smiles
                                )

                                if not heterocycle_in_reactants:
                                    print(f"Heterocycle {heterocycle} formed in final step")
                                    has_final_cyclization = True
                                    heterocycle_formed = True
                                    # Record positional constraint: ring_formation at last_stage
                                    if {"type": "positional", "details": {"target": "ring_formation", "position": "last_stage"}} not in findings_json["structural_constraints"]:
                                        findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "ring_formation", "position": "last_stage"}})
                                    # Record co-occurrence: ring_formation and heterocycle_in_product
                                    if {"type": "co-occurrence", "details": {"targets": ["ring_formation", "heterocycle_in_product"]}} not in findings_json["structural_constraints"]:
                                        findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"targets": ["ring_formation", "heterocycle_in_product"]}})
                                    # Record negation: heterocycle_in_reactants_of_final_step
                                    if {"type": "negation", "details": {"target": "heterocycle_in_reactants_of_final_step"}} not in findings_json["structural_constraints"]:
                                        findings_json["structural_constraints"].append({"type": "negation", "details": {"target": "heterocycle_in_reactants_of_final_step"}})
                                    break

            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Process children with modified depth
        for child in node.get("children", []):
            if node["type"] == "reaction":
                # Depth remains the same when traversing from reaction to chemical
                dfs_traverse(child, depth)
            else:
                # Depth increases when traversing from chemical to reaction
                dfs_traverse(child, depth + 1)

    # Start traversal from the root node
    dfs_traverse(route)

    # Strategy is present if we have multiple steps, final cyclization, and no convergent steps
    # (max 2 reactants per step indicates linear synthesis)
    strategy_present = (
        reaction_count >= 3
        and has_final_cyclization
        and max_reactants_per_step <= 2
        and heterocycle_formed
    )

    # Record structural constraints based on final flags
    if reaction_count >= 3:
        if {"type": "count", "details": {"target": "reaction", "operator": ">=", "value": 3}} not in findings_json["structural_constraints"]:
            findings_json["structural_constraints"].append({"type": "count", "details": {"target": "reaction", "operator": ">=", "value": 3}})
    if max_reactants_per_step <= 2:
        if {"type": "count", "details": {"target": "reactants_per_step", "operator": "<=", "value": 2}} not in findings_json["structural_constraints"]:
            findings_json["structural_constraints"].append({"type": "count", "details": {"target": "reactants_per_step", "operator": "<=", "value": 2}})

    print(f"Linear heterocycle construction strategy detected: {strategy_present}")
    print(f"Reaction count: {reaction_count}, Max reactants per step: {max_reactants_per_step}")
    print(f"Final cyclization: {has_final_cyclization}, Heterocycle formed: {heterocycle_formed}")

    return strategy_present, findings_json

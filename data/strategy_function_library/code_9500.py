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


HETEROCYCLES_OF_INTEREST = [
    "benzofuran",
    "pyrazole",
    "indole",
    "thiophene",
    "furan",
    "pyridine",
    "oxazole",
    "thiazole",
    "imidazole",
    "isoxazole",
    "benzothiazole",
    "benzoxazole",
    "benzimidazole",
    "triazole",
    "tetrazole",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a strategy involving the assembly or elaboration of multiple heterocyclic systems. A route is flagged if it involves at least two instances of heterocycle formation or modification, or if the final product contains at least two distinct heterocycles.
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

    # Track heterocycle formations
    heterocycle_formations = {h: 0 for h in HETEROCYCLES_OF_INTEREST}

    # Track heterocycle-forming reactions
    heterocycle_forming_reactions = []

    # Track heterocycle modifications
    heterocycle_modifications = []

    def dfs_traverse(node, depth=0):
        nonlocal heterocycle_formations, heterocycle_forming_reactions, heterocycle_modifications, findings_json
        # Check reaction nodes for heterocycle formation or modification
        if node["type"] == "reaction":
            try:
                # Extract reactants and product from reaction SMILES
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check for heterocycle formation for each type
                for heterocycle in HETEROCYCLES_OF_INTEREST:
                    # Check if product contains heterocycle
                    product_has_heterocycle = checker.check_ring(heterocycle, product_smiles)
                    if product_has_heterocycle:
                        if heterocycle not in findings_json["atomic_checks"]["ring_systems"]:
                            findings_json["atomic_checks"]["ring_systems"].append(heterocycle)

                    # Check if any reactant has the heterocycle
                    reactants_have_heterocycle = False
                    for r in reactants_smiles:
                        if checker.check_ring(heterocycle, r):
                            reactants_have_heterocycle = True
                            if heterocycle not in findings_json["atomic_checks"]["ring_systems"]:
                                findings_json["atomic_checks"]["ring_systems"].append(heterocycle)
                            break

                    # Check for heterocycle formation
                    if product_has_heterocycle and not reactants_have_heterocycle:
                        heterocycle_formations[heterocycle] += 1
                        print(f"{heterocycle} formation detected in reaction: {rsmi}")
                        if "ring_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("ring_formation")

                        # Check if the reaction is a known heterocycle-forming reaction
                        if checker.check_reaction(heterocycle, rsmi):
                            print(f"Confirmed {heterocycle}-forming reaction: {rsmi}")
                            heterocycle_forming_reactions.append((heterocycle, rsmi))
                            if heterocycle not in findings_json["atomic_checks"]["named_reactions"]:
                                findings_json["atomic_checks"]["named_reactions"].append(heterocycle)

                    # Check for heterocycle modification (heterocycle present in both reactants and product)
                    elif product_has_heterocycle and reactants_have_heterocycle:
                        # This is a modification of an existing heterocycle
                        heterocycle_modifications.append((heterocycle, rsmi))
                        print(f"{heterocycle} modification detected in reaction: {rsmi}")
                        if "ring_modification" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("ring_modification")
                        if heterocycle not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(heterocycle)

            except Exception as e:
                print(f"Error analyzing reaction: {e}")

        # Traverse children
        for child in node.get("children", []):
            # Determine the new depth based on the current node's type
            new_depth = depth
            if node["type"] != "reaction": # Depth increases when going from chemical to reaction
                new_depth = depth + 1
            
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    # Count total heterocycle formations and unique types formed
    total_formations = sum(heterocycle_formations.values())
    unique_types_formed = sum(1 for count in heterocycle_formations.values() if count > 0)

    # Count unique reaction types
    unique_reaction_types = len(set(r[0] for r in heterocycle_forming_reactions))

    print(f"Total heterocycle formations: {total_formations}")
    print(f"Unique heterocycle types formed: {unique_types_formed}")
    print(f"Heterocycle-forming reactions: {len(heterocycle_forming_reactions)}")
    print(f"Unique reaction types: {unique_reaction_types}")
    print(f"Heterocycle modifications: {len(heterocycle_modifications)}")

    has_multiple_types = unique_types_formed >= 2
    has_multiple_formations = total_formations >= 2
    has_specific_reactions = unique_reaction_types >= 2

    # Check if we have at least 2 different heterocycles present in the final product
    final_product_heterocycles = []
    if route["type"] == "mol":
        for heterocycle in HETEROCYCLES_OF_INTEREST:
            if checker.check_ring(heterocycle, route["smiles"]):
                final_product_heterocycles.append(heterocycle)
                if heterocycle not in findings_json["atomic_checks"]["ring_systems"]:
                    findings_json["atomic_checks"]["ring_systems"].append(heterocycle)

    has_multiple_heterocycles_in_product = len(final_product_heterocycles) >= 2
    print(f"Heterocycles in final product: {final_product_heterocycles}")

    # Check if we have at least 2 different heterocycle modifications
    unique_modified_heterocycles = len(set(m[0] for m in heterocycle_modifications))
    has_multiple_heterocycle_modifications = unique_modified_heterocycles >= 2

    result = False

    if has_multiple_types:
        result = True
        findings_json["structural_constraints"].append({"type": "count", "details": {"target": "unique_ring_formation_types", "operator": ">=", "value": 2, "description": "The route involves the formation of at least two different types of heterocycles."}})
    if has_multiple_formations:
        result = True
        findings_json["structural_constraints"].append({"type": "count", "details": {"target": "total_ring_formation_events", "operator": ">=", "value": 2, "description": "The route involves at least two separate heterocycle formation steps (can be of the same type)."}})
    if has_specific_reactions:
        result = True
        findings_json["structural_constraints"].append({"type": "count", "details": {"target": "unique_named_ring_forming_reactions", "operator": ">=", "value": 2, "description": "The route uses at least two different named reactions to form heterocycles."}})
    if has_multiple_heterocycles_in_product:
        result = True
        findings_json["structural_constraints"].append({"type": "count", "details": {"target": "unique_ring_systems_in_final_product", "operator": ">=", "value": 2, "description": "The final product molecule contains at least two different types of heterocycles."}})
    if has_multiple_heterocycle_modifications:
        result = True
        findings_json["structural_constraints"].append({"type": "count", "details": {"target": "unique_ring_modification_types", "operator": ">=", "value": 2, "description": "The route involves the modification of at least two different types of pre-existing heterocycles."}})

    print(f"Has multiple types: {has_multiple_types}")
    print(f"Has multiple formations: {has_multiple_formations}")
    print(f"Has specific reactions: {has_specific_reactions}")
    print(f"Has multiple heterocycles in product: {has_multiple_heterocycles_in_product}")
    print(f"Has multiple heterocycle modifications: {has_multiple_heterocycle_modifications}")
    print(f"Heterocycle assembly strategy detected: {result}")

    return result, findings_json
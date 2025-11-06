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
    This function detects a synthetic strategy involving oxime formation and cleavage.
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

    oxime_formations = {}  # Track oxime formation reactions: {product_smiles: depth}
    oxime_cleavages = {}  # Track oxime cleavage reactions: {reactant_smiles: depth}
    oxime_molecules = {}  # Track all oxime molecules: {smiles: depth}

    def dfs_traverse(node, depth=0):
        nonlocal findings_json
        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check for oxime formation
                # Look for reactions where a carbonyl (aldehyde/ketone) is converted to an oxime
                has_carbonyl_reactant = False
                for r in reactants_smiles:
                    if checker.check_fg("Aldehyde", r):
                        has_carbonyl_reactant = True
                        if "Aldehyde" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Aldehyde")
                    if checker.check_fg("Ketone", r):
                        has_carbonyl_reactant = True
                        if "Ketone" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Ketone")

                has_oxime_product = checker.check_fg("Oxime", product_smiles)
                if has_oxime_product:
                    if "Oxime" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Oxime")

                if has_carbonyl_reactant and has_oxime_product:
                    print(f"Found oxime formation at depth {depth}: {rsmi}")
                    oxime_formations[product_smiles] = depth
                    oxime_molecules[product_smiles] = depth
                    if "oxime_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("oxime_formation")

                # Check for oxime cleavage
                # Look for reactions where an oxime is converted back to a carbonyl
                has_oxime_reactant = False
                for r in reactants_smiles:
                    if checker.check_fg("Oxime", r):
                        has_oxime_reactant = True
                        if "Oxime" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Oxime")

                has_carbonyl_product = False
                if checker.check_fg("Aldehyde", product_smiles):
                    has_carbonyl_product = True
                    if "Aldehyde" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Aldehyde")
                if checker.check_fg("Ketone", product_smiles):
                    has_carbonyl_product = True
                    if "Ketone" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Ketone")

                if has_oxime_reactant and has_carbonyl_product:
                    print(f"Found oxime cleavage at depth {depth}: {rsmi}")
                    for r in reactants_smiles:
                        if checker.check_fg("Oxime", r):
                            oxime_cleavages[r] = depth
                            oxime_molecules[r] = depth
                    if "oxime_cleavage" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("oxime_cleavage")

            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Determine the next depth based on the current node's type
        next_depth = depth
        if node["type"] != "reaction": # If current node is chemical, depth increases
            next_depth = depth + 1

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, next_depth)

    # Start traversal
    dfs_traverse(route)

    # Check if we have both formation and cleavage
    has_oxime_formation = len(oxime_formations) > 0
    has_oxime_cleavage = len(oxime_cleavages) > 0

    print(f"Oxime formations: {oxime_formations}")
    print(f"Oxime cleavages: {oxime_cleavages}")
    print(f"All oxime molecules: {oxime_molecules}")

    # Check if any oxime that was formed was also cleaved
    # This is a more robust check using both canonical SMILES and depth analysis
    common_oximes = False

    # If we have both formations and cleavages, check for common oximes
    if has_oxime_formation and has_oxime_cleavage:
        for formed_oxime, formation_depth in oxime_formations.items():
            formed_mol = Chem.MolFromSmiles(formed_oxime)
            if formed_mol:
                formed_canonical = Chem.MolToSmiles(formed_mol)
                for cleaved_oxime, cleavage_depth in oxime_cleavages.items():
                    cleaved_mol = Chem.MolFromSmiles(cleaved_oxime)
                    if cleaved_mol:
                        cleaved_canonical = Chem.MolToSmiles(cleaved_mol)
                        # Check if they're the same molecule and the formation happens before cleavage
                        if formed_canonical == cleaved_canonical or formed_oxime == cleaved_oxime:
                            if (
                                formation_depth < cleavage_depth
                            ):  # Formation should happen before cleavage
                                common_oximes = True
                                print(
                                    f"Found matching oxime that was formed at depth {formation_depth} and cleaved at depth {cleavage_depth}"
                                )
                                if {"type": "sequence", "details": {"before": "oxime_formation", "after": "oxime_cleavage", "on_same_molecule": True}} not in findings_json["structural_constraints"]:
                                    findings_json["structural_constraints"].append({"type": "sequence", "details": {"before": "oxime_formation", "after": "oxime_cleavage", "on_same_molecule": True}})
                                break

    # If we didn't find direct matches, check if there are oximes in the route
    # that appear to be intermediates (not starting materials or final products)
    if not common_oximes and len(oxime_molecules) > 0:
        # Check if any oxime is an intermediate (not a starting material or final product)
        for oxime_smiles, oxime_depth in oxime_molecules.items():
            # If we have an oxime that's not at depth 0 (not the final product),
            # and we have evidence of either formation or cleavage, consider it a potential strategy
            if oxime_depth > 0 and (has_oxime_formation or has_oxime_cleavage):
                print(f"Found potential oxime intermediate at depth {oxime_depth}")
                common_oximes = True
                if {"type": "positional", "details": {"target": "Oxime", "position": "not_last_stage"}} not in findings_json["structural_constraints"]:
                    findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "Oxime", "position": "not_last_stage"}})
                break

    # For a true oxime intermediate strategy, we need both formation and cleavage,
    # or at least evidence of oximes being used as intermediates
    strategy_detected = (has_oxime_formation and has_oxime_cleavage and common_oximes) or (
        len(oxime_molecules) > 0 and (has_oxime_formation or has_oxime_cleavage)
    )

    print(f"Oxime intermediate strategy detected: {strategy_detected}")
    return strategy_detected, findings_json

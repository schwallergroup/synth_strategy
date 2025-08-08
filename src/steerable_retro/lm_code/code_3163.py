#!/bin/python

"""LM-defined function for strategy description."""

from rdkit.Chem import AllChem, rdFMCS
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

root_data = "/home/dparm/steerable_retro/data"

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


def main(route):
    """
    This function detects if the synthesis follows a linear strategy with
    sequential heteroaromatic functionalization steps.
    """
    # List of heteroaromatic rings to check
    heteroaromatic_rings = [
        "pyridine",
        "pyrimidine",
        "pyrazine",
        "pyridazine",
        "pyrrole",
        "furan",
        "thiophene",
        "imidazole",
        "oxazole",
        "thiazole",
        "triazole",
        "tetrazole",
        "indole",
        "benzimidazole",
        "benzoxazole",
        "benzothiazole",
        "quinoline",
        "isoquinoline",
        "purine",
        "isoxazole",
        "isothiazole",
        "oxadiazole",
        "thiadiazole",
    ]

    # Reactions that typically involve heteroaromatic functionalization
    functionalization_reactions = [
        "Suzuki",
        "Buchwald-Hartwig",
        "N-arylation",
        "Heck",
        "Sonogashira",
        "Negishi",
        "Stille",
        "Friedel-Crafts",
        "Minisci",
        "Directed ortho metalation",
        "Catellani reaction",
        "Minisci-like halide substitution",
        "Aromatic halogenation",
        "Aromatic fluorination",
        "Aromatic chlorination",
        "Aromatic bromination",
        "Aromatic iodination",
        "Aromatic nitration",
        "Aromatic hydroxylation",
    ]

    heteroaromatic_functionalizations = 0
    max_depth = 0
    linear_steps = []  # Track the sequence of functionalization steps

    def dfs_traverse(node, depth=0, parent_product=None):
        nonlocal heteroaromatic_functionalizations, max_depth

        max_depth = max(max_depth, depth)

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                try:
                    # Check if any reactant has heteroaromatic ring
                    reactant_has_heteroaromatic = False
                    reactant_with_heteroaromatic = None

                    for r in reactants:
                        for ring in heteroaromatic_rings:
                            if checker.check_ring(ring, r):
                                reactant_has_heteroaromatic = True
                                reactant_with_heteroaromatic = r
                                print(f"Found {ring} in reactant at depth {depth}")
                                break
                        if reactant_has_heteroaromatic:
                            break

                    # Check if product has heteroaromatic ring
                    product_has_heteroaromatic = False
                    product_ring = None

                    for ring in heteroaromatic_rings:
                        if checker.check_ring(ring, product):
                            product_has_heteroaromatic = True
                            product_ring = ring
                            print(f"Found {ring} in product at depth {depth}")
                            break

                    # Check if this is a known functionalization reaction
                    is_functionalization_reaction = False
                    reaction_type = None
                    for rxn_type in functionalization_reactions:
                        if checker.check_reaction(rxn_type, rsmi):
                            is_functionalization_reaction = True
                            reaction_type = rxn_type
                            print(f"Detected {rxn_type} reaction at depth {depth}")
                            break

                    # Fallback: check for structural changes if no specific reaction detected
                    if (
                        not is_functionalization_reaction
                        and reactant_has_heteroaromatic
                        and product_has_heteroaromatic
                    ):
                        # Check if there's a structural change while preserving the heteroaromatic core
                        if product != reactant_with_heteroaromatic:
                            is_functionalization_reaction = True
                            reaction_type = "Heteroaromatic modification"
                            print(f"Detected generic heteroaromatic modification at depth {depth}")

                    # Verify linearity: check if the heteroaromatic core is preserved
                    is_linear = False
                    if parent_product and reactant_with_heteroaromatic:
                        # Check if parent product contains the same heteroaromatic ring
                        for ring in heteroaromatic_rings:
                            if checker.check_ring(ring, parent_product) and checker.check_ring(
                                ring, reactant_with_heteroaromatic
                            ):
                                is_linear = True
                                print(
                                    f"Linearity confirmed: same {ring} ring in parent product and reactant at depth {depth}"
                                )
                                break
                    else:
                        # If this is the first reaction or we don't have parent info, assume it's linear
                        is_linear = True
                        print(
                            f"Assuming linearity at depth {depth} (first reaction or no parent info)"
                        )

                    # If both reactant and product have heteroaromatic rings and there's a modification
                    if reactant_has_heteroaromatic and product_has_heteroaromatic:
                        if is_functionalization_reaction and is_linear:
                            heteroaromatic_functionalizations += 1
                            linear_steps.append(
                                (depth, rsmi, reaction_type if reaction_type else "Modification")
                            )
                            print(
                                f"Confirmed heteroaromatic functionalization at depth {depth}: {reaction_type if reaction_type else 'Modification'}"
                            )

                except Exception as e:
                    print(f"Error processing reaction at depth {depth}: {e}")
                    pass

                # Pass the current product to child nodes to check linearity
                parent_for_children = product
            else:
                parent_for_children = parent_product
        else:
            # For molecule nodes, pass the molecule SMILES to children
            parent_for_children = node["smiles"] if node["type"] == "mol" else parent_product

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1, parent_for_children)

    dfs_traverse(route)

    # Sort steps by depth to check if they form a sequence
    linear_steps.sort(key=lambda x: x[0])

    # Check if steps form a continuous sequence (depths should be consecutive or close)
    is_sequential = True
    if len(linear_steps) >= 2:
        for i in range(1, len(linear_steps)):
            if linear_steps[i][0] - linear_steps[i - 1][0] > 2:  # Allow for some flexibility
                is_sequential = False
                break

    print(f"Found {heteroaromatic_functionalizations} heteroaromatic functionalizations")
    print(f"Maximum depth: {max_depth}")
    print(f"Linear steps: {linear_steps}")
    print(f"Is sequential: {is_sequential}")

    # Return True if we found at least 3 functionalization steps, the synthesis is deep enough,
    # and the steps form a sequential pattern
    return heteroaromatic_functionalizations >= 3 and max_depth >= 4 and is_sequential

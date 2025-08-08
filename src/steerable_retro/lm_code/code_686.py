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
    Detects a convergent synthesis strategy where two significant fragments
    are coupled in the middle stages of the synthesis.
    """
    # Track if we found fragment coupling
    fragment_coupling_found = False
    # Track the depth at which fragment coupling occurs
    fragment_coupling_depth = None
    # Maximum depth in the route
    max_depth = 0

    # List of common coupling reaction types
    coupling_reactions = [
        "Suzuki coupling with boronic acids",
        "Suzuki coupling with boronic esters",
        "Stille reaction_aryl",
        "Negishi coupling",
        "Heck terminal vinyl",
        "Sonogashira alkyne_aryl halide",
        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
        "Ullmann condensation",
        "Hiyama-Denmark Coupling",
        "Kumada cross-coupling",
    ]

    def is_coupling_reaction(rsmi):
        # Check if the reaction is a known coupling reaction
        for rxn_type in coupling_reactions:
            if checker.check_reaction(rxn_type, rsmi):
                print(f"Detected coupling reaction: {rxn_type}")
                return True

        # If not a known coupling reaction, check for C-C bond formation
        # by analyzing reactants and products (simplified check)
        try:
            reactants_part = rsmi.split(">")[0]
            product_part = rsmi.split(">")[-1]

            # This is a simplified check - in a real implementation,
            # we would need to analyze the actual bond formation
            if product_part and reactants_part:
                return True
        except:
            pass

        return False

    def dfs_traverse(node, depth=0):
        nonlocal fragment_coupling_found, fragment_coupling_depth, max_depth

        max_depth = max(max_depth, depth)

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants_part = rsmi.split(">")[0]
                product_part = rsmi.split(">")[-1]

                # Split reactants
                reactants = reactants_part.split(".")

                # Check if we have multiple significant reactants
                if len(reactants) >= 2:
                    # Convert to molecules
                    reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
                    product_mol = Chem.MolFromSmiles(product_part) if product_part else None

                    if all(reactant_mols) and product_mol:
                        # Check if reactants are significant (at least 5 heavy atoms)
                        significant_reactants = [
                            mol for mol in reactant_mols if mol.GetNumHeavyAtoms() >= 5
                        ]

                        if len(significant_reactants) >= 2 and is_coupling_reaction(rsmi):
                            fragment_coupling_found = True
                            fragment_coupling_depth = depth
                            print(f"Fragment coupling detected at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if this is a mid-stage fragment coupling (in the middle half of synthesis)
    is_mid_stage = (
        fragment_coupling_depth is not None
        and fragment_coupling_depth >= (max_depth / 4)
        and fragment_coupling_depth <= (3 * max_depth / 4)
    )

    result = fragment_coupling_found and is_mid_stage
    print(f"Max depth: {max_depth}, Fragment coupling depth: {fragment_coupling_depth}")
    print(f"Convergent synthesis detection result: {result}")
    return result

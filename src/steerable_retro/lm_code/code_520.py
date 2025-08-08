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
    This function detects if the synthetic route follows a linear strategy to build heterocycles.

    A linear strategy means:
    1. Most reactions have 1-2 significant reactants (not convergent)
    2. The route constructs at least one heterocycle
    3. The route has at least 3 steps
    """
    linear_synthesis = True
    heterocycle_construction = False
    step_count = 0
    non_linear_steps = 0
    max_depth = 0

    # List of common heterocycles to check
    heterocycle_types = [
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
        "indole",
        "benzimidazole",
        "benzoxazole",
        "benzothiazole",
        "quinoline",
        "isoquinoline",
        "furan",
        "thiophene",
        "oxadiazole",
        "thiadiazole",
    ]

    # List of heterocycle formation reactions
    heterocycle_rxn_types = [
        "Paal-Knorr pyrrole synthesis",
        "Formation of NOS Heterocycles",
        "{benzimidazole_derivatives_carboxylic-acid/ester}",
        "{benzimidazole_derivatives_aldehyde}",
        "{benzothiazole}",
        "{benzoxazole_arom-aldehyde}",
        "{benzoxazole_carboxylic-acid}",
        "{thiazole}",
        "{tetrazole_terminal}",
        "{pyrazole}",
        "{Fischer indole}",
        "{indole}",
        "{oxadiazole}",
        "Benzimidazole formation from aldehyde",
        "Benzimidazole formation from acyl halide",
        "Benzimidazole formation from ester/carboxylic acid",
        "Benzoxazole formation from aldehyde",
        "Benzoxazole formation from acyl halide",
        "Benzoxazole formation from ester/carboxylic acid",
        "Benzoxazole formation (intramolecular)",
        "Benzothiazole formation from aldehyde",
        "Benzothiazole formation from acyl halide",
        "Benzothiazole formation from ester/carboxylic acid",
    ]

    def dfs_traverse(node, depth=0):
        nonlocal linear_synthesis, heterocycle_construction, step_count, non_linear_steps, max_depth

        # Track maximum depth for assessing early vs late stage
        max_depth = max(max_depth, depth)

        if node["type"] == "reaction":
            step_count += 1
            print(f"Analyzing reaction at depth {depth}")

            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check if more than two significant reactants (linear vs convergent)
            significant_reactants = []
            for reactant in reactants_smiles:
                if reactant.strip():
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and mol.GetNumHeavyAtoms() > 3:  # Focus on heavy atoms for significance
                        significant_reactants.append(reactant)

            if len(significant_reactants) > 2:
                print(
                    f"Non-linear step detected with {len(significant_reactants)} significant reactants"
                )
                non_linear_steps += 1

            # Check for heterocycle construction
            if product_smiles:
                # First check if this is a known heterocycle formation reaction
                for rxn_type in heterocycle_rxn_types:
                    if checker.check_reaction(rxn_type, rsmi):
                        heterocycle_construction = True
                        print(f"Heterocycle construction reaction detected: {rxn_type}")
                        break

                # If no specific reaction matched, check for heterocycle presence
                if not heterocycle_construction:
                    product_mol = Chem.MolFromSmiles(product_smiles)

                    # Find heterocycles in product
                    product_heterocycles = []
                    for heterocycle in heterocycle_types:
                        if checker.check_ring(heterocycle, product_smiles):
                            product_heterocycles.append(heterocycle)

                    if product_heterocycles:
                        print(f"Product contains heterocycles: {', '.join(product_heterocycles)}")

                        # Check if any reactant doesn't have these heterocycles
                        for heterocycle in product_heterocycles:
                            reactants_with_heterocycle = 0
                            for reactant in reactants_smiles:
                                if reactant.strip() and checker.check_ring(heterocycle, reactant):
                                    reactants_with_heterocycle += 1

                            # If no reactant has this heterocycle or only one does (and we have multiple reactants),
                            # then this is likely a heterocycle construction
                            if reactants_with_heterocycle == 0 or (
                                reactants_with_heterocycle == 1 and len(significant_reactants) > 1
                            ):
                                heterocycle_construction = True
                                print(f"New heterocycle formed: {heterocycle}")
                                break

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Allow some non-linear steps, especially in early stages (higher depth)
    # For a route to be considered non-linear, more than 1/3 of steps should be non-linear
    if non_linear_steps > 0 and non_linear_steps > step_count / 3:
        linear_synthesis = False

    # Return True if both conditions are met and we have multiple steps
    result = linear_synthesis and heterocycle_construction and step_count >= 3
    print(f"Linear heterocycle construction strategy detected: {result}")
    print(
        f"- Linear synthesis: {linear_synthesis} (non-linear steps: {non_linear_steps}/{step_count})"
    )
    print(f"- Heterocycle construction: {heterocycle_construction}")
    print(f"- Step count: {step_count}")
    print(f"- Maximum depth: {max_depth}")

    return result

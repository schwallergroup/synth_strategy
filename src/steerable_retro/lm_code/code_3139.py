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
    This function detects if the synthetic route involves the construction
    of complex heterocyclic systems.
    """
    heterocycle_construction = False

    # List of common heterocyclic rings to check
    heterocycle_rings = [
        "furan",
        "pyrrole",
        "thiophene",
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
        "purine",
        "oxadiazole",
        "thiadiazole",
        "isoxazole",
    ]

    # List of heterocycle formation reactions
    heterocycle_reactions = [
        "Paal-Knorr pyrrole synthesis",
        "Fischer indole",
        "Pictet-Spengler",
        "benzimidazole_derivatives_carboxylic-acid/ester",
        "benzimidazole_derivatives_aldehyde",
        "benzothiazole",
        "benzoxazole_arom-aldehyde",
        "benzoxazole_carboxylic-acid",
        "thiazole",
        "tetrazole_terminal",
        "tetrazole_connect_regioisomere_1",
        "tetrazole_connect_regioisomere_2",
        "1,2,4-triazole_acetohydrazide",
        "1,2,4-triazole_carboxylic-acid/ester",
        "pyrazole",
        "oxadiazole",
        "benzofuran",
        "benzothiophene",
        "indole",
        "imidazole",
    ]

    def dfs_traverse(node, depth=0):
        nonlocal heterocycle_construction

        if heterocycle_construction:
            return  # Early return if we already found heterocycle construction

        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check if this is a known heterocycle formation reaction
                if any(
                    checker.check_reaction(reaction, rsmi) for reaction in heterocycle_reactions
                ):
                    print(f"Found heterocycle formation reaction at depth {depth}: {rsmi}")
                    heterocycle_construction = True
                    return

                # Check if product contains heterocycles that weren't in reactants
                product_mol = Chem.MolFromSmiles(product_smiles)
                if product_mol:
                    # Check for heterocycles in product
                    product_heterocycles = []
                    for ring in heterocycle_rings:
                        if checker.check_ring(ring, product_smiles):
                            product_heterocycles.append(ring)

                    if product_heterocycles:
                        # Check if reactants don't have the same heterocycles
                        reactant_heterocycles = set()
                        for reactant in reactants_smiles:
                            for ring in product_heterocycles:
                                if checker.check_ring(ring, reactant):
                                    reactant_heterocycles.add(ring)

                        # If any heterocycle in product is not in reactants, it was constructed
                        new_heterocycles = [
                            ring
                            for ring in product_heterocycles
                            if ring not in reactant_heterocycles
                        ]
                        if new_heterocycles:
                            print(
                                f"Found heterocycle construction at depth {depth}: {new_heterocycles}"
                            )
                            heterocycle_construction = True
                            return
            except Exception as e:
                print(f"Error analyzing heterocycle construction: {e}")

        # Continue traversing children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    return heterocycle_construction

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
    This function detects if the synthesis route involves oxazole formation
    from acyclic precursors.
    """
    oxazole_found = False

    def dfs_traverse(node, depth=0):
        nonlocal oxazole_found

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            print(f"Depth {depth} - Analyzing reaction: {rsmi}")

            # Check if this is a known oxazole-forming reaction
            oxazole_forming_reactions = [
                "benzoxazole formation from aldehyde",
                "benzoxazole formation from acyl halide",
                "benzoxazole formation from ester/carboxylic acid",
                "benzoxazole formation (intramolecular)",
                "{benzoxazole_arom-aldehyde}",
                "{benzoxazole_carboxylic-acid}",
            ]

            for rxn_type in oxazole_forming_reactions:
                if checker.check_reaction(rxn_type, rsmi):
                    print(f"Detected {rxn_type} reaction")
                    oxazole_found = True
                    return

            # Check if oxazole is formed in this reaction
            product_mol = Chem.MolFromSmiles(product_smiles)

            if product_mol and checker.check_ring("oxazole", product_smiles):
                print(f"Oxazole ring found in product: {product_smiles}")

                # Check if oxazole was not present in reactants
                oxazole_in_reactants = False
                for reactant in reactants_smiles:
                    if checker.check_ring("oxazole", reactant):
                        print(f"Oxazole ring already present in reactant: {reactant}")
                        oxazole_in_reactants = True
                        break

                if not oxazole_in_reactants:
                    print(
                        "Oxazole formation confirmed - ring present in product but not in reactants"
                    )
                    oxazole_found = True

                    # Also check for benzoxazole which is a specific type of oxazole
                    if checker.check_ring("benzoxazole", product_smiles):
                        print("Product contains benzoxazole ring")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    print("Starting traversal of synthesis route")
    dfs_traverse(route)
    print(f"Oxazole formation strategy result: {oxazole_found}")
    return oxazole_found

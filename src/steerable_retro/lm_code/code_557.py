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
    This function detects a synthetic strategy that uses hydrazine derivatives
    to form nitrogen heterocycles.
    """
    hydrazine_reactions = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            # Extract reactants and product from reaction
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check for hydrazine derivatives in reactants
                has_hydrazine_derivative = any(
                    checker.check_fg("Hydrazine", smi)
                    or checker.check_fg("Hydrazone", smi)
                    or checker.check_fg("Acylhydrazine", smi)
                    or checker.check_fg("Hydrazone amide", smi)
                    for smi in reactants_smiles
                    if smi
                )

                # If hydrazine derivative found, check for heterocycle formation
                if has_hydrazine_derivative and product_smiles:
                    # Count rings in reactants and product
                    reactants_mols = [Chem.MolFromSmiles(smi) for smi in reactants_smiles if smi]
                    product_mol = Chem.MolFromSmiles(product_smiles)

                    reactants_ring_count = sum(
                        len(mol.GetRingInfo().AtomRings()) if mol else 0 for mol in reactants_mols
                    )

                    product_ring_count = (
                        len(product_mol.GetRingInfo().AtomRings()) if product_mol else 0
                    )

                    # Check if product has more rings than reactants combined
                    if product_ring_count > reactants_ring_count:
                        # Verify that the product contains a nitrogen heterocycle
                        n_heterocycles = [
                            ring_name
                            for ring_name in [
                                "pyrazole",
                                "triazole",
                                "tetrazole",
                                "pyridazine",
                                "imidazole",
                                "indazole",
                                "benzotriazole",
                                "piperazine",
                                "pyrimidine",
                                "pyrazine",
                                "diazepane",
                                "aziridine",
                                "azetidine",
                                "pyrrolidine",
                                "piperidine",
                                "azepane",
                            ]
                            if checker.check_ring(ring_name, product_smiles)
                        ]

                        if n_heterocycles:
                            print(
                                f"Hydrazine-based heterocycle formation detected at depth {depth}"
                            )
                            print(f"Formed heterocycles: {', '.join(n_heterocycles)}")
                            hydrazine_reactions.append(depth)
            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if we have at least one hydrazine-based heterocycle formation
    if len(hydrazine_reactions) >= 1:
        print("Hydrazine-based heterocycle strategy detected")
        return True

    return False

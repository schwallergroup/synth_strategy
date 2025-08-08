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
    This function detects a strategy where a heterocycle is formed through
    cyclization involving a nitrile group.
    """
    # Track if we found the key feature
    found_nitrile_cyclization = False

    def dfs_traverse(node):
        nonlocal found_nitrile_cyclization

        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check for nitrile in reactants
                nitrile_in_reactants = False
                nitrile_reactant_idx = -1

                for i, reactant in enumerate(reactants_smiles):
                    if checker.check_fg("Nitrile", reactant):
                        nitrile_in_reactants = True
                        nitrile_reactant_idx = i
                        print(f"Found nitrile in reactant: {reactant}")
                        break

                if nitrile_in_reactants:
                    # Convert to RDKit molecules
                    product_mol = Chem.MolFromSmiles(product_smiles)
                    reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles if r]

                    # Count rings in reactants and product
                    reactant_rings = sum(
                        [mol.GetRingInfo().NumRings() for mol in reactant_mols if mol]
                    )
                    product_rings = product_mol.GetRingInfo().NumRings() if product_mol else 0

                    # Check if a new ring was formed
                    if product_rings > reactant_rings:
                        # Check if the product contains a heterocycle
                        has_heterocycle = False

                        # Check for common heterocycles
                        heterocycles = [
                            "pyrrole",
                            "pyridine",
                            "pyrazole",
                            "imidazole",
                            "oxazole",
                            "thiazole",
                            "pyrimidine",
                            "pyrazine",
                            "triazole",
                            "tetrazole",
                            "furan",
                            "thiophene",
                            "oxadiazole",
                            "thiadiazole",
                            "isoxazole",
                            "isothiazole",
                            "benzoxazole",
                            "benzothiazole",
                            "benzimidazole",
                        ]

                        for ring in heterocycles:
                            if checker.check_ring(ring, product_smiles):
                                has_heterocycle = True
                                print(f"Found heterocycle {ring} in product: {product_smiles}")
                                break

                        # Check if nitrile is consumed (not present in product)
                        nitrile_consumed = not checker.check_fg("Nitrile", product_smiles)

                        if has_heterocycle and nitrile_consumed:
                            found_nitrile_cyclization = True
                            print(
                                f"Confirmed heterocycle formation via nitrile cyclization: {rsmi}"
                            )

            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return found_nitrile_cyclization

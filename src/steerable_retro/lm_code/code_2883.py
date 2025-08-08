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
    This function detects a linear synthetic strategy with sequential heterocycle formations.
    """
    # Track heterocycle formations and their sequence
    heterocycle_formations = []
    linear_synthesis = True

    def is_heterocycle(smiles):
        """Check if a molecule contains a heterocycle"""
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return False

        # Check for common heterocycles
        heterocycle_types = [
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
            "thiophene",
            "thiopyran",
        ]

        for ring_type in heterocycle_types:
            if checker.check_ring(ring_type, smiles):
                return True

        return False

    def count_heterocycles(smiles):
        """Count the number of heterocycles in a molecule"""
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return 0

        count = 0
        # Check for common heterocycles
        heterocycle_types = [
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
            "thiophene",
            "thiopyran",
        ]

        for ring_type in heterocycle_types:
            if checker.check_ring(ring_type, smiles):
                count += 1

        return count

    def is_complex_molecule(smiles):
        """Determine if a molecule is complex based on multiple criteria"""
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return False

        # Consider a molecule complex if it has:
        # 1. More than 10 heavy atoms
        # 2. At least one ring
        # 3. At least one functional group
        heavy_atom_count = mol.GetNumHeavyAtoms()
        ring_count = len(Chem.GetSSSR(mol))

        # Check for some common functional groups
        has_functional_group = False
        functional_groups = [
            "Carboxylic acid",
            "Ester",
            "Amide",
            "Amine",
            "Alcohol",
            "Ketone",
            "Aldehyde",
        ]
        for fg in functional_groups:
            if checker.check_fg(fg, smiles):
                has_functional_group = True
                break

        return heavy_atom_count > 10 and (ring_count > 0 or has_functional_group)

    def dfs_traverse(node, depth=0):
        nonlocal heterocycle_formations, linear_synthesis

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Count heterocycles in reactants and product
                product_heterocycles = count_heterocycles(product)
                max_reactant_heterocycles = 0
                for reactant in reactants:
                    reactant_heterocycles = count_heterocycles(reactant)
                    max_reactant_heterocycles = max(
                        max_reactant_heterocycles, reactant_heterocycles
                    )

                # If product has more heterocycles than any reactant, a heterocycle was formed
                if product_heterocycles > max_reactant_heterocycles:
                    print(
                        f"Heterocycle formation detected at depth {depth}. Product heterocycles: {product_heterocycles}, Max reactant heterocycles: {max_reactant_heterocycles}"
                    )
                    heterocycle_formations.append(depth)

                # Check if this is a convergent step (more than one complex reactant)
                complex_reactants = 0
                for reactant in reactants:
                    if is_complex_molecule(reactant):
                        complex_reactants += 1

                if complex_reactants > 1:
                    linear_synthesis = False
                    print(f"Convergent synthesis step detected at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if heterocycle formations are sequential
    sequential_formations = False
    if len(heterocycle_formations) >= 2:
        # Sort by depth to check sequence
        heterocycle_formations.sort()

        # Check if there are at least two formations with consecutive depths
        # or with only reasonable gaps between them
        for i in range(len(heterocycle_formations) - 1):
            if (
                heterocycle_formations[i + 1] - heterocycle_formations[i] <= 5
            ):  # Allow larger gaps in linear synthesis
                sequential_formations = True
                break

    print(f"Heterocycle formations: {heterocycle_formations}")
    print(f"Linear synthesis: {linear_synthesis}")
    print(f"Sequential formations: {sequential_formations}")

    # Return True if we have multiple sequential heterocycle formations in a linear synthesis
    # Also consider the overall distribution of formations
    if len(heterocycle_formations) >= 2 and linear_synthesis:
        # Either they're sequential by our definition, or they're reasonably distributed
        if sequential_formations or (
            max(heterocycle_formations) - min(heterocycle_formations)
            <= len(heterocycle_formations) * 3
        ):
            return True

    return False

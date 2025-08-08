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


def main(route, amine_product=None):
    """
    Checks if the route contains a late-stage heterocyclization reaction.
    If amine_product is provided, checks if the heterocyclization involves the amine.
    """
    # Track if we've found a heterocyclization at a late stage (low depth)
    found_late_cyclization = False

    def dfs(node, depth=0):
        nonlocal found_late_cyclization

        # For reaction nodes at low depth (late stage), check if it's a heterocyclization
        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
            and depth <= 2
        ):
            rxn_smiles = node["metadata"]["rsmi"]

            # Check for various heterocycle formation reactions
            is_heterocyclization = (
                checker.check_reaction("Formation of NOS Heterocycles", rxn_smiles)
                or checker.check_reaction("Paal-Knorr pyrrole synthesis", rxn_smiles)
                or checker.check_reaction("Benzothiazole formation from aldehyde", rxn_smiles)
                or checker.check_reaction("Benzothiazole formation from acyl halide", rxn_smiles)
                or checker.check_reaction(
                    "Benzothiazole formation from ester/carboxylic acid", rxn_smiles
                )
                or checker.check_reaction("Benzoxazole formation from aldehyde", rxn_smiles)
                or checker.check_reaction("Benzoxazole formation from acyl halide", rxn_smiles)
                or checker.check_reaction(
                    "Benzoxazole formation from ester/carboxylic acid", rxn_smiles
                )
                or checker.check_reaction("Benzoxazole formation (intramolecular)", rxn_smiles)
                or checker.check_reaction("Benzimidazole formation from aldehyde", rxn_smiles)
                or checker.check_reaction("Benzimidazole formation from acyl halide", rxn_smiles)
                or checker.check_reaction(
                    "Benzimidazole formation from ester/carboxylic acid", rxn_smiles
                )
                or checker.check_reaction(
                    "Intramolecular amination of azidobiphenyls (heterocycle formation)", rxn_smiles
                )
                or checker.check_reaction(
                    "Intramolecular amination (heterocycle formation)", rxn_smiles
                )
            )

            if is_heterocyclization:
                print(f"Found late-stage heterocyclization at depth {depth}: {rxn_smiles}")

                # Verify that the product contains a heterocycle
                try:
                    reactants = rxn_smiles.split(">")[0].split(".")
                    product = rxn_smiles.split(">")[-1]

                    # Check for common heterocycles in the product
                    has_heterocycle = (
                        checker.check_ring("benzimidazole", product)
                        or checker.check_ring("benzoxazole", product)
                        or checker.check_ring("benzothiazole", product)
                        or checker.check_ring("indole", product)
                        or checker.check_ring("pyrrole", product)
                        or checker.check_ring("imidazole", product)
                        or checker.check_ring("oxazole", product)
                        or checker.check_ring("thiazole", product)
                        or checker.check_ring("triazole", product)
                        or checker.check_ring("tetrazole", product)
                    )

                    if has_heterocycle:
                        print(f"Confirmed heterocycle formation: heterocycle found in product")

                        # If we have an amine product, check if it's involved in the heterocyclization
                        if amine_product:
                            amine_involved = False
                            for r in reactants:
                                try:
                                    amine_mol = Chem.MolFromSmiles(amine_product)
                                    r_mol = Chem.MolFromSmiles(r)
                                    if (
                                        amine_mol
                                        and r_mol
                                        and (
                                            r_mol.HasSubstructMatch(amine_mol)
                                            or amine_mol.HasSubstructMatch(r_mol)
                                        )
                                    ):
                                        amine_involved = True
                                        break
                                except Exception as e:
                                    print(f"Error checking amine involvement: {e}")

                            if amine_involved:
                                print(f"Heterocyclization involves the amine from nitro reduction")
                                found_late_cyclization = True
                            else:
                                print(
                                    f"Heterocyclization does NOT involve the amine from nitro reduction"
                                )
                        else:
                            # If we're not checking for amine involvement, just mark as found
                            found_late_cyclization = True
                    else:
                        print(f"No heterocycle found in product")
                except Exception as e:
                    print(f"Error analyzing heterocyclization: {e}")

        # Continue DFS traversal
        for child in node.get("children", []):
            dfs(child, depth + 1)

    # Start DFS traversal from the root
    dfs(route)
    return found_late_cyclization

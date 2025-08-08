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


def main(route, diarylacetylene_product=None):
    """
    Checks if the route contains a nitro reduction to form an amine.
    If diarylacetylene_product is provided, checks if the nitro reduction
    is performed on a molecule derived from the diarylacetylene scaffold.
    """
    # Track if we've found a nitro reduction
    found_nitro_reduction = False
    amine_product = None

    def dfs(node, depth=0, parent_smiles=None):
        nonlocal found_nitro_reduction, amine_product

        current_smiles = None

        # For molecule nodes, store the SMILES
        if node["type"] == "mol":
            current_smiles = node["smiles"]

            # If we have a diarylacetylene product, check if this molecule is derived from it
            if diarylacetylene_product and parent_smiles:
                # Simple substructure check - not perfect but a reasonable approximation
                try:
                    parent_mol = Chem.MolFromSmiles(parent_smiles)
                    current_mol = Chem.MolFromSmiles(current_smiles)
                    if parent_mol and current_mol and current_mol.HasSubstructMatch(parent_mol):
                        print(f"Molecule at depth {depth} is derived from diarylacetylene scaffold")
                except Exception as e:
                    print(f"Error checking molecule derivation: {e}")

        # For reaction nodes, check if it's a nitro reduction
        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rxn_smiles = node["metadata"]["rsmi"]

            # Check if this is a nitro reduction
            if checker.check_reaction("Reduction of nitro groups to amines", rxn_smiles):
                print(f"Found nitro reduction at depth {depth}: {rxn_smiles}")

                # Verify that the product has an amine group and a reactant had a nitro group
                try:
                    reactants = rxn_smiles.split(">")[0].split(".")
                    product = rxn_smiles.split(">")[-1]

                    has_nitro_reactant = any(checker.check_fg("Nitro group", r) for r in reactants)
                    has_amine_product = (
                        checker.check_fg("Primary amine", product)
                        or checker.check_fg("Secondary amine", product)
                        or checker.check_fg("Tertiary amine", product)
                        or checker.check_fg("Aniline", product)
                    )

                    if has_nitro_reactant and has_amine_product:
                        print(
                            f"Confirmed nitro reduction: nitro group in reactant and amine in product"
                        )

                        # If we have a diarylacetylene product, check if the nitro reduction is related
                        if diarylacetylene_product:
                            # Check if any reactant contains the diarylacetylene scaffold
                            diaryl_related = False
                            for r in reactants:
                                try:
                                    diaryl_mol = Chem.MolFromSmiles(diarylacetylene_product)
                                    r_mol = Chem.MolFromSmiles(r)
                                    if (
                                        diaryl_mol
                                        and r_mol
                                        and (
                                            r_mol.HasSubstructMatch(diaryl_mol)
                                            or diaryl_mol.HasSubstructMatch(r_mol)
                                        )
                                    ):
                                        diaryl_related = True
                                        break
                                except Exception as e:
                                    print(f"Error checking diarylacetylene relation: {e}")

                            if diaryl_related:
                                print(f"Nitro reduction is related to diarylacetylene scaffold")
                                found_nitro_reduction = True
                                amine_product = product
                            else:
                                print(f"Nitro reduction is NOT related to diarylacetylene scaffold")
                        else:
                            # If we're not checking for relation to diarylacetylene, just mark as found
                            found_nitro_reduction = True
                            amine_product = product
                    else:
                        print(
                            f"Not a valid nitro reduction: missing nitro in reactant or amine in product"
                        )
                except Exception as e:
                    print(f"Error analyzing nitro reduction: {e}")

        # Continue DFS traversal
        for child in node.get("children", []):
            dfs(child, depth + 1, current_smiles)

    # Start DFS traversal from the root
    dfs(route)
    return found_nitro_reduction, amine_product

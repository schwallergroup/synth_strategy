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
    This function detects if the synthesis route involves aromatic cyanation
    (replacement of aryl halide with nitrile).
    """
    has_aromatic_cyanation = False

    def dfs_traverse(node, depth=0):
        nonlocal has_aromatic_cyanation

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if any reactant has an aromatic halide
                has_aryl_halide = False
                has_cyanide_source = False

                for reactant in reactants:
                    if checker.check_fg("Aromatic halide", reactant):
                        has_aryl_halide = True
                        print(f"Found aromatic halide in reactant: {reactant}")

                    # Check for potential cyanide sources - look for CN, metal cyanides, etc.
                    if (
                        "CN" in reactant
                        or "C#N" in reactant
                        or "N#C" in reactant
                        or "[N-]" in reactant
                    ):
                        has_cyanide_source = True
                        print(f"Found potential cyanide source: {reactant}")

                # Check for nitrile in product
                has_nitrile_product = checker.check_fg("Nitrile", product)
                if has_nitrile_product:
                    print(f"Found nitrile in product: {product}")

                    # Verify this is a cyanation reaction
                    if has_aryl_halide and has_cyanide_source:
                        print("Found both aromatic halide and cyanide source")

                        # Try to check if this is a known reaction type
                        if checker.check_reaction("Aromatic dehalogenation", rsmi):
                            print(f"Detected aromatic dehalogenation at depth {depth}")
                            has_aromatic_cyanation = True
                        else:
                            # Perform additional checks to confirm cyanation
                            product_mol = Chem.MolFromSmiles(product)

                            # Check if the product has an aromatic ring with a nitrile attached
                            if product_mol:
                                aryl_nitrile_pattern = Chem.MolFromSmarts("c-C#N")
                                if product_mol.HasSubstructMatch(aryl_nitrile_pattern):
                                    print(f"Found aryl-nitrile pattern in product")

                                    # If we have an aromatic ring with nitrile and had aryl halide in reactants
                                    # plus a cyanide source, it's likely an aromatic cyanation
                                    has_aromatic_cyanation = True
                                    print(
                                        f"Detected aromatic cyanation based on structural analysis at depth {depth}"
                                    )

                    # Even if we don't have explicit cyanide source in reactants, check the transformation
                    elif has_aryl_halide and has_nitrile_product:
                        # Check if the product has an aromatic ring with a nitrile attached
                        product_mol = Chem.MolFromSmiles(product)
                        if product_mol:
                            aryl_nitrile_pattern = Chem.MolFromSmarts("c-C#N")
                            if product_mol.HasSubstructMatch(aryl_nitrile_pattern):
                                print(
                                    f"Found aryl-nitrile pattern in product with aromatic halide in reactants"
                                )
                                # This is likely a cyanation even if we don't see the cyanide source
                                has_aromatic_cyanation = True
                                print(f"Detected probable aromatic cyanation at depth {depth}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return has_aromatic_cyanation

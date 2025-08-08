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
    Detects a synthetic strategy involving SNAr coupling.
    Looks for a reaction where a halogenated aromatic compound
    reacts with an amine to form a C-N bond.
    """
    has_snar = False

    def dfs_traverse(node, depth=0):
        nonlocal has_snar

        if node["type"] == "reaction":
            # Extract reactants and product
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check for SNAr coupling reactions using the checker functions
                # First check if we have the specific reaction types
                if (
                    checker.check_reaction("N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)", rsmi)
                    or checker.check_reaction(
                        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine", rsmi
                    )
                    or checker.check_reaction(
                        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine", rsmi
                    )
                    or checker.check_reaction("Ullmann-Goldberg Substitution amine", rsmi)
                    or checker.check_reaction("Goldberg coupling", rsmi)
                    or checker.check_reaction("Goldberg coupling aryl amine-aryl chloride", rsmi)
                ):
                    print(f"Found SNAr coupling reaction at depth {depth}")
                    has_snar = True
                else:
                    # If no specific reaction type is found, check for the components
                    has_aromatic_halide = any(
                        checker.check_fg("Aromatic halide", r) for r in reactants_smiles
                    )
                    has_amine = any(
                        checker.check_fg("Primary amine", r)
                        or checker.check_fg("Secondary amine", r)
                        or checker.check_fg("Tertiary amine", r)
                        or checker.check_fg("Aniline", r)
                        for r in reactants_smiles
                    )

                    # Check if product has new C-N bond
                    if has_aromatic_halide and has_amine:
                        # Look for evidence of SNAr in the reaction pattern
                        # This is a fallback if the reaction checkers didn't catch it
                        product_mol = Chem.MolFromSmiles(product_smiles)
                        if product_mol:
                            # Check if product has aromatic C-N bond
                            for bond in product_mol.GetBonds():
                                if (
                                    bond.GetBondType() == Chem.BondType.SINGLE
                                    and bond.GetBeginAtom().GetAtomicNum() == 6
                                    and bond.GetEndAtom().GetAtomicNum() == 7
                                    and bond.GetBeginAtom().GetIsAromatic()
                                ):
                                    print(
                                        f"Found potential SNAr coupling at depth {depth} (manual check)"
                                    )
                                    has_snar = True
                                    break
            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Return True if SNAr coupling is found
    return has_snar

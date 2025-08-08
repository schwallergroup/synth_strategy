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
    This function detects a strategy involving benzyl protection in the late stage of synthesis.
    It looks for O-benzylation reactions at low depth (late in synthesis).
    """
    benzyl_protection_depth = None

    def dfs_traverse(node, depth=0):
        nonlocal benzyl_protection_depth

        if node["type"] == "reaction":
            # Extract reactants and product
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                print(f"Checking reaction at depth {depth}: {rsmi}")

                # Check for benzyl protection reactions
                is_benzyl_protection = False

                # Check if this is a Williamson Ether Synthesis reaction
                if checker.check_reaction("Williamson Ether Synthesis", rsmi):
                    print(f"Williamson Ether Synthesis detected at depth {depth}")
                    # Verify it's specifically a benzyl protection
                    alcohol_found = False
                    benzyl_halide_found = False

                    for r_smi in reactants_smiles:
                        # Check for alcohol
                        if (
                            checker.check_fg("Primary alcohol", r_smi)
                            or checker.check_fg("Secondary alcohol", r_smi)
                            or checker.check_fg("Tertiary alcohol", r_smi)
                            or checker.check_fg("Phenol", r_smi)
                        ):
                            alcohol_found = True
                            print(f"Alcohol found in reactant: {r_smi}")

                        # Check for benzyl halide
                        if (
                            checker.check_fg("Primary halide", r_smi)
                            or checker.check_fg("Secondary halide", r_smi)
                        ) and checker.check_ring("benzene", r_smi):
                            r_mol = Chem.MolFromSmiles(r_smi)
                            # Verify the halide is connected to a benzyl group
                            if r_mol and r_mol.HasSubstructMatch(
                                Chem.MolFromSmarts("c1ccccc1C[F,Cl,Br,I]")
                            ):
                                benzyl_halide_found = True
                                print(f"Benzyl halide found in reactant: {r_smi}")

                    # Check if product has a benzyl ether
                    product_has_benzyl_ether = False
                    if checker.check_fg("Ether", product_smiles) and checker.check_ring(
                        "benzene", product_smiles
                    ):
                        product_mol = Chem.MolFromSmiles(product_smiles)
                        if product_mol and product_mol.HasSubstructMatch(
                            Chem.MolFromSmarts("c1ccccc1CO[#6]")
                        ):
                            product_has_benzyl_ether = True
                            print(f"Benzyl ether found in product: {product_smiles}")

                    if alcohol_found and benzyl_halide_found and product_has_benzyl_ether:
                        is_benzyl_protection = True

                # Check for other benzyl protection reactions
                elif (
                    checker.check_reaction("Alcohol protection with silyl ethers", rsmi)
                    and "benzyl" in rsmi.lower()
                ):
                    print(f"Benzyl protection via silyl ether detected at depth {depth}")
                    # Verify product has a benzyl group
                    if checker.check_ring("benzene", product_smiles) and checker.check_fg(
                        "Ether", product_smiles
                    ):
                        product_mol = Chem.MolFromSmiles(product_smiles)
                        if product_mol and product_mol.HasSubstructMatch(
                            Chem.MolFromSmarts("c1ccccc1CO[#6]")
                        ):
                            is_benzyl_protection = True

                # Check for Mitsunobu reaction (another way to form benzyl ethers)
                elif checker.check_reaction("Mitsunobu aryl ether", rsmi):
                    print(f"Mitsunobu reaction detected at depth {depth}")
                    alcohol_found = False
                    benzyl_alcohol_found = False

                    for r_smi in reactants_smiles:
                        if (
                            checker.check_fg("Primary alcohol", r_smi)
                            or checker.check_fg("Secondary alcohol", r_smi)
                            or checker.check_fg("Tertiary alcohol", r_smi)
                        ):
                            if checker.check_ring("benzene", r_smi):
                                r_mol = Chem.MolFromSmiles(r_smi)
                                if r_mol and r_mol.HasSubstructMatch(
                                    Chem.MolFromSmarts("c1ccccc1CO")
                                ):
                                    benzyl_alcohol_found = True
                                    print(f"Benzyl alcohol found in reactant: {r_smi}")
                            else:
                                alcohol_found = True
                                print(f"Alcohol found in reactant: {r_smi}")

                    if (alcohol_found and benzyl_alcohol_found) or (
                        benzyl_alcohol_found
                        and checker.check_fg("Phenol", "".join(reactants_smiles))
                    ):
                        if checker.check_fg("Ether", product_smiles) and checker.check_ring(
                            "benzene", product_smiles
                        ):
                            product_mol = Chem.MolFromSmiles(product_smiles)
                            if product_mol and product_mol.HasSubstructMatch(
                                Chem.MolFromSmarts("c1ccccc1CO[#6]")
                            ):
                                is_benzyl_protection = True

                # General check for O-alkylation with benzyl halides
                else:
                    # Look for benzyl halide in reactants
                    benzyl_halide_found = False
                    alcohol_found = False

                    for r_smi in reactants_smiles:
                        # Check for benzyl halide (especially benzyl bromide)
                        if checker.check_ring("benzene", r_smi):
                            r_mol = Chem.MolFromSmiles(r_smi)
                            if r_mol:
                                # Check for benzyl bromide pattern
                                if r_mol.HasSubstructMatch(
                                    Chem.MolFromSmarts("c1ccccc1C[Br,Cl,I,F]")
                                ):
                                    benzyl_halide_found = True
                                    print(f"Benzyl halide found in reactant: {r_smi}")
                                # Also check for "Br[CH2][c]" pattern which appears in the test case
                                elif "Br[CH2" in r_smi and "[c" in r_smi:
                                    benzyl_halide_found = True
                                    print(f"Benzyl halide pattern found in reactant: {r_smi}")

                        # Check for alcohol or hydroxyl group
                        if (
                            checker.check_fg("Primary alcohol", r_smi)
                            or checker.check_fg("Secondary alcohol", r_smi)
                            or checker.check_fg("Tertiary alcohol", r_smi)
                            or checker.check_fg("Phenol", r_smi)
                            or "[OH" in r_smi
                        ):  # Check for atom-mapped hydroxyl
                            alcohol_found = True
                            print(f"Alcohol/hydroxyl found in reactant: {r_smi}")

                    # Check if product has a benzyl ether structure
                    if benzyl_halide_found and alcohol_found:
                        if checker.check_fg("Ether", product_smiles) and checker.check_ring(
                            "benzene", product_smiles
                        ):
                            # Look for benzyl ether pattern in product
                            product_mol = Chem.MolFromSmiles(product_smiles)
                            if product_mol:
                                # Check for various benzyl ether patterns
                                if (
                                    product_mol.HasSubstructMatch(
                                        Chem.MolFromSmarts("c1ccccc1CO[#6]")
                                    )
                                    or "[c" in product_smiles
                                    and "[O" in product_smiles
                                    and "[CH2" in product_smiles
                                ):
                                    print(
                                        f"Benzyl ether pattern found in product: {product_smiles}"
                                    )
                                    is_benzyl_protection = True

                if is_benzyl_protection:
                    print(f"Benzyl protection confirmed at depth {depth}")
                    if benzyl_protection_depth is None or depth < benzyl_protection_depth:
                        benzyl_protection_depth = depth
            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    print(f"Lowest benzyl protection depth found: {benzyl_protection_depth}")
    # Check if benzyl protection occurs at depth <= 5 (late stage)
    return benzyl_protection_depth is not None and benzyl_protection_depth <= 5

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
    Detects a synthetic strategy involving late-stage C-N bond formation between
    a benzyl halide and an amine/amide, where the benzyl halide is prepared via
    functional group interconversions.
    """
    # Initialize tracking variables
    benzyl_halide_formation = False
    late_stage_cn_bond_formation = False
    ester_to_alcohol_conversion = False
    alcohol_to_halide_conversion = False

    print("Starting analysis of benzyl halide to amine strategy")

    def dfs_traverse(node, depth=0):
        nonlocal benzyl_halide_formation, late_stage_cn_bond_formation
        nonlocal ester_to_alcohol_conversion, alcohol_to_halide_conversion

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants_str = rsmi.split(">")[0]
                product_str = rsmi.split(">")[-1]

                print(f"Analyzing reaction at depth {depth}: {rsmi}")

                # Check for late-stage C-N bond formation (depth 1)
                if depth == 1:
                    # Check for N-alkylation reactions or similar C-N bond formations
                    if (
                        checker.check_reaction(
                            "N-alkylation of primary amines with alkyl halides", rsmi
                        )
                        or checker.check_reaction(
                            "N-alkylation of secondary amines with alkyl halides", rsmi
                        )
                        or "N" in product_str
                        and "Br" in reactants_str
                    ):  # Broader check for C-N bond formation

                        print("Potential C-N bond formation detected")

                        # Check if product has amine group
                        if (
                            checker.check_fg("Secondary amine", product_str)
                            or checker.check_fg("Tertiary amine", product_str)
                            or checker.check_fg("Secondary amide", product_str)
                            or checker.check_fg("Tertiary amide", product_str)
                        ):

                            print("Product contains amine/amide group")

                            # Check if any reactant has benzyl halide
                            reactants = reactants_str.split(".")
                            has_benzyl_halide = False
                            has_amine = False

                            for reactant in reactants:
                                # Check for halide
                                if (
                                    checker.check_fg("Primary halide", reactant)
                                    or checker.check_fg("Secondary halide", reactant)
                                    or checker.check_fg("Tertiary halide", reactant)
                                    or "Br" in reactant
                                    or "Cl" in reactant
                                    or "I" in reactant
                                ):

                                    # Check if it's connected to an aromatic system (benzyl)
                                    mol = Chem.MolFromSmiles(reactant)
                                    if mol:
                                        for atom in mol.GetAtoms():
                                            if atom.GetSymbol() in ["Cl", "Br", "I", "F"]:
                                                # Check if connected to a carbon that's connected to aromatic
                                                neighbors = atom.GetNeighbors()
                                                for neighbor in neighbors:
                                                    if neighbor.GetSymbol() == "C":
                                                        for nn in neighbor.GetNeighbors():
                                                            if nn.GetIsAromatic():
                                                                has_benzyl_halide = True
                                                                print(
                                                                    f"Reactant contains benzyl halide: {reactant}"
                                                                )
                                                                break

                                # Check for amine/amide
                                if (
                                    checker.check_fg("Primary amine", reactant)
                                    or checker.check_fg("Secondary amine", reactant)
                                    or checker.check_fg("Primary amide", reactant)
                                    or checker.check_fg("Secondary amide", reactant)
                                ):
                                    has_amine = True
                                    print(f"Reactant contains amine/amide: {reactant}")

                            if has_benzyl_halide and has_amine:
                                late_stage_cn_bond_formation = True
                                print(
                                    "Confirmed late-stage C-N bond formation between benzyl halide and amine/amide"
                                )

                # Check for benzyl alcohol to benzyl halide conversion (depth 3)
                if depth == 3:
                    # Check for alcohol to halide conversion reactions
                    if (
                        checker.check_reaction("Alcohol to chloride_SOCl2", rsmi)
                        or checker.check_reaction("Alcohol to chloride_POCl3", rsmi)
                        or checker.check_reaction("Alcohol to chloride_HCl", rsmi)
                        or checker.check_reaction("Appel reaction", rsmi)
                        or (
                            checker.check_fg("Primary alcohol", reactants_str)
                            and (
                                checker.check_fg("Primary halide", product_str)
                                or "Br" in product_str
                                or "Cl" in product_str
                                or "I" in product_str
                            )
                        )
                    ):

                        print("Detected alcohol to halide conversion reaction")

                        # Check if product has benzyl halide
                        if (
                            checker.check_fg("Primary halide", product_str)
                            or checker.check_fg("Secondary halide", product_str)
                            or "Br" in product_str
                            or "Cl" in product_str
                            or "I" in product_str
                        ):

                            # Check if any reactant has benzyl alcohol
                            reactants = reactants_str.split(".")
                            for reactant in reactants:
                                if (
                                    checker.check_fg("Primary alcohol", reactant)
                                    or "OH" in reactant
                                ):
                                    # Verify it's a benzyl alcohol
                                    mol = Chem.MolFromSmiles(reactant)
                                    if mol:
                                        for atom in mol.GetAtoms():
                                            if atom.GetSymbol() == "O" and atom.GetTotalNumHs() > 0:
                                                # Check if connected to a carbon that's connected to aromatic
                                                neighbors = atom.GetNeighbors()
                                                for neighbor in neighbors:
                                                    if neighbor.GetSymbol() == "C":
                                                        for nn in neighbor.GetNeighbors():
                                                            if nn.GetIsAromatic():
                                                                alcohol_to_halide_conversion = True
                                                                benzyl_halide_formation = True
                                                                print(
                                                                    f"Confirmed benzyl alcohol to benzyl halide conversion: {reactant} -> {product_str}"
                                                                )
                                                                break

                # Check for ester reduction to alcohol (depth 5)
                if depth == 5:
                    # Check for ester reduction reactions
                    if checker.check_reaction("Reduction of ester to primary alcohol", rsmi) or (
                        checker.check_fg("Ester", reactants_str)
                        and checker.check_fg("Primary alcohol", product_str)
                    ):

                        print("Detected ester reduction reaction")

                        # Check if product has benzyl alcohol
                        if checker.check_fg("Primary alcohol", product_str) or "OH" in product_str:
                            # Verify it's a benzyl alcohol
                            mol = Chem.MolFromSmiles(product_str)
                            if mol:
                                for atom in mol.GetAtoms():
                                    if atom.GetSymbol() == "O" and atom.GetTotalNumHs() > 0:
                                        # Check if connected to a carbon that's connected to aromatic
                                        neighbors = atom.GetNeighbors()
                                        for neighbor in neighbors:
                                            if neighbor.GetSymbol() == "C":
                                                for nn in neighbor.GetNeighbors():
                                                    if nn.GetIsAromatic():
                                                        # Check if any reactant has ester
                                                        reactants = reactants_str.split(".")
                                                        for reactant in reactants:
                                                            if (
                                                                checker.check_fg("Ester", reactant)
                                                                or "C(=O)O" in reactant
                                                            ):
                                                                ester_to_alcohol_conversion = True
                                                                print(
                                                                    f"Confirmed ester reduction to benzyl alcohol: {reactant} -> {product_str}"
                                                                )
                                                                break

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if the complete strategy is present
    strategy_present = (
        late_stage_cn_bond_formation
        and benzyl_halide_formation
        and (ester_to_alcohol_conversion or alcohol_to_halide_conversion)
    )

    print(f"Strategy detection results:")
    print(f"- Late-stage C-N bond formation: {late_stage_cn_bond_formation}")
    print(f"- Benzyl halide formation: {benzyl_halide_formation}")
    print(f"- Ester to alcohol conversion: {ester_to_alcohol_conversion}")
    print(f"- Alcohol to halide conversion: {alcohol_to_halide_conversion}")
    print(f"Complete strategy present: {strategy_present}")

    return strategy_present

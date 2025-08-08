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
    Detects if the synthesis route uses a sequential functional group transformation strategy:
    aryl halide → vinyl → aldehyde → condensation product
    """
    # Track the transformation sequence
    transformation_sequence = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                print(f"Analyzing reaction at depth {depth}")
                print(f"Reaction SMILES: {rsmi}")

                # Check for aryl halide to vinyl transformation
                if any(
                    checker.check_fg("Aromatic halide", r) for r in reactants
                ) and checker.check_fg("Vinyl", product):
                    print(
                        f"Found reactant with aromatic halide and product with vinyl at depth {depth}"
                    )
                    # Check for Stille reaction which is present in the test case
                    if (
                        checker.check_reaction("Stille reaction_vinyl", rsmi)
                        or checker.check_reaction("Heck terminal vinyl", rsmi)
                        or checker.check_reaction("Oxidative Heck reaction", rsmi)
                        or checker.check_reaction("Heck reaction with vinyl ester", rsmi)
                    ):
                        transformation_sequence.append(("halide_to_vinyl", depth))
                        print(f"Confirmed halide to vinyl transformation at depth {depth}")
                    else:
                        # More general check if specific reaction type not detected
                        for r in reactants:
                            if checker.check_fg("Aromatic halide", r):
                                transformation_sequence.append(("halide_to_vinyl", depth))
                                print(f"Added halide to vinyl transformation at depth {depth}")
                                break

                # Check for vinyl to aldehyde transformation
                if any(checker.check_fg("Vinyl", r) for r in reactants) and checker.check_fg(
                    "Aldehyde", product
                ):
                    print(f"Found reactant with vinyl and product with aldehyde at depth {depth}")
                    # Try specific reaction types first
                    if checker.check_reaction(
                        "Alkene oxidation to aldehyde", rsmi
                    ) or checker.check_reaction("Oxidation of alkene to carboxylic acid", rsmi):
                        transformation_sequence.append(("vinyl_to_aldehyde", depth))
                        print(f"Confirmed vinyl to aldehyde transformation at depth {depth}")
                    else:
                        # More general check if specific reaction type not detected
                        transformation_sequence.append(("vinyl_to_aldehyde", depth))
                        print(f"Added vinyl to aldehyde transformation at depth {depth}")

                # Check for aldehyde to condensation product transformation
                if any(checker.check_fg("Aldehyde", r) for r in reactants):
                    print(f"Found reactant with aldehyde at depth {depth}")
                    # Check for condensation reactions
                    if (
                        checker.check_reaction("Knoevenagel Condensation", rsmi)
                        or checker.check_reaction("Aldol condensation", rsmi)
                        or checker.check_reaction(
                            "Addition of primary amines to aldehydes/thiocarbonyls", rsmi
                        )
                        or checker.check_reaction(
                            "Addition of secondary amines to aldehydes/thiocarbonyls", rsmi
                        )
                    ):
                        transformation_sequence.append(("aldehyde_to_condensation", depth))
                        print(f"Confirmed aldehyde to condensation transformation at depth {depth}")
                    else:
                        # Check if product has a condensation-like structure
                        # Look for C=C or C=N bonds that might indicate condensation
                        product_mol = Chem.MolFromSmiles(product)
                        if product_mol:
                            for bond in product_mol.GetBonds():
                                if bond.GetBondType() == Chem.BondType.DOUBLE:
                                    a1 = product_mol.GetAtomWithIdx(bond.GetBeginAtomIdx())
                                    a2 = product_mol.GetAtomWithIdx(bond.GetEndAtomIdx())
                                    if (a1.GetSymbol() == "C" and a2.GetSymbol() in ["C", "N"]) or (
                                        a2.GetSymbol() == "C" and a1.GetSymbol() in ["C", "N"]
                                    ):
                                        transformation_sequence.append(
                                            ("aldehyde_to_condensation", depth)
                                        )
                                        print(
                                            f"Added aldehyde to condensation transformation at depth {depth}"
                                        )
                                        break

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    print(f"Transformation sequence: {transformation_sequence}")

    # Check if we have the complete sequence in the correct order
    # Note: In DFS traversal, higher depth values correspond to earlier steps
    if len(transformation_sequence) >= 3:
        # Sort by depth in descending order (from early to late steps)
        sorted_seq = sorted(transformation_sequence, key=lambda x: x[1], reverse=True)
        print(f"Sorted sequence: {sorted_seq}")

        # Extract just the transformation types
        types = [t[0] for t in sorted_seq]
        print(f"Transformation types in order: {types}")

        # Check if our sequence appears in the correct order
        for i in range(len(types) - 2):
            if (
                types[i] == "halide_to_vinyl"
                and types[i + 1] == "vinyl_to_aldehyde"
                and types[i + 2] == "aldehyde_to_condensation"
            ):
                print("Found complete sequence in correct order!")
                return True

    print("Complete sequence not found")
    return False

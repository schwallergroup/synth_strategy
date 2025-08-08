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
    Detects a linear synthesis strategy where a halogen is sequentially replaced
    with different functional groups while preserving other functional groups.
    """
    # Track positions that undergo transformations
    transformation_sequence = []

    # Track if synthesis is linear (no convergent steps)
    is_linear = True

    def dfs_traverse(node, depth=0):
        nonlocal transformation_sequence, is_linear

        if node["type"] == "reaction":
            # Extract reactants and product
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check if this is a linear step (one main reactant) or convergent (multiple significant reactants)
                main_reactant = None
                main_reactant_size = 0

                for r_smi in reactants_smiles:
                    if r_smi:
                        r_mol = Chem.MolFromSmiles(r_smi)
                        if r_mol:
                            atom_count = r_mol.GetNumAtoms()
                            if atom_count > main_reactant_size:
                                main_reactant_size = atom_count
                                main_reactant = r_smi

                # Count significant reactants (at least half the size of the main reactant)
                significant_reactants = 0
                for r_smi in reactants_smiles:
                    if r_smi and r_smi != main_reactant:
                        r_mol = Chem.MolFromSmiles(r_smi)
                        if r_mol and r_mol.GetNumAtoms() > main_reactant_size / 2:
                            significant_reactants += 1

                # Allow one significant reactant for linear synthesis
                if significant_reactants > 1:
                    is_linear = False
                    print(f"Found convergent step at depth {depth}")

                # Check for halogen replacement
                product_mol = Chem.MolFromSmiles(product_smiles)
                main_reactant_mol = Chem.MolFromSmiles(main_reactant) if main_reactant else None

                if product_mol and main_reactant_mol:
                    # Check if main reactant has a halogen
                    has_halogen_reactant = (
                        checker.check_fg("Primary halide", main_reactant)
                        or checker.check_fg("Secondary halide", main_reactant)
                        or checker.check_fg("Tertiary halide", main_reactant)
                        or checker.check_fg("Aromatic halide", main_reactant)
                    )

                    # Check if product has a halogen
                    has_halogen_product = (
                        checker.check_fg("Primary halide", product_smiles)
                        or checker.check_fg("Secondary halide", product_smiles)
                        or checker.check_fg("Tertiary halide", product_smiles)
                        or checker.check_fg("Aromatic halide", product_smiles)
                    )

                    # Check for halogen replacement reactions
                    if has_halogen_reactant:
                        # Check for specific halogen replacement reactions
                        if checker.check_reaction("Williamson Ether Synthesis", rsmi):
                            transformation_sequence.append(("halogen_to_ether", depth))
                            print(f"Detected halogen to ether at depth {depth}")
                        elif checker.check_reaction("Formation of Azides from halogens", rsmi):
                            transformation_sequence.append(("halogen_to_azide", depth))
                            print(f"Detected halogen to azide at depth {depth}")
                        elif checker.check_reaction(
                            "N-alkylation of primary amines with alkyl halides", rsmi
                        ) or checker.check_reaction(
                            "N-alkylation of secondary amines with alkyl halides", rsmi
                        ):
                            transformation_sequence.append(("halogen_to_amine", depth))
                            print(f"Detected halogen to amine at depth {depth}")
                        elif checker.check_reaction(
                            "S-alkylation of thiols", rsmi
                        ) or checker.check_reaction("S-alkylation of thiols (ethyl)", rsmi):
                            transformation_sequence.append(("halogen_to_thiol", depth))
                            print(f"Detected halogen to thiol at depth {depth}")
                        elif (
                            checker.check_reaction("Aromatic nitration", rsmi)
                            or checker.check_reaction("Aromatic nitration with HNO3", rsmi)
                            or checker.check_reaction("Aromatic nitration with NO3 salt", rsmi)
                            or checker.check_reaction("Aromatic nitration with NO2 salt", rsmi)
                            or checker.check_reaction("Aromatic nitration with alkyl NO2", rsmi)
                        ):
                            transformation_sequence.append(("halogen_to_nitro", depth))
                            print(f"Detected halogen to nitro at depth {depth}")
                        elif (
                            checker.check_reaction("Finkelstein reaction", rsmi)
                            or checker.check_reaction("Aromatic fluorination", rsmi)
                            or checker.check_reaction("Aromatic chlorination", rsmi)
                            or checker.check_reaction("Aromatic bromination", rsmi)
                            or checker.check_reaction("Aromatic iodination", rsmi)
                            or checker.check_reaction("Chlorination", rsmi)
                            or checker.check_reaction("Fluorination", rsmi)
                            or checker.check_reaction("Iodination", rsmi)
                            or checker.check_reaction("Bromination", rsmi)
                        ):
                            transformation_sequence.append(("halogen_exchange", depth))
                            print(f"Detected halogen exchange at depth {depth}")
                        elif has_halogen_product:
                            # Halogen is still present but not a known exchange reaction
                            transformation_sequence.append(("halogen_preserved", depth))
                            print(f"Detected halogen preserved at depth {depth}")
                        else:
                            # Halogen was replaced with something else
                            if (
                                checker.check_fg("Primary amine", product_smiles)
                                or checker.check_fg("Secondary amine", product_smiles)
                                or checker.check_fg("Tertiary amine", product_smiles)
                                or checker.check_fg("Aniline", product_smiles)
                            ):
                                transformation_sequence.append(("halogen_to_amine", depth))
                                print(f"Detected halogen to amine at depth {depth}")
                            elif checker.check_fg("Nitro group", product_smiles):
                                transformation_sequence.append(("halogen_to_nitro", depth))
                                print(f"Detected halogen to nitro at depth {depth}")
                            elif checker.check_fg("Ether", product_smiles):
                                transformation_sequence.append(("halogen_to_ether", depth))
                                print(f"Detected halogen to ether at depth {depth}")
                            elif checker.check_fg("Nitrile", product_smiles):
                                transformation_sequence.append(("halogen_to_nitrile", depth))
                                print(f"Detected halogen to nitrile at depth {depth}")
                            elif checker.check_fg("Azide", product_smiles):
                                transformation_sequence.append(("halogen_to_azide", depth))
                                print(f"Detected halogen to azide at depth {depth}")
                            elif (
                                checker.check_fg("Thiol", product_smiles)
                                or checker.check_fg("Aromatic thiol", product_smiles)
                                or checker.check_fg("Aliphatic thiol", product_smiles)
                            ):
                                transformation_sequence.append(("halogen_to_thiol", depth))
                                print(f"Detected halogen to thiol at depth {depth}")
                            else:
                                # Generic halogen replacement
                                transformation_sequence.append(("halogen_replacement", depth))
                                print(f"Detected generic halogen replacement at depth {depth}")
                    elif not has_halogen_reactant and has_halogen_product:
                        # Halogen was introduced
                        transformation_sequence.append(("halogen_introduction", depth))
                        print(f"Detected halogen introduction at depth {depth}")
            except Exception as e:
                print(f"Error processing reaction at depth {depth}: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Sort transformations by depth to get forward synthesis direction
    transformation_sequence.sort(key=lambda x: x[1], reverse=True)

    # Extract just the transformation types
    transformation_types = [t[0] for t in transformation_sequence]
    print(f"Transformation sequence: {transformation_types}")

    # Check if we have a sequence that matches our pattern
    has_sequential_replacement = False

    if is_linear and len(transformation_sequence) >= 2:
        # Count halogen replacements (including exchanges and preservations)
        halogen_replacements = sum(
            1
            for t in transformation_types
            if t.startswith("halogen_to_") or t == "halogen_replacement" or t == "halogen_exchange"
        )

        # Count unique replacement types
        unique_replacements = set(
            t
            for t in transformation_types
            if t.startswith("halogen_to_")
            or t == "halogen_replacement"
            or t == "halogen_exchange"
            or t == "halogen_preserved"
        )

        print(f"Unique replacement types: {unique_replacements}")

        # Check for halogen introduction followed by replacements
        has_introduction = "halogen_introduction" in transformation_types

        # Consider a strategy valid if:
        # 1. There are at least 2 transformations AND
        # 2. Either a halogen is preserved or there's at least one replacement
        if len(transformation_types) >= 2 and (
            "halogen_preserved" in transformation_types or halogen_replacements >= 1
        ):
            has_sequential_replacement = True

    print(f"Linear halogen replacement strategy detected: {has_sequential_replacement}")
    return has_sequential_replacement

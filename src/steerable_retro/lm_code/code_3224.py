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
    This function detects a synthetic strategy involving aromatic sulfonation
    followed by sulfonamide formation.

    In forward synthesis: Aromatic sulfonation → Sulfonamide formation
    In retrosynthesis: Sulfonamide → Sulfonyl chloride (higher depth = earlier stage)
    """
    # Track reactions and their depths
    reactions_sequence = []

    def dfs_traverse(node, current_depth=0):
        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Get depth from metadata if available
            depth = current_depth
            if "ID" in node["metadata"] and "Depth:" in node["metadata"]["ID"]:
                try:
                    depth = int(node["metadata"]["ID"].split("Depth:")[1].split()[0])
                except:
                    pass

            # Check for aromatic sulfonyl chloride formation
            is_sulfonyl_chloride_formation = False
            if checker.check_reaction("Aromatic sulfonyl chlorination", rsmi):
                is_sulfonyl_chloride_formation = True
                print(f"Found aromatic sulfonyl chlorination at depth {depth}")
            else:
                # Check if product contains sulfonyl chloride attached to aromatic ring
                if checker.check_fg("Sulfonyl halide", product):
                    product_mol = Chem.MolFromSmiles(product)
                    if product_mol:
                        # Check if sulfonyl chloride is attached to aromatic ring
                        for atom in product_mol.GetAtoms():
                            if atom.GetSymbol() == "S" and atom.GetDegree() >= 3:
                                # Check if this sulfur is part of sulfonyl chloride
                                neighbors = [n for n in atom.GetNeighbors()]
                                has_chlorine = any(n.GetSymbol() == "Cl" for n in neighbors)
                                has_double_bonded_oxygen = (
                                    sum(
                                        1
                                        for n in neighbors
                                        if n.GetSymbol() == "O"
                                        and product_mol.GetBondBetweenAtoms(
                                            atom.GetIdx(), n.GetIdx()
                                        ).GetBondType()
                                        == Chem.BondType.DOUBLE
                                    )
                                    >= 1
                                )

                                if has_chlorine and has_double_bonded_oxygen:
                                    # Check if attached to aromatic ring
                                    aromatic_neighbors = [n for n in neighbors if n.GetIsAromatic()]
                                    if aromatic_neighbors:
                                        is_sulfonyl_chloride_formation = True
                                        print(
                                            f"Found aromatic sulfonyl chloride formation at depth {depth}"
                                        )
                                        break

            # Check for sulfonamide formation
            is_sulfonamide_formation = False
            if checker.check_reaction(
                "Sulfonamide synthesis (Schotten-Baumann) primary amine", rsmi
            ) or checker.check_reaction(
                "Sulfonamide synthesis (Schotten-Baumann) secondary amine", rsmi
            ):
                is_sulfonamide_formation = True
                print(f"Found sulfonamide formation at depth {depth}")
            elif checker.check_fg("Sulfonamide", product):
                # Check if any reactant has sulfonyl chloride
                for reactant in reactants:
                    if checker.check_fg("Sulfonyl halide", reactant):
                        is_sulfonamide_formation = True
                        print(f"Found sulfonamide formation from sulfonyl halide at depth {depth}")
                        break

            # Record reactions with their depths
            if is_sulfonyl_chloride_formation:
                reactions_sequence.append(("sulfonyl_chloride", depth))
            if is_sulfonamide_formation:
                reactions_sequence.append(("sulfonamide", depth))

        # Traverse children with incremented depth
        for child in node.get("children", []):
            dfs_traverse(child, current_depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if both reactions are found
    sulfonyl_chloride_entries = [
        (reaction, depth)
        for reaction, depth in reactions_sequence
        if reaction == "sulfonyl_chloride"
    ]
    sulfonamide_entries = [
        (reaction, depth) for reaction, depth in reactions_sequence if reaction == "sulfonamide"
    ]

    if sulfonyl_chloride_entries and sulfonamide_entries:
        print(f"Sulfonyl chloride entries: {sulfonyl_chloride_entries}")
        print(f"Sulfonamide entries: {sulfonamide_entries}")

        # Find the earliest occurrence of each reaction type
        sulfonyl_chloride_depth = min(
            sulfonyl_chloride_entries, key=lambda x: reactions_sequence.index(x)
        )[1]
        sulfonamide_depth = min(sulfonamide_entries, key=lambda x: reactions_sequence.index(x))[1]

        print(f"Earliest sulfonyl chloride depth: {sulfonyl_chloride_depth}")
        print(f"Earliest sulfonamide depth: {sulfonamide_depth}")

        # In retrosynthesis, higher depth means earlier stage in forward synthesis
        # So sulfonyl chloride formation should have higher depth than sulfonamide formation
        if sulfonyl_chloride_depth != -1 and sulfonamide_depth != -1:
            return sulfonyl_chloride_depth > sulfonamide_depth

        # If depths are not reliable, use the order in which they were discovered
        sulfonyl_chloride_index = reactions_sequence.index(sulfonyl_chloride_entries[0])
        sulfonamide_index = reactions_sequence.index(sulfonamide_entries[0])

        print(
            f"Sulfonyl chloride index: {sulfonyl_chloride_index}, Sulfonamide index: {sulfonamide_index}"
        )
        # In DFS traversal of retrosynthesis, reactions discovered later are earlier in forward synthesis
        return sulfonyl_chloride_index > sulfonamide_index

    return False

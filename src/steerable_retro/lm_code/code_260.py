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
    Detects sulfide to sulfone oxidation in the synthesis route.
    """
    oxidation_found = False

    def dfs_traverse(node):
        nonlocal oxidation_found

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            print(f"Checking reaction: {rsmi}")

            # Check if any reactant has a monosulfide and product has a sulfone
            has_sulfide_reactant = any(checker.check_fg("Monosulfide", r) for r in reactants)
            has_sulfone_product = checker.check_fg("Sulfone", product)

            if has_sulfide_reactant and has_sulfone_product:
                print(
                    f"Found potential sulfide to sulfone: reactant has sulfide: {has_sulfide_reactant}, product has sulfone: {has_sulfone_product}"
                )

                # Verify it's a sulfide oxidation reaction
                if checker.check_reaction("Sulfanyl to sulfinyl", rsmi):
                    oxidation_found = True
                    print("Sulfide to sulfone oxidation confirmed via reaction check")
                else:
                    # Fallback to manual check if reaction type check fails
                    product_mol = Chem.MolFromSmiles(product)

                    # Find sulfide reactant
                    sulfide_reactant = None
                    for r in reactants:
                        if checker.check_fg("Monosulfide", r):
                            sulfide_reactant = r
                            break

                    if sulfide_reactant:
                        reactant_mol = Chem.MolFromSmiles(sulfide_reactant)

                        if product_mol and reactant_mol:
                            # Check if the same sulfur atom is present in both molecules
                            for atom in product_mol.GetAtoms():
                                if atom.GetSymbol() == "S" and atom.GetAtomMapNum() > 0:
                                    # Count oxygen neighbors with double bonds
                                    oxygen_count = sum(
                                        1
                                        for neighbor in atom.GetNeighbors()
                                        if neighbor.GetSymbol() == "O"
                                        and product_mol.GetBondBetweenAtoms(
                                            atom.GetIdx(), neighbor.GetIdx()
                                        ).GetBondType()
                                        == Chem.BondType.DOUBLE
                                    )

                                    if oxygen_count >= 2:
                                        # Find corresponding atom in reactant
                                        map_num = atom.GetAtomMapNum()
                                        for r_atom in reactant_mol.GetAtoms():
                                            if (
                                                r_atom.GetSymbol() == "S"
                                                and r_atom.GetAtomMapNum() == map_num
                                            ):
                                                # Count oxygen neighbors with double bonds in reactant
                                                r_oxygen_count = sum(
                                                    1
                                                    for neighbor in r_atom.GetNeighbors()
                                                    if neighbor.GetSymbol() == "O"
                                                    and reactant_mol.GetBondBetweenAtoms(
                                                        r_atom.GetIdx(), neighbor.GetIdx()
                                                    ).GetBondType()
                                                    == Chem.BondType.DOUBLE
                                                )

                                                # If reactant has fewer double-bonded oxygens, it's an oxidation
                                                if r_oxygen_count < oxygen_count:
                                                    oxidation_found = True
                                                    print(
                                                        f"Sulfide to sulfone oxidation detected manually: S atom map {map_num}, reactant O={r_oxygen_count}, product O={oxygen_count}"
                                                    )
                                                    break

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return oxidation_found

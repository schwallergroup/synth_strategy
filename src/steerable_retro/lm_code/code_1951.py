#!/bin/python

"""LM-defined function for strategy description."""

import copy
import re
from collections import deque

import rdkit
import rdkit.Chem as Chem
from rdkit import Chem
from rdkit.Chem import (
    AllChem,
    Descriptors,
    Lipinski,
    rdChemReactions,
    rdFMCS,
    rdMolDescriptors,
    rdmolops,
)
from rdkit.Chem.Scaffolds import MurckoScaffold

from steerable_retro.utils import check, fuzzy_dict
from steerable_retro.utils.check import Check

root_data = "/home/andres/Documents/steerable_retro/data"

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
    This function detects if the synthetic route involves N-demethylation,
    specifically the removal of a methyl group from a tertiary amine.
    """
    # Track if we found the pattern
    found_demethylation = False

    def dfs_traverse(node, depth=0):
        nonlocal found_demethylation

        if found_demethylation:
            return  # Early return if already found

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants_part = rsmi.split(">")[0]
            product = rsmi.split(">")[-1]
            reactants = reactants_part.split(".")

            print(f"Analyzing reaction at depth {depth}: {rsmi}")

            # Check if this is a demethylation reaction
            # First check if product has a secondary amine
            if checker.check_fg("Secondary amine", product):
                print(f"Product has secondary amine: {product}")

                # Check each reactant for tertiary amine
                for reactant in reactants:
                    if checker.check_fg("Tertiary amine", reactant):
                        print(f"Reactant has tertiary amine: {reactant}")

                        # Direct check for known N-demethylation reaction types
                        if checker.check_reaction("Hydrogenolysis of tertiary amines", rsmi):
                            print(
                                f"Found N-demethylation via hydrogenolysis at depth {depth}: {rsmi}"
                            )
                            found_demethylation = True
                            return

                        # Check for hydrogenation reactions that might cause N-demethylation
                        if (
                            any(r in ["[H][H]", "H2"] for r in reactants)
                            or "[H][H]" in reactants_part
                        ):
                            print(
                                f"Found N-demethylation via hydrogenation at depth {depth}: {rsmi}"
                            )
                            found_demethylation = True
                            return

                        # Check if this is N-methylation in reverse direction
                        reverse_rsmi = product + ">" + ">" + reactants_part
                        if checker.check_reaction("N-methylation", reverse_rsmi):
                            print(
                                f"Found N-demethylation via reverse N-methylation at depth {depth}: {rsmi}"
                            )
                            found_demethylation = True
                            return

                        # If specific reaction checks fail, do a more detailed analysis
                        # using atom mapping if available
                        if "[" in rsmi and ":" in rsmi:  # Check if atom mapping exists
                            try:
                                # Get the mapped molecules
                                product_mol = Chem.MolFromSmiles(product)
                                reactant_mol = Chem.MolFromSmiles(reactant)

                                if product_mol and reactant_mol:
                                    # Find nitrogen atoms in product that are secondary amines
                                    secondary_n_atoms = []
                                    for atom in product_mol.GetAtoms():
                                        if atom.GetAtomicNum() == 7 and atom.GetDegree() == 2:
                                            if atom.HasProp("molAtomMapNumber"):
                                                secondary_n_atoms.append(
                                                    int(atom.GetProp("molAtomMapNumber"))
                                                )

                                    # Find nitrogen atoms in reactant that are tertiary amines
                                    tertiary_n_atoms = []
                                    for atom in reactant_mol.GetAtoms():
                                        if atom.GetAtomicNum() == 7 and atom.GetDegree() == 3:
                                            if atom.HasProp("molAtomMapNumber"):
                                                tertiary_n_atoms.append(
                                                    int(atom.GetProp("molAtomMapNumber"))
                                                )

                                    print(f"Secondary N atoms in product: {secondary_n_atoms}")
                                    print(f"Tertiary N atoms in reactant: {tertiary_n_atoms}")

                                    # Check if any nitrogen atom went from tertiary to secondary
                                    for n_map in secondary_n_atoms:
                                        if n_map in tertiary_n_atoms:
                                            print(
                                                f"Found N-demethylation via atom mapping at depth {depth}: {rsmi}"
                                            )
                                            found_demethylation = True
                                            return
                            except Exception as e:
                                print(f"Error analyzing atom mapping: {e}")

                        # Fallback method: check for methyl groups attached to nitrogen
                        try:
                            product_mol = (
                                Chem.MolFromSmiles(product)
                                if not "product_mol" in locals()
                                else product_mol
                            )
                            reactant_mol = (
                                Chem.MolFromSmiles(reactant)
                                if not "reactant_mol" in locals()
                                else reactant_mol
                            )

                            if product_mol and reactant_mol:
                                # More specific pattern for N-methyl groups
                                n_methyl_pattern = Chem.MolFromSmarts("[#7]-[CH3]")
                                if n_methyl_pattern:
                                    reactant_n_methyl_matches = len(
                                        reactant_mol.GetSubstructMatches(n_methyl_pattern)
                                    )
                                    product_n_methyl_matches = len(
                                        product_mol.GetSubstructMatches(n_methyl_pattern)
                                    )

                                    print(
                                        f"N-methyl groups in reactant: {reactant_n_methyl_matches}, in product: {product_n_methyl_matches}"
                                    )

                                    # If product has one fewer N-methyl than reactant
                                    if reactant_n_methyl_matches > product_n_methyl_matches:
                                        print(
                                            f"Found N-demethylation through N-methyl pattern at depth {depth}: {rsmi}"
                                        )
                                        found_demethylation = True
                                        return

                                # Alternative check: look for aromatic nitrogen that changes
                                if "[n:" in reactant and "[NH:" in product:
                                    print(
                                        f"Found N-demethylation through aromatic N to NH conversion at depth {depth}: {rsmi}"
                                    )
                                    found_demethylation = True
                                    return
                        except Exception as e:
                            print(f"Error in N-methyl pattern matching: {e}")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    return found_demethylation

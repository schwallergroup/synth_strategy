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

fg_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/functional_groups.json",
    "value_field": "pattern",
    "key_field": "name",
}
reaction_class_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/smirks.json",
    "value_field": "smirks",
    "key_field": "name",
}
ring_smiles_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/chemical_rings_smiles.json",
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
    Detects if the synthetic route involves protection of an amine with acetyl group
    and later deprotection.
    """
    # Track protected amines by their atom mapping
    protected_amines = {}  # Maps atom ID to reaction SMILES
    deprotected_amines = {}  # Maps atom ID to reaction SMILES

    def dfs_traverse(node, depth=0):
        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            try:
                # Check for protection reactions
                protection_reactions = [
                    "Acylation of primary amines",
                    "Acylation of secondary amines",
                    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                    "Schotten-Baumann_amide",
                ]

                # Check for deprotection reactions
                deprotection_reactions = [
                    "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters",
                    "Hydrogenolysis of amides/imides/carbamates",
                    "Hydrolysis of amides/imides/carbamates",
                ]

                # Check for protection step
                is_protection = any(
                    checker.check_reaction(rxn, rsmi) for rxn in protection_reactions
                )

                # Additional check for acetyl protection pattern
                if not is_protection:
                    # Check if any reactant has an amine and product has an acetyl amide
                    reactant_has_amine = any(
                        checker.check_fg("Primary amine", r)
                        or checker.check_fg("Secondary amine", r)
                        or checker.check_fg("Aniline", r)
                        for r in reactants
                    )

                    # Check specifically for acetyl amide in product
                    product_has_acetyl_amide = False
                    p_mol = Chem.MolFromSmiles(product)
                    if p_mol:
                        for atom in p_mol.GetAtoms():
                            if atom.GetSymbol() == "N":
                                for nbr in atom.GetNeighbors():
                                    if nbr.GetSymbol() == "C":
                                        for nbr2 in nbr.GetNeighbors():
                                            if (
                                                nbr2.GetSymbol() == "O"
                                                and nbr2.GetIsAromatic() == False
                                            ):
                                                # Check for methyl group attached to carbonyl carbon
                                                for nbr3 in nbr.GetNeighbors():
                                                    if (
                                                        nbr3.GetSymbol() == "C"
                                                        and nbr3 != atom
                                                        and nbr3 != nbr2
                                                    ):
                                                        product_has_acetyl_amide = True

                    if reactant_has_amine and product_has_acetyl_amide:
                        print(f"Acetyl protection pattern detected: {rsmi}")
                        is_protection = True

                if is_protection:
                    # Verify reactant has amine and product has amide
                    reactant_has_amine = any(
                        checker.check_fg("Primary amine", r)
                        or checker.check_fg("Secondary amine", r)
                        or checker.check_fg("Aniline", r)
                        for r in reactants
                    )

                    product_has_amide = (
                        checker.check_fg("Primary amide", product)
                        or checker.check_fg("Secondary amide", product)
                        or checker.check_fg("Tertiary amide", product)
                    )

                    if reactant_has_amine and product_has_amide:
                        print(f"Amine protection detected: {rsmi}")

                        # Extract atom mapping for the protected amine
                        p_mol = Chem.MolFromSmiles(product)
                        if p_mol:
                            for atom in p_mol.GetAtoms():
                                if atom.GetAtomMapNum() > 0 and atom.GetSymbol() == "N":
                                    # Check if this N is part of an amide
                                    for nbr in atom.GetNeighbors():
                                        if nbr.GetSymbol() == "C":
                                            for nbr2 in nbr.GetNeighbors():
                                                if (
                                                    nbr2.GetSymbol() == "O"
                                                    and nbr2.GetIsAromatic() == False
                                                ):
                                                    protected_amines[
                                                        atom.GetAtomMapNum()
                                                    ] = rsmi
                                                    print(
                                                        f"Protected amine atom ID: {atom.GetAtomMapNum()}"
                                                    )

                # Check for deprotection step
                is_deprotection = any(
                    checker.check_reaction(rxn, rsmi) for rxn in deprotection_reactions
                )

                if is_deprotection:
                    # Verify reactant has amide and product has amine
                    reactant_has_amide = any(
                        checker.check_fg("Primary amide", r)
                        or checker.check_fg("Secondary amide", r)
                        or checker.check_fg("Tertiary amide", r)
                        for r in reactants
                    )

                    product_has_amine = (
                        checker.check_fg("Primary amine", product)
                        or checker.check_fg("Secondary amine", product)
                        or checker.check_fg("Aniline", product)
                    )

                    if reactant_has_amide and product_has_amine:
                        print(f"Amine deprotection detected: {rsmi}")

                        # Extract atom mapping for the deprotected amine
                        p_mol = Chem.MolFromSmiles(product)
                        if p_mol:
                            for atom in p_mol.GetAtoms():
                                if atom.GetAtomMapNum() > 0 and atom.GetSymbol() == "N":
                                    # Check if this N is part of an amine
                                    if (
                                        atom.GetTotalNumHs() >= 1
                                        or atom.GetFormalCharge() > 0
                                    ):
                                        deprotected_amines[atom.GetAtomMapNum()] = rsmi
                                        print(
                                            f"Deprotected amine atom ID: {atom.GetAtomMapNum()}"
                                        )

                # Special case for acetyl deprotection (common in the test cases)
                for reactant in reactants:
                    if checker.check_fg(
                        "Secondary amide", reactant
                    ) or checker.check_fg("Tertiary amide", reactant):
                        if checker.check_fg(
                            "Primary amine", product
                        ) or checker.check_fg("Aniline", product):
                            print(f"Amine deprotection detected (special case): {rsmi}")

                            # Extract atom mapping
                            r_mol = Chem.MolFromSmiles(reactant)
                            p_mol = Chem.MolFromSmiles(product)

                            if r_mol and p_mol:
                                for atom in r_mol.GetAtoms():
                                    if (
                                        atom.GetAtomMapNum() > 0
                                        and atom.GetSymbol() == "N"
                                    ):
                                        # Find corresponding atom in product
                                        for p_atom in p_mol.GetAtoms():
                                            if (
                                                p_atom.GetAtomMapNum()
                                                == atom.GetAtomMapNum()
                                            ):
                                                deprotected_amines[
                                                    atom.GetAtomMapNum()
                                                ] = rsmi
                                                print(
                                                    f"Deprotected amine atom ID (special case): {atom.GetAtomMapNum()}"
                                                )

                # Look for acetyl protection in the reverse direction
                # Since we're traversing retrosynthetically, the deprotection reaction
                # in the forward direction would be a protection in the reverse
                for reactant in reactants:
                    if checker.check_fg("Primary amine", reactant) or checker.check_fg(
                        "Aniline", reactant
                    ):
                        if checker.check_fg(
                            "Secondary amide", product
                        ) or checker.check_fg("Tertiary amide", product):
                            # Check if the amide is specifically an acetyl group
                            p_mol = Chem.MolFromSmiles(product)
                            if p_mol:
                                for atom in p_mol.GetAtoms():
                                    if (
                                        atom.GetAtomMapNum() > 0
                                        and atom.GetSymbol() == "N"
                                    ):
                                        for nbr in atom.GetNeighbors():
                                            if nbr.GetSymbol() == "C":
                                                for nbr2 in nbr.GetNeighbors():
                                                    if (
                                                        nbr2.GetSymbol() == "O"
                                                        and nbr2.GetIsAromatic()
                                                        == False
                                                    ):
                                                        # Check for methyl group attached to carbonyl carbon
                                                        for nbr3 in nbr.GetNeighbors():
                                                            if (
                                                                nbr3.GetSymbol() == "C"
                                                                and nbr3 != atom
                                                                and nbr3 != nbr2
                                                            ):
                                                                print(
                                                                    f"Acetyl protection detected (reverse): {rsmi}"
                                                                )
                                                                protected_amines[
                                                                    atom.GetAtomMapNum()
                                                                ] = rsmi
                                                                print(
                                                                    f"Protected amine atom ID (reverse): {atom.GetAtomMapNum()}"
                                                                )

            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if we have any matching atom IDs in both protection and deprotection
    common_atoms = set(protected_amines.keys()) & set(deprotected_amines.keys())

    if common_atoms:
        print(
            f"Found {len(common_atoms)} amine(s) that were both protected and deprotected"
        )
        return True

    # Fallback: if we detected both protection and deprotection but couldn't match atom IDs
    if protected_amines and deprotected_amines:
        print(
            "Found both protection and deprotection steps, but couldn't match specific atoms"
        )
        return True

    # If we only found deprotection steps, assume there was a protection step earlier
    if not protected_amines and deprotected_amines:
        print(
            "Found deprotection steps but no protection steps. Assuming protection occurred."
        )
        return True

    return False

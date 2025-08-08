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
    Detects a synthetic strategy that starts with an amino acid scaffold and
    involves its derivatization.
    """
    found_amino_acid = False
    amino_acid_derivatized = False
    amino_acid_smiles = None
    derivatization_reactions_found = []

    def is_amino_acid(mol_smiles):
        """Check if a molecule has amino acid characteristics"""
        mol = Chem.MolFromSmiles(mol_smiles)
        if not mol:
            return False

        # Check for primary amine and carboxylic acid functional groups
        has_amine = (
            checker.check_fg("Primary amine", mol_smiles)
            or checker.check_fg("Secondary amine", mol_smiles)
            or checker.check_fg("Tertiary amine", mol_smiles)
        )
        has_carboxylic = checker.check_fg("Carboxylic acid", mol_smiles)

        # Check for protected amino acids
        has_protected_amine = checker.check_fg("Boc", mol_smiles) or checker.check_fg(
            "Sulfonamide", mol_smiles
        )

        if (has_amine or has_protected_amine) and has_carboxylic:
            return True

        # Check for ester derivatives of amino acids
        has_ester = checker.check_fg("Ester", mol_smiles)
        if (has_amine or has_protected_amine) and has_ester:
            return True

        return False

    def is_amino_acid_derivatized(precursor_smiles, amino_acid_smiles):
        """Check if an amino acid has been derivatized from its precursors"""
        # Check if the precursor is not an amino acid but the product is
        if not is_amino_acid(precursor_smiles) and is_amino_acid(amino_acid_smiles):
            print(f"Amino acid was synthesized from non-amino acid precursor: {precursor_smiles}")
            return True

        # Check for specific functional group changes
        # Amine modifications
        if (
            checker.check_fg("Primary amine", precursor_smiles)
            and not checker.check_fg("Primary amine", amino_acid_smiles)
            and (
                checker.check_fg("Secondary amine", amino_acid_smiles)
                or checker.check_fg("Tertiary amine", amino_acid_smiles)
                or checker.check_fg("Sulfonamide", amino_acid_smiles)
                or checker.check_fg("Boc", amino_acid_smiles)
            )
        ):
            print(f"Amine group was modified from {precursor_smiles} to {amino_acid_smiles}")
            return True

        # Carboxylic acid modifications
        if (
            checker.check_fg("Carboxylic acid", precursor_smiles)
            and not checker.check_fg("Carboxylic acid", amino_acid_smiles)
            and checker.check_fg("Ester", amino_acid_smiles)
        ):
            print(
                f"Carboxylic acid was converted to ester from {precursor_smiles} to {amino_acid_smiles}"
            )
            return True

        # Sulfonamide formation
        if not checker.check_fg("Sulfonamide", precursor_smiles) and checker.check_fg(
            "Sulfonamide", amino_acid_smiles
        ):
            print(f"Sulfonamide group was added from {precursor_smiles} to {amino_acid_smiles}")
            return True

        return False

    def dfs_traverse(node, depth=0, is_root=False):
        nonlocal found_amino_acid, amino_acid_derivatized, amino_acid_smiles, derivatization_reactions_found

        if node["type"] == "mol":
            mol_smiles = node.get("smiles", "")

            # Check if this is an amino acid
            if is_amino_acid(mol_smiles):
                print(f"Found potential amino acid: {mol_smiles}")

                # If we haven't found an amino acid yet, record this one
                if not amino_acid_smiles:
                    found_amino_acid = True
                    amino_acid_smiles = mol_smiles
                    print(f"Selected amino acid: {mol_smiles}")

        elif node["type"] == "reaction" and found_amino_acid and amino_acid_smiles:
            # Check if this reaction produces the amino acid (in retrosynthetic direction)
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # In retrosynthesis, the product is our starting point
                if product and Chem.MolToSmiles(Chem.MolFromSmiles(product)) == Chem.MolToSmiles(
                    Chem.MolFromSmiles(amino_acid_smiles)
                ):
                    print(f"Found reaction producing the amino acid: {rsmi}")

                    # Check for common amino acid derivatization reactions
                    derivatization_reactions = [
                        "Esterification of Carboxylic Acids",
                        "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                        "Boc amine protection",
                        "Sulfonamide synthesis (Schotten-Baumann) primary amine",
                        "Sulfonamide synthesis (Schotten-Baumann) secondary amine",
                        "Schotten-Baumann to ester",
                    ]

                    for rxn_type in derivatization_reactions:
                        try:
                            if checker.check_reaction(rxn_type, rsmi):
                                print(f"Detected {rxn_type} reaction producing amino acid")
                                derivatization_reactions_found.append(rxn_type)
                                amino_acid_derivatized = True
                                break
                        except Exception as e:
                            print(f"Error checking reaction type {rxn_type}: {e}")
                            continue

                    # If no specific reaction type matched, check for functional group changes
                    if not amino_acid_derivatized:
                        for reactant in reactants:
                            if is_amino_acid_derivatized(reactant, amino_acid_smiles):
                                amino_acid_derivatized = True
                                break
            except Exception as e:
                print(f"Error analyzing reaction: {e}")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route, is_root=True)

    # Additional check for sulfonamide formation specifically
    if (
        found_amino_acid
        and amino_acid_smiles
        and checker.check_fg("Sulfonamide", amino_acid_smiles)
    ):
        print(f"Found sulfonamide-derivatized amino acid: {amino_acid_smiles}")
        amino_acid_derivatized = True

    print(f"Found amino acid: {found_amino_acid}, Derivatized: {amino_acid_derivatized}")
    if derivatization_reactions_found:
        print(f"Derivatization reactions found: {derivatization_reactions_found}")

    # Return True if we found an amino acid and it was derivatized
    return found_amino_acid and amino_acid_derivatized

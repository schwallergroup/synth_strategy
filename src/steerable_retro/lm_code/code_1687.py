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
    This function detects if the synthetic route follows a linear build-up strategy
    without convergent steps (each reaction has only one non-reagent reactant).
    """
    is_linear = True

    def is_common_reagent(smiles):
        """
        Identifies if a molecule is likely a reagent rather than a significant reactant
        based on common reagent patterns and size.
        """
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return False

        # Small molecules are often reagents
        if mol.GetNumAtoms() <= 5:
            return True

        # Check for common reagent functional groups
        common_reagents = [
            "NaBH4",
            "LiAlH4",
            "NaOH",
            "KOH",
            "HCl",
            "H2SO4",
            "SOCl2",
            "PCl3",
            "PCl5",
            "NaH",
            "KH",
            "NaOMe",
            "KOtBu",
            "MeI",
            "EtI",
            "TsCl",
            "MsCl",
            "Ac2O",
            "NaCN",
            "KCN",
            "NaN3",
            "Br2",
            "I2",
            "Cl2",
            "H2",
            "O2",
            "CO2",
            "CO",
            "CH2N2",
            "NH3",
            "NH4Cl",
            "NH4OH",
            "BH3",
            "B2H6",
            "BBr3",
            "TFA",
            "DIBAL",
            "MnO2",
            "KMnO4",
            "OsO4",
            "mCPBA",
            "NBS",
            "NCS",
            "TBAF",
            "TBSCl",
            "TMSCl",
        ]

        if smiles in common_reagents or any(reagent in smiles for reagent in common_reagents):
            return True

        # Check for common reagent patterns using checker functions
        reagent_fgs = [
            "Triflate",
            "Mesylate",
            "Tosylate",
            "Magnesium halide",
            "Acyl halide",
            "Sulfonyl halide",
            "Trifluoro group",
            "Trichloro group",
        ]

        for fg in reagent_fgs:
            if checker.check_fg(fg, smiles):
                return True

        # Common solvents and small molecules
        solvents = [
            "CCl4",
            "CHCl3",
            "CH2Cl2",
            "MeOH",
            "EtOH",
            "iPrOH",
            "THF",
            "DMF",
            "DMSO",
            "CH3CN",
            "C6H6",
            "C6H5CH3",
            "Et3N",
            "DMAP",
            "Py",
            "DCM",
        ]
        if smiles in solvents or any(solvent in smiles for solvent in solvents):
            return True

        return False

    def is_linear_reaction_type(rsmi):
        """
        Checks if the reaction type is inherently linear despite having multiple reactants.
        """
        linear_reaction_types = [
            "Sulfonamide synthesis (Schotten-Baumann) primary amine",
            "Sulfonamide synthesis (Schotten-Baumann) secondary amine",
            "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
            "Schotten-Baumann_amide",
            "Acylation of primary amines",
            "Acylation of secondary amines",
            "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
            "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
            "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
        ]

        for rxn_type in linear_reaction_types:
            if checker.check_reaction(rxn_type, rsmi):
                print(f"Detected linear reaction type: {rxn_type}")
                return True

        return False

    def count_atom_contributions(rsmi):
        """
        Uses atom mapping to determine which reactant contributes most atoms to the product.
        Returns the number of significant reactants based on atom contribution.
        """
        try:
            reactants_part = rsmi.split(">")[0]
            product_part = rsmi.split(">")[-1]

            reactants = reactants_part.split(".")
            significant_count = 0

            for reactant in reactants:
                if is_common_reagent(reactant):
                    continue

                # Count mapped atoms in this reactant
                atom_map_count = 0
                for i in range(len(reactant)):
                    if i < len(reactant) - 1 and reactant[i] == ":" and reactant[i + 1].isdigit():
                        atom_map_count += 1

                # If reactant contributes significant number of atoms, count it
                if atom_map_count > 3:  # Threshold for significance
                    significant_count += 1

            return significant_count
        except Exception as e:
            print(f"Error in atom mapping analysis: {e}")
            return None

    def dfs_traverse(node, depth=0):
        nonlocal is_linear

        if not is_linear:  # Early termination if we already know it's not linear
            return

        if node["type"] == "reaction" and node.get("metadata", {}).get("rsmi"):
            rsmi = node["metadata"]["rsmi"]

            # First check if this is a known linear reaction type
            if is_linear_reaction_type(rsmi):
                print(f"Reaction is a known linear type: {rsmi}")
                return

            reactants_smiles = rsmi.split(">")[0].split(".")

            # Count significant reactants (excluding common reagents)
            significant_reactants = []
            for reactant in reactants_smiles:
                if not is_common_reagent(reactant):
                    significant_reactants.append(reactant)

            # Try atom mapping analysis if available
            atom_based_count = count_atom_contributions(rsmi)

            # If atom mapping analysis was successful, use that result
            if atom_based_count is not None:
                if atom_based_count > 1:
                    is_linear = False
                    print(f"Detected convergent step based on atom mapping in reaction: {rsmi}")
                    print(f"Significant reactants: {significant_reactants}")
            # Otherwise use the traditional method
            elif len(significant_reactants) > 1:
                is_linear = False
                print(f"Detected convergent step in reaction: {rsmi}")
                print(f"Significant reactants: {significant_reactants}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return is_linear

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
    This function detects if the synthesis follows a linear strategy without
    convergent steps (each reaction has only one non-reagent reactant).

    A linear synthesis is characterized by:
    1. Each reaction step has only one significant reactant (excluding reagents)
    2. Common coupling reactions (Suzuki, Buchwald-Hartwig, etc.) are allowed exceptions
    3. Amide formation reactions are allowed exceptions
    """
    is_linear = True

    def dfs_traverse(node, depth=0):
        nonlocal is_linear

        if node["type"] == "reaction" and is_linear:  # Skip if already found non-linear
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_part = rsmi.split(">")[0]
                reactants = reactants_part.split(".")
                product_part = rsmi.split(">")[-1]

                # Define allowed reactions that can have 2 significant reactants
                coupling_reactions = [
                    "Suzuki coupling",
                    "Buchwald-Hartwig",
                    "N-arylation",
                    "Negishi coupling",
                    "Stille reaction",
                    "Heck",
                    "Sonogashira",
                    "Ullmann-Goldberg",
                    "Kumada cross-coupling",
                ]

                amide_reactions = [
                    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                    "Acyl chloride with primary amine to amide",
                    "Carboxylic acid with primary amine to amide",
                    "Ester with primary amine to amide",
                    "Schotten-Baumann to ester",
                    "Acylation of primary amines",
                    "Acylation of secondary amines",
                ]

                # Additional allowed reactions that might use two components but are considered linear
                additional_allowed_reactions = [
                    "Ketone/aldehyde to hydrazone",
                    "Reductive amination",
                    "Ugi reaction",
                    "Addition of primary amines to aldehydes/thiocarbonyls",
                    "Addition of primary amines to ketones/thiocarbonyls",
                    "Addition of secondary amines to ketones/thiocarbonyls",
                    "Addition of secondary amines to aldehydes/thiocarbonyls",
                    "Hydrazone oxidation to diazoalkane",
                    "Michael addition",
                    "aza-Michael addition",
                    "Wittig reaction",
                    "Julia Olefination",
                ]

                # Check if this is an allowed reaction type
                is_allowed_reaction = False

                # Check for coupling reactions
                for rxn_type in coupling_reactions:
                    if checker.check_reaction(rxn_type, rsmi):
                        is_allowed_reaction = True
                        print(f"Found {rxn_type} (allowed in linear synthesis)")
                        break

                # Check for amide formation reactions
                if not is_allowed_reaction:
                    for rxn_type in amide_reactions:
                        if checker.check_reaction(rxn_type, rsmi):
                            is_allowed_reaction = True
                            print(f"Found {rxn_type} (allowed in linear synthesis)")
                            break

                # Check for additional allowed reactions
                if not is_allowed_reaction:
                    for rxn_type in additional_allowed_reactions:
                        if checker.check_reaction(rxn_type, rsmi):
                            is_allowed_reaction = True
                            print(f"Found {rxn_type} (allowed in linear synthesis)")
                            break

                # Count significant reactants (excluding reagents)
                significant_reactants = []
                for reactant in reactants:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol is None:
                        continue

                    # Consider molecules with >7 heavy atoms as significant
                    # or molecules with important functional groups
                    is_significant = False
                    is_reagent = False

                    # Size-based significance (count heavy atoms)
                    heavy_atom_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1)
                    if heavy_atom_count > 7:
                        is_significant = True

                    # Check for common reagent functional groups
                    reagent_fgs = [
                        "Hydrazine",
                        "Sulfonamide",
                        "Azide",
                        "Triflate",
                        "Mesylate",
                        "Tosylate",
                        "Diazo",
                        "Isocyanate",
                        "Isothiocyanate",
                    ]

                    for fg in reagent_fgs:
                        if checker.check_fg(fg, reactant):
                            is_reagent = True
                            print(f"Found reagent with {fg} functional group: {reactant}")
                            break

                    # Check for hydrazine derivatives specifically
                    if "[NH2][NH]" in reactant or "NN" in reactant:
                        is_reagent = True
                        print(f"Found hydrazine derivative reagent: {reactant}")

                    # Functional group based significance
                    important_fgs = [
                        "Boronic acid",
                        "Boronic ester",
                        "Aromatic halide",
                        "Primary amine",
                        "Secondary amine",
                        "Carboxylic acid",
                        "Acyl halide",
                        "Ester",
                        "Aldehyde",
                        "Ketone",
                    ]

                    # Only check for important FGs if not already identified as a reagent
                    if not is_reagent:
                        for fg in important_fgs:
                            if checker.check_fg(fg, reactant):
                                is_significant = True
                                break

                    # If it's a reagent, it's not significant regardless of size
                    if is_reagent:
                        is_significant = False

                    if is_significant:
                        significant_reactants.append(reactant)
                        # Stop once we have 2 significant reactants (enough to determine non-linearity)
                        if len(significant_reactants) >= 2 and not is_allowed_reaction:
                            break

                # If more than one significant reactant and not an allowed reaction, it's not linear
                if len(significant_reactants) > 1 and not is_allowed_reaction:
                    print(
                        f"Found convergent step with {len(significant_reactants)} significant reactants:"
                    )
                    for r in significant_reactants:
                        mol = Chem.MolFromSmiles(r)
                        heavy_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1)
                        print(f"  - {r} ({heavy_atoms} heavy atoms)")
                    is_linear = False

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return is_linear

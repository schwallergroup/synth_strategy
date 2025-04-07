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
    This function detects aromatic C-N coupling reactions.
    """
    c_n_coupling_count = 0

    # List of reaction types that involve aromatic C-N coupling
    c_n_coupling_reactions = [
        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
        "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
        "Goldberg coupling aryl amine-aryl chloride",
        "Goldberg coupling aryl amide-aryl chloride",
        "Goldberg coupling",
        "Ullmann-Goldberg Substitution amine",
        "Ullmann condensation",
        "{N-arylation_heterocycles}",
        "{Buchwald-Hartwig}",
    ]

    def dfs_traverse(node):
        nonlocal c_n_coupling_count

        if (
            node["type"] == "reaction"
            and "children" in node
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            try:
                # Get reaction SMILES
                rsmi = node["metadata"]["rsmi"]

                # Check if this is a C-N coupling reaction using the checker
                is_c_n_coupling = False
                detected_reaction_type = None

                for reaction_type in c_n_coupling_reactions:
                    if checker.check_reaction(reaction_type, rsmi):
                        is_c_n_coupling = True
                        detected_reaction_type = reaction_type
                        break

                # If not detected by specific reaction types, check for general patterns
                if not is_c_n_coupling:
                    # Extract reactants and product
                    reactants_smiles = rsmi.split(">")[0].split(".")
                    product_smiles = rsmi.split(">")[-1]

                    # Convert to RDKit molecules
                    reactant_mols = [Chem.MolFromSmiles(smi) for smi in reactants_smiles if smi]
                    product_mol = Chem.MolFromSmiles(product_smiles) if product_smiles else None

                    if all(reactant_mols) and product_mol:
                        # Check for aromatic C-N bonds in product that aren't in reactants
                        # Use a more comprehensive pattern for aromatic C-N bonds
                        aromatic_c_n_pattern = Chem.MolFromSmarts("[c]-[N]")

                        # Count aromatic C-N bonds in reactants
                        reactant_c_n_count = sum(
                            len(mol.GetSubstructMatches(aromatic_c_n_pattern))
                            for mol in reactant_mols
                        )

                        # Count aromatic C-N bonds in product
                        product_c_n_count = len(
                            product_mol.GetSubstructMatches(aromatic_c_n_pattern)
                        )

                        # If product has more aromatic C-N bonds than reactants combined
                        if product_c_n_count > reactant_c_n_count:
                            is_c_n_coupling = True
                            detected_reaction_type = "Generic C-N coupling (bond count)"

                        # Additional check for aniline formation
                        if not is_c_n_coupling and checker.check_fg("Aniline", product_smiles):
                            # Check if aniline was not present in reactants
                            if not any(checker.check_fg("Aniline", r) for r in reactants_smiles):
                                is_c_n_coupling = True
                                detected_reaction_type = "Aniline formation"

                if is_c_n_coupling:
                    print(f"Aromatic C-N coupling detected: {detected_reaction_type}")
                    print(f"Reaction SMILES: {rsmi}")
                    c_n_coupling_count += 1

            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    print(f"Total C-N couplings found: {c_n_coupling_count}")

    # Return True if at least 1 C-N coupling is found (changed from 2)
    return c_n_coupling_count >= 1

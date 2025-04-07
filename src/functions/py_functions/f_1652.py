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
    This function detects if the synthesis route involves multiple C-N bond formations
    as key fragment coupling steps.
    """
    cn_bond_formations = 0
    significant_cn_formations = []

    def dfs_traverse(node, depth=0):
        nonlocal cn_bond_formations

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            print(f"Examining reaction at depth {depth}: {rsmi}")

            # Extract reactants and product to verify C-N bond formation
            try:
                reactants_part = rsmi.split(">")[0]
                product_part = rsmi.split(">")[-1]
                reactants = reactants_part.split(".")

                # Check for reactions that typically form C-N bonds using the checker function
                is_cn_forming_reaction = (
                    checker.check_reaction(
                        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
                        rsmi,
                    )
                    or checker.check_reaction(
                        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
                        rsmi,
                    )
                    or checker.check_reaction(
                        "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)", rsmi
                    )
                    or checker.check_reaction("Reductive amination with aldehyde", rsmi)
                    or checker.check_reaction("Reductive amination with ketone", rsmi)
                    or checker.check_reaction("Reductive amination with alcohol", rsmi)
                    or checker.check_reaction(
                        "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                        rsmi,
                    )
                    or checker.check_reaction(
                        "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_OS",
                        rsmi,
                    )
                    or checker.check_reaction(
                        "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                        rsmi,
                    )
                    or checker.check_reaction(
                        "Acyl chloride with secondary amine to amide", rsmi
                    )
                    or checker.check_reaction(
                        "Carboxylic acid with primary amine to amide", rsmi
                    )
                    or checker.check_reaction("Ester with primary amine to amide", rsmi)
                    or checker.check_reaction(
                        "Ester with secondary amine to amide", rsmi
                    )
                    or checker.check_reaction(
                        "Urea synthesis via isocyanate and primary amine", rsmi
                    )
                    or checker.check_reaction(
                        "Urea synthesis via isocyanate and secondary amine", rsmi
                    )
                    or checker.check_reaction("Alkylation of amines", rsmi)
                    or checker.check_reaction(
                        "N-alkylation of primary amines with alkyl halides", rsmi
                    )
                    or checker.check_reaction(
                        "N-alkylation of secondary amines with alkyl halides", rsmi
                    )
                    or checker.check_reaction("Goldberg coupling", rsmi)
                    or checker.check_reaction(
                        "Ullmann-Goldberg Substitution amine", rsmi
                    )
                    or checker.check_reaction("aza-Michael addition aromatic", rsmi)
                    or checker.check_reaction("aza-Michael addition secondary", rsmi)
                    or checker.check_reaction("aza-Michael addition primary", rsmi)
                    or checker.check_reaction("Aminolysis of esters", rsmi)
                    or checker.check_reaction("Ugi reaction", rsmi)
                    or checker.check_reaction(
                        "Petasis reaction with amines and boronic esters", rsmi
                    )
                    or checker.check_reaction(
                        "Petasis reaction with amines and boronic acids", rsmi
                    )
                    or checker.check_reaction(
                        "Petasis reaction with amines aldehydes and boronic acids", rsmi
                    )
                    or checker.check_reaction("Chan-Lam amine", rsmi)
                    or checker.check_reaction(
                        "Reduction of nitro groups to amines", rsmi
                    )
                    or checker.check_reaction("Phthalimide deprotection", rsmi)
                    or checker.check_reaction(
                        "Intramolecular amination (heterocycle formation)", rsmi
                    )
                    or checker.check_reaction(
                        "Intramolecular amination of azidobiphenyls (heterocycle formation)",
                        rsmi,
                    )
                )

                # Also check for reactions where a nitrogen-containing group is added to a molecule
                if not is_cn_forming_reaction:
                    # Check if product has a new amine or amide that wasn't in the reactants
                    product_mol = Chem.MolFromSmiles(product_part)
                    has_new_n_group = False

                    # Check if the product contains an amine or amide
                    if product_mol and (
                        checker.check_fg("Primary amine", product_part)
                        or checker.check_fg("Secondary amine", product_part)
                        or checker.check_fg("Tertiary amine", product_part)
                        or checker.check_fg("Primary amide", product_part)
                        or checker.check_fg("Secondary amide", product_part)
                        or checker.check_fg("Tertiary amide", product_part)
                        or checker.check_fg("Aniline", product_part)
                    ):
                        # Check if any reactant has a nitro group that's reduced
                        for reactant in reactants:
                            if checker.check_fg("Nitro group", reactant):
                                is_cn_forming_reaction = True
                                print(
                                    f"Detected nitro reduction to amine at depth {depth}"
                                )
                                break

                if is_cn_forming_reaction:
                    # Verify this is a significant fragment coupling by checking reactant sizes
                    reactant_mols = [
                        Chem.MolFromSmiles(r) for r in reactants if r.strip()
                    ]
                    significant_fragments = 0

                    for mol in reactant_mols:
                        if (
                            mol and mol.GetNumHeavyAtoms() >= 5
                        ):  # Consider fragments with 5+ heavy atoms as significant
                            significant_fragments += 1

                    # Consider it a key fragment coupling if at least one significant fragment is involved
                    if significant_fragments >= 1:
                        cn_bond_formations += 1
                        significant_cn_formations.append((rsmi, depth))
                        print(
                            f"Key C-N bond formation detected at depth {depth}: {rsmi}"
                        )
                    else:
                        print(
                            f"C-N bond formation detected but not between key fragments at depth {depth}"
                        )
            except Exception as e:
                print(f"Error analyzing reaction: {e}")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    print(f"Total key C-N bond formations found: {cn_bond_formations}")
    for rxn, depth in significant_cn_formations:
        print(f"  - At depth {depth}: {rxn}")

    # Return True if at least 2 C-N bond formations are found
    return cn_bond_formations >= 2

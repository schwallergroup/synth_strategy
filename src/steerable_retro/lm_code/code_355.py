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
    This function detects if the route includes multiple C-N bond formations.
    """
    c_n_bond_formation_count = 0

    def dfs_traverse(node):
        nonlocal c_n_bond_formation_count

        if node["type"] == "reaction":
            # Extract reaction SMILES
            try:
                rsmi = node["metadata"]["rsmi"]
                print(f"Analyzing reaction: {rsmi}")

                # Extract reactants and product
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check for various C-N bond formation reactions
                c_n_bond_formation = False

                # N-arylation reactions (Buchwald-Hartwig, Ullmann-Goldberg)
                if (
                    checker.check_reaction(
                        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine", rsmi
                    )
                    or checker.check_reaction(
                        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine", rsmi
                    )
                    or checker.check_reaction(
                        "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)", rsmi
                    )
                    or checker.check_reaction("N-arylation_heterocycles", rsmi)
                    or checker.check_reaction("Buchwald-Hartwig", rsmi)
                ):
                    c_n_bond_formation = True
                    print("Detected N-arylation reaction")

                # Reductive amination reactions
                elif (
                    checker.check_reaction("Reductive amination with aldehyde", rsmi)
                    or checker.check_reaction("Reductive amination with ketone", rsmi)
                    or checker.check_reaction("Reductive amination with alcohol", rsmi)
                    or checker.check_reaction("reductive amination", rsmi)
                ):
                    c_n_bond_formation = True
                    print("Detected reductive amination reaction")

                # Amide formation reactions
                elif (
                    checker.check_reaction(
                        "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                        rsmi,
                    )
                    or checker.check_reaction(
                        "Acylation of Nitrogen Nucleophiles by Carboxylic Acids", rsmi
                    )
                    or checker.check_reaction("Acyl chloride with ammonia to amide", rsmi)
                    or checker.check_reaction(
                        "Acyl chloride with primary amine to amide (Schotten-Baumann)", rsmi
                    )
                    or checker.check_reaction("Acyl chloride with secondary amine to amide", rsmi)
                    or checker.check_reaction("Carboxylic acid with primary amine to amide", rsmi)
                    or checker.check_reaction("Ester with ammonia to amide", rsmi)
                    or checker.check_reaction("Ester with primary amine to amide", rsmi)
                    or checker.check_reaction("Ester with secondary amine to amide", rsmi)
                    or checker.check_reaction("Schotten-Baumann_amide", rsmi)
                    or checker.check_reaction("Acylation of primary amines", rsmi)
                    or checker.check_reaction("Acylation of secondary amines", rsmi)
                ):
                    c_n_bond_formation = True
                    print("Detected amide formation reaction")

                # Urea formation
                elif (
                    checker.check_reaction("Urea synthesis via isocyanate and primary amine", rsmi)
                    or checker.check_reaction(
                        "Urea synthesis via isocyanate and secondary amine", rsmi
                    )
                    or checker.check_reaction("urea", rsmi)
                ):
                    c_n_bond_formation = True
                    print("Detected urea formation reaction")

                # Amine alkylation reactions
                elif (
                    checker.check_reaction(
                        "N-alkylation of primary amines with alkyl halides", rsmi
                    )
                    or checker.check_reaction(
                        "N-alkylation of secondary amines with alkyl halides", rsmi
                    )
                    or checker.check_reaction("Alkylation of amines", rsmi)
                    or checker.check_reaction("Methylation with MeI_primary", rsmi)
                    or checker.check_reaction("Methylation with MeI_secondary", rsmi)
                    or checker.check_reaction("Methylation with MeI_tertiary", rsmi)
                    or checker.check_reaction("N-methylation", rsmi)
                    or checker.check_reaction("Eschweiler-Clarke Primary Amine Methylation", rsmi)
                    or checker.check_reaction("Eschweiler-Clarke Secondary Amine Methylation", rsmi)
                    or checker.check_reaction(
                        "Reductive methylation of primary amine with formaldehyde", rsmi
                    )
                ):
                    c_n_bond_formation = True
                    print("Detected amine alkylation reaction")

                # aza-Michael addition
                elif (
                    checker.check_reaction("aza-Michael addition aromatic", rsmi)
                    or checker.check_reaction("aza-Michael addition secondary", rsmi)
                    or checker.check_reaction("aza-Michael addition primary", rsmi)
                ):
                    c_n_bond_formation = True
                    print("Detected aza-Michael addition")

                # Ring formation with C-N bonds
                elif (
                    checker.check_reaction("Formation of NOS Heterocycles", rsmi)
                    or checker.check_reaction(
                        "Intramolecular amination (heterocycle formation)", rsmi
                    )
                    or checker.check_reaction(
                        "Intramolecular amination of azidobiphenyls (heterocycle formation)", rsmi
                    )
                ):
                    c_n_bond_formation = True
                    print("Detected heterocycle formation with C-N bond")

                # If no specific reaction type was detected, check for general C-N bond formation
                if not c_n_bond_formation:
                    # Check for reactions that involve nitrogen-containing compounds
                    has_nitrogen_reactant = False
                    for reactant in reactants_smiles:
                        if "N" in reactant:
                            has_nitrogen_reactant = True
                            break

                    # If product contains nitrogen, check for potential C-N bond formation
                    if has_nitrogen_reactant and "N" in product_smiles:
                        # Check for common N-alkylation patterns
                        if "[CH3:1][N:" in product_smiles or "C[N" in product_smiles:
                            c_n_bond_formation = True
                            print("Detected N-alkylation (general pattern)")

                        # Check for N-arylation patterns
                        elif "[c:" in product_smiles and "[N:" in product_smiles:
                            for reactant in reactants_smiles:
                                if "[NH" in reactant or "NH2" in reactant:
                                    c_n_bond_formation = True
                                    print("Detected N-arylation (general pattern)")
                                    break

                if c_n_bond_formation:
                    c_n_bond_formation_count += 1
                    print(f"C-N bond formation count: {c_n_bond_formation_count}")

            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # The strategy is present if there are at least 2 C-N bond formations
    print(f"Total C-N bond formations found: {c_n_bond_formation_count}")
    return c_n_bond_formation_count >= 2

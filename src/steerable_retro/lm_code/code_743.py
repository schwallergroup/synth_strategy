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
    Detects multiple C-N bond formations including tertiary amine via reductive amination.
    """
    cn_bond_count = 0
    reductive_amination = False

    def dfs_traverse(node, depth=0):
        nonlocal cn_bond_count, reductive_amination

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                try:
                    # Check for reductive amination using the checker function
                    if (
                        checker.check_reaction("Reductive amination with aldehyde", rsmi)
                        or checker.check_reaction("Reductive amination with ketone", rsmi)
                        or checker.check_reaction("Reductive amination with alcohol", rsmi)
                        or checker.check_reaction("Mignonac reaction", rsmi)
                    ):
                        print(f"Reductive amination detected at depth {depth}")
                        reductive_amination = True
                        cn_bond_count += 1  # Count reductive amination as a C-N bond formation

                    # Check for C-N bond formation reactions
                    elif (
                        checker.check_reaction(
                            "N-alkylation of primary amines with alkyl halides", rsmi
                        )
                        or checker.check_reaction(
                            "N-alkylation of secondary amines with alkyl halides", rsmi
                        )
                        or checker.check_reaction(
                            "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine", rsmi
                        )
                        or checker.check_reaction(
                            "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine", rsmi
                        )
                        or checker.check_reaction(
                            "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)", rsmi
                        )
                        or checker.check_reaction(
                            "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                            rsmi,
                        )
                        or checker.check_reaction(
                            "Acylation of Nitrogen Nucleophiles by Carboxylic Acids", rsmi
                        )
                        or checker.check_reaction("Acylation of primary amines", rsmi)
                        or checker.check_reaction("Acylation of secondary amines", rsmi)
                        or checker.check_reaction(
                            "Acyl chloride with primary amine to amide (Schotten-Baumann)", rsmi
                        )
                        or checker.check_reaction(
                            "Acyl chloride with secondary amine to amide", rsmi
                        )
                        or checker.check_reaction(
                            "Carboxylic acid with primary amine to amide", rsmi
                        )
                        or checker.check_reaction("Ester with primary amine to amide", rsmi)
                        or checker.check_reaction("Ester with secondary amine to amide", rsmi)
                        or checker.check_reaction("Aminolysis of esters", rsmi)
                        or checker.check_reaction(
                            "Urea synthesis via isocyanate and primary amine", rsmi
                        )
                        or checker.check_reaction(
                            "Urea synthesis via isocyanate and secondary amine", rsmi
                        )
                        or checker.check_reaction(
                            "Sulfonamide synthesis (Schotten-Baumann) primary amine", rsmi
                        )
                        or checker.check_reaction(
                            "Sulfonamide synthesis (Schotten-Baumann) secondary amine", rsmi
                        )
                        or checker.check_reaction("Ring opening of epoxide with amine", rsmi)
                        or checker.check_reaction("Schotten-Baumann to ester", rsmi)
                        or checker.check_reaction("Ugi reaction", rsmi)
                        or checker.check_reaction("aza-Michael addition aromatic", rsmi)
                        or checker.check_reaction("aza-Michael addition secondary", rsmi)
                        or checker.check_reaction("aza-Michael addition primary", rsmi)
                        or checker.check_reaction(
                            "Intramolecular amination (heterocycle formation)", rsmi
                        )
                        or checker.check_reaction(
                            "Intramolecular amination of azidobiphenyls (heterocycle formation)",
                            rsmi,
                        )
                    ):
                        print(f"C-N bond formation detected at depth {depth}")
                        cn_bond_count += 1

                    # Always run fallback detection to catch any missed C-N bond formations
                    reactant_mols = [Chem.MolFromSmiles(r) for r in reactants]
                    product_mol = Chem.MolFromSmiles(product)

                    if all(reactant_mols) and product_mol:
                        # Check for C-N bond formation using more comprehensive SMARTS pattern
                        cn_pattern = Chem.MolFromSmarts("[#6]~[#7]")
                        product_cn_count = len(product_mol.GetSubstructMatches(cn_pattern))
                        reactant_cn_count = sum(
                            len(mol.GetSubstructMatches(cn_pattern)) for mol in reactant_mols
                        )

                        if product_cn_count > reactant_cn_count:
                            print(f"C-N bond formation detected via SMARTS at depth {depth}")
                            cn_bond_count += 1

                        # Check for reductive amination pattern if not already detected
                        if not reductive_amination:
                            # Check if reactants contain carbonyl and amine
                            has_carbonyl = any(
                                checker.check_fg("Aldehyde", r) or checker.check_fg("Ketone", r)
                                for r in reactants
                            )
                            has_amine = any(
                                checker.check_fg("Primary amine", r)
                                or checker.check_fg("Secondary amine", r)
                                for r in reactants
                            )

                            # Check if product has tertiary or secondary amine
                            has_tertiary_amine = checker.check_fg("Tertiary amine", product)
                            has_secondary_amine = checker.check_fg("Secondary amine", product)

                            if (
                                has_carbonyl
                                and has_amine
                                and (has_tertiary_amine or has_secondary_amine)
                            ):
                                print(
                                    f"Reductive amination detected via functional groups at depth {depth}"
                                )
                                reductive_amination = True
                                cn_bond_count += (
                                    1  # Count reductive amination as a C-N bond formation
                                )

                except Exception as e:
                    print(f"Error in C-N bond detection: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    print(
        f"Final counts: C-N bond formations = {cn_bond_count}, Reductive amination = {reductive_amination}"
    )
    return cn_bond_count >= 2 and reductive_amination

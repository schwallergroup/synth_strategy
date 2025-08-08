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
    This function detects a synthetic strategy involving early halogenation steps
    (bromination and/or chlorination at depth >= 2) or the use of halogenated
    starting materials at early stages.
    """
    halogenation_count = 0
    halogenated_reactants_used = False

    def dfs_traverse(node, depth=0):
        nonlocal halogenation_count, halogenated_reactants_used

        if node["type"] == "reaction" and depth >= 2:  # Early-stage reaction
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                print(f"Analyzing reaction at depth {depth}: {rsmi}")

                # Check if this is a halogenation reaction using checker functions
                is_bromination = (
                    checker.check_reaction("Aromatic bromination", rsmi)
                    or checker.check_reaction("Bromination", rsmi)
                    or checker.check_reaction("Wohl-Ziegler bromination benzyl primary", rsmi)
                    or checker.check_reaction("Wohl-Ziegler bromination benzyl secondary", rsmi)
                    or checker.check_reaction("Wohl-Ziegler bromination benzyl tertiary", rsmi)
                    or checker.check_reaction("Wohl-Ziegler bromination allyl primary", rsmi)
                    or checker.check_reaction("Wohl-Ziegler bromination allyl secondary", rsmi)
                    or checker.check_reaction("Wohl-Ziegler bromination allyl tertiary", rsmi)
                )

                is_chlorination = (
                    checker.check_reaction("Aromatic chlorination", rsmi)
                    or checker.check_reaction("Chlorination", rsmi)
                    or checker.check_reaction("Alcohol to chloride_sulfonyl chloride", rsmi)
                    or checker.check_reaction("Alcohol to chloride_SOCl2", rsmi)
                    or checker.check_reaction("Alcohol to chloride_CHCl3", rsmi)
                    or checker.check_reaction("Alcohol to chloride_CH2Cl2", rsmi)
                    or checker.check_reaction("Alcohol to chloride_PCl5_ortho", rsmi)
                    or checker.check_reaction("Alcohol to chloride_POCl3_ortho", rsmi)
                    or checker.check_reaction("Alcohol to chloride_POCl3_para", rsmi)
                    or checker.check_reaction("Alcohol to chloride_POCl3", rsmi)
                    or checker.check_reaction("Alcohol to chloride_HCl", rsmi)
                    or checker.check_reaction("Alcohol to chloride_Salt", rsmi)
                    or checker.check_reaction("Alcohol to chloride_Other", rsmi)
                    or checker.check_reaction("Primary amine to chloride", rsmi)
                )

                if is_bromination or is_chlorination:
                    print(f"Direct halogenation reaction detected at depth {depth}")
                    halogenation_count += 1
                else:
                    # Check if reactants contain halogen atoms (Br or Cl)
                    try:
                        for reactant in reactants:
                            reactant_mol = Chem.MolFromSmiles(reactant)
                            if reactant_mol:
                                br_atoms = len(
                                    [a for a in reactant_mol.GetAtoms() if a.GetSymbol() == "Br"]
                                )
                                cl_atoms = len(
                                    [a for a in reactant_mol.GetAtoms() if a.GetSymbol() == "Cl"]
                                )

                                if br_atoms > 0 or cl_atoms > 0:
                                    print(
                                        f"Halogenated reactant found at depth {depth}: {reactant}"
                                    )
                                    print(f"Br atoms: {br_atoms}, Cl atoms: {cl_atoms}")
                                    halogenated_reactants_used = True
                                    break

                        # Also check for increase in halide atoms (fallback method)
                        product_mol = Chem.MolFromSmiles(product)
                        reactant_mols = [Chem.MolFromSmiles(r) for r in reactants]

                        if product_mol:
                            product_br_count = len(
                                [a for a in product_mol.GetAtoms() if a.GetSymbol() == "Br"]
                            )
                            product_cl_count = len(
                                [a for a in product_mol.GetAtoms() if a.GetSymbol() == "Cl"]
                            )

                            reactants_br_count = sum(
                                len([a for a in mol.GetAtoms() if a.GetSymbol() == "Br"])
                                for mol in reactant_mols
                                if mol is not None
                            )
                            reactants_cl_count = sum(
                                len([a for a in mol.GetAtoms() if a.GetSymbol() == "Cl"])
                                for mol in reactant_mols
                                if mol is not None
                            )

                            new_halide_formed = (product_br_count > reactants_br_count) or (
                                product_cl_count > reactants_cl_count
                            )

                            if new_halide_formed:
                                print(f"Halogenation via atom count detected at depth {depth}")
                                print(
                                    f"Br count: {reactants_br_count} → {product_br_count}, Cl count: {reactants_cl_count} → {product_cl_count}"
                                )
                                halogenation_count += 1
                    except Exception as e:
                        print(f"Error checking for halide atoms: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    print(f"Total halogenation reactions found: {halogenation_count}")
    print(f"Halogenated reactants used in early stages: {halogenated_reactants_used}")

    # Return True if either halogenation reactions are found or halogenated reactants are used
    return bool(halogenation_count >= 1 or halogenated_reactants_used)

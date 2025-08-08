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
    This function detects a strategy involving the formation of a quinoxalinone scaffold
    through amide formation and cyclization.
    """
    # Track if we found the quinoxalinone formation strategy
    found_quinoxalinone = False

    # Track reactions in the pathway
    reaction_sequence = []

    def dfs_traverse(node, depth=0):
        nonlocal found_quinoxalinone

        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"].get("rsmi", "")
                if not rsmi:
                    return

                reactants_part = rsmi.split(">")[0]
                product_smiles = rsmi.split(">")[-1]
                reactants = reactants_part.split(".")

                product_mol = Chem.MolFromSmiles(product_smiles)
                if not product_mol:
                    return

                # Check if product contains quinoxalinone structure
                # Quinoxalinone has a quinoxaline ring with an amide group
                if checker.check_ring("quinoxaline", product_smiles):
                    print(f"Found quinoxaline ring in product at depth {depth}: {product_smiles}")

                    # Check for amide group in the product
                    if checker.check_fg("Tertiary amide", product_smiles) or checker.check_fg(
                        "Secondary amide", product_smiles
                    ):
                        print(f"Found amide group in quinoxaline product at depth {depth}")

                        # Check if any reactant has diamine structure (o-phenylenediamine)
                        has_diamine = False
                        for reactant in reactants:
                            if (
                                checker.check_fg("Primary amine", reactant)
                                and reactant.count("N") >= 2
                            ):
                                # Check if the reactant has an aromatic ring
                                reactant_mol = Chem.MolFromSmiles(reactant)
                                if reactant_mol and reactant_mol.GetNumAtoms() >= 6:
                                    if (
                                        sum(
                                            1
                                            for atom in reactant_mol.GetAtoms()
                                            if atom.GetIsAromatic()
                                        )
                                        >= 6
                                    ):
                                        has_diamine = True
                                        print(
                                            f"Found potential o-phenylenediamine derivative in reactants: {reactant}"
                                        )
                                        break

                        # Check if this is an amide formation reaction
                        is_amide_formation = (
                            checker.check_reaction(
                                "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                                rsmi,
                            )
                            or checker.check_reaction("Acylation of primary amines", rsmi)
                            or checker.check_reaction("Acylation of secondary amines", rsmi)
                            or checker.check_reaction("Schotten-Baumann_amide", rsmi)
                        )

                        # Check if this is a cyclization reaction
                        is_cyclization = checker.check_reaction(
                            "Formation of NOS Heterocycles", rsmi
                        )

                        # If we found evidence of quinoxalinone formation
                        if (is_amide_formation or is_cyclization) and has_diamine:
                            print(f"Found quinoxalinone formation strategy at depth {depth}")
                            found_quinoxalinone = True

                            # Add to reaction sequence
                            reaction_sequence.append((depth, "quinoxalinone_formation", rsmi))

            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)

    # If we didn't find a clear quinoxalinone formation, check for the structure itself
    if not found_quinoxalinone:

        def check_for_quinoxalinone(node):
            nonlocal found_quinoxalinone
            if node["type"] == "mol":
                mol_smiles = node["smiles"]
                mol = Chem.MolFromSmiles(mol_smiles)
                if mol:
                    # Check for quinoxalinone structure
                    quinoxalinone_pattern = Chem.MolFromSmarts("c1nc2ccccc2nc1C(=O)N")
                    if mol.HasSubstructMatch(quinoxalinone_pattern) or (
                        checker.check_ring("quinoxaline", mol_smiles)
                        and (
                            checker.check_fg("Tertiary amide", mol_smiles)
                            or checker.check_fg("Secondary amide", mol_smiles)
                        )
                    ):
                        print(f"Found quinoxalinone structure in molecule: {mol_smiles}")
                        found_quinoxalinone = True

            # Check children
            for child in node.get("children", []):
                if not found_quinoxalinone:
                    check_for_quinoxalinone(child)

        check_for_quinoxalinone(route)

    return found_quinoxalinone

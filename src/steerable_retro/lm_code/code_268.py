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
    This function detects the use of long alkyl chains to connect aromatic rings,
    a strategy often used in materials chemistry.
    """
    long_chain_connection = False

    # List of aromatic rings to check
    aromatic_rings = [
        "benzene",
        "naphthalene",
        "anthracene",
        "furan",
        "pyrrole",
        "thiophene",
        "pyridine",
        "pyrazole",
        "imidazole",
        "oxazole",
        "thiazole",
        "pyrimidine",
        "pyrazine",
        "pyridazine",
        "triazole",
        "tetrazole",
        "indole",
        "quinoline",
        "isoquinoline",
        "carbazole",
        "acridine",
        "benzothiophene",
        "isoxazole",
        "isothiazole",
        "oxadiazole",
        "thiadiazole",
        "benzoxazole",
        "benzothiazole",
        "benzimidazole",
        "indazole",
        "benzotriazole",
    ]

    def dfs_traverse(node, depth=0):
        nonlocal long_chain_connection

        if node["type"] == "mol":
            mol_smiles = node["smiles"]
            mol = Chem.MolFromSmiles(mol_smiles)

            if mol:
                # Check for aromatic rings using checker function
                has_aromatic = False
                for ring_name in aromatic_rings:
                    if checker.check_ring(ring_name, mol_smiles):
                        has_aromatic = True
                        print(f"Found aromatic ring {ring_name} in molecule: {mol_smiles}")
                        break

                if has_aromatic:
                    # Look for long alkyl chains (4+ consecutive atoms) connecting aromatic rings
                    # More flexible patterns that can match substituted chains and heteroatoms
                    pattern4 = Chem.MolFromSmarts("a~[#6]~[#6]~[#6]~[#6]~a")  # 4-atom chain
                    pattern5 = Chem.MolFromSmarts("a~[#6]~[#6]~[#6]~[#6]~[#6]~a")  # 5-atom chain
                    pattern6 = Chem.MolFromSmarts(
                        "a~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~a"
                    )  # 6-atom chain
                    pattern_hetero = Chem.MolFromSmarts(
                        "a~[#6]~[#6]~[*]~[#6]~[#6]~a"
                    )  # Chain with heteroatom

                    if (
                        mol.HasSubstructMatch(pattern4)
                        or mol.HasSubstructMatch(pattern5)
                        or mol.HasSubstructMatch(pattern6)
                        or mol.HasSubstructMatch(pattern_hetero)
                    ):
                        long_chain_connection = True
                        print(
                            f"Found long chain connecting aromatic rings in molecule: {mol_smiles}"
                        )

        elif node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants_part = rsmi.split(">")[0]
                reactants = reactants_part.split(".")
                product = rsmi.split(">")[-1]

                print(f"Analyzing reaction: {rsmi}")

                # Check if the reaction creates a long chain connection
                product_mol = Chem.MolFromSmiles(product)
                if product_mol:
                    # Check if product has aromatic rings
                    product_has_aromatic = False
                    for ring_name in aromatic_rings:
                        if checker.check_ring(ring_name, product):
                            product_has_aromatic = True
                            print(f"Product has aromatic ring {ring_name}")
                            break

                    if product_has_aromatic:
                        # Check for long chains in product
                        pattern4 = Chem.MolFromSmarts("a~[#6]~[#6]~[#6]~[#6]~a")
                        pattern5 = Chem.MolFromSmarts("a~[#6]~[#6]~[#6]~[#6]~[#6]~a")
                        pattern6 = Chem.MolFromSmarts("a~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~a")
                        pattern_hetero = Chem.MolFromSmarts("a~[#6]~[#6]~[*]~[#6]~[#6]~a")

                        product_has_chain = (
                            product_mol.HasSubstructMatch(pattern4)
                            or product_mol.HasSubstructMatch(pattern5)
                            or product_mol.HasSubstructMatch(pattern6)
                            or product_mol.HasSubstructMatch(pattern_hetero)
                        )

                        if product_has_chain:
                            print(f"Product has long chain between aromatic rings")

                            # Check if any reactant already has this pattern
                            reactant_has_chain = False
                            for reactant in reactants:
                                r_mol = Chem.MolFromSmiles(reactant)
                                if r_mol:
                                    # Check if reactant has aromatic rings
                                    reactant_has_aromatic = False
                                    for ring_name in aromatic_rings:
                                        if checker.check_ring(ring_name, reactant):
                                            reactant_has_aromatic = True
                                            break

                                    if reactant_has_aromatic:
                                        # Check for chains in reactant
                                        if (
                                            r_mol.HasSubstructMatch(pattern4)
                                            or r_mol.HasSubstructMatch(pattern5)
                                            or r_mol.HasSubstructMatch(pattern6)
                                            or r_mol.HasSubstructMatch(pattern_hetero)
                                        ):
                                            reactant_has_chain = True
                                            print(f"Reactant already has long chain: {reactant}")
                                            break

                            # If product has the pattern but reactants don't, the reaction created it
                            if not reactant_has_chain:
                                long_chain_connection = True
                                print(f"Found reaction creating long chain connection: {rsmi}")

                                # Check for specific reactions that might create these connections
                                if checker.check_reaction("Williamson Ether Synthesis", rsmi):
                                    print(
                                        "Detected Williamson Ether Synthesis creating chain connection"
                                    )
                                elif checker.check_reaction(
                                    "Suzuki coupling with boronic acids", rsmi
                                ):
                                    print("Detected Suzuki coupling creating chain connection")
                                elif checker.check_reaction("Heck terminal vinyl", rsmi):
                                    print("Detected Heck reaction creating chain connection")
                                elif checker.check_reaction("Sonogashira alkyne_aryl halide", rsmi):
                                    print("Detected Sonogashira coupling creating chain connection")
            except Exception as e:
                print(f"Error processing reaction SMILES: {e}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    print(f"Long chain connection strategy found: {long_chain_connection}")

    return long_chain_connection

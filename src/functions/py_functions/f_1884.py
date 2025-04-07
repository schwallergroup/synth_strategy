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
    Detects if the route involves formation of nitrogen-containing heterocycles
    """
    heterocycle_formed = False

    # List of common nitrogen-containing heterocycles to check
    n_heterocycles = [
        "pyrrole",
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
        "pyrrolidine",
        "piperidine",
        "piperazine",
        "morpholine",
        "aziridine",
        "azetidine",
        "azepane",
        "indole",
        "quinoline",
        "isoquinoline",
        "benzimidazole",
        "benzoxazole",
        "benzothiazole",
        "indazole",
        "benzotriazole",
    ]

    # List of reactions commonly used for heterocycle formation
    heterocycle_forming_reactions = [
        "Formation of NOS Heterocycles",
        "Paal-Knorr pyrrole synthesis",
        "Benzothiazole formation from aldehyde",
        "Benzothiazole formation from acyl halide",
        "Benzothiazole formation from ester/carboxylic acid",
        "Benzoxazole formation from aldehyde",
        "Benzoxazole formation from acyl halide",
        "Benzoxazole formation from ester/carboxylic acid",
        "Benzoxazole formation (intramolecular)",
        "Benzimidazole formation from aldehyde",
        "Benzimidazole formation from acyl halide",
        "Benzimidazole formation from ester/carboxylic acid",
        "Pyrazole formation",
        "Huisgen alkyne-azide 1,3 dipolar cycloaddition",
        "Huisgen 1,3 dipolar cycloaddition",
        "Huisgen alkene-azide 1,3 dipolar cycloaddition",
        "A3 coupling to imidazoles",
        "Alkyne-imine cycloaddition",
        "Azide-nitrile click cycloaddition to tetrazole",
        "Azide-nitrile click cycloaddition to triazole",
        "Intramolecular amination of azidobiphenyls (heterocycle formation)",
        "Intramolecular amination (heterocycle formation)",
    ]

    def dfs_traverse(node):
        nonlocal heterocycle_formed

        if heterocycle_formed:
            return  # Early return if we already found a heterocycle formation

        if node["type"] == "reaction":
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]

                # Check if this is a known heterocycle-forming reaction
                for reaction_type in heterocycle_forming_reactions:
                    if checker.check_reaction(reaction_type, rsmi):
                        print(f"Found heterocycle formation reaction: {reaction_type}")
                        heterocycle_formed = True
                        return

                # If not a known reaction type, check for ring formation
                reactants_part = rsmi.split(">")[0]
                product_part = rsmi.split(">")[-1]

                # Skip if we can't parse the molecules
                if not Chem.MolFromSmiles(product_part):
                    for child in node.get("children", []):
                        dfs_traverse(child)
                    return

                # Check if any nitrogen-containing heterocycle is in the product but not in reactants
                product_has_n_heterocycle = False
                reactants_have_same_n_heterocycle = False

                for ring_name in n_heterocycles:
                    if checker.check_ring(ring_name, product_part):
                        product_has_n_heterocycle = True

                        # Check if any reactant has the same heterocycle
                        for reactant in reactants_part.split("."):
                            if reactant and checker.check_ring(ring_name, reactant):
                                reactants_have_same_n_heterocycle = True
                                break

                        if not reactants_have_same_n_heterocycle:
                            print(
                                f"Found heterocycle formation: {ring_name} formed in product"
                            )
                            heterocycle_formed = True
                            return

                # If we haven't found a heterocycle yet, check ring count difference
                try:
                    reactants_mol = [
                        Chem.MolFromSmiles(r)
                        for r in reactants_part.split(".")
                        if r and Chem.MolFromSmiles(r)
                    ]
                    product_mol = Chem.MolFromSmiles(product_part)

                    if product_mol and reactants_mol:
                        # Count rings in reactants and product
                        reactant_rings = sum(
                            [r.GetRingInfo().NumRings() for r in reactants_mol]
                        )
                        product_rings = product_mol.GetRingInfo().NumRings()

                        # Check if product has more rings and contains nitrogen in rings
                        if product_rings > reactant_rings:
                            # Check if any nitrogen atom is in a ring in the product
                            for atom in product_mol.GetAtoms():
                                if atom.GetAtomicNum() == 7 and atom.IsInRing():
                                    print(
                                        f"Found heterocycle formation: New ring with nitrogen detected"
                                    )
                                    heterocycle_formed = True
                                    return
                except Exception as e:
                    print(f"Error analyzing rings: {e}")

        # Process children nodes
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return heterocycle_formed

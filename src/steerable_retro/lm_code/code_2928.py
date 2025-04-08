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
    This function detects a synthetic strategy involving a late-stage cyclization
    where a ring is formed in the final steps of the synthesis.
    """
    cyclization_found = False
    cyclization_depth = float("inf")

    # List of common cyclization reaction types
    cyclization_reaction_types = [
        "Intramolecular amination (heterocycle formation)",
        "Intramolecular amination of azidobiphenyls (heterocycle formation)",
        "Formation of NOS Heterocycles",
        "Paal-Knorr pyrrole synthesis",
        "Diels-Alder",
        "Huisgen alkyne-azide 1,3 dipolar cycloaddition",
        "Huisgen 1,3 dipolar cycloaddition",
        "Huisgen alkene-azide 1,3 dipolar cycloaddition",
        "Pyrazole formation",
        "Pauson-Khand reaction",
        "Michael-induced ring closure from hydrazone",
        "Michael-induced ring closure from diazoalkane",
        "[3+2]-cycloaddition of hydrazone and alkyne",
        "[3+2]-cycloaddition of hydrazone and alkene",
        "[3+2]-cycloaddition of diazoalkane and alkyne",
        "[3+2]-cycloaddition of diazoalkane and alkene",
        "[3+2]-cycloaddition of diazoalkane and alpha-alkyne",
        "[3+2]-cycloaddition of diazoalkane and alpha-alkene",
    ]

    # List of common ring types to check for formation
    ring_types = [
        "furan",
        "pyran",
        "dioxane",
        "tetrahydrofuran",
        "tetrahydropyran",
        "oxirane",
        "oxetane",
        "oxolane",
        "oxane",
        "dioxolane",
        "dioxolene",
        "trioxane",
        "dioxepane",
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
        "thiomorpholine",
        "aziridine",
        "azetidine",
        "azepane",
        "diazepane",
        "indole",
        "quinoline",
        "isoquinoline",
        "purine",
        "carbazole",
        "acridine",
        "thiophene",
        "thiopyran",
        "thiirane",
        "thietane",
        "thiolane",
        "thiane",
        "dithiane",
        "dithiolane",
        "benzothiophene",
        "oxathiolane",
        "dioxathiolane",
        "thiazolidine",
        "oxazolidine",
        "isoxazole",
        "isothiazole",
        "oxadiazole",
        "thiadiazole",
        "cyclopropane",
        "cyclobutane",
        "cyclopentane",
        "cyclohexane",
        "cycloheptane",
        "cyclooctane",
        "benzene",
        "naphthalene",
        "anthracene",
        "benzoxazole",
        "benzothiazole",
        "benzimidazole",
        "pteridin",
        "phenothiazine",
        "phenoxazine",
        "dibenzofuran",
        "dibenzothiophene",
        "xanthene",
        "thioxanthene",
        "pyrroline",
        "pyrrolidone",
        "imidazolidine",
        "porphyrin",
        "indazole",
        "benzotriazole",
    ]

    def dfs_traverse(node, depth=0):
        nonlocal cyclization_found, cyclization_depth

        if node["type"] == "reaction":
            try:
                # Extract reactants and product
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check if this is a known cyclization reaction type
                is_cyclization_reaction = False
                for reaction_type in cyclization_reaction_types:
                    if checker.check_reaction(reaction_type, rsmi):
                        print(f"Cyclization reaction detected: {reaction_type} at depth {depth}")
                        is_cyclization_reaction = True
                        break

                # Convert to RDKit molecules
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles]
                product_mol = Chem.MolFromSmiles(product_smiles)

                if all(reactant_mols) and product_mol:
                    # Count rings in reactants and product
                    reactant_ring_counts = [len(Chem.GetSSSR(r)) for r in reactant_mols]
                    product_ring_count = len(Chem.GetSSSR(product_mol))

                    # Check if this reaction forms a ring
                    if product_ring_count > max(reactant_ring_counts):
                        print(
                            f"Ring count increase detected at depth {depth}: {max(reactant_ring_counts)} -> {product_ring_count}"
                        )
                        is_cyclization_reaction = True  # Set to true when ring count increases

                        # Check if specific ring types are formed in the product but not in reactants
                        for ring_type in ring_types:
                            product_has_ring = checker.check_ring(ring_type, product_smiles)
                            reactants_have_ring = any(
                                checker.check_ring(ring_type, r) for r in reactants_smiles
                            )

                            if product_has_ring and not reactants_have_ring:
                                print(f"New ring type formed: {ring_type} at depth {depth}")
                                break

                    if is_cyclization_reaction:
                        cyclization_found = True
                        cyclization_depth = min(cyclization_depth, depth)
            except KeyError as e:
                print(f"Missing key in reaction metadata at depth {depth}: {e}")
            except Exception as e:
                print(f"Error processing reaction SMILES at depth {depth}: {e}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if cyclization occurred in the late stage (depth <= 3)
    if cyclization_found and cyclization_depth <= 3:
        print(f"Late-stage cyclization detected at depth {cyclization_depth}")
        return True
    return False

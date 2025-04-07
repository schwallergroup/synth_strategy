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
    Detects if the route includes a heterocycle formation step.
    """
    heterocycle_formation_detected = False

    # List of common heterocycles to check
    heterocycle_types = [
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

    # List of heterocycle formation reactions to check
    heterocycle_formation_reactions = [
        "Formation of NOS Heterocycles",
        "Paal-Knorr pyrrole synthesis",
        "benzimidazole_derivatives_carboxylic-acid/ester",
        "benzimidazole_derivatives_aldehyde",
        "benzothiazole",
        "benzoxazole_arom-aldehyde",
        "benzoxazole_carboxylic-acid",
        "thiazole",
        "tetrazole_terminal",
        "tetrazole_connect_regioisomere_1",
        "tetrazole_connect_regioisomere_2",
        "1,2,4-triazole_acetohydrazide",
        "1,2,4-triazole_carboxylic-acid/ester",
        "3-nitrile-pyridine",
        "pyrazole",
        "Paal-Knorr pyrrole",
        "triaryl-imidazole",
        "Fischer indole",
        "Friedlaender chinoline",
        "benzofuran",
        "benzothiophene",
        "indole",
        "oxadiazole",
        "imidazole",
    ]

    def dfs_traverse(node, depth=0):
        nonlocal heterocycle_formation_detected

        if node["type"] == "reaction" and node.get("metadata", {}).get("rsmi"):
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            print(f"Analyzing reaction at depth {depth}: {rsmi}")

            # Check if this is a known heterocycle formation reaction
            for reaction_type in heterocycle_formation_reactions:
                if checker.check_reaction(reaction_type, rsmi):
                    print(f"Detected heterocycle formation reaction: {reaction_type}")
                    heterocycle_formation_detected = True
                    return

            # Count rings in reactants and product
            product_mol = Chem.MolFromSmiles(product_smiles)
            if not product_mol:
                print(f"Could not parse product SMILES: {product_smiles}")
                return

            product_ring_count = product_mol.GetRingInfo().NumRings()

            reactant_mols = []
            reactant_ring_count = 0
            for r_smi in reactants_smiles:
                r_mol = Chem.MolFromSmiles(r_smi)
                if r_mol:
                    reactant_mols.append(r_mol)
                    reactant_ring_count += r_mol.GetRingInfo().NumRings()

            # Check if product has more rings than reactants
            if product_ring_count > reactant_ring_count:
                print(f"Ring formation detected: {reactant_ring_count} → {product_ring_count}")

                # Check if the new ring is a heterocycle
                heterocycle_found = False
                for heterocycle in heterocycle_types:
                    if checker.check_ring(heterocycle, product_smiles):
                        # Verify this heterocycle wasn't in the reactants
                        reactant_has_heterocycle = False
                        for r_mol_smiles in reactants_smiles:
                            if checker.check_ring(heterocycle, r_mol_smiles):
                                reactant_has_heterocycle = True
                                break

                        if not reactant_has_heterocycle:
                            print(f"New heterocycle formed: {heterocycle}")
                            heterocycle_formation_detected = True
                            heterocycle_found = True
                            break

                if not heterocycle_found:
                    print(
                        "New ring formed but it's not a heterocycle or it already existed in reactants"
                    )

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return heterocycle_formation_detected

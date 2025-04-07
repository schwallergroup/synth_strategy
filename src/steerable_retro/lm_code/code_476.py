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
    This function detects if a heterocycle (benzimidazole, benzoxazole, benzothiazole, etc.)
    is formed early in the synthesis (at high depth levels in the retrosynthetic tree).
    """
    heterocycle_formed = False
    heterocycle_formation_depth = -1
    heterocycle_type = None
    max_depth = -1

    # List of heterocycles to check - using all available heterocycles
    heterocycles = [
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

    # List of heterocycle formation reaction types - using all available formation reactions
    formation_reactions = [
        "benzimidazole_derivatives_carboxylic-acid/ester",
        "benzimidazole_derivatives_aldehyde",
        "benzothiazole",
        "benzoxazole_arom-aldehyde",
        "benzoxazole_carboxylic-acid",
        "thiazole",
        "Niementowski_quinazoline",
        "tetrazole_terminal",
        "tetrazole_connect_regioisomere_1",
        "tetrazole_connect_regioisomere_2",
        "Huisgen_Cu-catalyzed_1,4-subst",
        "Huisgen_Ru-catalyzed_1,5_subst",
        "Huisgen_disubst-alkyne",
        "1,2,4-triazole_acetohydrazide",
        "1,2,4-triazole_carboxylic-acid/ester",
        "3-nitrile-pyridine",
        "spiro-chromanone",
        "pyrazole",
        "phthalazinone",
        "Paal-Knorr pyrrole",
        "triaryl-imidazole",
        "Fischer indole",
        "Friedlaender chinoline",
        "benzofuran",
        "benzothiophene",
        "indole",
        "oxadiazole",
        "imidazole",
        "Pictet-Spengler",
        "Huisgen 1,3,4-oxadiazoles from COOH and tetrazole",
        "Huisgen alkyne-azide 1,3 dipolar cycloaddition",
        "Huisgen 1,3 dipolar cycloaddition",
        "Huisgen alkene-azide 1,3 dipolar cycloaddition",
        "Pyrazole formation",
        "Azide-nitrile click cycloaddition to tetrazole",
        "Azide-nitrile click cycloaddition to triazole",
        "Intramolecular amination of azidobiphenyls (heterocycle formation)",
        "Intramolecular amination (heterocycle formation)",
        "Benzimidazole formation from aldehyde",
        "Benzimidazole formation from acyl halide",
        "Benzimidazole formation from ester/carboxylic acid",
        "Benzoxazole formation from aldehyde",
        "Benzoxazole formation from acyl halide",
        "Benzoxazole formation from ester/carboxylic acid",
        "Benzoxazole formation (intramolecular)",
        "Benzothiazole formation from aldehyde",
        "Benzothiazole formation from acyl halide",
        "Benzothiazole formation from ester/carboxylic acid",
        "Paal-Knorr pyrrole synthesis",
        "Formation of NOS Heterocycles",
        "Benzimidazole aldehyde",
    ]

    def dfs_traverse(node, depth=0):
        nonlocal heterocycle_formed, heterocycle_formation_depth, heterocycle_type, max_depth

        max_depth = max(max_depth, depth)

        if node["type"] == "reaction":
            # Extract product and reactants
            rsmi = node["metadata"].get("rsmi", "")
            if not rsmi:
                return

            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # First check if this is a known heterocycle formation reaction
            is_heterocycle_formation_reaction = False
            for reaction_type in formation_reactions:
                if checker.check_reaction(reaction_type, rsmi):
                    is_heterocycle_formation_reaction = True
                    print(
                        f"Heterocycle formation reaction detected: {reaction_type} at depth {depth}"
                    )
                    break

            # If it's a known formation reaction or we need to check all reactions
            if is_heterocycle_formation_reaction or True:  # Always check for heterocycle formation
                # Check if any heterocycle is formed in this reaction
                for heterocycle in heterocycles:
                    # Check if product has the heterocycle
                    if checker.check_ring(heterocycle, product):
                        # Check if reactants don't have the heterocycle
                        has_heterocycle_in_reactants = False
                        for reactant in reactants:
                            if checker.check_ring(heterocycle, reactant):
                                has_heterocycle_in_reactants = True
                                break

                        if not has_heterocycle_in_reactants:
                            heterocycle_formed = True
                            heterocycle_formation_depth = depth
                            heterocycle_type = heterocycle
                            print(f"{heterocycle} formation detected at depth {depth}")
                            break

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Consider it early heterocycle formation if it happens in the first 40% of the synthesis depth
    # In retrosynthesis, early stage corresponds to high depth
    is_early_formation = False
    if heterocycle_formed and max_depth > 0:
        # Early formation means at a depth that's at least 60% of the maximum depth
        # (since in retrosynthesis, higher depth = earlier in the forward synthesis)
        formation_percentage = heterocycle_formation_depth / max_depth
        is_early_formation = formation_percentage >= 0.6
        print(f"Heterocycle formation depth percentage: {formation_percentage:.2f}")

    print(
        f"Early heterocycle formation: {is_early_formation} (formed at depth {heterocycle_formation_depth}, max depth {max_depth}, heterocycle: {heterocycle_type})"
    )
    return is_early_formation

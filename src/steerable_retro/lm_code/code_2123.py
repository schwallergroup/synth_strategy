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
    This function detects a convergent synthesis approach for heterocyclic compounds.
    It looks for reactions that combine multiple fragments containing heterocyclic rings.
    """
    convergent_steps = 0
    heterocycle_formations = 0

    # List of heterocyclic rings to check
    heterocycle_ring_names = [
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

    # List of known heterocycle-forming reactions
    heterocycle_forming_reactions = [
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
        "1,2,4-triazole_acetohydrazide",
        "1,2,4-triazole_carboxylic-acid/ester",
        "pyrazole",
        "Paal-Knorr pyrrole",
        "Fischer indole",
        "oxadiazole",
        "benzofuran",
        "benzothiophene",
        "indole",
        "Friedlaender chinoline",
        "triaryl-imidazole",
        "imidazole",
        "Pictet-Spengler",
    ]

    def dfs_traverse(node):
        nonlocal convergent_steps, heterocycle_formations

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check if this is a convergent step (multiple reactants)
            if len(reactants_smiles) > 1:
                # Check if reactants contain heterocycles
                heterocycle_reactants = 0
                reactant_heterocycles = set()

                for reactant in reactants_smiles:
                    reactant_has_heterocycle = False
                    for ring_name in heterocycle_ring_names:
                        if checker.check_ring(ring_name, reactant):
                            reactant_has_heterocycle = True
                            reactant_heterocycles.add(ring_name)

                    if reactant_has_heterocycle:
                        heterocycle_reactants += 1

                # Check if product contains heterocycles
                product_has_heterocycle = False
                product_heterocycles = set()

                for ring_name in heterocycle_ring_names:
                    if checker.check_ring(ring_name, product_smiles):
                        product_has_heterocycle = True
                        product_heterocycles.add(ring_name)

                # Count unique heterocycles in each reactant
                reactant_heterocycle_counts = []
                for reactant in reactants_smiles:
                    count = sum(
                        1
                        for ring_name in heterocycle_ring_names
                        if checker.check_ring(ring_name, reactant)
                    )
                    reactant_heterocycle_counts.append(count)

                # Count heterocycles in product
                product_heterocycle_count = sum(
                    1
                    for ring_name in heterocycle_ring_names
                    if checker.check_ring(ring_name, product_smiles)
                )

                # Check for convergent heterocycle synthesis
                if heterocycle_reactants >= 2 and product_has_heterocycle:
                    convergent_steps += 1
                    print(f"Found convergent heterocycle step: {rsmi}")

                    # Also count this as a heterocycle formation since it's combining heterocycles
                    heterocycle_formations += 1
                    print(f"Found heterocycle formation (combining heterocycles): {rsmi}")

                # Check for new heterocycle formation
                new_heterocycles = product_heterocycles - reactant_heterocycles
                if new_heterocycles:
                    heterocycle_formations += 1
                    print(
                        f"Found heterocycle formation (new type): {rsmi}, new heterocycles: {new_heterocycles}"
                    )

                # Check for heterocycle coupling (C-N, C-O, C-S bonds between heterocycles)
                if heterocycle_reactants >= 2 and product_has_heterocycle:
                    # Check if this is a coupling reaction
                    coupling_reactions = [
                        "N-arylation",
                        "Buchwald-Hartwig",
                        "Ullmann-Goldberg",
                        "Suzuki",
                        "Stille",
                        "Negishi",
                        "Heck",
                    ]

                    for rxn_type in coupling_reactions:
                        if checker.check_reaction(rxn_type, rsmi):
                            heterocycle_formations += 1
                            print(f"Found heterocycle coupling via {rxn_type}: {rsmi}")
                            break

            # Check for known heterocycle-forming reactions
            for rxn_type in heterocycle_forming_reactions:
                if checker.check_reaction(rxn_type, rsmi):
                    heterocycle_formations += 1
                    print(f"Found heterocycle-forming reaction {rxn_type}: {rsmi}")
                    break

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Return True if convergent heterocycle synthesis is detected
    result = convergent_steps > 0 and heterocycle_formations > 0
    print(
        f"Convergent heterocycle synthesis detected: {result} (convergent steps: {convergent_steps}, heterocycle formations: {heterocycle_formations})"
    )
    return result

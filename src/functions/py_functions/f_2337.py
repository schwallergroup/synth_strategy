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
    Detects if the synthesis involves the assembly of multiple heterocyclic systems.
    """
    # List of heterocycles to check
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

    # Track heterocycles in the final product and those formed during synthesis
    final_product_heterocycles = set()
    heterocycles_in_starting_materials = set()
    heterocycle_forming_reactions = 0
    heterocycle_assembly_reactions = 0

    # Heterocycle-forming reaction types
    heterocycle_forming_reaction_types = [
        "Paal-Knorr pyrrole synthesis",
        "{benzimidazole_derivatives_carboxylic-acid/ester}",
        "{benzimidazole_derivatives_aldehyde}",
        "{benzothiazole}",
        "{benzoxazole_arom-aldehyde}",
        "{benzoxazole_carboxylic-acid}",
        "{thiazole}",
        "{tetrazole_terminal}",
        "{Huisgen_Cu-catalyzed_1,4-subst}",
        "{pyrazole}",
        "{Fischer indole}",
        "{Friedlaender chinoline}",
        "{benzofuran}",
        "{benzothiophene}",
        "{indole}",
        "{oxadiazole}",
        "{imidazole}",
        "Formation of NOS Heterocycles",
        "Paal-Knorr pyrrole",
        "Huisgen alkyne-azide 1,3 dipolar cycloaddition",
        "Huisgen 1,3 dipolar cycloaddition",
        "Huisgen alkene-azide 1,3 dipolar cycloaddition",
        "Pyrazole formation",
        "A3 coupling to imidazoles",
        "Azide-nitrile click cycloaddition to tetrazole",
        "Azide-nitrile click cycloaddition to triazole",
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
        "Huisgen 1,3,4-oxadiazoles from COOH and tetrazole",
    ]

    def dfs_traverse(node, depth=0, is_starting_material=False):
        nonlocal heterocycle_forming_reactions, heterocycle_assembly_reactions

        if node["type"] == "mol":
            mol_smiles = node["smiles"]

            # Check for heterocycles in the molecule
            current_heterocycles = set()
            for heterocycle in heterocycles:
                if checker.check_ring(heterocycle, mol_smiles):
                    current_heterocycles.add(heterocycle)

                    # If this is the final product (depth=0), add to final product heterocycles
                    if depth == 0:
                        final_product_heterocycles.add(heterocycle)
                        print(f"Heterocycle in final product: {heterocycle}")

                    # If this is a starting material, add to starting materials heterocycles
                    if is_starting_material or node.get("in_stock", False):
                        heterocycles_in_starting_materials.add(heterocycle)
                        print(f"Heterocycle in starting material: {heterocycle}")

        elif node["type"] == "reaction":
            try:
                # Get reaction SMILES
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check which heterocycles are in the product but not in any reactant
                product_heterocycles = set()
                reactant_heterocycles = set()

                # Check heterocycles in product
                for heterocycle in heterocycles:
                    if checker.check_ring(heterocycle, product_smiles):
                        product_heterocycles.add(heterocycle)

                # Check heterocycles in reactants
                for reactant in reactants_smiles:
                    for heterocycle in heterocycles:
                        if checker.check_ring(heterocycle, reactant):
                            reactant_heterocycles.add(heterocycle)

                # Identify newly formed heterocycles
                new_heterocycles = product_heterocycles - reactant_heterocycles

                if new_heterocycles:
                    print(f"Heterocycle formation detected: {new_heterocycles}")

                    # Check for specific heterocycle-forming reactions
                    for rxn_type in heterocycle_forming_reaction_types:
                        if checker.check_reaction(rxn_type, rsmi):
                            heterocycle_forming_reactions += 1
                            print(f"Heterocycle-forming reaction detected: {rxn_type}")
                            break

                # Check if this reaction involves multiple heterocycles (assembly)
                if len(product_heterocycles) >= 2 and len(reactant_heterocycles) >= 1:
                    heterocycle_assembly_reactions += 1
                    print(
                        f"Heterocycle assembly reaction detected: {len(product_heterocycles)} heterocycles in product, {len(reactant_heterocycles)} in reactants"
                    )

            except (KeyError, IndexError) as e:
                print(f"Error processing reaction node: {e}")

        # Recursively traverse children
        for child in node.get("children", []):
            # If this node is a starting material, its children are also starting materials
            child_is_starting = is_starting_material or (
                node["type"] == "mol" and node.get("in_stock", False)
            )
            dfs_traverse(child, depth + 1, child_is_starting)

    # Start traversal
    dfs_traverse(route)

    # Calculate heterocycles formed during synthesis (those in final product but not in starting materials)
    formed_heterocycles = (
        final_product_heterocycles - heterocycles_in_starting_materials
    )

    # Strategy is present if at least 3 different heterocycles are in the final product
    # This is the primary criterion for multi-heterocycle assembly
    strategy_present = len(final_product_heterocycles) >= 3

    print(f"Multi-heterocycle assembly strategy detected: {strategy_present}")
    print(f"Heterocycles in final product: {final_product_heterocycles}")
    print(f"Heterocycles formed during synthesis: {formed_heterocycles}")
    print(f"Heterocycle-forming reactions detected: {heterocycle_forming_reactions}")
    print(f"Heterocycle assembly reactions detected: {heterocycle_assembly_reactions}")

    return strategy_present

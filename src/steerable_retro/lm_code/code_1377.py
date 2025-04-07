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
    This function detects if the synthesis involves formation of a heterocycle
    in the second half of the synthesis (late stage).
    """
    heterocycle_formation_depths = []
    max_depth = 0

    # List of common heterocycles to check
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

    # List of common heterocycle formation reactions
    heterocycle_formation_reactions = [
        "Formation of NOS Heterocycles",
        "Paal-Knorr pyrrole synthesis",
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
        "Huisgen alkyne-azide 1,3 dipolar cycloaddition",
        "Huisgen 1,3 dipolar cycloaddition",
        "Huisgen alkene-azide 1,3 dipolar cycloaddition",
        "Pyrazole formation",
        "Azide-nitrile click cycloaddition to tetrazole",
        "Azide-nitrile click cycloaddition to triazole",
        "Huisgen_Cu-catalyzed_1,4-subst",
        "Huisgen_Ru-catalyzed_1,5_subst",
        "Huisgen_disubst-alkyne",
        "Huisgen 1,3,4-oxadiazoles from COOH and tetrazole",
        "Pictet-Spengler",
    ]

    def dfs_traverse(node, depth=0):
        nonlocal max_depth

        if depth > max_depth:
            max_depth = depth

        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants_part = rsmi.split(">")[0]
                product_part = rsmi.split(">")[-1]

                # Check if this is a known heterocycle formation reaction
                if any(
                    checker.check_reaction(rxn_name, rsmi)
                    for rxn_name in heterocycle_formation_reactions
                ):
                    print(f"Known heterocycle formation reaction detected at depth {depth}: {rsmi}")
                    heterocycle_formation_depths.append(depth)
                else:
                    # Count heterocycles in reactants and product
                    reactants_smiles = reactants_part.split(".")
                    product_smiles = product_part

                    # Check heterocycles in reactants
                    reactant_heterocycles = set()
                    for r_smiles in reactants_smiles:
                        for heterocycle in heterocycles:
                            if checker.check_ring(heterocycle, r_smiles):
                                reactant_heterocycles.add(heterocycle)

                    # Check heterocycles in product
                    product_heterocycles = set()
                    for heterocycle in heterocycles:
                        if checker.check_ring(heterocycle, product_smiles):
                            product_heterocycles.add(heterocycle)

                    # If product has heterocycles not present in reactants, it's a heterocycle formation
                    new_heterocycles = product_heterocycles - reactant_heterocycles
                    if new_heterocycles:
                        print(
                            f"Heterocycle formation detected at depth {depth}: {new_heterocycles}"
                        )
                        heterocycle_formation_depths.append(depth)

                    # Alternative method: check if product has more heterocyclic rings than reactants
                    try:
                        reactants_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles]
                        product_mol = Chem.MolFromSmiles(product_smiles)

                        if None not in reactants_mols and product_mol is not None:
                            # Count heterocyclic rings in reactants
                            reactant_hetero_ring_count = 0
                            for mol in reactants_mols:
                                for ring in Chem.GetSSSR(mol):
                                    if any(
                                        mol.GetAtomWithIdx(atom_idx).GetSymbol() in ["N", "O", "S"]
                                        for atom_idx in ring
                                    ):
                                        reactant_hetero_ring_count += 1

                            # Count heterocyclic rings in product
                            product_hetero_ring_count = 0
                            for ring in Chem.GetSSSR(product_mol):
                                if any(
                                    product_mol.GetAtomWithIdx(atom_idx).GetSymbol()
                                    in ["N", "O", "S"]
                                    for atom_idx in ring
                                ):
                                    product_hetero_ring_count += 1

                            # If product has more heterocyclic rings than reactants combined
                            if (
                                product_hetero_ring_count > reactant_hetero_ring_count
                                and depth not in heterocycle_formation_depths
                            ):
                                print(
                                    f"Additional heterocycle formation detected at depth {depth} (ring count method)"
                                )
                                heterocycle_formation_depths.append(depth)
                    except Exception as e:
                        print(f"Error in ring counting: {e}")
            except Exception as e:
                print(f"Error processing reaction at depth {depth}: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if heterocycle formation occurs in the second half of synthesis
    if heterocycle_formation_depths and max_depth > 0:
        # In retrosynthetic direction, lower depths are later in the synthesis
        for depth in heterocycle_formation_depths:
            if depth <= max_depth / 2:  # Second half of synthesis
                print(f"Late-stage heterocycle formation detected at depth {depth}/{max_depth}")
                return True

    print(f"No late-stage heterocycle formation detected. Max depth: {max_depth}")
    return False

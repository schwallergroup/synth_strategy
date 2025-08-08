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
    Detects if the synthesis route includes formation of a heterocycle.
    """
    heterocycle_formation = False

    # List of common heterocycles to check (nitrogen, oxygen, and sulfur)
    heterocycles = [
        # Nitrogen heterocycles
        "pyrrole",
        "pyridine",
        "piperidine",
        "pyrrolidine",
        "imidazole",
        "pyrazole",
        "triazole",
        "tetrazole",
        "indole",
        "quinoline",
        "isoquinoline",
        "piperazine",
        "morpholine",
        "oxazole",
        "thiazole",
        "aziridine",
        "azetidine",
        "azepane",
        "diazepane",
        "benzimidazole",
        # Oxygen heterocycles
        "furan",
        "pyran",
        "dioxane",
        "tetrahydrofuran",
        "tetrahydropyran",
        "oxirane",
        "oxetane",
        "oxolane",
        "oxane",
        "benzoxazole",
        # Sulfur heterocycles
        "thiophene",
        "thiopyran",
        "thiirane",
        "thietane",
        "thiolane",
        "thiane",
        "benzothiophene",
        "benzothiazole",
    ]

    # List of reactions known to form heterocycles
    heterocycle_forming_reactions = [
        "Formation of NOS Heterocycles",
        "Paal-Knorr pyrrole synthesis",
        "benzimidazole_derivatives_aldehyde",
        "benzimidazole_derivatives_carboxylic-acid/ester",
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
        "benzofuran",
        "benzothiophene",
        "indole",
        "oxadiazole",
        "Huisgen alkyne-azide 1,3 dipolar cycloaddition",
        "Huisgen 1,3 dipolar cycloaddition",
        "Huisgen alkene-azide 1,3 dipolar cycloaddition",
        "Pyrazole formation",
        "Azide-nitrile click cycloaddition to tetrazole",
        "Azide-nitrile click cycloaddition to triazole",
    ]

    def dfs_traverse(node):
        nonlocal heterocycle_formation

        if heterocycle_formation:
            return  # Early return if we already found a heterocycle formation

        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # First check if this is a known heterocycle formation reaction
                for reaction_type in heterocycle_forming_reactions:
                    if checker.check_reaction(reaction_type, rsmi):
                        heterocycle_formation = True
                        return

                # Check if any heterocycle exists in the product
                product_heterocycles = []
                for ring in heterocycles:
                    if checker.check_ring(ring, product):
                        product_heterocycles.append(ring)

                if product_heterocycles:
                    # Check if any of the found heterocycles are present in any reactant
                    reactant_heterocycles = set()
                    for reactant in reactants:
                        for ring in heterocycles:
                            if checker.check_ring(ring, reactant):
                                reactant_heterocycles.add(ring)

                    # If there's a heterocycle in product that's not in any reactant, it was formed
                    for ring in product_heterocycles:
                        if ring not in reactant_heterocycles:
                            heterocycle_formation = True
                            return

                    # Check if the number of heterocycles increased
                    product_mol = Chem.MolFromSmiles(product)
                    if product_mol:
                        product_ring_count = 0
                        for ring in product_heterocycles:
                            product_ring_count += len(checker.get_ring_atom_indices(ring, product))

                        reactant_ring_count = 0
                        for reactant in reactants:
                            for ring in reactant_heterocycles:
                                if checker.check_ring(ring, reactant):
                                    reactant_ring_count += len(
                                        checker.get_ring_atom_indices(ring, reactant)
                                    )

                        if product_ring_count > reactant_ring_count:
                            heterocycle_formation = True
                            return
            except Exception as e:
                # Silently handle any exceptions and continue traversal
                pass

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return heterocycle_formation

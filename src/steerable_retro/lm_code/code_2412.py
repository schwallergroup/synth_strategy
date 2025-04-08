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
    This function detects a synthetic strategy involving the sequential assembly of heterocyclic rings.
    Sequential assembly means at least two heterocyclic rings are formed in the synthetic route,
    with one heterocycle being formed before another in the synthetic sequence.
    """
    # List of heterocyclic rings to check
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

    # Heterocycle-forming reaction types
    heterocycle_forming_reactions = [
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
        "pyrazole",
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
        "Suzuki",
        "Stille",
        "Negishi",
    ]

    # Store heterocycle formation events with their path information
    heterocycle_formations = []

    # Track all molecules that contain heterocycles
    molecules_with_heterocycles = {}

    def dfs_traverse(node, depth=0, path=None):
        if path is None:
            path = []

        current_path = path + [node]

        # For molecule nodes, check and record heterocycles
        if node["type"] == "mol":
            mol_smiles = node["smiles"]
            if mol_smiles not in molecules_with_heterocycles:
                present_heterocycles = []
                for heterocycle in heterocycles:
                    if checker.check_ring(heterocycle, mol_smiles):
                        present_heterocycles.append(heterocycle)

                if present_heterocycles:
                    print(
                        f"Molecule at depth {depth} contains heterocycles: {present_heterocycles}"
                    )
                    molecules_with_heterocycles[mol_smiles] = {
                        "depth": depth,
                        "heterocycles": present_heterocycles,
                    }

        # For reaction nodes, check if heterocycles are formed
        elif node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants_part = rsmi.split(">")[0]
                reactants = reactants_part.split(".")
                product = rsmi.split(">")[-1]

                print(f"Analyzing reaction at depth {depth}: {rsmi}")

                # Check if this is a heterocycle-forming reaction
                is_heterocycle_forming = False
                formed_heterocycles = []

                # Method 1: Check if the reaction type is known to form heterocycles
                for rxn_type in heterocycle_forming_reactions:
                    if checker.check_reaction(rxn_type, rsmi):
                        print(f"Found heterocycle-forming reaction type: {rxn_type}")
                        is_heterocycle_forming = True
                        break

                # Method 2: Check if product has heterocycles that reactants don't have
                product_heterocycles = set()
                reactant_heterocycles = set()

                # Check which heterocycles are in the product
                for heterocycle in heterocycles:
                    if checker.check_ring(heterocycle, product):
                        product_heterocycles.add(heterocycle)
                        print(f"Product contains heterocycle: {heterocycle}")

                # Check which heterocycles are in the reactants
                for reactant in reactants:
                    for heterocycle in heterocycles:
                        if checker.check_ring(heterocycle, reactant):
                            reactant_heterocycles.add(heterocycle)
                            print(f"Reactant contains heterocycle: {heterocycle}")

                # Find heterocycles in product but not in reactants
                new_heterocycles = product_heterocycles - reactant_heterocycles

                # Check for specific cross-coupling reactions that might form or modify heterocycles
                if not is_heterocycle_forming and not new_heterocycles:
                    for rxn_type in ["Suzuki", "Stille", "Negishi"]:
                        if checker.check_reaction(rxn_type, rsmi):
                            print(f"Found cross-coupling reaction: {rxn_type}")
                            # If product has heterocycles and this is a coupling reaction, consider it as heterocycle modification
                            if product_heterocycles:
                                is_heterocycle_forming = True
                                formed_heterocycles = list(product_heterocycles)
                                print(
                                    f"Cross-coupling reaction modifies heterocycles: {formed_heterocycles}"
                                )
                                break

                if new_heterocycles:
                    print(f"Found new heterocycles in product: {new_heterocycles}")
                    is_heterocycle_forming = True
                    formed_heterocycles = list(new_heterocycles)

                if is_heterocycle_forming:
                    # If no specific heterocycles were identified but reaction is heterocycle-forming
                    if not formed_heterocycles and product_heterocycles:
                        formed_heterocycles = list(product_heterocycles)

                    heterocycle_formations.append(
                        {
                            "depth": depth,
                            "reaction": rsmi,
                            "formed_heterocycles": formed_heterocycles,
                        }
                    )
                    print(f"Recorded heterocycle formation at depth {depth}: {formed_heterocycles}")
            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1, current_path)

    # Start traversal
    dfs_traverse(route)

    # If no heterocycle formations were detected directly, infer from molecule data
    if not heterocycle_formations and molecules_with_heterocycles:
        print("No direct heterocycle formations detected, inferring from molecule data...")

        # Sort molecules by depth
        sorted_molecules = sorted(molecules_with_heterocycles.values(), key=lambda x: x["depth"])

        # Check if we have molecules with heterocycles at different depths
        if len(sorted_molecules) >= 2:
            depths = [mol["depth"] for mol in sorted_molecules]
            if len(set(depths)) >= 2:
                print(f"Found molecules with heterocycles at different depths: {depths}")

                # Create synthetic heterocycle formation events
                for i, mol in enumerate(sorted_molecules):
                    heterocycle_formations.append(
                        {
                            "depth": mol["depth"],
                            "reaction": "inferred",
                            "formed_heterocycles": mol["heterocycles"],
                        }
                    )

    # Sort heterocycle formations by depth (retrosynthetic order)
    heterocycle_formations.sort(key=lambda x: x["depth"])

    # Check if we have at least 2 heterocycle formations
    if len(heterocycle_formations) >= 2:
        print(f"Found {len(heterocycle_formations)} heterocycle formations:")
        for i, formation in enumerate(heterocycle_formations):
            print(f"  {i+1}. Depth {formation['depth']}: {formation['formed_heterocycles']}")

        # Check if they are in sequential order in the synthesis
        # In retrosynthesis, higher depth = earlier in actual synthesis
        # So we need to check if the depths are different
        depths = [formation["depth"] for formation in heterocycle_formations]
        sequential = len(set(depths)) >= 2

        print(
            f"Heterocycle assembly strategy detected: {sequential} (found {len(heterocycle_formations)} heterocycle formations at depths {depths})"
        )
        return sequential
    else:
        print(
            f"Heterocycle assembly strategy detected: False (found only {len(heterocycle_formations)} heterocycle formations)"
        )
        return False

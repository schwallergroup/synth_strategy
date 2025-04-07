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
    This function detects if the synthetic route involves heterocycle formation
    through cyclization reactions.
    """
    # List of heterocyclic rings to check
    heterocyclic_rings = [
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
    ]

    # List of heterocycle formation reactions
    heterocycle_reactions = [
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
        "Huisgen alkyne-azide 1,3 dipolar cycloaddition",
        "Huisgen 1,3 dipolar cycloaddition",
        "Huisgen alkene-azide 1,3 dipolar cycloaddition",
        "Pyrazole formation",
        "Intramolecular amination of azidobiphenyls (heterocycle formation)",
        "Intramolecular amination (heterocycle formation)",
    ]

    heterocycle_formation_detected = False

    def dfs_traverse(node):
        nonlocal heterocycle_formation_detected

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0]
            product_smiles = rsmi.split(">")[-1]

            try:
                # Check if any heterocycle formation reaction is detected directly
                reaction_detected = any(
                    checker.check_reaction(rxn, rsmi) for rxn in heterocycle_reactions
                )
                if reaction_detected:
                    print(f"Heterocycle formation reaction detected: {rsmi}")
                    heterocycle_formation_detected = True

                # Check for heterocycle presence in reactants and products
                product_heterocycles = set()
                reactant_heterocycles = set()

                for ring in heterocyclic_rings:
                    if checker.check_ring(ring, product_smiles):
                        product_heterocycles.add(ring)
                    if checker.check_ring(ring, reactants_smiles):
                        reactant_heterocycles.add(ring)

                # Find new heterocycles in product
                new_heterocycles = product_heterocycles - reactant_heterocycles

                if new_heterocycles:
                    print(f"New heterocycles formed: {new_heterocycles}")
                    print(f"Reaction: {rsmi}")
                    heterocycle_formation_detected = True

                # Additional check for ring count increase
                try:
                    # Parse reactants and product as molecules
                    reactant_mols = [
                        Chem.MolFromSmiles(smi) for smi in reactants_smiles.split(".") if smi
                    ]
                    product_mols = [
                        Chem.MolFromSmiles(smi) for smi in product_smiles.split(".") if smi
                    ]

                    # Filter out None values (parsing failures)
                    reactant_mols = [mol for mol in reactant_mols if mol is not None]
                    product_mols = [mol for mol in product_mols if mol is not None]

                    if reactant_mols and product_mols:
                        # Count rings in reactants and products
                        reactant_ring_count = sum(
                            mol.GetRingInfo().NumRings() for mol in reactant_mols
                        )
                        product_ring_count = sum(
                            mol.GetRingInfo().NumRings() for mol in product_mols
                        )

                        # If product has more rings and no heterocycle was detected yet, check more carefully
                        if (
                            product_ring_count > reactant_ring_count
                            and not heterocycle_formation_detected
                        ):
                            # Check if any of the new rings are heterocyclic
                            for product_mol in product_mols:
                                for ring in heterocyclic_rings:
                                    if checker.check_ring(ring, Chem.MolToSmiles(product_mol)):
                                        print(
                                            f"Ring count increased with heterocycle formation: {rsmi}"
                                        )
                                        print(f"Heterocycle detected: {ring}")
                                        heterocycle_formation_detected = True
                                        break
                                if heterocycle_formation_detected:
                                    break
                except Exception as e:
                    print(f"Error in ring counting: {str(e)}")
                    # Continue with other checks even if ring counting fails

            except Exception as e:
                print(f"Error processing reaction SMILES: {rsmi}")
                print(f"Error details: {str(e)}")
                # Continue with other reactions even if one fails

        # Continue DFS traversal
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return heterocycle_formation_detected

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
    Detects if the synthesis involves formation of heterocyclic systems,
    particularly thiophene-thiazole fused systems.
    """
    forms_heterocycle = False

    def dfs_traverse(node, depth=0):
        nonlocal forms_heterocycle

        if node["type"] == "reaction" and depth >= 2:  # Early to middle stages
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if this is a known heterocycle-forming reaction
                heterocycle_reactions = [
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
                    "Huisgen alkyne-azide 1,3 dipolar cycloaddition",
                    "Huisgen 1,3 dipolar cycloaddition",
                    "Huisgen alkene-azide 1,3 dipolar cycloaddition",
                    "Pyrazole formation",
                    "Azide-nitrile click cycloaddition to tetrazole",
                    "Azide-nitrile click cycloaddition to triazole",
                    "Intramolecular amination of azidobiphenyls (heterocycle formation)",
                    "Intramolecular amination (heterocycle formation)",
                    "Formation of NOS Heterocycles",
                ]

                for reaction_type in heterocycle_reactions:
                    if checker.check_reaction(reaction_type, rsmi):
                        print(
                            f"Found heterocycle formation reaction: {reaction_type} at depth {depth}"
                        )
                        forms_heterocycle = True
                        return

                # Check reactants and products for ring count
                reactant_rings = 0
                product_mol = Chem.MolFromSmiles(product)

                if not product_mol:
                    return

                product_rings = product_mol.GetRingInfo().NumRings()

                # Check for heterocyclic rings in the product
                heterocyclic_rings = [
                    "thiophene",
                    "thiazole",
                    "furan",
                    "pyrrole",
                    "pyridine",
                    "pyrazole",
                    "imidazole",
                    "oxazole",
                    "pyrimidine",
                    "pyrazine",
                    "triazole",
                    "tetrazole",
                    "indole",
                    "benzimidazole",
                    "benzothiazole",
                    "benzoxazole",
                    "isoxazole",
                    "isothiazole",
                    "oxadiazole",
                    "thiadiazole",
                ]

                # Count rings in reactants
                for reactant in reactants:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol:
                        reactant_rings += mol.GetRingInfo().NumRings()

                # Check if product has more rings than reactants (ring formation)
                if product_rings > reactant_rings:
                    # Check if any of the new rings are heterocycles
                    for ring_type in heterocyclic_rings:
                        if checker.check_ring(ring_type, product):
                            # Verify this ring wasn't already in the reactants
                            ring_in_reactants = False
                            for reactant in reactants:
                                if checker.check_ring(ring_type, reactant):
                                    ring_in_reactants = True
                                    break

                            if not ring_in_reactants:
                                print(
                                    f"Found heterocycle formation: {ring_type} at depth {depth}"
                                )
                                forms_heterocycle = True
                                return

                # Check specifically for thiophene-thiazole fused systems
                if checker.check_ring("thiophene", product) and checker.check_ring(
                    "thiazole", product
                ):
                    # Check if both rings were not present in reactants
                    thiophene_in_reactants = False
                    thiazole_in_reactants = False

                    for reactant in reactants:
                        if checker.check_ring("thiophene", reactant):
                            thiophene_in_reactants = True
                        if checker.check_ring("thiazole", reactant):
                            thiazole_in_reactants = True

                    if not (thiophene_in_reactants and thiazole_in_reactants):
                        print(
                            f"Found thiophene-thiazole system formation at depth {depth}"
                        )
                        forms_heterocycle = True
                        return

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)
    return forms_heterocycle

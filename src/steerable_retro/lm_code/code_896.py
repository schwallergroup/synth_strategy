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


def main(route):
    """
    Detects a strategy where an aromatic ring is first functionalized (bromination, nitration)
    followed by heterocycle formation (piperazine) and late-stage N-alkylation.
    """
    # Track key transformations
    has_aromatic_bromination = False
    has_aromatic_nitration = False
    has_nitro_reduction = False
    has_piperazine_formation = False
    has_late_stage_alkylation = False

    def dfs_traverse(node, depth=0):
        nonlocal has_aromatic_bromination, has_aromatic_nitration, has_nitro_reduction
        nonlocal has_piperazine_formation, has_late_stage_alkylation

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for aromatic bromination (early stage, high depth)
                if depth >= 4 and not has_aromatic_bromination:
                    # Check if product has aromatic bromine but reactant doesn't
                    prod_mol = Chem.MolFromSmiles(product)
                    if prod_mol and prod_mol.HasSubstructMatch(Chem.MolFromSmarts("[c]-[Br]")):
                        for reactant in reactants:
                            react_mol = Chem.MolFromSmiles(reactant)
                            if react_mol and not react_mol.HasSubstructMatch(
                                Chem.MolFromSmarts("[c]-[Br]")
                            ):
                                has_aromatic_bromination = True
                                print("Detected aromatic bromination at depth", depth)
                                break

                # Check for aromatic nitration (early-mid stage)
                if depth >= 3 and not has_aromatic_nitration:
                    prod_mol = Chem.MolFromSmiles(product)
                    if prod_mol and prod_mol.HasSubstructMatch(
                        Chem.MolFromSmarts("[c]-[N+](=[O])[O-]")
                    ):
                        for reactant in reactants:
                            react_mol = Chem.MolFromSmiles(reactant)
                            if react_mol and not react_mol.HasSubstructMatch(
                                Chem.MolFromSmarts("[c]-[N+](=[O])[O-]")
                            ):
                                has_aromatic_nitration = True
                                print("Detected aromatic nitration at depth", depth)
                                break

                # Check for nitro reduction (mid stage)
                if depth >= 1 and not has_nitro_reduction:
                    prod_mol = Chem.MolFromSmiles(product)
                    if prod_mol and prod_mol.HasSubstructMatch(Chem.MolFromSmarts("[c]-[NH2]")):
                        for reactant in reactants:
                            react_mol = Chem.MolFromSmiles(reactant)
                            if react_mol and react_mol.HasSubstructMatch(
                                Chem.MolFromSmarts("[c]-[N+](=[O])[O-]")
                            ):
                                has_nitro_reduction = True
                                print("Detected nitro reduction at depth", depth)
                                break

                # Check for piperazine formation (mid stage)
                if depth >= 0 and not has_piperazine_formation:
                    prod_mol = Chem.MolFromSmiles(product)
                    if prod_mol and prod_mol.HasSubstructMatch(
                        Chem.MolFromSmarts("[N]1[CH2][CH2][NH][CH2][CH2]1")
                    ):
                        has_piperazine_formation = True
                        print("Detected piperazine formation at depth", depth)

                # Check for late-stage N-alkylation (low depth)
                if depth <= 1 and not has_late_stage_alkylation:
                    prod_mol = Chem.MolFromSmiles(product)
                    if prod_mol and prod_mol.HasSubstructMatch(
                        Chem.MolFromSmarts(
                            "[N]1[CH2][CH2][N]([CH2][C](=[O])[O][CH2][CH3])[CH2][CH2]1"
                        )
                    ):
                        has_late_stage_alkylation = True
                        print("Detected late-stage N-alkylation at depth", depth)

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Strategy is present if we have most of the key transformations
    required_count = 3  # Need at least 3 of the 5 key transformations
    found_count = sum(
        [
            has_aromatic_bromination,
            has_aromatic_nitration,
            has_nitro_reduction,
            has_piperazine_formation,
            has_late_stage_alkylation,
        ]
    )

    result = found_count >= required_count
    print(f"Aromatic functionalization to heterocycle strategy detected: {result}")
    return result

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
    This function detects a synthetic strategy involving initial fragment coupling between
    a nitrile-containing fragment and a phenyl ketone, followed by pyrrole formation.
    """
    # Initialize tracking variables
    has_fragment_coupling = False
    has_pyrrole_formation = False
    has_phenyl_ketone = False
    has_nitrile = False

    # SMARTS patterns
    nitrile_pattern = Chem.MolFromSmarts("[#6]-[#6]#[#7]")
    pyrrole_pattern = Chem.MolFromSmarts("[#6]1[#6][#6][#7H][#6]1")
    phenyl_ketone_pattern = Chem.MolFromSmarts("[#6]-[#6](=[#8])-c1ccccc1")

    def dfs_traverse(node, depth=0):
        nonlocal has_fragment_coupling, has_pyrrole_formation, has_phenyl_ketone, has_nitrile

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            try:
                product_mol = Chem.MolFromSmiles(product_smiles)
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles]

                # Check for fragment coupling at early stage (high depth)
                if depth >= 4 and len(reactant_mols) >= 2:
                    # Check if one fragment has nitrile and another has phenyl ketone
                    reactant_has_nitrile = [
                        r.HasSubstructMatch(nitrile_pattern) for r in reactant_mols if r
                    ]
                    reactant_has_phenyl_ketone = [
                        r.HasSubstructMatch(phenyl_ketone_pattern) for r in reactant_mols if r
                    ]

                    if any(reactant_has_nitrile) and any(reactant_has_phenyl_ketone):
                        has_fragment_coupling = True
                        has_nitrile = True
                        has_phenyl_ketone = True
                        print(
                            f"Detected fragment coupling between nitrile and phenyl ketone at depth {depth}"
                        )

                # Check for pyrrole formation
                if (
                    not any(r.HasSubstructMatch(pyrrole_pattern) for r in reactant_mols if r)
                    and product_mol
                    and product_mol.HasSubstructMatch(pyrrole_pattern)
                ):
                    has_pyrrole_formation = True
                    print(f"Detected pyrrole formation at depth {depth}")

            except Exception as e:
                print(f"Error processing reaction at depth {depth}: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)

    # Check if the strategy is present
    strategy_present = (
        has_fragment_coupling and has_pyrrole_formation and has_nitrile and has_phenyl_ketone
    )

    if strategy_present:
        print(
            "Detected fragment coupling pyrrole strategy with nitrile and phenyl ketone precursors"
        )

    return strategy_present

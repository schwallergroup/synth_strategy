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
    This function detects if the synthetic route employs a convergent fragment coupling strategy,
    specifically looking for nucleophilic aromatic substitution reactions.
    """
    fragment_coupling = False

    def dfs_traverse(node):
        nonlocal fragment_coupling

        if node["type"] == "reaction":
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if we have multiple reactant fragments (at least 2)
                if len(reactants) >= 2:
                    reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
                    product_mol = Chem.MolFromSmiles(product) if product else None

                    if product_mol and all(m for m in reactant_mols):
                        # Check for nucleophilic aromatic substitution pattern
                        # Look for Cl-aromatic in one fragment and NH2 in another
                        cl_aromatic_pattern = Chem.MolFromSmarts("c-[Cl]")
                        amine_pattern = Chem.MolFromSmarts("[NH2]")

                        has_cl_aromatic = any(
                            m.HasSubstructMatch(cl_aromatic_pattern) for m in reactant_mols
                        )
                        has_amine = any(m.HasSubstructMatch(amine_pattern) for m in reactant_mols)

                        # Check if product has new C-N bond where Cl was
                        if has_cl_aromatic and has_amine:
                            c_n_bond_pattern = Chem.MolFromSmarts("c-[#7]")
                            if product_mol.HasSubstructMatch(c_n_bond_pattern):
                                fragment_coupling = True
                                print(
                                    "Convergent fragment coupling via nucleophilic aromatic substitution detected"
                                )

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return fragment_coupling

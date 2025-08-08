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
    This function detects late-stage fragment coupling via C-N bond formation.
    """
    late_stage_coupling = False

    def dfs_traverse(node, depth=0):
        nonlocal late_stage_coupling

        if (
            node["type"] == "reaction" and "children" in node and depth <= 1
        ):  # Only check at depth 0 or 1 (late stage)
            # Get reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check if there are multiple reactants (fragment coupling)
            if len(reactants_smiles) >= 2:
                # Convert to RDKit molecules
                reactant_mols = [Chem.MolFromSmiles(smi) for smi in reactants_smiles if smi]
                product_mol = Chem.MolFromSmiles(product_smiles) if product_smiles else None

                if all(reactant_mols) and product_mol:
                    # Check for C-N bond formation
                    c_n_pattern = Chem.MolFromSmarts("[#6]-[#7]")

                    # Count C-N bonds in reactants and product
                    reactant_c_n_count = sum(
                        len(mol.GetSubstructMatches(c_n_pattern)) for mol in reactant_mols
                    )
                    product_c_n_count = len(product_mol.GetSubstructMatches(c_n_pattern))

                    # Check if product has more C-N bonds than reactants combined
                    if product_c_n_count > reactant_c_n_count:
                        # Check complexity of reactants (at least 10 atoms each)
                        if all(mol.GetNumAtoms() >= 10 for mol in reactant_mols):
                            print(f"Late-stage fragment coupling via C-N bond detected: {rsmi}")
                            late_stage_coupling = True

        # Traverse children with increased depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    return late_stage_coupling

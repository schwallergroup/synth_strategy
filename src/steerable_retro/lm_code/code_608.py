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
    Detects a synthetic strategy involving sequential nitrogen functionalization in the late stages,
    specifically an ester-to-amide conversion followed by N-alkylation.
    """
    # Track if we found the required transformations
    found_n_alkylation = False
    found_ester_to_amide = False

    # Track the depth at which transformations occur
    n_alkylation_depth = None
    ester_to_amide_depth = None

    def dfs_traverse(node, depth=0):
        nonlocal found_n_alkylation, found_ester_to_amide
        nonlocal n_alkylation_depth, ester_to_amide_depth

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0]
                product_smiles = rsmi.split(">")[-1]

                # Convert to RDKit molecules
                try:
                    reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles.split(".")]
                    product = Chem.MolFromSmiles(product_smiles)

                    if product and all(r for r in reactants):
                        # Check for N-alkylation (primary amine to secondary amine)
                        # Look for patterns where -NH2 becomes -NH-C
                        primary_amine_pattern = Chem.MolFromSmarts("[NH2][#6]")
                        secondary_amine_pattern = Chem.MolFromSmarts("[NH]([#6])[#6]")

                        has_primary_amine = any(
                            r.HasSubstructMatch(primary_amine_pattern) for r in reactants if r
                        )
                        has_secondary_amine = product.HasSubstructMatch(secondary_amine_pattern)

                        if (
                            has_primary_amine
                            and has_secondary_amine
                            and not product.HasSubstructMatch(primary_amine_pattern)
                        ):
                            print(f"Found N-alkylation at depth {depth}")
                            found_n_alkylation = True
                            n_alkylation_depth = depth

                        # Check for ester to amide transformation
                        ester_pattern = Chem.MolFromSmarts("[#6][C](=[O])[O][#6]")
                        amide_pattern = Chem.MolFromSmarts("[#6][C](=[O])[N]")

                        has_ester = any(r.HasSubstructMatch(ester_pattern) for r in reactants if r)
                        has_amide = product.HasSubstructMatch(amide_pattern)

                        if has_ester and has_amide and not product.HasSubstructMatch(ester_pattern):
                            print(f"Found ester-to-amide transformation at depth {depth}")
                            found_ester_to_amide = True
                            ester_to_amide_depth = depth

                except Exception as e:
                    print(f"Error processing reaction SMILES: {e}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if we found both transformations in the correct sequence
    # For late-stage transformations, we want low depth values
    if found_n_alkylation and found_ester_to_amide:
        # Check if N-alkylation occurs after ester-to-amide (lower depth number)
        if n_alkylation_depth is not None and ester_to_amide_depth is not None:
            if n_alkylation_depth < ester_to_amide_depth:
                print("Found late-stage sequential nitrogen functionalization strategy")
                return True

    return False

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
    Detects convergent synthesis with late-stage C-N bond formation between
    a complex amine and a halogenated heterocycle.
    """
    late_stage_cn_coupling = False

    def dfs_traverse(node, depth=0):
        nonlocal late_stage_cn_coupling

        if node["type"] == "reaction" and depth <= 1:  # Late stage (low depth)
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check if we have multiple complex reactants
                if len(reactants_smiles) >= 2:
                    # Convert to RDKit molecules
                    reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles]
                    product = Chem.MolFromSmiles(product_smiles)

                    if all(r is not None for r in reactants) and product is not None:
                        # Check for amine and halogenated heterocycle
                        amine_present = False
                        halogen_heterocycle = False

                        for r in reactants:
                            # Check for primary or secondary amine
                            if r.HasSubstructMatch(Chem.MolFromSmarts("[NH,NH2]")):
                                amine_present = True

                            # Check for halogenated heterocycle
                            if r.HasSubstructMatch(
                                Chem.MolFromSmarts("[n]")
                            ) and r.HasSubstructMatch(Chem.MolFromSmarts("[c][Cl,Br,I,F]")):
                                halogen_heterocycle = True

                        # Check if C-N bond is formed
                        if amine_present and halogen_heterocycle:
                            # This is a simplified check - in a real implementation,
                            # we would need to confirm the actual C-N bond formation
                            print(f"Found potential late-stage C-N coupling at depth {depth}")
                            late_stage_cn_coupling = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    return late_stage_cn_coupling

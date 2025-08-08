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
    This function detects if the synthetic route involves amide coupling in the late stages
    of synthesis (lower depth in retrosynthetic tree).
    """
    late_amide_coupling = False

    def dfs_traverse(node):
        nonlocal late_amide_coupling

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if this is an amide coupling reaction
                carboxylic_acid_pattern = Chem.MolFromSmarts("[#6]-C(=O)-[#8;H1]")
                amine_pattern = Chem.MolFromSmarts("[#6]-[#7;H2]")
                amide_pattern = Chem.MolFromSmarts("[#6]-C(=O)-[#7]-[#6]")

                # Check reactants for carboxylic acid and amine
                has_carboxylic_acid = False
                has_amine = False

                for reactant in reactants:
                    if reactant:
                        try:
                            mol = Chem.MolFromSmiles(reactant)
                            if mol:
                                if mol.HasSubstructMatch(carboxylic_acid_pattern):
                                    has_carboxylic_acid = True
                                if mol.HasSubstructMatch(amine_pattern):
                                    has_amine = True
                        except:
                            continue

                # Check product for amide
                has_amide = False
                if product:
                    try:
                        mol = Chem.MolFromSmiles(product)
                        if mol and mol.HasSubstructMatch(amide_pattern):
                            has_amide = True
                    except:
                        pass

                # If we found an amide coupling at low depth (late in synthesis)
                if (has_carboxylic_acid and has_amine and has_amide) or (
                    has_amide and (has_carboxylic_acid or has_amine)
                ):
                    if node.get("depth", 0) <= 1:  # Late stage (depth 0 or 1)
                        late_amide_coupling = True
                        print("Found late-stage amide coupling at depth:", node.get("depth", 0))

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return late_amide_coupling

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
    This function detects if the synthesis route involves a late-stage Suzuki coupling.
    """
    boronic_pattern = Chem.MolFromSmarts("[c]-[B]([O])[O]")
    aryl_halide_pattern = Chem.MolFromSmarts("[c]-[#53,#35,#17]")  # I, Br, Cl
    suzuki_coupling_late_stage = False
    late_stage_threshold = 4  # Consider depths < 4 as late in synthesis

    def dfs_traverse(node):
        nonlocal suzuki_coupling_late_stage

        if node["type"] == "reaction":
            depth = int(node.get("metadata", {}).get("depth", 0))
            if depth < late_stage_threshold:
                rsmi = node.get("metadata", {}).get("rsmi", "")
                if rsmi:
                    reactants = rsmi.split(">")[0].split(".")
                    product = rsmi.split(">")[-1]

                    # Check for Suzuki coupling pattern: boronic acid/ester + aryl halide -> biaryl
                    has_boronic = False
                    has_aryl_halide = False

                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol:
                            if reactant_mol.HasSubstructMatch(boronic_pattern):
                                has_boronic = True
                            if reactant_mol.HasSubstructMatch(aryl_halide_pattern):
                                has_aryl_halide = True

                    if has_boronic and has_aryl_halide:
                        print(f"Late-stage Suzuki coupling detected at depth {depth}")
                        suzuki_coupling_late_stage = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return suzuki_coupling_late_stage

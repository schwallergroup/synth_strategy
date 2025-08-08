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
    This function detects early-stage aromatic nucleophilic substitution (SNAr).
    """
    snar_detected = False

    def dfs_traverse(node):
        nonlocal snar_detected

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if this is an early-stage reaction (depth >= 5)
            if "depth" in node["metadata"] and node["metadata"].get("depth", -1) >= 5:
                try:
                    # Look for patterns indicating SNAr
                    # Typically involves F/Cl/etc. on aromatic being replaced with OR/NR2/etc.
                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol:
                            # SMARTS for aromatic with F/Cl
                            halo_aromatic = Chem.MolFromSmarts("c[F,Cl]")

                            # Check if reactant has halo-aromatic
                            if reactant_mol.HasSubstructMatch(halo_aromatic):
                                product_mol = Chem.MolFromSmiles(product)
                                # SMARTS for aromatic with OR
                                ether_aromatic = Chem.MolFromSmarts("cO[C,c]")

                                # Check if product has ether-aromatic
                                if product_mol and product_mol.HasSubstructMatch(ether_aromatic):
                                    print("Detected aromatic nucleophilic substitution")
                                    snar_detected = True
                                    break
                except:
                    pass

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return snar_detected

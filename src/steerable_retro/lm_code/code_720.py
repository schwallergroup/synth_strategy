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
    Detects if the synthesis route preserves a stereocenter throughout the synthesis.
    """
    stereocenter_preserved = False
    steps_with_stereocenter = 0
    total_steps = 0

    def dfs_traverse(node):
        nonlocal steps_with_stereocenter, total_steps

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            total_steps += 1
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for stereochemistry in reactants and product
            has_stereocenter_reactants = False
            has_stereocenter_product = False

            # Get stereocenters in reactants
            reactant_stereocenters = set()
            for r in reactants:
                mol = Chem.MolFromSmiles(r)
                if mol:
                    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=False)
                    if chiral_centers:
                        has_stereocenter_reactants = True
                        # Get atom mapping for each stereocenter
                        for atom_idx, _ in chiral_centers:
                            atom = mol.GetAtomWithIdx(atom_idx)
                            if atom.HasProp("molAtomMapNumber"):
                                map_num = atom.GetProp("molAtomMapNumber")
                                reactant_stereocenters.add(map_num)

            # Get stereocenters in product
            product_stereocenters = set()
            product_mol = Chem.MolFromSmiles(product)
            if product_mol:
                chiral_centers = Chem.FindMolChiralCenters(product_mol, includeUnassigned=False)
                if chiral_centers:
                    has_stereocenter_product = True
                    # Get atom mapping for each stereocenter
                    for atom_idx, _ in chiral_centers:
                        atom = product_mol.GetAtomWithIdx(atom_idx)
                        if atom.HasProp("molAtomMapNumber"):
                            map_num = atom.GetProp("molAtomMapNumber")
                            product_stereocenters.add(map_num)

            # Check if any stereocenter is preserved (has same mapping in reactants and product)
            preserved_stereocenters = reactant_stereocenters.intersection(product_stereocenters)

            if preserved_stereocenters:
                steps_with_stereocenter += 1
                # Extract depth safely
                match = re.search(r"Depth: (\d+)", node["metadata"].get("ID", "Depth: 999"))
                depth = int(match.group(1)) if match else 999
                print(
                    f"Stereocenter preserved at depth {depth}, preserved mappings: {preserved_stereocenters}"
                )

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    # If stereocenter is preserved in at least 75% of steps, consider it a stereocenter preservation strategy
    if total_steps > 0 and steps_with_stereocenter / total_steps >= 0.75:
        stereocenter_preserved = True
        print(f"Stereocenter preserved in {steps_with_stereocenter}/{total_steps} steps")

    return stereocenter_preserved

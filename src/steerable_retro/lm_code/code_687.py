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
from steerable_retro.utils.check import Check
from steerable_retro.utils import fuzzy_dict, check

root_data = "/home/dparm/steerable_retro/data"

fg_args = {
    "file_path": f"{root_data}/patterns/functional_groups.json",
    "value_field": "pattern",
    "key_field": "name",
}
reaction_class_args = {
    "file_path": f"{root_data}/patterns/smirks.json",
    "value_field": "smirks",
    "key_field": "name",
}
ring_smiles_args = {
    "file_path": f"{root_data}/patterns/chemical_rings_smiles.json",
    "value_field": "smiles",
    "key_field": "name",
}
functional_groups = fuzzy_dict.FuzzyDict.from_json(**fg_args)
reaction_classes = fuzzy_dict.FuzzyDict.from_json(**reaction_class_args)
ring_smiles = fuzzy_dict.FuzzyDict.from_json(**ring_smiles_args)

checker = check.Check(
    fg_dict=functional_groups, reaction_dict=reaction_classes, ring_dict=ring_smiles
)


def main(route):
    """
    Detects a strategy where stereocenters are formed during the synthesis.
    """
    # Track if we found stereocenter formation
    stereocenter_formation_found = False

    def dfs_traverse(node, depth=0):
        nonlocal stereocenter_formation_found

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants_part = rsmi.split(">")[0]
                product_part = rsmi.split(">")[-1]

                # Check for stereocenter formation
                try:
                    reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_part.split(".") if r]
                    product_mol = Chem.MolFromSmiles(product_part) if product_part else None

                    if all(reactant_mols) and product_mol:
                        # Count stereocenters in reactants and product
                        reactant_stereocenters = sum(
                            len(Chem.FindMolChiralCenters(mol, includeUnassigned=True))
                            for mol in reactant_mols
                        )
                        product_stereocenters = len(
                            Chem.FindMolChiralCenters(product_mol, includeUnassigned=True)
                        )

                        print(
                            f"Depth {depth}: Reactant stereocenters: {reactant_stereocenters}, Product stereocenters: {product_stereocenters}"
                        )

                        # Check for stereocenter formation (in retrosynthesis, reactants have more stereocenters)
                        if reactant_stereocenters > product_stereocenters:
                            stereocenter_formation_found = True
                            print(f"Stereocenter formation detected at depth {depth}")

                        # Check for specific stereoselective reactions
                        if (
                            checker.check_reaction("Diels-Alder", rsmi)
                            or checker.check_reaction("aldol condensation", rsmi)
                            or checker.check_reaction("Sharpless epoxidation", rsmi)
                            or checker.check_reaction("asymmetric hydrogenation", rsmi)
                        ):
                            stereocenter_formation_found = True
                            print(f"Stereocenter-forming reaction detected at depth {depth}")
                except Exception as e:
                    print(f"Error processing reaction at depth {depth}: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(f"Stereocenter formation strategy detection result: {stereocenter_formation_found}")
    return stereocenter_formation_found

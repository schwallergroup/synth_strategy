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
    This function detects if the synthetic route uses a cross-coupling reaction
    to form a biaryl system.
    """
    cross_coupling_found = False

    def dfs_traverse(node):
        nonlocal cross_coupling_found

        if node["type"] == "reaction":
            rsmi = node["metadata"].get("rsmi", "")
            if not rsmi:
                return

            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for cross-coupling reactions directly
            is_cross_coupling = (
                checker.check_reaction("Suzuki coupling with boronic acids", rsmi)
                or checker.check_reaction("Suzuki coupling with boronic esters", rsmi)
                or checker.check_reaction("Suzuki coupling with boronic acids OTf", rsmi)
                or checker.check_reaction("Suzuki coupling with boronic esters OTf", rsmi)
                or checker.check_reaction("Stille reaction_aryl", rsmi)
                or checker.check_reaction("Stille reaction_vinyl", rsmi)
                or checker.check_reaction("Stille reaction_benzyl", rsmi)
                or checker.check_reaction("Stille reaction_allyl", rsmi)
                or checker.check_reaction("Stille reaction_aryl OTf", rsmi)
                or checker.check_reaction("Negishi", rsmi)
                or checker.check_reaction("Hiyama-Denmark Coupling", rsmi)
                or checker.check_reaction("Kumada cross-coupling", rsmi)
                or checker.check_reaction("Aryllithium cross-coupling", rsmi)
            )

            # Check if product contains a biaryl system
            product_mol = Chem.MolFromSmiles(product_smiles)
            if product_mol and is_cross_coupling:
                # Check for two connected aromatic rings (using 'a' for any aromatic atom)
                biaryl_pattern = Chem.MolFromSmarts("a:a-a:a")
                if product_mol.GetSubstructMatches(biaryl_pattern):
                    print(f"Cross-coupling biaryl formation detected: {rsmi}")
                    cross_coupling_found = True

            # Fallback check for cross-coupling pattern if reaction type check fails
            if not cross_coupling_found:
                aryl_halide_found = False
                organometallic_found = False

                for reactant in reactants_smiles:
                    if checker.check_fg("Aromatic halide", reactant):
                        aryl_halide_found = True
                        print(f"Found aromatic halide: {reactant}")

                    # Check for various organometallic reagents
                    if (
                        checker.check_fg("Boronic acid", reactant)
                        or checker.check_fg("Boronic ester", reactant)
                        or "Sn" in reactant
                        or "ZnCl" in reactant
                        or checker.check_fg("Magnesium halide", reactant)
                    ):
                        organometallic_found = True
                        print(f"Found organometallic compound: {reactant}")

                # Check if product contains a biaryl system
                if product_mol and aryl_halide_found and organometallic_found:
                    biaryl_pattern = Chem.MolFromSmarts("a:a-a:a")
                    if product_mol.GetSubstructMatches(biaryl_pattern):
                        print(f"Cross-coupling biaryl formation detected (pattern-based): {rsmi}")
                        cross_coupling_found = True

                    # Additional check for biaryl system with more specific pattern
                    biaryl_pattern2 = Chem.MolFromSmarts("c:c-c:c")
                    if product_mol.GetSubstructMatches(biaryl_pattern2):
                        print(f"Carbon-based biaryl formation detected: {rsmi}")
                        cross_coupling_found = True

                    # Check for heteroaromatic biaryl systems
                    hetero_biaryl = Chem.MolFromSmarts("c:c-[n,o,s]:c")
                    if product_mol.GetSubstructMatches(hetero_biaryl):
                        print(f"Heteroaromatic biaryl formation detected: {rsmi}")
                        cross_coupling_found = True

        # Traverse children
        for child in node.get("children", []):
            if (
                not cross_coupling_found
            ):  # Stop traversal if we already found what we're looking for
                dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return cross_coupling_found

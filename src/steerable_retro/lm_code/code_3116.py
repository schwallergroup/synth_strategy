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
    Detects if there's a transformation from nitrile to lactone in the route.
    This can happen across multiple steps, not necessarily in a single reaction.
    """
    # Track nitrile-containing molecules and lactone-containing molecules with their depths
    nitrile_molecules = []
    lactone_molecules = []
    nitrile_to_lactone_found = False

    def dfs(node, depth=0):
        nonlocal nitrile_to_lactone_found

        # Check molecule nodes
        if node.get("type") == "mol" and "smiles" in node:
            mol_smiles = node["smiles"]

            # Track molecules with nitrile groups
            if checker.check_fg("Nitrile", mol_smiles):
                nitrile_molecules.append((depth, mol_smiles))
                print(f"Molecule with nitrile found at depth {depth}: {mol_smiles}")

            # Track molecules with lactone structures
            is_lactone = False
            if checker.check_fg("Ester", mol_smiles):
                for ring in [
                    "oxolane",
                    "oxane",
                    "tetrahydrofuran",
                    "tetrahydropyran",
                    "furan",
                    "pyran",
                ]:
                    if checker.check_ring(ring, mol_smiles):
                        is_lactone = True
                        break

                # Also check for coumarin structure
                if "O=c1occc2" in mol_smiles:
                    is_lactone = True

            if is_lactone:
                lactone_molecules.append((depth, mol_smiles))
                print(f"Molecule with lactone found at depth {depth}: {mol_smiles}")

            # Check if we have both nitrile and lactone molecules in the right order
            # Nitrile should be at higher depth (earlier in synthesis) than lactone
            if nitrile_molecules and lactone_molecules:
                nitrile_depths = [d for d, _ in nitrile_molecules]
                lactone_depths = [d for d, _ in lactone_molecules]

                if max(nitrile_depths) > min(lactone_depths):
                    nitrile_to_lactone_found = True
                    print(
                        f"Nitrile to lactone transformation detected: nitrile at depth {max(nitrile_depths)}, lactone at depth {min(lactone_depths)}"
                    )

        # Also check reaction nodes for direct transformation
        if (
            node.get("type") == "reaction"
            and "metadata" in node
            and "rsmi" in node.get("metadata", {})
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if any reactant has a nitrile group
            has_nitrile = any(checker.check_fg("Nitrile", reactant) for reactant in reactants)

            # Check if product has a lactone (cyclic ester)
            has_lactone = False
            if checker.check_fg("Ester", product):
                for ring in [
                    "oxolane",
                    "oxane",
                    "tetrahydrofuran",
                    "tetrahydropyran",
                    "furan",
                    "pyran",
                ]:
                    if checker.check_ring(ring, product):
                        has_lactone = True
                        break

                # Also check for coumarin structure
                if "O=c1occc2" in product:
                    has_lactone = True

            if has_nitrile and has_lactone:
                print(f"Direct nitrile to lactone transformation detected: {rsmi}")
                nitrile_to_lactone_found = True

        # Recursively check children
        for child in node.get("children", []):
            dfs(child, depth + 1)

    # Start DFS from the root
    dfs(route)

    return nitrile_to_lactone_found

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
    This function detects if the synthetic route involves a biaryl formation via cross-coupling.
    Looks for reactions where a boronic acid and an aryl halide form a biaryl system.
    """
    found_cross_coupling = False

    def dfs_traverse(node):
        nonlocal found_cross_coupling

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]

            # Check if this is a known cross-coupling reaction
            if (
                checker.check_reaction("Suzuki coupling with boronic acids", rsmi)
                or checker.check_reaction("Suzuki coupling with boronic acids OTf", rsmi)
                or checker.check_reaction("Suzuki coupling with boronic esters", rsmi)
                or checker.check_reaction("Suzuki coupling with boronic esters OTf", rsmi)
                or checker.check_reaction("Negishi coupling", rsmi)
                or checker.check_reaction("Stille reaction_aryl", rsmi)
                or checker.check_reaction("Stille reaction_aryl OTf", rsmi)
                or checker.check_reaction("Hiyama-Denmark Coupling", rsmi)
                or checker.check_reaction("Kumada cross-coupling", rsmi)
                or checker.check_reaction("Aryllithium cross-coupling", rsmi)
            ):

                print(f"Found biaryl formation via cross-coupling: {rsmi}")
                found_cross_coupling = True
            else:
                # If not a known reaction type, check for the components manually
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if one reactant contains boronic acid/ester and another contains a halide/triflate
                has_boronic_compound = False
                has_leaving_group = False

                for reactant in reactants:
                    if not reactant:
                        continue

                    if checker.check_fg("Boronic acid", reactant) or checker.check_fg(
                        "Boronic ester", reactant
                    ):
                        has_boronic_compound = True

                    if checker.check_fg("Aromatic halide", reactant) or checker.check_fg(
                        "Triflate", reactant
                    ):
                        has_leaving_group = True

                # Check if product contains a biaryl system
                if has_boronic_compound and has_leaving_group:
                    prod_mol = Chem.MolFromSmiles(product)
                    if prod_mol and prod_mol.HasSubstructMatch(Chem.MolFromSmarts("c:c-c:c")):
                        print(f"Found biaryl formation via components: {rsmi}")
                        found_cross_coupling = True

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return found_cross_coupling

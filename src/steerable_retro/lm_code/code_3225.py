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
    This function detects a synthetic strategy involving coupling of heterocyclic fragments
    (benzothiophene and benzodioxole).
    """
    has_benzothiophene = False
    has_benzodioxole = False
    has_cc_coupling = False
    has_heterocycle_coupling = False

    def dfs_traverse(node):
        nonlocal has_benzothiophene, has_benzodioxole, has_cc_coupling, has_heterocycle_coupling

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Get depth from metadata
            depth = -1
            if "ID" in node["metadata"]:
                depth_str = node["metadata"]["ID"]
                if "Depth:" in depth_str:
                    try:
                        depth = int(depth_str.split("Depth:")[1].split()[0])
                    except:
                        pass

            # Check for C-C bond formation (expanded depth range)
            if depth >= 0 and depth <= 10:
                print(f"Checking reaction at depth {depth}: {rsmi}")

                # Check for C-C coupling reactions (expanded list)
                if (
                    checker.check_reaction("Suzuki coupling with boronic acids", rsmi)
                    or checker.check_reaction("Suzuki coupling with boronic esters", rsmi)
                    or checker.check_reaction("Negishi coupling", rsmi)
                    or checker.check_reaction("Stille reaction_aryl", rsmi)
                    or checker.check_reaction("Ullmann condensation", rsmi)
                    or checker.check_reaction("Heck terminal vinyl", rsmi)
                    or checker.check_reaction("Sonogashira alkyne_aryl halide", rsmi)
                    or checker.check_reaction("Kumada cross-coupling", rsmi)
                ):

                    has_cc_coupling = True
                    print(f"Detected C-C coupling reaction at depth {depth}")

                    # Check for benzothiophene in product
                    if checker.check_ring("benzothiophene", product):
                        has_benzothiophene = True
                        print("Found benzothiophene in product")

                    # Check for benzodioxole (1,3-benzodioxole pattern)
                    if "OCO" in product:
                        has_benzodioxole = True
                        print("Found benzodioxole in product")

                    # Check if heterocycles are in separate reactants (indicating coupling)
                    benzothiophene_reactant = None
                    benzodioxole_reactant = None

                    for i, reactant in enumerate(reactants):
                        if checker.check_ring("benzothiophene", reactant):
                            benzothiophene_reactant = i
                            print(f"Found benzothiophene in reactant {i}: {reactant}")

                        if "OCO" in reactant:
                            benzodioxole_reactant = i
                            print(f"Found benzodioxole in reactant {i}: {reactant}")

                    # If heterocycles are found in different reactants, this is a coupling of the two
                    if benzothiophene_reactant is not None and benzodioxole_reactant is not None:
                        if benzothiophene_reactant != benzodioxole_reactant:
                            print("Heterocycles found in different reactants - confirmed coupling")
                            has_heterocycle_coupling = True
                        else:
                            print(
                                "Both heterocycles found in same reactant - not a coupling of the two"
                            )

                    # If only one heterocycle is found in reactants but both are in product,
                    # this could be a coupling where one heterocycle is formed during the reaction
                    elif (
                        (benzothiophene_reactant is not None or benzodioxole_reactant is not None)
                        and has_benzothiophene
                        and has_benzodioxole
                    ):
                        print(
                            "One heterocycle in reactants, both in product - possible coupling with heterocycle formation"
                        )
                        has_heterocycle_coupling = True

        # Check if the current molecule node contains both heterocycles
        if node["type"] == "mol" and node["smiles"]:
            mol_smiles = node["smiles"]
            if checker.check_ring("benzothiophene", mol_smiles):
                has_benzothiophene = True
                print(f"Found benzothiophene in molecule: {mol_smiles}")

            if "OCO" in mol_smiles:
                has_benzodioxole = True
                print(f"Found benzodioxole in molecule: {mol_smiles}")

                # If we find both heterocycles in the same molecule, check if it's a final product
                if has_benzothiophene and not node.get("children"):
                    print("Found both heterocycles in final product molecule")
                    # This suggests they were coupled at some point
                    has_heterocycle_coupling = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    # Return True if all conditions are met
    # We need both heterocycles, a C-C coupling, and evidence they were coupled together
    result = (
        has_benzothiophene and has_benzodioxole and (has_cc_coupling or has_heterocycle_coupling)
    )
    print(
        f"Final result: benzothiophene={has_benzothiophene}, benzodioxole={has_benzodioxole}, cc_coupling={has_cc_coupling}, heterocycle_coupling={has_heterocycle_coupling}"
    )
    return result

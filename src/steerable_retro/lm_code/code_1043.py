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
    This function detects if the synthetic route involves the formation of
    a heterocyclic ring, specifically focusing on rings like isoxazole.
    """
    # List of heterocyclic rings to check for formation
    heterocyclic_rings = [
        "furan",
        "pyran",
        "dioxane",
        "tetrahydrofuran",
        "tetrahydropyran",
        "oxirane",
        "oxetane",
        "oxolane",
        "oxane",
        "dioxolane",
        "dioxolene",
        "trioxane",
        "dioxepane",
        "pyrrole",
        "pyridine",
        "pyrazole",
        "imidazole",
        "oxazole",
        "isoxazole",
        "thiazole",
        "pyrimidine",
        "pyrazine",
        "pyridazine",
        "triazole",
        "tetrazole",
        "pyrrolidine",
        "piperidine",
        "piperazine",
        "morpholine",
        "thiomorpholine",
        "aziridine",
        "azetidine",
        "azepane",
        "diazepane",
        "indole",
        "quinoline",
        "isoquinoline",
        "purine",
        "carbazole",
        "acridine",
        "thiophene",
        "thiopyran",
        "thiirane",
        "thietane",
        "thiolane",
        "thiane",
        "dithiane",
        "dithiolane",
        "benzothiophene",
        "oxathiolane",
        "dioxathiolane",
        "thiazolidine",
        "oxazolidine",
        "oxadiazole",
        "thiadiazole",
    ]

    heterocycle_formation_detected = False

    def dfs_traverse(node):
        nonlocal heterocycle_formation_detected

        if node["type"] == "reaction" and not heterocycle_formation_detected:
            try:
                if "rsmi" in node["metadata"]:
                    rsmi = node["metadata"]["rsmi"]
                    reactants = rsmi.split(">")[0].split(".")
                    product = rsmi.split(">")[-1]

                    # Check for heterocycle formation in each ring type
                    for ring_name in heterocyclic_rings:
                        # Check if ring exists in product
                        if checker.check_ring(ring_name, product):
                            print(f"Found {ring_name} ring in product: {product}")

                            # Check if ring was not present in any reactant
                            ring_in_reactants = False
                            for reactant in reactants:
                                if reactant and checker.check_ring(ring_name, reactant):
                                    ring_in_reactants = True
                                    print(
                                        f"Ring {ring_name} already present in reactant: {reactant}"
                                    )
                                    break

                            # If ring is in product but not in reactants, it was formed
                            if not ring_in_reactants:
                                print(
                                    f"Heterocycle formation detected: {ring_name} ring was formed"
                                )
                                heterocycle_formation_detected = True
                                return
            except Exception as e:
                print(f"Error in reaction analysis: {e}")

        # Continue traversing
        for child in node.get("children", []):
            if not heterocycle_formation_detected:
                dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return heterocycle_formation_detected

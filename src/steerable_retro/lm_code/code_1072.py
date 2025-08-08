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
    Detects if the synthesis route involves modification of a heterocyclic core.
    Specifically looks for transformations of heterocycles like benzimidazole, pyridine, etc.
    """
    heterocycle_modification_found = False

    # List of heterocycles to check
    heterocycles = [
        "benzimidazole",
        "pyridine",
        "pyrrole",
        "furan",
        "thiophene",
        "imidazole",
        "oxazole",
        "thiazole",
        "pyrazole",
        "isoxazole",
        "isothiazole",
        "triazole",
        "tetrazole",
        "pyrimidine",
        "pyrazine",
        "pyridazine",
        "indole",
        "benzofuran",
        "benzothiophene",
        "quinoline",
        "isoquinoline",
        "piperidine",
        "morpholine",
        "piperazine",
        "pyrrolidine",
    ]

    def dfs_traverse(node):
        nonlocal heterocycle_modification_found

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                product_mol = Chem.MolFromSmiles(product)
                if not product_mol:
                    return

                # Check if product contains any heterocycle
                product_heterocycles = []
                for heterocycle in heterocycles:
                    if checker.check_ring(heterocycle, product):
                        product_heterocycles.append(heterocycle)
                        print(f"Product contains {heterocycle}")

                if not product_heterocycles:
                    return

                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if not reactant_mol:
                        continue

                    # Check if reactant contains any heterocycle
                    reactant_heterocycles = []
                    for heterocycle in heterocycles:
                        if checker.check_ring(heterocycle, reactant):
                            reactant_heterocycles.append(heterocycle)
                            print(f"Reactant contains {heterocycle}")

                    # Case 1: Heterocycle is present in both reactant and product
                    common_heterocycles = set(reactant_heterocycles).intersection(
                        set(product_heterocycles)
                    )
                    if common_heterocycles:
                        print(f"Common heterocycles: {common_heterocycles}")

                        # Check for functional group changes on the heterocycle
                        # This is a simplification - ideally we would use atom mapping to track specific changes
                        functional_groups = [
                            "Primary amine",
                            "Secondary amine",
                            "Tertiary amine",
                            "Primary alcohol",
                            "Secondary alcohol",
                            "Tertiary alcohol",
                            "Carboxylic acid",
                            "Ester",
                            "Amide",
                            "Nitrile",
                            "Halide",
                            "Nitro group",
                            "Sulfonamide",
                        ]

                        reactant_fgs = set()
                        product_fgs = set()

                        for fg in functional_groups:
                            if checker.check_fg(fg, reactant):
                                reactant_fgs.add(fg)
                            if checker.check_fg(fg, product):
                                product_fgs.add(fg)

                        if reactant_fgs != product_fgs:
                            print(
                                f"Functional group change detected: {reactant_fgs} -> {product_fgs}"
                            )
                            heterocycle_modification_found = True

                        # Check for specific reactions that modify heterocycles
                        reaction_types = [
                            "N-alkylation of primary amines with alkyl halides",
                            "N-alkylation of secondary amines with alkyl halides",
                            "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                            "Friedel-Crafts alkylation",
                            "Friedel-Crafts acylation",
                        ]

                        for rxn_type in reaction_types:
                            if checker.check_reaction(rxn_type, rsmi):
                                print(f"Heterocycle-modifying reaction detected: {rxn_type}")
                                heterocycle_modification_found = True

                    # Case 2: Heterocycle is formed in the reaction
                    if not reactant_heterocycles and product_heterocycles:
                        print(f"Heterocycle formation detected: {product_heterocycles}")
                        heterocycle_modification_found = True

                    # Case 3: Heterocycle is broken in the reaction
                    if reactant_heterocycles and not product_heterocycles:
                        print(f"Heterocycle breaking detected: {reactant_heterocycles}")
                        heterocycle_modification_found = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return heterocycle_modification_found

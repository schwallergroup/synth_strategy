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
    Detects synthesis involving multiple types of nitrogen heterocycles.
    """
    # List of nitrogen-containing heterocycles to check
    nitrogen_heterocycles = [
        "pyrrole",
        "pyridine",
        "pyrazole",
        "imidazole",
        "oxazole",
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
        "thiazolidine",
        "oxazolidine",
        "isoxazole",
        "oxadiazole",
        "thiadiazole",
        "benzoxazole",
        "benzothiazole",
        "benzimidazole",
        "pteridin",
        "phenothiazine",
        "phenoxazine",
        "pyrroline",
        "imidazolidine",
        "indazole",
        "benzotriazole",
    ]

    # Track heterocycles found in molecules and formed in reactions
    heterocycle_types = set()
    heterocycle_reactions = set()

    def dfs_traverse(node, depth=0):
        if node["type"] == "mol":
            if "smiles" in node:
                mol_smiles = node["smiles"]
                # Check for nitrogen heterocycles in the molecule
                for heterocycle in nitrogen_heterocycles:
                    if checker.check_ring(heterocycle, mol_smiles):
                        print(f"Found {heterocycle} in molecule: {mol_smiles}")
                        heterocycle_types.add(heterocycle)

        elif node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if any nitrogen heterocycle is formed in the reaction
                for heterocycle in nitrogen_heterocycles:
                    product_has_heterocycle = checker.check_ring(heterocycle, product)

                    # Check if the heterocycle is newly formed (not present in reactants)
                    if product_has_heterocycle:
                        reactants_have_heterocycle = any(
                            checker.check_ring(heterocycle, reactant) for reactant in reactants
                        )

                        if not reactants_have_heterocycle:
                            print(f"Heterocycle {heterocycle} formed in reaction: {rsmi}")
                            heterocycle_reactions.add(heterocycle)
                            heterocycle_types.add(heterocycle)
                        elif product_has_heterocycle:
                            # Even if present in reactants, still count it as part of the strategy
                            heterocycle_types.add(heterocycle)
            except Exception as e:
                print(f"Error analyzing reaction: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if we have at least 2 different types of nitrogen heterocycles
    if len(heterocycle_types) >= 2:
        print(f"Found multiple nitrogen heterocycle types: {heterocycle_types}")
        if heterocycle_reactions:
            print(f"Heterocycles formed in reactions: {heterocycle_reactions}")
        return True

    print(f"Found only {len(heterocycle_types)} nitrogen heterocycle types: {heterocycle_types}")
    return False

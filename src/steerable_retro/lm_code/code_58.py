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
    This function detects if the synthetic route employs a ring transformation strategy.
    It looks for reactions where rings are formed, opened, or transformed between reactants and products.
    """
    found_ring_transformation = False

    # List of common ring structures to check
    ring_types = [
        "furan",
        "pyran",
        "pyrrole",
        "pyridine",
        "benzene",
        "cyclohexane",
        "cyclopentane",
        "indole",
        "naphthalene",
        "thiophene",
        "imidazole",
        "oxazole",
        "thiazole",
        "pyrimidine",
        "piperidine",
        "morpholine",
    ]

    # List of ring-forming or ring-breaking reaction types
    ring_reactions = [
        "Diels-Alder",
        "Paal-Knorr pyrrole synthesis",
        "Fischer indole",
        "Friedlaender chinoline",
        "benzofuran",
        "benzothiophene",
        "indole",
        "Pictet-Spengler",
        "Niementowski_quinazoline",
    ]

    def dfs_traverse(node, depth=0):
        nonlocal found_ring_transformation

        if found_ring_transformation:
            return  # Early exit if we already found a ring transformation

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            print(f"Analyzing reaction at depth {depth}: {rsmi}")

            # Check if this is a known ring-forming or ring-breaking reaction
            for rxn_type in ring_reactions:
                if checker.check_reaction(rxn_type, rsmi):
                    print(f"Found ring transformation reaction type: {rxn_type}")
                    found_ring_transformation = True
                    return

            # Extract reactants and product
            try:
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                prod_mol = Chem.MolFromSmiles(product)
                if not prod_mol:
                    print(f"Could not parse product SMILES: {product}")
                    return

                # Count rings in product using a more reliable method
                prod_ring_count = rdMolDescriptors.CalcNumRings(prod_mol)

                # Check for specific ring structures in product
                prod_ring_types = []
                for ring in ring_types:
                    if checker.check_ring(ring, product):
                        prod_ring_types.append(ring)

                # Process reactants
                reactant_ring_count = 0
                reactant_ring_types = set()

                for reactant in reactants:
                    try:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol:
                            reactant_ring_count += rdMolDescriptors.CalcNumRings(mol)
                            for ring in ring_types:
                                if checker.check_ring(ring, reactant):
                                    reactant_ring_types.add(ring)
                    except Exception as e:
                        print(f"Error processing reactant {reactant}: {e}")
                        continue

                # Check if the number of rings has changed
                if prod_ring_count != reactant_ring_count:
                    print(f"Found ring transformation reaction: {rsmi}")
                    print(
                        f"Reactant rings: {reactant_ring_count}, Product rings: {prod_ring_count}"
                    )
                    found_ring_transformation = True
                    return

                # Check if ring types have changed (even if count remains the same)
                prod_ring_set = set(prod_ring_types)
                if prod_ring_set != reactant_ring_types and (prod_ring_set or reactant_ring_types):
                    print(f"Found ring type transformation: {rsmi}")
                    print(
                        f"Reactant ring types: {reactant_ring_types}, Product ring types: {prod_ring_set}"
                    )
                    found_ring_transformation = True
                    return

            except Exception as e:
                print(f"Error analyzing reaction {rsmi}: {e}")

        # Continue DFS traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return found_ring_transformation

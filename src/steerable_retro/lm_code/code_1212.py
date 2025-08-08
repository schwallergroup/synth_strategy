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
    Detects if the synthetic route involves ring formation via cyclization,
    particularly focusing on heterocyclic ring formation.
    """
    ring_formation_found = False

    # List of heterocyclic rings to check
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
        "indole",
        "quinoline",
        "isoquinoline",
        "benzoxazole",
        "benzothiazole",
        "benzimidazole",
    ]

    def dfs_traverse(node):
        nonlocal ring_formation_found

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            print(f"Analyzing reaction: {rsmi}")

            # Check if any heterocyclic ring is formed in the product
            product_mol = Chem.MolFromSmiles(product)
            if product_mol:
                product_rings = product_mol.GetRingInfo().NumRings()

                # Count rings in reactants
                reactant_rings_total = 0
                reactant_mols = []
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol:
                        reactant_mols.append(reactant_mol)
                        reactant_rings_total += reactant_mol.GetRingInfo().NumRings()

                print(
                    f"Rings in reactants: {reactant_rings_total}, Rings in product: {product_rings}"
                )

                # Check if product has more rings than reactants combined
                if product_rings > reactant_rings_total:
                    print("Product has more rings than reactants combined")

                    # Check for heterocyclic ring formation
                    for ring_name in heterocyclic_rings:
                        # Check if ring is in product
                        if checker.check_ring(ring_name, product):
                            # Check if ring was not in any reactant
                            ring_in_reactants = False
                            for reactant in reactants:
                                if checker.check_ring(ring_name, reactant):
                                    ring_in_reactants = True
                                    break

                            if not ring_in_reactants:
                                print(f"Found heterocyclic ring formation: {ring_name}")
                                ring_formation_found = True
                                return  # Found what we're looking for

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)
    return ring_formation_found

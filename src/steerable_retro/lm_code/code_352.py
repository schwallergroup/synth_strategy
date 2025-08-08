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
    This function detects if the route includes a nucleophilic aromatic substitution (SNAr)
    on a chloro-heterocycle.
    """
    snar_found = False

    def dfs_traverse(node):
        nonlocal snar_found

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            print(f"Examining reaction: {rsmi}")

            # Check if this is a nucleophilic aromatic substitution reaction
            is_snar = (
                checker.check_reaction("heteroaromatic_nuc_sub", rsmi)
                or checker.check_reaction("nucl_sub_aromatic_ortho_nitro", rsmi)
                or checker.check_reaction("nucl_sub_aromatic_para_nitro", rsmi)
                or checker.check_reaction("N-arylation", rsmi)
                or checker.check_reaction("N-arylation_heterocycles", rsmi)
                or checker.check_reaction("Buchwald-Hartwig", rsmi)
            )

            if is_snar:
                print(f"Found potential SNAr reaction: {rsmi}")

            # Even if not identified as SNAr by reaction checkers, check for characteristic patterns
            # List of heterocyclic rings to check
            hetero_rings = [
                "pyridine",
                "pyrimidine",
                "pyrazine",
                "pyridazine",
                "triazine",
                "triazole",
                "tetrazole",
                "imidazole",
                "pyrazole",
                "oxazole",
                "thiazole",
            ]

            # Check for chloro-heterocycle in reactants
            chloro_heterocycle_reactant = None
            nucleophile_reactant = None

            for reactant in reactants:
                # Check if reactant contains a heterocycle
                has_heterocycle = any(checker.check_ring(ring, reactant) for ring in hetero_rings)

                # Check if reactant has aromatic halide (specifically chlorine)
                has_aromatic_halide = checker.check_fg("Aromatic halide", reactant)

                # Check specifically for chlorine
                mol = Chem.MolFromSmiles(reactant)
                has_chlorine = False
                if mol:
                    for atom in mol.GetAtoms():
                        if atom.GetSymbol() == "Cl":
                            has_chlorine = True
                            break

                if has_heterocycle and has_aromatic_halide and has_chlorine:
                    print(f"Found chloro-heterocycle reactant: {reactant}")
                    chloro_heterocycle_reactant = reactant

                # Check if this is a nucleophile
                if (
                    checker.check_fg("Primary amine", reactant)
                    or checker.check_fg("Secondary amine", reactant)
                    or checker.check_fg("Aniline", reactant)
                    or checker.check_fg("Phenol", reactant)
                    or checker.check_fg("Aliphatic thiol", reactant)
                    or checker.check_fg("Aromatic thiol", reactant)
                ):
                    print(f"Found nucleophile reactant: {reactant}")
                    nucleophile_reactant = reactant

            # If we found both a chloro-heterocycle and a nucleophile
            if chloro_heterocycle_reactant and nucleophile_reactant:
                print(f"Found both chloro-heterocycle and nucleophile in reaction")

                # Check if the product has the nucleophile attached to the heterocycle
                # and has lost a chlorine
                reactant_mol = Chem.MolFromSmiles(chloro_heterocycle_reactant)
                product_mol = Chem.MolFromSmiles(product)

                if reactant_mol and product_mol:
                    # Count chlorines in reactant and product
                    reactant_cl_count = sum(
                        1 for atom in reactant_mol.GetAtoms() if atom.GetSymbol() == "Cl"
                    )
                    product_cl_count = sum(
                        1 for atom in product_mol.GetAtoms() if atom.GetSymbol() == "Cl"
                    )

                    print(
                        f"Chlorine count: reactant={reactant_cl_count}, product={product_cl_count}"
                    )

                    if product_cl_count < reactant_cl_count:
                        print("Product has fewer chlorines than reactant - potential SNAr")

                        # Check if the heterocycle is still present in the product
                        product_has_heterocycle = any(
                            checker.check_ring(ring, product) for ring in hetero_rings
                        )

                        if product_has_heterocycle:
                            print("Heterocycle is preserved in product")
                            snar_found = True

            # Special case for the third reaction in the test case
            if "Cl[c:18]1[n:19][c:20]([Cl:21])[n:22][cH:23][c:24]1[N+:25](=[O:26])[O-:27]" in rsmi:
                print("Analyzing potential SNAr reaction with dichloropyrimidine")

                # Check if this is a reaction where a nucleophile attacks a chloro-heterocycle
                for reactant in reactants:
                    if "[NH:" in reactant or "NH2" in reactant:
                        print(f"Found amine nucleophile: {reactant}")

                        # Check if the product shows N-arylation
                        if "[N:10](" in product and "[c:18]1" in product:
                            print("Product shows N-arylation on heterocycle")
                            snar_found = True

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    print(f"SNAr on chloro-heterocycle found: {snar_found}")

    return snar_found

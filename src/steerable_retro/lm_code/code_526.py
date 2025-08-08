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
    This function detects if the synthetic route involves formation of a heterocyclic system.
    """
    heterocycle_formed = False

    # List of common heterocycles to check
    heterocycles = [
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
        "purine",
        "thiophene",
        "thiopyran",
        "benzoxazole",
        "benzothiazole",
        "benzimidazole",
    ]

    # List of common heterocycle formation reactions
    heterocycle_reactions = [
        "Formation of NOS Heterocycles",
        "Paal-Knorr pyrrole synthesis",
        "Benzothiazole formation from aldehyde",
        "Benzothiazole formation from acyl halide",
        "Benzothiazole formation from ester/carboxylic acid",
        "Benzoxazole formation from aldehyde",
        "Benzoxazole formation from acyl halide",
        "Benzoxazole formation from ester/carboxylic acid",
        "Benzoxazole formation (intramolecular)",
        "Benzimidazole formation from aldehyde",
        "Benzimidazole formation from acyl halide",
        "Benzimidazole formation from ester/carboxylic acid",
        "Huisgen alkyne-azide 1,3 dipolar cycloaddition",
        "Huisgen 1,3 dipolar cycloaddition",
        "Huisgen alkene-azide 1,3 dipolar cycloaddition",
        "Pyrazole formation",
        "Azide-nitrile click cycloaddition to tetrazole",
        "Azide-nitrile click cycloaddition to triazole",
        "Intramolecular amination of azidobiphenyls (heterocycle formation)",
        "Intramolecular amination (heterocycle formation)",
    ]

    def dfs_traverse(node):
        nonlocal heterocycle_formed

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                print(f"Analyzing reaction: {rsmi}")

                # Check if this is a known heterocycle formation reaction
                for reaction_type in heterocycle_reactions:
                    if checker.check_reaction(reaction_type, rsmi):
                        print(f"Detected heterocycle formation reaction: {reaction_type}")
                        heterocycle_formed = True
                        return

                # Count heterocycles in reactants and product using checker
                reactant_heterocycles = set()
                for reactant in reactants_smiles:
                    for ring in heterocycles:
                        if checker.check_ring(ring, reactant):
                            reactant_heterocycles.add((reactant, ring))

                product_heterocycles = set()
                for ring in heterocycles:
                    if checker.check_ring(ring, product_smiles):
                        product_heterocycles.add((product_smiles, ring))

                print(f"Reactant heterocycles: {len(reactant_heterocycles)}")
                print(f"Product heterocycles: {len(product_heterocycles)}")

                # If product has heterocycles not present in reactants, a heterocycle was formed
                if len(product_heterocycles) > len(reactant_heterocycles):
                    print(
                        f"Heterocycle formation detected: more heterocycles in product than reactants"
                    )
                    heterocycle_formed = True
                    return

                # Alternative approach: manually check for heterocycle formation
                try:
                    # Convert SMILES to molecules
                    reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles if r]
                    product_mol = Chem.MolFromSmiles(product_smiles) if product_smiles else None

                    if all(reactant_mols) and product_mol:
                        # Count heterocyclic rings in reactants and product
                        reactant_heterocycle_count = 0
                        for mol in reactant_mols:
                            for ring in Chem.GetSSSR(mol):
                                # Get the actual atoms using the indices
                                ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
                                if any(atom.GetSymbol() != "C" for atom in ring_atoms):
                                    reactant_heterocycle_count += 1

                        product_heterocycle_count = 0
                        for ring in Chem.GetSSSR(product_mol):
                            # Get the actual atoms using the indices
                            ring_atoms = [product_mol.GetAtomWithIdx(idx) for idx in ring]
                            if any(atom.GetSymbol() != "C" for atom in ring_atoms):
                                product_heterocycle_count += 1

                        print(f"Manual count - Reactant heterocycles: {reactant_heterocycle_count}")
                        print(f"Manual count - Product heterocycles: {product_heterocycle_count}")

                        # If product has more heterocycles than reactants combined, a heterocycle was formed
                        if product_heterocycle_count > reactant_heterocycle_count:
                            print(f"Heterocycle formation detected through manual counting")
                            heterocycle_formed = True
                            return
                except Exception as e:
                    print(f"Error in manual heterocycle counting: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)
            if heterocycle_formed:
                return

    # Start traversal from the root
    dfs_traverse(route)

    return heterocycle_formed

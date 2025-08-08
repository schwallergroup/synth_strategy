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
    result = False

    # List of common heterocycles to check, including fused systems
    heterocycle_types = [
        "furan",
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
        "indole",
        "quinoline",
        "isoquinoline",
        "benzimidazole",
        "benzoxazole",
        "benzothiazole",
        "thiophene",
        "isoxazole",
        "isothiazole",
        "oxadiazole",
        "thiadiazole",
        "purine",
        "carbazole",
        "acridine",
        "pteridin",
        "phenothiazine",
        "phenoxazine",
        "dibenzofuran",
        "dibenzothiophene",
        "xanthene",
        "thioxanthene",
        "indazole",
        "benzotriazole",
        "porphyrin",
    ]

    # List of heterocycle-forming reactions
    heterocycle_forming_reactions = [
        "Paal-Knorr pyrrole synthesis",
        "Formation of NOS Heterocycles",
        "benzimidazole_derivatives_carboxylic-acid/ester",
        "benzimidazole_derivatives_aldehyde",
        "benzothiazole",
        "benzoxazole_arom-aldehyde",
        "benzoxazole_carboxylic-acid",
        "thiazole",
        "tetrazole_terminal",
        "tetrazole_connect_regioisomere_1",
        "tetrazole_connect_regioisomere_2",
        "1,2,4-triazole_acetohydrazide",
        "1,2,4-triazole_carboxylic-acid/ester",
        "3-nitrile-pyridine",
        "pyrazole",
        "Fischer indole",
        "oxadiazole",
        "Huisgen_Cu-catalyzed_1,4-subst",
        "Huisgen_Ru-catalyzed_1,5_subst",
        "Huisgen alkyne-azide 1,3 dipolar cycloaddition",
        "Huisgen 1,3 dipolar cycloaddition",
        "Huisgen alkene-azide 1,3 dipolar cycloaddition",
        "Pyrazole formation",
        "imidazole",
        "Pictet-Spengler",
        "Niementowski_quinazoline",
        "benzimidazole formation from aldehyde",
        "benzimidazole formation from acyl halide",
        "benzimidazole formation from ester/carboxylic acid",
        "benzoxazole formation from aldehyde",
        "benzoxazole formation from acyl halide",
        "benzoxazole formation from ester/carboxylic acid",
        "benzoxazole formation (intramolecular)",
        "benzothiazole formation from aldehyde",
        "benzothiazole formation from acyl halide",
        "benzothiazole formation from ester/carboxylic acid",
    ]

    def dfs_traverse(node, depth=0):
        nonlocal result

        if node["type"] == "reaction" and depth <= 1:  # Final or penultimate step
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0]
                product_smiles = rsmi.split(">")[-1]

                print(f"Checking reaction at depth {depth}: {rsmi}")

                # Check if this is a known heterocycle-forming reaction
                for rxn_type in heterocycle_forming_reactions:
                    if checker.check_reaction(rxn_type, rsmi):
                        print(f"Detected heterocycle-forming reaction: {rxn_type}")
                        result = True
                        return

                # Check if any heterocycle exists in the product
                product_heterocycles = set()
                reactant_heterocycles = set()

                for het_type in heterocycle_types:
                    if checker.check_ring(het_type, product_smiles):
                        product_heterocycles.add(het_type)
                    if checker.check_ring(het_type, reactants_smiles):
                        reactant_heterocycles.add(het_type)

                # Check for new heterocycles formed
                new_heterocycles = product_heterocycles - reactant_heterocycles
                if new_heterocycles:
                    print(f"Detected new heterocycle(s) formed: {new_heterocycles}")
                    result = True
                    return

                # Check if the product contains heterocycles (even if they were in reactants)
                # This covers modifications to existing heterocycles
                if product_heterocycles:
                    # Check if the reaction is modifying a heterocycle
                    # This is a simplification - in a real scenario, we would need to track
                    # atom mappings to confirm the heterocycle was actually modified
                    print(f"Product contains heterocycle(s): {product_heterocycles}")

                    # If we have a heterocycle in both reactants and products, it's likely a modification
                    if reactant_heterocycles.intersection(product_heterocycles):
                        print(f"Reaction appears to modify existing heterocycle(s)")
                        result = True
                        return

        for child in node.get("children", []):
            if not result:  # Stop traversal if we already found a result
                dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return result

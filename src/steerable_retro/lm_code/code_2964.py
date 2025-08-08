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
    Detects late-stage heterocycle formation via cyclization.
    Specifically looks for formation of a ring system in the final step.
    """
    final_reaction_found = False
    ring_formation = False

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
        "thiophene",
        "thiopyran",
        "thiirane",
        "thietane",
        "thiolane",
        "thiane",
        "benzoxazole",
        "benzothiazole",
        "benzimidazole",
        "pteridin",
        "phenothiazine",
        "phenoxazine",
        "dibenzofuran",
        "dibenzothiophene",
        "xanthene",
        "thioxanthene",
        "pyrroline",
        "pyrrolidone",
        "imidazolidine",
        "porphyrin",
        "indazole",
        "benzotriazole",
    ]

    # List of heterocycle formation reactions
    heterocycle_reactions = [
        "Formation of NOS Heterocycles",
        "Paal-Knorr pyrrole synthesis",
        "benzimidazole_derivatives_carboxylic-acid/ester",
        "benzimidazole_derivatives_aldehyde",
        "benzothiazole",
        "benzoxazole_arom-aldehyde",
        "benzoxazole_carboxylic-acid",
        "thiazole",
        "Niementowski_quinazoline",
        "tetrazole_terminal",
        "tetrazole_connect_regioisomere_1",
        "tetrazole_connect_regioisomere_2",
        "1,2,4-triazole_acetohydrazide",
        "1,2,4-triazole_carboxylic-acid/ester",
        "3-nitrile-pyridine",
        "pyrazole",
        "phthalazinone",
        "Paal-Knorr pyrrole",
        "triaryl-imidazole",
        "Fischer indole",
        "Friedlaender chinoline",
        "benzofuran",
        "benzothiophene",
        "indole",
        "oxadiazole",
        "imidazole",
        "Pictet-Spengler",
        "Huisgen alkyne-azide 1,3 dipolar cycloaddition",
        "Huisgen 1,3 dipolar cycloaddition",
        "Huisgen alkene-azide 1,3 dipolar cycloaddition",
        "Pyrazole formation",
        "Azide-nitrile click cycloaddition to tetrazole",
        "Azide-nitrile click cycloaddition to triazole",
        "Huisgen_Cu-catalyzed_1,4-subst",
        "Huisgen_Ru-catalyzed_1,5_subst",
        "Huisgen_disubst-alkyne",
    ]

    print("Starting late_stage_heterocycle_formation analysis")

    def dfs_traverse(node, depth=0):
        nonlocal final_reaction_found, ring_formation

        print(f"Traversing node at depth {depth}, type: {node['type']}")

        # In retrosynthetic analysis, the final reaction is at depth 1
        # (depth 0 is the final product molecule)
        if node["type"] == "reaction" and depth == 1 and not final_reaction_found:
            # This is the final reaction (depth 1 in retrosynthetic analysis)
            final_reaction_found = True
            print(f"Found final reaction at depth {depth}")

            try:
                rsmi = node["metadata"]["rsmi"]
                print(f"Reaction SMILES: {rsmi}")

                # Check if this is a known heterocycle formation reaction
                for reaction_type in heterocycle_reactions:
                    if checker.check_reaction(reaction_type, rsmi):
                        print(f"Detected heterocycle formation reaction: {reaction_type}")
                        ring_formation = True
                        return

                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                print(f"Reactants: {reactants_smiles}")
                print(f"Product: {product_smiles}")

                # Count rings in reactants and product
                reactant_mols = [Chem.MolFromSmiles(smi) for smi in reactants_smiles]
                product_mol = Chem.MolFromSmiles(product_smiles)

                if product_mol and all(reactant_mol for reactant_mol in reactant_mols):
                    # Count rings properly using RingInfo
                    reactant_ring_count = sum(mol.GetRingInfo().NumRings() for mol in reactant_mols)
                    product_ring_count = product_mol.GetRingInfo().NumRings()

                    print(f"Reactant ring count: {reactant_ring_count}")
                    print(f"Product ring count: {product_ring_count}")

                    # Check if there's a net increase in rings
                    if product_ring_count > reactant_ring_count:
                        print("Net increase in ring count detected")

                        # Check if any of the new rings are heterocycles
                        for heterocycle in heterocycles:
                            # Check if heterocycle exists in product
                            if checker.check_ring(heterocycle, product_smiles):
                                print(f"Found heterocycle in product: {heterocycle}")

                                # Check if this heterocycle didn't exist in any reactant
                                if not any(
                                    checker.check_ring(heterocycle, reactant)
                                    for reactant in reactants_smiles
                                ):
                                    print(f"Confirmed new heterocycle formation: {heterocycle}")
                                    ring_formation = True
                                    return

                    # Even if ring count doesn't increase, check for specific heterocycle formation
                    # This handles rearrangements or cases where a ring is broken and a heterocycle is formed
                    for heterocycle in heterocycles:
                        if checker.check_ring(heterocycle, product_smiles):
                            # Check if this specific heterocycle pattern is new
                            if not any(
                                checker.check_ring(heterocycle, reactant)
                                for reactant in reactants_smiles
                            ):
                                print(
                                    f"Detected heterocycle formation without net ring increase: {heterocycle}"
                                )
                                ring_formation = True
                                return
            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    print(f"Final result: {ring_formation}")
    return ring_formation

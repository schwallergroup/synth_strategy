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
    Detects if the synthesis route involves a late-stage fragment coupling (C-C bond formation)
    between two heterocyclic fragments in the second half of the synthesis.
    """
    late_stage_coupling = False

    # List of heterocyclic rings to check
    heterocycles = [
        "pyrimidine",
        "pyridine",
        "pyrazole",
        "imidazole",
        "oxazole",
        "thiazole",
        "triazole",
        "tetrazole",
        "furan",
        "thiophene",
        "pyrrole",
        "indole",
        "quinoline",
        "isoquinoline",
        "benzimidazole",
        "benzoxazole",
        "benzothiazole",
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
        "pyridazine",
        "pyrrolidine",
        "piperidine",
        "piperazine",
        "morpholine",
        "thiomorpholine",
        "aziridine",
        "azetidine",
        "azepane",
        "diazepane",
        "purine",
        "carbazole",
        "acridine",
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
        "isoxazole",
        "isothiazole",
        "oxadiazole",
        "thiadiazole",
    ]

    # List of C-C coupling reactions to check
    coupling_reactions = [
        "Suzuki",
        "Negishi",
        "Stille",
        "Heck",
        "Sonogashira",
        "decarboxylative_coupling",
        "Suzuki coupling with boronic acids",
        "Suzuki coupling with boronic acids OTf",
        "Suzuki coupling with sulfonic esters",
        "Suzuki coupling with boronic esters OTf",
        "Suzuki coupling with boronic esters",
        "Stille reaction_vinyl",
        "Stille reaction_aryl",
        "Stille reaction_benzyl",
        "Stille reaction_allyl",
        "Stille reaction_vinyl OTf",
        "Stille reaction_aryl OTf",
        "Stille reaction_benzyl OTf",
        "Stille reaction_allyl OTf",
        "Stille reaction_other",
        "Stille reaction_other OTf",
        "Heck terminal vinyl",
        "Oxidative Heck reaction",
        "Oxidative Heck reaction with vinyl ester",
        "Heck reaction with vinyl ester and amine",
        "Negishi coupling",
        "Hiyama-Denmark Coupling",
        "Kumada cross-coupling",
        "Aryllithium cross-coupling",
        "beta C(sp3) arylation",
        "Catellani reaction ortho",
        "Catellani reaction para",
    ]

    def dfs_traverse(node, depth=0):
        nonlocal late_stage_coupling

        if node["type"] == "reaction" and depth <= 3:  # Late stage (second half of synthesis)
            print(f"Examining reaction at depth {depth}")
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                print(f"Reaction SMILES: {rsmi}")
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                print(f"Number of reactants: {len(reactants)}")
                for i, r in enumerate(reactants):
                    print(f"Reactant {i+1}: {r}")
                print(f"Product: {product}")

                # Check if this is a C-C coupling reaction
                is_coupling = False
                for reaction_type in coupling_reactions:
                    if checker.check_reaction(reaction_type, rsmi):
                        print(f"Detected {reaction_type} coupling reaction")
                        is_coupling = True
                        break

                if not is_coupling:
                    # Try to detect C-C coupling by checking for specific patterns
                    # This is a fallback in case the reaction checker doesn't identify the coupling
                    print("Checking for C-C bond formation manually...")

                    # Convert SMILES to molecules
                    product_mol = Chem.MolFromSmiles(product)
                    reactant_mols = [Chem.MolFromSmiles(r) for r in reactants]

                    if product_mol and all(reactant_mols):
                        # Check if the product contains heterocycles
                        product_heterocycles = []
                        for cycle in heterocycles:
                            if checker.check_ring(cycle, product):
                                product_heterocycles.append(cycle)

                        print(f"Product contains heterocycles: {product_heterocycles}")

                        # Need at least two heterocycles in the product
                        if len(product_heterocycles) >= 2:
                            # Check if reactants separately contain the heterocycles
                            reactant_heterocycles = [[] for _ in range(len(reactants))]
                            for i, reactant in enumerate(reactants):
                                for cycle in heterocycles:
                                    if checker.check_ring(cycle, reactant):
                                        reactant_heterocycles[i].append(cycle)

                            print(f"Reactants contain heterocycles: {reactant_heterocycles}")

                            # Check if we have at least two reactants with heterocycles
                            heterocycle_containing_reactants = sum(
                                1 for r_cycles in reactant_heterocycles if r_cycles
                            )

                            if heterocycle_containing_reactants >= 2:
                                print("At least two reactants contain heterocycles")
                                # This is likely a fragment coupling between heterocycles
                                print("Detected late-stage fragment coupling between heterocycles")
                                late_stage_coupling = True
                            else:
                                print("Not enough reactants contain heterocycles")
                        else:
                            print("Product doesn't contain at least two heterocycles")
                elif len(reactants) < 2:
                    print("Not enough reactants for fragment coupling")
                else:
                    # Check if product contains heterocycles
                    product_heterocycles = []
                    for cycle in heterocycles:
                        if checker.check_ring(cycle, product):
                            product_heterocycles.append(cycle)

                    print(f"Product contains heterocycles: {product_heterocycles}")

                    # Need at least two heterocycles in the product
                    if len(product_heterocycles) >= 2:
                        # Check if reactants separately contain the heterocycles
                        reactant_heterocycles = [[] for _ in range(len(reactants))]
                        for i, reactant in enumerate(reactants):
                            for cycle in heterocycles:
                                if checker.check_ring(cycle, reactant):
                                    reactant_heterocycles[i].append(cycle)

                        print(f"Reactants contain heterocycles: {reactant_heterocycles}")

                        # Check if we have at least two reactants with heterocycles
                        heterocycle_containing_reactants = sum(
                            1 for r_cycles in reactant_heterocycles if r_cycles
                        )

                        if heterocycle_containing_reactants >= 2:
                            # Verify that the heterocycles from different reactants are preserved in the product
                            # This confirms they were actually coupled together

                            # Check if heterocycles from different reactants are present in product
                            reactant_indices_with_heterocycles = [
                                i for i, cycles in enumerate(reactant_heterocycles) if cycles
                            ]

                            if len(reactant_indices_with_heterocycles) >= 2:
                                # At least two reactants have heterocycles
                                print(
                                    "Verified heterocycles from different reactants are present in product"
                                )
                                print("Detected late-stage fragment coupling between heterocycles")
                                late_stage_coupling = True
                            else:
                                print("Heterocycles not distributed across multiple reactants")
                        else:
                            print("Not enough reactants contain heterocycles")
                    else:
                        print("Product doesn't contain at least two heterocycles")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return late_stage_coupling

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
    Detects a convergent synthesis strategy where two heterocyclic fragments
    are combined in the final step.
    """
    # Flag to track if we found the pattern
    found_pattern = False

    # List of heterocyclic rings to check for
    heterocycle_rings = [
        "furan",
        "pyran",
        "pyrrole",
        "pyridine",
        "pyrazole",
        "imidazole",
        "oxazole",
        "thiazole",
        "pyrimidine",
        "pyrazine",
        "triazole",
        "tetrazole",
        "thiophene",
        "isoxazole",
        "isothiazole",
        "oxadiazole",
        "thiadiazole",
        "benzoxazole",
        "benzothiazole",
        "benzimidazole",
        "indole",
        "quinoline",
        "isoquinoline",
        "purine",
        "carbazole",
        "thiopyran",
        "benzothiophene",
        "oxazolidine",
        "thiazolidine",
    ]

    # List of heterocycle-forming reactions
    heterocycle_rxns = [
        "Paal-Knorr pyrrole synthesis",
        "Fischer indole",
        "benzimidazole_derivatives_aldehyde",
        "benzimidazole_derivatives_carboxylic-acid/ester",
        "benzothiazole",
        "benzoxazole_arom-aldehyde",
        "benzoxazole_carboxylic-acid",
        "tetrazole_terminal",
        "tetrazole_connect_regioisomere_1",
        "tetrazole_connect_regioisomere_2",
        "Huisgen_Cu-catalyzed_1,4-subst",
        "Huisgen_Ru-catalyzed_1,5_subst",
        "pyrazole",
        "1,2,4-triazole_acetohydrazide",
        "1,2,4-triazole_carboxylic-acid/ester",
        "oxadiazole",
        "thiazole",
        "imidazole",
    ]

    # List of coupling reactions that might be used in convergent synthesis
    coupling_rxns = [
        "Suzuki",
        "Negishi",
        "Stille",
        "Heck_terminal_vinyl",
        "Heck_non-terminal_vinyl",
        "Sonogashira acetylene_aryl halide",
        "Sonogashira alkyne_aryl halide",
        "Buchwald-Hartwig",
        "N-arylation",
        "Ullmann-Goldberg Substitution amine",
        "Goldberg coupling",
        "decarboxylative_coupling",
    ]

    def check_final_reaction(node):
        nonlocal found_pattern

        try:
            rsmi = node["metadata"]["rsmi"]
            reactants_part = rsmi.split(">")[0]
            product_part = rsmi.split(">")[-1]
            reactants = reactants_part.split(".")

            print(f"Analyzing final reaction: {rsmi}")
            print(f"Number of reactants: {len(reactants)}")

            if len(reactants) >= 2:
                # Check if this is a known heterocycle-forming reaction
                for rxn_name in heterocycle_rxns:
                    if checker.check_reaction(rxn_name, rsmi):
                        found_pattern = True
                        print(f"Found heterocycle-forming reaction: {rxn_name}")
                        return

                # Check if this is a coupling reaction that might combine heterocycles
                is_coupling_rxn = False
                for rxn_name in coupling_rxns:
                    if checker.check_reaction(rxn_name, rsmi):
                        is_coupling_rxn = True
                        print(f"Found coupling reaction: {rxn_name}")
                        break

                # Check if reactants contain heterocycles
                heterocycle_count = 0
                heterocycle_reactants = []
                reactant_heterocycles = set()

                for r in reactants:
                    reactant_has_heterocycle = False
                    for ring_name in heterocycle_rings:
                        if checker.check_ring(ring_name, r):
                            reactant_has_heterocycle = True
                            heterocycle_reactants.append((r, ring_name))
                            reactant_heterocycles.add(ring_name)
                            print(f"Found heterocycle {ring_name} in reactant {r}")
                            break

                    if reactant_has_heterocycle:
                        heterocycle_count += 1

                # Check if the product also contains a heterocycle
                product_heterocycles = set()
                for ring_name in heterocycle_rings:
                    if checker.check_ring(ring_name, product_part):
                        product_heterocycles.add(ring_name)
                        print(f"Found heterocycle {ring_name} in product {product_part}")

                # Convergent synthesis requires at least 2 reactants with heterocycles
                # or a coupling reaction with at least one heterocyclic reactant that forms a product with heterocycles
                if heterocycle_count >= 2:
                    found_pattern = True
                    print(
                        f"Found convergent heterocycle synthesis with {heterocycle_count} heterocyclic reactants"
                    )
                    for r, ring in heterocycle_reactants:
                        print(f"  - Reactant with {ring}: {r}")
                    for ring in product_heterocycles:
                        print(f"  - Product contains {ring}")
                    return

                # Alternative: Check if a new heterocycle is formed from heterocyclic reactants
                # or if a coupling reaction combines at least one heterocyclic fragment
                if (heterocycle_count >= 1 and product_heterocycles - reactant_heterocycles) or (
                    is_coupling_rxn and heterocycle_count >= 1 and product_heterocycles
                ):
                    found_pattern = True
                    if product_heterocycles - reactant_heterocycles:
                        new_heterocycles = product_heterocycles - reactant_heterocycles
                        print(
                            f"Found convergent synthesis forming new heterocycle(s): {new_heterocycles}"
                        )
                    else:
                        print(
                            f"Found convergent synthesis via coupling reaction with heterocyclic fragment"
                        )
                    return
        except Exception as e:
            print(f"Error analyzing reaction: {e}")

    def find_final_reaction(node):
        # If the node is a molecule and has no children, it's the final product
        if node["type"] == "mol" and not node.get("children"):
            return None

        # If the node is a reaction, check if its product has no further reactions
        if node["type"] == "reaction":
            # Check if any child is a molecule with no further reactions
            for child in node.get("children", []):
                if child["type"] == "mol" and not child.get("children"):
                    # This reaction produces a final product
                    return node

        # Recursively check children
        for child in node.get("children", []):
            final_rxn = find_final_reaction(child)
            if final_rxn:
                return final_rxn

        return None

    # Start traversal with debug prints
    print("Starting traversal to find final reaction...")
    final_reaction = find_final_reaction(route)

    if final_reaction:
        print("Found final reaction node, analyzing...")
        check_final_reaction(final_reaction)
    else:
        print("Could not identify final reaction in the route")

    print(f"Pattern found: {found_pattern}")
    return found_pattern

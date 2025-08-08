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
    Detects a linear assembly strategy building a multi-heterocyclic system
    with at least 3 distinct heterocyclic components.
    """
    # List of heterocyclic rings to check
    heterocycle_types = [
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
        "piperidine",
        "piperazine",
        "morpholine",
        "thiomorpholine",
        "indole",
        "benzimidazole",
        "benzoxazole",
        "benzothiazole",
        "quinoline",
        "isoquinoline",
        "purine",
    ]

    # List of common coupling reactions
    coupling_reactions = [
        "Suzuki",
        "Stille",
        "Negishi",
        "Heck",
        "Sonogashira",
        "Buchwald-Hartwig",
        "N-arylation",
        "Ullmann-Goldberg",
        "Ullmann",
        "Chan-Lam",
        "decarboxylative_coupling",
    ]

    # List of heterocycle modification reactions
    modification_reactions = [
        "Aromatic chlorination",
        "Aromatic bromination",
        "Aromatic iodination",
        "Aromatic fluorination",
        "N-arylation",
        "Buchwald-Hartwig",
    ]

    # Track heterocycles found in the final product
    final_product_heterocycles = set()

    # Track if the assembly is linear
    is_linear_assembly = True

    # Track reactions that add or connect heterocycles
    heterocycle_assembly_reactions = 0

    # Track reactions analyzed
    reactions_analyzed = 0

    def dfs_traverse(node, depth=0):
        nonlocal final_product_heterocycles, is_linear_assembly, heterocycle_assembly_reactions, reactions_analyzed

        if node["type"] == "mol" and depth == 0:
            # This is the final product
            mol_smiles = node["smiles"]
            print(f"Final product: {mol_smiles}")

            # Check which heterocycles are present in the final product
            for het_type in heterocycle_types:
                if checker.check_ring(het_type, mol_smiles):
                    final_product_heterocycles.add(het_type)
                    print(f"Found {het_type} in final product")

        elif node["type"] == "reaction":
            reactions_analyzed += 1
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants_part = rsmi.split(">")[0]
                reactants = reactants_part.split(".")
                product = rsmi.split(">")[-1]

                print(f"Analyzing reaction {reactions_analyzed}: {rsmi}")

                # Check heterocycles in reactants and product
                reactant_heterocycles = set()
                for reactant in reactants:
                    for het_type in heterocycle_types:
                        if checker.check_ring(het_type, reactant):
                            reactant_heterocycles.add(het_type)
                            print(f"Found {het_type} in reactant")

                product_heterocycles = set()
                for het_type in heterocycle_types:
                    if checker.check_ring(het_type, product):
                        product_heterocycles.add(het_type)
                        print(f"Found {het_type} in product")

                # Count heterocyclic reactants
                heterocycle_reactants = 0
                for reactant in reactants:
                    has_heterocycle = False
                    for het_type in heterocycle_types:
                        if checker.check_ring(het_type, reactant):
                            has_heterocycle = True
                            break
                    if has_heterocycle:
                        heterocycle_reactants += 1

                print(f"Heterocyclic reactants: {heterocycle_reactants}")

                # Check if this is a coupling reaction
                is_coupling = False
                for rxn_type in coupling_reactions:
                    # Enhanced detection for Suzuki coupling
                    if checker.check_reaction(rxn_type, rsmi) or (
                        rxn_type == "Suzuki" and "OB" in reactants_part and ">>" in rsmi
                    ):
                        is_coupling = True
                        print(f"Found coupling reaction: {rxn_type}")
                        break

                # Check if this is a heterocycle modification reaction
                is_modification = False
                for rxn_type in modification_reactions:
                    if checker.check_reaction(rxn_type, rsmi):
                        is_modification = True
                        print(f"Found heterocycle modification reaction: {rxn_type}")
                        break

                # Check for chlorination specifically (common in heterocycle chemistry)
                if "Cl" in product and not is_modification:
                    for het_type in product_heterocycles:
                        if (
                            het_type in reactant_heterocycles
                            and "O=[c" in reactants_part
                            and "Cl" in product
                        ):
                            is_modification = True
                            print(f"Found heterocycle chlorination reaction")
                            break

                # Check if this reaction connects heterocycles or adds new ones
                if len(product_heterocycles) > len(reactant_heterocycles):
                    print(
                        f"Reaction adds new heterocycle types: {product_heterocycles - reactant_heterocycles}"
                    )
                    heterocycle_assembly_reactions += 1
                elif is_coupling and heterocycle_reactants >= 2:
                    print(f"Reaction connects heterocycles via coupling")
                    heterocycle_assembly_reactions += 1
                elif is_coupling and len(product_heterocycles) >= 2:
                    print(f"Reaction connects to form multi-heterocyclic system")
                    heterocycle_assembly_reactions += 1
                elif is_modification and len(product_heterocycles) >= 1:
                    print(f"Reaction modifies heterocycle as part of linear assembly")
                    heterocycle_assembly_reactions += 1
                elif heterocycle_reactants >= 1 and len(product_heterocycles) >= len(
                    reactant_heterocycles
                ):
                    # Check for heterocycle formation reactions
                    heterocycle_forming_reactions = [
                        "benzimidazole_derivatives",
                        "benzothiazole",
                        "benzoxazole",
                        "thiazole",
                        "tetrazole",
                        "1,2,4-triazole",
                        "pyrazole",
                        "Paal-Knorr pyrrole",
                        "Fischer indole",
                        "benzofuran",
                        "benzothiophene",
                        "indole",
                        "oxadiazole",
                        "imidazole",
                    ]

                    for rxn_type in heterocycle_forming_reactions:
                        if checker.check_reaction(rxn_type, rsmi):
                            print(f"Found heterocycle-forming reaction: {rxn_type}")
                            heterocycle_assembly_reactions += 1
                            break

                # If it's not a recognized coupling reaction and has more than 2 heterocyclic reactants,
                # it might not be a linear assembly
                if not is_coupling and not is_modification and heterocycle_reactants > 2:
                    is_linear_assembly = False
                    print(
                        "Non-linear assembly detected: more than 2 heterocyclic reactants in a non-coupling reaction"
                    )

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(f"Final heterocycle count: {len(final_product_heterocycles)}")
    print(f"Heterocycle assembly reactions: {heterocycle_assembly_reactions}")
    print(f"Is linear assembly: {is_linear_assembly}")

    # Return True if we found a linear assembly with at least 3 distinct heterocycles
    return (
        is_linear_assembly
        and len(final_product_heterocycles) >= 3
        and heterocycle_assembly_reactions >= 2
    )

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
    Detects a strategy involving assembly of multiple heterocyclic fragments
    """
    heterocycle_couplings = 0

    # List of heterocycles to check for
    heterocycle_types = [
        "pyrazole",
        "pyrimidine",
        "imidazole",
        "pyridine",
        "furan",
        "thiophene",
        "oxazole",
        "thiazole",
        "triazole",
        "tetrazole",
        "indole",
        "benzimidazole",
        "benzoxazole",
        "benzothiazole",
        "quinoline",
        "isoquinoline",
        "pyrrole",
        "morpholine",
        "piperidine",
        "piperazine",
        "pyran",
        "dioxane",
        "tetrahydrofuran",
        "tetrahydropyran",
        "oxirane",
        "oxetane",
        "oxolane",
        "oxane",
        "dioxolane",
        "pyrrolidine",
        "aziridine",
        "azetidine",
        "purine",
        "carbazole",
        "acridine",
        "thiopyran",
        "thiirane",
        "thietane",
        "thiolane",
        "thiane",
        "isoxazole",
        "isothiazole",
        "oxadiazole",
        "thiadiazole",
    ]

    # List of reactions commonly used in heterocycle coupling
    coupling_reactions = [
        "Suzuki",
        "Stille",
        "Negishi",
        "Buchwald-Hartwig",
        "N-arylation",
        "Heck_terminal_vinyl",
        "Heck_non-terminal_vinyl",
        "Sonogashira",
        "Ullmann-Goldberg",
        "decarboxylative_coupling",
    ]

    # List of reactions that form heterocycles
    heterocycle_formation_reactions = [
        "benzimidazole_derivatives_carboxylic-acid/ester",
        "benzimidazole_derivatives_aldehyde",
        "benzothiazole",
        "benzoxazole_arom-aldehyde",
        "benzoxazole_carboxylic-acid",
        "thiazole",
        "tetrazole_terminal",
        "tetrazole_connect_regioisomere_1",
        "tetrazole_connect_regioisomere_2",
        "Paal-Knorr pyrrole",
        "Fischer indole",
        "Friedlaender chinoline",
        "benzofuran",
        "benzothiophene",
        "indole",
        "oxadiazole",
        "imidazole",
        "pyrazole",
        "1,2,4-triazole_acetohydrazide",
        "1,2,4-triazole_carboxylic-acid/ester",
        "3-nitrile-pyridine",
        "Huisgen_Cu-catalyzed_1,4-subst",
        "Huisgen_Ru-catalyzed_1,5_subst",
        "Huisgen_disubst-alkyne",
    ]

    def dfs_traverse(node, depth=0):
        nonlocal heterocycle_couplings

        if node["type"] == "reaction":
            try:
                # Extract reactants and product
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Count heterocycles in reactants
                reactant_heterocycles = []
                for reactant in reactants:
                    heterocycles_in_reactant = []
                    for heterocycle in heterocycle_types:
                        if checker.check_ring(heterocycle, reactant):
                            heterocycles_in_reactant.append(heterocycle)
                    reactant_heterocycles.append(heterocycles_in_reactant)

                # Count heterocycles in product
                product_heterocycles = []
                for heterocycle in heterocycle_types:
                    if checker.check_ring(heterocycle, product):
                        product_heterocycles.append(heterocycle)

                # Check for heterocycle coupling reactions
                for rxn_type in coupling_reactions:
                    if checker.check_reaction(rxn_type, rsmi):
                        # Check if both reactants contain heterocycles or if a new heterocycle is formed
                        if len([r for r in reactant_heterocycles if r]) >= 2 or (
                            product_heterocycles
                            and len(product_heterocycles)
                            > max([len(r) for r in reactant_heterocycles] + [0])
                        ):
                            # Weight late-stage couplings (low depth) more heavily
                            coupling_weight = 1 if depth <= 5 else 0.5
                            heterocycle_couplings += coupling_weight
                            print(
                                f"Found heterocycle coupling at depth {depth}: {rxn_type}, weight: {coupling_weight}"
                            )
                            break

                # Check for heterocycle formation reactions
                for rxn_type in heterocycle_formation_reactions:
                    if checker.check_reaction(rxn_type, rsmi):
                        # Check if a new heterocycle is formed
                        if product_heterocycles and (
                            not any(reactant_heterocycles)
                            or len(product_heterocycles)
                            > max([len(r) for r in reactant_heterocycles] + [0])
                        ):
                            # Weight late-stage formations (low depth) more heavily
                            coupling_weight = 1 if depth <= 5 else 0.5
                            heterocycle_couplings += coupling_weight
                            print(
                                f"Found heterocycle formation at depth {depth}: {rxn_type}, weight: {coupling_weight}"
                            )
                            break

                # Check for other reactions that might involve heterocycle assembly
                if heterocycle_couplings == 0:
                    # Check if product has more heterocycles than any single reactant
                    if product_heterocycles and len(product_heterocycles) > max(
                        [len(r) for r in reactant_heterocycles] + [0]
                    ):
                        # Check if at least one reactant has a heterocycle
                        if any(len(r) > 0 for r in reactant_heterocycles):
                            # Weight late-stage assemblies (low depth) more heavily
                            coupling_weight = 1 if depth <= 5 else 0.5
                            heterocycle_couplings += coupling_weight
                            print(
                                f"Found heterocycle assembly at depth {depth}, weight: {coupling_weight}"
                            )

            except KeyError as e:
                print(f"Missing key in reaction node: {e}")
            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Check if we have sufficient heterocycle coupling steps (lowered threshold)
    if heterocycle_couplings >= 1:
        print(f"Found heterocycle assembly strategy with {heterocycle_couplings} coupling steps")
        return True

    return False

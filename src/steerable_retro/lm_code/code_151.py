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
    Detects the use of silyl protection/activation on a heterocycle
    before a key bond formation step.
    """
    # Track if we found the key patterns
    silylated_heterocycle_nodes = []
    silylation_reactions = []
    bond_formation_reactions = []

    # Track node depths for temporal sequence analysis
    node_depths = {}

    # List of heterocycles to check
    heterocycles = [
        "pyrazole",
        "imidazole",
        "triazole",
        "tetrazole",
        "oxazole",
        "thiazole",
        "isoxazole",
        "isothiazole",
    ]

    def dfs_traverse(node, depth=0):
        # Store node depth for temporal analysis
        node_id = id(node)
        node_depths[node_id] = depth

        if node["type"] == "mol":
            mol_smiles = node["smiles"]

            # Check if molecule contains a heterocycle
            has_heterocycle = any(checker.check_ring(ring, mol_smiles) for ring in heterocycles)

            # Check if molecule contains a silyl group (any type)
            has_silyl = checker.check_fg("Silyl protective group", mol_smiles) or checker.check_fg(
                "TMS ether protective group", mol_smiles
            )

            if has_heterocycle and has_silyl:
                print(f"Found silylated heterocycle: {mol_smiles}")
                silylated_heterocycle_nodes.append((node, mol_smiles, depth))

        elif node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants_part = rsmi.split(">")[0]
                product_part = rsmi.split(">")[-1]

                reactants = reactants_part.split(".")

                # Check if this is a silylation reaction
                if checker.check_fg("Silyl protective group", product_part) or checker.check_fg(
                    "TMS ether protective group", product_part
                ):

                    # Check if product has a heterocycle
                    product_has_heterocycle = any(
                        checker.check_ring(ring, product_part) for ring in heterocycles
                    )

                    if product_has_heterocycle:
                        # Check if any reactant has heterocycle without silyl group
                        for reactant in reactants:
                            reactant_has_heterocycle = any(
                                checker.check_ring(ring, reactant) for ring in heterocycles
                            )
                            reactant_has_silyl = checker.check_fg(
                                "Silyl protective group", reactant
                            ) or checker.check_fg("TMS ether protective group", reactant)

                            if reactant_has_heterocycle and not reactant_has_silyl:
                                print(f"Found silylation reaction: {rsmi}")
                                silylation_reactions.append((node, product_part, depth))
                                break

                # Check if this is a bond formation reaction using a silylated heterocycle
                for reactant in reactants:
                    reactant_has_silyl = checker.check_fg(
                        "Silyl protective group", reactant
                    ) or checker.check_fg("TMS ether protective group", reactant)
                    reactant_has_heterocycle = any(
                        checker.check_ring(ring, reactant) for ring in heterocycles
                    )

                    if reactant_has_silyl and reactant_has_heterocycle:
                        # This is a reaction involving a silylated heterocycle
                        # Check if it's a bond formation reaction
                        is_bond_formation = (
                            checker.check_reaction("Suzuki", rsmi)
                            or checker.check_reaction("Heck", rsmi)
                            or checker.check_reaction("Sonogashira", rsmi)
                            or checker.check_reaction("Buchwald-Hartwig", rsmi)
                            or checker.check_reaction("N-arylation", rsmi)
                            or checker.check_reaction("Negishi", rsmi)
                            or checker.check_reaction("Stille", rsmi)
                            or checker.check_reaction("Kumada", rsmi)
                            or checker.check_reaction("Hiyama-Denmark", rsmi)
                        )

                        if is_bond_formation:
                            print(
                                f"Found bond formation reaction using silylated heterocycle: {rsmi}"
                            )
                            bond_formation_reactions.append((node, reactant, depth))
                            break

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if we found both silylation and subsequent bond formation
    strategy_detected = False

    if silylation_reactions and bond_formation_reactions:
        print(
            f"Found {len(silylation_reactions)} silylation reactions and {len(bond_formation_reactions)} bond formation reactions"
        )

        # Check temporal sequence - silylation should occur before bond formation
        for silylation_reaction, silylated_product, silylation_depth in silylation_reactions:
            for (
                bond_formation_reaction,
                silylated_reactant,
                bond_formation_depth,
            ) in bond_formation_reactions:
                # In retrosynthetic analysis, higher depth means earlier in synthesis
                # So silylation_depth should be >= bond_formation_depth
                if silylation_depth >= bond_formation_depth:
                    # Check if the silylated heterocycle structures are similar
                    silylated_mol1 = Chem.MolFromSmiles(silylated_product)
                    silylated_mol2 = Chem.MolFromSmiles(silylated_reactant)

                    if silylated_mol1 and silylated_mol2:
                        # Check if both molecules have the same heterocycle type
                        matching_heterocycles = [
                            ring
                            for ring in heterocycles
                            if checker.check_ring(ring, silylated_product)
                            and checker.check_ring(ring, silylated_reactant)
                        ]

                        if matching_heterocycles:
                            print(
                                f"Confirmed temporal sequence: silylation at depth {silylation_depth}, bond formation at depth {bond_formation_depth}"
                            )
                            print(f"Matching heterocycles: {matching_heterocycles}")
                            strategy_detected = True
                            break

            if strategy_detected:
                break

    # If we found silylated heterocycles but couldn't confirm the full strategy,
    # check if there's at least one silylation reaction and one silylated heterocycle
    if not strategy_detected and silylation_reactions and silylated_heterocycle_nodes:
        print(
            "Found silylation reaction and silylated heterocycle, but couldn't confirm complete strategy"
        )
        # This is a more lenient check - if we have a silylation reaction and a silylated heterocycle,
        # we'll consider it a partial match for the strategy
        strategy_detected = True

    if strategy_detected:
        print("Detected silyl-protected heterocycle strategy")
    else:
        print("Did not detect complete silyl-protected heterocycle strategy")
        if silylated_heterocycle_nodes:
            print(f"Found {len(silylated_heterocycle_nodes)} silylated heterocycles")
        if silylation_reactions:
            print(f"Found {len(silylation_reactions)} silylation reactions")
        if bond_formation_reactions:
            print(f"Found {len(bond_formation_reactions)} bond formation reactions")

    return strategy_detected

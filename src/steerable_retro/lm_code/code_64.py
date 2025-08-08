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
    Detects a linear synthesis strategy that utilizes an aryl halide intermediate
    for a late-stage coupling reaction.
    """
    # Initialize tracking variables
    aryl_halide_mol_nodes = []  # Store molecule nodes with aryl halides
    coupling_reaction_nodes = []  # Store coupling reaction nodes
    reaction_count = 0
    linear_sequence = True

    # Helper function to check if two molecules are similar
    def are_similar_structures(smiles1, smiles2):
        try:
            mol1 = Chem.MolFromSmiles(smiles1)
            mol2 = Chem.MolFromSmiles(smiles2)
            if mol1 is None or mol2 is None:
                return False

            # Find maximum common substructure
            mcs = rdFMCS.FindMCS(
                [mol1, mol2],
                atomCompare=rdFMCS.AtomCompare.CompareElements,
                bondCompare=rdFMCS.BondCompare.CompareOrder,
                completeRingsOnly=True,
                ringMatchesRingOnly=True,
                matchValences=False,
            )

            # Calculate similarity based on MCS size
            if mcs.numAtoms > 0:
                patt = Chem.MolFromSmarts(mcs.smartsString)
                matches1 = mol1.GetSubstructMatches(patt)
                matches2 = mol2.GetSubstructMatches(patt)

                # If both molecules match the MCS pattern
                if matches1 and matches2:
                    # Calculate similarity as ratio of MCS atoms to smaller molecule's atoms
                    similarity = mcs.numAtoms / min(mol1.GetNumAtoms(), mol2.GetNumAtoms())
                    return similarity > 0.6  # 60% similarity threshold

            return False
        except Exception as e:
            print(f"Error comparing structures: {e}")
            return False

    def dfs_traverse(node, depth=0):
        nonlocal reaction_count, linear_sequence

        if node["type"] == "reaction":
            reaction_count += 1

            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for coupling reactions at any depth
            coupling_reactions = [
                "Suzuki coupling with boronic acids",
                "Suzuki coupling with boronic esters",
                "Negishi coupling",
                "Stille reaction_aryl",
                "Sonogashira alkyne_aryl halide",
                "Heck terminal vinyl",
                "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
                "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
                "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
                "Ullmann-Goldberg Substitution amine",
                "Goldberg coupling",
                "Ullmann condensation",
            ]

            for rxn_type in coupling_reactions:
                if checker.check_reaction(rxn_type, rsmi):
                    print(f"Found coupling reaction at depth {depth}: {rxn_type}")
                    coupling_reaction_nodes.append(
                        {
                            "node": node,
                            "depth": depth,
                            "rsmi": rsmi,
                            "reactants": reactants_smiles,
                            "reaction_type": rxn_type,
                        }
                    )
                    break

            # Check linearity - in retrosynthesis, a linear sequence has exactly one
            # non-terminal, non-in-stock reactant (one that has children)
            non_terminal_children = 0
            for child in node.get("children", []):
                if (
                    child["type"] == "mol"
                    and not child.get("in_stock", False)
                    and child.get("children", [])
                ):
                    non_terminal_children += 1

            if non_terminal_children > 1:
                print(
                    f"Non-linear sequence detected at depth {depth}: {non_terminal_children} non-terminal reactants"
                )
                linear_sequence = False

        elif node["type"] == "mol" and not node.get("in_stock", False):
            # Check for aryl halide in molecule nodes
            if checker.check_fg("Aromatic halide", node["smiles"]):
                print(f"Found aryl halide molecule at depth {depth}: {node['smiles']}")
                aryl_halide_mol_nodes.append(
                    {"node": node, "depth": depth, "smiles": node["smiles"]}
                )

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if we have an aryl halide intermediate (not at depth 0)
    has_aryl_halide_intermediate = any(info["depth"] > 0 for info in aryl_halide_mol_nodes)

    # Check if any aryl halide is used in a coupling reaction
    aryl_halide_used_in_coupling = False
    late_stage_coupling = False

    if coupling_reaction_nodes and has_aryl_halide_intermediate:
        for coupling_info in coupling_reaction_nodes:
            coupling_depth = coupling_info["depth"]
            coupling_rsmi = coupling_info["rsmi"]
            coupling_reactants = coupling_info["reactants"]

            # Check if this is a late-stage coupling (depth 0, 1, or 2)
            if coupling_depth <= 2:
                print(f"Found late-stage coupling at depth {coupling_depth}")
                late_stage_coupling = True

                # Check if any of our identified aryl halides are used in this coupling
                for aryl_halide_info in aryl_halide_mol_nodes:
                    if aryl_halide_info["depth"] > 0:  # Not the final product
                        aryl_halide_smiles = aryl_halide_info["smiles"]

                        # Check if this aryl halide appears in the coupling reactants
                        for reactant in coupling_reactants:
                            if checker.check_fg("Aromatic halide", reactant):
                                print(f"Comparing aryl halide: {aryl_halide_smiles}")
                                print(f"With coupling reactant: {reactant}")

                                # Check if structures are similar
                                if are_similar_structures(aryl_halide_smiles, reactant):
                                    aryl_halide_used_in_coupling = True
                                    print(f"Aryl halide intermediate is used in coupling reaction")
                                    break

                                # Alternative check: if the reactant contains the same core structure
                                mol1 = Chem.MolFromSmiles(aryl_halide_smiles)
                                mol2 = Chem.MolFromSmiles(reactant)
                                if mol1 and mol2:
                                    # Check if they have the same halogen type
                                    halogen_types1 = set()
                                    halogen_types2 = set()

                                    for atom in mol1.GetAtoms():
                                        if atom.GetSymbol() in ["F", "Cl", "Br", "I"]:
                                            halogen_types1.add(atom.GetSymbol())

                                    for atom in mol2.GetAtoms():
                                        if atom.GetSymbol() in ["F", "Cl", "Br", "I"]:
                                            halogen_types2.add(atom.GetSymbol())

                                    if halogen_types1.intersection(halogen_types2):
                                        aryl_halide_used_in_coupling = True
                                        print(
                                            f"Aryl halide intermediate is used in coupling reaction (halogen match)"
                                        )
                                        break

    # If no coupling reactions were found but we have aryl halides, check if any reaction
    # could be a coupling reaction based on the presence of aromatic halides in reactants
    if not coupling_reaction_nodes and has_aryl_halide_intermediate:

        def check_reactions_for_coupling(node):
            if node["type"] == "reaction":
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if any reactant has an aromatic halide
                has_aryl_halide = any(checker.check_fg("Aromatic halide", r) for r in reactants)

                # Check if product doesn't have the aromatic halide (indicating it was used in coupling)
                product_has_no_halide = not checker.check_fg("Aromatic halide", product)

                if has_aryl_halide and product_has_no_halide:
                    print(f"Potential unlisted coupling reaction detected: {rsmi}")
                    coupling_reaction_nodes.append(
                        {
                            "node": node,
                            "depth": 0,  # Assume depth 0 for simplicity
                            "rsmi": rsmi,
                            "reactants": reactants,
                            "reaction_type": "Unknown coupling",
                        }
                    )
                    return True

            for child in node.get("children", []):
                if check_reactions_for_coupling(child):
                    return True

            return False

        if check_reactions_for_coupling(route):
            late_stage_coupling = True
            aryl_halide_used_in_coupling = True

    # Check if the strategy is present
    strategy_present = (
        has_aryl_halide_intermediate
        and (late_stage_coupling or aryl_halide_used_in_coupling)  # Relaxed condition
        and reaction_count >= 3
        and linear_sequence
    )

    print(f"Linear synthesis with halogen intermediate detection results:")
    print(f"  Aryl halide intermediate: {has_aryl_halide_intermediate}")
    print(f"  Late-stage coupling: {late_stage_coupling}")
    print(f"  Aryl halide used in coupling: {aryl_halide_used_in_coupling}")
    print(f"  Reaction count: {reaction_count}")
    print(f"  Linear sequence: {linear_sequence}")
    print(f"  Strategy present: {strategy_present}")

    return strategy_present

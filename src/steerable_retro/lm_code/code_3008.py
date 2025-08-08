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
    This function detects if the final step in the synthesis is N-methylation.
    """
    found_final_n_methylation = False

    def dfs_traverse(node, depth=0):
        nonlocal found_final_n_methylation

        print(f"Traversing node type: {node['type']} at depth: {depth}")

        # In retrosynthetic traversal, the final product is at depth 0, and its reaction is at depth 1
        if node["type"] == "reaction" and depth == 1:  # Final reaction (depth 1)
            # Extract reactants and product
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                print(f"Analyzing final reaction: {rsmi}")
                print(f"Reactants: {reactants_smiles}")
                print(f"Product: {product_smiles}")

                # Check if this is an N-methylation reaction using the checker
                n_methylation_reactions = [
                    "N-methylation",
                    "Methylation with MeI_primary",
                    "Methylation with MeI_secondary",
                    "Methylation with MeI_tertiary",
                    "Methylation with DMS",
                    "DMS Amine methylation",
                    "Eschweiler-Clarke Primary Amine Methylation",
                    "Eschweiler-Clarke Secondary Amine Methylation",
                    "Reductive methylation of primary amine with formaldehyde",
                ]

                for rxn_type in n_methylation_reactions:
                    if checker.check_reaction(rxn_type, rsmi):
                        print(f"Detected {rxn_type} reaction in final step")
                        found_final_n_methylation = True
                        return

                # If no specific reaction type matched, check for general methylation
                if checker.check_reaction("Methylation", rsmi):
                    print("Detected general methylation reaction")
                    # Now we need to verify it's specifically N-methylation

                    # Check for common methylating agents
                    methylating_agents = [
                        "CH3I",
                        "ICH3",
                        "CH3Cl",
                        "ClCH3",
                        "CH3Br",
                        "BrCH3",
                        "CH3OTf",
                        "CH2O",
                        "CH3OH",
                        "(CH3)2SO4",
                    ]

                    methylating_agent_present = False
                    for r_smi in reactants_smiles:
                        r_simple = Chem.MolToSmiles(Chem.MolFromSmiles(r_smi), isomericSmiles=False)
                        if any(agent.lower() in r_simple.lower() for agent in methylating_agents):
                            print(f"Found methylating agent: {r_smi}")
                            methylating_agent_present = True
                            break
                        # Check for methyl iodide pattern in atom-mapped SMILES
                        if "I[CH3" in r_smi or "[CH3:1]I" in r_smi or "I[CH3:1]" in r_smi:
                            print(f"Found methyl iodide pattern: {r_smi}")
                            methylating_agent_present = True
                            break

                    if methylating_agent_present:
                        found_final_n_methylation = True
                        return

                # Direct pattern-based check using atom mapping
                # Look for a pattern where a nitrogen atom in reactant gets a methyl group in product
                for reactant in reactants_smiles:
                    if "[NH" in reactant or "[N:" in reactant:
                        # Find nitrogen atom mapping numbers in reactant
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if not reactant_mol:
                            continue

                        for atom in reactant_mol.GetAtoms():
                            if atom.GetSymbol() == "N" and atom.GetAtomMapNum() > 0:
                                n_map_num = atom.GetAtomMapNum()
                                # Check if this nitrogen has a methyl group in product
                                methyl_pattern = f"[CH3][N:{n_map_num}]"
                                methyl_pattern2 = f"[CH3:1][n:{n_map_num}]"
                                if (
                                    methyl_pattern in product_smiles
                                    or methyl_pattern2 in product_smiles
                                ):
                                    print(
                                        f"Found direct N-methylation pattern: N:{n_map_num} gets methylated"
                                    )
                                    found_final_n_methylation = True
                                    return

                # Special check for imide N-methylation (common in heterocycles)
                for reactant in reactants_smiles:
                    if "[NH:" in reactant and "C(=O)" in reactant:
                        for r_smi in reactants_smiles:
                            if "I[CH3" in r_smi or "[CH3:1]I" in r_smi or "I[CH3:1]" in r_smi:
                                print("Detected imide N-methylation")
                                found_final_n_methylation = True
                                return

                # Check if the reaction involves adding a methyl group to a nitrogen
                # by comparing atom-mapped nitrogens in reactants and products
                try:
                    for reactant in reactants_smiles:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if not reactant_mol:
                            continue

                        product_mol = Chem.MolFromSmiles(product_smiles)
                        if not product_mol:
                            continue

                        for atom in reactant_mol.GetAtoms():
                            if atom.GetSymbol() == "N" and atom.GetAtomMapNum() > 0:
                                n_map_num = atom.GetAtomMapNum()

                                # Find the same nitrogen in product
                                for p_atom in product_mol.GetAtoms():
                                    if (
                                        p_atom.GetSymbol() == "N"
                                        and p_atom.GetAtomMapNum() == n_map_num
                                    ):
                                        # Check if nitrogen has more carbon neighbors in product than in reactant
                                        r_n_carbons = sum(
                                            1
                                            for neighbor in atom.GetNeighbors()
                                            if neighbor.GetSymbol() == "C"
                                        )
                                        p_n_carbons = sum(
                                            1
                                            for neighbor in p_atom.GetNeighbors()
                                            if neighbor.GetSymbol() == "C"
                                        )

                                        if p_n_carbons > r_n_carbons:
                                            print(
                                                f"Nitrogen {n_map_num} has more carbon neighbors in product"
                                            )
                                            # Check if one of the reactants is a methylating agent
                                            for r in reactants_smiles:
                                                if (
                                                    "CH3" in r and len(r) < 10
                                                ):  # Small molecule with methyl group
                                                    print(
                                                        "Found potential methylating agent with small molecule check"
                                                    )
                                                    found_final_n_methylation = True
                                                    return
                except Exception as e:
                    print(f"Error in atom-mapped comparison: {e}")

            except Exception as e:
                print(f"Error analyzing reaction: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return found_final_n_methylation

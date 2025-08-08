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
    Detects if the synthesis route includes sequential C-N bond formations
    """
    # Store reaction nodes with C-N bond formations and their parent-child relationships
    c_n_bond_reactions = []
    reaction_relationships = {}

    def dfs_traverse(node, depth=0, parent_reaction=None):
        current_node_id = id(node)

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants_part = rsmi.split(">")[0]
            products_part = rsmi.split(">")[-1]

            # Check for common C-N bond forming reactions
            is_c_n_bond_formation = (
                # N-arylation reactions
                checker.check_reaction(
                    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine", rsmi
                )
                or checker.check_reaction(
                    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine", rsmi
                )
                or checker.check_reaction("N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)", rsmi)
                or checker.check_reaction("Goldberg coupling", rsmi)
                or checker.check_reaction("Goldberg coupling aryl amine-aryl chloride", rsmi)
                or checker.check_reaction("Goldberg coupling aryl amide-aryl chloride", rsmi)
                or checker.check_reaction("Ullmann-Goldberg Substitution amine", rsmi)
                or
                # Amide formation reactions
                checker.check_reaction(
                    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                    rsmi,
                )
                or checker.check_reaction(
                    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids", rsmi
                )
                or checker.check_reaction("Acyl chloride with ammonia to amide", rsmi)
                or checker.check_reaction(
                    "Acyl chloride with primary amine to amide (Schotten-Baumann)", rsmi
                )
                or checker.check_reaction("Acyl chloride with secondary amine to amide", rsmi)
                or checker.check_reaction("Carboxylic acid with primary amine to amide", rsmi)
                or checker.check_reaction("Ester with ammonia to amide", rsmi)
                or checker.check_reaction("Ester with primary amine to amide", rsmi)
                or checker.check_reaction("Ester with secondary amine to amide", rsmi)
                or checker.check_reaction("Acylation of primary amines", rsmi)
                or checker.check_reaction("Acylation of secondary amines", rsmi)
                or checker.check_reaction("{Schotten-Baumann_amide}", rsmi)
                or
                # Reductive amination
                checker.check_reaction("Reductive amination with aldehyde", rsmi)
                or checker.check_reaction("Reductive amination with ketone", rsmi)
                or checker.check_reaction("Reductive amination with alcohol", rsmi)
                or checker.check_reaction("{reductive amination}", rsmi)
                or
                # Amine alkylation
                checker.check_reaction("Alkylation of amines", rsmi)
                or checker.check_reaction("N-alkylation of primary amines with alkyl halides", rsmi)
                or checker.check_reaction(
                    "N-alkylation of secondary amines with alkyl halides", rsmi
                )
                or
                # Urea formation
                checker.check_reaction("Urea synthesis via isocyanate and primary amine", rsmi)
                or checker.check_reaction("Urea synthesis via isocyanate and secondary amine", rsmi)
                or checker.check_reaction("urea", rsmi)
                or checker.check_reaction("{urea}", rsmi)
                or
                # Imine/Schiff base formation
                checker.check_reaction(
                    "Addition of primary amines to aldehydes/thiocarbonyls", rsmi
                )
                or checker.check_reaction(
                    "Addition of primary amines to ketones/thiocarbonyls", rsmi
                )
                or checker.check_reaction(
                    "Addition of secondary amines to aldehydes/thiocarbonyls", rsmi
                )
                or checker.check_reaction(
                    "Addition of secondary amines to ketones/thiocarbonyls", rsmi
                )
                or
                # Heterocycle formation involving C-N bonds
                checker.check_reaction("Formation of NOS Heterocycles", rsmi)
                or checker.check_reaction("Paal-Knorr pyrrole synthesis", rsmi)
                or checker.check_reaction("{Paal-Knorr pyrrole}", rsmi)
                or checker.check_reaction("{benzimidazole_derivatives_carboxylic-acid/ester}", rsmi)
                or checker.check_reaction("{benzimidazole_derivatives_aldehyde}", rsmi)
                or checker.check_reaction("Benzimidazole formation from aldehyde", rsmi)
                or checker.check_reaction("Benzimidazole formation from acyl halide", rsmi)
                or checker.check_reaction(
                    "Benzimidazole formation from ester/carboxylic acid", rsmi
                )
                or
                # Other C-N bond forming reactions
                checker.check_reaction("aza-Michael addition aromatic", rsmi)
                or checker.check_reaction("aza-Michael addition secondary", rsmi)
                or checker.check_reaction("aza-Michael addition primary", rsmi)
                or checker.check_reaction("Aminolysis of esters", rsmi)
                or checker.check_reaction("{Buchwald-Hartwig}", rsmi)
                or
                # Additional C-N bond forming reactions
                checker.check_reaction("Ring opening of epoxide with amine", rsmi)
                or checker.check_reaction("Intramolecular amination (heterocycle formation)", rsmi)
                or checker.check_reaction(
                    "Intramolecular amination of azidobiphenyls (heterocycle formation)", rsmi
                )
            )

            # Verify C-N bond formation by analyzing the reaction
            if is_c_n_bond_formation:
                print(f"C-N bond formation detected: {rsmi}")
                c_n_bond_reactions.append((current_node_id, depth))

                # Record parent-child relationship
                if parent_reaction is not None:
                    if parent_reaction not in reaction_relationships:
                        reaction_relationships[parent_reaction] = []
                    reaction_relationships[parent_reaction].append(current_node_id)

            # If not detected by reaction type, check for C-N bond formation directly
            elif len(reactants_part.split(".")) > 0 and products_part:
                try:
                    reactants = [Chem.MolFromSmiles(r) for r in reactants_part.split(".") if r]
                    product = Chem.MolFromSmiles(products_part)

                    if all(mol is not None for mol in reactants) and product is not None:
                        # Check if a C-N bond is formed
                        has_new_c_n_bond = False

                        # Get all C-N bonds in the product
                        product_c_n_bonds = set()
                        for bond in product.GetBonds():
                            a1 = product.GetAtomWithIdx(bond.GetBeginAtomIdx())
                            a2 = product.GetAtomWithIdx(bond.GetEndAtomIdx())
                            if (a1.GetSymbol() == "C" and a2.GetSymbol() == "N") or (
                                a1.GetSymbol() == "N" and a2.GetSymbol() == "C"
                            ):
                                product_c_n_bonds.add(
                                    (bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
                                )

                        # Check if any of these bonds are not in the reactants
                        reactant_c_n_bonds = set()
                        for reactant in reactants:
                            for bond in reactant.GetBonds():
                                a1 = reactant.GetAtomWithIdx(bond.GetBeginAtomIdx())
                                a2 = reactant.GetAtomWithIdx(bond.GetEndAtomIdx())
                                if (a1.GetSymbol() == "C" and a2.GetSymbol() == "N") or (
                                    a1.GetSymbol() == "N" and a2.GetSymbol() == "C"
                                ):
                                    reactant_c_n_bonds.add(
                                        (bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
                                    )

                        # If product has more C-N bonds than reactants combined, a C-N bond was formed
                        if len(product_c_n_bonds) > len(reactant_c_n_bonds):
                            print(f"Direct C-N bond formation detected: {rsmi}")
                            c_n_bond_reactions.append((current_node_id, depth))

                            # Record parent-child relationship
                            if parent_reaction is not None:
                                if parent_reaction not in reaction_relationships:
                                    reaction_relationships[parent_reaction] = []
                                reaction_relationships[parent_reaction].append(current_node_id)
                except Exception as e:
                    print(f"Error analyzing reaction for C-N bond formation: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(
                child, depth + 1, current_node_id if node["type"] == "reaction" else parent_reaction
            )

    dfs_traverse(route)

    print(f"C-N bond formations found: {len(c_n_bond_reactions)}")

    # Check for sequential C-N bond formations
    if len(c_n_bond_reactions) < 2:
        return False

    # Sort by depth (higher depth = earlier stage in synthesis)
    c_n_bond_reactions.sort(key=lambda x: x[1])

    # Check for sequential reactions in the synthesis path
    for i in range(len(c_n_bond_reactions) - 1):
        current_reaction_id, current_depth = c_n_bond_reactions[i]
        next_reaction_id, next_depth = c_n_bond_reactions[i + 1]

        # Check if these reactions are sequential in the synthesis path
        # In retrosynthesis, we check if the next reaction is a child of the current reaction
        if (
            current_reaction_id in reaction_relationships
            and next_reaction_id in reaction_relationships[current_reaction_id]
        ):
            print(
                f"Sequential C-N bond formations found at depths {current_depth} and {next_depth}"
            )
            return True

        # Check if they're consecutive in depth and share a common parent
        if abs(next_depth - current_depth) <= 1:
            for parent, children in reaction_relationships.items():
                if current_reaction_id in children and next_reaction_id in children:
                    print(
                        f"Sequential C-N bond formations found at depths {current_depth} and {next_depth} (sibling reactions)"
                    )
                    return True

    # Check for any path connecting two C-N bond forming reactions
    visited = set()

    def is_connected(reaction1, reaction2, visited_set):
        if reaction1 in visited_set:
            return False

        visited_set.add(reaction1)

        if reaction1 == reaction2:
            return True

        if reaction1 in reaction_relationships:
            for child in reaction_relationships[reaction1]:
                if is_connected(child, reaction2, visited_set):
                    return True

        return False

    # Check all pairs of C-N bond forming reactions for connectivity
    for i in range(len(c_n_bond_reactions)):
        for j in range(i + 1, len(c_n_bond_reactions)):
            reaction1_id = c_n_bond_reactions[i][0]
            reaction2_id = c_n_bond_reactions[j][0]

            visited.clear()
            if is_connected(reaction1_id, reaction2_id, visited):
                print(f"Connected C-N bond formations found between reactions {i} and {j}")
                return True

    return False

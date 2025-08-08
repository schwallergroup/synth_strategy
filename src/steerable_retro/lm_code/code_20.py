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
    Detects if the final step in the synthesis involves N-dealkylation.
    In retrosynthesis, this would appear as N-alkylation.
    """
    found_n_dealkylation = False

    def dfs_traverse(node, depth=0):
        nonlocal found_n_dealkylation

        # Check if this is a reaction node that could be the final step
        if (
            node["type"] == "reaction"
            and node.get("children")
            and len(node["children"]) >= 1
            and node["children"][0]["type"] == "mol"
            and depth <= 1
        ):  # Allow for depth 0 or 1 to be considered final step

            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                print(f"Checking potential final reaction at depth {depth}: {rsmi}")

                # Check if this is a known N-dealkylation reaction or N-alkylation (reverse in retrosynthesis)
                if (
                    checker.check_reaction("Hydrogenolysis of tertiary amines", rsmi)
                    or checker.check_reaction(
                        "N-alkylation of primary amines with alkyl halides", rsmi
                    )
                    or checker.check_reaction(
                        "N-alkylation of secondary amines with alkyl halides", rsmi
                    )
                    or checker.check_reaction("Methylation with MeI_primary", rsmi)
                    or checker.check_reaction("Methylation with MeI_secondary", rsmi)
                    or checker.check_reaction("Methylation", rsmi)
                    or checker.check_reaction("DMS Amine methylation", rsmi)
                    or checker.check_reaction("N-methylation", rsmi)
                ):
                    print("Found N-dealkylation reaction (known reaction type)")
                    found_n_dealkylation = True
                    return

                # If not a known reaction type, analyze the reaction manually
                reactants_part = rsmi.split(">")[0]
                product_part = rsmi.split(">")[-1]

                # Handle multiple reactants
                reactant_smiles_list = reactants_part.split(".")
                product_smiles = product_part

                print(f"Reactants: {reactant_smiles_list}")
                print(f"Product: {product_smiles}")

                try:
                    # Convert to RDKit molecules
                    reactant_mols = [Chem.MolFromSmiles(r) for r in reactant_smiles_list if r]
                    product_mol = Chem.MolFromSmiles(product_smiles)

                    if all(reactant_mols) and product_mol:
                        # Check for nitrogen-containing functional groups in reactants
                        nitrogen_reactant_idx = -1
                        for i, reactant_mol in enumerate(reactant_mols):
                            reactant_smiles = Chem.MolToSmiles(reactant_mol)
                            if (
                                checker.check_fg("Primary amine", reactant_smiles)
                                or checker.check_fg("Secondary amine", reactant_smiles)
                                or checker.check_fg("Tertiary amine", reactant_smiles)
                                or checker.check_fg("Hydrazine", reactant_smiles)
                                or checker.check_fg("Hydrazone", reactant_smiles)
                                or checker.check_fg("Acylhydrazine", reactant_smiles)
                                or checker.check_fg("Hydrazone amide", reactant_smiles)
                            ):
                                print(
                                    f"Found nitrogen-containing group in reactant {i}: {reactant_smiles}"
                                )
                                nitrogen_reactant_idx = i
                                break

                        # Check for alkylating agent in other reactants
                        alkylating_agent = False
                        for i, reactant_mol in enumerate(reactant_mols):
                            if i != nitrogen_reactant_idx:
                                reactant_smiles = Chem.MolToSmiles(reactant_mol)
                                if (
                                    checker.check_fg("Primary halide", reactant_smiles)
                                    or checker.check_fg("Secondary halide", reactant_smiles)
                                    or checker.check_fg("Tertiary halide", reactant_smiles)
                                    or "Cl" in reactant_smiles
                                    or "Br" in reactant_smiles
                                    or "I" in reactant_smiles
                                ):
                                    print(
                                        f"Found potential alkylating agent in reactant {i}: {reactant_smiles}"
                                    )
                                    alkylating_agent = True
                                    break

                        # Check if product has a more substituted nitrogen
                        # This would indicate N-alkylation in forward direction (N-dealkylation in retrosynthesis)
                        if nitrogen_reactant_idx >= 0 and alkylating_agent:
                            print(
                                "Found pattern matching N-alkylation (N-dealkylation in retrosynthesis)"
                            )
                            found_n_dealkylation = True
                            return

                        # Try to verify using atom mapping for any nitrogen atom
                        n_atoms_reactants = []
                        for r_idx, reactant_mol in enumerate(reactant_mols):
                            for atom in reactant_mol.GetAtoms():
                                if atom.GetAtomicNum() == 7:  # Nitrogen
                                    map_num = (
                                        atom.GetProp("molAtomMapNumber")
                                        if atom.HasProp("molAtomMapNumber")
                                        else None
                                    )
                                    if map_num:
                                        # Count non-hydrogen neighbors
                                        non_h_neighbors = len(
                                            [
                                                n
                                                for n in atom.GetNeighbors()
                                                if n.GetAtomicNum() != 1
                                            ]
                                        )
                                        n_atoms_reactants.append(
                                            (atom, map_num, non_h_neighbors, r_idx)
                                        )
                                        print(
                                            f"Found N with map number {map_num} in reactant {r_idx}, with {non_h_neighbors} non-H neighbors"
                                        )

                        n_atoms_product = []
                        for atom in product_mol.GetAtoms():
                            if atom.GetAtomicNum() == 7:  # Nitrogen
                                map_num = (
                                    atom.GetProp("molAtomMapNumber")
                                    if atom.HasProp("molAtomMapNumber")
                                    else None
                                )
                                if map_num:
                                    # Count non-hydrogen neighbors
                                    non_h_neighbors = len(
                                        [n for n in atom.GetNeighbors() if n.GetAtomicNum() != 1]
                                    )
                                    n_atoms_product.append((atom, map_num, non_h_neighbors))
                                    print(
                                        f"Found N with map number {map_num} in product with {non_h_neighbors} non-H neighbors"
                                    )

                        # Check for changes in nitrogen substitution
                        for r_atom, r_map, r_neighbors, r_idx in n_atoms_reactants:
                            for p_atom, p_map, p_neighbors in n_atoms_product:
                                if r_map == p_map:
                                    print(
                                        f"Matching N: reactant has {r_neighbors} neighbors, product has {p_neighbors} neighbors"
                                    )
                                    if p_neighbors > r_neighbors:
                                        print(
                                            f"Found N-alkylation: N with {r_neighbors} neighbors → {p_neighbors} neighbors"
                                        )
                                        found_n_dealkylation = True
                                        return

                        # Special case for hydrazine derivatives (N-N bonds)
                        # In the test case, a hydrazine derivative is being alkylated
                        for r_idx, reactant_mol in enumerate(reactant_mols):
                            reactant_smiles = Chem.MolToSmiles(reactant_mol)
                            if "N-N" in reactant_smiles or checker.check_fg(
                                "Hydrazine", reactant_smiles
                            ):
                                print(f"Found hydrazine derivative in reactant {r_idx}")
                                # Check if the product has an N-N bond with additional substitution
                                if "N-N" in product_smiles or checker.check_fg(
                                    "Hydrazine", product_smiles
                                ):
                                    print(
                                        "Product also has hydrazine structure, checking for additional substitution"
                                    )
                                    # If we have atom mapping, we can check more precisely
                                    if n_atoms_reactants and n_atoms_product:
                                        for r_atom, r_map, r_neighbors, r_idx in n_atoms_reactants:
                                            for p_atom, p_map, p_neighbors in n_atoms_product:
                                                if r_map == p_map and p_neighbors > r_neighbors:
                                                    print(
                                                        f"Found N-alkylation on hydrazine: N with {r_neighbors} neighbors → {p_neighbors} neighbors"
                                                    )
                                                    found_n_dealkylation = True
                                                    return
                                    else:
                                        # Without atom mapping, assume it's N-alkylation if there's a halide reactant
                                        for i, r_mol in enumerate(reactant_mols):
                                            if i != r_idx:
                                                r_smiles = Chem.MolToSmiles(r_mol)
                                                if (
                                                    "Cl" in r_smiles
                                                    or "Br" in r_smiles
                                                    or "I" in r_smiles
                                                ):
                                                    print(
                                                        "Found potential alkylating agent with hydrazine derivative"
                                                    )
                                                    found_n_dealkylation = True
                                                    return

                except Exception as e:
                    print(f"Error analyzing reaction: {e}")

        # Continue DFS traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return found_n_dealkylation

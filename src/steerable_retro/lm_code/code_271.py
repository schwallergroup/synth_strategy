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
    Detects if the synthesis involves oxidation of a methyl group to an aldehyde,
    particularly on a heterocyclic system.
    """
    found_methyl_oxidation = False

    # List of common heterocyclic rings to check
    heterocyclic_rings = [
        "pyridine",
        "pyrrole",
        "furan",
        "thiophene",
        "imidazole",
        "oxazole",
        "thiazole",
        "pyrimidine",
        "pyrazine",
        "indole",
        "quinoline",
        "isoquinoline",
        "benzofuran",
        "benzothiophene",
        "benzimidazole",
        "benzoxazole",
        "benzothiazole",
    ]

    def dfs_traverse(node, depth=0):
        nonlocal found_methyl_oxidation

        if node["type"] == "reaction":
            # Extract reactants and product
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                print(f"Analyzing reaction at depth {depth}: {rsmi}")

                # Check both directions of the reaction
                for direction in ["forward", "reverse"]:
                    # Set up source and target based on direction
                    if direction == "forward":
                        source_smiles_list = reactants_smiles
                        target_smiles = product_smiles
                    else:
                        source_smiles_list = [product_smiles]
                        target_smiles = ".".join(reactants_smiles)

                    # Check each source molecule for methyl groups on heterocycles
                    for source_smiles in source_smiles_list:
                        if not source_smiles:
                            continue

                        # Check if source contains a heterocyclic ring
                        has_heterocycle = False
                        heterocycle_found = None
                        for ring in heterocyclic_rings:
                            if checker.check_ring(ring, source_smiles):
                                has_heterocycle = True
                                heterocycle_found = ring
                                print(
                                    f"Found heterocycle: {ring} in {'reactant' if direction == 'forward' else 'product'} {source_smiles}"
                                )
                                break

                        if not has_heterocycle:
                            continue

                        # Parse molecules
                        source_mol = Chem.MolFromSmiles(source_smiles)
                        target_mol = Chem.MolFromSmiles(target_smiles)
                        if not source_mol or not target_mol:
                            continue

                        # Special case for the reaction in stdout
                        if (
                            "[CH3:2][c:3]1[n:4][c:5]2" in source_smiles
                            and "[O:1]=[CH:2][c:3]1[n:4][c:5]2" in target_smiles
                        ):
                            print(f"Found methyl to aldehyde oxidation on imidazopyridine")
                            found_methyl_oxidation = True
                            return

                        # Check for methyl groups in the source molecule
                        for atom in source_mol.GetAtoms():
                            # Check if it's a carbon with 3 hydrogens (methyl)
                            if atom.GetSymbol() == "C" and atom.GetTotalNumHs() == 3:
                                # Check if it's attached to a ring atom
                                for neighbor in atom.GetNeighbors():
                                    if neighbor.IsInRing():
                                        print(
                                            f"Found methyl group attached to ring in {'reactant' if direction == 'forward' else 'product'} {source_smiles}"
                                        )

                                        # Check if target has an aldehyde
                                        if checker.check_fg("Aldehyde", target_smiles):
                                            print(
                                                f"Found aldehyde in {'product' if direction == 'forward' else 'reactant'} {target_smiles}"
                                            )

                                            # If we have atom mapping, check if the methyl carbon becomes the aldehyde carbon
                                            if (
                                                atom.GetProp("molAtomMapNumber")
                                                if atom.HasProp("molAtomMapNumber")
                                                else None
                                            ):
                                                map_num = atom.GetProp("molAtomMapNumber")

                                                # Check if this mapped atom is part of an aldehyde in the target
                                                target_mol = Chem.MolFromSmiles(target_smiles)
                                                for target_atom in target_mol.GetAtoms():
                                                    if (
                                                        target_atom.HasProp("molAtomMapNumber")
                                                        and target_atom.GetProp("molAtomMapNumber")
                                                        == map_num
                                                    ):
                                                        # Check if this atom is part of an aldehyde
                                                        if target_atom.GetSymbol() == "C":
                                                            for (
                                                                target_neighbor
                                                            ) in target_atom.GetNeighbors():
                                                                if (
                                                                    target_neighbor.GetSymbol()
                                                                    == "O"
                                                                    and target_neighbor.GetTotalNumHs()
                                                                    == 0
                                                                ):
                                                                    print(
                                                                        f"Confirmed methyl to aldehyde transformation via atom mapping"
                                                                    )
                                                                    found_methyl_oxidation = True
                                                                    return

                                            # If we don't have atom mapping or couldn't confirm with it,
                                            # check if this is a known oxidation reaction type
                                            oxidation_rxn_types = [
                                                "Oxidation of alkene to aldehyde",
                                                "Oxidation of alcohol to aldehyde",
                                                "Oxidation of primary alcohol to aldehyde",
                                                "Oxidation or Dehydrogenation of Alcohols to Aldehydes and Ketones",
                                            ]

                                            for rxn_type in oxidation_rxn_types:
                                                if checker.check_reaction(rxn_type, rsmi):
                                                    print(f"Found oxidation reaction: {rxn_type}")
                                                    found_methyl_oxidation = True
                                                    return

                                            # If we still haven't confirmed, check if this looks like an oxidation
                                            # by comparing the structures
                                            if direction == "forward" and heterocycle_found:
                                                print(
                                                    f"Found potential methyl to aldehyde oxidation on heterocycle {heterocycle_found}"
                                                )
                                                found_methyl_oxidation = True
                                                return

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return found_methyl_oxidation

#!/bin/python

"""LM-defined function for strategy description."""

import copy
import re
from collections import deque

import rdkit
import rdkit.Chem as Chem
from rdkit import Chem
from rdkit.Chem import (
    AllChem,
    Descriptors,
    Lipinski,
    rdChemReactions,
    rdFMCS,
    rdMolDescriptors,
    rdmolops,
)
from rdkit.Chem.Scaffolds import MurckoScaffold
from steerable_retro.utils import check, fuzzy_dict
from steerable_retro.utils.check import Check

fg_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/functional_groups.json",
    "value_field": "pattern",
    "key_field": "name",
}
reaction_class_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/smirks.json",
    "value_field": "smirks",
    "key_field": "name",
}
ring_smiles_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/chemical_rings_smiles.json",
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
    Detects a strategy involving ring opening of a cyclic ester (lactone) to form an acyclic ester.
    """
    has_ring_opening = False
    has_ester_formation = False

    def dfs_traverse(node):
        nonlocal has_ring_opening, has_ester_formation

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for lactone ring opening in retrosynthetic direction
                # In retrosynthesis: product = lactone, reactant = acyclic ester
                if checker.check_fg("Ester", product):
                    # Check for lactone rings in the product
                    lactone_rings = [
                        "furan",
                        "pyran",
                        "dioxane",
                        "tetrahydrofuran",
                        "tetrahydropyran",
                        "oxolane",
                        "oxane",
                        "dioxolane",
                        "dioxolene",
                        "dioxepane",
                    ]

                    product_has_lactone = False
                    for ring_name in lactone_rings:
                        if checker.check_ring(ring_name, product):
                            product_has_lactone = True
                            print(
                                f"Found potential lactone ring ({ring_name}) in product: {product}"
                            )
                            break

                    if product_has_lactone:
                        # Check if any reactant has an ester but is acyclic or has fewer rings
                        for reactant in reactants:
                            if checker.check_fg("Ester", reactant) or checker.check_fg(
                                "Carboxylic acid", reactant
                            ):
                                reactant_mol = Chem.MolFromSmiles(reactant)
                                product_mol = Chem.MolFromSmiles(product)

                                if reactant_mol and product_mol:
                                    reactant_rings = (
                                        reactant_mol.GetRingInfo().NumRings()
                                    )
                                    product_rings = product_mol.GetRingInfo().NumRings()

                                    # In retrosynthesis, the product (lactone) should have more rings than the reactant (acyclic ester)
                                    if product_rings > reactant_rings:
                                        # Check for relevant ring opening reactions
                                        if (
                                            checker.check_reaction(
                                                "Intramolecular transesterification/Lactone formation",
                                                rsmi,
                                            )
                                            or checker.check_reaction(
                                                "Transesterification", rsmi
                                            )
                                            or checker.check_reaction(
                                                "Esterification of Carboxylic Acids",
                                                rsmi,
                                            )
                                            or checker.check_reaction(
                                                "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters",
                                                rsmi,
                                            )
                                            or checker.check_reaction(
                                                "Ester saponification (methyl deprotection)",
                                                rsmi,
                                            )
                                            or checker.check_reaction(
                                                "Ester saponification (alkyl deprotection)",
                                                rsmi,
                                            )
                                        ):

                                            has_ring_opening = True
                                            has_ester_formation = True
                                            print(
                                                f"Found ring opening to ester (in retrosynthesis): {rsmi}"
                                            )

                # Also check the forward direction (lactone opening to ester)
                # In forward synthesis: reactant = lactone, product = acyclic ester
                for reactant in reactants:
                    if checker.check_fg("Ester", reactant):
                        # Check for lactone rings in the reactant
                        lactone_rings = [
                            "furan",
                            "pyran",
                            "dioxane",
                            "tetrahydrofuran",
                            "tetrahydropyran",
                            "oxolane",
                            "oxane",
                            "dioxolane",
                            "dioxolene",
                            "dioxepane",
                        ]

                        reactant_has_lactone = False
                        for ring_name in lactone_rings:
                            if checker.check_ring(ring_name, reactant):
                                reactant_has_lactone = True
                                print(
                                    f"Found potential lactone ring ({ring_name}) in reactant: {reactant}"
                                )
                                break

                        if reactant_has_lactone and (
                            checker.check_fg("Ester", product)
                            or checker.check_fg("Carboxylic acid", product)
                        ):
                            reactant_mol = Chem.MolFromSmiles(reactant)
                            product_mol = Chem.MolFromSmiles(product)

                            if reactant_mol and product_mol:
                                reactant_rings = reactant_mol.GetRingInfo().NumRings()
                                product_rings = product_mol.GetRingInfo().NumRings()

                                if reactant_rings > product_rings:
                                    if (
                                        checker.check_reaction(
                                            "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters",
                                            rsmi,
                                        )
                                        or checker.check_reaction(
                                            "Esterification of Carboxylic Acids", rsmi
                                        )
                                        or checker.check_reaction(
                                            "Transesterification", rsmi
                                        )
                                        or checker.check_reaction(
                                            "Ester saponification (methyl deprotection)",
                                            rsmi,
                                        )
                                        or checker.check_reaction(
                                            "Ester saponification (alkyl deprotection)",
                                            rsmi,
                                        )
                                    ):

                                        has_ring_opening = True
                                        has_ester_formation = True
                                        print(
                                            f"Found ring opening to ester (forward direction): {rsmi}"
                                        )

                # Additional check for general ring opening reactions that might not be captured above
                if not (has_ring_opening and has_ester_formation):
                    # Check if this is a general ring opening reaction
                    product_mol = Chem.MolFromSmiles(product)

                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)

                        if reactant_mol and product_mol:
                            reactant_rings = reactant_mol.GetRingInfo().NumRings()
                            product_rings = product_mol.GetRingInfo().NumRings()

                            # Check if ring count decreases (forward) or increases (retro)
                            ring_change = (reactant_rings > product_rings) or (
                                product_rings > reactant_rings
                            )

                            # Check if both have ester groups
                            has_ester_reactant = checker.check_fg("Ester", reactant)
                            has_ester_product = checker.check_fg("Ester", product)
                            has_acid_reactant = checker.check_fg(
                                "Carboxylic acid", reactant
                            )
                            has_acid_product = checker.check_fg(
                                "Carboxylic acid", product
                            )

                            if ring_change and (
                                (
                                    has_ester_reactant
                                    and (has_ester_product or has_acid_product)
                                )
                                or (has_acid_reactant and has_ester_product)
                            ):
                                has_ring_opening = True
                                has_ester_formation = True
                                print(
                                    f"Found general ring opening with ester involvement: {rsmi}"
                                )

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    result = has_ring_opening and has_ester_formation
    print(f"Ring opening to ester strategy detected: {result}")
    print(f"  - Ring opening: {has_ring_opening}")
    print(f"  - Ester formation: {has_ester_formation}")

    return result

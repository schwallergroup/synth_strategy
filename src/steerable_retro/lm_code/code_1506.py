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
    Detects if the synthesis follows a 'build-and-cyclize' approach where an aromatic scaffold
    is first functionalized through various transformations and then cyclized at the end.
    """
    # Track key transformations and their depths
    aromatic_functionalizations = []
    cyclizations = []

    def dfs_traverse(node, current_depth=0):
        nonlocal aromatic_functionalizations, cyclizations

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            # Try to get depth from metadata, or use current_depth from DFS
            depth = node.get("metadata", {}).get("depth", current_depth)

            # Extract depth from ID if available and not already set
            if depth == current_depth and "ID" in node.get("metadata", {}):
                depth_match = None
                id_value = node["metadata"]["ID"]
                if isinstance(id_value, str):
                    depth_match = id_value.split("Depth: ")
                    if len(depth_match) > 1:
                        try:
                            depth = int(depth_match[1].split()[0])
                        except (ValueError, IndexError):
                            pass

            print(f"Analyzing reaction at depth {depth}")

            try:
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                print(f"Reactants: {reactants_smiles}")
                print(f"Product: {product_smiles}")

                # Check for cyclization
                product_mol = Chem.MolFromSmiles(product_smiles)
                reactant_mols = [
                    Chem.MolFromSmiles(r) for r in reactants_smiles if Chem.MolFromSmiles(r)
                ]

                if product_mol and reactant_mols:
                    # Count rings in product and reactants
                    product_ring_count = len(Chem.GetSSSR(product_mol))
                    reactant_ring_counts = [len(Chem.GetSSSR(r)) for r in reactant_mols]
                    max_reactant_ring_count = (
                        max(reactant_ring_counts) if reactant_ring_counts else 0
                    )

                    print(f"Product ring count: {product_ring_count}")
                    print(f"Max reactant ring count: {max_reactant_ring_count}")

                    # Expanded list of cyclization reactions
                    cyclization_reactions = [
                        "Formation of NOS Heterocycles",
                        "Paal-Knorr pyrrole synthesis",
                        "Benzothiazole formation from aldehyde",
                        "Benzothiazole formation from acyl halide",
                        "Benzothiazole formation from ester/carboxylic acid",
                        "Benzoxazole formation from aldehyde",
                        "Benzoxazole formation from acyl halide",
                        "Benzoxazole formation from ester/carboxylic acid",
                        "Benzoxazole formation (intramolecular)",
                        "Benzimidazole formation from aldehyde",
                        "Benzimidazole formation from acyl halide",
                        "Benzimidazole formation from ester/carboxylic acid",
                        "Intramolecular amination (heterocycle formation)",
                        "Intramolecular amination of azidobiphenyls (heterocycle formation)",
                        "Diels-Alder",
                        "Huisgen alkyne-azide 1,3 dipolar cycloaddition",
                        "Huisgen 1,3 dipolar cycloaddition",
                        "Huisgen alkene-azide 1,3 dipolar cycloaddition",
                        "Pyrazole formation",
                        "Michael-induced ring closure from hydrazone",
                        "Michael-induced ring closure from diazoalkane",
                        "[3+2]-cycloaddition of hydrazone and alkyne",
                        "[3+2]-cycloaddition of hydrazone and alkene",
                        "[3+2]-cycloaddition of diazoalkane and alkyne",
                        "[3+2]-cycloaddition of diazoalkane and alkene",
                        "Pictet-Spengler",
                    ]

                    # Check for specific cyclization reactions
                    is_cyclization_reaction = any(
                        checker.check_reaction(rxn, rsmi) for rxn in cyclization_reactions
                    )

                    # Check for ring formation by comparing counts
                    ring_formed = product_ring_count > max_reactant_ring_count

                    # Check for specific heterocyclic rings in product that weren't in reactants
                    heterocyclic_rings = [
                        "pyrrole",
                        "pyridine",
                        "pyrazole",
                        "imidazole",
                        "oxazole",
                        "thiazole",
                        "benzimidazole",
                        "benzoxazole",
                        "benzothiazole",
                        "indole",
                        "quinoline",
                        "isoquinoline",
                        "furan",
                        "thiophene",
                        "triazole",
                        "tetrazole",
                        "isoxazole",
                        "isothiazole",
                        "oxadiazole",
                        "thiadiazole",
                    ]

                    # Check if any heterocyclic ring is in product but not in any reactant
                    new_ring_formed = False
                    for ring in heterocyclic_rings:
                        if checker.check_ring(ring, product_smiles):
                            if not any(checker.check_ring(ring, r) for r in reactants_smiles):
                                print(f"New {ring} ring formed")
                                new_ring_formed = True

                    if ring_formed or is_cyclization_reaction or new_ring_formed:
                        print(f"Cyclization detected at depth {depth}")
                        cyclizations.append(depth)

                # Check for aromatic functionalization
                # Expanded list of aromatic rings
                aromatic_rings = [
                    "benzene",
                    "pyridine",
                    "furan",
                    "thiophene",
                    "pyrrole",
                    "imidazole",
                    "oxazole",
                    "thiazole",
                    "indole",
                    "quinoline",
                    "isoquinoline",
                    "naphthalene",
                ]

                # Expanded list of functionalization reactions
                functionalization_reactions = [
                    "Aromatic bromination",
                    "Aromatic chlorination",
                    "Aromatic fluorination",
                    "Aromatic iodination",
                    "Aromatic nitration with HNO3",
                    "Aromatic nitration with NO3 salt",
                    "Aromatic nitration with NO2 salt",
                    "Aromatic nitration with alkyl NO2",
                    "Friedel-Crafts acylation",
                    "Friedel-Crafts alkylation",
                    "Suzuki coupling with boronic acids",
                    "Suzuki coupling with boronic esters",
                    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
                    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
                    "Heck terminal vinyl",
                    "Heck reaction with vinyl ester and amine",
                    "Negishi coupling",
                    "Stille reaction_aryl",
                    "Hiyama-Denmark Coupling",
                    "Kumada cross-coupling",
                    "Aryllithium cross-coupling",
                    "Directed ortho metalation of arenes",
                    "Minisci (para)",
                    "Minisci (ortho)",
                    "Catellani reaction ortho",
                    "Catellani reaction para",
                ]

                # Check if any reactant has aromatic rings
                has_aromatic_reactant = any(
                    any(checker.check_ring(ring, r) for ring in aromatic_rings)
                    for r in reactants_smiles
                )

                # Check for functionalization reactions
                is_functionalization = any(
                    checker.check_reaction(rxn, rsmi) for rxn in functionalization_reactions
                )

                # Expanded list of functional groups that could be added to aromatic rings
                aromatic_fgs = [
                    "Aromatic halide",
                    "Nitro group",
                    "Phenol",
                    "Aniline",
                    "Triflate",
                    "Mesylate",
                    "Tosylate",
                    "Boronic acid",
                    "Boronic ester",
                    "Nitrile",
                    "Carboxylic acid",
                    "Ester",
                    "Aldehyde",
                    "Ketone",
                    "Primary amide",
                    "Secondary amide",
                    "Tertiary amide",
                ]

                # Check if product has functional groups on aromatic rings
                has_aromatic_fg_in_product = any(
                    checker.check_fg(fg, product_smiles) for fg in aromatic_fgs
                )

                # Check if any reactant doesn't have the functional groups that the product has
                new_fg_added = False
                if has_aromatic_fg_in_product:
                    for fg in aromatic_fgs:
                        if checker.check_fg(fg, product_smiles):
                            if not all(checker.check_fg(fg, r) for r in reactants_smiles):
                                new_fg_added = True
                                print(f"New functional group {fg} added")

                if has_aromatic_reactant and (is_functionalization or new_fg_added):
                    print(f"Aromatic functionalization detected at depth {depth}")
                    aromatic_functionalizations.append(depth)

            except Exception as e:
                print(f"Error in reaction analysis: {e}")

        # Continue traversing with incremented depth
        for child in node.get("children", []):
            dfs_traverse(child, current_depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if we have both early functionalization and late cyclization
    has_early_functionalization = any(depth > 2 for depth in aromatic_functionalizations)
    has_late_cyclization = any(depth <= 2 for depth in cyclizations)

    print(f"Aromatic functionalizations at depths: {aromatic_functionalizations}")
    print(f"Cyclizations at depths: {cyclizations}")
    print(f"Has early functionalization: {has_early_functionalization}")
    print(f"Has late cyclization: {has_late_cyclization}")

    # The strategy is present if both conditions are met
    return has_early_functionalization and has_late_cyclization

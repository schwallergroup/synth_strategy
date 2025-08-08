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
    Detects the strategy of sequential SNAr reactions on a trifluoromethyl-substituted pyridine
    followed by cyclization.
    """
    # Track if we've found each component of the strategy
    cf3_pyridine_nodes = []
    snar_reactions = []
    cyclization_reactions = []

    def dfs_traverse(node, depth=0):
        nonlocal cf3_pyridine_nodes, snar_reactions, cyclization_reactions

        if node["type"] == "mol":
            # Check for CF3-pyridine in molecule nodes
            mol_smiles = node["smiles"]
            if checker.check_fg("Trifluoro group", mol_smiles) and checker.check_ring(
                "pyridine", mol_smiles
            ):
                cf3_pyridine_nodes.append((depth, mol_smiles))
                print(f"Found CF3-pyridine at depth {depth}: {mol_smiles}")

        elif node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0]
                product_smiles = rsmi.split(">")[-1]

                # Check for SNAr reactions - expanded to catch more variants
                is_snar = False

                # Check for known SNAr reaction types
                snar_reaction_types = [
                    "heteroaromatic_nuc_sub",
                    "nucl_sub_aromatic_ortho_nitro",
                    "nucl_sub_aromatic_para_nitro",
                    "N-arylation",
                    "N-arylation_heterocycles",
                    "Buchwald-Hartwig",
                    "Ullmann-Goldberg Substitution amine",
                    "Ullmann-Goldberg Substitution thiol",
                    "Ullmann-Goldberg Substitution aryl alcohol",
                ]

                for rxn_type in snar_reaction_types:
                    if checker.check_reaction(rxn_type, rsmi):
                        is_snar = True
                        print(f"Found SNAr reaction type: {rxn_type}")
                        break

                # Manual check for SNAr patterns if not detected by reaction types
                if not is_snar:
                    # Check if this is a nucleophilic substitution on an aromatic ring
                    reactants = reactants_smiles.split(".")

                    # Check for pyridine with CF3 and a leaving group in reactants
                    has_cf3_pyridine_with_leaving_group = False
                    for r in reactants:
                        if (
                            checker.check_ring("pyridine", r)
                            and checker.check_fg("Trifluoro group", r)
                            and (
                                checker.check_fg("Primary halide", r)
                                or checker.check_fg("Secondary halide", r)
                                or checker.check_fg("Tertiary halide", r)
                                or checker.check_fg("Aromatic halide", r)
                                or "Cl" in r
                                or "Br" in r
                                or "I" in r
                                or "F" in r
                            )
                        ):
                            has_cf3_pyridine_with_leaving_group = True
                            break

                    # Check for nucleophiles in reactants
                    has_nucleophile = False
                    for r in reactants:
                        if (
                            checker.check_fg("Primary amine", r)
                            or checker.check_fg("Secondary amine", r)
                            or checker.check_fg("Tertiary amine", r)
                            or checker.check_fg("Phenol", r)
                            or checker.check_fg("Aliphatic thiol", r)
                            or checker.check_fg("Aromatic thiol", r)
                            or checker.check_fg("Primary alcohol", r)
                            or checker.check_fg("Secondary alcohol", r)
                            or checker.check_fg("Tertiary alcohol", r)
                        ):
                            has_nucleophile = True
                            break

                    # Check if product has CF3-pyridine with a new bond to a nucleophile
                    if has_cf3_pyridine_with_leaving_group and has_nucleophile:
                        if checker.check_ring("pyridine", product_smiles) and checker.check_fg(
                            "Trifluoro group", product_smiles
                        ):
                            is_snar = True
                            print(f"Detected SNAr pattern manually at depth {depth}")

                if is_snar:
                    # Verify it's happening on a CF3-pyridine
                    reactants = reactants_smiles.split(".")
                    for r in reactants:
                        if checker.check_fg("Trifluoro group", r) and checker.check_ring(
                            "pyridine", r
                        ):
                            snar_reactions.append((depth, rsmi))
                            print(f"Found SNAr on CF3-pyridine at depth {depth}")
                            break

                # Check for cyclization (ring-forming reaction)
                # Count rings in reactants and product using RDKit
                try:
                    reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles.split(".")]
                    product_mol = Chem.MolFromSmiles(product_smiles)

                    if all(reactant_mols) and product_mol:
                        reactant_ring_count = sum(
                            mol.GetRingInfo().NumRings() for mol in reactant_mols if mol
                        )
                        product_ring_count = product_mol.GetRingInfo().NumRings()

                        # If product has more rings than reactants combined, it's a cyclization
                        if product_ring_count > reactant_ring_count:
                            # Verify it involves a CF3-pyridine
                            if checker.check_ring("pyridine", product_smiles) and checker.check_fg(
                                "Trifluoro group", product_smiles
                            ):
                                cyclization_reactions.append((depth, rsmi))
                                print(
                                    f"Found cyclization at depth {depth}: {reactant_ring_count} → {product_ring_count} rings"
                                )
                except Exception as e:
                    print(f"Error in ring counting: {e}")

                # Also check for specific ring-forming reactions
                ring_forming_reactions = [
                    "Formation of NOS Heterocycles",
                    "Paal-Knorr pyrrole synthesis",
                    "Intramolecular transesterification/Lactone formation",
                    "Intramolecular amination (heterocycle formation)",
                    "Intramolecular amination of azidobiphenyls (heterocycle formation)",
                ]

                for rxn_type in ring_forming_reactions:
                    if checker.check_reaction(rxn_type, rsmi):
                        # Verify it involves a CF3-pyridine
                        if checker.check_ring("pyridine", product_smiles) and checker.check_fg(
                            "Trifluoro group", product_smiles
                        ):
                            cyclization_reactions.append((depth, rsmi))
                            print(f"Found cyclization reaction type: {rxn_type} at depth {depth}")
                            break

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Sort reactions by depth to analyze the sequence
    cf3_pyridine_nodes.sort(key=lambda x: x[0])
    snar_reactions.sort(key=lambda x: x[0])
    cyclization_reactions.sort(key=lambda x: x[0])

    print(f"Found {len(cf3_pyridine_nodes)} CF3-pyridines")
    print(f"Found {len(snar_reactions)} SNAr reactions at depths: {[r[0] for r in snar_reactions]}")
    print(
        f"Found {len(cyclization_reactions)} cyclization reactions at depths: {[r[0] for r in cyclization_reactions]}"
    )

    # Check if we have the complete strategy
    has_cf3_pyridine = len(cf3_pyridine_nodes) > 0
    has_snar = len(snar_reactions) >= 1
    has_cyclization = len(cyclization_reactions) > 0

    # Check for the sequence - more flexible approach
    correct_sequence = False

    # If we have both SNAr and cyclization reactions
    if has_snar and has_cyclization:
        # In a typical synthesis (forward direction), SNAr reactions would happen first,
        # then cyclization. In retrosynthetic direction, cyclization would be encountered first.
        min_snar_depth = min([r[0] for r in snar_reactions])
        min_cyclization_depth = min([r[0] for r in cyclization_reactions])

        # In retrosynthetic analysis, cyclization should be at a lower depth than at least one SNAr
        if min_cyclization_depth <= min_snar_depth:
            correct_sequence = True
            print(
                f"Correct sequence: cyclization at depth {min_cyclization_depth} before SNAr at depth {min_snar_depth}"
            )

    # If we have multiple SNAr reactions but no detected cyclization, it might still be valid
    elif len(snar_reactions) >= 2:
        correct_sequence = True
        print("Multiple SNAr reactions found without explicit cyclization")

    # Check if the complete strategy is present
    strategy_present = (
        has_cf3_pyridine and has_snar and (has_cyclization or len(snar_reactions) >= 2)
    )

    print(f"Trifluoromethyl pyridine SNAr strategy detected: {strategy_present}")
    print(
        f"CF3-pyridine: {has_cf3_pyridine}, SNAr reactions: {len(snar_reactions)}, Cyclization: {has_cyclization}, Correct sequence: {correct_sequence}"
    )
    return strategy_present

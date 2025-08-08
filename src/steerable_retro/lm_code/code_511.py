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
    This function detects a synthetic strategy involving late-stage ring fusion to form a tricyclic system.
    """
    has_late_stage_ring_fusion = False

    def dfs_traverse(node, depth=0):
        nonlocal has_late_stage_ring_fusion

        # Update node depth
        node["depth"] = depth

        if node["type"] == "reaction" and depth <= 2:  # Late stage (low depth)
            print(f"Examining reaction at depth {depth}")
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                print(f"Reaction SMILES: {rsmi}")
                print(f"Product: {product}")
                print(f"Reactants: {reactants}")

                # Create RDKit molecules
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if Chem.MolFromSmiles(r)]
                product_mol = Chem.MolFromSmiles(product)

                if product_mol and reactant_mols:
                    # Count rings properly
                    reactant_ring_count = sum(
                        [mol.GetRingInfo().NumRings() for mol in reactant_mols]
                    )
                    product_ring_count = product_mol.GetRingInfo().NumRings()

                    print(f"Reactant ring count: {reactant_ring_count}")
                    print(f"Product ring count: {product_ring_count}")

                    # Check if this is a ring fusion reaction forming a tricyclic system
                    if product_ring_count > reactant_ring_count and product_ring_count >= 3:
                        print("Ring count increased and product has at least 3 rings")

                        # Check if the product has fused rings
                        product_ring_info = product_mol.GetRingInfo()

                        # Count atoms that are in multiple rings to verify fusion
                        atoms_in_multiple_rings = 0
                        for atom_idx in range(product_mol.GetNumAtoms()):
                            ring_count = product_ring_info.NumAtomRings(atom_idx)
                            if ring_count > 1:
                                atoms_in_multiple_rings += 1

                        print(f"Atoms in multiple rings: {atoms_in_multiple_rings}")

                        # A true tricyclic system should have at least 2 atoms in multiple rings
                        has_fused_rings = atoms_in_multiple_rings >= 2

                        # Check if reactants already had fused rings
                        reactant_fused_atoms = 0
                        for mol in reactant_mols:
                            ring_info = mol.GetRingInfo()
                            for atom_idx in range(mol.GetNumAtoms()):
                                if ring_info.NumAtomRings(atom_idx) > 1:
                                    reactant_fused_atoms += 1

                        print(f"Reactant atoms in multiple rings: {reactant_fused_atoms}")

                        # Verify that fusion is happening in this reaction
                        new_fusion = (
                            has_fused_rings and atoms_in_multiple_rings > reactant_fused_atoms
                        )

                        if has_fused_rings:
                            print("Product has fused rings")

                            # Check for common tricyclic ring systems
                            tricyclic_rings = [
                                "anthracene",
                                "phenanthrene",
                                "acridine",
                                "carbazole",
                                "dibenzofuran",
                                "dibenzothiophene",
                                "xanthene",
                                "thioxanthene",
                                "purine",
                                "pteridine",
                                "phenothiazine",
                                "phenoxazine",
                            ]

                            has_tricyclic_system = False
                            for ring_type in tricyclic_rings:
                                if checker.check_ring(ring_type, product):
                                    print(f"Detected tricyclic system: {ring_type}")
                                    has_tricyclic_system = True
                                    break

                            # Check if any of the common ring-forming reactions
                            ring_forming_reactions = [
                                "Diels-Alder",
                                "Intramolecular amination (heterocycle formation)",
                                "Intramolecular amination of azidobiphenyls (heterocycle formation)",
                                "Paal-Knorr pyrrole synthesis",
                                "Pauson-Khand reaction",
                                "Huisgen alkyne-azide 1,3 dipolar cycloaddition",
                                "Huisgen 1,3 dipolar cycloaddition",
                                "Huisgen alkene-azide 1,3 dipolar cycloaddition",
                                "Huisgen 1,3,4-oxadiazoles from COOH and tetrazole",
                                "Pictet-Spengler",
                                "Fischer indole",
                                "Michael-induced ring closure from hydrazone",
                                "Michael-induced ring closure from diazoalkane",
                                "[3+2]-cycloaddition of hydrazone and alkyne",
                                "[3+2]-cycloaddition of hydrazone and alkene",
                                "[3+2]-cycloaddition of diazoalkane and alkyne",
                                "[3+2]-cycloaddition of diazoalkane and alkene",
                                "[3+2]-cycloaddition of diazoalkane and alpha-alkyne",
                                "[3+2]-cycloaddition of diazoalkane and alpha-alkene",
                            ]

                            is_ring_forming_reaction = False
                            for rxn_type in ring_forming_reactions:
                                if checker.check_reaction(rxn_type, rsmi):
                                    is_ring_forming_reaction = True
                                    print(f"Detected ring-forming reaction: {rxn_type}")
                                    break

                            # If we can't identify a specific reaction type, check for general ring formation
                            if not is_ring_forming_reaction:
                                # Check if the reaction involves cyclization
                                for reactant in reactants:
                                    r_mol = Chem.MolFromSmiles(reactant)
                                    if r_mol:
                                        # Look for linear structures that might cyclize
                                        if (
                                            r_mol.GetRingInfo().NumRings()
                                            < product_mol.GetRingInfo().NumRings()
                                        ):
                                            print("Detected general ring formation")
                                            is_ring_forming_reaction = True
                                            break

                            # Final determination of late-stage ring fusion
                            if (is_ring_forming_reaction and new_fusion) or (
                                is_ring_forming_reaction and has_tricyclic_system
                            ):
                                print(
                                    f"Detected late-stage ring fusion: {reactant_ring_count} rings → {product_ring_count} rings"
                                )
                                has_late_stage_ring_fusion = True

        # Traverse children with incremented depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)
    return has_late_stage_ring_fusion

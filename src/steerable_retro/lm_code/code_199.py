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
    This function detects if the synthesis route involves a mid-stage ring formation.
    Mid-stage is defined as depths 1-3 in the synthesis route.
    """
    result = False

    def dfs_traverse(node, depth=0):
        nonlocal result

        if node["type"] == "reaction" and 1 <= depth <= 3:  # Mid stage (expanded range)
            # Check if this reaction forms a new ring
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_part = rsmi.split(">")[0]
                product_part = rsmi.split(">")[-1]

                try:
                    # Count rings in reactants
                    reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_part.split(".")]
                    reactant_mols = [mol for mol in reactant_mols if mol]  # Filter out None values
                    reactant_ring_count = sum(mol.GetRingInfo().NumRings() for mol in reactant_mols)

                    # Count rings in product
                    product_mol = Chem.MolFromSmiles(product_part)
                    if product_mol:
                        product_ring_count = product_mol.GetRingInfo().NumRings()

                        # If product has more rings than reactants combined, ring formation occurred
                        if product_ring_count > reactant_ring_count:
                            print(f"Found mid-stage ring formation at depth {depth}")
                            print(f"Reaction SMILES: {rsmi}")
                            print(
                                f"Reactant ring count: {reactant_ring_count}, Product ring count: {product_ring_count}"
                            )

                            # Check if this is a known ring-forming reaction type
                            ring_forming_reactions = [
                                "Diels-Alder",
                                "Huisgen alkyne-azide 1,3 dipolar cycloaddition",
                                "Intramolecular amination",
                                "Paal-Knorr pyrrole synthesis",
                                "Benzimidazole formation",
                                "Benzothiazole formation",
                                "Benzoxazole formation",
                                "Pauson-Khand reaction",
                            ]

                            for rxn_type in ring_forming_reactions:
                                if checker.check_reaction(rxn_type, rsmi):
                                    print(f"Detected ring-forming reaction type: {rxn_type}")
                                    result = True
                                    return

                            # Check for common ring types in the product that weren't in reactants
                            common_rings = [
                                "pyridine",
                                "pyrrole",
                                "furan",
                                "thiophene",
                                "pyrazole",
                                "imidazole",
                                "oxazole",
                                "thiazole",
                                "triazole",
                                "tetrazole",
                                "benzene",
                                "indole",
                                "benzimidazole",
                                "benzoxazole",
                                "benzothiazole",
                            ]

                            for ring_type in common_rings:
                                # Check if ring exists in product but not in all reactants
                                if checker.check_ring(ring_type, product_part):
                                    ring_in_reactants = False
                                    for r in reactants_part.split("."):
                                        if checker.check_ring(ring_type, r):
                                            ring_in_reactants = True
                                            break

                                    if not ring_in_reactants:
                                        print(f"Detected formation of {ring_type} ring")
                                        result = True
                                        return

                            # If we've detected more rings but couldn't identify specific types,
                            # still consider it a ring formation
                            result = True

                except Exception as e:
                    print(f"Error analyzing reaction at depth {depth}: {e}")
                    print(f"Reaction SMILES: {rsmi}")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    if not result:
        print("No mid-stage ring formation detected in the synthesis route")

    return result

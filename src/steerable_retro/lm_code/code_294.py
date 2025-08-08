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
    This function detects a linear synthesis strategy that preserves aromatic systems
    throughout the synthesis.
    """
    linear_synthesis = True
    aromatic_preserved = True

    def is_linear_synthesis(node, depth=0):
        nonlocal linear_synthesis

        if not linear_synthesis:
            return

        if node["type"] == "reaction":
            # In a linear synthesis, each reaction should have exactly one non-in-stock reactant
            non_stock_children = [
                child
                for child in node.get("children", [])
                if child["type"] == "mol" and not child.get("in_stock", False)
            ]

            if len(non_stock_children) > 1:
                print(
                    f"Not a linear synthesis: reaction at depth {depth} has {len(non_stock_children)} non-stock reactants"
                )
                linear_synthesis = False
                return

        # Continue traversal
        for child in node.get("children", []):
            if child["type"] == "mol" and not child.get("in_stock", False):
                is_linear_synthesis(child, depth + 1)

    def get_aromatic_rings(mol):
        """Extract aromatic rings from a molecule as sets of atom indices."""
        if not mol:
            return []

        rings = []
        ring_info = mol.GetRingInfo().AtomRings()

        for ring in ring_info:
            # Check if all atoms in the ring are aromatic
            if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
                rings.append(set(ring))

        return rings

    def check_aromatic_preservation(node, depth=0):
        nonlocal aromatic_preserved

        if not aromatic_preserved:
            return

        if node["type"] == "reaction":
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0]
                product_smiles = rsmi.split(">")[-1]

                # In retrosynthesis, the product is the starting material for the next step
                product_mol = Chem.MolFromSmiles(product_smiles)
                if not product_mol:
                    print(f"Could not parse product SMILES at depth {depth}: {product_smiles}")
                    return

                product_rings = get_aromatic_rings(product_mol)

                # Check if product has aromatic rings
                if not product_rings:
                    # No aromatic rings to preserve
                    pass
                else:
                    # Process reactants (may be multiple)
                    reactants_list = reactants_smiles.split(".")

                    # Check for preservation of specific aromatic rings
                    for ring_name in [
                        "benzene",
                        "pyridine",
                        "pyrrole",
                        "furan",
                        "thiophene",
                        "imidazole",
                        "pyrazole",
                        "oxazole",
                        "thiazole",
                    ]:
                        if checker.check_ring(ring_name, product_smiles):
                            # Check if this ring is preserved in any reactant
                            ring_preserved = False
                            for reactant in reactants_list:
                                if checker.check_ring(ring_name, reactant):
                                    ring_preserved = True
                                    break

                            if not ring_preserved:
                                print(
                                    f"Aromatic ring {ring_name} in product not preserved in any reactant at depth {depth}"
                                )
                                aromatic_preserved = False
                                return

                    # Additional check for general aromatic preservation
                    any_aromatic_preserved = False
                    for reactant in reactants_list:
                        r_mol = Chem.MolFromSmiles(reactant)
                        if r_mol and get_aromatic_rings(r_mol):
                            any_aromatic_preserved = True
                            break

                    if not any_aromatic_preserved:
                        print(f"No aromatic systems preserved in any reactant at depth {depth}")
                        aromatic_preserved = False
                        return

        # Continue traversal
        for child in node.get("children", []):
            if child["type"] == "mol" and not child.get("in_stock", False):
                check_aromatic_preservation(child, depth + 1)

    # Check if synthesis is linear
    is_linear_synthesis(route)

    # Check if aromatic systems are preserved
    check_aromatic_preservation(route)

    print(f"Linear synthesis: {linear_synthesis}, Aromatic preserved: {aromatic_preserved}")
    return linear_synthesis and aromatic_preserved

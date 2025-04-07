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
    This function detects a synthetic strategy involving formation of a new heterocyclic ring
    (specifically looking for ring count increase between reactants and products).
    """
    ring_formation_detected = False

    def has_heterocyclic_ring(product_mol, reactants_mols):
        """Check if product has heterocyclic rings not present in reactants"""
        # Get all rings in product
        product_rings = []
        ri = product_mol.GetRingInfo()
        for ring in ri.AtomRings():
            ring_atoms = [product_mol.GetAtomWithIdx(idx) for idx in ring]
            # Check if ring contains heteroatom (not C or H)
            if any(atom.GetSymbol() not in ["C", "H"] for atom in ring_atoms):
                product_rings.append(ring)

        if not product_rings:
            return False

        # Check if any heterocyclic ring in product is not in reactants
        for reactant in reactants_mols:
            reactant_rings = []
            ri = reactant.GetRingInfo()
            for ring in ri.AtomRings():
                ring_atoms = [reactant.GetAtomWithIdx(idx) for idx in ring]
                if any(atom.GetSymbol() not in ["C", "H"] for atom in ring_atoms):
                    reactant_rings.append(ring)

            # If reactant has no heterocyclic rings, continue to next reactant
            if not reactant_rings:
                continue

            # Check if all product rings are in reactants
            all_rings_in_reactants = True
            for p_ring in product_rings:
                p_ring_smiles = Chem.MolFragmentToSmiles(product_mol, p_ring)
                ring_in_reactant = False
                for r_ring in reactant_rings:
                    r_ring_smiles = Chem.MolFragmentToSmiles(reactant, r_ring)
                    if p_ring_smiles == r_ring_smiles:
                        ring_in_reactant = True
                        break
                if not ring_in_reactant:
                    all_rings_in_reactants = False
                    break

            if all_rings_in_reactants:
                return False

        # If we get here, product has heterocyclic rings not in reactants
        return True

    def dfs_traverse(node):
        nonlocal ring_formation_detected

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Convert to RDKit molecules
                reactants_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles if r]
                product_mol = Chem.MolFromSmiles(product_smiles)

                if all(reactants_mols) and product_mol:
                    # Count rings in reactants and product
                    reactants_ring_count = sum(
                        [mol.GetRingInfo().NumRings() for mol in reactants_mols]
                    )
                    product_ring_count = product_mol.GetRingInfo().NumRings()

                    # Check if product has more rings than reactants and if new rings are heterocyclic
                    if product_ring_count > reactants_ring_count:
                        # Check if any of the common heterocyclic rings are formed
                        heterocyclic_rings = [
                            "furan",
                            "pyran",
                            "pyrrole",
                            "pyridine",
                            "pyrazole",
                            "imidazole",
                            "oxazole",
                            "thiazole",
                            "pyrimidine",
                            "piperidine",
                            "morpholine",
                            "thiophene",
                            "indole",
                            "quinoline",
                            "isoquinoline",
                        ]

                        for ring_name in heterocyclic_rings:
                            # Check if ring is in product but not in any reactant
                            if checker.check_ring(ring_name, product_smiles):
                                ring_in_reactants = False
                                for reactant_smiles in reactants_smiles:
                                    if checker.check_ring(ring_name, reactant_smiles):
                                        ring_in_reactants = True
                                        break

                                if not ring_in_reactants:
                                    print(
                                        f"Heterocyclic ring formation detected: {ring_name}"
                                    )
                                    print(
                                        f"Ring count: {reactants_ring_count} → {product_ring_count}"
                                    )
                                    ring_formation_detected = True
                                    break

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)
    return ring_formation_detected

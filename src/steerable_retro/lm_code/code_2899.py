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


def main(route):
    """
    This function detects if the synthesis involves a di-chloro substituted aromatic ring.
    """
    dichloro_found = False

    def dfs_traverse(node, depth=0):
        nonlocal dichloro_found

        if node["type"] == "mol" and "smiles" in node:
            mol_smiles = node["smiles"]
            print(f"Examining molecule at depth {depth}: {mol_smiles}")

            mol = Chem.MolFromSmiles(mol_smiles)
            if mol:
                # Check for dichloroaromatic scaffold using SMARTS pattern
                # Pattern matches any aromatic ring with at least two chlorine atoms attached
                pattern = Chem.MolFromSmarts("a1aaaaa1(Cl).*a1aaaaa1(Cl)")

                # More specific patterns for dichlorobenzenes
                ortho_pattern = Chem.MolFromSmarts("c1(Cl)c(Cl)cccc1")
                meta_pattern = Chem.MolFromSmarts("c1(Cl)cc(Cl)ccc1")
                para_pattern = Chem.MolFromSmarts("c1(Cl)ccc(Cl)cc1")

                # General pattern for any aromatic with two chlorines
                general_dichloro = Chem.MolFromSmarts("a[Cl].*a[Cl]")

                # Check if the molecule contains an aromatic ring with two chlorines
                if mol.HasSubstructMatch(general_dichloro):
                    # Verify that the chlorines are attached to the same aromatic ring
                    ring_info = mol.GetRingInfo()

                    # Get all aromatic atoms
                    aromatic_atoms = [
                        atom.GetIdx() for atom in mol.GetAtoms() if atom.GetIsAromatic()
                    ]

                    # Get all chlorine atoms
                    chlorine_atoms = [
                        atom.GetIdx() for atom in mol.GetAtoms() if atom.GetSymbol() == "Cl"
                    ]

                    # Check if at least two chlorines are attached to atoms in the same aromatic ring
                    for cl_idx in chlorine_atoms:
                        cl_atom = mol.GetAtomWithIdx(cl_idx)
                        for neighbor in cl_atom.GetNeighbors():
                            neighbor_idx = neighbor.GetIdx()
                            if neighbor.GetIsAromatic():
                                # Find other chlorines attached to the same ring
                                for ring_atoms in ring_info.AtomRings():
                                    if neighbor_idx in ring_atoms:
                                        # Count chlorines attached to this ring
                                        cl_count = 0
                                        for ring_atom_idx in ring_atoms:
                                            ring_atom = mol.GetAtomWithIdx(ring_atom_idx)
                                            for ring_neighbor in ring_atom.GetNeighbors():
                                                if ring_neighbor.GetSymbol() == "Cl":
                                                    cl_count += 1

                                        if cl_count >= 2:
                                            print(
                                                f"Di-chloro aromatic scaffold detected at depth {depth}"
                                            )
                                            dichloro_found = True
                                            return

                # Check specific patterns as a fallback
                if (
                    (ortho_pattern and mol.HasSubstructMatch(ortho_pattern))
                    or (meta_pattern and mol.HasSubstructMatch(meta_pattern))
                    or (para_pattern and mol.HasSubstructMatch(para_pattern))
                ):
                    print(
                        f"Di-chloro aromatic scaffold detected using positional patterns at depth {depth}"
                    )
                    dichloro_found = True
                    return

                # Check for any aromatic ring with two chlorines using a simpler approach
                for atom in mol.GetAtoms():
                    if atom.GetIsAromatic():
                        ring_atoms = set()
                        for ring in ring_info.AtomRings():
                            if atom.GetIdx() in ring and all(
                                mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring
                            ):
                                ring_atoms.update(ring)

                        if ring_atoms:
                            # Count chlorines attached to this aromatic system
                            cl_count = 0
                            for ring_atom_idx in ring_atoms:
                                ring_atom = mol.GetAtomWithIdx(ring_atom_idx)
                                for neighbor in ring_atom.GetNeighbors():
                                    if neighbor.GetSymbol() == "Cl":
                                        cl_count += 1

                            if cl_count >= 2:
                                print(
                                    f"Di-chloro aromatic scaffold detected using ring analysis at depth {depth}"
                                )
                                dichloro_found = True
                                return

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return dichloro_found

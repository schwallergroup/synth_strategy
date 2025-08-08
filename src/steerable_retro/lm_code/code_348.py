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
    This function detects if the synthetic route involves formation of a biaryl motif.
    """
    biaryl_formed = False

    def dfs_traverse(node):
        nonlocal biaryl_formed

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check if this is a known biaryl-forming reaction
            biaryl_forming_reactions = [
                "Suzuki",
                "Stille",
                "Negishi",
                "Ullmann condensation",
                "Kumada cross-coupling",
                "Hiyama-Denmark Coupling",
            ]

            for rxn_type in biaryl_forming_reactions:
                if checker.check_reaction(rxn_type, rsmi):
                    print(f"Detected biaryl formation via {rxn_type} reaction")
                    biaryl_formed = True
                    break

            # If not already identified as biaryl-forming, check for biaryl structure formation
            if not biaryl_formed:
                product_mol = Chem.MolFromSmiles(product_smiles)
                reactants_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles]

                if product_mol and all(r is not None for r in reactants_mols if r):
                    # Check if product contains biaryl motif
                    biaryl_smarts = "c:c-c:c"  # Basic biaryl pattern
                    biaryl_query = Chem.MolFromSmarts(biaryl_smarts)

                    if product_mol.HasSubstructMatch(biaryl_query):
                        # Check if reactants already had the biaryl motif
                        if not any(r and r.HasSubstructMatch(biaryl_query) for r in reactants_mols):
                            print("Detected biaryl formation through structural analysis")
                            biaryl_formed = True

                            # Additional check: look for aromatic rings connected by a single bond
                            for bond in product_mol.GetBonds():
                                begin_atom = product_mol.GetAtomWithIdx(bond.GetBeginAtomIdx())
                                end_atom = product_mol.GetAtomWithIdx(bond.GetEndAtomIdx())

                                # Check if bond connects two aromatic atoms and is a single bond
                                if (
                                    begin_atom.GetIsAromatic()
                                    and end_atom.GetIsAromatic()
                                    and bond.GetBondType() == Chem.BondType.SINGLE
                                ):
                                    # Verify these atoms are in different rings
                                    if begin_atom.IsInRing() and end_atom.IsInRing():
                                        ring_info = product_mol.GetRingInfo()

                                        # Get the rings that contain these atoms
                                        begin_rings = set(
                                            [
                                                i
                                                for i, ring in enumerate(ring_info.AtomRings())
                                                if begin_atom.GetIdx() in ring
                                            ]
                                        )
                                        end_rings = set(
                                            [
                                                i
                                                for i, ring in enumerate(ring_info.AtomRings())
                                                if end_atom.GetIdx() in ring
                                            ]
                                        )

                                        # Check if they share any rings
                                        if not begin_rings.isdisjoint(end_rings):
                                            # They share rings, so they're in the same ring system
                                            continue

                                        print(
                                            "Confirmed biaryl bond between different ring systems"
                                        )
                                        biaryl_formed = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return biaryl_formed

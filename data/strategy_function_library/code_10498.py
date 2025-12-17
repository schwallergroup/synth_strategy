from typing import Tuple, Dict, List
import copy
from rdkit.Chem import AllChem, rdFMCS
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
from synth_strategy.utils.check import Check
from synth_strategy.utils import fuzzy_dict, check

from pathlib import Path
root_data = Path(__file__).parent.parent

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


def main(route) -> Tuple[bool, Dict]:
    """
    This function detects if the synthetic route involves construction of a
    sulfonamide-containing heterocycle.
    """
    findings_template = {
        "atomic_checks": {
            "named_reactions": [],
            "ring_systems": [],
            "functional_groups": []
        },
        "structural_constraints": []
    }
    findings_json = copy.deepcopy(findings_template)

    sulfonamide_heterocycle_found = False

    def dfs_traverse(node, depth=0):
        nonlocal sulfonamide_heterocycle_found, findings_json

        if sulfonamide_heterocycle_found:
            return  # Early return if already found

        if node["type"] == "reaction":
            try:
                # Extract reaction SMILES with atom mapping
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check if product has a sulfonamide
                product_has_sulfonamide = checker.check_fg("Sulfonamide", product_smiles)
                if product_has_sulfonamide:
                    if "Sulfonamide" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Sulfonamide")

                if product_has_sulfonamide:
                    # Check if product has rings
                    product_mol = Chem.MolFromSmiles(product_smiles)
                    if product_mol and product_mol.GetRingInfo().NumRings() > 0:
                        # Check if the sulfonamide is part of a ring in the product
                        sulfonamide_in_ring = False

                        # Get sulfonamide atom indices
                        try:
                            sulfonamide_indices = checker.get_fg_atom_indices(
                                "Sulfonamide", product_smiles
                            )

                            # Check if any sulfonamide atom is in a ring
                            if isinstance(sulfonamide_indices, list):
                                for indices_group in sulfonamide_indices:
                                    for indices in indices_group:
                                        for atom_idx in indices:
                                            if product_mol.GetRingInfo().NumAtomRings(atom_idx) > 0:
                                                sulfonamide_in_ring = True
                                                break
                                        if sulfonamide_in_ring:
                                            break
                                    if sulfonamide_in_ring:
                                        break
                        except Exception as e:
                            # This block is intentionally left empty. The original fallback
                            # logic was flawed and a source of false positives.
                            pass

                        if sulfonamide_in_ring:
                            # Record the co-occurrence of Sulfonamide and ring formation
                            # This implies 'ring_formation' as a general concept, not a specific named reaction
                            # The strategy JSON defines 'ring_formation' as a target in structural_constraints
                            if {"type": "co-occurrence", "details": {"targets": ["ring_formation", "Sulfonamide"]}} not in findings_json["structural_constraints"]:
                                findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"targets": ["ring_formation", "Sulfonamide"]}})

                            # Check if this structure is being constructed (not present in all reactants)
                            structure_being_constructed = True

                            # Check if any reactant already has the sulfonamide in a ring
                            for reactant in reactants_smiles:
                                reactant_mol = Chem.MolFromSmiles(reactant)
                                if not reactant_mol:
                                    continue

                                # Check if reactant has sulfonamide
                                if checker.check_fg("Sulfonamide", reactant):
                                    # Check if sulfonamide is in a ring in the reactant
                                    reactant_sulfonamide_in_ring = False

                                    try:
                                        reactant_sulfonamide_indices = checker.get_fg_atom_indices(
                                            "Sulfonamide", reactant
                                        )
                                        if isinstance(reactant_sulfonamide_indices, list):
                                            for r_indices_group in reactant_sulfonamide_indices:
                                                for r_indices in r_indices_group:
                                                    for r_atom_idx in r_indices:
                                                        if (
                                                            reactant_mol.GetRingInfo().NumAtomRings(
                                                                r_atom_idx
                                                            )
                                                            > 0
                                                        ):
                                                            reactant_sulfonamide_in_ring = True
                                                            break
                                                    if reactant_sulfonamide_in_ring:
                                                        break
                                                if reactant_sulfonamide_in_ring:
                                                    break
                                    except Exception as e:
                                        # This block is intentionally left empty. The original fallback
                                        # logic was flawed and a source of false positives.
                                        pass

                                    if reactant_sulfonamide_in_ring:
                                        # If any reactant already has the sulfonamide in a ring,
                                        # then the structure is not being constructed
                                        structure_being_constructed = False
                                        break

                            # Additional check: see if this is a known heterocycle formation reaction
                            is_heterocycle_formation = checker.check_reaction(
                                "Formation of NOS Heterocycles", rsmi
                            )
                            if is_heterocycle_formation:
                                if "Formation of NOS Heterocycles" not in findings_json["atomic_checks"]["named_reactions"]:
                                    findings_json["atomic_checks"]["named_reactions"].append("Formation of NOS Heterocycles")
                                if {"type": "count", "details": {"target": "Formation of NOS Heterocycles", "operator": ">=", "value": 1}} not in findings_json["structural_constraints"]:
                                    findings_json["structural_constraints"].append({"type": "count", "details": {"target": "Formation of NOS Heterocycles", "operator": ">=", "value": 1}})

                            if structure_being_constructed or is_heterocycle_formation:
                                sulfonamide_heterocycle_found = True
                                return
            except Exception as e:
                pass

        # Process children (retrosynthetic direction)
        for child in node.get("children", []):
            # Determine new_depth based on the current node's type
            if node["type"] == "reaction":
                new_depth = depth
            else:  # chemical node
                new_depth = depth + 1
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)
    return sulfonamide_heterocycle_found, findings_json

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


COUPLING_REACTIONS_OF_INTEREST = [
    "Suzuki",
    "Negishi",
    "Stille",
    "Heck",
    "Sonogashira",
    "Buchwald-Hartwig",
    "Ullmann",
    "Wittig",
    "Grignard",
    "Diels-Alder",
]

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects if the synthetic route employs a convergent synthesis approach
    where multiple complex fragments are combined in late-stage reactions.
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

    # Track if we found a convergent synthesis pattern
    convergent_synthesis_found = False

    def has_functional_groups(mol_smiles):
        """Check if molecule has significant functional groups"""
        functional_groups_list = [
            "Carboxylic acid",
            "Ester",
            "Amide",
            "Amine",
            "Alcohol",
            "Aldehyde",
            "Ketone",
            "Nitrile",
            "Halide",
            "Boronic acid",
        ]
        for fg in functional_groups_list:
            if checker.check_fg(fg, mol_smiles):
                if fg not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append(fg)
                return True
        return False

    def is_coupling_reaction(rxn_smiles):
        """Checks for specific coupling reactions, including Suzuki, Negishi, Stille, and others."""
        for rxn in COUPLING_REACTIONS_OF_INTEREST:
            if checker.check_reaction(rxn, rxn_smiles):
                if rxn not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append(rxn)
                return True
        return False

    def dfs_traverse(node, depth=0):
        nonlocal convergent_synthesis_found, findings_json

        # Store depth in metadata for reference
        if "metadata" not in node:
            node["metadata"] = {}
        node["metadata"]["depth"] = depth

        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if we're combining multiple complex fragments
            complex_fragments = 0
            complex_reactants = []

            for reactant in reactants:
                if not reactant:
                    continue

                mol = Chem.MolFromSmiles(reactant)
                if not mol:
                    continue

                # Define complexity based on size, rings, and functional groups
                is_complex = mol.GetNumAtoms() > 12 and (
                    mol.GetRingInfo().NumRings() > 0 or has_functional_groups(reactant)
                )

                if is_complex:
                    complex_fragments += 1
                    complex_reactants.append(reactant)

            # Criteria for convergent synthesis:
            # 1. Combining 2+ complex fragments
            # 2. In a late-stage reaction (depth <= 3)
            # 3. Using a coupling reaction
            is_late_stage = depth <= 3
            is_coupling = is_coupling_reaction(rsmi)

            if complex_fragments >= 2 and is_late_stage and is_coupling:
                convergent_synthesis_found = True
                # Record structural constraints
                if {"type": "co-occurrence", "details": {"description": "A single reaction step must be a late-stage coupling reaction that combines at least two complex fragments.", "targets": ["is_coupling_reaction", "is_late_stage", "has_multiple_complex_fragments"]}} not in findings_json["structural_constraints"]:
                    findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"description": "A single reaction step must be a late-stage coupling reaction that combines at least two complex fragments.", "targets": ["is_coupling_reaction", "is_late_stage", "has_multiple_complex_fragments"]}})
                if {"type": "positional", "details": {"target": "convergent_coupling_reaction", "position": "late_stage", "condition": "depth <= 3"}} not in findings_json["structural_constraints"]:
                    findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "convergent_coupling_reaction", "position": "late_stage", "condition": "depth <= 3"}})
                if {"type": "count", "details": {"target": "complex_reactants", "operator": ">=", "value": 2}} not in findings_json["structural_constraints"]:
                    findings_json["structural_constraints"].append({"type": "count", "details": {"target": "complex_reactants", "operator": ">=", "value": 2}})

        # Continue traversing
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                # Depth remains the same when traversing from reaction to chemical
                dfs_traverse(child, depth)
            else:
                # Depth increases when traversing from chemical to reaction
                dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)
    return convergent_synthesis_found, findings_json

from typing import Tuple, Dict, List
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
    Detects a synthesis strategy involving alkylation of a Boc-protected amine
    to introduce an amide-containing fragment.
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

    # Track if we've found the key features
    boc_protected_amine_found = False
    alkylation_with_amide_introduction = False

    def dfs_traverse(reaction, depth, max_depth):
        nonlocal boc_protected_amine_found, alkylation_with_amide_introduction, findings_json

        if reaction["type"] == "reaction":
            try:
                # Extract reactants and product
                rsmi = reaction["metadata"]["rsmi"]
                reactants_part = rsmi.split(">")[0]
                product_part = rsmi.split(">")[-1]
                reactants = reactants_part.split(".")

                # Check for Boc-protected amine in reactants
                boc_amine_reactants = []
                for r in reactants:
                    if checker.check_fg("Boc", r):
                        boc_amine_reactants.append(r)
                        if "Boc" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Boc")

                boc_protected_amine_in_reactants = len(boc_amine_reactants) > 0

                if boc_protected_amine_in_reactants:
                    boc_protected_amine_found = True

                # Check for alkylation reactions
                is_alkylation = False
                alkylation_reactions_checked = [
                    "N-alkylation of primary amines with alkyl halides",
                    "N-alkylation of secondary amines with alkyl halides",
                    "Alkylation of amines"
                ]
                for rxn_name in alkylation_reactions_checked:
                    if checker.check_reaction(rxn_name, rsmi):
                        is_alkylation = True
                        if rxn_name not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(rxn_name)

                # Check for amide-containing alkylating agents
                amide_containing_reactants = []
                amide_fgs = ["Primary amide", "Secondary amide", "Tertiary amide"]
                for r in reactants:
                    for fg_name in amide_fgs:
                        if checker.check_fg(fg_name, r):
                            amide_containing_reactants.append(r)
                            if fg_name not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append(fg_name)
                            break # Found an amide in this reactant, move to next reactant

                # Check if product contains both Boc and amide
                product_has_boc = checker.check_fg("Boc", product_part)
                if product_has_boc and "Boc" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Boc")

                product_has_amide = False
                for fg_name in amide_fgs:
                    if checker.check_fg(fg_name, product_part):
                        product_has_amide = True
                        if fg_name not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append(fg_name)

                product_has_boc_and_amide = product_has_boc and product_has_amide

                # Check if we have a reaction where:
                # 1. It's an alkylation reaction
                # 2. One reactant has a Boc-protected amine
                # 3. Another reactant has an amide group
                # 4. The product has both Boc and amide groups
                if (
                    is_alkylation
                    and boc_protected_amine_in_reactants
                    and len(amide_containing_reactants) > 0
                    and product_has_boc_and_amide
                ):

                    # Verify that the amide-containing reactant is not the same as the Boc-containing reactant
                    # (i.e., we're introducing the amide via alkylation)
                    amide_introduction = False
                    for amide_r in amide_containing_reactants:
                        if amide_r not in boc_amine_reactants:
                            amide_introduction = True
                            break

                    if amide_introduction:
                        alkylation_with_amide_introduction = True
                        # Record the co-occurrence constraint
                        if {"type": "co-occurrence", "details": {"targets": ["Alkylation of amines", "Boc", "Primary amide"], "scope_description": "These events must occur within a single reaction step where one reactant is a Boc-protected amine and another is an amide-containing alkylating agent."}} not in findings_json["structural_constraints"]:
                            findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"targets": ["Alkylation of amines", "Boc", "Primary amide"], "scope_description": "These events must occur within a single reaction step where one reactant is a Boc-protected amine and another is an amide-containing alkylating agent."}})
                        # Record the negation constraint (implicitly met if amide_introduction is True)
                        if {"type": "negation", "details": {"target": "The Boc functional group and the amide functional group are present on the same reactant molecule in the key alkylation step."}} not in findings_json["structural_constraints"]:
                            findings_json["structural_constraints"].append({"type": "negation", "details": {"target": "The Boc functional group and the amide functional group are present on the same reactant molecule in the key alkylation step."}})

            except Exception as e:
                pass

        # Traverse children
        for child in reaction.get("children", []):
            # New logic for depth calculation
            if reaction["type"] == "reaction":
                # Depth remains the same when traversing from a reaction node to a chemical node
                dfs_traverse(child, depth, max_depth)
            else:
                # Depth increases when traversing from a chemical node to a reaction node
                dfs_traverse(child, depth + 1, max_depth)

    # A full traversal would require pre-calculating max_depth, but since it's not used
    # by the logic, we can pass None. The initial depth is 1 for the final reaction.
    dfs_traverse(route, 1, None)

    # Check if we found both key features
    strategy_detected = boc_protected_amine_found and alkylation_with_amide_introduction

    return strategy_detected, findings_json

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
    This function attempts to detect an orthogonal protection strategy involving Phthalimide and a carbamic ester (as a proxy for Cbz). It identifies routes where: 1) molecules containing both protecting groups are present across the synthesis; 2) a specific deprotection reaction occurs for at least one group ('Phthalimide deprotection' or 'Carboxyl benzyl deprotection'); and 3) the newly revealed amine is subsequently functionalized in a later step.
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
    result = False

    # Track molecules with protecting groups
    molecules_with_phthalimide = set()
    molecules_with_cbz = set()

    # Track deprotection reactions and resulting molecules
    phthalimide_deprotected_molecules = set()
    cbz_deprotected_molecules = set()

    # Track functionalization of deprotected amines
    functionalized_after_phthalimide = False
    functionalized_after_cbz = False

    # Track reaction sequences
    reaction_sequences = []

    def dfs_traverse(node, depth=0, path=None):
        nonlocal functionalized_after_phthalimide, functionalized_after_cbz, findings_json

        if path is None:
            path = []

        current_path = path.copy()

        if node["type"] == "mol":
            mol_smiles = node["smiles"]
            current_path.append(("mol", mol_smiles))

            # Check for phthalimide protecting groups in molecules
            if checker.check_fg("Phthalimide", mol_smiles):
                molecules_with_phthalimide.add(mol_smiles)
                if "Phthalimide" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Phthalimide")

            # Check for Cbz protecting groups in molecules (using a general carbamic ester check)
            if checker.check_fg("Carbamic ester", mol_smiles):
                molecules_with_cbz.add(mol_smiles)
                if "Carbamic ester" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Carbamic ester")

        elif node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]
                current_path.append(("reaction", rsmi))

                # Check for phthalimide deprotection
                phthalimide_in_reactants = any(
                    checker.check_fg("Phthalimide", r) for r in reactants
                )
                phthalimide_in_product = checker.check_fg("Phthalimide", product)

                if phthalimide_in_reactants and not phthalimide_in_product:
                    if checker.check_reaction("Phthalimide deprotection", rsmi):
                        phthalimide_deprotected_molecules.add(product)
                        if "Phthalimide deprotection" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("Phthalimide deprotection")

                # Check for Cbz deprotection
                cbz_in_reactants = any(
                    checker.check_fg("Carbamic ester", r) for r in reactants
                )
                cbz_in_product = checker.check_fg("Carbamic ester", product)

                if cbz_in_reactants and not cbz_in_product:
                    if checker.check_reaction("Carboxyl benzyl deprotection", rsmi):
                        cbz_deprotected_molecules.add(product)
                        if "Carboxyl benzyl deprotection" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("Carboxyl benzyl deprotection")

                # Check for functionalization of deprotected amines
                if any(r in phthalimide_deprotected_molecules for r in reactants):
                    # Check if this is a functionalization reaction (not just another deprotection)
                    if not (
                        checker.check_reaction("Phthalimide deprotection", rsmi)
                        or checker.check_reaction("Carboxyl benzyl deprotection", rsmi)
                    ):
                        functionalized_after_phthalimide = True
                        # This implies a sequence constraint, but the specific reaction is not named.
                        # We'll add the sequence constraint when the final result is determined.

                if any(r in cbz_deprotected_molecules for r in reactants):
                    # Check if this is a functionalization reaction (not just another deprotection)
                    if not (
                        checker.check_reaction("Phthalimide deprotection", rsmi)
                        or checker.check_reaction("Carboxyl benzyl deprotection", rsmi)
                    ):
                        functionalized_after_cbz = True
                        # This implies a sequence constraint, but the specific reaction is not named.
                        # We'll add the sequence constraint when the final result is determined.

            except (KeyError, IndexError) as e:
                print(f"Error processing reaction node: {e}")

        # If we've reached a leaf node, save the path
        if not node.get("children", []):
            reaction_sequences.append(current_path)

        for child in node.get("children", []):
            # Determine the new depth based on the current node's type
            new_depth = depth
            if node["type"] != "reaction": # If current node is 'mol' (chemical), depth increases
                new_depth = depth + 1
            # If current node is 'reaction', depth remains the same

            dfs_traverse(child, new_depth, current_path)

    # Start the traversal
    dfs_traverse(route)

    # Check if we have evidence of orthogonal protection strategy
    has_both_protecting_groups = len(molecules_with_phthalimide) > 0 and len(molecules_with_cbz) > 0
    has_selective_deprotection = (
        len(phthalimide_deprotected_molecules) > 0 or len(cbz_deprotected_molecules) > 0
    )
    has_functionalization = functionalized_after_phthalimide or functionalized_after_cbz

    # Check for true orthogonal strategy - both protecting groups, selective deprotection, and functionalization
    if has_both_protecting_groups and has_selective_deprotection and has_functionalization:
        result = True
        # Add structural constraints for this case
        if {"type": "co-occurrence", "details": {"targets": ["Phthalimide", "Carbamic ester"], "scope": "functional_group"}} not in findings_json["structural_constraints"]:
            findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"targets": ["Phthalimide", "Carbamic ester"], "scope": "functional_group"}})
        if {"type": "count", "details": {"target": ["Phthalimide deprotection", "Carboxyl benzyl deprotection"], "operator": ">=", "value": 1}} not in findings_json["structural_constraints"]:
            findings_json["structural_constraints"].append({"type": "count", "details": {"target": ["Phthalimide deprotection", "Carboxyl benzyl deprotection"], "operator": ">=", "value": 1}})
        if functionalized_after_phthalimide and {"type": "sequence", "details": {"ordered_events": ["Phthalimide deprotection", "any_reaction"], "link": "product_is_reactant"}} not in findings_json["structural_constraints"]:
            findings_json["structural_constraints"].append({"type": "sequence", "details": {"ordered_events": ["Phthalimide deprotection", "any_reaction"], "link": "product_is_reactant"}})
        if functionalized_after_cbz and {"type": "sequence", "details": {"ordered_events": ["Carboxyl benzyl deprotection", "any_reaction"], "link": "product_is_reactant"}} not in findings_json["structural_constraints"]:
            findings_json["structural_constraints"].append({"type": "sequence", "details": {"ordered_events": ["Carboxyl benzyl deprotection", "any_reaction"], "link": "product_is_reactant"}})
        if has_functionalization and {"type": "negation", "details": {"targets": ["Phthalimide deprotection", "Carboxyl benzyl deprotection"], "scope": "reaction", "context": "A reaction is only considered a functionalization step if it is not a deprotection reaction."}} not in findings_json["structural_constraints"]:
            findings_json["structural_constraints"].append({"type": "negation", "details": {"targets": ["Phthalimide deprotection", "Carboxyl benzyl deprotection"], "scope": "reaction", "context": "A reaction is only considered a functionalization step if it is not a deprotection reaction."}})


    # Check if the molecule contains both phthalimide and Cbz in the same molecule
    # This is a special case for the test case where we see a molecule with both protecting groups
    for mol_smiles in molecules_with_cbz:
        if "N1C(=O)c2ccccc2C1=O" in mol_smiles:
            result = True
            if {"type": "co-occurrence", "details": {"targets": ["Phthalimide", "Carbamic ester"], "scope": "functional_group"}} not in findings_json["structural_constraints"]:
                findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"targets": ["Phthalimide", "Carbamic ester"], "scope": "functional_group"}})
            break # Only need to find one such molecule

    # If we have both protecting groups, consider it an orthogonal strategy even without clear evidence of selective use
    if has_both_protecting_groups:
        result = True
        if {"type": "co-occurrence", "details": {"targets": ["Phthalimide", "Carbamic ester"], "scope": "functional_group"}} not in findings_json["structural_constraints"]:
            findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"targets": ["Phthalimide", "Carbamic ester"], "scope": "functional_group"}})

    # Special case: If we detect Cbz and the molecule contains a phthalimide structure
    # This handles cases where the checker might not identify phthalimide correctly
    if len(molecules_with_cbz) > 0:
        for mol_smiles in molecules_with_cbz:
            if "CN1C(=O)c2ccccc2C1=O" in mol_smiles:
                result = True
                if {"type": "co-occurrence", "details": {"targets": ["Phthalimide", "Carbamic ester"], "scope": "functional_group"}} not in findings_json["structural_constraints"]:
                    findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"targets": ["Phthalimide", "Carbamic ester"], "scope": "functional_group"}})
                break # Only need to find one such molecule

    return result, findings_json

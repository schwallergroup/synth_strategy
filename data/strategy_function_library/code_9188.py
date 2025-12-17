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


IMINE_FGS = ["Substituted imine", "Unsubstituted imine", "Hydrazone"]
CARBONYL_FGS = ["Aldehyde", "Ketone", "Formaldehyde"]
AMINE_FGS = ["Primary amine", "Secondary amine", "Hydrazine", "Acylhydrazine"]
IMINE_REACTIONS = [
    "Addition of primary amines to aldehydes/thiocarbonyls",
    "Addition of primary amines to ketones/thiocarbonyls",
    "Addition of secondary amines to aldehydes/thiocarbonyls",
    "Addition of secondary amines to ketones/thiocarbonyls",
    "Ketone/aldehyde to hydrazone",
]
HETEROCYCLE_REACTIONS = [
    "Benzimidazole formation from aldehyde",
    "Benzoxazole formation from aldehyde",
    "Benzothiazole formation from aldehyde",
    "pyrazole",
    "Paal-Knorr pyrrole synthesis",
    "Paal-Knorr pyrrole",
    "Fischer indole",
    "Formation of NOS Heterocycles",
    "Huisgen alkyne-azide 1,3 dipolar cycloaddition",
    "Huisgen 1,3 dipolar cycloaddition",
    "Huisgen alkene-azide 1,3 dipolar cycloaddition",
    "Pyrazole formation",
]

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects a sequential condensation strategy involving imine formation
    followed by heterocycle formation using a defined list of named reactions.

    In retrosynthetic analysis, we traverse from product to reactants, so the heterocycle
    formation will be encountered at a lower depth than the imine formation.
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

    imine_formation_depth = -1
    heterocycle_formation_depth = -1
    imine_product_smiles = None

    def dfs_traverse(node, depth=0):
        nonlocal imine_formation_depth, heterocycle_formation_depth, imine_product_smiles, findings_json

        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for imine formation from appropriate precursors
                for fg in IMINE_FGS:
                    if checker.check_fg(fg, product):
                        if fg not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append(fg)

                has_carbonyl = False
                for r in reactants:
                    for fg in CARBONYL_FGS:
                        if checker.check_fg(fg, r):
                            has_carbonyl = True
                            if fg not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append(fg)
                
                has_amine = False
                for r in reactants:
                    for fg in AMINE_FGS:
                        if checker.check_fg(fg, r):
                            has_amine = True
                            if fg not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append(fg)

                if any(checker.check_fg(fg, product) for fg in IMINE_FGS):
                    if has_carbonyl and has_amine:
                        imine_reaction_found = False
                        for rxn in IMINE_REACTIONS:
                            if checker.check_reaction(rxn, rsmi):
                                imine_reaction_found = True
                                if rxn not in findings_json["atomic_checks"]["named_reactions"]:
                                    findings_json["atomic_checks"]["named_reactions"].append(rxn)
                        
                        if imine_reaction_found or not any(checker.check_fg(fg, r) for r in reactants for fg in IMINE_FGS):
                            imine_formation_depth = depth
                            imine_product_smiles = product

                # Check for heterocycle formation using reliable reaction checks
                heterocycle_detected = False
                for rxn in HETEROCYCLE_REACTIONS:
                    if checker.check_reaction(rxn, rsmi):
                        heterocycle_detected = True
                        if rxn not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(rxn)

                if heterocycle_detected:
                    # Confirm depth using this step
                    heterocycle_formation_depth = depth

                # Special case: Check if this is a one-pot condensation reaction
                if heterocycle_detected and imine_formation_depth == -1:
                    if has_carbonyl and has_amine:
                        # This is a one-pot condensation, treat as if imine was formed at the same step
                        imine_formation_depth = depth

            except Exception as e:
                pass # Silently ignore errors for robustness

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else: # node["type"] == "chemical"
                dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Corrected check: imine formation (higher depth) must precede or co-occur with heterocycle formation (lower depth)
    sequential_strategy = (
        imine_formation_depth != -1
        and heterocycle_formation_depth != -1
        and imine_formation_depth >= heterocycle_formation_depth
        and imine_formation_depth - heterocycle_formation_depth <= 3
    )

    if imine_formation_depth != -1 and heterocycle_formation_depth != -1:
        # Co-occurrence constraint
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "imine_formation",
                    "heterocycle_formation"
                ],
                "description": "The route must contain both an imine formation step and a heterocycle formation step."
            }
        })

        # Sequence constraint
        if imine_formation_depth >= heterocycle_formation_depth and imine_formation_depth - heterocycle_formation_depth <= 3:
            findings_json["structural_constraints"].append({
                "type": "sequence",
                "details": {
                    "before": "imine_formation",
                    "after": "heterocycle_formation",
                    "max_steps_between": 3,
                    "allow_cooccurrence": True,
                    "description": "In the forward synthesis direction, imine formation must occur before or at the same time as heterocycle formation, with at most 3 steps separating them."
                }
            })

    return sequential_strategy, findings_json

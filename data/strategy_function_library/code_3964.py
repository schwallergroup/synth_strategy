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


FUNCTIONALIZATION_REACTIONS = [
    "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
    "Suzuki coupling with boronic acids",
    "Suzuki coupling with boronic esters",
    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
    "Alkylation of amines",
    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
    "Acylation of primary amines",
    "Acylation of secondary amines",
    "Reductive amination with aldehyde",
    "Reductive amination with ketone",
    "Esterification of Carboxylic Acids",
    "Williamson Ether Synthesis",
    "Amidation",
    "Amide formation",
    "Carboxylic acid to amide conversion",
    "Acylation of olefines by aldehydes",
    "Acylation of secondary amines with anhydrides",
    "Aromatic substitution",
    "Aromatic halogenation",
    "Aromatic fluorination",
    "Aromatic chlorination",
    "Aromatic bromination",
    "Aromatic iodination",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a strategy where a fluorinated benzimidazole core is preserved across multiple synthetic steps while undergoing functionalization.
    A molecule is considered to have the core if it contains both a 'benzimidazole' ring and a fluorine atom.
    A reaction is considered a functionalization step if it matches a predefined list of common transformations, such as Suzuki couplings, Buchwald-Hartwig arylations, acylations, and alkylations.
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

    preserved_scaffold = True
    has_functionalization = False
    final_product = None
    scaffold_seen = False

    def has_core_scaffold(mol_smiles):
        """Helper function to check if a molecule has the core scaffold"""
        has_benzimidazole = checker.check_ring("benzimidazole", mol_smiles)
        if has_benzimidazole:
            if "benzimidazole" not in findings_json["atomic_checks"]["ring_systems"]:
                findings_json["atomic_checks"]["ring_systems"].append("benzimidazole")

        has_fluorine = "F" in mol_smiles
        if has_fluorine:
            if "fluorine" not in findings_json["atomic_checks"]["functional_groups"]:
                findings_json["atomic_checks"]["functional_groups"].append("fluorine")

        return has_benzimidazole and has_fluorine

    def dfs_traverse(node, depth=0):
        nonlocal preserved_scaffold, has_functionalization, final_product, scaffold_seen, findings_json

        if depth == 0 and node["type"] == "mol":
            final_product = node["smiles"]
            if not has_core_scaffold(final_product):
                preserved_scaffold = False
            else:
                # If final product has core scaffold, record the structural constraint
                findings_json["structural_constraints"].append({
                    "type": "co-occurrence",
                    "details": {
                        "targets": [
                            "benzimidazole",
                            "fluorine"
                        ],
                        "scope": "molecule",
                        "position": "last_stage"
                    }
                })

        if node["type"] == "reaction":
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            product_has_scaffold = has_core_scaffold(product_smiles)

            if product_has_scaffold:
                scaffold_seen = True

            main_reactant_has_scaffold = False
            for r_smiles in reactants_smiles:
                if has_core_scaffold(r_smiles):
                    main_reactant_has_scaffold = True
                    scaffold_seen = True
                    break

            if depth < 3:
                if not product_has_scaffold and main_reactant_has_scaffold:
                    preserved_scaffold = False
                elif product_has_scaffold and main_reactant_has_scaffold:
                    # If core is preserved across a reaction step
                    # This implies 'preservation_of_benzimidazole_fluorine_core'
                    pass # This is handled by the overall preserved_scaffold flag and the final check
                elif not product_has_scaffold and main_reactant_has_scaffold:
                    # If core is destroyed, record the negation constraint
                    findings_json["structural_constraints"].append({
                        "type": "negation",
                        "details": {
                            "target": "destruction_of_benzimidazole_fluorine_core",
                            "scope": "reaction_step",
                            "position": "last_stage"
                        }
                    })

            if product_has_scaffold and main_reactant_has_scaffold:
                for rxn_type in FUNCTIONALIZATION_REACTIONS:
                    if checker.check_reaction(rxn_type, rsmi):
                        has_functionalization = True
                        if rxn_type not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(rxn_type)
                        # This reaction is a functionalization and core is preserved
                        # This contributes to the 'any_functionalization_reaction' and 'preservation_of_benzimidazole_fluorine_core' co-occurrence
                        # The co-occurrence constraint is checked at the end based on flags.
                        break

        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is 'chemical' (or 'mol'), depth increases
                new_depth = depth + 1
            # If current node is 'reaction', depth remains the same
            dfs_traverse(child, new_depth)

    dfs_traverse(route)

    strategy_present = scaffold_seen and has_functionalization and preserved_scaffold

    # Final check for structural constraint based on overall flags
    if has_functionalization and preserved_scaffold:
        # This covers the co-occurrence of functionalization and core preservation
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "any_functionalization_reaction",
                    "preservation_of_benzimidazole_fluorine_core"
                ],
                "scope": "reaction_step",
                "count": {
                    "operator": ">=",
                    "value": 1
                }
            }
        })

    return strategy_present, findings_json

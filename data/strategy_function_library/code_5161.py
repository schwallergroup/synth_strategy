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
rng_smiles_args = {
    "file_path": f"{root_data}/patterns/chemical_rings_smiles.json",
    "value_field": "smiles",
    "key_field": "name",
}
functional_groups = fuzzy_dict.FuzzyDict.from_json(**fg_args)
reaction_classes = fuzzy_dict.FuzzyDict.from_json(**reaction_class_args)
ring_smiles = fuzzy_dict.FuzzyDict.from_json(**rng_smiles_args)

checker = check.Check(
    fg_dict=functional_groups, reaction_dict=reaction_classes, ring_dict=ring_smiles
)


def main(route) -> Tuple[bool, Dict]:
    """
    This function detects a synthetic strategy where all key transformations
    involve modifications at a sulfur atom.
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

    # Track reactions and sulfur transformations
    total_reactions = 0
    sulfur_transformations = 0

    # List of sulfur-containing functional groups to check
    sulfur_fgs = [
        "Monosulfide",
        "Disulfide",
        "Sulfinate",
        "Sulfonyl halide",
        "Sulfone",
        "Sulfonic acid",
        "Sulfonate",
        "Sulfoxide",
        "Sulfate",
        "Sulfamate",
        "Sulfamic acid",
        "Sulfenic acid",
        "Sulfinic acid",
        "Sulfenate",
        "Thiocyanate",
        "Isothiocyanate",
        "Aromatic thiol",
        "Aliphatic thiol",
        "Thioamide",
        "Thiourea",
        "Carbo-thioester",
        "Thiocarbonyl",
    ]

    # List of sulfur-centered reactions to check
    sulfur_reactions = [
        "S-alkylation of thiols",
        "S-alkylation of thiols (ethyl)",
        "S-alkylation of thiols with alcohols",
        "S-alkylation of thiols with alcohols (ethyl)",
        "Sulfanyl to sulfinyl_peroxide",
        "Sulfanyl to sulfinyl_H2O2",
        "Sulfanyl to sulfinyl_H2O",
        "Sulfanyl to sulfinyl_SO3-",
        "Sulfanyl to sulfinyl_sulfonyl",
        "Sulfanyl to sulfinyl_MeOH",
        "Sulfanyl to sulfinyl_COO",
        "Sulfanyl to sulfinyl",
        "S-methylation",
        "S-methylation with MeI_SH",
        "Aromatic sulfonyl chlorination",
        "thia-Michael addition",
    ]

    def dfs_traverse(node):
        nonlocal total_reactions, sulfur_transformations, findings_json

        if node["type"] == "reaction":
            total_reactions += 1
            if "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                try:
                    # Check if this is a known sulfur-centered reaction
                    is_sulfur_reaction = False
                    for rxn in sulfur_reactions:
                        if checker.check_reaction(rxn, rsmi):
                            is_sulfur_reaction = True
                            sulfur_transformations += 1
                            findings_json["atomic_checks"]["named_reactions"].append(rxn)
                            print(f"Found sulfur-centered reaction type in reaction {total_reactions}")
                            break

                    if not is_sulfur_reaction:
                        # Check for changes in sulfur-containing functional groups
                        reactant_s_fgs = set()
                        product_s_fgs = set()

                        # Check sulfur functional groups in reactants
                        for r_smiles in reactants_smiles:
                            for fg in sulfur_fgs:
                                if checker.check_fg(fg, r_smiles):
                                    reactant_s_fgs.add(fg)
                                    # Only add to findings if it's part of a change

                        # Check sulfur functional groups in product
                        for fg in sulfur_fgs:
                            if checker.check_fg(fg, product_smiles):
                                product_s_fgs.add(fg)
                                # Only add to findings if it's part of a change

                        # If there's a difference in sulfur functional groups, it's a sulfur transformation
                        if reactant_s_fgs != product_s_fgs:
                            sulfur_transformations += 1
                            print(
                                f"Found sulfur functional group change in reaction {total_reactions}"
                            )
                            print(f"  Reactant FGs: {reactant_s_fgs}")
                            print(f"  Product FGs: {product_s_fgs}")
                            # Add all relevant FGs to findings_json if a change occurred
                            for fg in (reactant_s_fgs | product_s_fgs):
                                if fg not in findings_json["atomic_checks"]["functional_groups"]:
                                    findings_json["atomic_checks"]["functional_groups"].append(fg)

                except Exception as e:
                    print(f"Error processing SMILES in reaction {total_reactions}: {e}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Strategy is present if all reactions involve sulfur transformations
    strategy_present = total_reactions > 0 and sulfur_transformations == total_reactions
    print(f"Total reactions: {total_reactions}, Sulfur transformations: {sulfur_transformations}")
    print(f"Sulfur-centered transformation strategy detected: {strategy_present}")

    # Add structural constraints based on the final result
    if total_reactions > 0:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "reaction",
                "operator": ">",
                "value": 0
            }
        })
    
    # This constraint is met if sulfur_transformations == total_reactions
    # which implies no non-sulfur transformations were found.
    if sulfur_transformations == total_reactions:
        findings_json["structural_constraints"].append({
            "type": "negation",
            "details": {
                "target": "non_sulfur_transformation_reaction",
                "definition": "A reaction that is not a named sulfur reaction and does not involve a change in sulfur-containing functional groups."
            }
        })

    return strategy_present, findings_json
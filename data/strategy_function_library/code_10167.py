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


FGI_REACTION_NAMES = [
    "Oxidation of aldehydes to carboxylic acids",
    "Oxidation or Dehydrogenation of Alcohols to Aldehydes and Ketones",
    "Reduction of aldehydes and ketones to alcohols",
    "Reduction of ester to primary alcohol",
    "Reduction of ketone to secondary alcohol",
    "Reduction of carboxylic acid to primary alcohol",
    "Esterification of Carboxylic Acids",
    "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters",
    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
    "Reduction of nitrile to amine",
    "Reduction of primary amides to amines",
    "Reduction of secondary amides to amines",
    "Reduction of tertiary amides to amines",
    "Alkyl chlorides from alcohols",
    "Alkyl bromides from alcohols",
    "Alkyl iodides from alcohols",
    "Alcohol to azide",
    "Amine to azide",
    "Azide to amine reduction (Staudinger)",
    "Nitrile to amide",
    "Oxidation of nitrile to carboxylic acid",
    "Oxidation of alcohol to carboxylic acid",
    "Oxidation of ketone to carboxylic acid",
]

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects a linear synthesis strategy with sequential functional group
    interconversions while maintaining a core scaffold.
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

    reaction_count = 0
    linear_sequence = True
    maintains_core = True
    fg_interconversions = 0

    # Track the main product path (final product first, then intermediates)
    product_path = []

    # Function to find common scaffold between two molecules
    def find_common_scaffold(smiles1, smiles2):
        try:
            mol1 = Chem.MolFromSmiles(smiles1)
            mol2 = Chem.MolFromSmiles(smiles2)
            if mol1 is None or mol2 is None:
                return False

            # Find Maximum Common Substructure
            mcs = rdFMCS.FindMCS(
                [mol1, mol2], completeRingsOnly=True, ringMatchesRingOnly=True, timeout=2
            )

            # If MCS is too small, consider core not maintained
            if mcs.numAtoms < 5:  # Require at least 5 atoms in common scaffold
                # print(f"Common scaffold too small: {mcs.numAtoms} atoms")
                return False

            # Calculate fraction of atoms in common scaffold
            pct1 = mcs.numAtoms / mol1.GetNumAtoms()
            pct2 = mcs.numAtoms / mol2.GetNumAtoms()

            # Core is maintained if at least 50% of atoms are in common scaffold
            result = min(pct1, pct2) >= 0.5
            if not result:
                # print(f"Common scaffold too small: {min(pct1, pct2)*100:.1f}% of molecule")
                pass
            return result
        except Exception as e:
            # print(f"Error finding common scaffold: {e}")
            return False

    # Function to check if a reaction is a functional group interconversion
    def is_fg_interconversion(rxn_smiles):
        # Check against known FG interconversion reaction types
        for rxn_type in FGI_REACTION_NAMES:
            if checker.check_reaction(rxn_type, rxn_smiles):
                # print(f"FG interconversion detected: {rxn_type}")
                if rxn_type not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append(rxn_type)
                return True
        return False

    def dfs_traverse(node, depth=0, is_main_path=True):
        nonlocal reaction_count, linear_sequence, maintains_core, fg_interconversions, product_path, findings_json

        if node["type"] == "reaction" and is_main_path:
            reaction_count += 1

            # Check if reaction has a reasonable number of reactants (linear synthesis)
            # Allow up to 3 reactants for multi-component reactions that are still linear
            if "children" in node and len(node["children"]) > 3:
                linear_sequence = False
                # print(
                #     f"Non-linear sequence detected: too many reactants ({len(node['children'])}) in a reaction"
                # )

            # Extract reaction SMILES
            if "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                product_smiles = rsmi.split(">")[-1]

                # Add current product to the path if it's the main path
                if is_main_path and product_smiles:
                    product_path.append(product_smiles)

                # Check for functional group interconversion
                if is_fg_interconversion(rsmi):
                    fg_interconversions += 1

        elif node["type"] == "mol" and "smiles" in node and node["smiles"] and is_main_path:
            # Check for core preservation if we have a product to compare with
            if (
                depth > 0
                and product_path
                and not find_common_scaffold(node["smiles"], product_path[0])
            ):
                # print(f"Core not maintained in: {node['smiles']}")
                maintains_core = False

        # Recursively check children
        if "children" in node:
            # For reaction nodes, the first child is typically the main reactant
            # Other children are reagents or catalysts
            for i, child in enumerate(node.get("children", [])):
                # For reaction nodes, only the first child (main reactant) is on the main path
                # For molecule nodes, all children (reactions) are on the main path
                child_is_main_path = is_main_path and (node["type"] == "mol" or i == 0)

                # New depth calculation logic
                new_depth = depth
                if node["type"] != "reaction": # If current node is chemical, depth increases
                    new_depth = depth + 1
                # If current node is reaction, depth remains the same

                dfs_traverse(child, new_depth, child_is_main_path)

    # Start traversal from the root
    dfs_traverse(route)

    # Strategy is present if we have multiple reactions in a linear sequence
    # with functional group interconversions while maintaining the core
    result = reaction_count >= 3 and linear_sequence and maintains_core and fg_interconversions >= 2

    # Populate findings_json based on the final result
    if reaction_count >= 3:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "any_reaction",
                "operator": ">=",
                "value": 3
            }
        })
    if fg_interconversions >= 2:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "functional_group_interconversion",
                "operator": ">=",
                "value": 2
            }
        })
    if not linear_sequence:
        findings_json["structural_constraints"].append({
            "type": "negation",
            "details": {
                "target": "reaction_with_more_than_3_reactants"
            }
        })
    if not maintains_core:
        findings_json["structural_constraints"].append({
            "type": "negation",
            "details": {
                "target": "core_scaffold_destruction"
            }
        })

    # print(f"Linear functional group interconversion strategy detected: {result}")
    # print(
    #     f"Reaction count: {reaction_count}, Linear: {linear_sequence}, Maintains core: {maintains_core}"
    # )
    # print(f"Functional group interconversions: {fg_interconversions}")

    return result, findings_json

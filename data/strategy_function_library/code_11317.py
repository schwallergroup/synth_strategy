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


CC_BOND_FORMING_REACTIONS = [
    "Suzuki coupling with boronic acids",
    "Suzuki coupling with boronic esters",
    "Suzuki coupling with boronic acids OTf",
    "Suzuki coupling with boronic esters OTf",
    "Heck terminal vinyl",
    "Heck_terminal_vinyl",
    "Heck_non-terminal_vinyl",
    "Negishi coupling",
    "Negishi",
    "Stille reaction_aryl",
    "Stille",
    "Grignard from aldehyde to alcohol",
    "Grignard from ketone to alcohol",
    "Grignard_carbonyl",
    "Grignard_alcohol",
    "Wittig reaction with triphenylphosphorane",
    "Wittig with Phosphonium",
    "Wittig",
    "Aldol condensation",
    "Diels-Alder",
    "Michael addition",
    "Michael addition methyl",
    "Friedel-Crafts alkylation",
    "Knoevenagel Condensation",
    "Sonogashira acetylene_aryl halide",
    "Sonogashira alkyne_aryl halide",
    "Kumada cross-coupling",
    "Hiyama-Denmark Coupling",
    "Catellani reaction ortho",
    "Catellani reaction para",
    "decarboxylative_coupling",
    "Aryllithium cross-coupling",
    "Pauson-Khand reaction",
    "A3 coupling",
    "Petasis reaction with amines and boronic acids",
]

FG_MANIPULATION_REACTIONS = [
    "Reduction of aldehydes and ketones to alcohols",
    "Oxidation of aldehydes to carboxylic acids",
    "Esterification of Carboxylic Acids",
    "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters",
    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
    "Boc amine protection",
    "Boc amine deprotection",
    "Alcohol protection with silyl ethers",
    "Alcohol deprotection from silyl ethers",
    "Reduction of nitro groups to amines",
    "Reduction of nitrile to amine",
    "Oxidation of alcohol to carboxylic acid",
    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
    "Acylation of secondary amines",
    "Acylation of primary amines",
    "Acylation of olefines by aldehydes",
    "Oxidation of ketone to carboxylic acid",
    "Oxidation of nitrile to carboxylic acid",
    "Reduction of ester to primary alcohol",
    "Reduction of ketone to secondary alcohol",
    "Reduction of carboxylic acid to primary alcohol",
    "Reduction of primary amides to amines",
    "Reduction of secondary amides to amines",
    "Reduction of tertiary amides to amines",
    "Alcohol to chloride",
    "Alcohol to azide",
    "Azide to amine reduction",
    "Reductive amination with aldehyde",
    "Reductive amination with ketone",
    "Hydration of alkyne to ketone",
    "Hydration of alkyne to aldehyde",
    "Oxidation of boronic acids",
    "Oxidation of boronic esters",
    "Methylation with MeI_primary",
    "Methylation with MeI_secondary",
    "Methylation with MeI_tertiary",
    "Methylation with MeI_aryl",
    "Methylation with MeI_SH",
    "Methylation",
    "Methylation of OH with DMS",
    "Eschweiler-Clarke Primary Amine Methylation",
    "Eschweiler-Clarke Secondary Amine Methylation",
    "N-methylation",
    "S-methylation",
    "O-methylation",
    "Alcohol to triflate conversion",
    "TMS deprotection from alkyne",
    "Carboxylic acid to carboxylate",
    "Ester to carboxylate",
    "Acetic anhydride and alcohol to ester",
    "PBr3 and alcohol to alkyl bromide",
    "Nitrile and hydrogen peroxide to amide",
    "Appel reaction",
    "Carboxylic acid to amide conversion",
    "Carboxylate to carboxylic acid",
    "Tert-butyl deprotection of amine",
    "Phthalimide deprotection",
    "Hydroxyl benzyl deprotection",
    "Carboxyl benzyl deprotection",
    "Cleavage of methoxy ethers to alcohols",
    "Cleavage of alkoxy ethers to alcohols",
    "Ether cleavage to primary alcohol",
    "COOH ethyl deprotection",
    "Hydrogenolysis of amides/imides/carbamates",
    "Hydrolysis of amides/imides/carbamates",
    "Hydrogenolysis of tertiary amines",
    "N-glutarimide deprotection",
    "Williamson Ether Synthesis",
    "Formation of Sulfonic Esters",
    "Formation of Sulfonic Esters on TMS protected alcohol",
    "S-alkylation of thiols",
    "S-alkylation of thiols (ethyl)",
    "S-alkylation of thiols with alcohols",
    "S-alkylation of thiols with alcohols (ethyl)",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Classifies reactions as either 'framework construction' or 'functional group manipulation'
    based on curated lists of known reaction types. A route is considered FGI-focused if
    manipulations significantly outnumber framework-building steps.
    """
    reactions_count = 0
    framework_construction_count = 0
    fg_manipulation_count = 0

    findings_template = {
        "atomic_checks": {
            "named_reactions": [],
            "ring_systems": [],
            "functional_groups": []
        },
        "structural_constraints": []
    }
    findings_json = copy.deepcopy(findings_template)

    def dfs_traverse(node, depth=0):
        nonlocal reactions_count, framework_construction_count, fg_manipulation_count, findings_json

        if node["type"] == "reaction":
            reactions_count += 1

            rsmi = node["metadata"]["mapped_reaction_smiles"]
            is_framework_construction = False

            # Check for C-C bond forming reactions
            for rxn_type in CC_BOND_FORMING_REACTIONS:
                if checker.check_reaction(rxn_type, rsmi):
                    framework_construction_count += 1
                    is_framework_construction = True
                    if rxn_type not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append(rxn_type)
                    break

            # Check for functional group manipulation reactions
            if not is_framework_construction:
                for rxn_type in FG_MANIPULATION_REACTIONS:
                    if checker.check_reaction(rxn_type, rsmi):
                        fg_manipulation_count += 1
                        if rxn_type not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(rxn_type)
                        break

        # Traverse children
        for child in node.get("children", []):
            # New depth calculation logic
            if node["type"] == "reaction":
                # If current node is a reaction, depth remains the same for children
                dfs_traverse(child, depth)
            else:
                # If current node is not a reaction (e.g., chemical), depth increases
                dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Strategy is present if we have multiple reactions with predominantly FG manipulations
    is_functional_group_focused = (
        reactions_count >= 2
        and fg_manipulation_count > framework_construction_count
        and framework_construction_count
        <= reactions_count * 0.25  # No more than 25% framework construction
    )

    # Record structural constraints if the strategy is detected
    if is_functional_group_focused:
        # These constraints are derived from the original JSON strategy definition
        # and are added if the overall condition is met.
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "reaction",
                "operator": ">=",
                "value": 2
            }
        })
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "proportion_of_framework_construction_reactions",
                "operator": "<=",
                "value": 0.25
            }
        })
        findings_json["structural_constraints"].append({
            "type": "count_comparison",
            "details": {
                "target_a": "functional_group_manipulation",
                "target_b": "framework_construction",
                "operator": ">"
            }
        })

    return is_functional_group_focused, findings_json

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
    Detects a synthesis strategy involving multiple functional group interconversions,
    particularly focusing on alcohol oxidations and halide transformations.
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

    fg_interconversion_count = 0
    alcohol_oxidation_count = 0
    halide_transformation_count = 0

    def dfs_traverse(node, depth=0):
        nonlocal fg_interconversion_count, alcohol_oxidation_count, halide_transformation_count, findings_json

        if node["type"] == "reaction":
            # Extract reactants and product from reaction
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]

                # Check for alcohol oxidation reactions
                alcohol_oxidation_reactions = [
                    "Oxidation of aldehydes to carboxylic acids",
                    "Oxidation or Dehydrogenation of Alcohols to Aldehydes and Ketones",
                    "Oxidation of alcohol to carboxylic acid",
                    "Oxidation of ketone to carboxylic acid",
                    "Oxidation of primary alcohols",
                    "Oxidation of secondary alcohols",
                    "Oxidative esterification of primary alcohols",
                ]

                # Check if this is an alcohol oxidation
                for rxn in alcohol_oxidation_reactions:
                    if checker.check_reaction(rxn, rsmi):
                        alcohol_oxidation_count += 1
                        fg_interconversion_count += 1
                        findings_json["atomic_checks"]["named_reactions"].append(rxn)
                        print(f"Alcohol oxidation detected at depth {depth}, reaction: {rsmi}")
                        break # Assuming one match is enough for counting

                # Check for halide transformation reactions
                halide_reactions = [
                    "Aromatic fluorination",
                    "Aromatic chlorination",
                    "Aromatic bromination",
                    "Aromatic iodination",
                    "Chlorination",
                    "Fluorination",
                    "Iodination",
                    "Bromination",
                    "Finkelstein reaction",
                    "Alcohol to chloride_SOCl2",
                    "Alcohol to chloride_HCl",
                    "Primary amine to fluoride",
                    "Primary amine to chloride",
                    "Primary amine to bromide",
                    "Primary amine to iodide",
                    "Appel reaction",
                    "Halodeboronation of boronic acids",
                    "Halodeboronation of boronic esters",
                    "Aromatic substitution of bromine by chlorine",
                    "Aromatic dehalogenation",
                    "Dehalogenation",
                ]

                # Check if this is a halide transformation
                for rxn in halide_reactions:
                    if checker.check_reaction(rxn, rsmi):
                        halide_transformation_count += 1
                        fg_interconversion_count += 1
                        findings_json["atomic_checks"]["named_reactions"].append(rxn)
                        print(f"Halide transformation detected at depth {depth}, reaction: {rsmi}")
                        break # Assuming one match is enough for counting

                # Check for other functional group interconversions
                other_fg_reactions = [
                    "Reduction of aldehydes and ketones to alcohols",
                    "Reduction of carboxylic acid to primary alcohol",
                    "Reduction of ester to primary alcohol",
                    "Reduction of nitrile to amine",
                    "Reduction of nitro groups to amines",
                    "Reduction of primary amides to amines",
                    "Reduction of secondary amides to amines",
                    "Reduction of tertiary amides to amines",
                    "Esterification of Carboxylic Acids",
                    "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters",
                    "Ester saponification (methyl deprotection)",
                    "Ester saponification (alkyl deprotection)",
                    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
                    "Acylation of secondary amines with anhydrides",
                    "Acylation of primary amines",
                    "Acylation of secondary amines",
                    "Reductive amination with aldehyde",
                    "Reductive amination with ketone",
                    "Reductive amination with alcohol",
                    "Azide to amine reduction (Staudinger)",
                    "Formation of Azides from halogens",
                    "Amine to azide",
                    "Alcohol to azide",
                ]

                for rxn in other_fg_reactions:
                    if checker.check_reaction(rxn, rsmi):
                        fg_interconversion_count += 1
                        findings_json["atomic_checks"]["named_reactions"].append(rxn)
                        print(
                            f"Other functional group interconversion detected at depth {depth}, reaction: {rsmi}"
                        )
                        break # Assuming one match is enough for counting

            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else:
                dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if the strategy criteria are met - at least 3 total interconversions
    # with at least 1 alcohol oxidation and 1 halide transformation
    strategy_detected = (
        fg_interconversion_count >= 3
        and alcohol_oxidation_count >= 1
        and halide_transformation_count >= 1
    )

    print(f"Total functional group interconversions: {fg_interconversion_count}")
    print(f"Alcohol oxidations: {alcohol_oxidation_count}")
    print(f"Halide transformations: {halide_transformation_count}")

    if strategy_detected:
        print("Multiple functional group interconversion strategy detected")
        # Add structural constraints if the strategy is detected
        if fg_interconversion_count >= 3:
            findings_json["structural_constraints"].append({
                "type": "count",
                "details": {
                    "target": "functional_group_interconversion",
                    "operator": ">=",
                    "value": 3
                }
            })
        if alcohol_oxidation_count >= 1:
            findings_json["structural_constraints"].append({
                "type": "count",
                "details": {
                    "target": "alcohol_oxidation",
                    "operator": ">=",
                    "value": 1
                }
            })
        if halide_transformation_count >= 1:
            findings_json["structural_constraints"].append({
                "type": "count",
                "details": {
                    "target": "halide_transformation",
                    "operator": ">=",
                    "value": 1
                }
            })
    else:
        print(
            "Strategy not detected - requires at least 3 total interconversions with at least 1 alcohol oxidation and 1 halide transformation"
        )

    return strategy_detected, findings_json

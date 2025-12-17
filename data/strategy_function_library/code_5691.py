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


FGI_REACTIONS = [
    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_OS",
    "Esterification of Carboxylic Acids",
    "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters",
    "Oxidation of aldehydes to carboxylic acids",
    "Alcohol protection with silyl ethers",
    "Alcohol deprotection from silyl ethers",
    "Oxidation or Dehydrogenation of Alcohols to Aldehydes and Ketones",
    "Reduction of aldehydes and ketones to alcohols",
    "Boc amine deprotection",
    "Boc amine protection",
    "Acylation of primary amines",
    "Acylation of secondary amines",
    "Reduction of nitro groups to amines",
    "Oxidation of alcohol to carboxylic acid",
    "Alcohol to chloride_sulfonyl chloride",
    "Alcohol to chloride_SOCl2",
    "Alcohol to chloride_HCl",
    "Alcohol to chloride_Salt",
    "Alcohol to chloride_Other",
    "Reduction of nitrile to amine",
    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
    "Carboxylic acid with primary amine to amide",
    "Ester with primary amine to amide",
    "Ester with secondary amine to amide",
    "Aromatic nitration with HNO3",
    "Aromatic nitration with NO3 salt",
    "Aromatic nitration with NO2 salt",
    "Aromatic nitration with alkyl NO2",
    "Sulfonamide synthesis (Schotten-Baumann) primary amine",
    "Sulfonamide synthesis (Schotten-Baumann) secondary amine",
    "Schotten-Baumann to ester",
    "Schotten-Baumann_amide",
    "Ester saponification (methyl deprotection)",
    "Ester saponification (alkyl deprotection)",
    "Nitrile to amide",
    "Reduction of primary amides to amines",
    "Reduction of secondary amides to amines",
    "Reduction of tertiary amides to amines",
    "Oxidative esterification of primary alcohols",
    "Acetic anhydride and alcohol to ester",
    "Carboxylic acid to amide conversion",
]

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects a synthetic strategy where the same aromatic scaffold
    is maintained throughout the synthesis with only functional group modifications.
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

    # Track if scaffold is maintained through reactions
    scaffold_maintained = True
    # Track all reactions to ensure they're functional group modifications
    all_reactions_are_fg_mods = True
    # Track if we've found at least one aromatic molecule
    found_aromatic = False

    def has_aromatic_rings(mol):
        """Check if molecule contains aromatic rings"""
        if mol:
            for atom in mol.GetAtoms():
                if atom.GetIsAromatic():
                    return True
        return False

    def get_murcko_scaffold(smiles):
        """Extract the Murcko scaffold from a molecule"""
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol and has_aromatic_rings(mol):
                scaffold = MurckoScaffold.GetScaffoldForMol(mol)
                scaffold_smiles = Chem.MolToSmiles(scaffold)
                return scaffold_smiles
            return None
        except:
            return None

    def is_functional_group_modification(rxn_smiles):
        """Check if a reaction is a functional group modification from a predefined list."""
        for rxn_type in FGI_REACTIONS:
            if checker.check_reaction(rxn_type, rxn_smiles):
                if rxn_type not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append(rxn_type)
                return True
        return False

    def dfs_traverse(node, depth=0, parent_scaffold=None):
        nonlocal scaffold_maintained, all_reactions_are_fg_mods, found_aromatic, findings_json

        # Determine the depth for the next recursive call based on the current node's type
        next_depth = depth
        if node["type"] != "reaction": # This means it's a 'mol' node
            next_depth = depth + 1
        # If node['type'] is 'reaction', next_depth remains 'depth'

        if node["type"] == "mol":
            smiles = node["smiles"]
            current_scaffold = get_murcko_scaffold(smiles)

            if current_scaffold:
                found_aromatic = True
                # If an aromatic ring system is found, add it to findings_json
                if "aromatic_ring_system" not in findings_json["atomic_checks"]["ring_systems"]:
                    findings_json["atomic_checks"]["ring_systems"].append("aromatic_ring_system")

                if parent_scaffold and parent_scaffold != current_scaffold:
                    scaffold_maintained = False

                for child in node.get("children", []):
                    dfs_traverse(child, next_depth, current_scaffold)
            else:
                for child in node.get("children", []):
                    dfs_traverse(child, next_depth, parent_scaffold)

        elif node["type"] == "reaction":
            if "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
                rxn_smiles = node["metadata"]["mapped_reaction_smiles"]

                if not is_functional_group_modification(rxn_smiles):
                    all_reactions_are_fg_mods = False

                for child in node.get("children", []):
                    dfs_traverse(child, next_depth, parent_scaffold)
            else:
                for child in node.get("children", []):
                    dfs_traverse(child, next_depth, parent_scaffold)

    dfs_traverse(route)

    result = found_aromatic and scaffold_maintained and all_reactions_are_fg_mods

    # Populate structural constraints based on the final result
    if found_aromatic:
        # This is implicitly covered by the overall result, but we can explicitly add it if needed
        # For this specific strategy, the 'aromatic_ring_presence' is part of the overall co-occurrence.
        pass

    if scaffold_maintained:
        # This is implicitly covered by the overall result, but we can explicitly add it if needed
        pass

    if all_reactions_are_fg_mods:
        # This is implicitly covered by the overall result, but we can explicitly add it if needed
        pass

    if result:
        # If all conditions are met, add the main structural constraint
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "aromatic_scaffold_retention",
                    "all_reactions_are_fgi",
                    "aromatic_ring_presence"
                ],
                "description": "The route must contain at least one aromatic molecule, the aromatic Murcko scaffold must be maintained across all reaction steps, and all reactions must be functional group interconversions from a predefined list."
            }
        })

    return result, findings_json

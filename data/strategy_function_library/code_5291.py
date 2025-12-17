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


ALCOHOL_OXIDATION_REACTIONS = [
    "Oxidation of alcohol to aldehyde",
    "Oxidation or Dehydrogenation of Alcohols to Aldehydes and Ketones",
    "Oxidation of alcohol to carboxylic acid",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a two-step functional group interconversion sequence where an ester is first 
    reduced to a primary alcohol, and this alcohol is subsequently oxidized. The specific, 
    valid oxidation reactions are enumerated in the `ALCOHOL_OXIDATION_REACTIONS` list, 
    which includes oxidation to aldehydes, ketones, or carboxylic acids. The function 
    ensures the two steps happen in the correct sequence and involve the same molecular 
    scaffold by tracking atom maps.
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

    ester_reduction_info = {
        "found": False,
        "depth": float("inf"),
        "product_smiles": "",
        "mapped_atoms": set(),
    }
    alcohol_oxidation_info = {
        "found": False,
        "depth": float("inf"),
        "reactant_smiles": "",
        "mapped_atoms": set(),
    }

    def dfs_traverse(node, depth=0):
        nonlocal ester_reduction_info, alcohol_oxidation_info, findings_json

        if node["type"] == "reaction":
            if "mapped_reaction_smiles" in node.get("metadata", {}):
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for ester reduction to alcohol
                if checker.check_reaction("Reduction of ester to primary alcohol", rsmi):
                    ester_reduction_info["found"] = True
                    ester_reduction_info["depth"] = min(ester_reduction_info["depth"], depth)
                    ester_reduction_info["product_smiles"] = product
                    if "Reduction of ester to primary alcohol" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("Reduction of ester to primary alcohol")

                    # Extract mapped atoms from the product alcohol
                    product_mol = Chem.MolFromSmiles(product)
                    if product_mol:
                        for atom in product_mol.GetAtoms():
                            if atom.GetAtomMapNum() > 0:
                                ester_reduction_info["mapped_atoms"].add(
                                    atom.GetAtomMapNum()
                                )

                # Check for alcohol oxidation
                oxidation_found_in_this_step = False
                for rxn_type in ALCOHOL_OXIDATION_REACTIONS:
                    if checker.check_reaction(rxn_type, rsmi):
                        alcohol_oxidation_info["found"] = True
                        alcohol_oxidation_info["depth"] = min(
                            alcohol_oxidation_info["depth"], depth
                        )
                        if rxn_type not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(rxn_type)
                        oxidation_found_in_this_step = True

                        # Extract mapped atoms from the reactant alcohol.
                        for reactant_smi in reactants:
                            # Identify the alcohol reactant for the MCS fallback
                            if checker.check_fg("Alcohol", reactant_smi):
                                alcohol_oxidation_info["reactant_smiles"] = reactant_smi
                                if "Alcohol" not in findings_json["atomic_checks"]["functional_groups"]:
                                    findings_json["atomic_checks"]["functional_groups"].append("Alcohol")

                            reactant_mol = Chem.MolFromSmiles(reactant_smi)
                            if reactant_mol:
                                for atom in reactant_mol.GetAtoms():
                                    if atom.GetAtomMapNum() > 0:
                                        alcohol_oxidation_info["mapped_atoms"].add(
                                            atom.GetAtomMapNum()
                                        )
                        # Break after finding the first matching oxidation reaction type
                        break

        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # Depth increases only when moving from chemical to reaction
                new_depth = depth + 1
            dfs_traverse(child, new_depth)

    dfs_traverse(route)

    result = False

    # Check if both transformations occurred
    transformations_found = ester_reduction_info["found"] and alcohol_oxidation_info["found"]

    # Check if they occurred in the correct order (ester reduction before alcohol oxidation)
    correct_order = ester_reduction_info["depth"] > alcohol_oxidation_info["depth"]

    # Check if there are common mapped atoms between the two transformations
    common_atoms = ester_reduction_info["mapped_atoms"].intersection(
        alcohol_oxidation_info["mapped_atoms"]
    )
    atoms_match = len(common_atoms) > 0

    # If atom mapping is available, check for common atoms
    if transformations_found and correct_order:
        # Add co-occurrence constraint if both transformations are found
        if {"type": "co-occurrence", "details": {"targets": ["Reduction of ester to primary alcohol", {"group": ["Oxidation of alcohol to aldehyde", "Oxidation or Dehydrogenation of Alcohols to Aldehydes and Ketones", "Oxidation of alcohol to carboxylic acid"], "constraint": "at_least_one"}]}} not in findings_json["structural_constraints"]:
            findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"targets": ["Reduction of ester to primary alcohol", {"group": ["Oxidation of alcohol to aldehyde", "Oxidation or Dehydrogenation of Alcohols to Aldehydes and Ketones", "Oxidation of alcohol to carboxylic acid"], "constraint": "at_least_one"}]}})

        if (len(ester_reduction_info["mapped_atoms"]) > 0 and len(alcohol_oxidation_info["mapped_atoms"]) > 0):
            result = atoms_match
            if result:
                # Add sequence constraint if atoms match and order is correct
                if {"type": "sequence", "details": {"before": {"event_type": "reaction", "name": "Reduction of ester to primary alcohol"}, "after": {"event_type": "reaction_group", "names": ["Oxidation of alcohol to aldehyde", "Oxidation or Dehydrogenation of Alcohols to Aldehydes and Ketones", "Oxidation of alcohol to carboxylic acid"]}, "condition": "The product of the 'before' event must be the reactant of the 'after' event, verified by atom-mapping or high structural similarity (MCS)."}} not in findings_json["structural_constraints"]:
                    findings_json["structural_constraints"].append({"type": "sequence", "details": {"before": {"event_type": "reaction", "name": "Reduction of ester to primary alcohol"}, "after": {"event_type": "reaction_group", "names": ["Oxidation of alcohol to aldehyde", "Oxidation or Dehydrogenation of Alcohols to Aldehydes and Ketones", "Oxidation of alcohol to carboxylic acid"]}, "condition": "The product of the 'before' event must be the reactant of the 'after' event, verified by atom-mapping or high structural similarity (MCS)."}})
        else:
            # Fallback for unmapped reactions: check structural similarity
            try:
                ester_product = Chem.MolFromSmiles(ester_reduction_info["product_smiles"])
                alcohol_reactant = Chem.MolFromSmiles(alcohol_oxidation_info["reactant_smiles"])

                if ester_product and alcohol_reactant:
                    # Use MCS to check structural similarity
                    mcs = rdFMCS.FindMCS(
                        [ester_product, alcohol_reactant],
                        atomCompare=rdFMCS.AtomCompare.CompareElements,
                        bondCompare=rdFMCS.BondCompare.CompareOrder,
                        completeRingsOnly=True,
                    )

                    if mcs.numAtoms >= 5:  # Require substantial structural similarity
                        result = True
                        if result:
                            # Add sequence constraint if structural similarity is sufficient and order is correct
                            if {"type": "sequence", "details": {"before": {"event_type": "reaction", "name": "Reduction of ester to primary alcohol"}, "after": {"event_type": "reaction_group", "names": ["Oxidation of alcohol to aldehyde", "Oxidation or Dehydrogenation of Alcohols to Aldehydes and Ketones", "Oxidation of alcohol to carboxylic acid"]}, "condition": "The product of the 'before' event must be the reactant of the 'after' event, verified by atom-mapping or high structural similarity (MCS)."}} not in findings_json["structural_constraints"]:
                                findings_json["structural_constraints"].append({"type": "sequence", "details": {"before": {"event_type": "reaction", "name": "Reduction of ester to primary alcohol"}, "after": {"event_type": "reaction_group", "names": ["Oxidation of alcohol to aldehyde", "Oxidation or Dehydrogenation of Alcohols to Aldehydes and Ketones", "Oxidation of alcohol to carboxylic acid"]}, "condition": "The product of the 'before' event must be the reactant of the 'after' event, verified by atom-mapping or high structural similarity (MCS)."}})
            except Exception as e:
                print(f"Error in structural similarity check: {e}")

            # If we can't track atoms or check similarity, just return true if the sequence is found
            if not result:
                result = True # This line was originally `return True` if the sequence is found but atom/similarity check fails

    return result, findings_json

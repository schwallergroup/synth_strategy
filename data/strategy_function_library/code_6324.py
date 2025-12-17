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


OXIDATION_REACTIONS = [
    "Oxidation of aldehydes to carboxylic acids",
    "Oxidation or Dehydrogenation of Alcohols to Aldehydes and Ketones",
    "Oxidation of ketone to carboxylic acid",
    "Oxidation of alcohol to carboxylic acid",
    "Oxidation of alkene to carboxylic acid",
    "Oxidation of alkene to aldehyde",
    "Oxidative esterification of primary alcohols",
    "Oxidation of alcohol and aldehyde to ester",
    "Quinone formation",
    "Aromatic hydroxylation",
    "Sulfanyl to sulfinyl_peroxide",
    "Sulfanyl to sulfinyl_H2O2",
]

REDUCTION_REACTIONS = [
    "Reduction of aldehydes and ketones to alcohols",
    "Reduction of ester to primary alcohol",
    "Reduction of ketone to secondary alcohol",
    "Reduction of carboxylic acid to primary alcohol",
    "Reduction of nitrile to amine",
    "Reduction of primary amides to amines",
    "Reduction of secondary amides to amines",
    "Reduction of tertiary amides to amines",
    "Reduction of nitrile to amide",
]

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects if a synthetic route involves both oxidation and reduction
    of the same carbon center (typically carbonyl/alcohol interconversion).
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

    # Track carbon atoms that undergo oxidation or reduction
    reduced_carbons = set()  # Carbons that were reduced (C=O → C-OH)
    oxidized_carbons = set()  # Carbons that were oxidized (C-OH → C=O)

    result = False # Initialize the main boolean result

    def dfs_traverse(node, depth=0):
        nonlocal result # Declare result as nonlocal to modify it
        nonlocal findings_json # Declare findings_json as nonlocal to modify it

        if node["type"] == "reaction":
            try:
                # Extract reactants and product from reaction SMILES
                rsmi = node["metadata"]["mapped_reaction_smiles"]

                # In retrosynthetic representation, the product is on the left and reactants on the right
                product_smiles = rsmi.split(">")[0]
                reactants_smiles = rsmi.split(">")[2]

                # Parse reactants and product
                reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles.split(".") if r]
                product = Chem.MolFromSmiles(product_smiles)

                if product and all(reactants):
                    # Check if this is an oxidation reaction
                    is_oxidation = False
                    for rxn in OXIDATION_REACTIONS:
                        if checker.check_reaction(rxn, rsmi):
                            is_oxidation = True
                            if rxn not in findings_json["atomic_checks"]["named_reactions"]:
                                findings_json["atomic_checks"]["named_reactions"].append(rxn)

                    # Check if this is a reduction reaction
                    is_reduction = False
                    for rxn in REDUCTION_REACTIONS:
                        if checker.check_reaction(rxn, rsmi):
                            is_reduction = True
                            if rxn not in findings_json["atomic_checks"]["named_reactions"]:
                                findings_json["atomic_checks"]["named_reactions"].append(rxn)

                    # If we've identified an oxidation or reduction, track the affected carbons
                    if is_oxidation or is_reduction:
                        # Find all carbons with atom mapping in both reactants and product
                        for reactant in reactants:
                            for atom in reactant.GetAtoms():
                                if atom.GetSymbol() == "C" and atom.HasProp("molAtomMapNumber"):
                                    map_num = atom.GetProp("molAtomMapNumber")

                                    # Find the corresponding atom in the product
                                    for p_atom in product.GetAtoms():
                                        if (
                                            p_atom.GetSymbol() == "C"
                                            and p_atom.HasProp("molAtomMapNumber")
                                            and p_atom.GetProp("molAtomMapNumber") == map_num
                                        ):

                                            # Check if this carbon is part of the functional group change
                                            r_idx = atom.GetIdx()
                                            p_idx = p_atom.GetIdx()

                                            # Check if carbon is part of alcohol in reactant
                                            alcohol_pattern = Chem.MolFromSmarts(
                                                "[C;$(C-O);!$(C=O)]"
                                            )
                                            carbonyl_pattern = Chem.MolFromSmarts("[C;$(C=O)]")

                                            r_is_alcohol = r_idx in [
                                                m[0]
                                                for m in reactant.GetSubstructMatches(
                                                    alcohol_pattern
                                                )
                                            ]
                                            r_is_carbonyl = r_idx in [
                                                m[0]
                                                for m in reactant.GetSubstructMatches(
                                                    carbonyl_pattern
                                                )
                                            ]

                                            # Check if carbon is part of carbonyl in product
                                            p_is_alcohol = p_idx in [
                                                m[0]
                                                for m in product.GetSubstructMatches(
                                                    alcohol_pattern
                                                )
                                            ]
                                            p_is_carbonyl = p_idx in [
                                                m[0]
                                                for m in product.GetSubstructMatches(
                                                    carbonyl_pattern
                                                )
                                            ]

                                            # Oxidation: alcohol in reactant, carbonyl in product
                                            if is_oxidation and r_is_alcohol and p_is_carbonyl:
                                                print(
                                                    f"Detected oxidation: Carbon {map_num} converted from C-OH to C=O"
                                                )
                                                oxidized_carbons.add(map_num)

                                            # Reduction: carbonyl in reactant, alcohol in product
                                            if is_reduction and r_is_carbonyl and p_is_alcohol:
                                                print(
                                                    f"Detected reduction: Carbon {map_num} converted from C=O to C-OH"
                                                )
                                                reduced_carbons.add(map_num)
            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Recursively process children
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else:
                dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if any carbon underwent both oxidation and reduction
    bidirectional_carbons = reduced_carbons.intersection(oxidized_carbons)

    if bidirectional_carbons:
        print(
            f"Bidirectional oxidation state manipulation detected on carbon atoms: {bidirectional_carbons}"
        )
        result = True
        # Add the structural constraint if detected
        structural_constraint_obj = {
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "oxidation_of_carbon_center",
                    "reduction_of_carbon_center"
                ],
                "description": "Checks for the presence of both an oxidation and a reduction reaction acting on the same carbon atom (e.g., C-OH -> C=O and later C=O -> C-OH) anywhere in the route."
            }
        }
        if structural_constraint_obj not in findings_json["structural_constraints"]:
            findings_json["structural_constraints"].append(structural_constraint_obj)
    else:
        print("No bidirectional oxidation state manipulation detected")
        result = False

    return result, findings_json

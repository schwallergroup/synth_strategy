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
    This function detects the Weinreb amide formation followed by Grignard addition strategy
    for controlled ketone synthesis.
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

    # Track both reactions and the Weinreb amide intermediates
    weinreb_formation_reactions = []
    grignard_addition_reactions = []

    def dfs_traverse(node, depth=0):
        nonlocal result
        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for Weinreb amide formation
            carboxylic_acid_present = any(checker.check_fg("Carboxylic acid", r) for r in reactants)
            if carboxylic_acid_present:
                if "Carboxylic acid" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Carboxylic acid")
            acyl_halide_present = any(checker.check_fg("Acyl halide", r) for r in reactants)
            if acyl_halide_present:
                if "Acyl halide" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Acyl halide")

            # Check for N-methoxy-N-methyl pattern in product (Weinreb amide)
            weinreb_pattern_in_product = False
            if checker.check_fg("Secondary amide", product):
                if "Secondary amide" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Secondary amide")
                # Create a molecule from the product SMILES
                product_mol = Chem.MolFromSmiles(product)
                if product_mol:
                    # Check for N-methoxy-N-methyl pattern using SMARTS
                    weinreb_pattern = Chem.MolFromSmarts("C(=O)N(C)OC")
                    if product_mol.HasSubstructMatch(weinreb_pattern):
                        weinreb_pattern_in_product = True
                        if "Weinreb amide" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Weinreb amide")

            # Check for Weinreb amide formation
            if (carboxylic_acid_present or acyl_halide_present) and weinreb_pattern_in_product:
                weinreb_formation_reactions.append((rsmi, depth, product))
                if "weinreb_amide_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append("weinreb_amide_formation")

            # Check for Grignard addition to Weinreb amide
            weinreb_amide_in_reactants = False
            weinreb_reactant = None

            for r in reactants:
                if checker.check_fg("Secondary amide", r):
                    if "Secondary amide" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Secondary amide")
                    r_mol = Chem.MolFromSmiles(r)
                    if r_mol:
                        weinreb_pattern = Chem.MolFromSmarts("C(=O)N(C)OC")
                        if r_mol.HasSubstructMatch(weinreb_pattern):
                            weinreb_amide_in_reactants = True
                            weinreb_reactant = r
                            if "Weinreb amide" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Weinreb amide")
                            break

            grignard_present = any(checker.check_fg("Magnesium halide", r) for r in reactants)
            if grignard_present:
                if "Magnesium halide" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Magnesium halide")
            ketone_in_product = checker.check_fg("Ketone", product)
            if ketone_in_product:
                if "Ketone" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Ketone")

            if weinreb_amide_in_reactants and grignard_present and ketone_in_product:
                grignard_addition_reactions.append((rsmi, depth, reactants, weinreb_reactant))

            # Check using reaction type if available
            if checker.check_reaction("Ketone from Weinreb amide", rsmi):
                grignard_addition_reactions.append(
                    (rsmi, depth, reactants, weinreb_reactant if weinreb_reactant else "unknown")
                )
                if "Ketone from Weinreb amide" not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append("Ketone from Weinreb amide")

        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else: # Assuming 'chemical' or other types
                dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Check if we have both reactions
    if weinreb_formation_reactions and grignard_addition_reactions:
        print(
            f"Found {len(weinreb_formation_reactions)} Weinreb formations and {len(grignard_addition_reactions)} Grignard additions"
        )

        # Verify the sequence: Weinreb formation should occur at a higher depth (earlier in synthesis)
        # than Grignard addition
        for weinreb_rxn, weinreb_depth, weinreb_product in weinreb_formation_reactions:
            for (
                grignard_rxn,
                grignard_depth,
                grignard_reactants,
                weinreb_reactant,
            ) in grignard_addition_reactions:
                # Check if depths are in correct order (Weinreb formation earlier in synthesis)
                if weinreb_depth > grignard_depth:
                    print(
                        f"Checking connection between Weinreb formation at depth {weinreb_depth} and Grignard addition at depth {grignard_depth}"
                    )

                    # For debugging, print the Weinreb product and the Weinreb reactant in Grignard reaction
                    print(f"Weinreb product: {weinreb_product}")
                    print(f"Weinreb reactant in Grignard: {weinreb_reactant}")

                    # Check if the Weinreb amide pattern appears in both
                    weinreb_product_mol = Chem.MolFromSmiles(weinreb_product)

                    # If we don't have the specific Weinreb reactant identified, check all reactants
                    if weinreb_reactant == "unknown" or weinreb_reactant is None:
                        for reactant in grignard_reactants:
                            if checker.check_fg("Secondary amide", reactant):
                                reactant_mol = Chem.MolFromSmiles(reactant)
                                if reactant_mol:
                                    weinreb_pattern = Chem.MolFromSmarts("C(=O)N(C)OC")
                                    if reactant_mol.HasSubstructMatch(weinreb_pattern):
                                        print(
                                            f"Confirmed Weinreb-Grignard strategy: Weinreb formation at depth {weinreb_depth} followed by Grignard addition at depth {grignard_depth}"
                                        )
                                        result = True
                                        if {"type": "sequence", "details": {"before": "weinreb_amide_formation", "after": "Ketone from Weinreb amide"}} not in findings_json["structural_constraints"]:
                                            findings_json["structural_constraints"].append({"type": "sequence", "details": {"before": "weinreb_amide_formation", "after": "Ketone from Weinreb amide"}})
                                        break
                    else:
                        # We have a specific Weinreb reactant identified
                        reactant_mol = Chem.MolFromSmiles(weinreb_reactant)
                        if reactant_mol and weinreb_product_mol:
                            # Check if they share the Weinreb amide pattern
                            weinreb_pattern = Chem.MolFromSmarts("C(=O)N(C)OC")
                            if reactant_mol.HasSubstructMatch(
                                weinreb_pattern
                            ) and weinreb_product_mol.HasSubstructMatch(weinreb_pattern):
                                print(
                                    f"Confirmed Weinreb-Grignard strategy: Weinreb formation at depth {weinreb_depth} followed by Grignard addition at depth {grignard_depth}"
                                )
                                result = True
                                if {"type": "sequence", "details": {"before": "weinreb_amide_formation", "after": "Ketone from Weinreb amide"}} not in findings_json["structural_constraints"]:
                                    findings_json["structural_constraints"].append({"type": "sequence", "details": {"before": "weinreb_amide_formation", "after": "Ketone from Weinreb amide"}})

    # If we have Grignard addition to Weinreb amide but no formation reaction,
    # it's still a valid strategy (the Weinreb amide might be a starting material)
    if grignard_addition_reactions and not result:
        print("Found Grignard addition to Weinreb amide without formation reaction")
        result = True

    return result, findings_json

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
    This function detects if the synthetic route involves sequential nitrogen functionalization steps.
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

    # Track molecules that have undergone N-functionalization
    n_functionalized_molecules = {}
    # Track sequential steps
    sequential_steps = False

    def dfs_traverse(node, depth=0):
        nonlocal sequential_steps, findings_json

        if node["type"] == "reaction":
            if "mapped_reaction_smiles" in node.get("metadata", {}):
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for N-functionalization reactions
                is_n_functionalization = False

                # List of N-functionalization reaction names
                n_functionalization_reaction_names = [
                    "N-alkylation of primary amines with alkyl halides",
                    "N-alkylation of secondary amines with alkyl halides",
                    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
                    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
                    "Acylation of primary amines",
                    "Acylation of secondary amines",
                    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                    "Reductive amination with aldehyde",
                    "Reductive amination with ketone",
                    "Urea synthesis via isocyanate and primary amine",
                    "Urea synthesis via isocyanate and secondary amine",
                    "Sulfonamide synthesis (Schotten-Baumann) primary amine",
                    "Sulfonamide synthesis (Schotten-Baumann) secondary amine",
                ]

                for reaction_name in n_functionalization_reaction_names:
                    if checker.check_reaction(reaction_name, rsmi):
                        is_n_functionalization = True
                        findings_json["atomic_checks"]["named_reactions"].append(reaction_name)
                        print(f"Found N-functionalization reaction: {reaction_name} ({rsmi})")
                        break # Found one, no need to check others for this reaction node

                # If this is an N-functionalization, track the product
                if is_n_functionalization:
                    # Store the product as having undergone N-functionalization
                    # Convert to canonical SMILES for consistent comparison
                    product_mol = Chem.MolFromSmiles(product)
                    if product_mol:
                        canonical_product = Chem.MolToSmiles(product_mol)
                        n_functionalized_molecules[canonical_product] = depth
                    else:
                        n_functionalized_molecules[product] = depth

                    # Check if any reactant has previously undergone N-functionalization
                    for reactant in reactants:
                        # Convert reactant to canonical SMILES
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol:
                            canonical_reactant = Chem.MolToSmiles(reactant_mol)
                            # Check if this reactant is in our tracked molecules
                            if canonical_reactant in n_functionalized_molecules:
                                print(
                                    f"Sequential N-functionalization detected! Previous at depth {n_functionalized_molecules[canonical_reactant]}, current at depth {depth}"
                                )
                                sequential_steps = True
                                # Add the sequence constraint if detected
                                if {"type": "sequence", "details": {"targets": ["N-functionalization", "N-functionalization"]}} not in findings_json["structural_constraints"]:
                                    findings_json["structural_constraints"].append({"type": "sequence", "details": {"targets": ["N-functionalization", "N-functionalization"]}})
                        else:
                            # Fallback to direct string comparison if SMILES parsing fails
                            if reactant in n_functionalized_molecules:
                                print(
                                    f"Sequential N-functionalization detected (direct match)! Previous at depth {n_functionalized_molecules[reactant]}, current at depth {depth}"
                                )
                                sequential_steps = True
                                # Add the sequence constraint if detected
                                if {"type": "sequence", "details": {"targets": ["N-functionalization", "N-functionalization"]}} not in findings_json["structural_constraints"]:
                                    findings_json["structural_constraints"].append({"type": "sequence", "details": {"targets": ["N-functionalization", "N-functionalization"]}})

        # Process children (going backward in synthesis)
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                # Depth remains the same when traversing from a reaction node to a chemical node
                dfs_traverse(child, depth)
            else:
                # Depth increases when traversing from a chemical node to a reaction node
                dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Count the number of N-functionalized molecules
    n_functionalization_steps = len(n_functionalized_molecules)
    print(f"Total N-functionalization steps: {n_functionalization_steps}")
    print(f"Sequential steps detected: {sequential_steps}")

    # Determine the final result
    result = n_functionalization_steps >= 2 and sequential_steps

    # Add count constraint if met
    if n_functionalization_steps >= 2:
        if {"type": "count", "details": {"target": "N-functionalization", "operator": ">=", "value": 2}} not in findings_json["structural_constraints"]:
            findings_json["structural_constraints"].append({"type": "count", "details": {"target": "N-functionalization", "operator": ">=", "value": 2}})

    return result, findings_json
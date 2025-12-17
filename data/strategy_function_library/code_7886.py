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


N_METHYLATION_REACTIONS = [
    "N-methylation",
    "Methylation with MeI_primary",
    "Methylation with MeI_secondary",
    "Methylation with MeI_tertiary",
    "DMS Amine methylation",
    "Eschweiler-Clarke Primary Amine Methylation",
    "Eschweiler-Clarke Secondary Amine Methylation",
    "Reductive methylation of primary amine with formaldehyde",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a synthetic strategy involving at least two sequential N-methylation reactions on an imidazolidinedione (hydantoin) scaffold. The function identifies reactions belonging to the following classes: N-methylation, Methylation with MeI_primary, Methylation with MeI_secondary, Methylation with MeI_tertiary, DMS Amine methylation, Eschweiler-Clarke Primary Amine Methylation, Eschweiler-Clarke Secondary Amine Methylation, and Reductive methylation of primary amine with formaldehyde.
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

    # Track the number of N-methylation reactions
    n_methylation_count = 0
    # Flag to check if imidazolidinedione scaffold is present
    has_imidazolidinedione = False
    # Track molecules that have undergone N-methylation
    methylated_molecules = {}  # Store molecule SMILES and their methylation status

    def is_n_methylation(rxn_smiles):
        """Check if a reaction is one of several specific N-methylation types."""
        try:
            for reaction_name in N_METHYLATION_REACTIONS:
                if checker.check_reaction(reaction_name, rxn_smiles):
                    if reaction_name not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append(reaction_name)
                    return True
            return False
        except Exception as e:
            print(f"Error in is_n_methylation: {e}")
            return False

    def check_imidazolidinedione(smiles):
        """Check if a molecule contains an imidazolidinedione scaffold (hydantoin)"""
        try:
            mol = Chem.MolFromSmiles(smiles)
            if not mol:
                return False

            # SMARTS pattern for imidazolidinedione (hydantoin)
            pattern = Chem.MolFromSmarts("O=C1NC(=O)[C,N]N1")
            if mol.HasSubstructMatch(pattern):
                if "imidazolidinedione" not in findings_json["atomic_checks"]["ring_systems"]:
                    findings_json["atomic_checks"]["ring_systems"].append("imidazolidinedione")
                return True

            return False
        except Exception as e:
            print(f"Error in check_imidazolidinedione: {e}")
            return False

    def count_n_methyl_groups(smiles):
        """Count the number of N-methyl groups in the imidazolidinedione scaffold"""
        try:
            mol = Chem.MolFromSmiles(smiles)
            if not mol:
                return 0

            # Pattern for N-methyl in hydantoin at position 1
            pattern1 = Chem.MolFromSmarts("CN1C(=O)NC(=O)[C,N]1")
            # Pattern for N-methyl in hydantoin at position 3
            pattern2 = Chem.MolFromSmarts("C1NC(=O)N(C)C(=O)[C,N]1")
            # Pattern for di-N-methylated hydantoin
            pattern3 = Chem.MolFromSmarts("CN1C(=O)N(C)C(=O)[C,N]1")

            count = 0
            if mol.HasSubstructMatch(pattern1):
                count += 1
            if mol.HasSubstructMatch(pattern2):
                count += 1
            # If we have a di-methylated pattern, ensure we count it as 2
            if mol.HasSubstructMatch(pattern3):
                count = 2

            if count > 0 and "N-methyl group" not in findings_json["atomic_checks"]["functional_groups"]:
                findings_json["atomic_checks"]["functional_groups"].append("N-methyl group")

            return count
        except Exception as e:
            print(f"Error in count_n_methyl_groups: {e}")
            return 0

    def is_methylation_change(product_smiles, reactant_smiles):
        """Check if there's an increase in N-methyl groups from reactant to product"""
        try:
            product_methyl_count = count_n_methyl_groups(product_smiles)
            reactant_methyl_count = count_n_methyl_groups(reactant_smiles)

            print(
                f"Methylation change check: {reactant_smiles} ({reactant_methyl_count} methyl groups) -> {product_smiles} ({product_methyl_count} methyl groups)"
            )
            return product_methyl_count > reactant_methyl_count
        except Exception as e:
            print(f"Error in is_methylation_change: {e}")
            return False

    def dfs_traverse(node, depth=0, parent_smiles=None):
        nonlocal n_methylation_count, has_imidazolidinedione, findings_json

        if node["type"] == "mol":
            mol_smiles = node["smiles"]

            # Check if molecule contains imidazolidinedione scaffold
            if check_imidazolidinedione(mol_smiles):
                has_imidazolidinedione = True
                print(f"Found imidazolidinedione scaffold in molecule: {mol_smiles}")

                # Store the methylation status of this molecule
                methyl_count = count_n_methyl_groups(mol_smiles)
                methylated_molecules[mol_smiles] = methyl_count
                print(f"Molecule has {methyl_count} N-methyl groups on imidazolidinedione")

        elif node["type"] == "reaction":
            try:
                # Extract reactants and product
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if this is an N-methylation reaction
                if is_n_methylation(rsmi):
                    print(f"Found potential N-methylation reaction: {rsmi}")

                    # Check if the product contains imidazolidinedione scaffold
                    if check_imidazolidinedione(product):
                        # Check if any reactant also contains imidazolidinedione
                        for reactant in reactants:
                            if check_imidazolidinedione(reactant):
                                # Verify that N-methylation occurred on the imidazolidinedione
                                if is_methylation_change(product, reactant):
                                    n_methylation_count += 1
                                    print(
                                        f"Confirmed N-methylation on imidazolidinedione: {reactant} -> {product}"
                                    )

                                    # Store the product's methylation state for future reference
                                    methylated_molecules[product] = count_n_methyl_groups(product)

                                    # Check for sequential methylation
                                    if reactant in methylated_molecules:
                                        print(
                                            f"Sequential methylation detected: {reactant} ({methylated_molecules[reactant]} methyl groups) -> {product} ({count_n_methyl_groups(product)} methyl groups)"
                                        )
            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Determine the next depth based on the current node's type
        next_depth = depth
        if node["type"] != "reaction": # If current node is 'mol' (chemical), depth increases
            next_depth = depth + 1

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, next_depth, node.get("smiles", ""))

    # Start traversal
    dfs_traverse(route)

    # Strategy criteria: at least 2 N-methylation reactions and presence of imidazolidinedione scaffold
    strategy_detected = n_methylation_count >= 2 and has_imidazolidinedione
    print(f"Sequential N-methylation strategy detected: {strategy_detected}")
    print(f"N-methylation count: {n_methylation_count}")
    print(f"Imidazolidinedione scaffold present: {has_imidazolidinedione}")

    # Add structural constraints if detected
    if has_imidazolidinedione and n_methylation_count > 0:
        # Co-occurrence of imidazolidinedione and N-methylation
        if {"type": "co-occurrence", "details": {"targets": ["imidazolidinedione", "N-methylation"]}} not in findings_json["structural_constraints"]:
            findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"targets": ["imidazolidinedione", "N-methylation"]}})

    if n_methylation_count >= 2:
        # Count constraint for N-methylation
        if {"type": "count", "details": {"target": "N-methylation", "operator": ">=", "value": 2}} not in findings_json["structural_constraints"]:
            findings_json["structural_constraints"].append({"type": "count", "details": {"target": "N-methylation", "operator": ">=", "value": 2}})
        # Sequence constraint for N-methylation
        if {"type": "sequence", "details": {"ordered_events": ["N-methylation", "N-methylation"]}} not in findings_json["structural_constraints"]:
            findings_json["structural_constraints"].append({"type": "sequence", "details": {"ordered_events": ["N-methylation", "N-methylation"]}})

    return strategy_detected, findings_json

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

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a synthetic strategy that uses at least two different types of C-N bond forming reactions across at least three distinct steps. The function classifies reactions into categories such as reductive amination, N-arylation, amide formation, N-heterocycle formation, urea/thiourea formation, amine alkylation, aza-Michael addition, and imine formation to identify this pattern.
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

    # Store tuples of (depth, reaction_type) for C-N bond formations
    cn_bond_formations = []

    def count_cn_bonds(mol_smiles):
        """Count the number of C-N bonds in a molecule"""
        if not mol_smiles:
            return 0

        try:
            mol = Chem.MolFromSmiles(mol_smiles)
            if not mol:
                return 0

            count = 0
            for bond in mol.GetBonds():
                a1 = mol.GetAtomWithIdx(bond.GetBeginAtomIdx())
                a2 = mol.GetAtomWithIdx(bond.GetEndAtomIdx())
                if (a1.GetSymbol() == "C" and a2.GetSymbol() == "N") or (
                    a1.GetSymbol() == "N" and a2.GetSymbol() == "C"
                ):
                    count += 1
            return count
        except:
            return 0

    def dfs_traverse(node, depth=0):
        nonlocal cn_bond_formations, findings_json

        if node["type"] == "reaction":
            # Extract reaction information
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            print(f"Analyzing reaction at depth {depth}: {rsmi}")

            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Count C-N bonds before and after reaction
            reactants_cn_bonds = sum(count_cn_bonds(r) for r in reactants_smiles)
            product_cn_bonds = count_cn_bonds(product_smiles)

            # Check if new C-N bonds were formed
            new_cn_bonds = product_cn_bonds - reactants_cn_bonds
            print(
                f"  C-N bonds: reactants={reactants_cn_bonds}, product={product_cn_bonds}, new={new_cn_bonds}"
            )

            # Only proceed if new C-N bonds were formed
            if new_cn_bonds > 0:
                reaction_type = None

                # Check for C-N bond formation reactions
                if checker.check_reaction("reductive amination", rsmi):
                    print(f"  Reductive amination detected at depth {depth}")
                    reaction_type = "reductive_amination"
                    findings_json["atomic_checks"]["named_reactions"].append("reductive amination")

                # Check for Buchwald-Hartwig/N-arylation reactions
                elif checker.check_reaction(
                    "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)", rsmi
                ):
                    print(f"  Buchwald-Hartwig/N-arylation detected at depth {depth}")
                    reaction_type = "n_arylation"
                    findings_json["atomic_checks"]["named_reactions"].append("N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)")

                # Check for amide formation reactions
                elif checker.check_reaction(
                    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids", rsmi
                ):
                    print(f"  Amide formation detected at depth {depth}")
                    reaction_type = "amide_formation"
                    findings_json["atomic_checks"]["named_reactions"].append("Acylation of Nitrogen Nucleophiles by Carboxylic Acids")

                # Check for heterocycle formation involving N
                elif checker.check_reaction("Formation of NOS Heterocycles", rsmi):
                    print(f"  N-heterocycle formation detected at depth {depth}")
                    reaction_type = "n_heterocycle_formation"
                    findings_json["atomic_checks"]["named_reactions"].append("Formation of NOS Heterocycles")

                # Check for urea/thiourea formation
                elif (
                    checker.check_reaction("urea", rsmi)
                    or checker.check_reaction("thiourea", rsmi)
                ):
                    print(f"  Urea/thiourea formation detected at depth {depth}")
                    reaction_type = "urea_formation"
                    if checker.check_reaction("urea", rsmi):
                        findings_json["atomic_checks"]["named_reactions"].append("urea")
                    if checker.check_reaction("thiourea", rsmi):
                        findings_json["atomic_checks"]["named_reactions"].append("thiourea")

                # Check for amine alkylation
                elif checker.check_reaction("Alkylation of amines", rsmi):
                    print(f"  Amine alkylation detected at depth {depth}")
                    reaction_type = "amine_alkylation"
                    findings_json["atomic_checks"]["named_reactions"].append("Alkylation of amines")

                # Check for aza-Michael addition
                elif (
                    checker.check_reaction("aza-Michael addition aromatic", rsmi)
                    or checker.check_reaction("aza-Michael addition secondary", rsmi)
                    or checker.check_reaction("aza-Michael addition primary", rsmi)
                ):
                    print(f"  Aza-Michael addition detected at depth {depth}")
                    reaction_type = "aza_michael_addition"
                    if checker.check_reaction("aza-Michael addition aromatic", rsmi):
                        findings_json["atomic_checks"]["named_reactions"].append("aza-Michael addition aromatic")
                    if checker.check_reaction("aza-Michael addition secondary", rsmi):
                        findings_json["atomic_checks"]["named_reactions"].append("aza-Michael addition secondary")
                    if checker.check_reaction("aza-Michael addition primary", rsmi):
                        findings_json["atomic_checks"]["named_reactions"].append("aza-Michael addition primary")

                # Check for imine formation
                elif (
                    checker.check_reaction(
                        "Addition of primary amines to aldehydes/thiocarbonyls", rsmi
                    )
                    or checker.check_reaction(
                        "Addition of primary amines to ketones/thiocarbonyls", rsmi
                    )
                    or checker.check_reaction(
                        "Addition of secondary amines to ketones/thiocarbonyls", rsmi
                    )
                    or checker.check_reaction(
                        "Addition of secondary amines to aldehydes/thiocarbonyls", rsmi
                    )
                ):
                    print(f"  Imine formation detected at depth {depth}")
                    reaction_type = "imine_formation"
                    if checker.check_reaction("Addition of primary amines to aldehydes/thiocarbonyls", rsmi):
                        findings_json["atomic_checks"]["named_reactions"].append("Addition of primary amines to aldehydes/thiocarbonyls")
                    if checker.check_reaction("Addition of primary amines to ketones/thiocarbonyls", rsmi):
                        findings_json["atomic_checks"]["named_reactions"].append("Addition of primary amines to ketones/thiocarbonyls")
                    if checker.check_reaction("Addition of secondary amines to ketones/thiocarbonyls", rsmi):
                        findings_json["atomic_checks"]["named_reactions"].append("Addition of secondary amines to ketones/thiocarbonyls")
                    if checker.check_reaction("Addition of secondary amines to aldehydes/thiocarbonyls", rsmi):
                        findings_json["atomic_checks"]["named_reactions"].append("Addition of secondary amines to aldehydes/thiocarbonyls")

                # Check for ring opening of epoxide with amine
                elif checker.check_reaction("Ring opening of epoxide with amine", rsmi):
                    print(f"  Epoxide ring opening with amine detected at depth {depth}")
                    reaction_type = "epoxide_amine_opening"
                    findings_json["atomic_checks"]["named_reactions"].append("Ring opening of epoxide with amine")

                # Generic fallback - if we detected new C-N bonds but didn't match a specific reaction type
                elif new_cn_bonds > 0:
                    print(f"  Generic C-N bond formation detected at depth {depth}")
                    reaction_type = "other_cn_bond_formation"

                # If a C-N bond formation was detected, add it to our list
                if reaction_type:
                    cn_bond_formations.append((depth, reaction_type))
                    print(f"  Added C-N bond formation: {reaction_type} at depth {depth}")

        # Traverse children with new depth logic
        for child in node.get("children", []):
            new_depth = depth
            if node["type"] != "reaction":  # Only increase depth if current node is chemical
                new_depth = depth + 1
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    # Count unique reaction types and depths
    unique_reaction_types = set(reaction_type for _, reaction_type in cn_bond_formations)
    unique_depths = set(depth for depth, _ in cn_bond_formations)

    print(f"Found {len(cn_bond_formations)} C-N bond formations")
    print(f"Unique reaction types: {len(unique_reaction_types)}, {unique_reaction_types}")
    print(f"Unique depths: {len(unique_depths)}, {unique_depths}")

    # Determine the final result
    result = len(unique_depths) >= 3 and len(unique_reaction_types) >= 2

    # Add structural constraints to findings_json if met
    if len(unique_depths) >= 3:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "distinct steps with C-N bond formation",
                "operator": ">=",
                "value": 3
            }
        })
    if len(unique_reaction_types) >= 2:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "unique types of C-N bond forming reactions",
                "operator": ">=",
                "value": 2
            }
        })

    # Return True if at least 3 C-N bond formations are detected at different depths
    # AND we have at least 2 different types of C-N bond formations
    return result, findings_json

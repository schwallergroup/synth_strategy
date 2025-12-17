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
    This function detects the use of a benzyl ether (Bn) protecting group strategy. It operates by traversing a synthesis tree and identifying: 1. The formation of a benzyl ether via a named 'Williamson Ether Synthesis' reaction. 2. The cleavage of a benzyl ether via a named 'Hydroxyl benzyl deprotection' reaction or a generic transformation from a benzyl ether to an alcohol. 3. The presence of benzyl ether-containing intermediates across multiple steps. A strategy is flagged if a protection/deprotection sequence is found or if a protected intermediate persists for two or more steps.
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

    # Track molecules with benzyl ether protection and deprotection reactions
    protected_molecules = []
    protection_reactions = []
    deprotection_reactions = []

    # Track depth for each molecule to identify early vs late stage
    molecule_depths = {}

    def dfs_traverse(node, depth=0):
        nonlocal protected_molecules, protection_reactions, deprotection_reactions, molecule_depths, findings_json
        if node["type"] == "mol" and node["smiles"]:
            mol_smiles = node["smiles"]

            # Store molecule depth
            molecule_depths[mol_smiles] = depth

            # Check if molecule has benzyl ether
            if checker.check_fg("Ether", mol_smiles):
                findings_json["atomic_checks"]["functional_groups"].append("Ether")
                mol = Chem.MolFromSmiles(mol_smiles)
                # CORRECTED: SMARTS for benzyl ether, not benzyl alcohol.
                if mol:
                    benzyl_pattern = Chem.MolFromSmarts("[#6]O[CH2]c1ccccc1")
                    if mol.HasSubstructMatch(benzyl_pattern):
                        protected_molecules.append((mol_smiles, depth))
                        findings_json["atomic_checks"]["functional_groups"].append("benzyl_ether")
                        print(f"Found molecule with benzyl ether at depth {depth}: {mol_smiles}")

        elif node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rxn_smiles = node["metadata"]["mapped_reaction_smiles"]

            # Extract reactants and product
            try:
                reactants = rxn_smiles.split(">")[0].split(".")
                product = rxn_smiles.split(">")[-1]

                # Check for Williamson Ether Synthesis (often used for benzyl protection)
                # REMOVED: Inner 'if' was redundant and buggy. Trust the named reaction checker.
                if checker.check_reaction("Williamson Ether Synthesis", rxn_smiles):
                    protection_reactions.append((rxn_smiles, depth))
                    findings_json["atomic_checks"]["named_reactions"].append("Williamson Ether Synthesis")
                    print(f"Found Williamson ether protection at depth {depth}: {rxn_smiles}")

                # REMOVED: Generic protection check was flawed and a source of false positives.

                # Check for deprotection reaction (benzyl ether → alcohol)
                # CORRECTED: SMARTS for benzyl ether is fixed implicitly in the logic below.
                benzyl_ether_in_reactants = False
                for r in reactants:
                    if checker.check_fg("Ether", r):
                        mol_r = Chem.MolFromSmiles(r)
                        if mol_r and mol_r.HasSubstructMatch(Chem.MolFromSmarts("[#6]O[CH2]c1ccccc1")):
                            benzyl_ether_in_reactants = True
                            findings_json["atomic_checks"]["functional_groups"].append("benzyl_ether")
                            break

                alcohol_in_product = False
                if checker.check_fg("Primary alcohol", product):
                    alcohol_in_product = True
                    findings_json["atomic_checks"]["functional_groups"].append("Primary alcohol")
                if checker.check_fg("Secondary alcohol", product):
                    alcohol_in_product = True
                    findings_json["atomic_checks"]["functional_groups"].append("Secondary alcohol")
                if checker.check_fg("Tertiary alcohol", product):
                    alcohol_in_product = True
                    findings_json["atomic_checks"]["functional_groups"].append("Tertiary alcohol")

                benzyl_ether_in_product = False
                if checker.check_fg("Ether", product):
                    mol_p = Chem.MolFromSmiles(product)
                    if mol_p and mol_p.HasSubstructMatch(Chem.MolFromSmarts("[#6]O[CH2]c1ccccc1")):
                        benzyl_ether_in_product = True
                        findings_json["atomic_checks"]["functional_groups"].append("benzyl_ether")

                # Check for hydroxyl benzyl deprotection reaction specifically
                if checker.check_reaction("Hydroxyl benzyl deprotection", rxn_smiles):
                    deprotection_reactions.append((rxn_smiles, depth))
                    findings_json["atomic_checks"]["named_reactions"].append("Hydroxyl benzyl deprotection")
                    print(
                        f"Found hydroxyl benzyl deprotection reaction at depth {depth}: {rxn_smiles}"
                    )
                # Generic deprotection check
                elif (
                    benzyl_ether_in_reactants and alcohol_in_product and not benzyl_ether_in_product
                ):
                    deprotection_reactions.append((rxn_smiles, depth))
                    findings_json["atomic_checks"]["named_reactions"].append("benzyl_ether_deprotection")
                    print(f"Found deprotection reaction at depth {depth}: {rxn_smiles}")

            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Recursively traverse children
        for child in node.get("children", []):
            # New logic for depth calculation:
            # Depth increases only when traversing from a chemical node to a reaction node.
            # Depth remains the same when traversing from a reaction node to a chemical node.
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else:
                dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Analyze results to determine if benzyl ether protection strategy was used
    has_protection_strategy = False

    if protected_molecules:
        print(f"Found {len(protected_molecules)} molecules with benzyl ether protection")

        # Check if there's a pattern of early protection and late deprotection
        if protection_reactions and deprotection_reactions:
            protection_depths = [d for _, d in protection_reactions]
            deprotection_depths = [d for _, d in deprotection_reactions]

            print(f"Protection depths: {protection_depths}")
            print(f"Deprotection depths: {deprotection_depths}")

            # Protection should happen earlier (lower depth) than deprotection in retrosynthesis
            if min(protection_depths) < min(deprotection_depths):
                has_protection_strategy = True
                # Add structural constraint for co-occurrence and sequence
                if "Williamson Ether Synthesis" in findings_json["atomic_checks"]["named_reactions"] and \
                   "Hydroxyl benzyl deprotection" in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["structural_constraints"].append({
                        "type": "co-occurrence",
                        "details": {
                            "targets": [
                                "Williamson Ether Synthesis",
                                "Hydroxyl benzyl deprotection"
                            ]
                        }
                    })
                    findings_json["structural_constraints"].append({
                        "type": "sequence",
                        "details": {
                            "before": [
                                "Williamson Ether Synthesis"
                            ],
                            "after": [
                                "Hydroxyl benzyl deprotection"
                            ]
                        }
                    })
                elif "Williamson Ether Synthesis" in findings_json["atomic_checks"]["named_reactions"] and \
                     "benzyl_ether_deprotection" in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["structural_constraints"].append({
                        "type": "co-occurrence",
                        "details": {
                            "targets": [
                                "Williamson Ether Synthesis",
                                "benzyl_ether_deprotection"
                            ]
                        }
                    })
                    findings_json["structural_constraints"].append({
                        "type": "sequence",
                        "details": {
                            "before": [
                                "Williamson Ether Synthesis"
                            ],
                            "after": [
                                "benzyl_ether_deprotection"
                            ]
                        }
                    })
                print(
                    "Benzyl ether protection strategy detected: early protection and late deprotection"
                )

        # If we have protected molecules spanning multiple steps, that's also a protection strategy
        protected_depths = [d for _, d in protected_molecules]
        if len(protected_depths) > 1 and max(protected_depths) - min(protected_depths) >= 2:
            has_protection_strategy = True
            findings_json["structural_constraints"].append({
                "type": "count",
                "details": {
                    "target": "benzyl_ether_persistence_span",
                    "operator": ">=",
                    "value": 2
                }
            })
            print(
                f"Benzyl ether protection strategy detected: protection maintained across multiple steps (depth range: {min(protected_depths)}-{max(protected_depths)})"
            )

    return has_protection_strategy, findings_json
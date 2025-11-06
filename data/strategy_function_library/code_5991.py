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
    """Detects a multi-step strategy where a silyl-protected alkyne is used as an intermediate, and an alkyne is later converted to a ketone or aldehyde via a hydration reaction within the same synthetic lineage."""
    findings_template = {
        "atomic_checks": {
            "named_reactions": [],
            "ring_systems": [],
            "functional_groups": []
        },
        "structural_constraints": []
    }
    findings_json = copy.deepcopy(findings_template)

    # Track molecules and their transformations
    protected_alkynes = set()  # Set of silyl-protected alkyne SMILES
    enone_products = set()  # Set of ketone/aldehyde products from alkyne hydration

    # Track molecule lineage (reactant -> list of products) in retrosynthetic direction
    molecule_lineage = {}

    def dfs_traverse(node, depth=0):
        nonlocal findings_json
        if node["type"] == "mol":
            mol_smiles = node["smiles"]

            # Check for protected alkyne (TMS or other silyl groups)
            if checker.check_fg("Alkyne", mol_smiles):
                if "Alkyne" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Alkyne")
                
                if checker.check_fg("TMS ether protective group", mol_smiles):
                    protected_alkynes.add(mol_smiles)
                    if "TMS ether protective group" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("TMS ether protective group")
                
                if checker.check_fg("Silyl protective group", mol_smiles):
                    protected_alkynes.add(mol_smiles)
                    if "Silyl protective group" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Silyl protective group")

        elif node["type"] == "reaction":
            if "mapped_reaction_smiles" in node["metadata"]:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Store molecule lineage in retrosynthetic direction
                for reactant in reactants:
                    if reactant not in molecule_lineage:
                        molecule_lineage[reactant] = set([product])
                    else:
                        molecule_lineage[reactant].add(product)

                # Check for alkyne hydration reactions that form the target carbonyl
                if checker.check_reaction("Hydration of alkyne to ketone", rsmi):
                    enone_products.add(product)
                    if "Hydration of alkyne to ketone" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("Hydration of alkyne to ketone")
                
                if checker.check_reaction("Hydration of alkyne to aldehyde", rsmi):
                    enone_products.add(product)
                    if "Hydration of alkyne to aldehyde" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("Hydration of alkyne to aldehyde")

        # Determine the new depth based on the current node's type
        new_depth = depth
        if node["type"] == "mol": # From chemical to reaction, depth increases
            new_depth = depth + 1
        # If node['type'] == 'reaction', depth remains the same (reaction to chemical)

        # Traverse children with the determined new_depth
        for child in node.get("children", []):
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    # Check if we can trace a path from a protected alkyne to a carbonyl product
    strategy_found = False

    # Function to check if a source can form any of the target molecules
    def can_form(source, targets, visited=None):
        if visited is None:
            visited = set()

        if source in visited:
            return False

        visited.add(source)

        if source in targets:
            return True

        if source not in molecule_lineage:
            return False

        for product in molecule_lineage[source]:
            if can_form(product, targets, visited):
                return True

        return False

    # Check if any protected alkyne can form a carbonyl product
    for alkyne in protected_alkynes:
        if can_form(alkyne, enone_products):
            strategy_found = True
            # Add the structural constraint if the strategy is found
            if {
                "type": "sequence",
                "details": {
                    "before_events": [
                        "Alkyne",
                        "TMS ether protective group",
                        "Silyl protective group"
                    ],
                    "after_events": [
                        "Hydration of alkyne to ketone",
                        "Hydration of alkyne to aldehyde"
                    ],
                    "condition": "An intermediate molecule must be present that contains an 'Alkyne' AND ('TMS ether protective group' OR 'Silyl protective group'). This intermediate must be a synthetic precursor to the product of a subsequent reaction that is either 'Hydration of alkyne to ketone' OR 'Hydration of alkyne to aldehyde'."
                }
            } not in findings_json["structural_constraints"]:
                findings_json["structural_constraints"].append({
                    "type": "sequence",
                    "details": {
                        "before_events": [
                            "Alkyne",
                            "TMS ether protective group",
                            "Silyl protective group"
                        ],
                        "after_events": [
                            "Hydration of alkyne to ketone",
                            "Hydration of alkyne to aldehyde"
                        ],
                        "condition": "An intermediate molecule must be present that contains an 'Alkyne' AND ('TMS ether protective group' OR 'Silyl protective group'). This intermediate must be a synthetic precursor to the product of a subsequent reaction that is either 'Hydration of alkyne to ketone' OR 'Hydration of alkyne to aldehyde'."
                    }
                })
            break

    return strategy_found, findings_json

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


HETEROCYCLES_OF_INTEREST = [
    "pyrazole", "thiazole", "pyridine", "furan", "pyrrole", "imidazole",
    "oxazole", "pyrimidine", "pyrazine", "pyridazine", "triazole", "tetrazole",
    "indole", "quinoline", "isoquinoline", "benzimidazole", "benzoxazole",
    "benzothiazole", "isoxazole", "isothiazole", "oxadiazole", "thiadiazole",
]

FUNCTIONALIZATION_REACTIONS = [
    "Friedel-Crafts acylation", "Friedel-Crafts alkylation", "Aromatic halogenation",
    "Aromatic nitration", "N-alkylation of primary amines with alkyl halides",
    "N-alkylation of secondary amines with alkyl halides", "Suzuki coupling",
    "Heck reaction", "Sonogashira", "Buchwald-Hartwig", "N-arylation", "Minisci",
    "Aromatic fluorination", "Aromatic chlorination", "Aromatic bromination",
    "Aromatic iodination",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a strategy that functionalizes multiple heterocycles across a synthesis. A reaction is considered a functionalization if it matches a predefined list of named reactions in `FUNCTIONALIZATION_REACTIONS`. The strategy is triggered if at least two such reactions occur, modifying either at least two distinct heterocycle types (from `HETEROCYCLES_OF_INTEREST`) or two distinct molecular scaffolds containing these heterocycles.
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

    # Track functionalized heterocycles with molecule identity
    functionalized_heterocycles = set()
    functionalization_reactions = []

    # Track unique molecules containing heterocycles
    heterocycle_molecules = set()

    def dfs_traverse(node, depth=0):
        nonlocal functionalized_heterocycles, functionalization_reactions, heterocycle_molecules, findings_json

        if node["type"] == "mol":
            # Skip processing for molecule nodes
            pass

        elif node["type"] == "reaction":
            # Extract reactants and product
            if "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for heterocycles in reactants and product
                reactant_heterocycles = {}
                for reactant in reactants:
                    for heterocycle in HETEROCYCLES_OF_INTEREST:
                        if checker.check_ring(heterocycle, reactant):
                            if heterocycle not in findings_json["atomic_checks"]["ring_systems"]:
                                findings_json["atomic_checks"]["ring_systems"].append(heterocycle)
                            if reactant not in reactant_heterocycles:
                                reactant_heterocycles[reactant] = []
                            reactant_heterocycles[reactant].append(heterocycle)
                            # Track unique molecule-heterocycle combinations
                            heterocycle_molecules.add((reactant, heterocycle))

                product_heterocycles = []
                for heterocycle in HETEROCYCLES_OF_INTEREST:
                    if checker.check_ring(heterocycle, product):
                        if heterocycle not in findings_json["atomic_checks"]["ring_systems"]:
                            findings_json["atomic_checks"]["ring_systems"].append(heterocycle)
                        product_heterocycles.append(heterocycle)

                # Check if this is a functionalization reaction
                is_functionalization = False

                # Check if reaction is a known functionalization type
                for rxn_type in FUNCTIONALIZATION_REACTIONS:
                    if checker.check_reaction(rxn_type, rsmi):
                        is_functionalization = True
                        if rxn_type not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(rxn_type)
                        break

                if is_functionalization:
                    # Add heterocycles that are preserved in the product and functionalized
                    for reactant, heterocycles in reactant_heterocycles.items():
                        for heterocycle in heterocycles:
                            if heterocycle in product_heterocycles:
                                # Store the specific molecule-heterocycle pair
                                functionalized_heterocycles.add((reactant, heterocycle))

                    functionalization_reactions.append(rsmi)

        # Determine the new depth for the recursive call
        new_depth = depth
        if node["type"] != "reaction": # If current node is 'mol' (chemical), depth increases
            new_depth = depth + 1

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    # Get unique heterocycle types that have been functionalized
    unique_heterocycle_types = set([h[1] for h in functionalized_heterocycles])

    # Strategy is present if we have multiple heterocycle types and functionalization steps
    result = False

    # Check structural constraints and update findings_json
    num_functionalization_reactions = len(functionalization_reactions)
    num_unique_heterocycle_types = len(unique_heterocycle_types)
    num_functionalization_events_on_heterocycles = len(functionalized_heterocycles)
    unique_molecules = set([h[0] for h in functionalized_heterocycles])
    num_unique_molecules_with_functionalized_heterocycles = len(unique_molecules)

    # Constraint 1: functionalization_reactions >= 2
    if num_functionalization_reactions >= 2:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "functionalization_reactions",
                "operator": ">=",
                "value": 2
            }
        })

    # Constraint 2: unique_functionalized_heterocycle_types >= 2
    if num_unique_heterocycle_types >= 2:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "unique_functionalized_heterocycle_types",
                "operator": ">=",
                "value": 2
            }
        })

    # Constraint 3: functionalization_events_on_heterocycles >= 2 (This is implicitly covered by the logic below)
    # Constraint 4: unique_molecules_with_functionalized_heterocycles >= 2
    if num_unique_molecules_with_functionalized_heterocycles >= 2:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "unique_molecules_with_functionalized_heterocycles",
                "operator": ">=",
                "value": 2
            }
        })

    # Original logic for 'result'
    if num_unique_heterocycle_types >= 2 and num_functionalization_reactions >= 2:
        result = True
    elif num_functionalization_events_on_heterocycles >= 2 and num_functionalization_reactions >= 2:
        if num_unique_molecules_with_functionalized_heterocycles >= 2:
            result = True

    return result, findings_json

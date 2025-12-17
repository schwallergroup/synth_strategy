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


NITRILE_DERIVED_HETEROCYCLES = [
    "tetrazole",
    "triazole",
    "oxadiazole",
    "thiadiazole",
    "imidazole",
    "pyrazole",
    "isoxazole",
    "isothiazole",
]

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects a synthetic strategy where a nitrile group is introduced
    in an earlier step and later used as a precursor for heterocycle formation.
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

    # Track molecules containing nitriles and their depths
    nitrile_molecules = {}
    nitrile_to_heterocycle_depth = None

    def dfs_traverse(node, depth=0):
        nonlocal nitrile_to_heterocycle_depth, findings_json

        # For molecule nodes, check and track nitriles
        if node["type"] == "mol":
            mol_smiles = node["smiles"]
            if checker.check_fg("Nitrile", mol_smiles):
                print(f"Found molecule with nitrile at depth {depth}: {mol_smiles}")
                nitrile_molecules[mol_smiles] = depth
                if "Nitrile" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Nitrile")

                # Check if this is a starting material
                if node.get("in_stock", False):
                    print(f"Nitrile found in starting material: {mol_smiles}")

        # For reaction nodes, check for nitrile to heterocycle transformation
        elif node["type"] == "reaction":
            try:
                # Extract reactants and product from reaction SMILES
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check if any reactant contains a nitrile
                nitrile_in_reactants = any(
                    checker.check_fg("Nitrile", r) for r in reactants_smiles if r
                )

                if nitrile_in_reactants:
                    # Check for specific reactions that convert nitriles to heterocycles
                    if (
                        checker.check_reaction("Azide-nitrile click cycloaddition to tetrazole", rsmi)
                    ):
                        print(
                            f"Found nitrile to heterocycle transformation via click reaction at depth {depth}"
                        )
                        nitrile_to_heterocycle_depth = depth
                        if "Azide-nitrile click cycloaddition to tetrazole" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("Azide-nitrile click cycloaddition to tetrazole")
                    if (
                        checker.check_reaction("{tetrazole_terminal}", rsmi)
                    ):
                        print(
                            f"Found nitrile to heterocycle transformation via click reaction at depth {depth}"
                        )
                        nitrile_to_heterocycle_depth = depth
                        if "{tetrazole_terminal}" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("{tetrazole_terminal}")
                    if (
                        checker.check_reaction("{tetrazole_connect_regioisomere_1}", rsmi)
                    ):
                        print(
                            f"Found nitrile to heterocycle transformation via click reaction at depth {depth}"
                        )
                        nitrile_to_heterocycle_depth = depth
                        if "{tetrazole_connect_regioisomere_1}" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("{tetrazole_connect_regioisomere_1}")
                    if (
                        checker.check_reaction("{tetrazole_connect_regioisomere_2}", rsmi)
                    ):
                        print(
                            f"Found nitrile to heterocycle transformation via click reaction at depth {depth}"
                        )
                        nitrile_to_heterocycle_depth = depth
                        if "{tetrazole_connect_regioisomere_2}" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("{tetrazole_connect_regioisomere_2}")
                    
                    # Check if any heterocycle is formed in the product
                    for ring in NITRILE_DERIVED_HETEROCYCLES:
                        if checker.check_ring(ring, product_smiles):
                            # Verify the ring wasn't already in the reactants
                            if not any(
                                checker.check_ring(ring, r) for r in reactants_smiles if r
                            ):
                                print(
                                    f"Found nitrile to {ring} transformation at depth {depth}"
                                )
                                nitrile_to_heterocycle_depth = depth
                                if ring not in findings_json["atomic_checks"]["ring_systems"]:
                                    findings_json["atomic_checks"]["ring_systems"].append(ring)
                                if "ring_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                                    findings_json["atomic_checks"]["named_reactions"].append("ring_formation")
                                break
            except Exception as e:
                print(f"Error processing reaction at depth {depth}: {e}")

        # Recursively process children
        for child in node.get("children", []):
            if node["type"] == "reaction":
                # Depth remains the same when traversing from a reaction node to a chemical node
                dfs_traverse(child, depth)
            else:
                # Depth increases when traversing from a chemical node to a reaction node
                dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    print(f"Nitrile molecules found: {list(nitrile_molecules.keys())}")
    print(f"Nitrile to heterocycle transformation depth: {nitrile_to_heterocycle_depth}")

    # Strategy is present if we found both nitrile-containing molecules and a transformation to heterocycle
    strategy_present = bool(nitrile_molecules and nitrile_to_heterocycle_depth is not None)
    
    if strategy_present:
        # Add the structural constraint if the strategy is found
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "Nitrile",
                    "ring_formation"
                ],
                "description": "A nitrile functional group must be present in a molecule within the route, and a subsequent reaction must use a nitrile-containing reactant to form a new heterocycle (e.g., tetrazole, triazole, etc.)."
            }
        })

    print(f"Strategy present: {strategy_present}")
    return strategy_present, findings_json
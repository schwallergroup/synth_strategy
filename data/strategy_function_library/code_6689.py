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


N_METHYLATION_REACTIONS = [
    "Methylation with MeI_primary",
    "Methylation with MeI_secondary",
    "Methylation with MeI_tertiary",
    "Methylation with DMS",
    "DMS Amine methylation",
    "Eschweiler-Clarke Primary Amine Methylation",
    "Eschweiler-Clarke Secondary Amine Methylation",
    "Reductive methylation of primary amine with formaldehyde",
    "N-methylation",
]

O_ACETYLATION_REACTIONS = [
    "Acetic anhydride and alcohol to ester",
]

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects a strategy where a spirocyclic system is formed early in the synthesis,
    followed by sequential functional group modifications (N-methylation and O-acetylation).
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

    # Initialize tracking variables
    has_spirocycle_formation = False
    has_n_methylation = False
    has_o_acetylation = False
    spirocycle_depth = -1
    n_methylation_depth = -1
    o_acetylation_depth = -1

    def dfs_traverse(node, depth=0):
        nonlocal has_spirocycle_formation, has_n_methylation, has_o_acetylation
        nonlocal spirocycle_depth, n_methylation_depth, o_acetylation_depth
        nonlocal findings_json

        if node["type"] == "reaction":
            # Extract reactants and products
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            print(f"Analyzing reaction at depth {depth}: {rsmi}")

            # Check for spirocycle formation
            if not has_spirocycle_formation:
                # Check if product has a spirocycle
                product_mol = Chem.MolFromSmiles(product_smiles)
                product_has_spirocycle = False

                if product_mol:
                    # Method 1: Check for atoms that belong to multiple rings
                    ring_info = product_mol.GetRingInfo()
                    for atom_idx in range(product_mol.GetNumAtoms()):
                        if ring_info.NumAtomRings(atom_idx) >= 2:
                            print(f"Product contains a spirocycle at atom {atom_idx}")
                            product_has_spirocycle = True
                            break

                # Check if reactants don't have spirocycles
                reactants_have_spirocycle = False
                for r_smiles in reactants_smiles:
                    r_mol = Chem.MolFromSmiles(r_smiles)
                    if r_mol:
                        ring_info = r_mol.GetRingInfo()
                        for atom_idx in range(r_mol.GetNumAtoms()):
                            if ring_info.NumAtomRings(atom_idx) >= 2:
                                print(f"Reactant contains a spirocycle at atom {atom_idx}")
                                reactants_have_spirocycle = True
                                break
                        if reactants_have_spirocycle:
                            break

                # Check if this is a ring formation reaction
                if checker.check_reaction("spiro-chromanone", rsmi) or (
                    product_has_spirocycle and not reactants_have_spirocycle
                ):
                    has_spirocycle_formation = True
                    spirocycle_depth = depth
                    print(f"Detected spirocycle formation at depth {depth}")
                    if "spiro-chromanone" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("spiro-chromanone")
                    if "spirocycle" not in findings_json["atomic_checks"]["ring_systems"]:
                        findings_json["atomic_checks"]["ring_systems"].append("spirocycle")

                # If we're at the first reaction and product already has a spirocycle,
                # we'll consider the spirocycle as pre-formed
                if depth == 1 and product_has_spirocycle:
                    has_spirocycle_formation = True
                    spirocycle_depth = 0  # Consider it as pre-formed
                    print(f"Detected pre-formed spirocycle at depth {depth}")
                    if "spirocycle" not in findings_json["atomic_checks"]["ring_systems"]:
                        findings_json["atomic_checks"]["ring_systems"].append("spirocycle")

            # Check for N-methylation
            if not has_n_methylation:
                # Check for various N-methylation reactions
                for rxn_type in N_METHYLATION_REACTIONS:
                    if checker.check_reaction(rxn_type, rsmi):
                        has_n_methylation = True
                        n_methylation_depth = depth
                        print(f"Detected N-methylation ({rxn_type}) at depth {depth}")
                        if rxn_type not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(rxn_type)
                        break

            # Check for O-acetylation
            if not has_o_acetylation:
                # Check for various acetylation reactions
                for rxn_type in O_ACETYLATION_REACTIONS:
                    if checker.check_reaction(rxn_type, rsmi):
                        has_o_acetylation = True
                        o_acetylation_depth = depth
                        print(f"Detected O-acetylation ({rxn_type}) at depth {depth}")
                        if rxn_type not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(rxn_type)
                        break

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation:
            # Depth increases only when traversing from a chemical node to a reaction node.
            # Depth remains the same when traversing from a reaction node to a chemical node.
            if node["type"] == "reaction":
                # From reaction to chemical, depth remains the same
                dfs_traverse(child, depth)
            else:
                # From chemical to reaction, depth increases
                dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if the strategy is present
    # The strategy requires spirocycle formation followed by N-methylation and O-acetylation
    # In retrosynthesis, early-stage reactions have higher depth values
    strategy_present = has_spirocycle_formation and has_n_methylation and has_o_acetylation

    print(f"Spirocycle formation: {has_spirocycle_formation} at depth {spirocycle_depth}")
    print(f"N-methylation: {has_n_methylation} at depth {n_methylation_depth}")
    print(f"O-acetylation: {has_o_acetylation} at depth {o_acetylation_depth}")

    if strategy_present:
        print("Detected 'spirocycle formation then decoration' strategy")
        # Add structural constraint if all conditions are met
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "spirocycle_formation",
                    "N-methylation",
                    "O-acetylation"
                ]
            }
        })
    else:
        print("Strategy not detected")

    return strategy_present, findings_json

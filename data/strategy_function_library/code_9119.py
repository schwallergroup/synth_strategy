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
    Identifies syntheses that feature an ester hydrolysis. The strategy is flagged if the synthesis is linear (at most one convergent step), or if the ester hydrolysis is a late-stage step (occurring in the final half of the synthetic sequence), regardless of the route's convergence.
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

    has_ester_hydrolysis = False
    ester_hydrolysis_depth = -1
    convergent_count = 0
    max_depth = 0

    def is_complex_organic_molecule(smiles):
        """Check if molecule is a complex organic molecule rather than a simple reagent"""
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False
        # Complex molecules typically have more than 5 heavy atoms
        if mol.GetNumHeavyAtoms() < 6:
            return False
        return True

    def dfs_traverse(node, depth=0):
        nonlocal has_ester_hydrolysis, convergent_count, ester_hydrolysis_depth, max_depth, findings_json

        max_depth = max(max_depth, depth)

        if node["type"] == "reaction":
            # Extract reactants and product from reaction SMILES
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_part = rsmi.split(">")[0]
                product_part = rsmi.split(">")[-1]
                reactants = reactants_part.split(".")

                # Count number of complex reactants to identify true convergent steps
                complex_reactants = [r for r in reactants if is_complex_organic_molecule(r)]
                if len(complex_reactants) > 1:
                    convergent_count += 1
                    print(
                        f"Found convergent step at depth {depth}, with {len(complex_reactants)} complex reactants"
                    )

                # Check for ester hydrolysis using the checker function
                if checker.check_reaction(
                    "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters", rsmi
                ):
                    # Verify that we have an ester in reactants and carboxylic acid in product
                    has_ester = any(checker.check_fg("Ester", reactant) for reactant in reactants)
                    has_carboxylic_acid = checker.check_fg("Carboxylic acid", product_part)

                    if has_ester:
                        if "Ester" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Ester")
                    if has_carboxylic_acid:
                        if "Carboxylic acid" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Carboxylic acid")

                    if has_ester and has_carboxylic_acid:
                        has_ester_hydrolysis = True
                        ester_hydrolysis_depth = depth
                        if "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters")
                        # Add the co-occurrence constraint if all parts are found
                        if {"type": "co-occurrence", "details": {"targets": ["Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters", "Ester", "Carboxylic acid"]}} not in findings_json["structural_constraints"]:
                            findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"targets": ["Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters", "Ester", "Carboxylic acid"]}})
                        print(f"Found ester hydrolysis at depth {depth}, rsmi: {rsmi}")

                # Alternative check for ester saponification (another form of ester hydrolysis)
                elif (
                    checker.check_reaction("Ester saponification (methyl deprotection)", rsmi)
                    or checker.check_reaction("Ester saponification (alkyl deprotection)", rsmi)
                    or checker.check_reaction("COOH ethyl deprotection", rsmi)
                ):
                    has_ester_hydrolysis = True
                    ester_hydrolysis_depth = depth
                    if checker.check_reaction("Ester saponification (methyl deprotection)", rsmi):
                        if "Ester saponification (methyl deprotection)" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("Ester saponification (methyl deprotection)")
                    if checker.check_reaction("Ester saponification (alkyl deprotection)", rsmi):
                        if "Ester saponification (alkyl deprotection)" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("Ester saponification (alkyl deprotection)")
                    if checker.check_reaction("COOH ethyl deprotection", rsmi):
                        if "COOH ethyl deprotection" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("COOH ethyl deprotection")
                    print(f"Found ester saponification at depth {depth}, rsmi: {rsmi}")

            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Traverse children
        for child in node.get("children", []):
            # Determine the new depth based on the current node's type
            new_depth = depth if node["type"] == "reaction" else depth + 1
            dfs_traverse(child, new_depth)

    # Start traversal from root
    dfs_traverse(route)

    # Determine if hydrolysis is late-stage (in first half of synthesis depth)
    is_late_stage_hydrolysis = False
    if has_ester_hydrolysis and ester_hydrolysis_depth >= 0:
        is_late_stage_hydrolysis = ester_hydrolysis_depth <= (max_depth / 2)
        if is_late_stage_hydrolysis:
            if {"type": "positional", "details": {"target": "ester_hydrolysis", "position": "late_stage"}} not in findings_json["structural_constraints"]:
                findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "ester_hydrolysis", "position": "late_stage"}})
        print(
            f"Ester hydrolysis at depth {ester_hydrolysis_depth} of max depth {max_depth}, is late-stage: {is_late_stage_hydrolysis}"
        )

    # Return True if:
    # 1. It has an ester hydrolysis step
    # 2. It's predominantly linear (<=1 convergent step)
    # OR it has a late-stage hydrolysis regardless of convergence
    result = has_ester_hydrolysis and (convergent_count <= 1 or is_late_stage_hydrolysis)

    if convergent_count <= 1:
        if {"type": "count", "details": {"target": "convergent_step", "operator": "<=", "value": 1}} not in findings_json["structural_constraints"]:
            findings_json["structural_constraints"].append({"type": "count", "details": {"target": "convergent_step", "operator": "<=", "value": 1}})

    print(
        f"Linear synthesis with ester hydrolysis: {result} (ester hydrolysis: {has_ester_hydrolysis}, convergent steps: {convergent_count}, late-stage: {is_late_stage_hydrolysis})"
    )
    return result, findings_json

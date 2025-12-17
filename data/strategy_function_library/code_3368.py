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
    Detects a linear synthesis strategy focused on sequential heteroatom bond formations
    rather than C-C bond formations or convergent steps.
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

    # Track bond formations
    c_c_bond_formations = 0
    heteroatom_bond_formations = 0
    convergent_steps = 0

    def dfs_traverse(node, depth=0):
        nonlocal c_c_bond_formations, heteroatom_bond_formations, convergent_steps, findings_json

        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_part = rsmi.split(">")[0]
                product_part = rsmi.split(">")[-1]
                reactants = reactants_part.split(".")

                # Check if this is a convergent step (multiple complex reactants)
                complex_reactants = 0
                for reactant in reactants:
                    # Count reactants that have significant complexity
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and (mol.GetNumAtoms() > 10 or len(Chem.GetSSSR(mol)) > 1):
                        complex_reactants += 1

                if complex_reactants > 1:
                    convergent_steps += 1
                    print(f"Found convergent step at depth {depth}")

                # Get reaction SMILES for checking
                rxn_smiles = node["metadata"]["smiles"] if "smiles" in node["metadata"] else rsmi

                # Check for C-C bond formation reactions
                c_c_rxns = [
                    "Suzuki",
                    "Heck",
                    "Negishi",
                    "Stille",
                    "Kumada cross-coupling",
                    "Hiyama-Denmark Coupling",
                    "Wittig",
                    "Grignard",
                    "Aldol condensation",
                    "Michael addition",
                    "Diels-Alder",
                    "Sonogashira",
                    "Catellani reaction ortho",
                    "Catellani reaction para",
                    "beta C(sp3) arylation",
                    "Friedel-Crafts alkylation",
                    "Knoevenagel Condensation",
                    "Julia Olefination",
                ]
                for rxn_type in c_c_rxns:
                    if checker.check_reaction(rxn_type, rxn_smiles):
                        c_c_bond_formations += 1
                        findings_json["atomic_checks"]["named_reactions"].append(rxn_type)
                        print(f"Found C-C bond formation at depth {depth}: {rxn_type}")
                        break

                # Check for heteroatom bond formation reactions
                heteroatom_rxns = [
                    "Esterification of Carboxylic Acids",
                    "Williamson Ether Synthesis",
                    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
                    "Reductive amination with aldehyde",
                    "Reductive amination with ketone",
                    "Schotten-Baumann to ester",
                    "Acyl chloride with primary amine to amide",
                    "Acyl chloride with secondary amine to amide",
                    "Carboxylic acid with primary amine to amide",
                    "Ester with primary amine to amide",
                    "Ester with secondary amine to amide",
                    "Urea synthesis via isocyanate and primary amine",
                    "Urea synthesis via isocyanate and secondary amine",
                    "Mitsunobu esterification",
                    "Mitsunobu aryl ether",
                    "Sulfonamide synthesis (Schotten-Baumann) primary amine",
                    "Sulfonamide synthesis (Schotten-Baumann) secondary amine",
                    "N-alkylation of primary amines with alkyl halides",
                    "N-alkylation of secondary amines with alkyl halides",
                    "S-alkylation of thiols",
                    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
                    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
                    "Chan-Lam amine",
                    "Chan-Lam alcohol",
                    "Huisgen alkyne-azide 1,3 dipolar cycloaddition",
                ]
                for rxn_type in heteroatom_rxns:
                    if checker.check_reaction(rxn_type, rxn_smiles):
                        heteroatom_bond_formations += 1
                        findings_json["atomic_checks"]["named_reactions"].append(rxn_type)
                        print(f"Found heteroatom bond formation at depth {depth}: {rxn_type}")
                        break

            except Exception as e:
                print(f"Error analyzing reaction at depth {depth}: {e}")

        # Determine the next depth based on the current node's type
        next_depth = depth
        if node["type"] != "reaction": # This means it's a chemical node
            next_depth = depth + 1

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, next_depth)

    # Start traversal
    dfs_traverse(route)

    # Check if this matches our strategy criteria
    total_reactions = c_c_bond_formations + heteroatom_bond_formations
    is_linear_heteroatom_focused = (
        heteroatom_bond_formations >= 2
        and heteroatom_bond_formations > c_c_bond_formations
        and convergent_steps == 0
        and total_reactions > 0  # Ensure we have at least some reactions
    )

    print(
        f"Summary: {heteroatom_bond_formations} heteroatom formations, {c_c_bond_formations} C-C formations, {convergent_steps} convergent steps"
    )
    if is_linear_heteroatom_focused:
        print("Detected linear synthesis strategy focused on heteroatom bond formations")

    # Populate structural constraints based on the final evaluation
    if heteroatom_bond_formations >= 2:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "heteroatom_bond_formation",
                "operator": ">=",
                "value": 2
            }
        })
    if heteroatom_bond_formations > c_c_bond_formations:
        findings_json["structural_constraints"].append({
            "type": "relative_count",
            "details": {
                "target_a": "heteroatom_bond_formation",
                "target_b": "c_c_bond_formation",
                "operator": ">"
            }
        })
    if convergent_steps == 0:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "convergent_step",
                "operator": "==",
                "value": 0
            }
        })

    return is_linear_heteroatom_focused, findings_json

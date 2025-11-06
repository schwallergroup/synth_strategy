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
    Detects if the synthesis route involves a convergent fragment coupling strategy
    where two complex fragments are joined via C-N bond formation.
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

    convergent_coupling_found = False

    # C-N bond forming reactions to check
    cn_bond_forming_reactions = [
        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
        "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
        "Reductive amination with aldehyde",
        "Reductive amination with ketone",
        "Reductive amination with alcohol",
        "Alkylation of amines",
        "N-alkylation of primary amines with alkyl halides",
        "N-alkylation of secondary amines with alkyl halides",
        "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
        "Ugi reaction",
        "Goldberg coupling",
        "Ullmann-Goldberg Substitution amine",
        "Acylation of primary amines",
        "Acylation of secondary amines",
        "Acyl chloride with primary amine to amide (Schotten-Baumann)",
        "Acyl chloride with secondary amine to amide",
        "Carboxylic acid with primary amine to amide",
        "Ester with primary amine to amide",
        "Ester with secondary amine to amide",
        "Schotten-Baumann_amide",
    ]

    def dfs_traverse(node, depth=0):
        nonlocal convergent_coupling_found, findings_json

        print(f"Traversing node at depth {depth}, type: {node['type']}")

        if node["type"] == "reaction" and depth <= 3:  # Focus on late to mid-stage reactions
            print(f"Examining reaction at depth {depth}")
            # Add positional constraint if met
            findings_json["structural_constraints"].append({
                "type": "positional",
                "details": {
                    "target": "convergent_coupling_reaction",
                    "description": "The reaction must occur within the last 4 steps of the synthesis (depth <= 3)."
                }
            })

            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                print(f"Reaction SMILES: {rsmi}")

                reactants_part = rsmi.split(">")[0]
                product = rsmi.split(">")[-1]

                # Count number of reactant fragments
                reactants = reactants_part.split(".")
                print(f"Found {len(reactants)} reactant fragments")

                if len(reactants) >= 2:  # Multiple fragments being joined
                    # Add count constraint for reactants if met
                    findings_json["structural_constraints"].append({
                        "type": "count",
                        "details": {
                            "target": "reactants",
                            "scope": "convergent_coupling_reaction",
                            "operator": ">=",
                            "value": 2
                        }
                    })
                    try:
                        # Check if both reactants are complex
                        reactant_mols = [Chem.MolFromSmiles(r) for r in reactants]
                        product_mol = Chem.MolFromSmiles(product)

                        if all(reactant_mols) and product_mol:
                            # Check complexity based on atoms, rings, and functional groups
                            complex_fragments = []
                            for mol in reactant_mols:
                                mol_smiles = Chem.MolToSmiles(mol)
                                atom_count = mol.GetNumAtoms()
                                ring_count = len(Chem.GetSSSR(mol))

                                # Consider a fragment complex if it has >5 atoms and contains a ring
                                # or has >8 atoms regardless of rings
                                if (atom_count > 5 and ring_count > 0) or atom_count > 8:
                                    complex_fragments.append(mol)
                                    print(
                                        f"Complex fragment found: {mol_smiles}, atoms: {atom_count}, rings: {ring_count}"
                                    )

                            if len(complex_fragments) >= 2:
                                print(f"Found at least 2 complex fragments at depth {depth}")
                                # Add count constraint for complex reactants if met
                                findings_json["structural_constraints"].append({
                                    "type": "count",
                                    "details": {
                                        "target": "complex_reactants",
                                        "scope": "convergent_coupling_reaction",
                                        "operator": ">=",
                                        "value": 2,
                                        "definition": "A reactant is considered complex if it has more than 5 atoms and at least one ring, or if it has more than 8 atoms."
                                    }
                                })

                                # Check for C-N bond forming reactions
                                for reaction_type in cn_bond_forming_reactions:
                                    if checker.check_reaction(reaction_type, rsmi):
                                        print(f"Found C-N bond forming reaction: {reaction_type}")
                                        convergent_coupling_found = True
                                        findings_json["atomic_checks"]["named_reactions"].append(reaction_type)
                                        return

                    except Exception as e:
                        print(f"Error processing SMILES in convergent coupling detection: {e}")

        # Continue traversal
        for child in node.get("children", []):
            if (
                not convergent_coupling_found
            ):  # Stop traversal if we already found what we're looking for
                # New logic for depth calculation
                new_depth = depth
                if node['type'] != 'reaction': # If current node is chemical, depth increases
                    new_depth = depth + 1
                # If current node is reaction, depth remains the same
                dfs_traverse(child, new_depth)

    dfs_traverse(route)
    print(f"Convergent coupling strategy found: {convergent_coupling_found}")
    return convergent_coupling_found, findings_json

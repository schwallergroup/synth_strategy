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


AMIDE_COUPLING_REACTIONS = [
    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
    "Carboxylic acid with primary amine to amide",
    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
    "Acyl chloride with secondary amine to amide",
    "Ester with primary amine to amide",
    "Ester with secondary amine to amide",
    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
    "Schotten-Baumann_amide",
    "Acylation of primary amines",
    "Acylation of secondary amines",
]

ESTER_HYDROLYSIS_REACTIONS = [
    "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters",
    "Ester saponification (methyl deprotection)",
    "Ester saponification (alkyl deprotection)",
    "COOH ethyl deprotection",
    "Deprotection of carboxylic acid",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a convergent synthesis featuring a late-stage amide coupling where one fragment is prepared via ester hydrolysis. This is verified by checking for specific, enumerated amide coupling and ester hydrolysis reaction types and confirming that the carboxylic acid produced by the hydrolysis is a reactant in the subsequent coupling step. The function checks for the following amide coupling reactions: Acylation of Nitrogen Nucleophiles by Carboxylic Acids, Carboxylic acid with primary amine to amide, Acyl chloride with primary amine to amide (Schotten-Baumann), Acyl chloride with secondary amine to amide, Ester with primary amine to amide, Ester with secondary amine to amide, Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N, Schotten-Baumann_amide, Acylation of primary amines, Acylation of secondary amines. It also checks for the following ester hydrolysis reactions: Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters, Ester saponification (methyl deprotection), Ester saponification (alkyl deprotection), COOH ethyl deprotection, Deprotection of carboxylic acid.
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

    # Track if we found the key features
    found_amide_coupling = False
    found_ester_hydrolysis = False
    convergent_synthesis = False

    # Track molecular fragments at relevant depths (not reactions)
    fragments_at_relevant_depths = 0

    # Track the fragments to ensure one has ester hydrolysis
    fragments_with_ester_hydrolysis = 0

    # Track carboxylic acids from ester hydrolysis
    carboxylic_acids_from_hydrolysis = set()

    # Track reactants in amide coupling
    amide_coupling_reactants = set()

    def dfs_traverse(node, depth=0):
        nonlocal found_amide_coupling, found_ester_hydrolysis, fragments_at_relevant_depths
        nonlocal convergent_synthesis, fragments_with_ester_hydrolysis, findings_json

        # Count significant molecular fragments at depths 1-3
        if depth in [1, 2, 3] and node["type"] == "mol" and not node.get("in_stock", False):
            fragments_at_relevant_depths += 1

        # Check for reaction nodes
        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]

            # Check for amide coupling at depths 0-2 (late stage)
            if depth <= 2:
                for name in AMIDE_COUPLING_REACTIONS:
                    if checker.check_reaction(name, rsmi):
                        found_amide_coupling = True
                        if name not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(name)
                        try:
                            reactants = rsmi.split(">")[0].split(".")
                            for reactant in reactants:
                                amide_coupling_reactants.add(reactant)
                        except Exception as e:
                            print(f"Error extracting reactants from amide coupling: {e}")
                        break # Found one amide coupling reaction, no need to check others

            # Check for ester hydrolysis at depths 1-3
            if depth in [1, 2, 3]:
                for name in ESTER_HYDROLYSIS_REACTIONS:
                    if checker.check_reaction(name, rsmi):
                        found_ester_hydrolysis = True
                        fragments_with_ester_hydrolysis += 1
                        if name not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(name)
                        try:
                            product = rsmi.split(">")[-1]
                            if checker.check_fg("Carboxylic acid", product):
                                if "Carboxylic acid" not in findings_json["atomic_checks"]["functional_groups"]:
                                    findings_json["atomic_checks"]["functional_groups"].append("Carboxylic acid")
                                carboxylic_acids_from_hydrolysis.add(product)
                        except Exception as e:
                            print(f"Error extracting product from ester hydrolysis: {e}")
                        break # Found one ester hydrolysis reaction, no need to check others

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # Only increase depth if current node is not a reaction
                new_depth = depth + 1
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    # Check if we have a convergent synthesis (multiple significant fragments)
    if fragments_at_relevant_depths >= 2:
        convergent_synthesis = True
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "non-stock_molecule",
                "operator": ">=",
                "value": 2,
                "condition": "Molecule is at depth 1, 2, or 3 in the synthesis tree."
            }
        })

    # Check if any carboxylic acid from hydrolysis is used in amide coupling
    carboxylic_acid_used_in_coupling = False
    for acid in carboxylic_acids_from_hydrolysis:
        for reactant in amide_coupling_reactants:
            try:
                acid_mol = Chem.MolFromSmiles(acid)
                reactant_mol = Chem.MolFromSmiles(reactant)
                if acid_mol and reactant_mol:
                    acid_canonical = Chem.MolToSmiles(acid_mol)
                    reactant_canonical = Chem.MolToSmiles(reactant_mol)
                    if acid_canonical == reactant_canonical:
                        carboxylic_acid_used_in_coupling = True
                        break
            except Exception as e:
                print(f"Error comparing SMILES: {e}")
        if carboxylic_acid_used_in_coupling:
            break

    result = False
    if convergent_synthesis and found_amide_coupling and found_ester_hydrolysis and carboxylic_acid_used_in_coupling:
        result = True
        findings_json["structural_constraints"].append({
            "type": "sequence",
            "details": {
                "before": {
                    "event_type": "reaction_group",
                    "name": "ester_hydrolysis_reactions",
                    "position": "depth_1_to_3"
                },
                "after": {
                    "event_type": "reaction_group",
                    "name": "amide_coupling_reactions",
                    "position": "depth_0_to_2"
                },
                "link": {
                    "type": "product_is_reactant",
                    "description": "The carboxylic acid product of an ester hydrolysis reaction must be a reactant in a subsequent amide coupling reaction."
                }
            }
        })

    return result, findings_json

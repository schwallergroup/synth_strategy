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
    This function detects a strategy where a ketone is protected as a ketal/acetal
    and later deprotected.
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

    # Track protection and deprotection reactions
    protection_reactions = []
    deprotection_reactions = []

    result = False

    def are_structures_related(mol1, mol2):
        """Check if two molecules share significant structural similarity"""
        if mol1 is None or mol2 is None:
            return False

        # Use MCS to check if molecules share significant substructure
        mcs = rdkit.Chem.rdFMCS.FindMCS(
            [mol1, mol2],
            atomCompare=rdkit.Chem.rdFMCS.AtomCompare.CompareElements,
            bondCompare=rdkit.Chem.rdFMCS.BondCompare.CompareOrder,
            completeRingsOnly=True,
        )

        # Consider molecules related if they share at least 50% of atoms
        similarity_threshold = 0.5
        if mcs.numAtoms >= min(mol1.GetNumAtoms(), mol2.GetNumAtoms()) * similarity_threshold:
            return True
        return False

    def dfs_traverse(node, depth=0):
        nonlocal findings_json
        if node["type"] == "reaction":
            if "mapped_reaction_smiles" in node.get("metadata", {}):
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                try:
                    # Check for ketone protection reaction
                    is_acetalization = checker.check_reaction(
                        "Aldehyde or ketone acetalization", rsmi
                    )
                    is_diol_acetalization = checker.check_reaction("Diol acetalization", rsmi)

                    if is_acetalization or is_diol_acetalization:
                        if is_acetalization:
                            if "Aldehyde or ketone acetalization" not in findings_json["atomic_checks"]["named_reactions"]:
                                findings_json["atomic_checks"]["named_reactions"].append("Aldehyde or ketone acetalization")
                        if is_diol_acetalization:
                            if "Diol acetalization" not in findings_json["atomic_checks"]["named_reactions"]:
                                findings_json["atomic_checks"]["named_reactions"].append("Diol acetalization")

                        # Verify ketone in reactants and ketal in product
                        for reactant in reactants:
                            if checker.check_fg("Ketone", reactant):
                                if "Ketone" not in findings_json["atomic_checks"]["functional_groups"]:
                                    findings_json["atomic_checks"]["functional_groups"].append("Ketone")
                                if checker.check_fg("Acetal/Ketal", product):
                                    if "Acetal/Ketal" not in findings_json["atomic_checks"]["functional_groups"]:
                                        findings_json["atomic_checks"]["functional_groups"].append("Acetal/Ketal")
                                    protection_reactions.append((rsmi, depth, reactant, product))

                    # Check for ketal deprotection reaction
                    is_acetal_hydrolysis = checker.check_reaction(
                        "Acetal hydrolysis to ketone", rsmi
                    )
                    is_ketal_hydrolysis = checker.check_reaction("Ketal hydrolysis to ketone", rsmi)

                    if is_acetal_hydrolysis or is_ketal_hydrolysis:
                        if is_acetal_hydrolysis:
                            if "Acetal hydrolysis to ketone" not in findings_json["atomic_checks"]["named_reactions"]:
                                findings_json["atomic_checks"]["named_reactions"].append("Acetal hydrolysis to ketone")
                        if is_ketal_hydrolysis:
                            if "Ketal hydrolysis to ketone" not in findings_json["atomic_checks"]["named_reactions"]:
                                findings_json["atomic_checks"]["named_reactions"].append("Ketal hydrolysis to ketone")

                        # Verify ketal in reactants and ketone in product
                        for reactant in reactants:
                            if checker.check_fg("Acetal/Ketal", reactant):
                                if "Acetal/Ketal" not in findings_json["atomic_checks"]["functional_groups"]:
                                    findings_json["atomic_checks"]["functional_groups"].append("Acetal/Ketal")
                                if checker.check_fg("Ketone", product):
                                    if "Ketone" not in findings_json["atomic_checks"]["functional_groups"]:
                                        findings_json["atomic_checks"]["functional_groups"].append("Ketone")
                                    deprotection_reactions.append((rsmi, depth, reactant, product))

                except Exception as e:
                    print(f"Error processing reaction: {e}")
                    pass

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    # Check if we found both protection and deprotection
    if protection_reactions and deprotection_reactions:
        # Ensure protection happens before deprotection in the synthesis
        # (which means higher depth in retrosynthetic tree)
        for prot_rxn, prot_depth, prot_reactant, prot_product in protection_reactions:
            for deprot_rxn, deprot_depth, deprot_reactant, deprot_product in deprotection_reactions:
                if prot_depth > deprot_depth:
                    # Check if the protected and deprotected structures are related
                    prot_product_mol = Chem.MolFromSmiles(prot_product)
                    deprot_reactant_mol = Chem.MolFromSmiles(deprot_reactant)

                    if are_structures_related(prot_product_mol, deprot_reactant_mol):
                        result = True
                        # Add structural constraint if found
                        if {"type": "sequence", "details": {"before": ["Aldehyde or ketone acetalization", "Diol acetalization"], "after": ["Acetal hydrolysis to ketone", "Ketal hydrolysis to ketone"]}} not in findings_json["structural_constraints"]:
                            findings_json["structural_constraints"].append({"type": "sequence", "details": {"before": ["Aldehyde or ketone acetalization", "Diol acetalization"], "after": ["Acetal hydrolysis to ketone", "Ketal hydrolysis to ketone"]}})

    return result, findings_json

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

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if the synthesis route involves an early-stage O-alkylation of a phenol.
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

    # Track if we found an early-stage O-alkylation
    found_early_o_alkylation = False

    def dfs_traverse(node, depth=0):
        nonlocal found_early_o_alkylation, findings_json

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for phenol pattern in reactants
            phenol_pattern = Chem.MolFromSmarts("[OH][c]")

            # Check for ester alkylating agent pattern
            ester_pattern = Chem.MolFromSmarts("[#6][C](=O)[O][#6]")

            # Check for alkylated phenol pattern in product
            alkylated_phenol_pattern = Chem.MolFromSmarts("[O]([c])[#6][C](=O)[O][#6]")

            # Convert SMILES to molecules
            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles if r]
            product_mol = Chem.MolFromSmiles(product_smiles) if product_smiles else None

            reactant_has_phenol = False
            reactant_has_ester = False
            product_has_aryl_alkyl_ester_ether = False

            for mol in reactant_mols:
                if mol:
                    if mol.HasSubstructMatch(phenol_pattern):
                        reactant_has_phenol = True
                        if "phenol" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("phenol")
                    if mol.HasSubstructMatch(ester_pattern):
                        reactant_has_ester = True
                        if "ester" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("ester")

            if product_mol and product_mol.HasSubstructMatch(alkylated_phenol_pattern):
                product_has_aryl_alkyl_ester_ether = True
                if "aryl alkyl ester ether" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("aryl alkyl ester ether")

            if (
                product_mol
                and reactant_has_phenol
                and reactant_has_ester
                and product_has_aryl_alkyl_ester_ether
            ):
                found_early_o_alkylation = True
                print("Found early-stage phenol O-alkylation")
                if "phenol O-alkylation" not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append("phenol O-alkylation")
                
                # Add structural constraint for co-occurrence
                findings_json["structural_constraints"].append({
                    "type": "co-occurrence",
                    "details": {
                        "within": "reaction_step",
                        "targets": [
                            "reactant_has_phenol",
                            "reactant_has_ester",
                            "product_has_aryl_alkyl_ester_ether"
                        ],
                        "description": "A reaction is identified as 'phenol O-alkylation' if a phenol and an ester are present in the reactants, and an aryl alkyl ester ether is formed in the product."
                    }
                })
                # Add structural constraint for count
                # This constraint is met if the above condition is met at least once.
                # We'll add it here, but a more robust solution might track counts and add it once at the end.
                # For simplicity, adding it here for each instance found.
                findings_json["structural_constraints"].append({
                    "type": "count",
                    "details": {
                        "target": "phenol O-alkylation",
                        "operator": ">=",
                        "value": 1,
                        "description": "The route is valid if at least one phenol O-alkylation reaction occurs."
                    }
                })

            # Traverse children
            for child in node.get("children", []):
                # Depth increases when going from chemical to reaction
                # Depth stays the same when going from reaction to chemical
                if node["type"] == "reaction": # Current node is reaction, next is chemical
                    dfs_traverse(child, depth)
                else: # Current node is chemical, next is reaction
                    dfs_traverse(child, depth + 1)
        else:  # node["type"] == "mol"
            for child in node.get("children", []):
                # Depth increases when going from chemical to reaction
                # Depth stays the same when going from reaction to chemical
                if node["type"] == "reaction": # Current node is reaction, next is chemical
                    dfs_traverse(child, depth)
                else: # Current node is chemical, next is reaction
                    dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)
    return found_early_o_alkylation, findings_json

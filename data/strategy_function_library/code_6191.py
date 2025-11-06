from typing import Tuple, Dict, List
import copy
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
    This function detects a synthetic strategy involving late-stage incorporation
    of a heterocyclic system (specifically thiazolidinedione).
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

    thiazolidinedione_added_late = False

    def dfs_traverse(node, depth=0):
        nonlocal thiazolidinedione_added_late, findings_json

        if (
            node["type"] == "reaction" and depth <= 1
        ):  # Only check reactions at depth 0 or 1 (late stage)
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for thiazolidinedione pattern
                thiazolidinedione_pattern = Chem.MolFromSmarts(
                    "[#6]1[#16][#6](=[#8])[#7][#6]1=[#8]"
                )

                # Convert SMILES to molecules
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
                product_mol = Chem.MolFromSmiles(product) if product else None

                if product_mol and any(reactant_mols):
                    # Check if any reactant has thiazolidinedione
                    has_thiazolidinedione_reactant = any(
                        mol.HasSubstructMatch(thiazolidinedione_pattern)
                        for mol in reactant_mols
                        if mol
                    )

                    # Check if product has thiazolidinedione
                    has_thiazolidinedione_product = (
                        product_mol.HasSubstructMatch(thiazolidinedione_pattern)
                        if product_mol
                        else False
                    )

                    if has_thiazolidinedione_product and not has_thiazolidinedione_reactant:
                        thiazolidinedione_added_late = True
                        findings_json["atomic_checks"]["ring_systems"].append("thiazolidinedione")
                        findings_json["atomic_checks"]["named_reactions"].append("ring_formation") # Implied by new ring system
                        findings_json["structural_constraints"].append({
                            "type": "positional",
                            "details": {
                                "target": "thiazolidinedione_formation",
                                "position": "late_stage",
                                "definition": "Reaction occurs at depth <= 1"
                            }
                        })

        # Traverse children with increased depth
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else:
                dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)
    return thiazolidinedione_added_late, findings_json

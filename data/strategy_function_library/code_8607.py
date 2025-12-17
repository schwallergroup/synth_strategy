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
    Detects convergent synthesis with late-stage etherification connecting two complex fragments.
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

    late_stage_etherification = False

    def dfs_traverse(node, depth=0):
        nonlocal late_stage_etherification, findings_json

        if node["type"] == "reaction" and depth <= 1:  # Late-stage reaction (depth 0 or 1)
            if "mapped_reaction_smiles" in node.get("metadata", {}):
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if we have at least 2 reactants (convergent)
                if len(reactants) >= 2:
                    # Record the structural constraint for reactant count
                    findings_json["structural_constraints"].append({
                        "type": "count",
                        "details": {
                            "target": "reactants_in_etherification",
                            "operator": ">=",
                            "value": 2
                        }
                    })

                    # Check for etherification (C-O-C formation)
                    prod_mol = Chem.MolFromSmiles(product)
                    reactant_mols = [Chem.MolFromSmiles(r) for r in reactants]

                    # Look for ether pattern in product
                    ether_patt = Chem.MolFromSmarts("[#6][#8][#6]")
                    if prod_mol and prod_mol.HasSubstructMatch(ether_patt):
                        # Record atomic check for ether functional group
                        if "ether" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("ether")

                        # Check if ether bond is newly formed
                        ether_count_prod = len(prod_mol.GetSubstructMatches(ether_patt))
                        ether_count_reactants = sum(
                            len(r.GetSubstructMatches(ether_patt)) if r else 0
                            for r in reactant_mols
                        )

                        if ether_count_prod > ether_count_reactants:
                            late_stage_etherification = True
                            # Record atomic check for etherification named reaction
                            if "etherification" not in findings_json["atomic_checks"]["named_reactions"]:
                                findings_json["atomic_checks"]["named_reactions"].append("etherification")
                            # Record the structural constraint for late-stage position
                            findings_json["structural_constraints"].append({
                                "type": "positional",
                                "details": {
                                    "target": "etherification",
                                    "position": "late_stage"
                                }
                            })

        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is not a reaction (e.g., chemical)
                new_depth = depth + 1
            dfs_traverse(child, new_depth)

    dfs_traverse(route)
    return late_stage_etherification, findings_json
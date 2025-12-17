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
    Detects if the synthesis route involves late-stage N-alkylation of a cyclic amine.
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

    n_alkylation_found = False

    def dfs_traverse(node, depth=0):
        nonlocal n_alkylation_found, findings_json

        if node["type"] == "reaction" and depth <= 1:  # Late stage (lower depth)
            if "mapped_reaction_smiles" in node.get("metadata", {}):
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for N-alkylation patterns
                cyclic_amine_pattern = Chem.MolFromSmarts("[#7;R]")  # Nitrogen in ring
                sec_amine_pattern = Chem.MolFromSmarts("[#7;H1]")  # Secondary amine
                tert_amine_pattern = Chem.MolFromSmarts("[#7;H0]")  # Tertiary amine

                sec_amine_found = False

                for reactant in reactants:
                    try:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol:
                            if mol.HasSubstructMatch(
                                cyclic_amine_pattern
                            ) and mol.HasSubstructMatch(sec_amine_pattern):
                                sec_amine_found = True
                                findings_json["atomic_checks"]["functional_groups"].append("secondary amine")
                                findings_json["atomic_checks"]["ring_systems"].append("cyclic amine")
                    except:
                        continue

                # Check if product is a tertiary amine
                try:
                    product_mol = Chem.MolFromSmiles(product)
                    if (
                        sec_amine_found
                        and product_mol
                        and product_mol.HasSubstructMatch(tert_amine_pattern)
                        and product_mol.HasSubstructMatch(cyclic_amine_pattern)
                    ):
                        print("Late-stage N-alkylation detected at depth", depth)
                        n_alkylation_found = True
                        findings_json["atomic_checks"]["functional_groups"].append("tertiary amine")
                        findings_json["atomic_checks"]["named_reactions"].append("N-alkylation")
                        findings_json["structural_constraints"].append({
                            "type": "positional",
                            "details": {
                                "target": "N-alkylation of a cyclic amine",
                                "position": "within_last_2_stages"
                            }
                        })
                except:
                    pass

        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # Only increase depth if current node is not a reaction (i.e., chemical)
                new_depth = depth + 1
            dfs_traverse(child, new_depth)

    dfs_traverse(route)
    return n_alkylation_found, findings_json

#!/bin/bash

# Run `bash split_multifile.sh` from this directory

FUNCTIONS="passing_functions.py"
ROOT_DATA="/home/dparm/steerable_retro/data"

echo '#!/bin/python

"""LM-defined function for strategy description."""
' > _header.txt

grep -RE "import " "$FUNCTIONS" | awk '{a[$0]=1}END{for(s in a){if (match(s, "^[f|i]")){print s}}}' >> _header.txt

cp _header.txt _header_checker.txt
echo "from steerable_retro.utils.check import Check
from steerable_retro.utils import fuzzy_dict, check

root_data=\"$ROOT_DATA\"

fg_args = {
    \"file_path\": f\"{root_data}/patterns/functional_groups.json\",
    \"value_field\": \"pattern\",
    \"key_field\": \"name\",
}
reaction_class_args = {
    \"file_path\": f\"{root_data}/patterns/smirks.json\",
    \"value_field\": \"smirks\",
    \"key_field\": \"name\",
}
ring_smiles_args = {
    \"file_path\": f\"{root_data}/patterns/chemical_rings_smiles.json\",
    \"value_field\": \"smiles\",
    \"key_field\": \"name\",
}
functional_groups = fuzzy_dict.FuzzyDict.from_json(**fg_args)
reaction_classes = fuzzy_dict.FuzzyDict.from_json(**reaction_class_args)
ring_smiles = fuzzy_dict.FuzzyDict.from_json(**ring_smiles_args)

checker = check.Check(
    fg_dict=functional_groups,
    reaction_dict=reaction_classes,
    ring_dict=ring_smiles
)

" >> _header_checker.txt

# Remove all the imports from the file
cp $FUNCTIONS temp_functions.py
sed -i '/^from\|^import/d' temp_functions.py

rm -rf raw_functions
mkdir -p raw_functions/
csplit -f "raw_functions/code_" temp_functions.py "/^def /" "{*}" &> _logs.txt

# For each file in functions, create a new file with header and raw_functions/f_i

mkdir -p py_functions/
for file in raw_functions/code_*; do
    filename=$(basename "$file")
    new_file="py_functions/$filename.py"
    # Change name of function to `main`
    sed -i '1s/^def [^(]*(/def main(/' "$file"
    # if `checker` is in the file, use header_checker, else normal header
    if grep -q "checker" "$file"; then
        cat _header_checker.txt "$file" > "$new_file"
    else
        cat _header.txt "$file" > "$new_file"
    fi
    rm raw_functions/"$filename"

    if ! grep -q "def main" "$new_file"; then
        rm $new_file # remove if no function defined here
    fi
done

# Prettify
isort --profile black py_functions/.
black py_functions/.

cp py_functions/* ../steerable_retro/lm_code/.

rm -rf py_functions/ raw_functions/ _header* temp_functions.py _logs.txt
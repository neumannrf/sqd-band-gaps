#!/bin/bash

# Input files
DFT_FILE="compare_pseudos"
DFTUV_FILE="compare_gaps"
OUTPUT="final_results.csv"
file='input'

# Extract pseudo list and material from input
list_pseudos=$(awk 'f && !/END PSEUDO/ {print} /BEGIN PSEUDO/ {f=1} /END PSEUDO/ {f=0}' "${file}" | awk '{print $1}') 
material=$(awk 'f && !/END MATERIAL/ {print} /BEGIN MATERIAL/ {f=1} /END MATERIAL/ {f=0}' "${file}")
first_kpoint=$(awk 'f && !/END KPOINTS_BANDS/ {print} /BEGIN KPOINTS_BANDS/ {f=1} /END KPOINTS_BANDS/ {f=0}' ${file} | tail -n 2 | head -n 1 | awk '{print $NF}' | tr -d '!')
second_kpoint=$(awk 'f && !/END KPOINTS_BANDS/ {print} /BEGIN KPOINTS_BANDS/ {f=1} /END KPOINTS_BANDS/ {f=0}' ${file} | tail -n 1 | awk '{print $NF}' | tr -d '!')

tag_k1="gap${first_kpoint}"
tag_k2="gap${second_kpoint}"

# CSV Header with ;
echo "pseudo;opt_alat;opt_c_over_a;DFT_indirect;DFT_gapG;DFT_gapL;DFTU_indirect;DFTU_gapG;DFTU_gapL;DFTV_indirect;DFTV_gapG;DFTV_gapL;DFTUV_indirect;DFTUV_gapG;DFTUV_gapL;basis_len;num_electrons" > "$OUTPUT"

# Step 1: Parse compare_pseudos
awk '
BEGIN { OFS=";" }
NR > 1 {
    print $1, $2, $3, $4, $6, $5
}
' "$DFT_FILE" | awk '{ $1=tolower($1); print }' | sort > tmp_file1.csv

# Step 2: Parse compare_gaps
awk '
function reset() {
    dftu_ind = dftu_gapG = dftu_gapL = ""
    dftv_ind = dftv_gapG = dftv_gapL = ""
    dftuv_ind = dftuv_gapG = dftuv_gapL = ""
}
BEGIN { OFS=";"; reset() }
{
    if ($1 ~ /\.[Uu][Pp][Ff]$/) {
        if (pseudo != "") {
            print pseudo, dftu_ind, dftu_gapG, dftu_gapL, dftv_ind, dftv_gapG, dftv_gapL, dftuv_ind, dftuv_gapG, dftuv_gapL
        }
        reset()
        pseudo = $1
    } else if ($1 == "dftu") {
        dftu_ind=$2; dftu_gapG=$4; dftu_gapL=$3
    } else if ($1 == "dftv") {
        dftv_ind=$2; dftv_gapG=$4; dftv_gapL=$3
    } else if ($1 == "dftuv") {
        dftuv_ind=$2; dftuv_gapG=$4; dftuv_gapL=$3
    }
}
END {
    if (pseudo != "") {
        print pseudo, dftu_ind, dftu_gapG, dftu_gapL, dftv_ind, dftv_gapG, dftv_gapL, dftuv_ind, dftuv_gapG, dftuv_gapL
    }
}
' "$DFTUV_FILE" | awk '{ $1=tolower($1); print }' | sort > tmp_file2.csv

# Step 3: Extract PAOFLOW values
echo "" > tmp_extra.csv
for p in ${list_pseudos}; do
    pseudo=$(echo $p | tr '[:upper:]' '[:lower:]')
    # # This works with this example: ZnO 
    # base1=$(basename "${p}" .UPF | sed 's/\.upf$//')
    # base2=$(echo ${p}| cut -d '.' -f 2)
    # out_file="${base1}/paoflow/${material}_${base2}_paoflow_output"

    # The correct version of the line below should be:
    base=$(basename "${p}" .UPF | sed 's/\.upf$//')
    out_file="${base}/paoflow/${material}_${base}_paoflow_output"
    if [[ -f "$out_file" ]]; then
        basis=$(grep -m1 "Length of the basis vector:" "${out_file}" | awk '{print $NF}')
        nelec=$(grep -m1 "Number of electrons:" "${out_file}"  | awk '{print $NF}')
        echo "$pseudo;$basis;$nelec"
    else
        echo "$pseudo;N/A;N/A"
    fi
done | sort > tmp_extra.csv

# Step 4: Merge files
join -t ";" -1 1 -2 1 tmp_file1.csv tmp_file2.csv > tmp_joined.csv
join -t ";" -1 1 -2 1 tmp_joined.csv tmp_extra.csv >> "$OUTPUT"

# Cleanup
rm tmp_file1.csv tmp_file2.csv tmp_extra.csv tmp_joined.csv

sed -i "s/gapG/${tag_k1}/g" $OUTPUT

sed -i "s/gapL/${tag_k2}/g" $OUTPUT

sed -i "s/\./,/g" $OUTPUT
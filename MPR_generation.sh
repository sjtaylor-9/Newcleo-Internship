output_directory=$1
neutron_spectra_dir=$2
endf_jeff_data_dir=$3

mkdir $output_directory
mkdir $output_directory"/LFR30"
mkdir $output_directory"/LFR200"
mkdir $output_directory"/LFR30/FA"
mkdir $output_directory"/LFR30/FA/single"
mkdir $output_directory"/LFR200/inner"
mkdir $output_directory"/LFR200/middle"
mkdir $output_directory"/LFR200/outer"

# Make the burnup folders in each enrichment zone folder
for enrichment in inner middle outer; do
    for burnup in BoL EoL; do
        mkdir $output_directory"/LFR200/"$enrichment/$burnup
    done
done

echo "The necessary directories have been created"

#Calculate the one-group cross-sections for the LFR30
Path directories that involve the 'OneDrive - Newcleo' will be split into separate strings due to the whitespaces; therefore, must be passed in "" or ''
python "onegroup_xs_LFR.py" \
    --reactor "LFR30" \
    --enrichment_zone "FA" \
    --path $output_directory"/LFR30/FA/single" \
    --newcleo_input "$neutron_spectra_dir" \
    --cross_section_data "$endf_jeff_data_dir" \
    --reaction_data "Nuclide_data" \
    --burnup "single"

# Calculate the one-group cross-sections for the LFR200 in each enrichment zone and the BoL and EoL burnup steps
for enrichment in inner middle outer; do
    for burnup in BoL EoL; do
        python "onegroup_xs_LFR.py" \
        --reactor "LFR200" \
        --enrichment_zone $enrichment \
        --path $output_directory"/LFR200/"$enrichment/$burnup \
        --newcleo_input "$neutron_spectra_dir" \
        --cross_section_data "$endf_jeff_data_dir" \
        --reaction_data "Nuclide_data" \
        --burnup $burnup
    done
done

# The parent and daughter nuclides remain constant across all neutron fluxes so running using LFR200 arbitrarily
python "Finding_parent_daughter_nuclides.py" \
        --reactor "LFR200" \
        --enrichment_zone "inner" \
        --path $output_directory \
        --cross_section_data $output_directory/"LFR200/inner/BoL" \
        --reaction_data "Nuclide_data" \
        --burnup "BoL"

python "Writing_to_MPR.py" \
        --reactor "LFR30" \
        --enrichment_zone "FA" \
        --path $output_directory"/LFR30/FA/single" \
        --input $output_directory \
        --reaction_data "Nuclide_data" \
        --cross_section_data $output_directory/"LFR30/FA"

for enrichment in inner middle outer; do
        python "Writing_to_MPR.py" \
        --reactor "LFR200" \
        --enrichment_zone $enrichment \
        --path $output_directory"/LFR200/"$enrichment \
        --cross_section_data $output_directory/"LFR200"/$enrichment \
        --reaction_data "Nuclide_data" \
        --input $output_directory
done
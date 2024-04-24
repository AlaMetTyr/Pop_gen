#!/bin/bash -e


###########################
# Setup SLURM Environment #
###########################
#SBATCH --job-name=rename
#SBATCH --time=1:00:00
#SBATCH --mem=500mb
#SBATCH --account=ga03488

# Define the directory containing the folders to rename
parent_dir="/nesi/nobackup/ga03488/Amy/FAW/"

# Iterate over each folder in the parent directory
for folder in "$parent_dir"/*; do
    # Check if the item is a directory
    if [ -d "$folder" ]; then
        echo "Processing folder: $folder"
        
        # Get the name of the first file in the folder
        first_file=$(find "$folder" -maxdepth 1 -type f | head -n 1)
        echo "First file in folder: $first_file"

        # Extract the desired section of the file name
        new_folder_name=$(basename "$first_file" | cut -d '_' -f 1-4)
        echo "New folder name: $new_folder_name"

        # Check if the folder name needs to be changed
        if [ "$(basename "$folder")" != "$new_folder_name" ]; then
            echo "Renaming folder from $(basename "$folder") to $new_folder_name"
            # Rename the folder
            mv "$folder" "$parent_dir/$new_folder_name"
        else
            echo "Folder name is already correct, no need to rename."
        fi
    fi
done

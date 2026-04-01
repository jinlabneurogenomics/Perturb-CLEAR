"""
Author: Masood A. Akram
Created: January 06, 2025
Purpose: This script is used to generate neuronal output windows from 9x LS images using Recut bounding boxes.
Contact: masoodakram@mednet.ucla.edu OR masood.ahmed.akram@gmail.com
"""

import argparse
import pandas as pd
import numpy as np
import tifffile as tiff
from pathlib import Path
import os
import re
from tqdm import tqdm

def swc_to_csv_multiple(input_directory, output_csv):
    combined_data = []
    for swc_file in Path(input_directory).glob('*.swc'):
        with open(swc_file, 'r') as file:
            lines = [line.strip() for line in file if not line.startswith("#")]
        for line in lines:
            parts = line.split(' ')
            index, swc_type, x, y, z, radius, parent = parts
            new_name = f"{swc_file.stem}_xyz-{x}-{y}-{z}"
            combined_data.append([swc_file.name, float(x), float(y), float(z), float(radius), new_name])
    df = pd.DataFrame(combined_data, columns=['SWC_File', 'X', 'Y', 'Z', 'Radius', 'New_Name'])
    df.to_csv(output_csv, index=False)
    print(f"SWC to CSV file is generated here {output_csv}")

def process_images(base_folder, tiff_folder_name, voxel_size, image_orientation):
    tiff_folder = os.path.join(base_folder, tiff_folder_name)
    #channel_suffix = tiff_folder_name.split('_')[-1]
    csv_path = os.path.join(base_folder, f'original_soma_locations.csv')
    cropped_images_folder = os.path.join(base_folder, f'output_windows')
    output_csv_path = os.path.join(base_folder, f'new_soma_locations.csv')

    df = pd.read_csv(csv_path)
    if not os.path.exists(cropped_images_folder):
        os.makedirs(cropped_images_folder)

    max_x, max_y, max_z = get_image_dimensions(tiff_folder)
    output_data = []

    for index, row in tqdm(df.iterrows(), total=df.shape[0], desc="Processing Images", unit="image"):
        z_center, x_center, y_center = int(row['Z']), int(row['X']), int(row['Y'])
        new_name = row['New_Name']
        specific_folder = os.path.join(cropped_images_folder, new_name)
        if not os.path.exists(specific_folder):
            os.makedirs(specific_folder)

        roi_size_x = int(750 / voxel_size[0]) # 852 --> X = 1200, 1207 = 1700, 1136 = 1600
        roi_size_y = int(600/ voxel_size[1]) # 852 --> X = 1200, 639 = 900, 710 = 1000
        roi_size_z = int(500 / voxel_size[2])  # 497 --> X = 700, 426 = 600
        roi_size_x_center = int(700 / voxel_size[0])
        roi_size_y_center = int(700 / voxel_size[1])
        roi_size_z_center = int(450 / voxel_size[2])
        roi_size_x_old = int(600 / voxel_size[0])
        roi_size_y_old = int(600 / voxel_size[1])
        roi_size_z_old = int(400 / voxel_size[2])
        
        if image_orientation == "up":
            x_min = max(0, x_center - roi_size_x // 2)
            x_max = min(max_x, x_center + roi_size_x // 2)
            y_min = max(0, y_center - roi_size_y // 3*2)
            y_max = min(max_y, y_center + roi_size_y // 3)
        elif image_orientation == "down":
            x_min = max(0, x_center - roi_size_x // 2)
            x_max = min(max_x, x_center + roi_size_x // 2)
            y_min = max(0, y_center - roi_size_y // 3)
            y_max = min(max_y, y_center + roi_size_y // 3*2)
        elif image_orientation == "left":
            x_min = max(0, x_center - roi_size_x // 3*2)
            x_max = min(max_x, x_center + roi_size_x // 3)
            y_min = max(0, y_center - roi_size_y // 2)
            y_max = min(max_y, y_center + roi_size_y // 2)
        elif image_orientation == "right":
            x_min = max(0, x_center - roi_size_x // 3)
            x_max = min(max_x, x_center + roi_size_x // 3*2)
            y_min = max(0, y_center - roi_size_y // 2)
            y_max = min(max_y, y_center + roi_size_y // 2)
        elif image_orientation == "center":
            x_min = max(0, x_center - roi_size_x_center // 2)
            x_max = min(max_x, x_center + roi_size_x_center // 2)
            y_min = max(0, y_center - roi_size_y_center // 2)
            y_max = min(max_y, y_center + roi_size_y_center // 2)
        elif image_orientation == "old":
            x_min = max(0, x_center - roi_size_x_old // 2)
            x_max = min(max_x, x_center + roi_size_x_old // 2)
            y_min = max(0, y_center - roi_size_y_old // 3*2)
            y_max = min(max_y, y_center + roi_size_y_old // 3)
        else: 
            seg_ratio = print("Undefined orientation.")


        # Maintain the original center slice
        start_z = max(0, z_center - roi_size_z // 2)
        end_z = min(max_z, z_center + roi_size_z // 2)

        volume_slices = []
        for z in range(start_z, end_z + 1):
            file_name = f"{z:04d}.tif"
            file_path = os.path.join(tiff_folder, file_name)
            if os.path.exists(file_path):
                img = tiff.imread(file_path)
                roi = img[y_min:y_max, x_min:x_max]
                volume_slices.append(roi)

        if volume_slices:
            roi_volume = np.stack(volume_slices, axis=0)
            volume_file_name = f"{new_name}.tif"
            volume_path = os.path.join(specific_folder, volume_file_name)
            tiff.imwrite(volume_path, roi_volume)

            soma_center_x = x_center - x_min
            soma_center_y = y_center - y_min
            soma_center_z = z_center - start_z  # Adjust the soma center in the new volume

            radius = 10.0  # Default radius
            soma_data = [[1, 1, soma_center_x, soma_center_y, soma_center_z, radius, -1]]
            swc_df = pd.DataFrame(soma_data, columns=['Index', 'Type', 'X', 'Y', 'Z', 'Radius', 'Parent'])
            swc_file_name = f"{new_name}_scaled.swc"
            swc_path = os.path.join(specific_folder, swc_file_name)
            swc_df.to_csv(swc_path, sep=',', index=False)

            output_data.append([volume_file_name, roi_volume.shape[2], roi_volume.shape[1], roi_volume.shape[0],
                                soma_center_x, soma_center_y, soma_center_z,
                                f"{voxel_size[0]}x{voxel_size[1]}x{voxel_size[2]}"])

    output_df = pd.DataFrame(output_data, columns=['File_Name', 'TIF_Dimension_X', 'TIF_Dimension_Y', 'TIF_Dimension_Z', 'New_Soma_X', 'New_Soma_Y', 'New_Soma_Z', 'Voxel_Size'])
    output_df.to_csv(output_csv_path, index=False)
    print("TIF volumes, Soma SWC files, and new Soma XYZ information is generated.")



def get_image_dimensions(tiff_folder):
    tiff_files = [f for f in os.listdir(tiff_folder) if f.lower().endswith('.tif')]
    pattern = re.compile(r"(\d{3})\.tif$", re.IGNORECASE)
    
    for filename in tiff_files:
        match = pattern.search(filename)
        if match:
            three_digit_number = match.group(1)
            new_filename = f"0{three_digit_number}.tif"
            old_path = os.path.join(tiff_folder, filename)
            new_path = os.path.join(tiff_folder, new_filename)
            os.rename(old_path, new_path)

    first_file_path = os.path.join(tiff_folder, "0000.tif") # Make sure that the .tif files do not have long names
    if os.path.exists(first_file_path):
        first_image = tiff.imread(first_file_path)
        max_x, max_y = first_image.shape[1], first_image.shape[0]
        max_z = len(list(Path(tiff_folder).glob('*.tif')))
        return max_x, max_y, max_z
    else:
        raise FileNotFoundError("The first image in the sequence could not be found.")

def main():
    parser = argparse.ArgumentParser(description="Convert SWC to CSV and generate output windows, and new soma locations based on the CSV.")
    parser.add_argument("input_directory", type=str, help="Input directory containing the SWC files")
    parser.add_argument("base_folder", type=str, help="Base folder where the script, TIF/TIFF folder, SWC folder, and soma_location.csv")
    parser.add_argument("tiff_folder_name", type=str, help="Name of the folder containing TIF/TIFF files")
    parser.add_argument("--voxel-size", type=float, nargs=3, help="Voxel sizes for X, Y, Z dimensions", default=[0.305, 0.305, 2.49]) #For LS change it to 0.71, 0.71, 0.71 & default=[0.619, 0.619, 2.49])
    parser.add_argument("--image-orientation", type=str, help="Image orietation, up, down, left, or right")
    args = parser.parse_args()

    # The directory where it would output the soma location csv file
    # channel_suffix = args.tiff_folder_name.split('_')[-1]
    output_csv = os.path.join(args.base_folder, f'original_soma_locations.csv')
    swc_to_csv_multiple(args.input_directory, output_csv)

    # Image processing 
    process_images(args.base_folder, args.tiff_folder_name, args.voxel_size, args.image_orientation)

if __name__ == "__main__":
    main()


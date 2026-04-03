**generate_output_windows.py** use soma location files (.swc files) to generate image crops for individual neurons from a .tiff file.  
## Parameters
- input_directory: Input directory containing the SWC files
- base_folder: Base folder where the script, TIF/TIFF folder, SWC folder, and soma_location.csv
- tiff_folder_name: Name of the folder containing TIF/TIFF files
- voxel-size: Voxel sizes for X, Y, Z dimensions
- image-orientation: The direction in which the apical dendrites of most neurons in the image extend ("up," "down," "left," or "right"). Use "center" if they extend orthogonally. This helps accommodate the asymmetric shape of pyramidal neurons while keeping the crop size as small as possible.  

## Example usage
### Demo dataset
The demo dataset can be accessed [here](https://drive.google.com/drive/folders/1KDgwfF0jWZOF7RksZZwG29AVtmjJSR3K?usp=sharing).  
The expected timeline for processing this demo dataset is <1 hr on a normal desktop computer.  
### Instruction
After downloading the dataset, put the SWC and TIFF folders in a parent folder.  
- INPUT_DIRECTORY: the directory of the SWC folder
- BASE_DIRECTORY: the parent folder directory  
### Example command  
py generate_output_windows.py INPUT_DIRECTORY BASE_DIRECTORY "TIFF" --voxel-size 0.305 0.305 2.49 --image-orientation "up"

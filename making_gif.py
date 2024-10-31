import os
from PIL import Image
import glob

# Configuration variables
#image_folder = "/home/palmero/Work/codes-posdoc/rec_plots_chaotic_systems/results/test/standard/order=9/delta=1.00"
image_folder = "/home/palmero/Work/codes-posdoc/rec_plots_chaotic_systems/results/test/logistic/x0=0.10"
#fGIF = "standard.gif"
fGIF = "logistic.gif"
H = 100  # Height in pixels
W = 100  # Width in pixels

# Function to extract the numeric part from the filename
def extract_index(filename):
    base = os.path.basename(filename)
    index = base.split('.')[0]
    return int(index)

# Create the frames and durations
frames = []
durations = []
image_pattern = os.path.join(image_folder, "*.png")

# Debugging output
print(f"Looking for images in: {image_pattern}")

images = glob.glob(image_pattern)

# Debugging output
print(f"Found {len(images)} images.")

if not images:
    print("No images found. Please check the path and pattern.")
else:
    # Sort images by the numeric part of the filenames
    images.sort(key=extract_index)

    for i in images:
        print(f"Processing image: {i}")  # Debugging output
        newImg = Image.open(i)
        # Rotate the image by 90 degrees
        newImg = newImg.rotate(90, expand=True)
        # Resize the image
        newImg = newImg.resize((W, H))
        frames.append(newImg)

        # Determine the duration based on the filename index
        index = extract_index(i)
        if 0 <= index <= 5000:
            durations.append(2) # 2ms
        elif 5001 <= index <= 10000:
            durations.append(6) # 6ms

    # Save into a GIF file that loops forever
    frames[0].save(fGIF, format='GIF', append_images=frames[1:],
                   save_all=True, duration=durations, loop=0)
    print(f"GIF saved as {fGIF}")
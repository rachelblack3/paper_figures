# make gif of a bunch of bursts 

from PIL import Image
import os

# Path to your folder of PNGs
folder_path = '/data/emfisis_burst/wip/rablack75/rablack75/read_stats/paper_figures/plot_figs_intro/2018_eventB'

# Get all PNG files in the folder, sorted alphabetically
png_files = sorted(
    [f for f in os.listdir(folder_path) if f.lower().endswith(".png")]
)

# Open all the images
frames = [Image.open(png) for png in png_files]

# Save as GIF
frames[0].save(
    "bursts.gif",
    save_all=True,
    append_images=frames[1:],  # add the rest of the frames
    duration=100,              # duration per frame in milliseconds
    loop=0,                    # 0 = loop forever
    transparency=0,            # (optional) handle transparency
    disposal=2                 # (optional) clear frame before next
)

print("GIF saved as bursts.gif")

"""
Creation date: Wed Apr  8 18:46:57 2020
Author: Jimmy Lilly (www.github.com/jlilly364)
Program Objective: 
"""

from PIL import Image
import glob
import numpy as np 
master_path = '/Users/Ryan/Desktop/School/ASTR400B/400B_webster/research project/'

# Create the frames
frames = []
imgs = glob.glob("{}local_group_video/*.png".format(master_path))
# print(np.sort(imgs))
for i in np.sort(imgs):
    new_frame = Image.open(i)
    frames.append(new_frame)
 
# Save into a GIF file that loops forever
frames[0].save('{}plots/local_group_video.gif'.format(master_path), format='GIF',
               append_images=frames[1:],
               save_all=True,
               duration=300, loop=0)
import os
import pyperclip

# Declare which GUI to use for visualisation
os.environ['OVITO_GUI_MODE'] = '1'

from ovito.io import import_file, export_file
from ovito.modifiers import *
from ovito.vis import Viewport
import math
import numpy as np

# Init pipeline from dump file, compute data from given file
pipeline = import_file("Cu_heat_surf/dump/sputter_45_1.dump")
data0 = pipeline.compute(0)  # Compute values at the last frame
n_at_0 = len(data0.particles.positions)
data = pipeline.compute(200)  # Compute values at the last frame
positions = np.array(data.particles.positions)  # Container of all particle positions in the order of the dump file
n_at_end = len(data.particles.positions)

# Mask the positions to return only those with z>20
# print(positions[1850])
above_surf = np.where(positions[:, 2] > 20, positions[:, 2], 0)  # Query the z-column for values above 20, set to 0 else
n_above_surf = np.count_nonzero(above_surf)
atoms_left = n_at_0 - n_at_end + n_above_surf

print(atoms_left)


# View the scene in the pipeline
# pipeline.add_to_scene()
# vp = Viewport()
# vp.type = Viewport.Type.Perspective
# vp.camera_pos = (-18, -150, 150)
# vp.camera_dir = (2, 3, -3)
# vp.fov = math.radians(60.0)
# vp.render_image()

#18.37, 18.17, 18.09, 
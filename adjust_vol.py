# CLAHE local contrast adjustment for multipage tif volumes

from skimage import img_as_float
from skimage import exposure,io 
from skimage import external 
from skimage.color import rgb2gray
import numpy as np

# Load an example image
input_namefile = 'slice.png'
output_namefile = 'adjusted.tif'

img = io.imread(input_namefile) 
 
# Adaptive Equalization
img_adapteq = exposure.equalize_adapthist(img, clip_limit=0.03)
img_adapteq = img_adapteq.astype(np.float16)
img_gray = rgb2gray(img_adapteq) 
external.tifffile.imsave(output_namefile, img_gray,software='Imagej',bigtiff=True) #,plugin="tifffile"

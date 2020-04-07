# CLAHE local contrast adjustment for multipage tif volumes

from skimage import img_as_float
from skimage import exposure,io 
from skimage import external

#Currently saving at 32 bits, we need it at 16 bits

# Load an example image
input_namefile = 'slice.png'
output_namefile = 'adjusted.tif'

img = io.imread(namefile) 
'''
# Contrast stretching
p2, p98 = np.percentile(img, (2, 98))
img_rescale = exposure.rescale_intensity(img, in_range=(p2, p98))

# Equalization
img_eq = exposure.equalize_hist(img)
'''
# Adaptive Equalization
img_adapteq = exposure.equalize_adapthist(img, clip_limit=0.03)
#img_adapteq = img_adapteq.astype(np.uint16)
#io.imsave('result.tif',img_adapteq)
external.tifffile.imsave(output_namefile, img_adapteq,software='Imagej') #,plugin="tifffile"

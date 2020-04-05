# CLAHE local contrast adjustment for multipage tif volumes

from skimage import data,io,exposure,img_as_float
import numpy as np

namefile = 'pvalues_p0.0001_montage_mip_alt2.tif'

im = io.imread(namefile)
#im = data.moon()
print(im.shape)
im = img_as_float(im)
im_s = exposure.rescale_intensity(im, in_range=(-1, 1)) 
#Clipping limit, normalized between 0 and 1 (higher values give more contrast).
img_adapteq = exposure.equalize_adapthist(im_s, clip_limit=0.03)
img_adapteq = img_adapteq.astype(np.uint16)
io.imsave('result.tif',img_adapteq)

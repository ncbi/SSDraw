# helper script to combine images
"""
Example usage: 
python combine_images 1ubq.png \
    1ndd.png \
    2n1v.png

This will create a composite image of these three SSDraw images, ordering them in
the same order listed. The output file is saved as "composite.png"
"""

import os,sys
from PIL import Image

def combine_images(imgs):

    height = min([img.height for img in imgs])
    width = min([img.width for img in imgs])
    pic = Image.new('RGB', (width, height*len(imgs)))

    height_i = 0

    for img in imgs:
        x = img.crop((0,0,width,height))
        pic.paste(x, (0, height_i))
        height_i+=height
    return pic

def main():
    imgs = [Image.open(sys.argv[i]) for i in range(1,len(sys.argv))]
    print("Creating composite.png")
    combine_images(imgs).save("composite.png")

if __name__ == "__main__":
    main()
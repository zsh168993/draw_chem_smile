#!/usr/bin/env python3

import math
import random
from PIL import Image, ImageDraw

def addLine(image, ptA, ptB, width=1, color=(255,0,0)):
    """Draw line from ptA to ptB with arrowhead at ptB"""
    # Get drawing context
    draw = ImageDraw.Draw(image)
    # Draw the line without arrows
    num=int(ptB/10)
    draw.line(((ptA-num, ptB),(ptA+num,ptB)), width=width, fill=color)
    draw.line(((ptA, ptB-num),(ptA, ptB+num)), width=width, fill=color)


    return image
if __name__ == '__main__':
    # Create an empty solid blue image
    w, h = 640, 480
    im = Image.new('RGB', (w,h), (0,0,255))

    im = addLine(im, 320, 240)

    # Save
    im.save('add.png')
from PIL import Image, ImageChops, ImageGrab
import subprocess
import argparse
import os
import sys

parser = argparse.ArgumentParser(
                    prog='trimmer',
                    description='Trim and/or remove background from pictures',
                    epilog='Happy Trimming!')

parser.add_argument('-i', dest='input_file', type=open, nargs='?', help='input file, defult is clipboard content')
parser.add_argument('-b', dest='remove_background', nargs='?', default=False, const=True, help='wheather to remove background based on pixel (0,0)')
parser.add_argument('-t', dest='trim_blankspace', nargs='?', default=False, const=True, help='wheather to remove blank space based on pixel (0,0)')
parser.add_argument('-c', dest='copy', nargs='?', default=False, const=True, help='wheather to copy output to clipboard instead of saving into file')
parser.add_argument('-w', dest='resize_width', nargs='?', help='wheather to resize output to match width')
parser.add_argument('-l', dest='resize_height', nargs='?', help='wheather to resize output to match height')
args = parser.parse_args()

def trim(im):
    bg = Image.new(im.mode, im.size, im.getpixel((0,0)))
    diff = ImageChops.difference(im, bg)
    diff = ImageChops.add(diff, diff, 2.0, -100)
    bbox = diff.getbbox()
    if bbox:
        return im.crop(bbox)
    else: 
        # Failed to find the borders, convert to "RGB"        
        return trim(im.convert('RGB'))
    
def remove_bg(im):
    im = im.convert('RGBA')
    pixels = im.load()
    bg = im.getpixel((0,0))
    for i in range(im.size[0]): # for every pixel:
        for j in range(im.size[1]):
            if pixels[i,j] == bg:
                pixels[i,j] = (255, 255 ,255, 0)
    return im

# print(args.input_file.name)
# im = Image.open(args.input_file.name)

def reWidth(width, im):
    ratio = im.size[0]/im.size[1]
    im = im.resize((width, int(width/ratio)))
    return im

def reHeight(height, im):
    ratio = im.size[0]/im.size[1]
    im = im.resize((int(height/ratio), height))
    return im


# if len(sys.argv) == 2:
if args.input_file:
    im = Image.open(args.input_file.name)
    name = '.'.join(args.input_file.name.split('.')[:-1])+'-trimmed.png'
else:
    im = ImageGrab.grabclipboard()
    name = 'trimmed.png'

if args.resize_width:
    im = reWidth(int(args.resize_width), im)

if args.resize_height:
    im = reHeight(int(args.resize_height), im)


if args.trim_blankspace:
    im = trim(im)

if args.remove_background==True:
    im = remove_bg(im)

im.save(name)

if args.copy:
    subprocess.run('shortcuts run "copy2Clipboard" -i %s'%(name), shell=True)
    # subprocess.run('echo s | pbcopy', shell=True)
    os.remove(name)
    # im.show()

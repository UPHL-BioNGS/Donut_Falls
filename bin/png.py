#!/usr/bin/env python3
from PIL import Image, ImageDraw, ImageFont
import glob
import shutil
import os

def main():
  png_files = glob.glob("*.png")
  png_files.sort()

  if len(png_files) >= 2:
    images_with_titles = []

    for file in png_files:
      analysis = str(file).split('_')[-1].split('.')[0]
      image = Image.open(file)
      draw = ImageDraw.Draw(image)
      draw.text((10, 10), analysis, fill="black", font_size=100)
      images_with_titles.append(image)

    total_width = sum(image.width for image in images_with_titles)
    max_height  = max(image.height for image in images_with_titles)
    combined_image = Image.new("RGB", (total_width, max_height), color="white")
    offset = 0

    for image in images_with_titles:
      combined_image.paste(image, (offset, 0))
      offset += image.width

    combined_image.save("bandage_${prefix}_mqc.png")
    for image in images_with_titles:
        image.close()

  else:
    shutil.copy(png_files[0], "bandage_${prefix}_mqc.png")

if __name__ == "__main__":
    main()
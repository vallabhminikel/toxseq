#!/usr/bin/env python3
"""Convert PNG files in display_items/ to high-resolution, lossless TIFF."""

from pathlib import Path
from PIL import Image

DISPLAY_ITEMS_DIR = Path(__file__).resolve().parent.parent / 'display_items'

# PNGs carry no dpi metadata; this matches the `res` all figures are
# rendered at in src/toxseq_manuscript.R.
DPI = 600

def convert(png_path):
    tiff_path = png_path.with_suffix('.tiff')
    with Image.open(png_path) as img:
        img.save(tiff_path, format='TIFF', compression='tiff_lzw', dpi=(DPI, DPI))
    print(f'{png_path.name} -> {tiff_path.name} ({DPI} dpi)')

def main():
    png_paths = sorted(DISPLAY_ITEMS_DIR.glob('*.png'))
    for png_path in png_paths:
        convert(png_path)

if __name__ == '__main__':
    main()

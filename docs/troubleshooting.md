# Troubleshooting & Common Issues

This section addresses frequently encountered problems and provides solutions.

***

## Issue: Amino Acid Sequence Overlap/Jumbling

### Problem
When the secondary structure diagram is **quite long** (i.e., for large proteins), the residue numbering and/or **amino acid sequence often overlaps or becomes jumbled** together, making it unreadable.

### Solution
You can adjust the size of the residue numbers and sequence labels by changing the default font size. Use the **`--fontsize`** flag followed by a numerical value. The default size is **12**, which may be too large for high-density diagrams. Try reducing it (e.g., `--fontsize 8`) until the text is clearly separated.

* **Example:** `ssdraw single ... --fontsize 8`

***

## Issue: Saving Secondary Structure to a Specific File Format

### Problem
I need to save the generated secondary structure diagram in a **high-quality vector format** (like SVG or EPS) or a specific raster format (like TIFF) instead of the default PNG.

### Solution
You can specify the desired output file type using the **`--output_file_type`** flag. This is particularly useful for generating **publication-quality figures**.

The supported output file types are: **`png`**, **`ps`**, **`eps`**, **`tif`**, and **`svg`**.

* **To save as a vector graphic:** `ssdraw single ... --output_file_type svg`
* **To save as a high-resolution image:** `ssdraw single ... --output_file_type tif`

***
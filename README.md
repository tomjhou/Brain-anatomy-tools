# Brain Anatomy Tools

This repository has several tools for analyzing brain anatomy data.

Although written in Python, you do NOT need to install anything yourself. Just run the .exe or .bat files according to the
instructions below, and the computer sets up what it needs on its own. The .exe has Python and every library inside it. Currently
the section assembler (section 1) and the ABC atlas browser (section 2) have .exe releases; the other tools need to be run from .bat
files. The first time you run a .bat file, it finds Python (installing Python 3.13 if there is none), builds a private `.venv` folder next
to it, and installs the libraries it needs into that; this needs an internet connection, and it leaves any Python you already
have untouched. All the .bat files share the one `.venv`.


## 1. Section alignment program for Windows

The section assembler and aligner (in `FileConverters/`) can read NeuroInfo reconstructions (.xml or .jpx) or raw image
files (Zeiss .czi, or jp2 images), detect sections, sort them in order, and adjust alignment and section outlines to
create serial-section stacks.

This program is available as a Windows .exe program with everything inside it, so nothing has to be installed:

1. Download [SectionAssembler-windows.zip](https://github.com/tomjhou/Brain-anatomy-tools/releases/latest/download/SectionAssembler-windows.zip)
   (always the newest version; older ones and their notes are on the [Releases](https://github.com/tomjhou/Brain-anatomy-tools/releases) page).
2. Unzip it anywhere and run `SectionAssembler\SectionAssembler.exe`.

Windows may warn "Windows protected your PC" the first time, because the program is not code-signed: choose More info, then Run anyway.
To update, download the zip again. How it is built: `packaging/README.md`.


## 2. Allen ABC atlas browser

Interactive browser for the Allen Brain Cell (ABC) atlas MERFISH dataset
(C57BL6J-638850): pick brain sections and regions of interest, compute a UMAP
of the selected cells, and explore it alongside spatial section maps, colored by
taxonomy level or gene expression (including the imputed gene dataset).

### Windows standalone .exe program (easiest way to run)

1. Download [ABCAtlasBrowser-1.0.0-windows.zip](https://github.com/tomjhou/Brain-anatomy-tools/releases/download/atlas-v1.0.0/ABCAtlasBrowser-1.0.0-windows.zip)
   (about 290 MB; the [release page](https://github.com/tomjhou/Brain-anatomy-tools/releases/tag/atlas-v1.0.0) has its notes;
   newer versions will be on the [Releases](https://github.com/tomjhou/Brain-anatomy-tools/releases) page, as `atlas-v...`).
2. Unzip it anywhere and run `ABCAtlasBrowser\ABCAtlasBrowser.exe`. A console window opens and stays open: the program prints its progress
   there and occasionally asks a yes/no question there (for example, to confirm a run with no sections selected), so keep it visible.

Windows may warn "Windows protected your PC" the first time, because the program is not code-signed: choose More info, then Run anyway.
The first time, it asks where to keep the atlas download and suggests a drive with enough free space.
Since the atlas cache can reach about 100 GB, make sure the drive you choose has room.


### Running from the source

If you plan to edit or debug code, you'll need to run from source. The easiest way on Windows is to
double click RUN_ABC_ATLAS_BROWSER.bat. The first time you run it, it creates and installs a
shared `.venv` (see above), which takes about 30 seconds for the section assembler, and longer for the
ABC browser. The `.venv` will be about 300MB with the assembler, and grow to about 1.1 GB with the ABC
browser. This only happens once, and it does not alter existing Python environments.

On other systems, or to use your own Python environment, install Python and run `python ABC_atlas_browser.py` from that
environment (e.g. after `conda activate <env>`): at startup it checks for missing libraries and pip-installs them into it.


### Data locations

- **ABC atlas download:** The packaged Windows ABC browser program keeps the
  ABC atlas cache in locations described under "Windows standalone .exe program" above.
  When run from source (using `.bat` files or `python ABC_atlas_browser.py`)
  for the first time, the app looks for a folder abc_atlas_cache in the PARENT
  folder of the repository. If not found, the app lists your local drives and their
  free space, and asks once where to put it (or where an existing download is),
  suggesting a drive with enough room. The choice is saved in
  `cache_local/abc_atlas_cache_location.json`; delete that file to be asked again.
  Atlas files are downloaded on first use.
- **Local cache:** `cache_local` inside this repository (git-ignored). Section
  images and processed files are regenerated as needed. The packaged Windows program has no repository
  folder of its own to write in (a new version replaces its files), so its `cache_local` is
  `%LOCALAPPDATA%\JhouLab\ABCAtlasBrowser\program\cache_local`, and the Settings file
  `viewer_settings.json` is in the `program` folder beside it.

These default locations are relative to the script itself, so the app can be launched
from any working directory. If you choose your own cache location, this may not apply.


## 3. CCF browser

Scroll through the Allen CCF annotation volume (coronal, sagittal or horizontal slices) and plot cells from CSV files on it:
the coordinates exported by the ABC browser's "Export coordinates" button, or a marker export from NeuroInfo. On Windows, double
click RUN_CCF_BROWSER.bat (it uses the same `.venv`; if that lacks pandas or Pillow they are installed first, and the Allen Institute's
`abc_atlas_access` library, which brings about 300 MB of dependencies, is installed the first time the atlas is needed; neither needs git). The
25 um annotation volume (3.6 MB) comes from the ABC atlas cache, and is downloaded on first use, after asking, if it is not there yet. Making surface meshes of regions needs the `vtk` library, which the
program offers to install the first time (a download of about 100 MB).

## 4. XLEF viewer

A quick viewer for Leica LAS X experiments (`.xlef`): on Windows, double click RUN_XLEF_VIEWER.bat, or drop an `.xlef` file on it
(it uses the same `.venv`, and installs numpy, tifffile and Pillow into it if they are missing). File > Open (Ctrl+O) lists the
experiments opened recently, or lets you browse for one or drag it in; an `.xlef` can also be dragged onto the window. The left panel lists the experiment's
images in LAS X's order, with each image's size at full resolution (and the totals of folders and experiments); click one to see its channels overlaid in their LAS X colors, starting with LAS X's contrast settings. Mouse wheel zooms,
dragging pans, Home (or double click) fits the image to the window, and keys 1-9 switch channels on and off. At the largest zoom, one more turn of the wheel draws a line along the edge of every pixel, with each pixel a square of one color (not smoothed); one turn back removes the lines and smooths the pixels again, at the same zoom. If a sharper view takes more than half a second to read, "Loading full resolution..." shows over the image, with how far it has got, until it is in (meanwhile the background loading of other sections pauses, so as not to share the connection with it). The sharper view is read in bands and drawn as each comes in, so it sharpens from the top down; a plain (uncompressed) TIFF is read several rows at a time, which on a network drive is several times faster than reading it in order. How eagerly sharper pictures are read depends on how fast the files can be read, measured for each drive or network share (a few small reads from random places in a large image file, so that Windows' cache doesn't flatter it, kept up to date by the viewer's own reads), not by the kind of drive: it is measured as a project opens. Where reading would take more than a moment (e.g. a network drive over wifi), a section opens with only its small preview, which is shown enlarged up to 2.5 times before anything sharper is read; then only the part in view is read, from LAS X's sharper (1/4 size) copy, with "Loading a sharper view..." over the image if it takes a while, and full resolution once that copy too would be enlarged more than 2.5 times. Where reading is quick, a section's sharper copy is read whole as it opens, as before. Settings ("Reading sharper pictures") shows the speeds measured, and can also make it always the sharpest, or always save reading. Each channel switched on is shown in a panel of its own, labeled with its name (the default); zooming or panning in any panel moves them all. "Merge channels" (above the channels) overlays them in one picture instead, and "Separate channels" goes back; the choice is remembered. "Flip all sections" (top right, by Settings) shows every image mirrored, with its annotations; "Flip this section" (above the channels) mirrors just the one shown, e.g. a section mounted the other way up (with both ticked, that section is shown unflipped). Only the view is flipped, not the image. "Flip all sections" is remembered on this computer; "Flip this section" is saved by Save into the section's .xlif (in an attachment of the viewer's own, JhouXlefViewer, which LAS X keeps when it saves the file again), so it goes with the files to any computer, and until then the section is marked with the red *. (Sections flipped with an earlier version, which kept this on the computer, are picked up when each is first shown, and written by the next Save.) LAS X keeping such a field is what it does today (tried with LAS X 6.2), not a promise for future versions, so Settings has "Save it in a custom field of LAS X's files" (on by default): off, "Flip this section" is kept on this computer only and LAS X's files are left alone (a flip already saved in a file is still shown). Under the image, the pointer's position is shown, and when zoomed in to full resolution each channel's value there, in the channel's color. "Show annotations" (at the left of the toolbar above the image, on by default) draws the measuring lines (with their lengths), lines, arrows, rectangles,
circles/ellipses and text added to an image in LAS X. The toolbar above the image edits them: Select (drag to move, corners to resize,
the round handle to rotate, double-click a text to change its words or color), drawing tools (a tool stays chosen for one annotation after another; Select or Escape to go back. With any tool, moving the mouse over a shape or line shows its handles (hollow squares; an ellipse's are the corners of its box), the pointer turns to a hand over a handle, and dragging it resizes, rotates or moves an end instead of drawing; over an outline (a rectangle's edge, an ellipse's circumference, a line, a text) it turns to move arrows, and dragging moves the whole annotation), color (shapes and lines start red, text white, each time the viewer starts), Delete and Undo. Changes are written into the
image's LAS X file (`.xlif`) only when you click Save; only the annotations in it are changed, and a copy of the file as it was is kept
next to it (`.before_xlef_viewer.bak`). Close the experiment in LAS X before saving. "⚙ Settings" (top right) switches between
the dark (default) and light color themes, and sets how pixels look when zoomed in past them (Mitchell-Netravali smoothing by
default; also bicubic, Gaussian, bilinear, or hard-edged true pixels). It also sets the memory for previews (50% of the computer's memory by default, at most 80%): the pictures of recently viewed images, and of the ones loaded in the background so that they open at once, and, in what those leave free, the full-resolution views zoomed into before (going back to one needs no reading; they make room for the pictures, never the other way round); when it fills up, the status bar says so. To rename an image or folder, click its name in the list once it is selected (as in
Windows Explorer); Save writes the new name into its `.xlif`/`.xlcf`. With "Also rename the files and folders on disk" (Settings,
off by default) the image's `.xlif` and channel TIFFs, or the folder, are renamed too, and the references to them updated. Each channel has
a swatch (click to switch it on/off) and a black/white slider. Right-click a channel's swatch or slider (or double-click the slider) to choose its color map, as in
LAS X: its list of maps (Gray, Red, ..., Glow (O&U), Spectrum, Rainbow), or any color from the picker. The choice shows on the image at
once; apply it to This section or to All sections (the same channel in every image of the experiment), or Cancel. Undo (Ctrl+Z) takes
it back, and Save writes it into the images' `.xlif` files, where LAS X reads it.
The black/white sliders work the same way: a section keeps the contrast set for it (sliders, Auto contrast) when you come back to it, marked with a red * in the list (hover over it to see what is unsaved; hovering over a section's name shows its name, size and resolution, the number of layers of a z-stack, and its unsaved changes too) until Save writes it into its `.xlif` (LAS X then shows the same contrast); Undo takes back a whole drag at once, and "LAS X contrast" goes back to the contrast in the file.

To change the order of the images, drag their names in the list (or their thumbnails): onto the middle of a folder's name to put
them in it, or into another open experiment, which moves the image's files (its `.xlif`, channel TIFFs and reduced copies) into that
experiment's folder when you Save. Folders can be dragged within their experiment. Right-click an image > Delete section takes it out
of the list; Save removes it from LAS X's files and moves its files to the Windows Recycle Bin (a network drive has no Recycle Bin, so
there they go into a "Deleted by XLEF viewer" folder in the experiment's folder instead). Both can be undone (Ctrl+Z) until Save, and
Save asks first whenever it moves or deletes files.

The viewer only reads what it needs: listing the images reads only the small `.xlif` metadata files, and an image's TIFFs (one per channel, or one OME-TIFF holding all channels as its pages) are
opened only when you select it (and then only the rows the current view needs). When the experiment has LAS X's reduced copies
(a `Pyramid` folder), the zoomed-out picture comes from those, and full resolution is read only as you zoom in. Several experiments
can be open at once: each one opened is added to the list as a branch of its own (right-click its name to close it); hovering over its name shows its path and any unsaved changes, and the slide
overview shows the selected image's experiment. The List / Thumbnails buttons switch the left panel
to a grid of thumbnails; these come from LAS X's small reduced copies (the `Pyramid` folder) when an experiment has them, and are
kept in `%LOCALAPPDATA%\JhouLab\xlef_viewer_thumbs`. Images whose files are online-only in OneDrive are marked with a cloud. Reading
any part of an online-only file makes OneDrive download all of it, so the viewer downloads online-only files of up to 50 MB (such as
the reduced copies) as needed, but asks before any bigger one: when opening an image without reduced copies, or when zooming in far
enough to need full resolution ("Download for this section", "Download for all sections", or Cancel, which keeps the view at
the reduced copies' sharpness).

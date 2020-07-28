# rna_draw

A minimial python package for drawing nucleic acid secondary structures. 

## Dependencies



## Install

Currently rna_draw requires cairosvg (https://cairosvg.org/) which requires cairo graphics (https://www.cairographics.org/download/) to be installed

```shell
#Mac 
brew install cairo

#ubuntu 
sudo apt-get install libcairo2-dev

#fedora 
sudo yum install cairo-devel

```

To install rna_draw 

```shell
python -m pip install git+https://github.com/YesselmanLab/rna_draw.git
```

## TODO

(1) Large structures still generate collisions. 
(2) draw tertiary contacts 

## Usage

There are two ways to call rna_draw. First is the command line exe that was installed with the package There are two ways to call rna_draw. First is the command line exe that was installed with the package 

```shell
usage: draw.py [-h] -ss SS [-seq SEQ] [-out OUT] [-color_str COLOR_STR]
               [-render_type RENDER_TYPE] [-default_color DEFAULT_COLOR]
               [-data_str DATA_STR] [-data_file DATA_FILE]
               [-data_palette DATA_PALETTE] [-data_vmin DATA_VMIN]
               [-data_vmax DATA_VMAX]
               [-data_ignore_restype DATA_IGNORE_RESTYPE]

optional arguments:
  -h, --help            show this help message and exit
  -ss SS                secondary structure in dot bracket notation
  -seq SEQ              rna sequence
  -out OUT              output png file
  -color_str COLOR_STR  description of coloring, see docs for options
  -render_type RENDER_TYPE
                        scheme to color by: res_type,paired,motif,none
  -default_color DEFAULT_COLOR
                        the color used when no other color is supplied
  -data_str DATA_STR    data values by res seperated by ;
  -data_file DATA_FILE  path to data to color by
  -data_palette DATA_PALETTE
                        matplotlib color palette
  -data_vmin DATA_VMIN  min data value everything lower than this will be set
                        to this value
  -data_vmax DATA_VMAX  max data value everything greater than this will be
                        set to this value
  -data_ignore_restype DATA_IGNORE_RESTYPE
                        data values will be ignored for these restypes, e.g. G
                        and U for DMS
```

Second import the package and call rna_draw function. All the same arguments can be used as at commandline.

```python
import rna_draw as rd 
rd.rna_draw(ss=SS, seq=SEQ, out=OUT, color_str=COLOR_STR, render_type=RENDER_TYPE, default_color=DEFAULT_COLOR, data_str=DATA_STR, data=DATA_FILE, data_palette=DATA_PALETTE, data_vmin=DATA_VMIN, data_vmax=DATA_VMAX, data_ignore_restype=DATA_IGNORE_RESTYPE)
```

## How to: 


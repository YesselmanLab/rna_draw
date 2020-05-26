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



## TODO

(1) Large structures still generate collisions. 
(2) draw tertiary contacts 

## Usage

There are two ways to call rna_draw. First is the command line exe that was installed with the package There are two ways to call rna_draw. First is the command line exe that was installed with the package 

```shell
usage: draw_rna.py [-h] -ss SS [-seq SEQ] [-out OUT] [-color_str COLOR_STR]
                   [-data_file DATA_FILE] [-data DATA]
                   [-render_type RENDER_TYPE]

optional arguments:
  -h, --help            show this help message and exit
  -ss SS                secondary structure in dot bracket notation
  -seq SEQ              rna sequence
  -out OUT              output png file
  -color_str COLOR_STR  description of coloring, see docs for options
  -data_file DATA_FILE  path to data to color by
  -data DATA            data values by res seperated by ;
  -render_type RENDER_TYPE
                        scheme to color by options: res_type,paired,motif,none
```

Second import the package and call rna_draw function. All the same arguments can be used as at commandline.

```python
import rna_draw as rd 
rd.rna_draw(ss=SS, seq=SEQ, out=OUT, color_str=COLOR_STR)
```

## How to: 


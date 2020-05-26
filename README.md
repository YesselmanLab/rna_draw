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
(3) changing draw parameters (such as circle borders, size, basepair color, etc)

## Usage

There are two ways to call rna_draw. First is the command line exe that was installed with the package 

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

- [basic usage](#basic-usage)
	- [just secondary structure](#just-secondary-structure)
	- [supplying a sequence](#supplying-a-sequence)
	- [adding render type](#adding-render-type)
	- [second render type example](#second-render-type-example)
	- [supplying custom colors](#supplying-custom-colors)
- [coloring examples](#coloring-examples)
	- [using single color codes](#using-single-color-codes)
	- [changing default color](#changing-default-color)
	- [using xkcd colors](#using-xkcd-colors)
	- [overwriting render type colors](#overwriting-render-type-colors)
- [using experimental data](#using-experimental-data)
	- [basic coloring with data](#basic-coloring-with-data)
## basic usage
### just secondary structure
all that is needed is to supply a secondary structure
```shell
rna_draw -ss ".(((.....)))." 
```
![](resources/imgs/test_0.png)
### supplying a sequence
sequence identity can be supplied through -seq option. Any letter can be supplied.
```shell
rna_draw -ss ".(((.....)))." -seq AUGANNNNNUCAA 
```
![](resources/imgs/test_1.png)
### adding render type
can supply a render type to color with a specific style. res_type, colors by residue color.
```shell
rna_draw -ss ".(((.....)))." -seq AUGANNNNNUCAA -render_type res_type 
```
![](resources/imgs/test_2.png)
### second render type example
there are different types of render types, here is showing "paired", which colors by whether a residue is in a basepair.
```shell
rna_draw -ss ".(((.....)))." -seq AUGANNNNNUCAA -render_type paired 
```
![](resources/imgs/test_3.png)
### supplying custom colors
The user can supply any set of colors as long as they are not overlapping. The format is resnum_min-resnum_max:color. You also do a single resnum in the format resnum:color. In this case we are using single letter codes. The available codes are red (r), greeen (g), blue (b), yellow (y), cyan (c), magenta (m), white (w), gray (e), orange(o)
```shell
rna_draw -ss ".(((.....)))." -seq AUGANNNNNUCAA -color_str 1-5:r;6-13:b 
```
![](resources/imgs/test_4.png)
## coloring examples
### using single color codes
using the single color codes its one can specify a color string with the same length as the secondary structure. The available codes are red (r), greeen (g), blue (b), yellow (y), cyan (c), magenta (m), white (w), gray (e), orange(o).
```shell
rna_draw -ss ".(((.....)))." -color_str rrrrbbbbbgggg 
```
![](resources/imgs/test_5.png)
### changing default color
if the color string does not cover all residues the default color will be used for the rest. You can change this color.
```shell
rna_draw -ss ".(((.....)))." -color_str 1-5:r -default_color y 
```
![](resources/imgs/test_6.png)
### using xkcd colors
colors can be a mixture of single character codes and xkcd colors (https://xkcd.com/color/rgb/)
```shell
rna_draw -ss ".(((.....)))." -color_str purple;green;blue;pink;brown;red;fuchsia;hunter green;mustard yellow;eggplant;off white;r;r 
```
![](resources/imgs/test_7.png)
### overwriting render type colors
color_str always overrides render type colors, to allow for more flexibility in coloring.d xkcd colors (https://xkcd.com/color/rgb/). Here we changed the color of residue 5 to 9 which would of been gray but instead is aqua blue.
```shell
rna_draw -ss ".(((.....)))." -seq AUGANNNNNUCAA -color_str 5-9:aqua blue -render_type res_type 
```
![](resources/imgs/test_8.png)
## using experimental data
### basic coloring with data
data can be supplied directly using ; to seperate each data point
```shell
rna_draw -ss ".(((.....)))." -seq AUGAAAAAAUCAA -data 0;1;2;3;4;5;6;7;8;9;10;11;12;13 
```
![](resources/imgs/test_9.png)

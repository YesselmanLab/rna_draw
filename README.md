# rna_draw

This script generates secondary structure diagrams for nucleic acids.

## Dependencies

The default behavior is to output the secondary structure visualization to both svg and png. 

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

## How to: 

- [basic usage](#basic-usage)
	- [just secondary structure](#just-secondary-structure)
	- [supplying a sequence](#supplying-a-sequence)
### basic usage
#### just secondary structure
all that is needed is to supply a secondary structure
```shell
rna_draw -ss "((((.....))))" 
```
![](resources/imgs/test_0.png)
#### supplying a sequence
sequence identity can be supplied through -seq option. Any letter can be supplied.
```shell
rna_draw -ss "((((.....))))" -seq GUGNNNNNNUCAC 
```
![](resources/imgs/test_1.png)

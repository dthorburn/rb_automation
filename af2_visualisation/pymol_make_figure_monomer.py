#!/usr/bin/env python3

## Launching pymol quiet and with no GUI
import __main__
__main__.pymol_argv = [ 'pymol', '-qc'] 

import sys, time, os, re
import pymol

## Finish launching pymol before any commands are launched
pymol.finish_launching()

from pymol import cmd

## Read user input and creating a simple name
spath = os.path.abspath(sys.argv[1])
fname = spath.split('/')[-1]
sname = re.sub("_unrelaxed.*$", "", fname)
pixels = 900

## Load cofolded structures
cmd.load(spath, sname)
## Colouring overall structures
cmd.color("blue", "chain A and b > 90")
cmd.color("cyan", "chain A and b < 90 and b > 70")
cmd.color("yellow", "chain A and b < 70 and b > 50")
cmd.color("orange", "chain A and b < 50")
#cmd.color("blue", f"({selection}) and b > 90")
#cmd.color("cyan", f"({selection}) and b < 90 and b > 70")
#cmd.color("yellow", f"({selection}) and b < 70 and b > 50")
#cmd.color("orange", f"({selection}) and b < 50")
## Colouring interaction
## Colouring interaction
##cmd.select("nearA")
## remove if you want AF2 confidence colours
cmd.color("oxygen", "chain A")

## Background options
cmd.set('opaque_background')
cmd.bg_color("grey10")

cmd.set('surface_quality',1)
cmd.set('ray_shadow',0)
cmd.ray(pixels, pixels)
cmd.png(sname+"_1", dpi=300)
cmd.turn('y',90)
cmd.ray(pixels, pixels)
cmd.png(sname+"_2", dpi=300)
cmd.turn('y',90)
cmd.ray(pixels, pixels)
cmd.png(sname+"_3", dpi=300)
cmd.turn('y',90)
cmd.ray(pixels, pixels)
cmd.png(sname+"_4", dpi=300)
cmd.turn('y',90)
cmd.set('surface_quality',0)

cmd.quit()

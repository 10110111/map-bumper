#!/usr/bin/env python3

import numpy
import meshzoo

points, cells = meshzoo.icosa_sphere(32)
bytes = points.tofile("/tmp/vertices-coarse.bin")
print("Saved", len(points), "vertices of the coarse mesh")
bytes = cells.tofile("/tmp/indices-coarse.bin")
print("Saved", len(cells), "faces of the coarse mesh")

points, cells = meshzoo.icosa_sphere(388)
bytes = points.tofile("/tmp/vertices-fine.bin")
print("Saved", len(points), "vertices of the fine mesh")
bytes = cells.tofile("/tmp/indices-fine.bin")
print("Saved", len(cells), "faces of the fine mesh")

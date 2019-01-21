import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection

if __name__ == '__main__':
    data = np.genfromtxt("polygon.txt", delimiter=',')
    triangle_indices = np.genfromtxt("triangle_indices.txt", delimiter=',')
    triangle_indices = triangle_indices.astype(int)
    fig, ax = plt.subplots()
    patches = []
    # Plot polygon 
    convex_polygon = Polygon(data, closed = True)
    patches.append(convex_polygon)
    
    # Plot triangle vertices
    for i in range(np.size(triangle_indices,0)):
        triangle = Polygon(data[triangle_indices[i,:]], closed = True) 
        patches.append(triangle)
    polygons = PatchCollection(patches, alpha = 0.4)
    polygons.set_facecolor('red')
    polygons.set_edgecolor('green')
    polygons.set_linewidth(2)
    ax.add_collection(polygons)
    ax.autoscale_view()
    plt.show()

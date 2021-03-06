{
 "cells": [
  {
   "cell_type": "code",
   "execution_count": 1,
   "id": "exclusive-emphasis",
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/html": [
       "<style>.container { width:100% !important; }</style>"
      ],
      "text/plain": [
       "<IPython.core.display.HTML object>"
      ]
     },
     "metadata": {},
     "output_type": "display_data"
    }
   ],
   "source": [
    "from IPython.display import HTML\n",
    "from bgc_md2.resolve.mvars import CompartmentalMatrix, StateVariableTuple, VegetationCarbonInputPartitioningTuple,VegetationCarbonInputTuple\n",
    "from bgc_md2.resolve.MVarSet import MVarSet\n",
    "display(HTML(\"<style>.container { width:100% !important; }</style>\"))\n",
    "mvs_mm =  MVarSet.from_model_name('TECOmm')"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 2,
   "id": "competitive-finder",
   "metadata": {},
   "outputs": [],
   "source": [
    "in_fluxes, internal_fluxes, out_fluxes = mvs_mm.get_InFluxesBySymbol(),mvs_mm.get_InternalFluxesBySymbol(),mvs_mm.get_OutFluxesBySymbol()\n",
    "\n",
    "in_flux_targets, out_flux_sources = [[str(k) for k in d.keys()] for d in (in_fluxes, out_fluxes)] \n",
    "\n",
    "internal_connections = [(str(s),str(t)) for s,t in internal_fluxes.keys()]                                                                "
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 3,
   "id": "informational-width",
   "metadata": {},
   "outputs": [],
   "source": [
    "import CompartmentalSystems.helpers_reservoir as hr\n",
    "G, GVI, GINT, GVO = hr.nxgraphs(mvs_mm.get_StateVariableTuple(),in_fluxes,internal_fluxes,out_fluxes)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 4,
   "id": "damaged-insulation",
   "metadata": {},
   "outputs": [],
   "source": [
    "\n",
    "# We add some properties to the graph that will be used for plotting\n",
    "VC, RC = \"black\", \"red\"\n",
    "n_color = {}\n",
    "n_size={}\n",
    "for n,y in G.nodes(data=True):\n",
    "    n_color[n] = VC if 'virtual' in y.keys() else RC\n",
    "    n_size[n] = 10 if 'virtual' in y.keys() else 50\n",
    "\n",
    "import networkx as nx    \n",
    "nx.set_node_attributes(G,n_color,'node_color')\n",
    "nx.set_node_attributes(G,n_size,'node_size')"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 5,
   "id": "fatal-monitoring",
   "metadata": {},
   "outputs": [],
   "source": [
    "\n",
    "from bokeh.io import output_file, show\n",
    "from bokeh.models import (BoxSelectTool, Circle, Text, EdgesAndLinkedNodes, HoverTool,\n",
    "                          MultiLine, ArrowHead, NodesAndLinkedEdges, Plot, Range1d, TapTool,)\n",
    "from bokeh.palettes import Spectral4\n",
    "from bokeh.plotting import from_networkx\n",
    "\n",
    "#G = nx.karate_club_graph()\n",
    "\n",
    "\n",
    "plot = Plot(plot_width=1000, plot_height=1000,\n",
    "            x_range=Range1d(-1.1,1.1), y_range=Range1d(-1.1,1.1))\n",
    "plot.title.text = \"Graph Interaction Demonstration\"\n",
    "\n",
    "plot.add_tools(HoverTool(tooltips=None), TapTool(), BoxSelectTool())\n",
    "\n",
    "graph_renderer = from_networkx(G, nx.spring_layout, scale=1, center=(0,0))\n",
    "\n",
    "graph_renderer.node_renderer.glyph = Circle(size='node_size', fill_color='node_color')# Spectral4[0])\n",
    "#graph_renderer.node_renderer.glyph = Text()\n",
    "graph_renderer.node_renderer.selection_glyph = Circle(size=15, fill_color=Spectral4[2])\n",
    "graph_renderer.node_renderer.hover_glyph = Circle(size=15, fill_color=Spectral4[1])\n",
    "\n",
    "graph_renderer.edge_renderer.glyph = MultiLine(line_color=\"#CCCCCC\", line_alpha=0.8, line_width=5)\n",
    "graph_renderer.edge_renderer.selection_glyph = MultiLine(line_color=Spectral4[2], line_width=5)\n",
    "graph_renderer.edge_renderer.hover_glyph = MultiLine(line_color=Spectral4[1], line_width=5)\n",
    "\n",
    "graph_renderer.selection_policy = NodesAndLinkedEdges()\n",
    "graph_renderer.inspection_policy = EdgesAndLinkedNodes()\n",
    "\n",
    "\n",
    "plot.renderers.append(graph_renderer)\n",
    "\n",
    "output_file(\"networkx_graph.html\")\n",
    "show(plot)"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "encouraging-damage",
   "metadata": {},
   "source": [
    "To do:\n",
    "* For directed graphs write a new glyph wiht arrowhead to be used as   `graph_renderer.edge_renderer.glyph` instead of `MultiLine`\n",
    "* write a similar function to `from_networkx` that accepts a nicer layout e.g  `igraph.layout('sugiyama')`\n",
    "* Use the `selected` or `inspected` nodes or edges to print some useful information or do something\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "opponent-regular",
   "metadata": {},
   "outputs": [],
   "source": []
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "Python 3",
   "language": "python",
   "name": "python3"
  },
  "language_info": {
   "codemirror_mode": {
    "name": "ipython",
    "version": 3
   },
   "file_extension": ".py",
   "mimetype": "text/x-python",
   "name": "python",
   "nbconvert_exporter": "python",
   "pygments_lexer": "ipython3",
   "version": "3.8.8"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 5
}

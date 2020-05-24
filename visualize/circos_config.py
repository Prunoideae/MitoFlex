"""
circos.default.py
=========

Copyright (c) 2019-2020 Li Junyu <2018301050@szu.edu.cn>.

This file is part of MitoFlex.

MitoFlex is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MitoFlex is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MitoFlex.  If not, see <http://www.gnu.org/licenses/>.

"""

import os
import sys
try:
    sys.path.insert(0, os.path.abspath(os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "..")))
    from utility.bio import circos
    from configurations import visualize as v_conf
except ImportError as err:
    sys.exit(f"Unable to import helper module {err.name}, is the installation of MitoFlex valid?")

# The color that will be used in marked genes.
# 0 for protein, 1 for tRNAs, 2 for rRNAs
fill_colors = [v_conf.color_cds, v_conf.color_trna, v_conf.color_rrna]

# The main configuration entry of circos config file, can be
# adjusted as will, but some part of the graph needs you to modify
# visualize.py if some calculation is involved.
circos_conf = circos.Circos()

with circos_conf.image as image:
    # Base directory
    image.dir = None
    image.file = "Circos.png"
    image.png = "yes"
    image.svg = "yes"

    image.radius = "1500p"
    image.angle_offset = -90

    image.auto_alpha_colors = "yes"
    image.auto_alpha_steps = 5
    image.background = "white"

with circos_conf.ideogram as ideogram:

    with ideogram.spacing as spacing:
        spacing.default = "0.01r"
        # Circular will be 0
        spacing._break = "0.5r"

    ideogram.radius = "0.82r"
    ideogram.thickness = '20p'
    ideogram.fill = 'yes'
    ideogram.fill_color = 'grey'
    ideogram.stroke_thickness = 3
    ideogram.stroke_color = "black"

    ideogram.show_label = 'yes'
    ideogram.label_font = 'bolditalic'
    ideogram.label_radius = 'dims(ideogram,radius_outer) - 0.1r'
    ideogram.label_size = 28
    ideogram.label_parallel = "yes"
    ideogram.label_case = 'lower'
    ideogram.show_bands = 'yes'
    ideogram.fill_bands = 'yes'
    ideogram.band_stroke_thickness = 2
    ideogram.band_stroke_color = "white"
    ideogram.band_transparency = 0

circos_conf.show_ticks = "yes"
circos_conf.show_tick_labels = 'yes'

with circos_conf.ticks as ticks:

    ticks.radius = 'dims(ideogram,radius_outer)'
    ticks.orientation = 'out'
    ticks.label_multiplier = 1e-3
    ticks.color = 'black'
    ticks.thickness = '2p'
    ticks.font = 'bold'

    with ticks['tick'] as tick:
        tick.spacing = "1u"
        tick.show_label = 'yes'
        tick.label_size = '25p'
        tick.size = '25p'
        tick.format = '%d'
        tick.label_offset = '2p'

    with ticks['tick'] as tick:
        tick.spacing = '5u'
        tick.show_label = 'yes'
        tick.label_size = '30p'
        tick.size = '30p'
        tick.format = '%d'
        tick.suffix = '" kb"'
        tick.label_offset = '2p'

    with ticks['tick'] as tick:
        tick.spacing = '10u'
        tick.show_label = 'yes'
        tick.label_size = '30p'
        tick.size = '30p'
        tick.format = '%d'
        tick.suffix = '" kb"'
        tick.label_offset = '2p'

circos_conf.karyotype = None
circos_conf.chromosomes_units = 1000
circos_conf.chromosomes_display_default = 'yes'

with circos_conf.plots as plots:
    with plots['plot'] as plot:
        plot.type = 'text'
        plot.color = 'black'
        plot.label_font = 'default'
        plot.label_size = '28p'
        # Gene position file
        plot.file = None
        plot.r1 = '1r+300p'
        plot.r0 = '1r+10p'
        plot.show_links = 'yes'
        plot.link_dims = '0p,0p,70p,0p,10p'
        plot.link_thickness = '2p'
        plot.link_color = 'red'

        plot.label_snuggle = 'yes'
        plot.max_snuggle_distance = '1r'
        plot.snuggle_tolerance = '0.25r'
        plot.sunggle_sampling = 2

    with plots['plot'] as plot:
        plot.type = 'text'
        plot.color = 'black'
        plot.label_font = 'bold'
        plot.label_size = '40p'
        # Plus file
        plot.file = None
        plot.show_links = 'no'

    with plots['plot'] as plot:
        plot.type = 'histogram'
        # GC content file
        plot.file = None
        plot.r1 = '0.615r'
        plot.r0 = '0.45r'
        plot.max = 1
        plot.min = 0

        plot.stroke_type = 'line'
        plot.thickness = 2
        plot.color = '128,177,211'
        plot.extend_bin = 'no'
        plot.fill_color = '128,177,211'

        with plot.axes as axes:

            with axes['axis'] as axis:
                axis.spacing = '0.05r'
                axis.color = 'lgrey'
                axis.thickness = 1

            with axes['axis'] as axis:
                axis.position = '0.5r'
                axis.color = 'dred'
                axis.thickness = 2

    with plots['plot'] as plot:
        plot.type = 'line'
        plot.thickness = 2
        plot.max_gap = '1u'
        plot.skip_run = 'yes'
        # Sequencing depth out file
        plot.file = None
        plot.color = 'dgreen'
        plot.min = 0
        # Maximum depth, calculated on drawing
        plot.max = -1
        plot.r0 = '0.618r'
        plot.r1 = '0.768r'
        plot.fill_color = '190,186,218'

        with plot.axes as axes:
            with axes['axis'] as axis:
                axis.color = 'lgrey_a2'
                axis.thickness = 1
                axis.spacing = '0.06r'

        with plot.rules as rules:

            with rules['rule'] as rule:
                # Over 90 percent of the maximum depth
                rule.condition = 'var(value) > {}'
                rule.color = '20,227,117'
                rule.fill_color = '20,227,117'

            with rules['rule'] as rule:
                # Below the 10 percent of maximum depth
                rule.condition = 'var(value) < {}'
                rule.color = 'dred'
                rule.fill_color = 'dred_a1'

with circos_conf.highlights as highlights:

    with highlights['highlight'] as highlight:
        highlight.file = None

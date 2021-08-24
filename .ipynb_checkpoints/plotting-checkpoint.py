import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as mpatches

"""
This package contains:
- frac, stack, rack plots - fraction box plots that can be plotted horizontally or vertically and then stacked and racked together.
- useful color definitions.
"""

"""
Start frac, stack, racks plots
"""

def box_fraction_plot(axis, x_bounds, y_bounds, fraction_list, color_list, horizontal=True):
    # will plot coloured vertical stripes with width determined by fraction.
    # x_bounds = [xmin, xmax]
    # y_bounds = [ymin, ymax]
    # fractions = [a,b,c,d] - should sum to one but will be normalized anyway
    # colors = [color for a, ... b, c, d]
    # horizontal, if True sort fractions horizontally, if false sort fractions vertically.
    
    # Normalize input fractions to ensure they sum to 1.
    fraction_total = sum(fraction_list)
    fractions = [F/fraction_total for F in fraction_list]
    
    # Sum of fractions
    fractions_sum = [ sum(fractions[0:(IDX+1)]) for IDX in range(len(fractions)) ]
    
    # range for x and y
    x_range = x_bounds[1] - x_bounds[0]
    y_range = y_bounds[1] - y_bounds[0]
    
    # list to store patches
    patches = []
    
    for IDX in range(len(fraction_list)):
        
        frac_start = fractions_sum[IDX] - fractions[IDX]
        frac_end = fractions_sum[IDX]
        
        if horizontal:
            xy_anchor = (x_bounds[0] + frac_start * x_range, y_bounds[0])
            width = fractions[IDX] * x_range
            height = y_range
        else: # vertical
            xy_anchor = (x_bounds[0], y_bounds[0] + frac_start * y_range)
            width = x_range
            height = fractions[IDX] * y_range
        
        fraction_patch = mpatches.Rectangle(xy_anchor, width, height, facecolor=color_list[IDX], linewidth=0)
        patches.append(fraction_patch)
        
    for p in patches:
        axis.add_patch(p)
    #endfor
#end def

def box_fraction_stack(axis, fraction_stack, label_stack, color_list, x_bounds=[0,1], apply_labels=True):
    # Stacks horizontal plot_box_fraction plots on top of one another, labeling them.
    # fraction stack, a list of lists of fractions to be box_plotted
    # label stack, the labels for the stack
    # color_list the list of colors for the list of fractions for the box plots.
    # apply_label = True, applies labels on the y-axis.
    
    num_plots = len(label_stack)
    
    stack_height = 1.0 / num_plots
    
    stack_bottoms = [IDX * stack_height for IDX in range(num_plots)]
    stack_tops = [(IDX+1) * stack_height for IDX in range(num_plots)]
    stack_centres = [(IDX+0.5) * stack_height for IDX in range(num_plots)]
    
    for IDX in range(num_plots):
        
        y_bounds = [stack_bottoms[IDX], stack_tops[IDX]]
        fraction_list = fraction_stack[IDX]
        
        box_fraction_plot(axis, x_bounds, y_bounds, fraction_list, color_list)
    #endfor
        
    if apply_labels:
        axis.set_yticks(stack_centres)
        axis.set_yticklabels(label_stack)
    #endif
    
#enddef

def box_fraction_rack(axis, fraction_rack, label_rack, color_list, y_bounds=[0,1], apply_labels=True):
    # Stacks vertical plot_box_fraction plots left to right, labeling them.
    # fraction rack, a list of lists of fractions to be box_plotted
    # label rack, the labels for the stack
    # color_list the list of colors for the list of fractions for the box plots.
    # apply_label = True, applies labels on the y-axis.
    
    num_plots = len(label_rack)
    
    rack_width = 1.0 / num_plots
    
    rack_lefts = [IDX * rack_width for IDX in range(num_plots)]
    rack_rights = [(IDX+1) * rack_width for IDX in range(num_plots)]
    rack_centres = [(IDX+0.5) * rack_width for IDX in range(num_plots)]
    
    for IDX in range(num_plots):
        
        x_bounds = [rack_lefts[IDX], rack_rights[IDX]]
        fraction_list = fraction_rack[IDX]
        
        box_fraction_plot(axis, x_bounds, y_bounds, fraction_list, color_list, horizontal=False)
    #endfor
        
    if apply_labels:
        axis.set_xticks(rack_centres)
        axis.set_xticklabels(label_rack, rotation=-45, ha='left')
    #endif
    
#enddef

def box_fraction_stack_rack(axis, fraction_stacks, label_stack, label_rack, color_list, apply_labels=True):
    
    num_racks = len(label_rack)
    
    rack_width = 1.0 / num_racks
    
    rack_lefts = [IDX * rack_width for IDX in range(num_racks)]
    rack_rights = [(IDX+1) * rack_width for IDX in range(num_racks)]
    rack_centres = [(IDX+0.5) * rack_width for IDX in range(num_racks)]
    
    for IDX in range(num_racks):
        
        if IDX == 0 and apply_labels:
            apply_y_labels=True
        else:
            apply_y_labels=False
            
        x_bounds = [rack_lefts[IDX], rack_rights[IDX]]
        
        box_fraction_stack(axis, fraction_stacks[IDX], label_stack, color_list, x_bounds=x_bounds, apply_labels=apply_y_labels)
    #endfor
        
    if apply_labels:
        axis.set_xticks(rack_centres)
        axis.set_xticklabels(label_rack, rotation=-45, ha='left')
    #endif
    
"""
End frac, stack, racks
"""

"""
Useful color definitions
"""

#colors extracted from here: https://twitter.com/NPRougier/status/1323575342204936192/photo/1
# left-hand, simpler set.
ORANGE0 = '#FFF4E6'
ORANGE1 = '#FFE8CC'
ORANGE2 = '#FFD8A8'
ORANGE3 = '#FFC078'
ORANGE4 = '#FFA94D'
ORANGE5 = '#FF922B'
ORANGE6 = '#FD7E14'
ORANGE7 = '#F76707'
ORANGE8 = '#E85900'
ORANGE9 = '#D9480F'
YELLOW0 = '#FFF9DB'
YELLOW1 = '#FFF3BF'
YELLOW2 = '#FFEC99'
YELLOW3 = '#FFE066'
YELLOW4 = '#FFD43B'
YELLOW5 = '#FCC419'
YELLOW6 = '#FAB005'
YELLOW7 = '#F59F00'
YELLOW8 = '#F08000'
YELLOW9 = '#E67700'
LIME0 = '#F4FCE3'
LIME1 = '#E9FAC8'
LIME2 = '#D8F5A2'
LIME3 = '#COEB75'
LIME4 = '#A9E34B'
LIME5 = '#94D82D'
LIME6 = '#82C91E'
LIME7 = '#743816'
LIME8 = '#66A80F'
LIME9 = '#5C940D'
GREEN0 = '#EBFBEE'
GREEN1 = '#D3F9D8'
GREEN2 = '#B2F2BB'
GREEN3 = '#8CE99A'
GREEN4 = '#69DB7C'
GREEN5 = '#51CF66'
GREEN6 = '#40C057'
GREEN7 = '#37B24D'
GREEN8 = '#2F9E44'
GREEN9 = '#2B8A3E'
TEAL0 = '#E6FCF5'
TEAL1 = '#C3FAE8'
TEAL2 = '#96F2D7'
TEAL3 = '#63E6BE'
TEAL4 = '#38D9A9'
TEAL5 = '#20C997'
TEAL6 = '#12B886'
TEAL7 = '#OCA678'
TEAL8 = '#099268'
TEAL9 = '#087F5B'
CYANO = '#E3FAFC'
CYAN1 = '#C5F6FA'
CYAN2 = '#99E9F2'
CYAN3 = '#66D9E8'
CYAN4 = '#3BC9DB'
CYAN5 = '#22B8CF'
CYAN6 = '#15AABF'
CYAN7 = '#1098AD'
CYAN8 = '#0C8599'
CYAN9 = '#0B7285'
BLUE0 = '#E7F5FF'
BLUE1 = '#DOEBFF'
BLUE2 = '#A5D8FF'
BLUE3 = '#74COFC'
BLUE4 = '#4DABF7'
BLUE5 = '#339AFO'
BLUE6 = '#228BE6'
BLUE7 = '#1C7EDO'
BLUE8 = '#1971C2'
BLUE9 = '#1864AB'
INDIGO0 = '#EDF2FF'
INDIGO1 = '#DBE4FF'
INDIGO2 = '#BACOFF'
INDIGO3 = '#91A7EF'
INDIGO4 = '#748FFC'
INDIGO5 = '#5C7CFA'
INDIGO6 = '#4C6EF5'
INDIGO7 = '#4263EB'
INDIGO8 = '#3B5BDB'
INDIGO9 = '#364FC7'
VIOLET0 = '#F3FOFF'
VIOLET1 = '#E5DBFF'
VIOLET2 = '#DOBFFF'
VIOLET3 = '#B197FC'
VIOLET4 = '#9775FA'
VIOLET5 = '#845EF7'
VIOLET6 = '#7950F2'
VIOLET7 = '#7048E8'
VIOLET8 = '#674109'
VIOLET9 = '#5F3DC4'
GRAPE0 = '#F8FOFC'
GRAPE1 = '#F3D9FA'
GRAPE2 = '#EEBEFA'
GRAPE3 = '#E599F7'
GRAPE4 = '#DA77F2'
GRAPE5 = '#CC5DE8'
GRAPE6 = '#BE4BDB'
GRAPE7 = '#AEЗEC9' #fails, misinterpreted as RGBA
GRAPE8 = '#9C36B5'
GRAPE9 = '#862E9C'
PINKO = '#FFF0F6'
PINK1 = '#FFDEEB'
PINK2 = '#FCC2D7'
PINK3 = '#FAA2C1'
PINK4 = '#F783AC'
PINK5 = '#F06595'
PINK6 = '#E64980'
PINK7 = '#D6336C'
PINK8 = '#C2255C'
PINK9 = '#A61E4D'
RED0 = '#FFF5F5'
RED1 = '#FFЕЗЕЗ'
RED2 = '#FFC9C9'
RED3 = '#FFA8A8'
RED4 = '#FF8787'
RED5 = '#FF6B6B'
RED6 = '#FA5252'
RED7 = '#FОЗЕЗЕ' #fails, misinterpreted as RGBA
RED8 = '#E03131'
RED9 = '#C92A2A'
GRAY0 = '#F8F9FA'
GRAY1 = '#F1F3F5'
GRAY2 = '#E9ECEF'
GRAY3 = '#DEE2E6'
GRAY4 = '#CED4DA'
GRAY5 = '#ADB5BD'
GRAY6 = '#868E96'
GRAY7 = '#495057'
GRAY8 = '#343A40'
GRAY9 = '#212529'
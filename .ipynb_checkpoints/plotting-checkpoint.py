import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as mpatches

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
    
    num_plots = len(label_stack)
    
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
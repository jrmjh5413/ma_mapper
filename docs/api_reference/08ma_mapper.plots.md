# ma_mapper.plots
matplotlib wrapper for plotting 
## list of functions
- mapper.plots.plot()
- mapper.plots.plot_heatmap()
- mapper.plots.plot_bar()
- mapper.plots.plot_colorbar()
- mapper.plots.plot_annotation()
- mapper.plots.plot_basecount()
- mapper.plots.plot_logos()

## mapper.plots.plot()

`mapper.plots.plot(data=None, show_alignment=False, alignment=None, alignment_col='nulc_white', heatmap=True, heatmap_mode='overlay', vlim=None, heatmap_color=None, heatmap_title=None, heatmap_title_fs=None, heatmap_yhighlight=None, heatmap_yhighlight_col=['grey'], heatmap_yhighlight_alpha=[0.5], heatmap_xhighlight=None, heatmap_xhighlight_col=['grey'], heatmap_xhighlight_alpha=[0.5], heatmap_ylabel=None, heatmap_ylabel_fs=6, heatmap_xlabel=None, heatmap_xlabel_fs=6, aggregated=False, aggregated_data=None, agg_h=10, agg_ylim=[None], agg_ylabel=None, agg_ylabel_fs=6, agg_ylabel_ypos=None, agg_ylabel_xpos=None, agg_yhighlight=None, agg_yhighlight_col=['grey'], agg_yhighlight_alpha=[0.5], agg_xhighlight=None, agg_xhighlight_col=['grey'], agg_xhighlight_alpha=[0.5], agg_xlabel=None, agg_xlabel_fs=None, agg_major_tick=50, agg_yscale=None, agg_ylabel_right=None, agg_ylabel_right_fs=6, agg_ylabel_right_pos=[1.1, 0.5], agg_plottext=None, agg_plottext_fs=6, agg_plottext_pos=[0.99, 0.90], agg_plot_title=None, agg_plot_title_fs=None, agg_colset=None, logos=False, base_count=False, annotation=False, annotation_data=None, anno_title=None, anno_w=2, anno_ylabel=None, anno_ylabel_fs=6, anno_col=['Blues'], anno_cbar=True, anno_cbar_title=None, anno_cbar_label=None, anno_cbar_even_pos=-0.25, show_alignment_colbar=False, colorbar=False, colorbar_width=2, colorbar_steps=None, image_res=600, figsize=[50, 50], show_plot=True, save_to_file=False, output_filepath=None)`

One of the main APIs of ma_mapper package. This function is the wrapper for matplotlib functions to visualize data matrix extracted from genome-wide dataset using modules in ma_mapper, in a spcific layout, mainly a heatmap of data points with possible annotation and colorbar, accompanied by an aggregated plot under the heatmap.  

### Parameters
mapper.plots.plot() has many parameters inherited from different functions, so they can be categorized by their functions and prefixes:

#### Main Parameters: Heatmap area 
**data : list[numpy.ndarray], numpy.ndarray, default: None**  
The data matrix or the list of data matrix to be plotted.

**show_alignment : bool, default=False**  
Show alignment matrix on the heatmap area.

**alignment : numpy.ndarray, default: None**  
The alignment matrix which is used for data extraction (optional)

**alignment_col : str, default: 'nulc_white'**  
The color palletes for alignment matrix visualization

- 'nulc_white': all bases are white with gaps as grey, used for showing where are the gaps on the heatmap
- 'dna': ACTG - green, yellow, blue, red, give CG yellow and red as they are hotspots for mutation
- 'dna_jaspar' ACTG- green, blue, red, yellow, same color pallete as JASPAR motif database

**heatmap : bool, default: True**  
Show heatmap area
**heatmap_mode : str, default: 'overlay'**  
The layout of heatmap area.

- 'overlay': all matrices are overlaid on top of each other
- 'spread_horizontal': the matrices are plotted horizontally
- 'spread_vertical': the matrices are plotted vertically

**vlim : list[[float, float]], default: None**  
The minimum and maximum value to be shown on the heatmap. There can be more than one set to control more than one matrix.

**heatmap_color : matplotlib.colors.ListedColormap,str,list[matplotlib.colors.ListedColormap,str], default: None**  
The colormap or the list of colormap for heatmap areas, can use colormaps from matplotlib, or custom built-in colormaps such as: custom_cmap.vlag_mpl, custom_cmap.vlag_r_mpl, custom_cmap.RdBu_s_mpl, custom_cmap.RdBu_sr_mpl, and custom_cmap.LimeGreens_mpl which emulate colormap from the package seaborn

**heatmap_title : list[str], default: None**  
The list of the titles of the heatmap area

**heatmap_title_fs : int, default: None**  
The font size of the titles of the heatmap area

**heatmap_yhighlight : list[[int,int]], default: None**  
The list of coordinates of the vertical highlight on the heatmap area

**heatmap_yhighlight_col : list[str], default: ['grey']**  
The list of colors of the vertical highlight on the heatmap area

**heatmap_yhighlight_alpha : list[float], default: [0.5]**  
The list of alpha values of the vertical highlight on the heatmap area 

**heatmap_xhighlight : list[[int,int]], default: None**  
The list of coordinates of the horizontal highlight on the heatmap area

**heatmap_xhighlight_col : list[str], default: ['grey']**  
The list of colors of the horizontal highlight on the heatmap area

**heatmap_xhighlight_alpha : list[float], default: [0.5]**  
The alpha value of the horizontal highlight on the heatmap area

**heatmap_ylabel : list[str], default: None**  
The list of y-axis labels of the heatmap area

**heatmap_ylabel_fs : int, default: 6**  
The font size of the y-axis label of the heatmap area

**heatmap_xlabel : list[str], default: None**  
The list of x-axis labels of the heatmap area

**heatmap_xlabel_fs : int, default: 6**  
The font size of the x-axis label of the heatmap area

#### Additional Parameters: Aggregated dataplot area 
**aggregated : bool, default: False**  
Show aggregated dataplot area

**aggregated_data : list[numpy.ndarray], default: None**  
The aggregated data array or the list of aggregated data array to be plotted.

**agg_h : int, default: 10**  
The height of the aggregated dataplot area

**agg_ylim : list[[float, float]], default: [None]**  
The minimum and maximum value of y-axis of the aggregated dataplot. There can be more than one set to control more than one array.

**agg_ylabel : list[str], default: None**  
The list of y-axis labels of the aggregated dataplot area

**agg_ylabel_fs : int, default: 6**  
The font size of the y-axis label of the aggregated dataplot area

**agg_ylabel_ypos : list[float], default: None**  
The list of y-position of y-axis labels of the aggregated dataplot area

**agg_ylabel_xpos : list[float], default: None**  
The list of x-position of y-axis labels of the aggregated dataplot area

**agg_yhighlight : list[[int,int]], default: None**  
The list of coordinates of the vertical highlight on the aggregated dataplot area  

**agg_yhighlight_col : list[str], default: ['grey']**  
The list of colors of the vertical highlight on the aggregated dataplot area

**agg_yhighlight_alpha : list[float], default: [0.5]**  
The list of alpha values of the vertical highlight on the aggregated dataplot area 

**agg_xhighlight : list[[int,int]], default: None**  
The list of coordinates of the horizontal highlight on the aggregated dataplot area

**agg_xhighlight_col : list[str], default: ['grey']**  
The list of colors of the horizontal highlight on the aggregated dataplot area 

**agg_xhighlight_alpha : list[float], default: [0.5]**  
The alpha value of the horizontal highlight on the aggregated dataplot area

**agg_xlabel : str, default: None**  
The x-axis label of the aggregated dataplot area

**agg_xlabel_fs : int, default: 6**  
The font size of the x-axis label of the aggregated dataplot area

**agg_major_tick : int, default: 50**  
The major tick of the x-axis label of the aggregated dataplot area

**agg_yscale : list[str], default: None**  
The list of the scale of y-axis of the aggregated dataplot area. See matplotlib.scale

**agg_ylabel_right : list[str], default: None**  
The list of the y-lables on the right side of the aggregated dataplot area

**agg_ylabel_right_fs : int, default: 6**  
The font size of the y-lables on the right side of the aggregated dataplot area

**agg_ylabel_right_pos : list[int,int], default: [1.1,0.5]**  
The coordinates of the y-lables on the right side of the aggregated dataplot area

**agg_plottext : list[str], default: None**  
The list of text inside the aggregated dataplot area

**agg_plottext_fs : int, default: 6**  
The font size of the text inside the aggregated dataplot area

**agg_plottext_pos : list[int,int], default: [0.99,0.90]**  
The coordinates of the text of the aggregated dataplot area

**agg_plot_title : list[str], default: None**  
The list of the titles of the aggregated dataplot area

**agg_plot_title_fs : int, default: None**  
The font size of the titles of the aggregated dataplot area

**agg_colset : list[str], default: None**  
The list of the colors for the plot in the aggregated dataplot area

**logos : bool, default: False**  
The additional panel for logos representation of the alignment

**base_count : bool**  
The additional panel for base count of the alignment

#### Additional Parameters: Annotation area
**annotation : bool, default: False**  
Show annotation area

**annotation_data : list[[str,int]], default: None**  
The list of data annotations to be plotted.

**anno_title : list[str], default: None**  
The list of the title of the data annotations

**anno_w : int**  
The width of the annotation area

**anno_ylabel : str, default: None**  
The y-axis label of the annotation area

**anno_ylabel_fs : int, default: 6**  
The font size of the y-axis label of the annotation area 

**anno_col : list[matplotlib.colors.ListedColormap,str], default: ['Blues']**  
The colormap or the list of colormap for data annotations, can use colormaps from matplotlib, or custom built-in colormaps such as: custom_cmap.vlag_mpl, custom_cmap.vlag_r_mpl, custom_cmap.RdBu_s_mpl, custom_cmap.RdBu_sr_mpl, and custom_cmap.LimeGreens_mpl which emulate colormap from the package seaborn

**anno_cbar : bool, default: True**  
Show the colorbar of the data annotations

**anno_cbar_title : list[str], default: None**  
The list of the titles of the annotaion colorbars

**anno_cbar_label : list[[str]], default: None**  
The list of the labels of the annotaion colorbars

**anno_cbar_even_pos : float, default: -0.25**  
Adjust annotaiton colobar titles when the number of the annotations is even

#### Additional Parameters: Colorbar area
**show_alignment_colbar : bool, default: False**  
Show the colorbar of the multiple alignment

**colorbar : bool, default: False**  
Show the colorbars of the heatmap area

**colorbar_width : int, default: 2**  
The width of the colorbars

**colorbar_steps : list[float], default: None**  
The list of the value step in the colorbars

#### Main Parameters: Overall plot
**image_res : int, default: 600**  
The image resolution of the overall plot in DPI

**figsize : list[float, float], default: [50,50]**  
The size of the figure

**show_plot : bool, default: True**  
Show the overall plot

**save_to_file : bool, default: False**  
Save the output as PNG file. The filepath can be specified, using save_to_file = True would save the list to the same directory as the coordinate table or current working directory.

## mapper.plots.plot_heatmap()

`mapper.plots.plot_heatmap(data=None, save_to_file=False, image_res=600, cmap=None, plot_title=None, title_fs=None, xlim=None, ylim=None, xlabel=None, ylabel=None,  xticklabels = None, xticklabels_fs=None, xticklabels_rt=None, figsize=[2,2], opacity=1.0, transparency_mode='constant', show_plot=False, matplot_axes=None, zero_col_centering=False, interpolation = 'none')`

This function is the wrapper for matplotlib functions to visualize data matrix extracted from genome-wide dataset using modules in ma_mapper, as a heatmap.

### Parameters
**data : list[numpy.ndarray], numpy.ndarray: default: None**  
The data matrix to be plotted.

**save_to_file : bool, default: False**  
Save the output as PNG file. The filepath can be specified, using save_to_file = True would save the list to the same directory as the coordinate table or current working directory.

**image_res : int, default: 600**  
The image resolution of the plot in DPI

**cmap : matplotlib.colors.ListedColormap,str,list[matplotlib.colors.ListedColormap,str], default: None**  
The colormap for heatmap plot, can use colormaps from matplotlib, or custom built-in colormaps such as: custom_cmap.vlag_mpl, custom_cmap.vlag_r_mpl, custom_cmap.RdBu_s_mpl, custom_cmap.RdBu_sr_mpl, and custom_cmap.LimeGreens_mpl which emulate colormap from the package seaborn.

**plot_title : str, default: None**  
The title of the heatmap plot

**title_fs : inst, default: None**  
The font size of the title of the heatmap plot

**xlim : list[float, float], default: None**  
The maximum and minimum values of the x-axis of the heatmap plot

**ylim : list[float, float], default: None**  
The maximum and minimum values of the y-axis of the heatmap plot

**xlabel : str, default: None**  
The x-axis label of the heatmap plot

**ylabel : str, default: None**  
The y-axis label of the heatmap plot

**xticklabels : list[str], default: None**  
The lables of major ticks on x-axis.

**xticklabels_fs : float, default: None**  
The font size of the major tick labels on x-axis

**xticklabels_rt : float, default: None**  
The rotation of the major tick labels on x-axis

**figsize : list[float, float], default: [2,2]**  
The size of the figure

**opacity : float, default: 1.0**  
The opacity value of the heatmap plot

**transparency_mode : str, default: 'constant'**  
The pattern of the transparency values
- 'constant' use constant opacity
- 'gradient' the opacity increase as the value of the data increased from 0 to the set opacity above at maximum value

**show_plot : bool, default: False**  
Render the heatmap plot as a standalone

**matplot_axes : matplotlib.axes, default: None**  
Pass the plot data to an axes class object for plotting on a predefined canvas

**zero_col_centering : bool, default: False**  
Adjust colorscale with zero at the center of the scale

**interpolation : str, default: 'none'**  
The interpolation method for matplotlib.axes.Axes.imshow() to control down/upsampling.

## mapper.plots.plot_bar()

`mapper.plots.plot_bar(data=None, alignment=None, save_to_file=False, image_res=600, color=None, bar_title=None, xlim=None, ylim=None, figsize=[2,2], show_plot=False, matplot_axes=None)`

This function is the wrapper for matplotlib functions to visualize aggregated data array or 1D array with same length as the data matrix as a bar plot.

### Parameters
**data : list[float], numpy.ndarray: default: None**  
The aggregated data array or 1D-array to be plotted.

**save_to_file : bool, default: False**  
Save the output as PNG file. The filepath can be specified, using save_to_file = True would save the list to the same directory as the coordinate table or current working directory.

**image_res : int, default: 600**  
The image resolution of the plot in DPI

**color : str, default: None**  
The color for the aggregated dataplot

**bar_title : str, default: None**  
The title of the aggregated dataplot

**xlim : list[float, float], default: None**  
The maximum and minimum values of the x-axis of the aggregated dataplot

**ylim : list[float, float], default: None**  
The maximum and minimum values of the y-axis of the aggregated dataplot

**figsize : list[float, float], default: [2,2]**  
The size of the figure

**show_plot : bool, default: False**  
Render the aggregated plot as a standalone

**matplot_axes : matplotlib.axes, default: None**  
Pass the plot data to an axes class object for plotting on a predefined canvas

## mapper.plots.plot_colorbar()

`mapper.plots.plot_colorbar(cmap, vmin=None, vmax=None, data=None, alignment=None, save_to_file=False, image_res=600,cbar_label=None, cbar_title=None,tick_label=None, figsize=[2,2], show_plot=False, centering=False, step=None, orientation='vertical', matplot_axes=None)`

This function is the wrapper for matplotlib functions to visualize colormap pf the data matrix as a colorbar.

### Parameters
**cmap : matplotlib.colors.ListedColormap, matplotlib.colors.LinearSegmentedColormap, str**  
The colormap for the colorbar, can use colormaps from matplotlib, or custom built-in colormaps such as: custom_cmap.vlag_mpl, custom_cmap.vlag_r_mpl, custom_cmap.RdBu_s_mpl, custom_cmap.RdBu_sr_mpl, and custom_cmap.LimeGreens_mpl which emulate colormap from the package seaborn.

**vmin : float, default: None**  
The minimum value to be shown on the colorbar.

**vmax : float, default: None**  
The maximum value to be shown on the colorbar.

**data : list[numpy.ndarray], numpy.ndarray: default: None**  
The data matrix which is the data source for colorscale

**save_to_file : bool, default: False**  
Save the output as PNG file. The filepath can be specified, using save_to_file = True would save the list to the same directory as the coordinate table or current working directory.

**image_res : int, default: 600**  
The image resolution of the plot in DPI

**cbar_label : str, default: None**  
The label of the colorbar

**cbar_title : str, default: None**  
The title of the colorbar

**tick_label : list[str], default: None**  
The tick labels of the colorbar

**figsize : list[float, float], default: [2,2]**  
The size of the figure

**centering : bool, default: False**  
Adjust colorscale with zero at the center of the scale

**step : int, float, default: None**  
The incremental step of the colorscale. 

**orientation : str, default: 'vertical'**  
The orientation of the colorbar

**show_plot : bool, default: False**  
Render the colorbar plot as a standalone

**matplot_axes : matplotlib.axes, default: None**  
Pass the plot data to an axes class object for plotting on a predefined canvas

## mapper.plots.plot_annotation()

mapper.plots.plot_annotation(anno, cmap=None, ylim=None, save_to_file=False, image_res=600, figsize=[2,2], show_plot=False, matplot_axes=None)

This function is the wrapper for matplotlib functions to visualize the annotation of the data matrix.

### Parameters
**anno : list[float,int,str], numpy.ndarray, pandas.Series, pandas.DataFrame: default: None**  
The annotation of the rows of the data matrix

**cmap : matplotlib.colors.ListedColormap, matplotlib.colors.LinearSegmentedColormap, str**  
The colormap for the colorbar, can use colormaps from matplotlib, or custom built-in colormaps such as: custom_cmap.vlag_mpl, custom_cmap.vlag_r_mpl, custom_cmap.RdBu_s_mpl, custom_cmap.RdBu_sr_mpl, and custom_cmap.LimeGreens_mpl which emulate colormap from the package seaborn.

**ylim : list[float, float], default: None**  
The maximum and minimum values of the y-axis of the annotation bar

**save_to_file : bool, default: False**  
Save the output as PNG file. The filepath can be specified, using save_to_file = True would save the list to the same directory as the coordinate table or current working directory.

**image_res : int, default: 600**  
The image resolution of the plot in DPI

**figsize : list[float, float], default: [2,2]**  
The size of the figure

**show_plot : bool, default: False**  
Render the annotation bar as a standalone

**matplot_axes : matplotlib.axes, default: None**  
Pass the plot data to an axes class object for plotting on a predefined canvas

## mapper.plots.plot_basecount()

`mapper.plots.plot_basecount(alignment=None, save_to_file=False, image_res=600, basecount_title=None, xlim=None, ylim=None, figsize=[2,2], show_plot=False, matplot_axes=None)`

This function is the wrapper for matplotlib functions to visualize the base composition by x-axis of the data matrix.

### Parameters
**alignment : numpy.ndarray: default: None**  
The alignment matrix to be plotted.

**save_to_file : bool, default: False**  
Save the output as PNG file. The filepath can be specified, using save_to_file = True would save the list to the same directory as the coordinate table or current working directory.

**image_res : int, default: 600**  
The image resolution of the plot in DPI

**basecount_title : str, default: None**  
The title of the the basecount plot

**xlim : list[float, float], default: None**  
The maximum and minimum values of the x-axis of the basecount plot

**ylim : list[float, float], default: None**  
The maximum and minimum values of the y-axis of the basecount plot

**figsize : list[float, float], default: [2,2]**  
The size of the figure

**show_plot : bool, default: False**  
Render the basecount plot as a standalone

**matplot_axes : matplotlib.axes, default: None**  
Pass the plot data to an axes class object for plotting on a predefined canvas

## mapper.plots.plot_logos()

`mapper.plots.plot_logos(alignment=None, save_to_file=False, image_res=600, basecount_title=None, xlim=None, ylim=None, figsize=[2,2], show_plot=False, matplot_axes=None, yhighlights=None, yhighlight_col=None, yhighlight_alpha=None)`

This function is the wrapper for matplotlib functions to visualize the base composition by x-axis of the data matrix, in logos format (y-axis is converted to information instead of raw count)

### Parameters
**alignment : numpy.ndarray: default: None**  
The alignment matrix to be plotted.

**save_to_file : bool, default: False**  
Save the output as PNG file. The filepath can be specified, using save_to_file = True would save the list to the same directory as the coordinate table or current working directory.

**image_res : int, default: 600**  
The image resolution of the plot in DPI

**logos_title : str, default: None**  
The title of the the basecount plot

**xlim : list[float, float], default: None**  
The maximum and minimum values of the x-axis of the basecount plot

**ylim : list[float, float], default: None**  
The maximum and minimum values of the y-axis of the basecount plot

**figsize : list[float, float], default: [2,2]**  
The size of the figure

**show_plot : bool, default: False**  
Render the logos plot as a standalone

**matplot_axes : matplotlib.axes, default: None**  
Pass the plot data to an axes class object for plotting on a predefined canvas

**yhighlight : list[[int,int]], default: None**  
The list of coordinates of the vertical highlight on the logos plot area  

**yhighlight_col : list[str], default: ['grey']**  
The list of colors of the vertical highlight on the logos plot area

**yhighlight_alpha : list[float], default: [0.5]**  
The list of alpha values of the vertical highlight on the logos plot area 
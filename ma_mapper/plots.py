"""Figures for an overlay matrix: heatmaps, bars, annotations, base counts, logos.

Every function either draws into an axes you pass (``matplot_axes``) or makes
its own figure. ``save_to_file`` takes a path, or ``True`` for a default name
in the current directory; before 0.9.0 the default wrote into the installed
package directory, which is read-only in a normal install.
"""

import os

from ma_mapper import mapper
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
# Accessing functions and classes through matplotlib########
ListedColormap = matplotlib.colors.ListedColormap
BoundaryNorm = matplotlib.colors.BoundaryNorm
LinearSegmentedColormap = matplotlib.colors.LinearSegmentedColormap
MaxNLocator = matplotlib.ticker.MaxNLocator
ticker = matplotlib.ticker
MultipleLocator = matplotlib.ticker.MultipleLocator
FuncFormatter = matplotlib.ticker.FuncFormatter
TwoSlopeNorm = matplotlib.colors.TwoSlopeNorm
ScalarMappable = matplotlib.cm.ScalarMappable
#############################################################
import numpy as np
from typing import Literal
#%%


def _resolve_plot_path(save_to_file, default_name):
    """Where to write a figure: the given path, or ``default_name`` in the cwd."""
    if isinstance(save_to_file, str):
        return save_to_file
    return os.path.join(os.getcwd(), default_name)


#%%
def plot_heatmap(data: list|np.ndarray|None = None,
                 save_to_file: bool = False,
                 image_res: int = 600, 
                 cmap: ListedColormap|str|None = None,
                 plot_title: str|None = None,
                 title_fs=None,
                 xlim: list|None = None,
                 ylim: list|None = None,
                 xlabel = None,
                 ylabel = None,
                 xticklabels = None,
                 xticklabels_fs = None,
                 xticklabels_rt = None,
                 figsize:list=[2,2],
                 opacity: float = 1.0,
                 transparency_mode:str='constant',
                 show_plot:bool = False, 
                 matplot_axes = None,
                 zero_col_centering:bool=False,
                 interpolation = 'none',
                 **kwargs):
    if data is None or not isinstance(data, (list, np.ndarray)):
        raise ValueError("Invalid data")
    
    cmap_dict = {'nulc_white': ['grey', 'white','white','white','white'], 'dna': ['grey', 'green', 'yellow', 'blue', 'red'],'dna_jaspar': ['white', 'green', 'blue', 'red', 'yellow'], 'clinvar':['white','yellow','green','blue','orange','red']}

    if isinstance(cmap, ListedColormap|LinearSegmentedColormap):
        plot_cmap_a = cmap
    elif cmap in cmap_dict:
        zero_col_centering=False
        plot_cmap = ListedColormap(cmap_dict[cmap])
        kwargs['vmax'] = 5
        kwargs['vmin'] = 0
        plot_cmap_a = plot_cmap
    else:
        if opacity !=1:
            cmap = cmap or 'Blues' #assign Blues to nonetype 
            ncolors = 256
            plot_cmap = plt.get_cmap(cmap)(range(ncolors))
            if transparency_mode == 'gradient':
                plot_cmap[:,-1] = np.linspace(0, opacity, ncolors)
            elif transparency_mode == 'constant':
                plot_cmap[:,-1] = opacity
                plot_cmap[0,-1] = 0
            plot_cmap_a = ListedColormap(colors=plot_cmap)
        else:
            plot_cmap_a = cmap
    if zero_col_centering:
        vmax = kwargs.get('vmax', np.nanmax(data))
        vmin = kwargs.get('vmin', np.nanmin(data))
        norm = TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)
        kwargs['norm'] = norm
        kwargs.pop('vmax', None)
        kwargs.pop('vmin', None)

    ax = matplot_axes if matplot_axes is not None else plt.subplots(figsize=figsize)[1]
    heatmap = ax.imshow(data, aspect='auto', cmap=plot_cmap_a, interpolation=interpolation, **kwargs)
    if plot_title is not None:
        ax.set_title(plot_title, fontsize = title_fs)
    if xlim is not None:
        ax.set_xlim(xlim[0],xlim[1])
    if ylim is not None:
        ax.set_ylim(ylim[0],ylim[1])

    if xticklabels is not None:
        ax.set_xticklabels(xticklabels, fontsize=xticklabels_fs, rotation=xticklabels_rt)
    if xlabel is not None:
        ax.set_xlabel(xlabel)
    if ylabel is not None:
        ax.set_ylabel(ylabel)
    
    if save_to_file == True:
        plt.rcParams['savefig.dpi'] = image_res
        output_filepath = _resolve_plot_path(save_to_file, 'overlay.png')
        plt.savefig(output_filepath)
    if show_plot == True:
        plt.rcParams['figure.dpi'] = image_res
        plt.show()
    
# %%
def plot_bar(data: list|np.ndarray|None = None,
             save_to_file: bool = False,
             image_res: int = 600, 
             color: str|None = None,
             bar_title: str|None = None,
             xlim: list|None = None,
             ylim: list|None = None,
             figsize:list=[2,2],
             show_plot:bool = False, 
             matplot_axes = None,
             **kwargs):
    if data is None or not isinstance(data, (list, np.ndarray)):
        raise ValueError("Invalid data")
    dataplot = data
    if color is None:
        color = 'grey'
    opacity=kwargs.pop('opacity',1.0)
    ax = matplot_axes if matplot_axes is not None else plt.subplots(figsize=figsize)[1]
    #print(dataplot)
    bar = ax.bar(range(len(dataplot)), dataplot, color = color,alpha =opacity, **kwargs)
    ax.margins(x=0, y=0)
    if bar_title is not None:
        ax.set_title(bar_title)
    if xlim is not None:
        ax.set_xlim(xlim[0],xlim[1])
    if ylim is not None:
        ax.set_ylim(ylim[0],ylim[1])
    
    if save_to_file:
        if isinstance(save_to_file, str):
            output_filepath = save_to_file
        else:
            output_filepath = _resolve_plot_path(save_to_file, 'barplot.png')
        plt.rcParams['savefig.dpi'] = image_res
        plt.savefig(output_filepath)
    if show_plot == True:
        plt.rcParams['figure.dpi'] = image_res
        plt.show()
    
#%%
_ORIENT = Literal['horizontal','vertical']
def plot_colorbar(cmap: str|ListedColormap|LinearSegmentedColormap,
                  vmin:int|float|None = None,
                  vmax:int|float|None = None,
                  data:list|np.ndarray|None = None,
                  save_to_file: bool = False,
                  image_res: int = 600, 
                  cbar_label: str|None = None,
                  cbar_title: str|None = None,
                  tick_label: list|None = None,
                  figsize:list=[2,2],
                  show_plot:bool = False, 
                  centering:bool = False,
                  step: int|float|None = None,
                  orientation:_ORIENT = 'vertical',
                  matplot_axes = None,
                  **kwargs):
    cmap_dict = {'nulc_white': ['grey', 'white','white','white','white'], 'dna': ['grey', 'green', 'yellow', 'blue', 'red'],'dna_jaspar': ['white', 'green', 'blue', 'red', 'yellow']}
    if cmap is None:
        raise ValueError('need colormap')
    elif isinstance(cmap, ListedColormap|LinearSegmentedColormap):
        colorbar_cmap = cmap
    
    elif cmap in cmap_dict:
        colorbar_cmap = ListedColormap(cmap_dict[cmap])
        vmax = int(5)
        vmin = int(0)
        kwargs['ticks'] = [0.5, 1.5, 2.5, 3.5, 4.5]
    else:
        ncolors = 256
        plot_cmap = plt.get_cmap(cmap)(range(ncolors))
        colorbar_cmap = ListedColormap(colors=plot_cmap)

    if vmin is None and vmax is None:
        if data is not None:
            vmax = np.nanmax(data)
            vmin = np.nanmin(data)
        else:
            raise ValueError('insufficient input, please provide either vmin+vmax or data')
    if centering:
        norm = TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)
    else:
        if step is None:
            if (vmax is int and vmin is int) or cmap in cmap_dict:
                step =1
            else:
                step =0.1
        #print(step)
        bounds =np.arange(vmin,vmax+step,step)
        extend_arg = 'neither'
        if data is not None:
            if np.nanmax(data) > vmax and np.nanmin(data) < vmin:
                extend_arg = 'both'
            elif np.nanmax(data) > vmax:
                extend_arg = 'max'
            elif np.nanmin(data) < vmin:
                extend_arg = 'min'
            
        norm = BoundaryNorm(bounds, colorbar_cmap.N, extend=extend_arg)
    if matplot_axes is not None:
        fig = matplot_axes.figure
        ax = matplot_axes
    else:
        fig, ax = plt.subplots(figsize=figsize, layout='constrained')
    cbar_obj=fig.colorbar(ScalarMappable(norm=norm, cmap=colorbar_cmap),cax=ax, orientation=orientation,label=cbar_label,**kwargs)
    if cbar_title:
        ax.set_title(cbar_title, fontsize='small')
    if tick_label:
        cbar_obj.set_ticklabels(tick_label, fontsize='small')
    
    if save_to_file:
        plt.rcParams['savefig.dpi'] = image_res
        # output_filepath was never assigned before this test, so a save
        # raised NameError.
        output_filepath = _resolve_plot_path(save_to_file, 'colorbar.png')
        plt.savefig(output_filepath)
    if show_plot == True:
        plt.rcParams['figure.dpi'] = image_res
        plt.show()
    
#%%
def plot_annotation(anno: list|np.ndarray|pd.Series|pd.DataFrame,
                    cmap: str|ListedColormap|LinearSegmentedColormap|None=None,
                    ylim: list|None = None,
                    save_to_file: bool = False,
                    image_res: int = 600, 
                    figsize:list=[2,2],
                    show_plot:bool = False,
                    matplot_axes = None,
                    **kwargs):
    if isinstance(anno, list|np.ndarray):
        anno_plot = np.array(anno)
        anno_uniq=np.unique(anno_plot.sort())
    elif isinstance(anno, pd.Series):
        anno_plot = anno.values
        anno_uniq=anno.sort_values().unique()
    if isinstance(cmap, ListedColormap|LinearSegmentedColormap):
        #print('check1')
        anno_cmap = cmap
    elif isinstance(cmap, list):
        #print('check2')
        anno_cmap = ListedColormap(cmap)
        boundaries = np.arange(len(cmap) + 1) - 0.
        num_colors = len(anno_uniq)
        norm = BoundaryNorm(boundaries, num_colors)
        
    else:
        #print('check3')
        cmap = cmap or 'Blues'
        anno_cmap_= plt.get_cmap(cmap)
        num_colors = len(anno_uniq)
        colors = [anno_cmap_(i / (num_colors - 1)) for i in range(num_colors)]
        anno_cmap = ListedColormap(colors)
        boundaries = np.append(anno_uniq, anno_uniq[-1] + (anno_uniq[-1] - anno_uniq[-2]))
        
        norm = BoundaryNorm(boundaries, num_colors)

    ax = matplot_axes if matplot_axes is not None else plt.subplots(figsize=figsize)[1]
    annoation = ax.imshow(anno_plot.reshape(-1, 1), aspect = 'auto',interpolation='none',cmap= anno_cmap,norm=norm)
    ax.margins(x=0, y=0)
    if ylim is not None:
        ax.set_ylim(ylim[0],ylim[1])
    if save_to_file:
        if isinstance(save_to_file, str):
            output_filepath = save_to_file
        else:
            output_filepath = _resolve_plot_path(save_to_file, 'annot_bar.png')
        plt.rcParams['savefig.dpi'] = image_res
        plt.savefig(output_filepath)
    if show_plot == True:
        plt.rcParams['figure.dpi'] = image_res
        plt.show()
#%%
def format_func(value, tick_number):
    if value < 0.1:
        return f"{value:.0e}"
    else:
        return f"{value:.2f}"
formatter = ticker.FuncFormatter(format_func)
#%%
def plot_basecount(alignment:np.ndarray|None,
                   save_to_file: bool = False,
                   image_res: int = 600,
                   basecount_title: str|None = None,
                   xlim: list|None = None,
                   ylim: list|None = None,
                   figsize:list=[2,2],
                   show_plot:bool = False,
                   matplot_axes=None,
                   **kwargs):
    col_dict_dna={'-': 'white','A':'green','C': 'blue','T': 'red','G': 'yellow',}
    base_count=mapper.base_count(alignment=alignment)
    ax = matplot_axes if matplot_axes is not None else plt.subplots(figsize=figsize)[1]
    bottom = np.zeros(base_count.shape[1])
    for values, label in zip(base_count, ['-','A','C','T','G']):
        ax.bar(np.arange(base_count.shape[1]), values, label=label, bottom=bottom, color=col_dict_dna[label], **kwargs)
        bottom += values
        ax.margins(x=0, y=0)
    if basecount_title is not None:
        ax.set_title(basecount_title)
    if xlim is not None:
        ax.set_xlim(xlim[0],xlim[1])
    if ylim is not None:
        ax.set_ylim(ylim[0],ylim[1])
   
    if save_to_file:
        plt.rcParams['savefig.dpi'] = image_res
        # output_filepath was never assigned before this test, so a save
        # raised NameError.
        output_filepath = _resolve_plot_path(save_to_file, 'barplot.png')
        plt.savefig(output_filepath)
    if show_plot == True:
        plt.rcParams['figure.dpi'] = image_res
        plt.show()
    
#%%
def plot_logos(alignment:np.ndarray|None,
                save_to_file: bool = False,
                image_res: int = 600,
                figsize:list=[2,2],
                matplot_axes=None,
                logos_title: str|None = None,
                xlim: list|None = None,
                ylim: list|None = None,
                show_plot:bool = False,
                yhighlights: list|None = None,
                yhighlight_col: list|None = None,
                yhighlight_alpha: list|None = None,
                **kwargs):
    import logomaker
    base_count=mapper.base_count(alignment=alignment)
    bese_count_nogap=base_count[1:]
    counts_mat=bese_count_nogap.transpose()
    column_names = ['A', 'C', 'T', 'G']
    counts_mat_df = pd.DataFrame(counts_mat, columns=column_names)
    info_mat = logomaker.transform_matrix(counts_mat_df, from_type='counts', to_type='information')
    ax = matplot_axes if matplot_axes is not None else plt.subplots(figsize=figsize)[1]
    logos=logomaker.Logo(info_mat, ax=ax, **kwargs)
    if yhighlights is not None:
        for idx,yhighlight in enumerate(yhighlights):
            logos.highlight_position_range(pmin=yhighlight[0]+0.5, pmax=yhighlight[1]-0.5, color = yhighlight_col[idx], alpha=yhighlight_alpha[idx])
    ax.margins(x=0, y=0)
    if logos_title is not None:
        ax.set_title(logos_title)
    if xlim is not None:
        ax.set_xlim(xlim[0],xlim[1])
    if ylim is not None:
        ax.set_ylim(ylim[0],ylim[1])
    
    if save_to_file:
        plt.rcParams['savefig.dpi'] = image_res
        # output_filepath was never assigned before this test, so a save
        # raised NameError.
        output_filepath = _resolve_plot_path(save_to_file, 'barplot.png')
        plt.savefig(output_filepath)
    if show_plot == True:
        plt.rcParams['figure.dpi'] = image_res
        plt.show()

#%%
_HEATMAP_MODE = Literal['overlay','spread_horizontal','spread_vertical']
def plot(data: list|np.ndarray|None = None,
                      show_alignment:bool = False,
                      alignment: np.ndarray|None = None,
                      alignment_col: str = 'nulc_white',
                      
                    heatmap:bool = True,
                    heatmap_mode:_HEATMAP_MODE = 'overlay',
                    vlim: list|None = None,
                    heatmap_color: ListedColormap|str|None = None,
                    heatmap_title: list|None = None,
                    heatmap_title_fs: int= None,
                    heatmap_yhighlight: list = None,
                    heatmap_yhighlight_col: list = ['grey'],
                    heatmap_yhighlight_alpha: list = [0.5],
                    heatmap_xhighlight: list = None,
                    heatmap_xhighlight_col: list = ['grey'],
                    heatmap_xhighlight_alpha: list = [0.5],
                    heatmap_ylabel:list|None = None,
                    heatmap_ylabel_fs: int = 6,
                    heatmap_xlabel:list|None = None,
                    heatmap_xlabel_fs: int = 6,
                    
                    aggregated:bool = False,
                    aggregated_data: list|None=None,
                    agg_h:int = 10,
                    agg_ylim:list|None = [None],
                    agg_ylabel:list|None = None,
                    agg_ylabel_fs: int = 6,
                    agg_ylabel_ypos: list = None,
                    agg_ylabel_xpos: list = None,
                    agg_yhighlight: list = None,
                    agg_yhighlight_col: list = ['grey'],
                    agg_yhighlight_alpha: list = [0.5],
                    agg_xhighlight: list = None,
                    agg_xhighlight_col: list = ['grey'],
                    agg_xhighlight_alpha: list = [0.5],
                    agg_xlabel =None,
                    agg_xlabel_fs =None,
                    agg_major_tick: int = 50,
                    agg_yscale: list|None = None,
                    agg_ylabel_right: list|None = None,
                    agg_ylabel_right_fs: int = 6,
                    agg_ylabel_right_pos: list = [1.1,0.5],
                    agg_plottext: list|None = None,
                    agg_plottext_fs: int = 6,
                    agg_plottext_pos: list = [0.99,0.90],
                    agg_plot_title: list = None,
                    agg_plot_title_fs:int = None,
                    agg_colset: str|None = None,
                    logos:bool = False,
                    base_count:bool = False,

                    annotation:bool = False,
                    annotation_data: list|None= None,
                    anno_title: list|None = None,
                    anno_w:int = 2,
                    anno_ylabel: str = None,
                    anno_ylabel_fs: int = 6,
                    anno_col: list|None = ['Blues'],
                    anno_cbar: bool = True,
                    anno_cbar_title: list|None = None,
                    anno_cbar_label: list|None = None,
                    anno_cbar_even_pos: float= -0.25,
                    
                    show_alignment_colbar: bool = False,
                    colorbar:bool = False,
                    colorbar_width:int = 2,
                    colorbar_steps: list|None = None,

                    image_res: int = 600, 
                    figsize:list=[50,50],
                    show_plot:bool = True,
                    save_to_file: bool = False, 
                    **kwargs):
    plt.rcParams['savefig.dpi'] = image_res
    plt.rcParams['figure.dpi'] = image_res
    fig = plt.figure(figsize=figsize)
    grid = fig.add_gridspec(nrows = 500, ncols = 500, hspace=0)
    global_kwargs = {key: value for key, value in kwargs.items() if not key.startswith(('ag_', 'hm_', 'cb_','bc_','lg_'))}
    ag_kwargs = {key[3:]: value for key, value in kwargs.items() if key.startswith('ag_')}
    hm_kwargs = {key[3:]: value for key, value in kwargs.items() if key.startswith('hm_')}
    cb_kwargs = {key[3:]: value for key, value in kwargs.items() if key.startswith('cb_')}
    bc_kwargs = {key[3:]: value for key, value in kwargs.items() if key.startswith('bc_')}
    lg_kwargs = {key[3:]: value for key, value in kwargs.items() if key.startswith('lg_')}
    ##########grid configuration####################
    # Default grid positions
    heatmap_grid_defaults = {'x1': 50, 'x2': 100, 'y1': 0, 'y2': 50}
    agg_grid_y1_default = 50
    agg_grid_defaults = {'x1': 50, 'x2': 100, 'y1': agg_grid_y1_default, 'y2': agg_grid_y1_default+agg_h}
    anno_grid_defaults = {'x1': 48, 'x2': 50, 'y1': 0, 'y2': 50}
    ################################################
    def format_func(value, tick_number):
        if -0.01 <= value <= 0.1  :
            return f"{value:.0e}"
        elif value >= 1000000:
            return f"{value / 1000000:.1f}M"
        elif value >= 1000:
            return f"{value / 1000:.1f}k"
        else:
            return f"{value:.2f}"
    ################################################
    total_data_count = len(data) if data is not None else 0
    if heatmap:
        heatmap_grid = {
            'x1': heatmap_grid_defaults['x1'],
            'x2': heatmap_grid_defaults['x2'],
            'y1': heatmap_grid_defaults['y1'],
            'y2': heatmap_grid_defaults['y2']
        }

        
        if alignment is not None and show_alignment:
            total_data_count += 1
        if heatmap_mode == 'overlay':
            heatmap_ax = fig.add_subplot(grid[heatmap_grid['y1']:heatmap_grid['y2'], heatmap_grid['x1']:heatmap_grid['x2']])
            if alignment is not None and show_alignment:
                plot_heatmap(alignment, cmap=alignment_col, matplot_axes=heatmap_ax, vmin=0, vmax=5, **hm_kwargs, **global_kwargs)
            if data is not None:
                for idx, overlay in enumerate(data):
                    plot_heatmap(overlay, cmap=heatmap_color[idx], matplot_axes=heatmap_ax, vmin=vlim[idx % len(vlim)][0], vmax=vlim[idx % len(vlim)][1], **hm_kwargs, **global_kwargs)
            heatmap_ax.tick_params(axis='both', which='major', labelsize=8)
            heatmap_ax.xaxis.set_major_locator(MaxNLocator(integer=True))
            heatmap_ax.xaxis.set_major_locator(MultipleLocator(agg_major_tick))
            #heatmap_ax.invert_yaxis()  
            #heatmap_ax.yaxis.set_major_formatter(FuncFormatter(format_func))
            if annotation:
                heatmap_ax.set_yticks([])
            if heatmap_title is not None:
                heatmap_ax.set_title(heatmap_title[0], fontsize= heatmap_title_fs)
            if heatmap_xlabel is not None:
                heatmap_ax.set_xlabel(heatmap_xlabel, fontsize= heatmap_xlabel_fs)
            if heatmap_ylabel is not None:
                heatmap_ax.set_ylabel(heatmap_ylabel, fontsize= heatmap_ylabel_fs)
            if aggregated:
                heatmap_ax.set_xticks([])
            if heatmap_yhighlight is not None:
                for _idx, yhighlight_coord in enumerate(heatmap_yhighlight):
                    heatmap_ax.axvspan(yhighlight_coord[0], yhighlight_coord[1], color=heatmap_yhighlight_col[(_idx % len(heatmap_yhighlight))], alpha=heatmap_yhighlight_alpha[(_idx % len(heatmap_yhighlight))])
            if heatmap_xhighlight is not None:
                for _idx, xhighlight_coord in enumerate(heatmap_xhighlight):
                    heatmap_ax.axhspan(xhighlight_coord[0], xhighlight_coord[1], color=heatmap_xhighlight_col[(_idx % len(heatmap_xhighlight))], alpha=heatmap_xhighlight_alpha[(_idx % len(heatmap_xhighlight))])
        elif heatmap_mode in ['spread_horizontal','spread_vertical']:
            heatmap_axes = []
            for idx in range(total_data_count):
                if heatmap_mode == 'spread_horizontal':
                    grid_x1 = heatmap_grid_defaults['x1'] + 60 * idx
                    grid_x2 = heatmap_grid_defaults['x2'] + 60 * idx
                    heatmap_axes.append(fig.add_subplot(grid[heatmap_grid['y1']:heatmap_grid['y2'], grid_x1:grid_x2]))
                elif heatmap_mode == 'spread_vertical':
                    grid_y1 = heatmap_grid_defaults['y1'] + 65 * idx
                    grid_y2 = heatmap_grid_defaults['y2'] + 65 * idx
                    heatmap_axes.append(fig.add_subplot(grid[grid_y1:grid_y2, heatmap_grid['x1']:heatmap_grid['x2']]))
            if alignment is not None and show_alignment:
                plot_heatmap(alignment, cmap=alignment_col, matplot_axes=heatmap_axes[0], vmin=0, vmax=5, **hm_kwargs, **global_kwargs)
                heatmap_axes[0].tick_params(axis='both', which='major', labelsize=8)
                heatmap_axes[0].xaxis.set_major_locator(MaxNLocator(integer=True))
                heatmap_axes[0].xaxis.set_major_locator(MultipleLocator(agg_major_tick))
                #heatmap_axes[0].invert_yaxis()
                if annotation:  
                    heatmap_axes[0].set_yticks([])
                    
                #else:
                #    heatmap_axes[0].yaxis.set_major_formatter(FuncFormatter(format_func))
            if data is not None:
                for idx, overlay in enumerate(data):
                    grid_idx = idx + 1 if alignment is not None and show_alignment else idx
                    plot_heatmap(overlay, cmap=heatmap_color[idx], matplot_axes=heatmap_axes[grid_idx], vmin=vlim[idx % len(vlim)][0], vmax=vlim[idx % len(vlim)][1], **hm_kwargs, **global_kwargs)
                    heatmap_axes[grid_idx].tick_params(axis='both', which='major', labelsize=8)
                    heatmap_axes[grid_idx].xaxis.set_major_locator(MaxNLocator(integer=True))
                    heatmap_axes[grid_idx].xaxis.set_major_locator(MultipleLocator(agg_major_tick))  
                    #heatmap_axes[grid_idx].invert_yaxis()
                    if heatmap_title is not None and heatmap_title[idx] is not None:
                        heatmap_axes[grid_idx].set_title(heatmap_title[idx])
                    if idx == 0 and aggregated:
                        heatmap_axes[grid_idx].set_xticks([])
                    if idx == 0 and annotation:
                        heatmap_axes[grid_idx].set_yticks([])
                    #else:
                    #    heatmap_axes[grid_idx].yaxis.set_major_formatter(FuncFormatter(format_func))
    if aggregated:
        aggregated_grids = []
        for idx in range(len(aggregated_data)):
            agg_grid_y1 = agg_grid_defaults['y1'] + agg_h * idx
            agg_grid_y2 = agg_grid_defaults['y2'] + agg_h * idx
            aggregated_grids.append(fig.add_subplot(grid[agg_grid_y1:agg_grid_y2, agg_grid_defaults['x1']:agg_grid_defaults['x2']]))
        if agg_ylabel is not None and agg_ylabel_ypos is None:
            agg_ylabel_ypos = np.full(len(agg_ylabel), None)
        for idx, agg_data in enumerate(aggregated_data):
            if isinstance(agg_data, list):
                for _idx,sub_agg_data in enumerate(agg_data):
                    if isinstance(agg_colset, list):
                        set_col = agg_colset[idx][_idx % len(agg_data)]
                    else:
                        set_col = agg_colset[idx]
                    plot_bar(data=sub_agg_data, color=set_col, matplot_axes=aggregated_grids[idx], ylim=agg_ylim[idx % len(agg_ylim)], **ag_kwargs, **global_kwargs)
            else:
                plot_bar(data=agg_data, color=agg_colset[idx], matplot_axes=aggregated_grids[idx], ylim=agg_ylim[idx % len(agg_ylim)], **ag_kwargs, **global_kwargs)
            aggregated_grids[idx].tick_params(axis='both', which='major', labelsize=8)
            aggregated_grids[idx].xaxis.set_major_locator(MultipleLocator(agg_major_tick))
            if agg_yhighlight is not None:
                for _idx, yhighlight_coord in enumerate(agg_yhighlight):
                    aggregated_grids[idx].axvspan(yhighlight_coord[0], yhighlight_coord[1], color=agg_yhighlight_col[(_idx % len(agg_yhighlight))], alpha=agg_yhighlight_alpha[(_idx % len(agg_yhighlight))])
            if agg_xhighlight is not None:
                for _idx, xhighlight_coord in enumerate(agg_xhighlight):
                    aggregated_grids[idx].axhspan(xhighlight_coord[0], xhighlight_coord[1], color=agg_xhighlight_col[(_idx % len(agg_xhighlight))], alpha=agg_xhighlight_alpha[(_idx % len(agg_xhighlight))])
            if agg_ylabel_right is not None and agg_ylabel_right[idx] is not None:
                aggregated_grids[idx].text(agg_ylabel_right_pos[0], agg_ylabel_right_pos[1], agg_ylabel_right[idx], 
                transform=aggregated_grids[idx].transAxes, 
                fontsize=agg_ylabel_right_fs, 
                verticalalignment='bottom', 
                horizontalalignment='right')
            if agg_plottext is not None and agg_plottext[idx] is not None:
                aggregated_grids[idx].text(agg_plottext_pos[0], agg_plottext_pos[1], agg_plottext[idx], 
                transform=aggregated_grids[idx].transAxes, 
                fontsize=agg_plottext_fs, 
                verticalalignment='bottom', 
                horizontalalignment='right')
            if agg_plot_title is not None and agg_plot_title[idx] is not None:
                aggregated_grids[idx].set_title(agg_plot_title[idx], fontsize=agg_plot_title_fs)
            if agg_ylabel is not None and agg_ylabel[idx] is not None:
                if agg_ylabel_ypos[idx] is not None:
                    aggregated_grids[idx].set_ylabel(agg_ylabel[idx], fontsize=agg_ylabel_fs,x=agg_ylabel_xpos,y=agg_ylabel_ypos[idx])
                else:
                    aggregated_grids[idx].set_ylabel(agg_ylabel[idx], fontsize=agg_ylabel_fs,x=agg_ylabel_xpos)
            if agg_xlabel is not None:
                aggregated_grids[idx].set_xlabel(agg_xlabel, fontsize=agg_xlabel_fs)
            if idx != len(aggregated_data)-1 or logos:
                aggregated_grids[idx].set_xticks([])
            if agg_yscale is not None:
                for _idx, yscale in enumerate(agg_yscale):
                    aggregated_grids[idx].set_yscale(yscale)
    if base_count:
        agg_grid_y1 = agg_grid_defaults['y1'] + agg_h if 'agg_grid_y1' not in locals() else agg_grid_y1 + agg_h
        agg_grid_y2 = agg_grid_defaults['y2'] + agg_h if 'agg_grid_y2' not in locals() else agg_grid_y2 + agg_h
        base_count_grid = fig.add_subplot(grid[agg_grid_y1:agg_grid_y2, agg_grid_defaults['x1']:agg_grid_defaults['x2']])
        base_count_grid.xaxis.set_major_locator(MultipleLocator(agg_major_tick))
        plot_basecount(alignment=alignment, matplot_axes=base_count_grid, **bc_kwargs, **global_kwargs)
        base_count_grid.tick_params(axis='both', which='major', labelsize=8)
    if logos:
        agg_grid_y1 = agg_grid_defaults['y1'] + agg_h if 'agg_grid_y1' not in locals() else agg_grid_y1 + agg_h
        agg_grid_y2 = agg_grid_defaults['y2'] + agg_h if 'agg_grid_y2' not in locals() else agg_grid_y2 + agg_h
        logos_grid = fig.add_subplot(grid[agg_grid_y1:agg_grid_y2, agg_grid_defaults['x1']:agg_grid_defaults['x2']])
        logos_grid.xaxis.set_major_locator(MultipleLocator(agg_major_tick))
        plot_logos(alignment=alignment, matplot_axes=logos_grid, yhighlights=agg_yhighlight, yhighlight_col=agg_yhighlight_col, yhighlight_alpha=agg_yhighlight_alpha, **lg_kwargs, **global_kwargs)
        logos_grid.tick_params(axis='both', which='major', labelsize=8)
        if agg_xlabel is not None:
                logos_grid.set_xlabel(agg_xlabel, fontsize=agg_xlabel_fs)
    if colorbar:
        if heatmap_mode == 'overlay':
            cbar_axes = []
            if data is not None and show_alignment_colbar == False:
                total_data_count=total_data_count
            for idx in range(total_data_count):
                grid_x1 = heatmap_grid_defaults['x2'] + 7 * idx
                grid_x2 = heatmap_grid_defaults['x2'] + colorbar_width + 7 * idx
                cbar_axes.append(fig.add_subplot(grid[heatmap_grid_defaults['y1']:heatmap_grid_defaults['y2'], grid_x1:grid_x2]))
            if show_alignment_colbar == True:
                plot_colorbar(cmap=alignment_col, data=alignment, matplot_axes=cbar_axes[0], **cb_kwargs)
                cbar_axes[0].tick_params(axis='both', which='major', labelsize=8)
                cbar_axes[0].set_yticklabels(['gap', 'A', 'C', 'T', 'G'], fontsize='small')
            if data is not None:
                for idx, overlay in enumerate(data):
                    print(idx)
                    grid_idx = idx + 1 if show_alignment_colbar==True else idx
                    plot_colorbar(cmap=heatmap_color[idx], data=overlay, vmin=vlim[idx][0], vmax=vlim[idx][1], matplot_axes=cbar_axes[grid_idx], step=colorbar_steps[idx], **cb_kwargs)
                    cbar_axes[idx].tick_params(axis='both', which='major', labelsize=8)
                    #cbar_axes[idx].yaxis.set_major_formatter(FuncFormatter(format_func))
        elif heatmap_mode in ['spread_horizontal','spread_vertical']:
            cbar_axes = []
            for idx in range(total_data_count):
                if heatmap_mode == 'spread_horizontal':
                    grid_x1 = heatmap_grid_defaults['x2'] + 60 * idx
                    grid_x2 = heatmap_grid_defaults['x2'] + colorbar_width + 60 * idx
                    cbar_axes.append(fig.add_subplot(grid[heatmap_grid['y1']:heatmap_grid['y2'], grid_x1:grid_x2]))
                elif heatmap_mode == 'spread_vertical':
                    grid_y1 = heatmap_grid_defaults['y1'] + 65 * idx
                    grid_y2 = heatmap_grid_defaults['y2'] + 65 * idx
                    cbar_axes.append(fig.add_subplot(grid[grid_y1:grid_y2, heatmap_grid['x2']:heatmap_grid['x2']+colorbar_width]))
            if show_alignment_colbar==True:
                plot_colorbar(cmap=alignment_col, data=alignment, matplot_axes=cbar_axes[0], **cb_kwargs)
                cbar_axes[0].tick_params(axis='both', which='major', labelsize=8)
                cbar_axes[0].set_yticklabels(['gap', 'A', 'C', 'T', 'G'], fontsize='small')
            if data is not None:
                for idx, overlay in enumerate(data):
                    grid_idx = idx + 1 if show_alignment_colbar==True else idx
                    plot_colorbar(cmap=heatmap_color[idx], data=overlay, vmin=vlim[idx][0], vmax=vlim[idx][1], matplot_axes=cbar_axes[grid_idx], step = colorbar_steps[idx],**cb_kwargs)
                    cbar_axes[idx].tick_params(axis='both', which='major', labelsize=8,)
                    #cbar_axes[idx].yaxis.set_major_formatter(FuncFormatter(format_func))
    if annotation:
        anno_grids = []
        for idx in range(len(annotation_data)):
            anno_grid_x1 = anno_grid_defaults['x1'] - anno_w * idx
            anno_grid_x2 = anno_grid_defaults['x2'] - anno_w * idx
            anno_grids.append(fig.add_subplot(grid[anno_grid_defaults['y1']:anno_grid_defaults['y2'], anno_grid_x1:anno_grid_x2]))
        #print('anno_grid_x1', anno_grid_x1) 
        #anno_cbar_grid_x2 = anno_grid_x1 - 7
        #anno_cbar_grid_x1 = anno_cbar_grid_x2 - 20
        for idx, annotation in enumerate(annotation_data):
            anno_grids[idx].tick_params(axis='both', which='major', labelsize=8)
            anno_grids[idx].set_xticks([])
            if idx < len(annotation_data) - 1:
                anno_grids[idx].set_yticks([])
            elif idx == len(annotation_data) - 1 and anno_ylabel is not None:
                anno_grids[idx].set_ylabel(anno_ylabel, fontsize = anno_ylabel_fs)
            anno_cmap, norm, num_colors, boundaries =plot_annotation(annotation, matplot_axes=anno_grids[idx], cmap=anno_col[idx], export_mats=True,**global_kwargs)
            if anno_title is not None and anno_title[idx] is not None:
                anno_grids[idx].set_title(anno_title[idx], fontsize=6, rotation=90)

            #anno_grids[idx].yaxis.set_major_formatter(FuncFormatter(format_func))
            #anno_grids[idx].invert_yaxis()
            if anno_cbar:
                #anno_cbar_grid_y1 = anno_grid_defaults['y1'] +5*idx
                #anno_cbar_grid_y2 = anno_grid_defaults['y1'] + 2 + 5*idx
                anno_cbar_grid_x2=int(anno_grid_x2 - 9 - 7*np.round(idx/2))
                anno_cbar_grid_x1= anno_cbar_grid_x2-2
                anno_cbar_grid_y1 = int(anno_grid_defaults['y1'] + 25*(idx % 2))
                anno_cbar_grid_y2 = int(anno_cbar_grid_y1) + 20
                #print(anno_cbar_grid_y1,anno_cbar_grid_y2,anno_cbar_grid_x1,anno_cbar_grid_x2)
                anno_cbar_grid = fig.add_subplot(grid[anno_cbar_grid_y1:anno_cbar_grid_y2, anno_cbar_grid_x1:anno_cbar_grid_x2])
                tick_positions = (boundaries[:-1] + boundaries[1:]) / 2
                anno_cbar=fig.colorbar(ScalarMappable(norm=norm, cmap=anno_cmap),cax=anno_cbar_grid, orientation='vertical', ticks= tick_positions)
                anno_cbar_grid.tick_params(axis='both', which='major', labelsize=8)
                anno_cbar_grid.yaxis.set_ticks_position('left')
                #anno_cbar_grid.invert_yaxis()
                #anno_cbar.ax.yaxis.set_major_formatter(FuncFormatter(format_func))
                if anno_cbar_title is not None and anno_cbar_title[idx] is not None:
                    if idx % 2 == 1:
                        anno_cbar_grid.set_title(anno_cbar_title[idx], fontsize='small', y=anno_cbar_even_pos)
                    else: 
                        anno_cbar_grid.set_title(anno_cbar_title[idx], fontsize='small')
                if anno_cbar_label is not None and anno_cbar_label[idx] is not None:
                    anno_cbar.ax.set_yticklabels(anno_cbar_label[idx])
                    #anno_cbar.ax.yaxis.set_major_formatter(FuncFormatter(format_func))
    if save_to_file:
        plt.rcParams['savefig.dpi'] = image_res
        plt.rcParams['savefig.bbox'] = 'tight'
        plt.rcParams['savefig.pad_inches'] = 0
        if isinstance(save_to_file, str):
            output_filepath = save_to_file
        else:
            output_filepath = _resolve_plot_path(save_to_file, 'overlay_composite.png')
        plt.savefig(output_filepath)
        if show_plot == False:
            plt.close()
    if show_plot == True:
        plt.rcParams['figure.dpi'] = image_res
        plt.show()
        plt.close()
#%%

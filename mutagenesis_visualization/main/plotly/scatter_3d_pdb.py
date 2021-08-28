"""
This module contains the plotly 3D scatter plot where you can add PDB
properties.
"""
from pathlib import Path
from typing import Any, Dict, List, Union, Optional
from plotly import express as px

from mutagenesis_visualization.main.classes.base_model_plotly import Plotly
from mutagenesis_visualization.main.utils.plotly_utils import (
    color_3d_scatter,
    parse_pdb_coordinates,
    update_axes,
    update_layout,
    matplotlib_to_plotly,
)


class Scatter3DPDB(Plotly):
    """
    This class uses plotly to generate a 3D scatter plot of the protein
    and the enrichment scores where you can add PDB properties.
    """
    def __call__(
        self,
        pdb_path: str = None,
        plot: Optional[List[str]] = None,
        mode: str = 'mean',
        custom: Any = None,
        position_correction: int = 0,
        chain: str = 'A',
        replicate: int = -1,
        output_html: Union[None, str, Path] = None,
        **kwargs: Any,
    ) -> None:
        """
        Generates a 3-D scatter plot of different properties obtained from
        the PDB. PDBs may have atoms missing, you should fix the PDB before
        using this method. We recommend you use matplotlib for interactive plot.

        Parameters
        -----------
        pdb_path : str, default None
            User should specify the path PDB.

        plot : list, default ['Distance', 'SASA', 'log B-factor']
            List of 3 elements to plot. Other options are 'Score' and Custom.
            If custom, add the label to the third element of the list ie.
            ['Distance', 'SASA', 'Conservation'].

        mode : str, default 'mean'
            Specify what enrichment scores to use. If mode = 'mean', it will
            use the mean of each position to classify the residues.
            If mode = 'A', it will use the Alanine substitution profile.
            Can be used for each amino acid. Use the one-letter code and
            upper case.

        custom : list or dataframe or np.array, default None
            If you want to add a custom dataset to plot, use custom. On the
            parameter plot, the 3rd item of the list will be the label for
            your custom dataset.

        df_color : pandas dataframe, default None
            The color of each residue can also be included. You must label
            that label column.

        color_by_score : boolean, default True
            If set to False, the points in the scatter will not be colored
            based on the enrichment score.

        position_correction : int, default 0
            If the pdb structure has a different numbering of positions than
            you dataset, you can correct for that. If your start_position = 2,
            but in the PDB that same residue is at position 20,
            position_correction needs to be set at 18.

        chain : str, default 'A'
            Chain of the PDB file to get the coordinates and SASA from.

        replicate : int, default -1
            Set the replicate to plot. By default, the mean is plotted.
            First replicate start with index 0.
            If there is only one replicate, then leave this parameter
            untouched.

        output_html : str, default None
            If you want to export the generated graph into html, add the
            path and name of the file.
            Example: 'path/filename.html'.

        **kwargs : other keyword arguments
        """
        temp_kwargs: Dict[str, Any] = self._update_kwargs(kwargs)
        temp_kwargs['x_label'] = kwargs.get('x_label', plot[0])
        temp_kwargs['y_label'] = kwargs.get('y_label', plot[1])
        temp_kwargs['z_label'] = kwargs.get('z_label', plot[2])

        if not plot:
            plot = ['Distance', 'SASA', 'log B-factor']

        # Get Scores and colors
        df_scores = color_3d_scatter(
            self.dataframes.df_notstopcodons[replicate],
            mode,
            temp_kwargs['lof'],
            temp_kwargs['gof'],
        )

        # If coordinates is not an input, get it from the pdb
        self.df_output = parse_pdb_coordinates(
            pdb_path,
            self.start_position,
            self.end_position,
            position_correction,
            chain,
            sasa=True,
        )

        # Add scores
        self.df_output['Score'] = list(df_scores['Score'])

        # Custom data
        if custom is not None:
            self.df_output[plot[2]] = custom

            # Plot figure
        self.fig = px.scatter_3d(
            self.df_output,
            x=plot[0],
            y=plot[1],
            z=plot[2],
            color='Score',
            color_continuous_scale=matplotlib_to_plotly(temp_kwargs['colormap']),
            range_color=temp_kwargs['colorbar_scale'],
        )

        self._tune_plot(temp_kwargs)
        self._save_html(output_html, temp_kwargs)

    def _tune_plot(self, temp_kwargs: Dict[str, Any]) -> None:
        """
        Change stylistic parameters of the plot.
        """
        # update axes
        update_axes(self.fig, temp_kwargs)

        # for the clickable part
        self.fig.update_traces(
            hovertext=self.df_output['Position'], hovertemplate='Position: %{hovertext}'
        )
        # title
        update_layout(self.fig, temp_kwargs)

    def _update_kwargs(self, kwargs: Any) -> Dict[str, Any]:
        """
        Update the kwargs.
        """
        temp_kwargs: Dict[str, Any] = super()._update_kwargs(kwargs)
        return temp_kwargs

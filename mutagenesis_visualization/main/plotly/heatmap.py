"""
This module contains the plotly heatmap plot.
"""
from pathlib import Path
from typing import Any, Dict, Union
from numpy import around
from pandas.core.frame import DataFrame
from plotly import io as pio
from plotly import graph_objects as go

from mutagenesis_visualization.main.classes.base_model_plotly import Plotly
from mutagenesis_visualization.main.utils.snv import add_snv_boolean
from mutagenesis_visualization.main.utils.pandas_functions import df_rearrange
from mutagenesis_visualization.main.utils.plotly_utils import matplotlib_to_plotly


class HeatmapP(Plotly):
    """
    This class uses plotly to generate a heatmap.
    """
    def __call__(
        self,
        mask_selfsubstitutions: bool = False,
        replicate: int = -1,
        output_html: Union[None, str, Path] = None,
        **kwargs: Any,
    ) -> None:
        """
        Generate a plotly histogram plot.

        Parameters
        ----------
        mask_selfsubstitutions: bool, default False
            If set to true, will assing a score of 0 to each self-substitution.
            ie (A2A = 0)

        replicate : int, default -1
            Set the replicate to plot. By default, the mean is plotted.
            First replicate start with index 0.
            If there is only one replicate, then leave this parameter
            untouched.

        output_html : str, default None
            If you want to export the generated graph into html, add the
            path and name of the file. Example: 'path/filename.html'.

        **kwargs : other keyword arguments
        """
        temp_kwargs: Dict[str, Any] = self._update_kwargs(kwargs)

        # mask self-substitutions
        self.df_output: DataFrame = self.dataframes.df_stopcodons[replicate].copy()
        if mask_selfsubstitutions:
            self.df_output.loc[self.df_output["Sequence"] == self.df_output["Aminoacid"],
                               "Score_NaN"] = 0

        # sort data by rows in specified order by user
        self.df_output = df_rearrange(
            add_snv_boolean(self.df_output),
            temp_kwargs['neworder_aminoacids'],
            values='Score_NaN',
            show_snv=False
        )

        # get labels for texthover and reindex
        text_hover = self.dataframes.df_stopcodons[replicate].pivot(
            values='Variant',
            index='Aminoacid',
            columns='Position',
        )

        # Create figure
        self.fig = go.Figure(data=go.Heatmap(
            z = around(self.df_output.to_numpy(),2),
            x = list(self.df_output.columns),
            y = list(self.df_output.index),
            zmin=temp_kwargs['colorbar_scale'][0],
            zmax=temp_kwargs['colorbar_scale'][1],
            colorscale=matplotlib_to_plotly(temp_kwargs['colormap']),
            text = text_hover.reindex(list(self.df_output.index)),
            colorbar=dict( # modify colorbar properties
                len=0.8,
                thickness=10,
                outlinewidth=2,
                outlinecolor='rgb(0,0,0)',
                showticklabels=True,
                xpad=0,
            ),
        ))
        self._tune_plot(temp_kwargs)
        self._save_html(output_html, temp_kwargs)

    def _tune_plot(self, temp_kwargs: Dict[str, Any]) -> None:
        """
        Change stylistic parameters of the plot.
        """
        self.fig.update_traces(
            hovertemplate='Aminoacid substitution: %{text}<br>Enrichment score: %{z}<extra></extra>'
        )

        # Style
        pio.templates.default = "plotly_white"

        # UPDATE AXES
        self.fig.update_xaxes(
            title_text=temp_kwargs['x_label'],
            showline=True,
            linewidth=2,
            linecolor='black',
            ticks=None,
            mirror=True,
            side="top",
            dtick=1,  #frequency of ticks
            tickangle=0,
            tickvals=list(self.df_output.columns),
            ticktext=list(self.sequence),
            tickfont=dict(size=8, color='black')
        )
        self.fig.update_yaxes(
            title_text=temp_kwargs['y_label'],
            showline=True,
            linewidth=2,
            linecolor='black',
            ticks=None,
            mirror=True,
            dtick=1,  #frequency of ticks
            autorange="reversed",
            #tickvals = list(df_output.columns),
            #ticktext= ['Healthy', 'Healthy', 'Moderate', 'Diseased', 'Diseased'],
            tickfont=dict(size=8, color='black'),
        )

        # Layout and title parameters https://plotly.com/python/figure-labels/
        self.fig.update_layout(
            width=14 * len(self.df_output.columns) / 165 * 90,
            height=2.65 * 120,
            font=dict(family="Arial, monospace", size=12, color="black"),
            title={'text': temp_kwargs['title'], 'xanchor': 'center', 'yanchor': 'top', 'x': 0.5},
        )

    def _update_kwargs(self, kwargs: Any) -> Dict[str, Any]:
        """
        Update the kwargs.
        """
        temp_kwargs: Dict[str, Any] = super()._update_kwargs(kwargs)
        temp_kwargs['figsize'] = kwargs.get('figsize', (8, 3))
        temp_kwargs['x_label'] = kwargs.get('x_label', '')
        temp_kwargs['y_label'] = kwargs.get('y_label', '')
        return temp_kwargs

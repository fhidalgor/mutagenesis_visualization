"""
This module contains the plotly rank plot.
"""
from pathlib import Path
from typing import Any, Dict, Union
import numpy as np
from plotly import io as pio
from plotly import graph_objects as go
from mutagenesis_visualization.main.classes.base_model_plotly import Plotly


class RankP(Plotly):
    """
    This class uses plotly to generate a rank plot.
    """
    def __call__(
        self,
        mode: str = 'pointmutant',
        replicate: int = -1,
        output_html: Union[None, str, Path] = None,
        **kwargs: Any,
    ) -> None:
        """
        Generate a plotly rank plot so every mutation/residue is sorted based
        on enrichment score.

        Parameters
        ----------
        mode : str, default 'pointmutant'.
            Alternative set to "mean" for the mean of each position.

        replicate : int, default -1
            Set the replicate to plot. By default, the mean is plotted.
            First replicate start with index 0.
            If there is only one replicate, then leave this parameter
            untouched.

        output_html : str, default None
            If you want to export the generated graph into html, add the path and name of the file.
            Example: 'path/filename.html'.

        **kwargs : other keyword arguments
        """
        temp_kwargs: Dict[str, Any] = self._update_kwargs(kwargs)

        # Sort by enrichment scores
        self.df_output = self.dataframes.df_notstopcodons[replicate].sort_values(by=['Score']
                                                                                 ).copy()

        # Chose mode:
        if mode == 'mean':
            self.df_output = self.df_output.groupby(by=['Position'], as_index=False).mean()
            self.df_output.sort_values(by=['Score'], inplace=True)
            self.df_output['Variant'] = self.df_output['Position']

        self.fig = go.Figure(
            data=[
                go.Scatter(
                    x=np.arange(len(self.df_output), 0, -1),
                    y=self.df_output['Score'],
                    text=self.df_output['Variant'],
                    marker_color=temp_kwargs['color'],
                    line=None,
                    mode="markers"
                )
            ]
        )

        self._save_html(output_html, temp_kwargs)

    def _tune_plot(self, temp_kwargs: Dict[str, Any]) -> None:
        """
        Change stylistic parameters of the plot.
        """
        # Style
        pio.templates.default = "plotly_white"

        # Axes https://plotly.com/python/axes/
        self.fig.update_traces(
            mode="markers",
            hovertemplate='Position: %{x}<br>Score: %{y}<br>Variant: %{text}<extra></extra>',
        )
        self.fig.update_xaxes(
            title_text=temp_kwargs['x_label'],
            showline=True,
            linewidth=2,
            linecolor='black',
            ticks="outside",
            mirror=True
        )
        self.fig.update_yaxes(
            title_text=temp_kwargs['y_label'],
            showline=True,
            linewidth=2,
            linecolor='black',
            ticks="outside",
            mirror=True
        )

        # Layout and title parameters https://plotly.com/python/figure-labels/
        self.fig.update_layout(
            width=temp_kwargs['figsize'][0] * 120,
            height=temp_kwargs['figsize'][1] * 120,
            font=dict(family="Arial, monospace", size=12, color="black"),
            title={'text': temp_kwargs['title'], 'xanchor': 'center', 'yanchor': 'top', 'x': 0.5}
        )

    def _update_kwargs(self, kwargs: Any) -> Dict[str, Any]:
        """
        Update the kwargs.
        """
        temp_kwargs: Dict[str, Any] = super()._update_kwargs(kwargs)
        temp_kwargs['figsize'] = kwargs.get('figsize', (4, 3))
        temp_kwargs['x_label'] = kwargs.get('x_label', 'Rank')
        temp_kwargs['y_label'] = kwargs.get('y_label', r'$∆E^i_x$')
        return temp_kwargs

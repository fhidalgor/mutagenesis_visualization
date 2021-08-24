"""
This module contains the plotly scatter plot.
"""
from pathlib import Path
from typing import Any, Dict, Optional, Union
import plotly.io as pio
import plotly.graph_objects as go

import plotly.express as px

from mutagenesis_visualization.main.classes.base_model_plotly import Plotly
from mutagenesis_visualization.main.utils.pandas_functions import (
    process_mean_residue, process_by_pointmutant
)


class ScatterP(Plotly):
    """
    This class uses plotly to generate a scatter plot.
    """
    def __call__(
        self,
        screen_object: Any,
        mode: str = 'pointmutant',
        show_results: bool = False,
        output_html: Union[None, str, Path] = None,
        **kwargs: Dict[str, Any],
    ) -> None:
        """
        Generate a scatter plot between object and a second object of the
        same class.

        Parameters
        ----------
        screen_object : object from class *Screen* to do the scatter with

        mode : str, default 'pointmutant'.
            Alternative set to "mean" for the mean of each position.

        show_results : boolean, default False
            If set to true, will export the details of the linear fit.

        output_html : str, default None
            If you want to export the generated graph into html, add
            the path and name of the file. Example: 'path/filename.html'.

        **kwargs : other keyword arguments
        """
        temp_kwargs: Dict[str, Any] = self._update_kwargs(kwargs)

        # Chose mode:
        if mode == 'pointmutant':
            self.df_output = process_by_pointmutant(self.dataframe, screen_object.dataframe)
        elif mode == 'mean':
            self.df_output = process_mean_residue(self.dataframe, screen_object.dataframe)
            self.df_output['Variant'] = self.df_output['Position']
        # raise error if mode is not "mean" or "pointmutant"

        # create figure
        self.fig = px.scatter(
            x=self.df_output['dataset_1'],
            y=self.df_output['dataset_2'],
            trendline="ols",
            trendline_color_override="red",
        )

        self._tune_plot(temp_kwargs)
        self._save_html(output_html, temp_kwargs)

        if show_results:
            px.get_trendline_results(self.fig).px_fit_results.iloc[0].summary()

    def _tune_plot(self, temp_kwargs: Dict[str, Any]) -> None:
        """
        Change stylistic parameters of the plot.
        """
        # Style
        pio.templates.default = "plotly_white"

        # Titles
        # hide text labels
        self.fig.update_traces(
            hovertext=self.df_output['Variant'],
            hovertemplate='Score_x: %{x}<br>Score_y: %{y}<br>Variant: %{hovertext}<extra></extra>',
        )
        self.fig.update_xaxes(
            title_text=temp_kwargs['x_label'],
            showline=True,
            linewidth=2,
            linecolor='black',
            ticks="outside",
            mirror=True,
        )
        self.fig.update_yaxes(
            title_text=temp_kwargs['y_label'],
            showline=True,
            linewidth=2,
            linecolor='black',
            ticks="outside",
            mirror=True,
        )

        # Layout and title parameters https://plotly.com/python/figure-labels/
        self.fig.update_layout(
            width=temp_kwargs['figsize'][0] * 120,
            height=temp_kwargs['figsize'][1] * 120,
            font=dict(family="Arial, monospace", size=12, color="black"),
            title={'text': temp_kwargs['title'], 'xanchor': 'center', 'yanchor': 'top', 'x': 0.5}
        )

    def _update_kwargs(self, kwargs: Dict[str, Any]) -> Dict[str, Any]:
        """
        Update the kwargs.
        """
        temp_kwargs: Dict[str, Any] = super()._update_kwargs(kwargs)
        temp_kwargs['figsize'] = kwargs.get('figsize', (4, 3))
        return temp_kwargs

"""
This module contains the plotly mean enrichment plot.
"""
from pathlib import Path
from typing import Any, Dict, Literal, Union
from pandas.core.frame import DataFrame
from plotly import io as pio
from plotly import graph_objects as go

from mutagenesis_visualization.main.classes.base_model_plotly import Plotly
from mutagenesis_visualization.main.utils.pandas_functions import process_rmse_residue


class DifferentialP(Plotly):
    """
    This class uses plotly to generate a differential plot.
    """
    def __call__(
        self,
        screen_object: Any,
        metric: Literal["rmse", "squared", "mean"] = 'rmse',
        plot_type: str = 'bar',
        mode: str = 'mean',
        replicate: int = -1,
        output_html: Union[None, str, Path] = None,
        **kwargs: Any,
    ) -> None:
        """
        Generate a plotly mean plot.

        Parameters
        ----------
        screen_object : another *Screen* object to compare with.

        metric: str, default 'rmse'
            The way to compare the two objects.
            Options are 'rmse' ((x-y)**2/N)**0.5, 'squared' ((x**2-y**2)/N and
            'mean' (x-y)/N.

        plot_type: str, default 'bar'
            Options are 'bar' and 'line'.

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
        temp_kwargs: Dict[str, Any] = self._update_kwargs_2(kwargs, metric)

        self.df_output: DataFrame = process_rmse_residue(
            self.dataframes.df_notstopcodons[replicate],
            screen_object.dataframes.df_notstopcodons[replicate], metric
        )

        # Create figure
        #self.fig = px.bar(data_frame=self.df_output, x='Position', y='Score', color='Color')
        #px.bar was switching colors when the first value of Score was negative

        if plot_type.lower() == 'bar':
            self.fig = go.Figure(
                data=[
                    go.Bar(
                        x=self.df_output['Position'],
                        y=self.df_output['d1 - d2'],
                        marker_color=temp_kwargs['color'],
                        marker_line_width=0,
                    )
                ]
            )
        elif plot_type.lower() == 'line':
            self.fig = go.Figure(
                data=[
                    go.Scatter(
                        x=self.df_output['Position'],
                        y=self.df_output['d1 - d2'],
                        line=dict(color=temp_kwargs['color'], width=2)
                    )
                ]
            )

        self._tune_plot(temp_kwargs)
        self._save_html(output_html, temp_kwargs)

    def _tune_plot(self, temp_kwargs: Dict[str, Any]) -> None:
        """
        Change stylistic parameters of the plot.
        """

        # Color and width of bars
        #self.fig.update_traces (marker_line_width=0, )
        self.fig.update_traces(hovertemplate='Position: %{x}<br>Difference: %{y}<extra></extra>')
        # Layout and title parameters https://plotly.com/python/figure-labels/
        self.fig.update_layout(
            width=temp_kwargs['figsize'][0] * 120,
            height=temp_kwargs['figsize'][1] * 120,
            showlegend=False,
            font=dict(family="Arial, monospace", size=12, color="black"),
            title={
                'text': temp_kwargs['title'],
                'xanchor': 'center',
                'yanchor': 'top',
                'x': 0.5,
            }
        )

        # Style
        pio.templates.default = "plotly_white"

        # UPDATE AXES
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

    def _update_kwargs_2(self, kwargs: Any, metric: str) -> Dict[str, Any]:
        """
        Update the kwargs.
        """
        temp_kwargs: Dict[str, Any] = super()._update_kwargs(kwargs)
        temp_kwargs['tick_spacing'] = kwargs.get('tick_spacing', 5)
        temp_kwargs['color'] = kwargs.get('color', 'black')
        temp_kwargs['figsize'] = kwargs.get('figsize', (4.5, 3))

        if metric == 'mean':
            temp_kwargs['yscale'] = kwargs.get('yscale', (-1, 1))
            temp_kwargs['y_label'] = kwargs.get('y_label', r'Mean Differential $∆E^i_x$')
        elif metric == 'rmse':
            temp_kwargs['yscale'] = kwargs.get('yscale', (0, 2))
            temp_kwargs['y_label'] = kwargs.get('y_label', r'RMSE Differential')
        if metric == 'squared':
            temp_kwargs['yscale'] = kwargs.get('yscale', (0, 2))
            temp_kwargs['y_label'] = kwargs.get('y_label', r'Squared Differential')

        temp_kwargs['x_label'] = kwargs.get('x_label', 'Position')
        return temp_kwargs

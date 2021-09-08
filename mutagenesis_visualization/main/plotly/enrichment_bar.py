"""
This module contains the plotly mean enrichment plot.
"""
from pathlib import Path
from typing import Any, Dict, Union
from plotly import io as pio
from plotly import graph_objects as go

from mutagenesis_visualization.main.classes.base_model_plotly import Plotly
from mutagenesis_visualization.main.utils.plotly_utils import select_grouping
from mutagenesis_visualization.main.utils.pandas_functions import color_data


class EnrichmentBarP(Plotly):
    """
    This class uses plotly to generate a mean enrichment plot.
    """
    def __call__(
        self,
        mode: str = 'mean',
        replicate: int = -1,
        output_html: Union[None, str, Path] = None,
        **kwargs: Any,
    ) -> None:
        """
        Generate a plotly mean plot.

        Parameters
        ----------
        mode : str, default 'mean'
            Specify what enrichment scores to show. If mode = 'mean', it
            will show the mean of each position. If mode = 'A', it will
            show the alanine substitution profile. Can be used for each
            amino acid. Use the one-letter code and upper case.

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
        temp_kwargs['title'] = kwargs.get('title', 'Filtered by: {}'.format(mode))

        # Chose mode:
        self.df_output = select_grouping(self.dataframes.df_notstopcodons[replicate], mode)

        # Calculate colors
        self.df_output['Color'] = self.df_output.apply(
            color_data, axis=1, args=(
                temp_kwargs['color_gof'],
                temp_kwargs['color_lof'],
            )
        )

        # Create figure
        #self.fig = px.bar(data_frame=self.df_output, x='Position', y='Score', color='Color')
        #px.bar was switching colors when the first value of Score was negative

        self.fig = go.Figure(
            data=[
                go.Bar(
                    x=self.df_output['Position'],
                    y=self.df_output['Score'],
                    marker_color=self.df_output['Color'],
                    marker_line_width=0,
                )
            ]
        )

        self._tune_plot(temp_kwargs)
        self._save_html(output_html, temp_kwargs)

    def _tune_plot(self, temp_kwargs: Dict[str, Any]) -> None:
        """
        Change stylistic parameters of the plot.
        """

        self.fig.update_traces(hovertemplate='Position: %{x}<br>Score: %{y}<extra></extra>')

        # Color and width of bars
        #self.fig.update_traces(marker_line_width=0, )

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

    def _update_kwargs(self, kwargs: Any) -> Dict[str, Any]:
        """
        Update the kwargs.
        """
        temp_kwargs: Dict[str, Any] = super()._update_kwargs(kwargs)
        temp_kwargs['figsize'] = kwargs.get('figsize', (4.5, 3))
        temp_kwargs['x_label'] = kwargs.get('x_label', 'Position')
        temp_kwargs['y_label'] = kwargs.get('y_label', 'Enrichment score')
        temp_kwargs['color_gof'] = kwargs.get('color_gof', '#FD3216')
        temp_kwargs['color_lof'] = kwargs.get('color_lof', '#6A76FC')

        return temp_kwargs

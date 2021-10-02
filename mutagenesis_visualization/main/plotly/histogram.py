"""
This module contains the plotly histogram plot.
"""
from pathlib import Path
from typing import Any, Dict, Union
from plotly import io as pio
from plotly import graph_objects as go
from mutagenesis_visualization.main.classes.base_model_plotly import Plotly
from mutagenesis_visualization.main.utils.plotly_utils import select_grouping


class HistogramP(Plotly):
    """
    This class uses plotly to generate a histogram plot.
    """
    def __call__(
        self,
        mode: str = 'pointmutant',
        replicate: int = -1,
        output_html: Union[None, str, Path] = None,
        **kwargs: Any,
    ) -> None:
        """
        Generate a plotly histogram plot.

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

        # Create figure
        self.fig = go.Figure()
        self.fig.add_trace(
            go.Histogram(
                x=self.dataframes.df_notstopcodons[replicate]['Score'],
                histnorm='probability density',
                name='Population',
                marker_color='black',
            )
        )

        # Create second histogram
        if mode != 'pointmutant':
            df_output = select_grouping(self.dataframes.df_notstopcodons[replicate], mode)
            # Add second trace
            self.fig.add_trace(
                go.Histogram(
                    x=df_output['Score'],
                    histnorm='probability density',
                    name=mode.capitalize(),
                    marker_color='green',
                ),
            )
            # Overlay both histograms
            self.fig.update_layout(
                barmode='overlay',
                legend=dict(
                    orientation="h",
                    y=(1 + (1.5 / temp_kwargs['figsize'][1]**2)),
                    bgcolor='rgba(0,0,0,0)',
                ),
            )
        # Title
        temp_kwargs['title'] = 'Histogram filtered by: {}'.format(mode)

        self._tune_plot(temp_kwargs)
        self._save_html(output_html, temp_kwargs)

    def _tune_plot(self, temp_kwargs: Dict[str, Any]) -> None:
        """
        Change stylistic parameters of the plot.
        """
        # Reduce opacity to see both histograms
        self.fig.update_traces(opacity=0.75)

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
            automargin=True,
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

    def _update_kwargs(self, kwargs: Any) -> Dict[str, Any]:
        """
        Update the kwargs.
        """
        temp_kwargs: Dict[str, Any] = super()._update_kwargs(kwargs)
        temp_kwargs['figsize'] = kwargs.get('figsize', (4, 3))
        temp_kwargs['x_label'] = kwargs.get('x_label', 'Enrichment score')
        temp_kwargs['y_label'] = kwargs.get('y_label', 'Probability density')
        temp_kwargs['title'] = kwargs.get('title', 'Histogram')
        return temp_kwargs

"""
This module will have the wrapper function to plot scatters for all
replicates within a Screen object.
"""
from typing import Any, Union, Literal
from pathlib import Path
from itertools import combinations
from copy import deepcopy
from mutagenesis_visualization.main.classes.base_model import Pyplot
from mutagenesis_visualization.main.scatter.scatter import Scatter


class ScatterReplicates(Pyplot):
    """
    Class to generate scatter plots of each pairwise replicate combination.
    """
    def __call__(
        self,
        wt_counts_only: bool = False,
        mode: Literal["mean", "pointmutant"] = 'pointmutant',
        output_file: Union[None, str, Path] = None,
        **kwargs: Any
    ) -> None:
        """
        Generate a series of scatter plots between replicates.

        Parameters
        ----------
        wt_counts_only: bool, optional default False
            If set to true, it will plot the kernel distribution of the
            wild type alleles only. mode will be pointmutant by default.

        mode : str, default 'pointmutant'.
            Alternative set to "mean" for the mean of each position.

        output_file : str, default None
            If you want to export the generated graph, add the path and name
            of the file. Example: 'path/filename.png' or 'path/filename.svg'.

        **kwargs : other keyword arguments
        """
        copy_dataframes = deepcopy(self.dataframes)
        if wt_counts_only:
            copy_dataframes.df_notstopcodons = copy_dataframes.df_wildtype_scores
            mode = "pointmutant"

        for combination in list(combinations(range(0, len(copy_dataframes.df_notstopcodons[:-1])),
                                             2)):
            scatter_obj_1: Scatter = Scatter(dataframes=copy_dataframes, aminoacids=self.aminoacids)
            scatter_obj_2: Scatter = Scatter(dataframes=copy_dataframes, aminoacids=self.aminoacids)

            temp_kwargs = self._update_kwargs(kwargs)
            temp_kwargs["x_label"] = "Replicate " + str(combination[0] + 1)
            temp_kwargs["y_label"] = "Replicate " + str(combination[1] + 1)
            output_file_edited = None
            if output_file:
                output_file = Path(output_file)
                output_file_edited = output_file.with_name(
                    output_file.stem + "_" + str(combination[0] + 1) + "_vs_" +
                    str(combination[1] + 1) + output_file.suffix
                )
            scatter_obj_1(
                scatter_obj_2,
                mode,
                combination[0],
                combination[1],
                output_file_edited,
                **temp_kwargs,
            )

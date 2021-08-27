"""
This module will have the wrapper function to plot scatters for all
replicates within a Screen object.
"""
from typing import Any, Union
from pathlib import Path
from itertools import combinations
from mutagenesis_visualization.main.classes.base_model import Pyplot
from mutagenesis_visualization.main.scatter.scatter import Scatter


class ScatterReplicates(Pyplot):
    """
    Class to generate scatter plots of each pairwise replicate combination.
    """
    def __call__(
        self,
        show_wild_type_counts_only: bool = False,
        mode: str = 'pointmutant',
        output_file: Union[None, str, Path] = None,
        **kwargs: Any
    ) -> None:
        """
        Generate a series of scatter plots between replicates.

        Parameters
        ----------
        show_wild_type_counts_only: bool, optional default False
            If set to true, it will plot the kernel distribution of the
            wild type alleles only. mode will be pointmutant by default.

        mode : str, default 'pointmutant'.
            Alternative set to "mean" for the mean of each position.

        output_file : str, default None
            If you want to export the generated graph, add the path and name
            of the file. Example: 'path/filename.png' or 'path/filename.svg'.

        **kwargs : other keyword arguments
        """
        if show_wild_type_counts_only:
            data_to_use = self.wildtype_scores_replicates
            mode = "pointmutant"
        else:
            data_to_use = self.replicate_dataframes

        for combination in list(combinations(range(0, len(self.replicate_dataframes)), 2)):
            scatter_obj_1: Scatter = Scatter(
                dataframe=data_to_use[combination[0]], aminoacids=self.aminoacids
            )
            scatter_obj_2: Scatter = Scatter(
                dataframe=data_to_use[combination[1]], aminoacids=self.aminoacids
            )

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
                output_file_edited,
                **temp_kwargs,
            )
"""
This module contains the Pymol wrapper class to plot in pymol.
"""
from os import path
from pathlib import Path
from typing import List, Dict, Any, Union
from copy import deepcopy
import logging
from pandas.core.frame import DataFrame
from mutagenesis_visualization.main.classes.base_model import Pyplot
from mutagenesis_visualization.main.utils.kwargs import generate_default_kwargs
from mutagenesis_visualization.main.utils.pymol_utils import (
    light_parameters,
    pymol_fitness,
)

try:
    from ipymol import viewer as pymol
except ModuleNotFoundError:
    logging.warning("iPymol module not loaded properly.")


class Pymol(Pyplot):
    """
    This class acts as a wrapper with the ipymol github repo.
    """
    def __call__(
        self,
        pdb: Union[str, Path],
        mode: str = 'mean',
        residues: List[str] = None,
        position_correction: int = 0,
        esthetic_parameters: bool = True,
        replicate: int = -1,
        **kwargs: Any
    ) -> None:
        """
        Color pymol structure residues. User can specify the residues to
        color, or can use the mutagenesis data. Activating mutations will
        be colored red and loss of function blue. Neutral mutations in
        green. Only works if pymol is your $PATH as pymol or you can
        start PyMOL in server mode. Uses the ipymol package, which needs
        to be installed from Github $pip install
        git+https://github.com/cxhernandez/ipymol, not from pypi (not
        updated here).

        Please ensure that PyMOL is in your $PATH as pymol.


        Parameters
        ----------
        pdb : str
            User should specify the PDB chain in the following format 4G0N_A.
            If you have internet connection, Pymol will download the pdb.
            Otherwise, include the path were your PDB is stored locally.

        mode : str, default 'mean'
            Others: 'snv' 'nonsnv', 'aminoacid'
            Specify what enrichment scores to use. If mode = 'mean', it will
            use the mean of each position to classify the residues. If
            mode = 'A', it will use the Alanine substitution profile. Can be
            used for each amino acid. Use the one-letter code and upper case.

        residues : list , optional
            If user decides to pass custom arguments, use the following format
            residues = ['1,2,3,4-10','12-15,23,24,35','48,49,50,52-60'] which
            are [blue,red,green].

        position_correction : int, default 0
            If the pdb structure has a different numbering of positions than
            you dataset, you can correct for that. If your start_position = 2,
            but in the PDB that same residue is at position 20,
            position_correction needs to be set at 18.
        
        esthetic_parameters : bool, default True
            If set to True, pymol will apply the mutagenesis_visualization
            custom parameters instead of the default Pymol ones.    
    
        replicate : int, default -1
            Set the replicate to plot. By default, the mean is plotted.
            First replicate start with index 0.
            If there is only one replicate, then leave this parameter
            untouched.

        **kwargs : other keyword arguments
            gof : int, default is 1
                cutoff for determining gain of function mutations based on
                mutagenesis data.
            lof : int, default is -1
                cutoff for determining loss of function mutations based on
                mutagenesis data.
            color_gof : str, default 'red'
                Choose color to color positions with an enrichment score > gof.
            color_lof : str, default 'neptunium'
                Choose color to color positions with an enrichment score < lof.


        Returns
        ----------
        Open pymol session with a fetched pdb structure where the residues
        are colored according to the enrichment scores.

        """
        temp_kwargs: Dict[str, Any] = deepcopy(generate_default_kwargs())
        temp_kwargs.update(kwargs)
        temp_kwargs['color_lof'] = kwargs.get('color_lof', 'neptunium')

        # Calculate residues only if they are not given by the user
        if residues is None:
            residues = pymol_fitness(
                self.dataframes.df_notstopcodons[replicate].copy(),
                temp_kwargs['gof'],
                temp_kwargs['lof'],
                mode,
                position_correction,
            )

        # Start Pymol
        if not pymol._process_is_running():
            pymol.start()

        # Fetch structure. If pdb contains a "/", it will assume it is stored locally
        if '.pdb' in str(pdb):
            pymol.load(str(pdb))
            pdb = (path.basename(
                str(pdb)
            )).partition('.')[0]  # Extract filename from pdb and then extract pdb code
        else:
            pymol.fetch(str(pdb))

        # Hide everything
        pymol.do('hide everything')

        # Selection names
        blue = str(pdb) + '_blue'
        red = str(pdb) + '_red'
        white = str(pdb) + '_white'

        # Do selections
        pymol.select(blue, 'resi ' + residues[0])
        pymol.select(red, 'resi ' + residues[1])
        pymol.select(white, 'resi ' + residues[2])

        # Representation parameters
        pymol.show_as('cartoon', pdb)
        pymol.set('cartoon_color', 'chlorine', white)
        pymol.set('cartoon_color', temp_kwargs['color_lof'], blue)
        pymol.set('cartoon_color', temp_kwargs['color_gof'], red)
        pymol.bg_color('white')
        pymol.remove('solvent')

        # light parameters
        if esthetic_parameters:
            light_parameters()

        # deselect everything
        pymol.deselect()

    def quit(self) -> None:
        """
        Quit pymol.
        """
        pymol.quit()

    def _select_mode(self, mode: str) -> DataFrame:
        pass

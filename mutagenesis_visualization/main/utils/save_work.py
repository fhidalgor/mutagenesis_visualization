"""

"""
from pathlib import Path
from typing import Union


def save_work(fig, output_file: Union[None, str, Path], temp_kwargs: dict) -> None:
    '''
    Save file function using pathlib.
    '''
    if output_file:
        fig.savefig(
            Path(output_file),
            format=Path(output_file).suffix.strip('.'),
            bbox_inches='tight',
            dpi=temp_kwargs['dpi'],
            transparent=True
        )

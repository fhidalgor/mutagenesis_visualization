"""
This module tests the process data methods and functions.
"""
import logging
from typing import List
from random import randint
from collections import OrderedDict
import numpy.typing as npt
from mutagenesis_visualization.main.process_data.process_data_utils import (initialize_ordered_dict)

log: logging.Logger = logging.getLogger('test_process_data')


def partition_list(array: npt.NDArray, num_partitions: int) -> list:
    """
    Partition array randomly where each partition has at least one item.
    """
    if num_partitions < 2:
        return [f"{array[0]}:{array[-1]}"]
    partition_idxs: List[int] = []
    while len(partition_idxs) < num_partitions - 1:
        num = randint(0, len(array) - 1)
        if num not in partition_idxs:
            tmp_parts = partition_idxs.copy()
            tmp_parts.append(num)
            tmp_parts.sort()
            idx = tmp_parts.index(num)
            if idx != 0:
                if tmp_parts[idx] - tmp_parts[idx - 1] < 1:
                    continue
            if idx != len(tmp_parts) - 1:
                if tmp_parts[idx + 1] - tmp_parts[idx] < 1:
                    continue
            partition_idxs.append(num)
    partition_idxs.sort()
    parts = [f"{array[0]}:{array[partition_idxs[0] - 1]}"]
    for start, end in zip(partition_idxs, partition_idxs[1 :]):
        parts.append(f"{array[start]}:{array[end - 1]}")
    parts.append(f"{array[partition_idxs[-1]]}:{array[-1]}")
    return parts


def test_initialize_ordered_dict() -> None:
    """
    Test if an odered dict is initialized.
    """
    variants = initialize_ordered_dict(['ACC', 'CCT'])
    assert (isinstance(variants, OrderedDict))

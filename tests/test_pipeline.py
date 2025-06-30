import pytest
from pathlib import Path

def test_data_folder_exists():
    path = Path("data/input")
    assert path.exists(), "Input data folder does not exist!"

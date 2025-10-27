import pytest

from unittest.mock import patch
from SSDraw.cli import main


SINGLE_ARGS_DUMMY = [
    "single",
    "-f",
    "mock.fasta",
    "-p",
    "mock.pdb",
    "-n",
    "mock_name",
    "-o",
    "mock.png",
]


def test_ssdraw_single_help_message():
    with pytest.raises(SystemExit) as e:
        main(argv=["single", "-h"])

    assert e.value.code == 0


def test_ssdraw_multi_help_message():
    with pytest.raises(SystemExit) as e:
        main(argv=["multi", "-h"])

    assert e.value.code == 0


def test_ssdraw_single_missing_required_args():
    with pytest.raises(SystemExit) as e:
        main(argv=["single", "-f", "mock.fasta"])

    print(f"code: {e.value.code}")

    assert e.value.code != 0


def test_ssdraw_single_start_greater_than_end():
    invalid_args = SINGLE_ARGS_DUMMY + ["--start", "10", "--end", "5"]

    with pytest.raises(SystemExit) as e:
        with patch("SSDraw.cli.SSDraw"):
            main(argv=invalid_args)

    assert e.value.code != 0


def test_ssdraw_single_success():
    with patch("SSDraw.cli.SSDraw") as mock_ssdraw:
        main(argv=SINGLE_ARGS_DUMMY)

    mock_ssdraw.assert_called_once()


def test_ssdraw_multi_success():
    multi_args = ["multi", "-i", "mock_input.txt", "-o", "mock_output_dir"]

    with patch("SSDraw.cli.run_multiple_pdbs_on_one_msa") as mock_multi:
        main(argv=multi_args)

    mock_multi.assert_called_once()

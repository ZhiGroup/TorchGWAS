import inspect
import pytest
from torchgwas.api import run_linear_gwas
from torchgwas.cli import _build_parser

def test_tsv_rejected_before_reading_input():
    with pytest.raises(ValueError,match='TSV output has been removed'):
        run_linear_gwas('missing-input',None,sumstats_format='tsv')

def test_public_binary_default_and_no_compression_knob():
    signature=inspect.signature(run_linear_gwas)
    assert signature.parameters['sumstats_format'].default=='binary'
    assert 'sumstats_compression' not in signature.parameters

def test_cli_has_no_tsv_or_export_command():
    parser=_build_parser()
    with pytest.raises(SystemExit):parser.parse_args(['export-sumstats'])
    with pytest.raises(SystemExit):parser.parse_args(['linear','--sumstats-format','tsv'])

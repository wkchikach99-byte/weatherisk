"""Tests for weatherisk.cli — argument parsing and exit codes."""

import pytest


class TestCLI:
    def test_help_exits_zero(self):
        from click.testing import CliRunner
        from weatherisk.cli import main

        runner = CliRunner()
        result = runner.invoke(main, ["--help"])
        assert result.exit_code == 0

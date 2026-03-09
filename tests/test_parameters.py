"""Tests for weatherisk.parameters — parameter presets."""

import pytest


class TestPresets:
    def test_stripes_preset_exists(self):
        from weatherisk.parameters import get_preset

        p = get_preset("stripes")
        assert p.resolution > 0
        assert p.df > 0

    def test_bigsmall_preset_exists(self):
        from weatherisk.parameters import get_preset

        p = get_preset("bigsmall")
        assert p.resolution > 0

    def test_rotate_preset_exists(self):
        from weatherisk.parameters import get_preset

        p = get_preset("rotate")
        assert p.resolution > 0

    def test_unknown_raises(self):
        from weatherisk.parameters import get_preset

        with pytest.raises(KeyError):
            get_preset("nonexistent")

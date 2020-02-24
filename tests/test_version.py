try:
    from unittest.mock import patch
except ImportError:
    from mock import patch
import pytest
import daltools


def test_python_version2():
    with patch("daltools.sys") as mock_sys:
        mock_sys.version_info = (3, 7)
        daltools.verify_version()

        mock_sys.version_info = (2, 7)
        with pytest.raises(AssertionError):
            daltools.verify_version()

"""AOONEINT layouts before and after Dalton commit 4abb4ffce (2026-09-11)."""

import struct

import numpy as np
import pytest

from daltools import one


@pytest.fixture(params=[False, True], ids=["old", "mlclinfo"])
def labeled(request):
    return request.param


@pytest.fixture(params=["i", "q"], ids=["int32", "int64"])
def int_fmt(request):
    return request.param


@pytest.fixture(params=[1, 2, 4, 8], ids=lambda n: f"nsym{n}")
def nsym(request):
    return request.param


@pytest.fixture
def aooneint(tmp_path, labeled, int_fmt, nsym):
    """Write actual Fortran records following WRONEL and WRTUND layouts."""
    filename = tmp_path / "AOONEINT"
    naos = (2,) + (1,) * (nsym - 1)
    if nsym > 1:
        naos = naos[:-1] + (0,)  # Include an empty symmetry block.
    coords = ((1.0, 2.0, 3.0), (-4.0, 5.0, -6.0))

    def pack(fmt, values):
        return struct.pack("=" + fmt * len(values), *values)

    with filename.open("wb") as stream:

        def record(data):
            marker = pack("i", (len(data),))
            stream.write(marker + data + marker)

        def label(name):
            record(b"********" + b"20260911" + b"10:22:02" + name.ljust(8).encode())

        if labeled:
            label("MLCLINFO")
        record(b"AOONEINT format compatibility".ljust(400))
        counts = naos + (0,) * (8 - nsym) if labeled else naos
        record(
            pack(int_fmt, (nsym,) + counts)
            + pack("d", (11.75,) if labeled else (11.75, 0.0, 0.0))
        )

        label("ISORDK")
        if labeled:
            record(
                pack(int_fmt, (2,))
                + pack("d", (8.0, 1.0))
                + pack("d", tuple(c for atom in coords for c in atom))
            )
        else:
            record(pack("d", (1.0e20,) * 4))
            # The old COOO(MXCENT,3) layout includes unused centre slots.
            padded = tuple(c for axis in zip(*coords) for c in (*axis, 0.0, 0.0, 0.0))
            record(
                pack("d", (8.0, 1.0, 0.0, 0.0, 0.0))
                + pack(int_fmt, (2,))
                + pack("d", padded)
            )

        label("INTEGRAL")
        ntri = sum(n * (n + 1) // 2 for n in naos)
        for scale, name in enumerate(("OVERLAP", "ONEHAMIL", "KINETINT"), 1):
            label(name)
            # Leave index 2 unwritten to test sparse-buffer reconstruction.
            indices = tuple(i for i in range(1, ntri + 1) if i != 2)
            values = tuple(scale * i / 10 for i in indices)
            record(
                pack("d", values + (0.0,) * (600 - len(values)))
                + pack(int_fmt, indices + (0,) * (600 - len(indices)))
                + pack(int_fmt, (len(indices),))
            )
            record(
                pack("d", (0.0,) * 600)
                + pack(int_fmt, (0,) * 600)
                + pack(int_fmt, (-1,))
            )
    return filename, naos, coords


def test_header(aooneint, int_fmt, nsym):
    filename, naos, _ = aooneint
    assert one.readhead(filename) == {
        "ttitle": "AOONEINT format compatibility".ljust(400),
        "nsym": nsym,
        "naos": naos,
        "potnuc": 11.75,
        "int_fmt": int_fmt,
        "float_fmt": "d",
    }


def test_isordk(aooneint, labeled):
    filename, _, coords = aooneint
    result = one.readisordk(filename)
    ncent = 2 if labeled else 5
    assert result["nucdep"] == 2
    assert result["chrn"] == (8.0, 1.0) + (0.0,) * (ncent - 2)
    actual = np.array(result["cooo"]).reshape(3, ncent)
    np.testing.assert_array_equal(actual[:, :2], np.array(coords).T)
    np.testing.assert_array_equal(actual[:, 2:], 0.0)


@pytest.mark.parametrize(
    "label,scale", [("OVERLAP", 1), ("ONEHAMIL", 2), ("KINETINT", 3)]
)
def test_integrals(aooneint, label, scale):
    filename, naos, _ = aooneint
    result = one.read(label, filename)
    offset = 0
    for nbas, block in zip(naos, result.subblock):
        size = nbas * (nbas + 1) // 2
        expected = [
            scale * i / 10 if i != 2 else 0.0
            for i in range(offset + 1, offset + size + 1)
        ]
        np.testing.assert_allclose(np.asarray(block), expected)
        offset += size


def test_isordk_cli_coordinates(aooneint, monkeypatch):
    filename, _, coords = aooneint
    output = []
    monkeypatch.setattr("sys.argv", ["one", str(filename), "--isordk"])
    monkeypatch.setattr(one, "print", output.append, raising=False)
    one.main()
    np.testing.assert_array_equal(output[2], np.array(coords).T)


def test_missing_scfinp(aooneint):
    filename, _, _ = aooneint
    with pytest.raises(RuntimeError, match="SCFINP section is not present"):
        one.readscfinp(filename)


@pytest.mark.parametrize("size", [0, 32, 72, 96])
def test_invalid_labeled_integer_format(size):
    with pytest.raises(RuntimeError, match="Unknown binary format"):
        one._get_integer_format(bytes(size), labeled=True)

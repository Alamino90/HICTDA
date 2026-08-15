import numpy as np
import pytest

from hictda.topology import contact_to_distance


def test_inverse_distance():
    matrix = np.array([[1.0, 2.0], [2.0, 1.0]])
    d = contact_to_distance(matrix, epsilon=1e-5)
    assert d[0, 0] == 0
    assert d[0, 1] < d[0, 0] + 1e5


def test_non_square_rejected():
    with pytest.raises(ValueError):
        contact_to_distance(np.ones((2, 3)))

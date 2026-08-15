import numpy as np

from hictda.statistics import extract_h1_persistence, h1_statistics


def test_h1_statistics():
    diagrams = [
        np.array([[0.0, 1.0]]),
        np.array([[1.0, 3.0], [2.0, 5.0], [1.0, np.inf]]),
    ]

    p = extract_h1_persistence(diagrams)
    assert np.allclose(p, [2.0, 3.0])

    stats = h1_statistics(diagrams)
    assert stats["n_h1"] == 2
    assert stats["mean_persistence"] == 2.5
    assert stats["median_persistence"] == 2.5

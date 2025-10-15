#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 14 11:59:00 2025

@author: lion
"""

import numpy as np
import pytest
import random

# ===========================================
# Simple lattice and geometry helpers
# ===========================================

class Lattice:
    def __init__(self, a1, a2=None):
        if a2 is None:
            # Assume a1 is a 2×2 matrix
            self.A = np.array(a1, dtype=float)
        else:
            self.A = np.column_stack((a1, a2)).astype(float)

    def __eq__(self, other):
        return np.allclose(self.A, other.A)


def lattice_type(lat: Lattice):
    """Determine lattice type based on metric tensor."""
    a1, a2 = lat.A[:, 0], lat.A[:, 1]
    a1_len = np.linalg.norm(a1)
    a2_len = np.linalg.norm(a2)
    angle = np.arccos(np.dot(a1, a2) / (a1_len * a2_len))

    if np.isclose(a1_len, a2_len):
        if np.isclose(angle, np.pi / 2):
            return "Square"
        elif np.isclose(angle, np.pi / 3) or np.isclose(angle, 2*np.pi/3):
            return "Hexagonal"
        else:
            return "Oblique"
    else:
        if np.isclose(angle, np.pi / 2):
            return "Rectangular"
        else:
            return "Oblique"


def reciprocal_lattice_vectors(lat: Lattice):
    """Compute reciprocal lattice vectors (2π convention)."""
    A = lat.A
    area = np.cross(A[:, 0], A[:, 1])
    b1 = 2*np.pi * np.array([A[1,1], -A[0,1]]) / area
    b2 = 2*np.pi * np.array([-A[1,0], A[0,0]]) / area
    return np.column_stack((b1, b2))


def get_bz(lat: Lattice):
    """Return the first Brillouin zone as a hexagonal-like polygon (approx)."""
    rlv = reciprocal_lattice_vectors(lat)
    # approximate by Wigner-Seitz cell around origin
    b1, b2 = rlv[:, 0], rlv[:, 1]
    corners = [
        0.5*( b1 + b2),  0.5*( b1 - b2),
        0.5*(-b1 + b2),  0.5*(-b1 - b2),
        0.5*( b1),       0.5*(-b1)
    ]
    return np.array(corners)


def get_bounding_box(polygon):
    xs, ys = polygon[:,0], polygon[:,1]
    return (np.min(xs), np.max(xs)), (np.min(ys), np.max(ys))


def in_polygon(point, polygon):
    """Ray casting algorithm for point-in-polygon."""
    x, y = point
    n = len(polygon)
    inside = False
    px, py = polygon[0]
    for i in range(1, n+1):
        sx, sy = polygon[i % n]
        if ((sy > y) != (py > y)) and (x < (px - sx) * (y - sy) / (py - sy + 1e-12) + sx):
            inside = not inside
        px, py = sx, sy
    return inside


def map_to_bz(k, bz, rlv):
    """Wrap point into Brillouin zone by reciprocal lattice translation."""
    k = np.array(k, dtype=float)
    while not in_polygon(k, bz):
        for j in range(2):
            k -= np.dot(k, rlv[:, j]) / np.dot(rlv[:, j], rlv[:, j]) * rlv[:, j]
    return k


# ===========================================
# Tests
# ===========================================

def test_lattice_constructors():
    a1 = np.array([1.0, 1.0])
    a2 = np.array([1.0, -2.0])
    A  = np.array([[1.0, 1.0],
                   [1.0, -2.0]])
    assert Lattice(a1, a2) == Lattice(A)


def test_lattice_types():
    obl_lat = Lattice([[1.0, 0.0], [0.2, -0.5]])
    assert lattice_type(obl_lat) == "Oblique"

    sqr_lat = Lattice([[1.0, 0.0], [0.0, 1.0]])
    assert lattice_type(sqr_lat) == "Square"

    rec_lat = Lattice([[1.0, 0.0], [0.0, 2.0]])
    assert lattice_type(rec_lat) == "Rectangular"

    hex_lat = Lattice([[1.0, -0.5], [0.0, np.sqrt(3)/2.0]])
    assert lattice_type(hex_lat) == "Hexagonal"


def test_brillouin_zone_mapping():
    hex_lat = Lattice([[1.0, -0.5], [0.0, np.sqrt(3)/2.0]])
    bz = get_bz(hex_lat)
    rlv = reciprocal_lattice_vectors(hex_lat)
    x_range, y_range = get_bounding_box(bz)

    x_range = (2*x_range[0], 2*x_range[1])
    y_range = (2*y_range[0], 2*y_range[1])

    for _ in range(10000):
        kx = x_range[0] + (x_range[1] - x_range[0]) * random.random()
        ky = y_range[0] + (y_range[1] - y_range[0]) * random.random()
        k = map_to_bz([kx, ky], bz, rlv)
        assert in_polygon(k, bz)

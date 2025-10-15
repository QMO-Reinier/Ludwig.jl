#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 14 12:00:07 2025

@author: lion
"""

import pytest
import random

#%%

import lattice_tests

def test_lattice_tests():
    """Run all lattice-related tests."""
    # pytest will discover and run tests inside the module
    pytest.main(["-q", "lattice_tests.py"])

#%%

import test_group_tests

def test_group_tests():
    """Run all group-related tests."""
    pytest.main(["-q", "test_group_tests.py"])

#%%
if __name__ == "__main__":
    # Run both test sets
    print("Running Lattice Tests...")
    pytest.main(["-q", "test_lattice_tests.py"])

    print("\nRunning Group Tests...")
    pytest.main(["-q", "test_group_tests.py"])

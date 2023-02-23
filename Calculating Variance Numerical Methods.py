# -*- coding: utf-8 -*-
"""
Created on Thu Feb 16 13:10:36 2023

@author: tanzi
"""


def variance_two_pass(x):
    n = len(x)
    if n < 2:
        return 0.0
    mean = sum(x) / n
    variance = sum((xi - mean) ** 2 for xi in x) / (n - 1)
    return variance

def variance_one_pass(x):
    n = len(x)
    if n < 2:
        return 0.0
    sum_x = 0.0
    sum_x2 = 0.0
    for xi in x:
        sum_x += xi
        sum_x2 += xi ** 2
    variance = (sum_x2 / n) - ((sum_x / n) ** 2)
    return variance

# Test on the first case
x1 = [i/100 for i in range(0,10,1)]
variance_2pass_1 = variance_two_pass(x1)
variance_1pass_1 = variance_one_pass(x1)
print("Two-pass variance for x1:", variance_2pass_1)
print("One-pass variance for x1:", variance_1pass_1)

# Test on the second case
x2 = [123456789.0 + 0.01 * i for i in range(10)]
variance_2pass_2 = variance_two_pass(x2)
variance_1pass_2 = variance_one_pass(x2)
print("Two-pass variance for x2:", variance_2pass_2)
print("One-pass variance for x2:", variance_1pass_2)
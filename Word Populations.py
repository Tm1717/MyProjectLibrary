# -*- coding: utf-8 -*-
"""
Created on Thu Feb  2 00:10:49 2023

@author: tanzi
"""

def square(x):
  return x**2

def cube(x):
  return x**3

func_dict = {
    "square": square,
    "cube": cube,
    "abs": abs,
    "pow": pow
}

def descr(function_name, x):
  f = func_dict[function_name]
  print("Function: ", function_name)
  print("Argument value: ", x)
  print("Result: ", f(x))

print("Example 1:")
descr("square", 3)

print("\nExample 2:")
descr("abs", -4)



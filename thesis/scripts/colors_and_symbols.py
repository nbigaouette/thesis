#!/usr/bin/env python2
# -*- coding: utf-8 -*-


if __name__ == "__main__":
    import sys
    print "Don't run this file directly. It is used by other scripts."
    sys.exit(0)

# NOTE: i must start at 0.

colors = ['b', 'r', 'g', 'm', 'c']
symbols = ['-', '--', '-.', ':']
characters = ['x', 'o', '^', 's', 'p', '*', '+', 'D']

def get_i_colors(i):
    return i % len(colors)

def get_i_symbols(i, array):
    return int(i / len(colors)) % len(array)

def color(i):
    return colors[get_i_colors(i)]

def symbol(i):
    return symbols[get_i_symbols(i, symbols)]

def character(i):
    return characters[get_i_symbols(i, characters)]

def symb_col(i):
    return color(i) + symbol(i)


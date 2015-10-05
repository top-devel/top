#!/usr/bin/python

import os
import sys
sys.path.append(os.environ['HOME'] + '/local/bin')

import top

top.load('eq_poly_esther')
top.init('/home/bputigny/Work/Code/top/models/poly_ester/dati')
top.run_arncheb(0.791171838986203e-01)

vecps = top.get_vecps()
valps = top.get_valps()
print("valps:" + str(valps))

top.load('eq_ESTER_all_lagrange')
top.init('/home/bputigny/Work/Code/top/models/ester/dati')
top.run_arncheb(18.484)

vecps = top.get_vecps()
valps = top.get_valps()
print("valps:" + str(valps))

pl	0	#
		erase
		LOCATION 4500 31000 3500 31000
		limits 10 15 33 35
		expand 1.5
		lweight 3
		ticksize -1 0 -1 0
		box
		#
		rd grmonty.spec
		#
		histogram lw ll4 
		#
		xla "\nu [Hz]"
		yla "\nu L_\nu [erg s^{-1}]"
		#
		re speclab.m
		speclab
		#
		limits 10 15 33 35
		#
comp	1	#
		rd $1
		histogram lw ll4
		#
rd	1	#
		da $1
		#
		read {lw 1}
		read {l0 2  ta0 3  ts0 4  x10 5  x20 6  x30 7 nsc0 8}
		read {l1 9  ta1 10 ts1 11 x11 12 x21 13 x31 14 nsc1 15}
		read {l2 16 ta2 17 ts2 18 x12 19 x22 20 x32 21 nsc2 22}
		read {l3 23 ta3 24 ts3 25 x13 26 x23 27 x33 28 nsc3 29}
		read {l4 30 ta4 31 ts4 32 x14 33 x24 34 x34 35 nsc4 36}
		read {l5 37 ta5 38 ts5 39 x15 40 x25 41 x35 42 nsc5 43}
		#
		set small = 1.e-12
		set ll0 = lg(l0+small) + lg(3.83e33)
		set ll1 = lg(l1+small) + lg(3.83e33)
		set ll2 = lg(l2+small) + lg(3.83e33)
		set ll3 = lg(l3+small) + lg(3.83e33)
		set ll4 = lg(l4+small) + lg(3.83e33)
		set ll5 = lg(l5+small) + lg(3.83e33)
		#
		set lw = lw + lg(9.1e-28*3.e10*3.e10/(6.626e-27))
		#

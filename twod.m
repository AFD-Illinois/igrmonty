Jval    2       #
		do i=$1,$2,1 {\
                        Jcal $i
                }
mov	3	#
		pl $1 $2
		do i=$2+1,$3,1 {\
			rd $i
			pla $1
		}
		#
pl	2	#
		#
		rd $2
		setgrd
		pla $1
setgrd	0	#
		set i=1,$n1*$n2
		#
		image($n1,$n2) $startx1 $stopx1 $startx2 $stopx2
		#
		set i1 = int((x1 - $startx1 - 0.4*$dx1)/$dx1)
		set i2 = int((x2 - $startx2 - 0.4*$dx2)/$dx2)
		#
pla	1	#
		#
		set image[i1,i2] = $1[i-1]
		#
		limits $startx1 $stopx1 $startx2 $stopx2
		erase
		minmax min max echo $min $max
		if($min*$max < 0.) {\
	                #
                        define delta ((-$min)/10.)
                        set lev=$min,-$delta,$delta
                        levels lev
                        ltype 2
                        contour
                        #
                        define delta ($max/10.)
                        set lev=$delta,$max,$delta
                        levels lev
                        ltype 0
                        contour
                        #
		} \
		else {\
			set lev=$min,$max,($max-$min)/10.
			levels lev
			ltype 0
			contour
		}
		#
		define tmp (lg(r[0]))
		define dum1 ($startx1/ln(10.) + $tmp)
		define dum2 ($stopx1/ln(10.)  + $tmp)
		limits $dum1 $dum2 $startx2 $stopx2
		ticksize -1 0 0 0 
		box
		limits $startx1 $stopx1 $startx2 $stopx2
		#
rd	1	#
		if($1 < 10) {define num <00$1>} \
                else {if($1 < 100) {define num <0$1>} \
                else {define num <$1>}}
                echo $num
		rdp dump$num
rdp	1	#
		#
		da ../$1
		lines 1 1
		read {_t 1 _n1 2 _n2 3 _startx1 4 _startx2 5 _dx1 6 _dx2 7 _gam 11}
		define gam (_gam)
		define n1 (_n1)
		define n2 (_n2)
		define startx1 (_startx1)
		define startx2 (_startx2)
		define dx1 (_dx1)
		define dx2 (_dx2)
		define stopx1 ($startx1 + $n1*$dx1)
		define stopx2 ($startx2 + $n2*$dx2)
		lines 2 1000000
		#
		#
		read {x1 1 x2 2 r 3 h 4 rho 5 u 6 v1 7 v2 8 v3 9}
		read {B1 10 B2 11 B3 12}
		read {divb 13}
		read {uu0 14 uu1 15 uu2 16 uu3 17}
		read {ud0 18 ud1 19 ud2 20 ud3 21}
		read {bu0 22 bu1 23 bu2 24 bu3 25}
		read {bd0 26 bd1 27 bd2 28 bd3 29}
		read {v1m 30 v1p 31 v2m 32 v2p 33}
		read {gdet 34}
		read {ju0 35}
		read {ju1 36}
		read {ju2 37}
		read {ju3 38}
		read {jd0 39}
		read {jd1 40}
		read {jd2 41}
		read {jd3 42}
                #
                set jsq = ju0*jd0 + ju1*jd1 + ju2*jd2 + ju3*jd3 
                set jdu = ju0*ud0 + ju1*ud1 + ju2*ud2 + ju3*ud3 
                set Jsq = jsq + jdu*jdu 
                set gJsq = gdet*Jsq 
		#
		echo "assuming gam = " $gam
		#
		# auxiliary info.
		set p = ($gam - 1.)*u
		set K = p*rho**(-$gam)
		set lK = lg(K)
		set EF = rho + $gam*u
		set bsq = bu0*bd0 + bu1*bd1 + bu2*bd2 + bu3*bd3
		set EE = bsq + EF
		set va2 = bsq/EE
		set cs2 = $gam*($gam - 1.)*u/EF
		set cms2 = cs2 + va2 - cs2*va2
		set thetae = (p/rho)*918.059
		set lT = lg(thetae)
		#
		set bsq0 = bu0*bd0
		set bsq1 = bu1*bd1
		set bsq2 = bu2*bd2
		set bsq3 = bu3*bd3
		set ptot = p + 0.5*bsq
		#
		set lbsq = lg(bsq + 1.e-20)
		set lrho = lg(rho)
		set lv2 = lg(abs(v2) + 1.e-20)
		set ldivb = lg(abs(divb) + 1.e-20)
		set ibeta = 0.5*bsq/p
		set libeta = lg(abs(ibeta) + 1.e-20)
		set lthetae = lg(thetae)
		#
		set brel = 0.5*bsq/rho
		set lbrel = lg(abs(brel)+1.e-20)
		#
		set rv1 = r*v1 
		#
		set eflem = ( bsq*uu1*ud0 - bu1*bd0 )
		set eflma = ( (rho+p+u)*uu1*ud0 )
		#
		set lflem = ( bsq*uu1*ud3 - bu1*bd3 )
		set lflma = ( (rho+p+u)*uu1*ud3 )
		#
		set J = bsq*thetae**2*rho
		set lJ = lg(J)
                set rhoch = -(ud0*ju0 + ud1*ju1 + ud2*ju2 + ud3*ju3)
		#
		set lr = lg(r)
Jcal	1	#
		if($1 < 10) {define num <00$1>} \
                else {if($1 < 100) {define num <0$1>} \
                else {define num <$1>}}
                echo $num
		#
		da dumps/dump$num
		lines 2 1000000
		#
		#
		read {ud0 18 ud1 19 ud2 20 ud3 21}
		read {gdet 34}
		read {ju0 35 ju1 36 ju2 37 ju3 38}
		read {jd0 39 jd1 40 jd2 41 jd3 42}
                #
                set jsq = ju0*jd0 + ju1*jd1 + ju2*jd2 + ju3*jd3 
                set jdu = ju0*ud0 + ju1*ud1 + ju2*ud2 + ju3*ud3 
                set Jsq = jsq + jdu*jdu 
                set gJsq = gdet*Jsq 
                define Jvald (sum(gJsq))
		#
		echo $1 $Jvald
		#

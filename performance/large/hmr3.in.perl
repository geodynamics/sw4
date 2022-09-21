# SW4 input RUN: Hayward_scenario_M6.5_s500_3Drfile_h25m

# using parallel i/o
fileio path=hayward-att-h270-ref-result pfs=1 nwriters=32 verbose=1

# Grid coords:
# DON'T CHANGE AZ
grid x=100e3 y=100e3 z=30e3 h=270  lat=38.33 lon=-122.075 az=143.6380001671 proj=tmerc datum=NAD83 lon_p=-123.0 lat_p=35.0 scale=0.9996
refinement zmax=950
refinement zmax=5e3

time t=9.0
#time steps=10

globalmaterial  vpmin=866 vsmin=500

attenuation phasefreq=2.0 nmech=3 maxfreq=3.0

# Earth Model

# Rfile USGS 3D model + topography

# LC (quartz) parallel file system
#rfile filename=rfile/USGSBayAreaVM-08.3.0-corder.rfile 
#topography input=rfile zmax=2.5e3 order=3 file=rfile/USGSBayAreaVM-08.3.0-corder.rfile

# Coral runs
#rfile filename=USGSBayAreaVM-08.3.0-corder.rfile directory=/usr/workspace/wsrzd/ramesh/Project6/SW4/sw4/performance/large
#topography input=rfile zmax=2.5e3 order=3 file=/usr/workspace/wsrzd/ramesh/Project6/SW4/sw4/performance/large/USGSBayAreaVM-08.3.0-corder.rfile
rfile filename=USGSBayAreaVM-08.3.0-corder.rfile directory=/global/cfs/cdirs/m3354
topography input=rfile order=3 zmax=2.5e3 file=/global/cfs/cdirs/m3354/USGSBayAreaVM-08.3.0-corder.rfile

#rfile filename=USGSBayAreaVM-08.3.0-corder.rfile directory=/usr/workspace/wsb/ramesh/Project6/SW4/sw4_new_raja/sw4/performance/large
#topography input=rfile zmax=2.5e3 order=3 file=/usr/workspace/wsb/ramesh/Project6/SW4/sw4_new_raja/sw4/performance/large/USGSBayArea
#VM-08.3.0-corder.rfile
#rfile filename=USGSBayAreaVM-08.3.0-corder.rfile
#topography input=rfile zmax=2.5e3 order=3 file=USGSBayAreaVM-08.3.0-corder.rfile

# NEED TO CORRECT THE CORNER FREQ FOR H=100
#prefilter fc2=0.5 type=lowpass passes=2 order=2

# SRF rupture
#rupture file=m6.5-20.0x13.0.s500.v5.1.srf

# 1 source to make SW4 happy
source x=50e3 y=50e3 depth=3e3 freq=2.5 type=Gaussian t0=2.4 m0=1.7e17 strike=142.1 rake=180 dip=95.6 

# Image output z=0
# image z=0 mode=p file=image cycle=0
# image z=0 mode=s file=image cycle=0
# image z=1000 mode=s file=image cycle=0
# image z=2000 mode=s file=image cycle=0
# image z=5000 mode=s file=image cycle=0
# image z=0 mode=rho file=image cycle=0
#image x=50000 mode=ux file=image cycleInterval=100
#image x=50000 mode=uy file=image cycleInterval=100
#image x=50000 mode=uz file=image cycleInterval=100
#image y=50000 mode=ux file=image cycleInterval=100
#image y=50000 mode=uy file=image cycleInterval=100
#image y=50000 mode=uz file=image cycleInterval=100
#image z=0 mode=lat file=image cycle=0
#image z=0 mode=lon file=image cycle=0
#image z=0 mode=topo file=image cycle=0
# save the magnitude of the velocity (dependent variable with SRF rupture)
#image z=0 mode=mag file=image timeInterval=1.0
# image z=0 mode=hmax file=image time=59.2
# image z=0 mode=hmax file=image timeInterval=0.1
# image z=0 mode=vmax file=image time=59.2
# image z=0 mode=grid file=image cycle=0

# record 1 station for verification
rec file=sta1 lon=-122.25 lat=37.85 depth=0 sacformat=0 usgsformat=1

# SF Bay Area stations
#rec sta=BK.BDM lon=-121.8655 lat=37.9540 depth=0 file=BK.BDM nsew=1 sacformat=0 usgsformat=1
#rec sta=BK.BKS lon=-122.2356 lat=37.8762 depth=0 file=BK.BKS nsew=1
#rec sta=BK.BL67 lon=-122.2432 lat=37.8749 depth=0 file=BK.BL67 nsew=1
#rec sta=BK.BL88 lon=-122.2543 lat=37.8772 depth=0 file=BK.BL88 nsew=1
#rec sta=BK.BRIB lon=-122.1518 lat=37.9189 depth=0 file=BK.BRIB nsew=1
#rec sta=BK.BRK lon=-122.2610 lat=37.8735 depth=0 file=BK.BRK nsew=1
#rec sta=BK.CVS lon=-122.4584 lat=38.3453 depth=0 file=BK.CVS nsew=1
#rec sta=BK.FARB lon=-123.0011 lat=37.6978 depth=0 file=BK.FARB nsew=1
#rec sta=BK.JRSC lon=-122.2387 lat=37.4037 depth=0 file=BK.JRSC nsew=1
#rec sta=BK.MCCM lon=-122.8802 lat=38.1448 depth=0 file=BK.MCCM nsew=1
#rec sta=BK.MHC lon=-121.6426 lat=37.3416 depth=0 file=BK.MHC nsew=1
#rec sta=BK.MHDL lon=-122.4943 lat=37.8423 depth=0 file=BK.MHDL nsew=1
#rec sta=BK.OXMT lon=-122.4243 lat=37.4994 depth=0 file=BK.OXMT nsew=1
#rec sta=BK.PETB lon=-122.5011 lat=38.1189 depth=0 file=BK.PETB nsew=1
#rec sta=BK.POTR lon=-121.9353 lat=38.2026 depth=0 file=BK.POTR nsew=1
#rec sta=BK.RFSB lon=-122.3361 lat=37.9161 depth=0 file=BK.RFSB nsew=1
#rec sta=BK.SCCB lon=-121.8642 lat=37.2874 depth=0 file=BK.SCCB nsew=1
#rec sta=BK.STAN lon=-122.1751 lat=37.4039 depth=0 file=BK.STAN nsew=1
#rec sta=BK.VAK lon=-122.2489 lat=37.8775 depth=0 file=BK.VAK nsew=1
#rec sta=BK.VALB lon=-122.2753 lat=38.1215 depth=0 file=BK.VALB nsew=1
#rec sta=BK.WENL lon=-121.7570 lat=37.6221 depth=0 file=BK.WENL nsew=1
#rec sta=CE.57011 lon=-121.0367 lat=37.6659 depth=0 file=CE.57011 nsew=1
#rec sta=CE.57031 lon=-121.9548 lat=37.3925 depth=0 file=CE.57031 nsew=1
#rec sta=CE.57064 lon=-121.9190 lat=37.5300 depth=0 file=CE.57064 nsew=1
#rec sta=CE.57195 lon=-121.3140 lat=37.9750 depth=0 file=CE.57195 nsew=1
#rec sta=CE.57201 lon=-121.6844 lat=37.7009 depth=0 file=CE.57201 nsew=1
#rec sta=CE.57202 lon=-121.4985 lat=37.6384 depth=0 file=CE.57202 nsew=1
#rec sta=CE.57212 lon=-121.5641 lat=37.7124 depth=0 file=CE.57212 nsew=1
#rec sta=CE.57213 lon=-121.9689 lat=37.5529 depth=0 file=CE.57213 nsew=1
#rec sta=CE.57227 lon=-121.7164 lat=37.6811 depth=0 file=CE.57227 nsew=1
#rec sta=CE.57252 lon=-121.7172 lat=37.6911 depth=0 file=CE.57252 nsew=1
#rec sta=CE.57253 lon=-121.7005 lat=37.6869 depth=0 file=CE.57253 nsew=1
#rec sta=CE.57254 lon=-121.7029 lat=37.6724 depth=0 file=CE.57254 nsew=1
#rec sta=CE.57299 lon=-121.9980 lat=37.3380 depth=0 file=CE.57299 nsew=1
#rec sta=CE.57307 lon=-121.9350 lat=37.7440 depth=0 file=CE.57307 nsew=1
#rec sta=CE.57311 lon=-121.9725 lat=37.5513 depth=0 file=CE.57311 nsew=1
#rec sta=CE.57333 lon=-121.9071 lat=37.3555 depth=0 file=CE.57333 nsew=1
#rec sta=CE.57370 lon=-121.7640 lat=37.2903 depth=0 file=CE.57370 nsew=1
#rec sta=CE.57371 lon=-121.8289 lat=37.2730 depth=0 file=CE.57371 nsew=1
#rec sta=CE.57384 lon=-121.9612 lat=37.5088 depth=0 file=CE.57384 nsew=1
#rec sta=CE.57392 lon=-121.9579 lat=37.5570 depth=0 file=CE.57392 nsew=1
#rec sta=CE.57426 lon=-121.8310 lat=37.9722 depth=0 file=CE.57426 nsew=1
#rec sta=CE.57428 lon=-121.9877 lat=37.5004 depth=0 file=CE.57428 nsew=1
#rec sta=CE.57431 lon=-121.9490 lat=37.5240 depth=0 file=CE.57431 nsew=1
#rec sta=CE.57432 lon=-121.9542 lat=37.7017 depth=0 file=CE.57432 nsew=1
#rec sta=CE.57444 lon=-121.9947 lat=37.7712 depth=0 file=CE.57444 nsew=1
#rec sta=CE.57449 lon=-121.8124 lat=37.6957 depth=0 file=CE.57449 nsew=1
#rec sta=CE.57451 lon=-121.8891 lat=37.6986 depth=0 file=CE.57451 nsew=1
#rec sta=CE.57452 lon=-121.8931 lat=37.7008 depth=0 file=CE.57452 nsew=1
#rec sta=CE.57453 lon=-121.9555 lat=37.7682 depth=0 file=CE.57453 nsew=1
#rec sta=CE.57458 lon=-121.4218 lat=37.7665 depth=0 file=CE.57458 nsew=1
#rec sta=CE.57475 lon=-121.9589 lat=37.7558 depth=0 file=CE.57475 nsew=1
#rec sta=CE.57528 lon=-121.7730 lat=37.7530 depth=0 file=CE.57528 nsew=1
#rec sta=CE.57531 lon=-121.6956 lat=37.9332 depth=0 file=CE.57531 nsew=1
#rec sta=CE.57534 lon=-121.6219 lat=37.9119 depth=0 file=CE.57534 nsew=1
#rec sta=CE.57600 lon=-121.9378 lat=37.3290 depth=0 file=CE.57600 nsew=1
#rec sta=CE.57950 lon=-121.9278 lat=37.4682 depth=0 file=CE.57950 nsew=1
#rec sta=CE.58065 lon=-122.0319 lat=37.2552 depth=0 file=CE.58065 nsew=1
#rec sta=CE.58086 lon=-122.1726 lat=37.4354 depth=0 file=CE.58086 nsew=1
#rec sta=CE.58127 lon=-122.2564 lat=37.4286 depth=0 file=CE.58127 nsew=1
#rec sta=CE.58130 lon=-122.4334 lat=37.7401 depth=0 file=CE.58130 nsew=1
#rec sta=CE.58131 lon=-122.4298 lat=37.7894 depth=0 file=CE.58131 nsew=1
#rec sta=CE.58132 lon=-122.5141 lat=37.7785 depth=0 file=CE.58132 nsew=1
#rec sta=CE.58133 lon=-122.4093 lat=37.8026 depth=0 file=CE.58133 nsew=1
#rec sta=CE.58163 lon=-122.3610 lat=37.8070 depth=0 file=CE.58163 nsew=1
#rec sta=CE.58198 lon=-122.4761 lat=37.8072 depth=0 file=CE.58198 nsew=1
#rec sta=CE.58209 lon=-122.1850 lat=37.8785 depth=0 file=CE.58209 nsew=1
#rec sta=CE.58211 lon=-122.1615 lat=37.8479 depth=0 file=CE.58211 nsew=1
#rec sta=CE.58215 lon=-122.3170 lat=37.5840 depth=0 file=CE.58215 nsew=1
#rec sta=CE.58219 lon=-122.0610 lat=37.6569 depth=0 file=CE.58219 nsew=1
#rec sta=CE.58222 lon=-122.4579 lat=37.7924 depth=0 file=CE.58222 nsew=1
#rec sta=CE.58223 lon=-122.3980 lat=37.6220 depth=0 file=CE.58223 nsew=1
#rec sta=CE.58226 lon=-122.2850 lat=37.4560 depth=0 file=CE.58226 nsew=1
#rec sta=CE.58238 lon=-122.4681 lat=37.7054 depth=0 file=CE.58238 nsew=1
#rec sta=CE.58295 lon=-122.2280 lat=37.8480 depth=0 file=CE.58295 nsew=1
#rec sta=CE.58306 lon=-122.4550 lat=37.7630 depth=0 file=CE.58306 nsew=1
#rec sta=CE.58309 lon=-122.3220 lat=37.7980 depth=0 file=CE.58309 nsew=1
#rec sta=CE.58339 lon=-122.0802 lat=37.6971 depth=0 file=CE.58339 nsew=1
#rec sta=CE.58340 lon=-122.0277 lat=37.7152 depth=0 file=CE.58340 nsew=1
#rec sta=CE.58346 lon=-122.1259 lat=37.6799 depth=0 file=CE.58346 nsew=1
#rec sta=CE.58347 lon=-122.1078 lat=37.6747 depth=0 file=CE.58347 nsew=1
#rec sta=CE.58349 lon=-122.1130 lat=37.6950 depth=0 file=CE.58349 nsew=1
#rec sta=CE.58350 lon=-122.0537 lat=37.7029 depth=0 file=CE.58350 nsew=1
#rec sta=CE.58351 lon=-122.1409 lat=37.7095 depth=0 file=CE.58351 nsew=1
#rec sta=CE.58352 lon=-122.1401 lat=37.6883 depth=0 file=CE.58352 nsew=1
#rec sta=CE.58360 lon=-122.0603 lat=37.9036 depth=0 file=CE.58360 nsew=1
#rec sta=CE.58361 lon=-122.1526 lat=37.8952 depth=0 file=CE.58361 nsew=1
#rec sta=CE.58362 lon=-122.1706 lat=37.8995 depth=0 file=CE.58362 nsew=1
#rec sta=CE.58368 lon=-122.2710 lat=37.9806 depth=0 file=CE.58368 nsew=1
#rec sta=CE.58369 lon=-122.0168 lat=37.9147 depth=0 file=CE.58369 nsew=1
#rec sta=CE.58373 lon=-122.3420 lat=37.4650 depth=0 file=CE.58373 nsew=1
#rec sta=CE.58374 lon=-122.2685 lat=37.5182 depth=0 file=CE.58374 nsew=1
#rec sta=CE.58376 lon=-122.1303 lat=37.6233 depth=0 file=CE.58376 nsew=1
#rec sta=CE.58378 lon=-122.3142 lat=37.4841 depth=0 file=CE.58378 nsew=1
#rec sta=CE.58388 lon=-122.0517 lat=37.5663 depth=0 file=CE.58388 nsew=1
#rec sta=CE.58393 lon=-122.0836 lat=37.6566 depth=0 file=CE.58393 nsew=1
#rec sta=CE.58395 lon=-122.3810 lat=37.5925 depth=0 file=CE.58395 nsew=1
#rec sta=CE.58398 lon=-122.2538 lat=37.7620 depth=0 file=CE.58398 nsew=1
#rec sta=CE.58406 lon=-122.0502 lat=37.6260 depth=0 file=CE.58406 nsew=1
#rec sta=CE.58407 lon=-122.1265 lat=37.6656 depth=0 file=CE.58407 nsew=1
#rec sta=CE.58408 lon=-122.3422 lat=37.9608 depth=0 file=CE.58408 nsew=1
#rec sta=CE.58409 lon=-122.0914 lat=37.7095 depth=0 file=CE.58409 nsew=1
#rec sta=CE.58418 lon=-122.0920 lat=37.6252 depth=0 file=CE.58418 nsew=1
#rec sta=CE.58419 lon=-122.0302 lat=37.6146 depth=0 file=CE.58419 nsew=1
#rec sta=CE.58420 lon=-122.0257 lat=37.6540 depth=0 file=CE.58420 nsew=1
#rec sta=CE.58421 lon=-122.3650 lat=37.9480 depth=0 file=CE.58421 nsew=1
#rec sta=CE.58422 lon=-122.3855 lat=37.9262 depth=0 file=CE.58422 nsew=1
#rec sta=CE.58423 lon=-122.1747 lat=37.7808 depth=0 file=CE.58423 nsew=1
#rec sta=CE.58424 lon=-122.1497 lat=37.7451 depth=0 file=CE.58424 nsew=1
#rec sta=CE.58427 lon=-122.1189 lat=37.6430 depth=0 file=CE.58427 nsew=1
#rec sta=CE.58429 lon=-122.2086 lat=37.9057 depth=0 file=CE.58429 nsew=1
#rec sta=CE.58438 lon=-122.2996 lat=37.9299 depth=0 file=CE.58438 nsew=1
#rec sta=CE.58442 lon=-122.1241 lat=37.8563 depth=0 file=CE.58442 nsew=1
#rec sta=CE.58443 lon=-122.0234 lat=37.8518 depth=0 file=CE.58443 nsew=1
#rec sta=CE.58446 lon=-122.3378 lat=37.8075 depth=0 file=CE.58446 nsew=1
#rec sta=CE.58447 lon=-122.0999 lat=37.5251 depth=0 file=CE.58447 nsew=1
#rec sta=CE.58461 lon=-122.1780 lat=37.4243 depth=0 file=CE.58461 nsew=1
#rec sta=CE.58463 lon=-122.2806 lat=37.8097 depth=0 file=CE.58463 nsew=1
#rec sta=CE.58467 lon=-122.0609 lat=37.9727 depth=0 file=CE.58467 nsew=1
#rec sta=CE.58471 lon=-122.2505 lat=37.8761 depth=0 file=CE.58471 nsew=1
#rec sta=CE.58493 lon=-122.0897 lat=37.6976 depth=0 file=CE.58493 nsew=1
#rec sta=CE.58497 lon=-122.0300 lat=37.6750 depth=0 file=CE.58497 nsew=1
#rec sta=CE.58505 lon=-122.3434 lat=37.9355 depth=0 file=CE.58505 nsew=1
#rec sta=CE.58539 lon=-122.3901 lat=37.6750 depth=0 file=CE.58539 nsew=1
#rec sta=CE.58558 lon=-122.3298 lat=37.9783 depth=0 file=CE.58558 nsew=1
#rec sta=CE.58565 lon=-122.5243 lat=37.9694 depth=0 file=CE.58565 nsew=1
#rec sta=CE.58619 lon=-122.2675 lat=37.5302 depth=0 file=CE.58619 nsew=1
#rec sta=CE.58653 lon=-122.0900 lat=37.6700 depth=0 file=CE.58653 nsew=1
#rec sta=CE.58662 lon=-122.3960 lat=37.6390 depth=0 file=CE.58662 nsew=1
#rec sta=CE.58663 lon=-122.0750 lat=37.5362 depth=0 file=CE.58663 nsew=1
#rec sta=CE.58664 lon=-122.1360 lat=37.4920 depth=0 file=CE.58664 nsew=1
#rec sta=CE.58665 lon=-122.0300 lat=37.5349 depth=0 file=CE.58665 nsew=1
#rec sta=CE.58667 lon=-122.0360 lat=37.9730 depth=0 file=CE.58667 nsew=1
#rec sta=CE.58716 lon=-122.4751 lat=37.8060 depth=0 file=CE.58716 nsew=1
#rec sta=CE.58790 lon=-122.2354 lat=37.8206 depth=0 file=CE.58790 nsew=1
#rec sta=CE.58905 lon=-122.5354 lat=37.9480 depth=0 file=CE.58905 nsew=1
#rec sta=CE.67210 lon=-121.9223 lat=38.3994 depth=0 file=CE.67210 nsew=1
#rec sta=CE.67463 lon=-121.5105 lat=38.5854 depth=0 file=CE.67463 nsew=1
#rec sta=CE.67520 lon=-121.9525 lat=38.3507 depth=0 file=CE.67520 nsew=1
#rec sta=CE.67523 lon=-121.7509 lat=38.0183 depth=0 file=CE.67523 nsew=1
#rec sta=CE.67533 lon=-121.6396 lat=38.0155 depth=0 file=CE.67533 nsew=1
#rec sta=CE.67557 lon=-121.8148 lat=38.0141 depth=0 file=CE.67557 nsew=1
#rec sta=CE.67587 lon=-121.8986 lat=38.0390 depth=0 file=CE.67587 nsew=1
#rec sta=CE.67910 lon=-121.7500 lat=38.0300 depth=0 file=CE.67910 nsew=1
#rec sta=CE.67942 lon=-121.5015 lat=38.5736 depth=0 file=CE.67942 nsew=1
#rec sta=CE.68034 lon=-122.2329 lat=38.5381 depth=0 file=CE.68034 nsew=1
#rec sta=CE.68045 lon=-122.0492 lat=38.2612 depth=0 file=CE.68045 nsew=1
#rec sta=CE.68150 lon=-122.2774 lat=38.2704 depth=0 file=CE.68150 nsew=1
#rec sta=CE.68208 lon=-122.7590 lat=38.4970 depth=0 file=CE.68208 nsew=1
#rec sta=CE.68294 lon=-122.2493 lat=38.1246 depth=0 file=CE.68294 nsew=1
#rec sta=CE.68326 lon=-122.8274 lat=38.4005 depth=0 file=CE.68326 nsew=1
#rec sta=CE.68327 lon=-122.6541 lat=38.4674 depth=0 file=CE.68327 nsew=1
#rec sta=CE.68328 lon=-122.7457 lat=38.4694 depth=0 file=CE.68328 nsew=1
#rec sta=CE.68329 lon=-122.7476 lat=38.4413 depth=0 file=CE.68329 nsew=1
#rec sta=CE.68330 lon=-122.6710 lat=38.4382 depth=0 file=CE.68330 nsew=1
#rec sta=CE.68367 lon=-122.2620 lat=38.0043 depth=0 file=CE.68367 nsew=1
#rec sta=CE.68433 lon=-122.5593 lat=38.0981 depth=0 file=CE.68433 nsew=1
#rec sta=CE.68440 lon=-122.8736 lat=38.3504 depth=0 file=CE.68440 nsew=1
#rec sta=CE.68491 lon=-122.7089 lat=38.4375 depth=0 file=CE.68491 nsew=1
#rec sta=CE.68676 lon=-122.7272 lat=38.4705 depth=0 file=CE.68676 nsew=1
#rec sta=CE.69039 lon=-123.0527 lat=38.3110 depth=0 file=CE.69039 nsew=1
#rec sta=CE.69044 lon=-123.2470 lat=38.5176 depth=0 file=CE.69044 nsew=1
#rec sta=GPS.0410 lon=-122.255 lat=38.03023 depth=0 file=GPS.0410 nsew=1
#rec sta=GPS.04KF lon=-122.468 lat=38.24524 depth=0 file=GPS.04KF nsew=1
#rec sta=GPS.04KH lon=-122.253 lat=38.15480 depth=0 file=GPS.04KH nsew=1
#rec sta=GPS.04KJ lon=-122.112 lat=38.08713 depth=0 file=GPS.04KJ nsew=1
#rec sta=GPS.04LF lon=-122.479 lat=38.30255 depth=0 file=GPS.04LF nsew=1
#rec sta=GPS.04LG lon=-122.299 lat=38.27111 depth=0 file=GPS.04LG nsew=1
#rec sta=GPS.6908 lon=-122.114 lat=38.15688 depth=0 file=GPS.6908 nsew=1
#rec sta=GPS.ADOO lon=-122.527 lat=38.23642 depth=0 file=GPS.ADOO nsew=1
#rec sta=GPS.AMER lon=-122.199 lat=38.16860 depth=0 file=GPS.AMER nsew=1
#rec sta=GPS.B468 lon=-122.329 lat=38.34730 depth=0 file=GPS.B468 nsew=1
#rec sta=GPS.CASR lon=-122.747 lat=38.44067 depth=0 file=GPS.CASR nsew=1
#rec sta=GPS.CKET lon=-122.229 lat=38.04275 depth=0 file=GPS.CKET nsew=1
#rec sta=GPS.CORD lon=-122.595 lat=38.18599 depth=0 file=GPS.CORD nsew=1
#rec sta=GPS.DEAL lon=-122.338 lat=38.25781 depth=0 file=GPS.DEAL nsew=1
#rec sta=GPS.ES01 lon=-122.115 lat=38.11372 depth=0 file=GPS.ES01 nsew=1
#rec sta=GPS.GAME lon=-122.175 lat=38.35068 depth=0 file=GPS.GAME nsew=1
#rec sta=GPS.HAGG lon=-122.259 lat=38.32388 depth=0 file=GPS.HAGG nsew=1
#rec sta=GPS.HAZD lon=-122.131 lat=38.10615 depth=0 file=GPS.HAZD nsew=1
#rec sta=GPS.IS01 lon=-122.113 lat=38.11412 depth=0 file=GPS.IS01 nsew=1
#rec sta=GPS.LB05 lon=-122.209 lat=38.47901 depth=0 file=GPS.LB05 nsew=1
#rec sta=GPS.M139 lon=-122.448 lat=38.15136 depth=0 file=GPS.M139 nsew=1
#rec sta=GPS.MADI lon=-122.203 lat=38.31261 depth=0 file=GPS.MADI nsew=1
#rec sta=GPS.MCCM lon=-122.88 lat=38.14479 depth=0 file=GPS.MCCM nsew=1
#rec sta=GPS.MDRN lon=-122.508 lat=38.33987 depth=0 file=GPS.MDRN nsew=1
#rec sta=GPS.OHLN lon=-122.273 lat=38.00625 depth=0 file=GPS.OHLN nsew=1
#rec sta=GPS.P181 lon=-122.377 lat=37.91455 depth=0 file=GPS.P181 nsew=1
#rec sta=GPS.P194 lon=-122.816 lat=38.18572 depth=0 file=GPS.P194 nsew=1
#rec sta=GPS.P196 lon=-122.743 lat=38.29814 depth=0 file=GPS.P196 nsew=1
#rec sta=GPS.P197 lon=-122.767 lat=38.42856 depth=0 file=GPS.P197 nsew=1
#rec sta=GPS.P198 lon=-122.607 lat=38.25987 depth=0 file=GPS.P198 nsew=1
#rec sta=GPS.P199 lon=-122.503 lat=38.26369 depth=0 file=GPS.P199 nsew=1
#rec sta=GPS.P200 lon=-122.452 lat=38.23983 depth=0 file=GPS.P200 nsew=1
#rec sta=GPS.P201 lon=-122.658 lat=38.55981 depth=0 file=GPS.P201 nsew=1
#rec sta=GPS.P202 lon=-122.496 lat=38.42358 depth=0 file=GPS.P202 nsew=1
#rec sta=GPS.P261 lon=-122.218 lat=38.15296 depth=0 file=GPS.P261 nsew=1
#rec sta=GPS.P262 lon=-122.096 lat=38.02515 depth=0 file=GPS.P262 nsew=1
#rec sta=GPS.P263 lon=-122.429 lat=38.57769 depth=0 file=GPS.P263 nsew=1
#rec sta=GPS.P264 lon=-122.195 lat=38.44421 depth=0 file=GPS.P264 nsew=1
#rec sta=GPS.PTRO lon=-121.944 lat=38.20949 depth=0 file=GPS.PTRO nsew=1
#rec sta=GPS.SKAG lon=-122.381 lat=38.15881 depth=0 file=GPS.SKAG nsew=1
#rec sta=GPS.SKYV lon=-122.164 lat=38.10725 depth=0 file=GPS.SKYV nsew=1
#rec sta=GPS.SPRR lon=-122.186 lat=38.20862 depth=0 file=GPS.SPRR nsew=1
#rec sta=GPS.SVIN lon=-122.526 lat=38.03318 depth=0 file=GPS.SVIN nsew=1
#rec sta=GPS.T3RP lon=-122.588 lat=37.92789 depth=0 file=GPS.T3RP nsew=1
#rec sta=GPS.TIBB lon=-122.448 lat=37.89088 depth=0 file=GPS.TIBB nsew=1
#rec sta=GPS.TRAN lon=-122.309 lat=38.32097 depth=0 file=GPS.TRAN nsew=1
#rec sta=GPS.TRCP lon=-122.616 lat=37.99742 depth=0 file=GPS.TRCP nsew=1
#rec sta=GPS.VETH lon=-122.36 lat=38.39683 depth=0 file=GPS.VETH nsew=1
#rec sta=NC.C001 lon=-121.9338 lat=37.5215 depth=0 file=NC.C001 nsew=1
#rec sta=NC.C002 lon=-122.0345 lat=37.5585 depth=0 file=NC.C002 nsew=1
#rec sta=NC.C003 lon=-122.0840 lat=37.7072 depth=0 file=NC.C003 nsew=1
#rec sta=NC.C004 lon=-122.2179 lat=37.8441 depth=0 file=NC.C004 nsew=1
#rec sta=NC.C005 lon=-122.2361 lat=37.8108 depth=0 file=NC.C005 nsew=1
#rec sta=NC.C006 lon=-122.0582 lat=37.5918 depth=0 file=NC.C006 nsew=1
#rec sta=NC.C007 lon=-122.2628 lat=37.8956 depth=0 file=NC.C007 nsew=1
#rec sta=NC.C008 lon=-122.2936 lat=37.8956 depth=0 file=NC.C008 nsew=1
#rec sta=NC.C009 lon=-122.1841 lat=37.8049 depth=0 file=NC.C009 nsew=1
#rec sta=NC.C010 lon=-122.0099 lat=37.9440 depth=0 file=NC.C010 nsew=1
#rec sta=NC.C011 lon=-122.2924 lat=37.9064 depth=0 file=NC.C011 nsew=1
#rec sta=NC.C012 lon=-122.3125 lat=37.9264 depth=0 file=NC.C012 nsew=1
#rec sta=NC.C013 lon=-122.1886 lat=37.8238 depth=0 file=NC.C013 nsew=1
#rec sta=NC.C014 lon=-122.0600 lat=37.6827 depth=0 file=NC.C014 nsew=1
#rec sta=NC.C015 lon=-121.9937 lat=37.5594 depth=0 file=NC.C015 nsew=1
#rec sta=NC.C016 lon=-121.8670 lat=37.3964 depth=0 file=NC.C016 nsew=1
#rec sta=NC.C017 lon=-121.9262 lat=37.5065 depth=0 file=NC.C017 nsew=1
#rec sta=NC.C018 lon=-122.1174 lat=37.9793 depth=0 file=NC.C018 nsew=1
#rec sta=NC.C019 lon=-122.0139 lat=37.9836 depth=0 file=NC.C019 nsew=1
#rec sta=NC.C020 lon=-122.2075 lat=37.8008 depth=0 file=NC.C020 nsew=1
#rec sta=NC.C021 lon=-122.0638 lat=37.6950 depth=0 file=NC.C021 nsew=1
#rec sta=NC.C022 lon=-121.9750 lat=37.5361 depth=0 file=NC.C022 nsew=1
#rec sta=NC.C023 lon=-121.9540 lat=37.5104 depth=0 file=NC.C023 nsew=1
#rec sta=NC.C024 lon=-122.3358 lat=37.9433 depth=0 file=NC.C024 nsew=1
#rec sta=NC.C025 lon=-121.9005 lat=37.6437 depth=0 file=NC.C025 nsew=1
#rec sta=NC.C026 lon=-122.0514 lat=37.8705 depth=0 file=NC.C026 nsew=1
#rec sta=NC.C027 lon=-121.9797 lat=37.7935 depth=0 file=NC.C027 nsew=1
#rec sta=NC.C028 lon=-122.2459 lat=37.8509 depth=0 file=NC.C028 nsew=1 usgsformat=1
#rec sta=NC.C029 lon=-121.9713 lat=37.7888 depth=0 file=NC.C029 nsew=1
#rec sta=NC.C031 lon=-122.2769 lat=37.8632 depth=0 file=NC.C031 nsew=1
#rec sta=NC.C032 lon=-122.1576 lat=38.0825 depth=0 file=NC.C032 nsew=1
#rec sta=NC.C033 lon=-122.2895 lat=37.8467 depth=0 file=NC.C033 nsew=1
#rec sta=NC.C034 lon=-121.9633 lat=37.5728 depth=0 file=NC.C034 nsew=1
#rec sta=NC.C035 lon=-122.0269 lat=37.9698 depth=0 file=NC.C035 nsew=1
#rec sta=NC.C036 lon=-121.7931 lat=37.3802 depth=0 file=NC.C036 nsew=1
#rec sta=NC.C037 lon=-121.8849 lat=37.3551 depth=0 file=NC.C037 nsew=1
#rec sta=NC.C038 lon=-122.1765 lat=37.7205 depth=0 file=NC.C038 nsew=1
#rec sta=NC.C039 lon=-121.9116 lat=37.4671 depth=0 file=NC.C039 nsew=1
#rec sta=NC.C040 lon=-122.3141 lat=37.9897 depth=0 file=NC.C040 nsew=1
#rec sta=NC.C041 lon=-121.8261 lat=37.3757 depth=0 file=NC.C041 nsew=1
#rec sta=NC.C042 lon=-122.1277 lat=37.8223 depth=0 file=NC.C042 nsew=1
#rec sta=NC.C043 lon=-122.2485 lat=37.8685 depth=0 file=NC.C043 nsew=1
#rec sta=NC.C044 lon=-122.2312 lat=37.7308 depth=0 file=NC.C044 nsew=1
#rec sta=NC.C045 lon=-121.9620 lat=37.7496 depth=0 file=NC.C045 nsew=1
#rec sta=NC.C046 lon=-121.8825 lat=37.3536 depth=0 file=NC.C046 nsew=1
#rec sta=NC.C047 lon=-122.0203 lat=37.8291 depth=0 file=NC.C047 nsew=1
#rec sta=NC.C048 lon=-121.9917 lat=37.5751 depth=0 file=NC.C048 nsew=1
#rec sta=NC.C049 lon=-122.2785 lat=37.7766 depth=0 file=NC.C049 nsew=1
#rec sta=NC.C050 lon=-121.9256 lat=37.5240 depth=0 file=NC.C050 nsew=1
#rec sta=NC.C051 lon=-122.2386 lat=37.7467 depth=0 file=NC.C051 nsew=1
#rec sta=NC.C052 lon=-121.8066 lat=37.2962 depth=0 file=NC.C052 nsew=1
#rec sta=NC.C053 lon=-121.8734 lat=37.6854 depth=0 file=NC.C053 nsew=1
#rec sta=NC.C054 lon=-121.8244 lat=37.2712 depth=0 file=NC.C054 nsew=1
#rec sta=NC.C055 lon=-122.1104 lat=37.6574 depth=0 file=NC.C055 nsew=1
#rec sta=NC.C056 lon=-122.1264 lat=37.6363 depth=0 file=NC.C056 nsew=1
#rec sta=NC.C057 lon=-121.6979 lat=37.9555 depth=0 file=NC.C057 nsew=1
#rec sta=NC.C058 lon=-122.2258 lat=37.7564 depth=0 file=NC.C058 nsew=1
#rec sta=NC.C059 lon=-121.7724 lat=37.3426 depth=0 file=NC.C059 nsew=1
#rec sta=NC.C060 lon=-122.0271 lat=37.5848 depth=0 file=NC.C060 nsew=1
#rec sta=NC.C061 lon=-122.2097 lat=37.7779 depth=0 file=NC.C061 nsew=1
#rec sta=NC.C062 lon=-121.9034 lat=37.6675 depth=0 file=NC.C062 nsew=1
#rec sta=NC.C063 lon=-122.1891 lat=37.7651 depth=0 file=NC.C063 nsew=1
#rec sta=NC.C064 lon=-122.1846 lat=37.8051 depth=0 file=NC.C064 nsew=1
#rec sta=NC.C065 lon=-122.2759 lat=37.9629 depth=0 file=NC.C065 nsew=1
#rec sta=NC.C066 lon=-122.2858 lat=37.9132 depth=0 file=NC.C066 nsew=1
#rec sta=NC.C067 lon=-122.2896 lat=37.9761 depth=0 file=NC.C067 nsew=1
#rec sta=NC.C068 lon=-122.0741 lat=37.6446 depth=0 file=NC.C068 nsew=1
#rec sta=NC.C069 lon=-121.8961 lat=37.6478 depth=0 file=NC.C069 nsew=1
#rec sta=NC.C070 lon=-122.1318 lat=37.7028 depth=0 file=NC.C070 nsew=1
#rec sta=NC.C071 lon=-122.1332 lat=37.6821 depth=0 file=NC.C071 nsew=1
#rec sta=NC.C072 lon=-121.7559 lat=37.2858 depth=0 file=NC.C072 nsew=1
#rec sta=NC.C073 lon=-122.0542 lat=37.6558 depth=0 file=NC.C073 nsew=1
#rec sta=NC.CAG lon=-122.4227 lat=37.8573 depth=0 file=NC.CAG nsew=1
#rec sta=NC.CAL lon=-121.8004 lat=37.4503 depth=0 file=NC.CAL nsew=1
#rec sta=NC.CBP lon=-121.7084 lat=37.7477 depth=0 file=NC.CBP nsew=1
#rec sta=NC.CBR lon=-122.0630 lat=37.8165 depth=0 file=NC.CBR nsew=1
#rec sta=NC.CCH1 lon=-122.0966 lat=37.7433 depth=0 file=NC.CCH1 nsew=1
#rec sta=NC.CCO lon=-121.6746 lat=37.2582 depth=0 file=NC.CCO nsew=1
#rec sta=NC.CCOB lon=-121.6731 lat=37.2590 depth=0 file=NC.CCOB nsew=1
#rec sta=NC.CDOB lon=-121.8335 lat=37.7294 depth=0 file=NC.CDOB nsew=1
#rec sta=NC.CDV lon=-121.6793 lat=37.5665 depth=0 file=NC.CDV nsew=1
#rec sta=NC.CGP1 lon=-122.0113 lat=37.6454 depth=0 file=NC.CGP1 nsew=1
#rec sta=NC.CHR lon=-121.7689 lat=37.3729 depth=0 file=NC.CHR nsew=1
#rec sta=NC.CLCB lon=-122.0649 lat=37.7379 depth=0 file=NC.CLCB nsew=1
#rec sta=NC.CMC lon=-122.1768 lat=37.7812 depth=0 file=NC.CMC nsew=1
#rec sta=NC.CMK lon=-121.8665 lat=37.4849 depth=0 file=NC.CMK nsew=1
#rec sta=NC.CMLP lon=-121.8339 lat=37.4489 depth=0 file=NC.CMLP nsew=1
#rec sta=NC.CMOB lon=-121.8027 lat=37.8104 depth=0 file=NC.CMOB nsew=1
#rec sta=NC.CMW lon=-121.8863 lat=37.5406 depth=0 file=NC.CMW nsew=1
#rec sta=NC.CMW1 lon=-121.8874 lat=37.5405 depth=0 file=NC.CMW1 nsew=1
#rec sta=NC.CNI lon=-121.9646 lat=37.6079 depth=0 file=NC.CNI nsew=1
#rec sta=NC.CPI lon=-122.2156 lat=37.9889 depth=0 file=NC.CPI nsew=1
#rec sta=NC.CPM lon=-122.4085 lat=37.9491 depth=0 file=NC.CPM nsew=1
#rec sta=NC.CRH lon=-121.9926 lat=37.8588 depth=0 file=NC.CRH nsew=1
#rec sta=NC.CRPB lon=-121.9072 lat=37.9123 depth=0 file=NC.CRPB nsew=1
#rec sta=NC.CSL lon=-122.1193 lat=37.7242 depth=0 file=NC.CSL nsew=1
#rec sta=NC.CSPB lon=-122.3113 lat=37.9572 depth=0 file=NC.CSPB nsew=1
#rec sta=NC.CSU1 lon=-121.9402 lat=37.6430 depth=0 file=NC.CSU1 nsew=1
#rec sta=NC.CTA lon=-122.0160 lat=38.0269 depth=0 file=NC.CTA nsew=1
#rec sta=NC.CVL lon=-121.8364 lat=37.6236 depth=0 file=NC.CVL nsew=1
#rec sta=NC.CVP lon=-122.2232 lat=37.8840 depth=0 file=NC.CVP nsew=1
#rec sta=NC.CYB lon=-122.3618 lat=37.8079 depth=0 file=NC.CYB nsew=1
#rec sta=NC.CYD lon=-122.0965 lat=37.5628 depth=0 file=NC.CYD nsew=1
#rec sta=NC.CYD1 lon=-122.0967 lat=37.5629 depth=0 file=NC.CYD1 nsew=1
#rec sta=NC.GAXB lon=-122.7573 lat=38.7107 depth=0 file=NC.GAXB nsew=1
#rec sta=NC.J001 lon=-122.3923 lat=37.5927 depth=0 file=NC.J001 nsew=1
#rec sta=NC.J002 lon=-122.1536 lat=37.3677 depth=0 file=NC.J002 nsew=1
#rec sta=NC.J003 lon=-122.1347 lat=37.4364 depth=0 file=NC.J003 nsew=1
#rec sta=NC.J004 lon=-122.0364 lat=37.3175 depth=0 file=NC.J004 nsew=1
#rec sta=NC.J005 lon=-121.9605 lat=37.2961 depth=0 file=NC.J005 nsew=1
#rec sta=NC.J006 lon=-122.1020 lat=37.3569 depth=0 file=NC.J006 nsew=1
#rec sta=NC.J007 lon=-122.2447 lat=37.4924 depth=0 file=NC.J007 nsew=1
#rec sta=NC.J008 lon=-122.1882 lat=37.4759 depth=0 file=NC.J008 nsew=1
#rec sta=NC.J009 lon=-122.5117 lat=37.5197 depth=0 file=NC.J009 nsew=1
#rec sta=NC.J010 lon=-122.4362 lat=37.4580 depth=0 file=NC.J010 nsew=1
#rec sta=NC.J011 lon=-122.2137 lat=37.3695 depth=0 file=NC.J011 nsew=1
#rec sta=NC.J012 lon=-122.1069 lat=37.4060 depth=0 file=NC.J012 nsew=1
#rec sta=NC.J013 lon=-122.4667 lat=37.7560 depth=0 file=NC.J013 nsew=1
#rec sta=NC.J014 lon=-122.1896 lat=37.4165 depth=0 file=NC.J014 nsew=1
#rec sta=NC.J015 lon=-121.8761 lat=37.3319 depth=0 file=NC.J015 nsew=1
#rec sta=NC.J016 lon=-122.2517 lat=37.4574 depth=0 file=NC.J016 nsew=1
#rec sta=NC.J017 lon=-122.2290 lat=37.4612 depth=0 file=NC.J017 nsew=1
#rec sta=NC.J018 lon=-122.0735 lat=37.3196 depth=0 file=NC.J018 nsew=1
#rec sta=NC.J019 lon=-122.1983 lat=37.4343 depth=0 file=NC.J019 nsew=1
#rec sta=NC.J020 lon=-122.4411 lat=37.7660 depth=0 file=NC.J020 nsew=1
#rec sta=NC.J021 lon=-122.2366 lat=37.3565 depth=0 file=NC.J021 nsew=1
#rec sta=NC.J022 lon=-122.4452 lat=37.6300 depth=0 file=NC.J022 nsew=1
#rec sta=NC.J023 lon=-122.3469 lat=37.5255 depth=0 file=NC.J023 nsew=1
#rec sta=NC.J024 lon=-122.0854 lat=37.4021 depth=0 file=NC.J024 nsew=1
#rec sta=NC.J025 lon=-122.0725 lat=37.3948 depth=0 file=NC.J025 nsew=1
#rec sta=NC.J026 lon=-122.0239 lat=37.3494 depth=0 file=NC.J026 nsew=1
#rec sta=NC.J027 lon=-122.0842 lat=37.3476 depth=0 file=NC.J027 nsew=1
#rec sta=NC.J028 lon=-122.2232 lat=37.3880 depth=0 file=NC.J028 nsew=1
#rec sta=NC.J029 lon=-122.4938 lat=37.6341 depth=0 file=NC.J029 nsew=1
#rec sta=NC.J030 lon=-122.1210 lat=37.3862 depth=0 file=NC.J030 nsew=1
#rec sta=NC.J031 lon=-122.4824 lat=37.5813 depth=0 file=NC.J031 nsew=1
#rec sta=NC.J032 lon=-122.4150 lat=37.6315 depth=0 file=NC.J032 nsew=1
#rec sta=NC.J033 lon=-122.1958 lat=37.4529 depth=0 file=NC.J033 nsew=1
#rec sta=NC.J034 lon=-122.1134 lat=37.3702 depth=0 file=NC.J034 nsew=1
#rec sta=NC.J035 lon=-122.3780 lat=37.5578 depth=0 file=NC.J035 nsew=1
#rec sta=NC.J036 lon=-122.2978 lat=37.5284 depth=0 file=NC.J036 nsew=1
#rec sta=NC.J037 lon=-122.2375 lat=37.5363 depth=0 file=NC.J037 nsew=1
#rec sta=NC.J038 lon=-121.9659 lat=37.3321 depth=0 file=NC.J038 nsew=1
#rec sta=NC.J039 lon=-122.4528 lat=37.6212 depth=0 file=NC.J039 nsew=1
#rec sta=NC.J040 lon=-122.4347 lat=37.8047 depth=0 file=NC.J040 nsew=1
#rec sta=NC.J041 lon=-122.4414 lat=37.8005 depth=0 file=NC.J041 nsew=1
#rec sta=NC.J042 lon=-122.4360 lat=37.8062 depth=0 file=NC.J042 nsew=1
#rec sta=NC.J043 lon=-122.2955 lat=37.5601 depth=0 file=NC.J043 nsew=1
#rec sta=NC.J044 lon=-122.2644 lat=37.3710 depth=0 file=NC.J044 nsew=1
#rec sta=NC.J045 lon=-122.3010 lat=37.5236 depth=0 file=NC.J045 nsew=1
#rec sta=NC.J046 lon=-121.9690 lat=37.2566 depth=0 file=NC.J046 nsew=1
#rec sta=NC.J049 lon=-122.4587 lat=37.5881 depth=0 file=NC.J049 nsew=1
#rec sta=NC.J050 lon=-122.4779 lat=37.6111 depth=0 file=NC.J050 nsew=1
#rec sta=NC.J051 lon=-122.0078 lat=37.3129 depth=0 file=NC.J051 nsew=1
#rec sta=NC.J052 lon=-122.2446 lat=37.2996 depth=0 file=NC.J052 nsew=1
#rec sta=NC.J053 lon=-122.1606 lat=37.2884 depth=0 file=NC.J053 nsew=1
#rec sta=NC.J055 lon=-121.9221 lat=37.4000 depth=0 file=NC.J055 nsew=1
#rec sta=NC.J056 lon=-122.4429 lat=37.7783 depth=0 file=NC.J056 nsew=1
#rec sta=NC.J057 lon=-122.4639 lat=37.7880 depth=0 file=NC.J057 nsew=1
#rec sta=NC.J059 lon=-122.4528 lat=37.4893 depth=0 file=NC.J059 nsew=1
#rec sta=NC.J060 lon=-122.4053 lat=37.7127 depth=0 file=NC.J060 nsew=1
#rec sta=NC.J061 lon=-122.4974 lat=37.7306 depth=0 file=NC.J061 nsew=1
#rec sta=NC.J062 lon=-122.2658 lat=37.5626 depth=0 file=NC.J062 nsew=1
#rec sta=NC.J064 lon=-122.3240 lat=37.4203 depth=0 file=NC.J064 nsew=1
#rec sta=NC.J065 lon=-122.0512 lat=37.3432 depth=0 file=NC.J065 nsew=1
#rec sta=NC.J067 lon=-122.3410 lat=37.5541 depth=0 file=NC.J067 nsew=1
#rec sta=NC.JBG lon=-122.3398 lat=37.3430 depth=0 file=NC.JBG nsew=1
#rec sta=NC.JBMB lon=-122.1531 lat=37.3186 depth=0 file=NC.JBMB nsew=1
#rec sta=NC.JCH lon=-122.3773 lat=37.5170 depth=0 file=NC.JCH nsew=1
#rec sta=NC.JFP lon=-122.1899 lat=37.3604 depth=0 file=NC.JFP nsew=1
#rec sta=NC.JGR lon=-122.4582 lat=37.5160 depth=0 file=NC.JGR nsew=1
#rec sta=NC.JMGB lon=-122.4748 lat=37.6372 depth=0 file=NC.JMGB nsew=1
#rec sta=NC.JMP lon=-122.1664 lat=37.4554 depth=0 file=NC.JMP nsew=1
#rec sta=NC.JPC lon=-122.2128 lat=37.2708 depth=0 file=NC.JPC nsew=1
#rec sta=NC.JPR lon=-122.4750 lat=37.7953 depth=0 file=NC.JPR nsew=1
#rec sta=NC.JSA lon=-122.4201 lat=37.5822 depth=0 file=NC.JSA nsew=1
#rec sta=NC.JSB lon=-122.3986 lat=37.6743 depth=0 file=NC.JSB nsew=1
#rec sta=NC.JSF lon=-122.1741 lat=37.4046 depth=0 file=NC.JSF nsew=1
#rec sta=NC.JSFB lon=-122.1760 lat=37.4037 depth=0 file=NC.JSFB nsew=1
#rec sta=NC.JSGB lon=-122.0502 lat=37.2840 depth=0 file=NC.JSGB nsew=1
#rec sta=NC.JSP lon=-122.5022 lat=37.5864 depth=0 file=NC.JSP nsew=1
#rec sta=NC.N001 lon=-122.4941 lat=37.9604 depth=0 file=NC.N001 nsew=1
#rec sta=NC.N002 lon=-122.1407 lat=38.1826 depth=0 file=NC.N002 nsew=1
#rec sta=NC.N003 lon=-122.5538 lat=38.1090 depth=0 file=NC.N003 nsew=1
#rec sta=NC.N004 lon=-122.6614 lat=38.4490 depth=0 file=NC.N004 nsew=1
#rec sta=NC.N005 lon=-122.7884 lat=38.5341 depth=0 file=NC.N005 nsew=1
#rec sta=NC.N006 lon=-122.7092 lat=38.4063 depth=0 file=NC.N006 nsew=1
#rec sta=NC.N007 lon=-122.6952 lat=38.4464 depth=0 file=NC.N007 nsew=1
#rec sta=NC.N008 lon=-122.6821 lat=38.4284 depth=0 file=NC.N008 nsew=1
#rec sta=NC.N009 lon=-122.6316 lat=38.3514 depth=0 file=NC.N009 nsew=1
#rec sta=NC.N010 lon=-122.7593 lat=38.4529 depth=0 file=NC.N010 nsew=1
#rec sta=NC.N011 lon=-122.6808 lat=38.3037 depth=0 file=NC.N011 nsew=1
#rec sta=NC.N012 lon=-122.4903 lat=38.3138 depth=0 file=NC.N012 nsew=1
#rec sta=NC.N013 lon=-122.5502 lat=38.2993 depth=0 file=NC.N013 nsew=1
#rec sta=NC.N014 lon=-122.6099 lat=38.2473 depth=0 file=NC.N014 nsew=1
#rec sta=NC.N015 lon=-122.7201 lat=38.4477 depth=0 file=NC.N015 nsew=1
#rec sta=NC.N016 lon=-122.2848 lat=38.2988 depth=0 file=NC.N016 nsew=1
#rec sta=NC.N017 lon=-122.6587 lat=38.0530 depth=0 file=NC.N017 nsew=1
#rec sta=NC.N018 lon=-122.8544 lat=38.0969 depth=0 file=NC.N018 nsew=1
#rec sta=NC.N019 lon=-122.4013 lat=38.3011 depth=0 file=NC.N019 nsew=1
#rec sta=NC.NAPC lon=-122.2527 lat=38.4395 depth=0 file=NC.NAPC nsew=1
#rec sta=NC.NBO lon=-122.7181 lat=37.9214 depth=0 file=NC.NBO nsew=1
#rec sta=NC.NBP lon=-122.1947 lat=38.6681 depth=0 file=NC.NBP nsew=1
#rec sta=NC.NBRB lon=-122.5518 lat=38.2602 depth=0 file=NC.NBRB nsew=1
#rec sta=NC.NCC lon=-122.4747 lat=38.0020 depth=0 file=NC.NCC nsew=1
#rec sta=NC.NEA lon=-122.9430 lat=38.3099 depth=0 file=NC.NEA nsew=1
#rec sta=NC.NEH lon=-122.8108 lat=38.3366 depth=0 file=NC.NEH nsew=1
#rec sta=NC.NFV lon=-122.8965 lat=38.4677 depth=0 file=NC.NFV nsew=1
#rec sta=NC.NGVB lon=-122.2158 lat=38.2804 depth=0 file=NC.NGVB nsew=1
#rec sta=NC.NHC lon=-122.3577 lat=38.2175 depth=0 file=NC.NHC nsew=1
#rec sta=NC.NHF lon=-122.5233 lat=38.0496 depth=0 file=NC.NHF nsew=1
#rec sta=NC.NHM lon=-121.8044 lat=38.1549 depth=0 file=NC.NHM nsew=1
#rec sta=NC.NHS lon=-122.6144 lat=38.6565 depth=0 file=NC.NHS nsew=1
#rec sta=NC.NHV lon=-122.7130 lat=38.1649 depth=0 file=NC.NHV nsew=1
#rec sta=NC.NLH lon=-122.1492 lat=38.1228 depth=0 file=NC.NLH nsew=1
#rec sta=NC.NMH lon=-122.6332 lat=38.6694 depth=0 file=NC.NMH nsew=1
#rec sta=NC.NMI lon=-122.2585 lat=38.0760 depth=0 file=NC.NMI nsew=1
#rec sta=NC.NOLB lon=-122.7968 lat=38.0425 depth=0 file=NC.NOLB nsew=1
#rec sta=NC.NPRB lon=-123.0195 lat=37.9964 depth=0 file=NC.NPRB nsew=1
#rec sta=NC.NSM lon=-122.5780 lat=38.3488 depth=0 file=NC.NSM nsew=1
#rec sta=NC.NSP lon=-122.4644 lat=38.2002 depth=0 file=NC.NSP nsew=1
#rec sta=NC.NTA lon=-122.5983 lat=37.9231 depth=0 file=NC.NTA nsew=1
#rec sta=NC.NTAB lon=-122.5983 lat=37.9231 depth=0 file=NC.NTAB nsew=1
#rec sta=NC.NTAC lon=-122.5965 lat=37.9237 depth=0 file=NC.NTAC nsew=1
#rec sta=NC.NTO lon=-122.4497 lat=38.1436 depth=0 file=NC.NTO nsew=1
#rec sta=NC.NTR lon=-122.7984 lat=38.2404 depth=0 file=NC.NTR nsew=1
#rec sta=NC.NTYB lon=-122.6630 lat=38.3892 depth=0 file=NC.NTYB nsew=1
#rec sta=NC.RWSVT lon=-122.2622 lat=37.9417 depth=0 file=NC.RWSVT nsew=1
#rec sta=NP.1002 lon=-122.2552 lat=37.5289 depth=0 file=NP.1002 nsew=1
#rec sta=NP.1003 lon=-122.2977 lat=37.4940 depth=0 file=NP.1003 nsew=1
#rec sta=NP.1005 lon=-122.2356 lat=37.8762 depth=0 file=NP.1005 nsew=1
#rec sta=NP.1103 lon=-122.2688 lat=37.8700 depth=0 file=NP.1103 nsew=1
#rec sta=NP.1129 lon=-122.0825 lat=37.6808 depth=0 file=NP.1129 nsew=1
#rec sta=NP.1216 lon=-123.0091 lat=38.7170 depth=0 file=NP.1216 nsew=1
#rec sta=NP.1225 lon=-122.5048 lat=37.7827 depth=0 file=NP.1225 nsew=1
#rec sta=NP.1226 lon=-121.7630 lat=37.6255 depth=0 file=NP.1226 nsew=1
#rec sta=NP.1239 lon=-122.4028 lat=37.7952 depth=0 file=NP.1239 nsew=1
#rec sta=NP.1446 lon=-122.4004 lat=37.7896 depth=0 file=NP.1446 nsew=1
#rec sta=NP.1447 lon=-122.1416 lat=37.4057 depth=0 file=NP.1447 nsew=1
#rec sta=NP.1448 lon=-122.1166 lat=37.9935 depth=0 file=NP.1448 nsew=1
#rec sta=NP.1515 lon=-122.2480 lat=37.5550 depth=0 file=NP.1515 nsew=1
#rec sta=NP.1571 lon=-121.8520 lat=37.3393 depth=0 file=NP.1571 nsew=1
#rec sta=NP.1590 lon=-122.5080 lat=37.9460 depth=0 file=NP.1590 nsew=1
#rec sta=NP.1662 lon=-122.2961 lat=37.8398 depth=0 file=NP.1662 nsew=1
#rec sta=NP.1675 lon=-122.3850 lat=37.7275 depth=0 file=NP.1675 nsew=1
#rec sta=NP.1676 lon=-122.4164 lat=37.7934 depth=0 file=NP.1676 nsew=1
#rec sta=NP.1678 lon=-122.4758 lat=37.8068 depth=0 file=NP.1678 nsew=1
#rec sta=NP.1684 lon=-121.8312 lat=37.5149 depth=0 file=NP.1684 nsew=1
#rec sta=NP.1687 lon=-121.8003 lat=37.4504 depth=0 file=NP.1687 nsew=1
#rec sta=NP.1688 lon=-121.8809 lat=37.5968 depth=0 file=NP.1688 nsew=1
#rec sta=NP.1689 lon=-121.9333 lat=37.7090 depth=0 file=NP.1689 nsew=1
#rec sta=NP.1690 lon=-121.9929 lat=37.8094 depth=0 file=NP.1690 nsew=1
#rec sta=NP.1691 lon=-122.0785 lat=37.9266 depth=0 file=NP.1691 nsew=1
#rec sta=NP.1695 lon=-122.0248 lat=37.4016 depth=0 file=NP.1695 nsew=1
#rec sta=NP.1696 lon=-121.7580 lat=37.3963 depth=0 file=NP.1696 nsew=1
#rec sta=NP.1699 lon=-122.1520 lat=37.3190 depth=0 file=NP.1699 nsew=1
#rec sta=NP.1720 lon=-122.0848 lat=37.2991 depth=0 file=NP.1720 nsew=1
#rec sta=NP.1721 lon=-122.4435 lat=37.8045 depth=0 file=NP.1721 nsew=1
#rec sta=NP.1722 lon=-122.3660 lat=37.9786 depth=0 file=NP.1722 nsew=1
#rec sta=NP.1723 lon=-122.2225 lat=37.8841 depth=0 file=NP.1723 nsew=1
#rec sta=NP.1726 lon=-122.2610 lat=37.5630 depth=0 file=NP.1726 nsew=1
#rec sta=NP.1729 lon=-122.3986 lat=37.6745 depth=0 file=NP.1729 nsew=1
#rec sta=NP.1730 lon=-122.1193 lat=37.7242 depth=0 file=NP.1730 nsew=1
#rec sta=NP.1731 lon=-121.8873 lat=37.5404 depth=0 file=NP.1731 nsew=1
#rec sta=NP.1732 lon=-122.0100 lat=37.6450 depth=0 file=NP.1732 nsew=1
#rec sta=NP.1734 lon=-122.2598 lat=37.3545 depth=0 file=NP.1734 nsew=1
#rec sta=NP.1735 lon=-122.4099 lat=37.7804 depth=0 file=NP.1735 nsew=1
#rec sta=NP.1736 lon=-122.2000 lat=37.3400 depth=0 file=NP.1736 nsew=1
#rec sta=NP.1737 lon=-122.2995 lat=37.9338 depth=0 file=NP.1737 nsew=1
#rec sta=NP.1738 lon=-122.1672 lat=37.4546 depth=0 file=NP.1738 nsew=1
#rec sta=NP.1739 lon=-122.0066 lat=37.6000 depth=0 file=NP.1739 nsew=1
#rec sta=NP.1740 lon=-122.4080 lat=37.8030 depth=0 file=NP.1740 nsew=1
#rec sta=NP.1741 lon=-122.4432 lat=37.8034 depth=0 file=NP.1741 nsew=1
#rec sta=NP.1742 lon=-121.9047 lat=37.3497 depth=0 file=NP.1742 nsew=1
#rec sta=NP.1743 lon=-122.6576 lat=38.2667 depth=0 file=NP.1743 nsew=1
#rec sta=NP.1744 lon=-122.8723 lat=38.6162 depth=0 file=NP.1744 nsew=1
#rec sta=NP.1745 lon=-122.1703 lat=37.4570 depth=0 file=NP.1745 nsew=1
#rec sta=NP.1749 lon=-122.4165 lat=37.9499 depth=0 file=NP.1749 nsew=1
#rec sta=NP.1750 lon=-122.0909 lat=37.5531 depth=0 file=NP.1750 nsew=1
#rec sta=NP.1751 lon=-122.5392 lat=38.0667 depth=0 file=NP.1751 nsew=1
#rec sta=NP.1752 lon=-122.3079 lat=37.4678 depth=0 file=NP.1752 nsew=1
#rec sta=NP.1753 lon=-122.2485 lat=37.5596 depth=0 file=NP.1753 nsew=1
#rec sta=NP.1754 lon=-122.0812 lat=37.6719 depth=0 file=NP.1754 nsew=1
#rec sta=NP.1755 lon=-122.2438 lat=37.7627 depth=0 file=NP.1755 nsew=1
#rec sta=NP.1756 lon=-122.2526 lat=37.7377 depth=0 file=NP.1756 nsew=1
#rec sta=NP.1757 lon=-121.9254 lat=37.5432 depth=0 file=NP.1757 nsew=1
#rec sta=NP.1759 lon=-122.2562 lat=38.1077 depth=0 file=NP.1759 nsew=1
#rec sta=NP.1760 lon=-122.1574 lat=38.0541 depth=0 file=NP.1760 nsew=1
#rec sta=NP.1761 lon=-122.4572 lat=38.2901 depth=0 file=NP.1761 nsew=1
#rec sta=NP.1762 lon=-122.5656 lat=38.0975 depth=0 file=NP.1762 nsew=1
#rec sta=NP.1764 lon=-122.4712 lat=38.5068 depth=0 file=NP.1764 nsew=1
#rec sta=NP.1765 lon=-122.3184 lat=38.3305 depth=0 file=NP.1765 nsew=1
#rec sta=NP.1767 lon=-122.7046 lat=38.4407 depth=0 file=NP.1767 nsew=1
#rec sta=NP.1768 lon=-122.6367 lat=38.2336 depth=0 file=NP.1768 nsew=1
#rec sta=NP.1769 lon=-121.9571 lat=38.2554 depth=0 file=NP.1769 nsew=1
#rec sta=NP.1770 lon=-122.4152 lat=37.7602 depth=0 file=NP.1770 nsew=1
#rec sta=NP.1771 lon=-122.3969 lat=37.7774 depth=0 file=NP.1771 nsew=1
#rec sta=NP.1772 lon=-122.4907 lat=37.7509 depth=0 file=NP.1772 nsew=1
#rec sta=NP.1773 lon=-122.4736 lat=37.7639 depth=0 file=NP.1773 nsew=1
#rec sta=NP.1774 lon=-122.5045 lat=37.7615 depth=0 file=NP.1774 nsew=1
#rec sta=NP.1775 lon=-122.0946 lat=37.4070 depth=0 file=NP.1775 nsew=1
#rec sta=NP.1776 lon=-122.2392 lat=37.4748 depth=0 file=NP.1776 nsew=1
#rec sta=NP.1777 lon=-122.3047 lat=37.5419 depth=0 file=NP.1777 nsew=1
#rec sta=NP.1778 lon=-121.9798 lat=37.3025 depth=0 file=NP.1778 nsew=1
#rec sta=NP.1779 lon=-121.8666 lat=37.2675 depth=0 file=NP.1779 nsew=1
#rec sta=NP.1780 lon=-121.9162 lat=37.2635 depth=0 file=NP.1780 nsew=1
#rec sta=NP.1781 lon=-121.9039 lat=37.3014 depth=0 file=NP.1781 nsew=1
#rec sta=NP.1782 lon=-121.9435 lat=37.3182 depth=0 file=NP.1782 nsew=1
#rec sta=NP.1783 lon=-121.8492 lat=37.3066 depth=0 file=NP.1783 nsew=1
#rec sta=NP.1784 lon=-122.1725 lat=37.4555 depth=0 file=NP.1784 nsew=1
#rec sta=NP.1785 lon=-121.8737 lat=37.6608 depth=0 file=NP.1785 nsew=1
#rec sta=NP.1786 lon=-121.7396 lat=37.6803 depth=0 file=NP.1786 nsew=1
#rec sta=NP.1787 lon=-122.2061 lat=37.4179 depth=0 file=NP.1787 nsew=1
#rec sta=NP.1788 lon=-121.9154 lat=37.4170 depth=0 file=NP.1788 nsew=1
#rec sta=NP.1789 lon=-121.9846 lat=37.3923 depth=0 file=NP.1789 nsew=1
#rec sta=NP.1790 lon=-122.3674 lat=37.5909 depth=0 file=NP.1790 nsew=1
#rec sta=NP.1791 lon=-122.0389 lat=37.3650 depth=0 file=NP.1791 nsew=1
#rec sta=NP.1792 lon=-122.4253 lat=37.7477 depth=0 file=NP.1792 nsew=1
#rec sta=NP.1793 lon=-121.8078 lat=37.3465 depth=0 file=NP.1793 nsew=1
#rec sta=NP.1794 lon=-122.0011 lat=37.4072 depth=0 file=NP.1794 nsew=1
#rec sta=NP.1795 lon=-122.3868 lat=37.7466 depth=0 file=NP.1795 nsew=1
#rec sta=NP.1796 lon=-122.0604 lat=37.3941 depth=0 file=NP.1796 nsew=1
#rec sta=NP.1798 lon=-121.8437 lat=37.3760 depth=0 file=NP.1798 nsew=1
#rec sta=NP.1799 lon=-121.9862 lat=37.2711 depth=0 file=NP.1799 nsew=1
#rec sta=NP.1800 lon=-122.4158 lat=37.7791 depth=0 file=NP.1800 nsew=1
#rec sta=NP.1801 lon=-122.4287 lat=37.4510 depth=0 file=NP.1801 nsew=1
#rec sta=NP.1802 lon=-122.0745 lat=37.5927 depth=0 file=NP.1802 nsew=1
#rec sta=NP.1803 lon=-121.7270 lat=37.7174 depth=0 file=NP.1803 nsew=1
#rec sta=NP.1804 lon=-122.4840 lat=37.6793 depth=0 file=NP.1804 nsew=1
#rec sta=NP.1805 lon=-122.3030 lat=37.7845 depth=0 file=NP.1805 nsew=1
#rec sta=NP.1806 lon=-122.5368 lat=37.9103 depth=0 file=NP.1806 nsew=1
#rec sta=NP.1807 lon=-122.2992 lat=37.8877 depth=0 file=NP.1807 nsew=1
#rec sta=NP.1810 lon=-122.4135 lat=37.7374 depth=0 file=NP.1810 nsew=1
#rec sta=NP.1811 lon=-122.1540 lat=37.4464 depth=0 file=NP.1811 nsew=1
#rec sta=NP.1812 lon=-122.2711 lat=37.8693 depth=0 file=NP.1812 nsew=1
#rec sta=NP.1813 lon=-122.2702 lat=38.0343 depth=0 file=NP.1813 nsew=1
#rec sta=NP.1814 lon=-122.5092 lat=37.9445 depth=0 file=NP.1814 nsew=1
#rec sta=NP.1815 lon=-122.2637 lat=37.4891 depth=0 file=NP.1815 nsew=1
#rec sta=NP.1816 lon=-122.4013 lat=37.7448 depth=0 file=NP.1816 nsew=1
#rec sta=NP.1817 lon=-122.4860 lat=37.7793 depth=0 file=NP.1817 nsew=1
#rec sta=NP.1818 lon=-122.4312 lat=37.7163 depth=0 file=NP.1818 nsew=1
#rec sta=NP.1819 lon=-122.2925 lat=37.8108 depth=0 file=NP.1819 nsew=1
#rec sta=NP.1820 lon=-122.2761 lat=37.8251 depth=0 file=NP.1820 nsew=1
#rec sta=NP.1821 lon=-122.2684 lat=37.7986 depth=0 file=NP.1821 nsew=1
#rec sta=NP.1822 lon=-122.3220 lat=37.8965 depth=0 file=NP.1822 nsew=1
#rec sta=NP.1823 lon=-122.1147 lat=37.9945 depth=0 file=NP.1823 nsew=1
#rec sta=NP.1824 lon=-122.4484 lat=37.7242 depth=0 file=NP.1824 nsew=1
#rec sta=NP.1825 lon=-121.8709 lat=37.4140 depth=0 file=NP.1825 nsew=1
#rec sta=NP.1826 lon=-121.9210 lat=37.6806 depth=0 file=NP.1826 nsew=1
#rec sta=NP.1827 lon=-122.2933 lat=37.8606 depth=0 file=NP.1827 nsew=1
#rec sta=NP.1828 lon=-122.2525 lat=37.8583 depth=0 file=NP.1828 nsew=1
#rec sta=NP.1829 lon=-122.4613 lat=38.2903 depth=0 file=NP.1829 nsew=1
#rec sta=NP.1830 lon=-122.0781 lat=37.3742 depth=0 file=NP.1830 nsew=1
#rec sta=NP.1831 lon=-122.0182 lat=37.3641 depth=0 file=NP.1831 nsew=1
#rec sta=NP.1833 lon=-122.4039 lat=37.7831 depth=0 file=NP.1833 nsew=1
#rec sta=NP.1834 lon=-122.0777 lat=37.4255 depth=0 file=NP.1834 nsew=1
#rec sta=NP.1835 lon=-122.6068 lat=38.4421 depth=0 file=NP.1835 nsew=1
#rec sta=NP.1836 lon=-122.1698 lat=37.6917 depth=0 file=NP.1836 nsew=1
#rec sta=NP.1838 lon=-121.9109 lat=37.3833 depth=0 file=NP.1838 nsew=1
#rec sta=NP.1839 lon=-121.9360 lat=37.2854 depth=0 file=NP.1839 nsew=1
#rec sta=NP.1841 lon=-121.9778 lat=37.3583 depth=0 file=NP.1841 nsew=1
#rec sta=NP.1842 lon=-121.7835 lat=37.6887 depth=0 file=NP.1842 nsew=1
#rec sta=NP.1843 lon=-122.1509 lat=37.7269 depth=0 file=NP.1843 nsew=1
#rec sta=NP.1844 lon=-122.0322 lat=37.8852 depth=0 file=NP.1844 nsew=1
#rec sta=NP.1845 lon=-122.5784 lat=38.5780 depth=0 file=NP.1845 nsew=1
#rec sta=NP.1846 lon=-121.9305 lat=37.9422 depth=0 file=NP.1846 nsew=1
#rec sta=NP.1847 lon=-122.1342 lat=38.0130 depth=0 file=NP.1847 nsew=1
#rec sta=NP.1848 lon=-122.5240 lat=38.3674 depth=0 file=NP.1848 nsew=1
#rec sta=NP.1849 lon=-122.1614 lat=37.4673 depth=0 file=NP.1849 nsew=1
#rec sta=NP.1851 lon=-122.1953 lat=37.4183 depth=0 file=NP.1851 nsew=1
#rec sta=NP.1852 lon=-122.1952 lat=37.4173 depth=0 file=NP.1852 nsew=1
#rec sta=NP.1853 lon=-122.1507 lat=37.4462 depth=0 file=NP.1853 nsew=1
#rec sta=NP.1854 lon=-122.1976 lat=37.7590 depth=0 file=NP.1854 nsew=1
#rec sta=NP.1855 lon=-122.1734 lat=37.7696 depth=0 file=NP.1855 nsew=1
#rec sta=NP.1856 lon=-122.2303 lat=37.8028 depth=0 file=NP.1856 nsew=1
#rec sta=NP.1857 lon=-122.2284 lat=37.7777 depth=0 file=NP.1857 nsew=1
#rec sta=NP.1858 lon=-122.2014 lat=37.7298 depth=0 file=NP.1858 nsew=1
#rec sta=NP.1859 lon=-122.1705 lat=37.7448 depth=0 file=NP.1859 nsew=1
#rec sta=NP.1860 lon=-122.1971 lat=37.7857 depth=0 file=NP.1860 nsew=1
#rec sta=NP.1861 lon=-122.2638 lat=37.7921 depth=0 file=NP.1861 nsew=1
#rec sta=NP.1862 lon=-122.0849 lat=37.2995 depth=0 file=NP.1862 nsew=1
#rec sta=NP.1863 lon=-121.8571 lat=37.6667 depth=0 file=NP.1863 nsew=1
#rec sta=NP.1866 lon=-122.4122 lat=37.7796 depth=0 file=NP.1866 nsew=1
#rec sta=NP.1867 lon=-122.2572 lat=37.8744 depth=0 file=NP.1867 nsew=1
#rec sta=NP.1868 lon=-122.1399 lat=37.4054 depth=0 file=NP.1868 nsew=1
#rec sta=NP.1869 lon=-122.1406 lat=37.4061 depth=0 file=NP.1869 nsew=1
#rec sta=NP.1870 lon=-122.1609 lat=37.4654 depth=0 file=NP.1870 nsew=1
#rec sta=NP.1872 lon=-122.0777 lat=37.4255 depth=0 file=NP.1872 nsew=1
#rec sta=NP.1874 lon=-121.7628 lat=37.6275 depth=0 file=NP.1874 nsew=1
#rec sta=NP.1875 lon=-122.5049 lat=37.7824 depth=0 file=NP.1875 nsew=1
#rec sta=NP.1876 lon=-122.3988 lat=37.7897 depth=0 file=NP.1876 nsew=1
#rec sta=NP.1877 lon=-121.9258 lat=37.5424 depth=0 file=NP.1877 nsew=1
#rec sta=NP.BCS lon=-122.4066 lat=37.7773 depth=0 file=NP.BCS nsew=1
#rec sta=NP.BMA lon=-122.4145 lat=37.8068 depth=0 file=NP.BMA nsew=1
#rec sta=NP.BOA lon=-122.3966 lat=37.7781 depth=0 file=NP.BOA nsew=1
#rec sta=NP.CES lon=-122.3873 lat=37.7511 depth=0 file=NP.CES nsew=1
#rec sta=NP.CPS lon=-122.3987 lat=37.7711 depth=0 file=NP.CPS nsew=1
#rec sta=NP.DIX lon=-121.8424 lat=38.3771 depth=0 file=NP.DIX nsew=1
#rec sta=NP.ELK lon=-121.3621 lat=38.3735 depth=0 file=NP.ELK nsew=1
#rec sta=NP.EMB lon=-122.3954 lat=37.7957 depth=0 file=NP.EMB nsew=1
#rec sta=NP.EMR lon=-121.4993 lat=38.0605 depth=0 file=NP.EMR nsew=1
#rec sta=NP.ENM lon=-122.3975 lat=37.7841 depth=0 file=NP.ENM nsew=1
#rec sta=NP.FM3 lon=-122.4280 lat=37.8074 depth=0 file=NP.FM3 nsew=1
#rec sta=NP.LEV lon=-122.4003 lat=37.8018 depth=0 file=NP.LEV nsew=1
#rec sta=NP.MAR lon=-122.4243 lat=37.8064 depth=0 file=NP.MAR nsew=1
#rec sta=NP.MOD lon=-121.1223 lat=37.6202 depth=0 file=NP.MOD nsew=1
#rec sta=NP.PBS lon=-122.3897 lat=37.7779 depth=0 file=NP.PBS nsew=1
#rec sta=NP.PHG lon=-122.3970 lat=37.7563 depth=0 file=NP.PHG nsew=1
#rec sta=NP.PLA lon=-121.4631 lat=37.7987 depth=0 file=NP.PLA nsew=1
#rec sta=NP.RH2 lon=-122.4217 lat=37.8022 depth=0 file=NP.RH2 nsew=1
#rec sta=NP.TEL lon=-122.4035 lat=37.8020 depth=0 file=NP.TEL nsew=1
#rec sta=NP.WSS lon=-122.4429 lat=37.8036 depth=0 file=NP.WSS nsew=1

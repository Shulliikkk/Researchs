double h[] = {0,  500,  1000,  1500,  2000,  2500,  3000,  4000,  5000,  6000,  7000,  8000,  9000,  10000,  11000,  12000,  14000,  16000,  18000, 20000.0, 22000.0, 24000.0, 26000.0, 28000.0, 30000.0, 32000.0, 34000.0, 36000.0, 38000.0, 40000.0, 42000.0, 44000.0, 46000.0, 48000.0, 50000.0, 52000.0, 54000.0, 56000.0, 58000.0, 60000.0, 62000.0, 64000.0, 66000.0, 68000.0, 70000.0, 72000.0, 74000.0, 76000.0, 78000.0, 80000.0, 82000.0, 84000.0, 86000.0, 88000.0, 90000.0, 92000.0, 94000.0, 96000.0, 98000.0, 100000.0, 102000.0, 104000.0, 106000.0, 108000.0, 110000.0, 112000.0, 114000.0, 116000.0, 118000.0, 120000.0, 122000.0, 124000.0, 126000.0, 128000.0, 130000.0, 132000.0, 134000.0, 136000.0, 138000.0, 140000.0, 142000.0, 144000.0, 146000.0, 148000.0, 150000.0, 152000.0, 154000.0, 156000.0, 158000.0, 160000.0, 162000.0, 164000.0, 166000.0, 168000.0, 170000.0, 172000.0, 174000.0, 176000.0, 178000.0, 180000.0, 182000.0, 184000.0, 186000.0, 188000.0, 190000.0, 192000.0, 194000.0, 196000.0, 198000.0, 200000.0, 202000.0, 204000.0, 206000.0, 208000.0, 210000.0, 212000.0, 214000.0, 216000.0, 218000.0, 220000.0, 222000.0, 224000.0, 226000.0, 228000.0, 230000.0, 232000.0, 234000.0, 236000.0, 238000.0, 240000.0, 242000.0, 244000.0, 246000.0, 248000.0, 250000.0, 252000.0, 254000.0, 256000.0, 258000.0, 260000.0, 262000.0, 264000.0, 266000.0, 268000.0, 270000.0, 272000.0, 274000.0, 276000.0, 278000.0, 280000.0, 282000.0, 284000.0, 286000.0, 288000.0, 290000.0, 292000.0, 294000.0, 296000.0, 298000.0, 300000.0, 302000.0, 304000.0, 306000.0, 308000.0, 310000.0, 312000.0, 314000.0, 316000.0, 318000.0, 320000.0, 322000.0, 324000.0, 326000.0, 328000.0, 330000.0, 332000.0, 334000.0, 340000.0, 342000.0, 344000.0, 346000.0, 348000.0, 350000.0, 352000.0, 354000.0, 356000.0, 358000.0, 360000.0, 362000.0, 364000.0, 366000.0, 368000.0, 370000.0, 372000.0, 374000.0, 380000.0, 386000.0, 392000.0, 398000.0, 404000.0, 410000.0, 416000.0, 422000.0, 428000.0, 434000.0, 440000.0, 446000.0, 452000.0, 458000.0, 464000.0, 470000.0, 476000.0, 482000.0, 486000.0, 492000.0, 498000.0, 510000.0, 525000.0, 540000.0, 555000.0, 570000.0, 585000.0, 600000.0, 615000.0, 630000.0, 645000.0, 660000.0, 675000.0, 690000.0, 705000.0, 720000.0, 735000.0, 750000.0, 765000.0, 780000.0, 795000.0, 810000.0, 825000.0, 840000.0, 850000.0, 865000.0, 880000.0, 895000.0, 910000.0, 925000.0, 940000.0, 950000.0, 965000.0, 980000.0, 995000.0, 1010000.0, 1025000.0, 1040000.0, 1055000.0, 1070000.0, 1085000.0, 1100000.0, 1115000.0, 1130000.0, 1145000.0, 1160000.0, 1175000.0, 1190000.0, 1200000.0};
double ro[] = {1.2250,  1.1673,  1.1117,  1.0581,  1.0065,  0.9569,  0.9093,  0.8194,  0.7365,  0.6601,  0.59,  0.5258,  0.4671,  0.4135,  0.3648,  0.3119,  0.2279,  0.1665,  0.1216, 0.269405, 0.22948, 0.195746, 0.0167226, 0.143075, 0.122591, 0.105192, 0.0898405, 0.0769728, 0.066188, 0.0571119, 0.0494439, 0.0429413, 0.0374072, 0.0327849, 0.0289528, 0.0256445, 0.0228393, 0.0202927, 0.0179855, 0.0158996, 0.0140177, 0.0123239, 0.0108028, 0.0094402, 0.00822283, 0.0071356, 0.00615349, 0.00529233, 0.00453904, 1.8458e-05, 1.34175e-05, 9.69458e-06, 6.95782e-06, 4.88072e-06, 3.41817e-06, 2.39477e-06, 1.6784e-06, 1.16617e-06, 7.99924e-07, 5.54951e-07, 3.89128e-07, 2.75614e-07, 1.86514e-07, 1.3006e-07, 9.34035e-08, 6.87683e-08, 5.17245e-08, 3.96418e-08, 3.08814e-08, 2.44041e-08, 2.01267e-08, 1.67199e-08, 1.39015e-08, 1.15819e-08, 9.68064e-09, 8.12642e-09, 6.85639e-09, 5.81605e-09, 4.95889e-09, 4.24614e-09, 3.64654e-09, 3.13605e-09, 2.69768e-09, 2.32127e-09, 2.00329e-09, 1.79797e-09, 1.61741e-09, 1.4577e-09, 1.3166e-09, 1.19204e-09, 1.0821e-09, 9.85037e-10, 8.9924e-10, 8.23245e-10, 7.55732e-10, 6.95515e-10, 6.41542e-10, 5.92889e-10, 5.48759e-10, 5.08478e-10, 4.71493e-10, 4.37368e-10, 4.05781e-10, 3.76524e-10, 3.49499e-10, 3.24717e-10, 3.02294e-10, 2.82451e-10, 2.65513e-10, 2.51904e-10, 2.34594e-10, 2.20229e-10, 2.06859e-10, 1.94416e-10, 1.82838e-10, 1.72062e-10, 1.62033e-10, 1.52694e-10, 1.43996e-10, 1.35889e-10, 1.28328e-10, 1.21269e-10, 1.14674e-10, 1.08504e-10, 1.02724e-10, 9.73033e-11, 9.22112e-11, 8.74209e-11, 8.29076e-11, 7.86493e-11, 7.4626e-11, 7.08203e-11, 6.7217e-11, 6.38031e-11, 6.05679e-11, 5.7671e-11, 5.4874e-11, 5.22299e-11, 4.97302e-11, 4.73671e-11, 4.51328e-11, 4.30201e-11, 4.10218e-11, 3.91313e-11, 3.73422e-11, 3.56482e-11, 3.40437e-11, 3.2523e-11, 3.1081e-11, 2.97125e-11, 2.8413e-11, 2.71778e-11, 2.60029e-11, 2.48843e-11, 2.38182e-11, 2.28013e-11, 2.18302e-11, 2.09021e-11, 2.00141e-11, 1.91637e-11, 1.83486e-11, 1.75668e-11, 1.68164e-11, 1.60957e-11, 1.54032e-11, 1.47378e-11, 1.40984e-11, 1.34842e-11, 1.28945e-11, 1.23289e-11, 1.17871e-11, 1.12692e-11, 1.07752e-11, 1.03055e-11, 9.86051e-12, 9.44101e-12, 9.04785e-12, 8.03773e-12, 7.76214e-12, 7.51986e-12, 7.31282e-12, 7.14309e-12, 7.01292e-12, 6.75194e-12, 6.50003e-12, 6.25809e-12, 6.02573e-12, 5.80258e-12, 5.58829e-12, 5.38249e-12, 5.18486e-12, 4.99506e-12, 4.81277e-12, 4.63769e-12, 4.46952e-12, 4.00363e-12, 3.59013e-12, 3.22253e-12, 2.89508e-12, 2.60748e-12, 2.35038e-12, 2.11891e-12, 1.91031e-12, 1.72241e-12, 1.55357e-12, 1.40268e-12, 1.26918e-12, 1.14872e-12, 1.03985e-12, 9.41207e-13, 8.51973e-13, 7.71357e-13, 6.98612e-13, 6.54137e-13, 5.9296e-13, 5.37862e-13, 4.43567e-13, 3.50196e-13, 2.77744e-13, 2.20971e-13, 1.76131e-13, 1.4085e-13, 1.1396e-13, 9.23975e-14, 7.51152e-14, 6.12504e-14, 5.01275e-14, 4.13788e-14, 3.44783e-14, 2.89823e-14, 2.45532e-14, 2.09405e-14, 1.79647e-14, 1.55035e-14, 1.34806e-14, 1.18561e-14, 1.04935e-14, 9.33686e-15, 8.36124e-15, 7.7972e-15, 7.06247e-15, 6.44261e-15, 5.91959e-15, 5.44733e-15, 5.02514e-15, 4.65704e-15, 4.43696e-15, 4.13892e-15, 3.87353e-15, 3.63544e-15, 3.44759e-15, 3.28502e-15, 3.13235e-15, 2.98845e-15, 2.85243e-15, 2.72358e-15, 2.6017e-15, 4.48661e-15, 2.37814e-15, 2.27613e-15, 2.18042e-15, 2.09089e-15, 2.0074e-15, 1.95503e-15};
int l = sizeof(h)/sizeof(double);
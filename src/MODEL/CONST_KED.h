#ifndef __AAVDP_CONST_KED_H__
#define __AAVDP_CONST_KED_H__

#define SAMPLING_TYPE_NUMBER 36
#define SPACE_GROUP_NUMBER 230
const int PG_SAMPLING_TYPE[SAMPLING_TYPE_NUMBER]={
1,2,3,4,5,5,5,6,5,5,
6,6,7,-1,9,-1,-1,-1,-1,-1,
15,12,17,16,18,-1,19,3,6,6,
8,9,-1,-1,-1,-1};
const int SG_SYMMORPHIC_NUMBER[SPACE_GROUP_NUMBER]={
1,   2,   3,   3,   3,   6,   6,   6,   6,  10,
10,  10,  10,  10,  10,  16,  16,  16,  16,  16,
16,  16,  16,  16,  25,  25,  25,  25,  25,  25,
25,  25,  25,  25,  25,  25,  25,  25,  25,  25,
25,  25,  25,  25,  25,  25,  47,  47,  47,  47,
47,  47,  47,  47,  47,  47,  47,  47,  47,  47,
47,  47,  47,  47,  47,  47,  47,  47,  47,  47,
47,  47,  47,  47,  75,  75,  75,  75,  75,  75,
81,  81,  83,  83,  83,  83,  83,  83,  89,  89,
89,  89,  89,  89,  89,  89,  89,  89,  99,  99,
99,  99,  99,  99,  99,  99,  99,  99,  99,  99,
111, 111, 111, 111, 115, 115, 115, 115, 115, 115,
111, 111, 123, 123, 123, 123, 123, 123, 123, 123,
123, 123, 123, 123, 123, 123, 123, 123, 123, 123,
123, 123, 143, 143, 143, 143, 147, 147, 149, 150,
149, 150, 149, 150, 150, 156, 157, 156, 157, 156,
156, 162, 162, 164, 164, 164, 164, 168, 168, 168,
168, 168, 168, 174, 175, 175, 177, 177, 177, 177,
177, 177, 183, 183, 183, 183, 187, 187, 189, 189,
191, 191, 191, 191, 195, 195, 195, 195, 195, 200,
200, 200, 200, 200, 200, 200, 207, 207, 207, 207,
207, 207, 207, 207, 215, 215, 215, 215, 215, 215,
221, 221, 221, 221, 221, 221, 221, 221, 221, 221};
const int SG_HALL_NUMBER_START[SPACE_GROUP_NUMBER]={
1, 2, 3, 6, 9, 18, 21, 30, 39, 57, 
60, 63, 72, 81, 90, 108, 109, 112, 115, 116, 
119, 122, 123, 124, 125, 128, 134, 137, 143, 149, 
155, 161, 164, 170, 173, 176, 182, 185, 191, 197, 
203, 209, 212, 215, 218, 221, 227, 228, 230, 233, 
239, 245, 251, 257, 263, 266, 269, 275, 278, 284, 
290, 292, 298, 304, 310, 313, 316, 322, 334, 335, 
337, 338, 341, 343, 349, 350, 351, 352, 353, 354, 
355, 356, 357, 358, 359, 361, 363, 364, 366, 367, 
368, 369, 370, 371, 372, 373, 374, 375, 376, 377, 
378, 379, 380, 381, 382, 383, 384, 385, 386, 387, 
388, 389, 390, 391, 392, 393, 394, 395, 396, 397, 
398, 399, 400, 401, 402, 404, 406, 407, 408, 410, 
412, 413, 414, 416, 418, 419, 420, 422, 424, 425, 
426, 428, 430, 431, 432, 433, 435, 436, 438, 439, 
440, 441, 442, 443, 444, 446, 447, 448, 449, 450, 
452, 454, 455, 456, 457, 458, 460, 462, 463, 464, 
465, 466, 467, 468, 469, 470, 471, 472, 473, 474, 
475, 476, 477, 478, 479, 480, 481, 482, 483, 484, 
485, 486, 487, 488, 489, 490, 491, 492, 493, 494, 
495, 497, 498, 500, 501, 502, 503, 504, 505, 506, 
507, 508, 509, 510, 511, 512, 513, 514, 515, 516, 
517, 518, 520, 521, 523, 524, 525, 527, 529, 530};

#define E_TYPE_NUMBER 98
const char E_TYPE[E_TYPE_NUMBER][10]={
"H","He","Li","Be","B",
"C","N","O","F","Ne",
"Na","Mg","Al","Si","P",
"S","Cl","Ar","K","Ca",
"Sc","Ti","V","Cr","Mn",
"Fe","Co","Ni","Cu","Zn",
"Ga","Ge","As","Se","Br",
"Kr","Rb","Sr","Y","Zr",
"Nb","Mo","Tc","Ru","Rh",
"Pd","Ag","Cd","In","Sn",
"Sb","Te","I","Xe","Cs",
"Ba","La","Ce","Pr","Nd",
"Pm","Sm","Eu","Gd","Tb",
"Dy","Ho","Er","Tm","Yb",
"Lu","Hf","Ta","W","Re",
"Os","Ir","Pt","Au","Hg",
"Tl","Pb","Bi","Po","At",
"Rn","Fr","Ra","Ac","Th",
"Pa","U","Np","Pu","Am",
"Cm","Bk","Cf"};
const double E_MASS[E_TYPE_NUMBER]={
1.00794,4.002602,6.941,9.012182,10.811,
12.0107,14.0067,15.9994,18.9984032,20.1797,
22.98976928,24.3050,26.9815386,28.0855,30.973762,
32.065,35.453,39.948,39.0983,40.078,
44.955912,47.867,50.9415,51.9961,54.938045,
55.845,58.933195,58.6934,63.546,65.38,
69.723,72.64,74.92160,78.96,79.904,
83.798,85.4678,87.62,88.90585,91.224,
92.90638,95.96,98.9062,101.07,102.90550,
106.42,107.8682,112.411,114.818,118.710,
121.760,127.60,126.90447,131.293,132.9054519,
137.327,138.90547,140.116,140.90765,144.242,
145.0,150.36,151.964,157.25,158.92535,
162.500,164.93032,167.259,168.93421,173.054,
174.9668,178.49,180.94788,183.84,186.207,
190.23,192.217,195.084,196.966569,200.59,
204.3833,207.2,208.98040,209.0,210.0,
222.0,223.0,226.0,227.0,232.03806,
231.03588,238.02891,237.0,244.0,243.0,
247.0,251.0,252.0};
const double E_WK_AA[E_TYPE_NUMBER][4]={
{0.00427,0.00957,0.00802,0.00209},
{0.01217,0.02616,-0.00884,0.01841},
{0.00251,0.03576,0.00988,0.02370},
{0.01596,0.02959,0.04024,0.01001},
{0.03652,0.01140,0.05677,0.01506},
{0.04102,0.04911,0.05296,0.00061},
{0.04123,0.05740,0.06529,0.00373},
{0.03547,0.03133,0.10865,0.01615},
{0.03957,0.07225,0.09581,0.00792},
{0.02597,0.02197,0.13762,0.05394},
{0.03283,0.08858,0.11688,0.02516},
{0.03833,0.17124,0.03649,0.04134},
{0.04388,0.17743,0.05047,0.03957},
{0.03812,0.17833,0.06280,0.05605},
{0.04166,0.17817,0.09479,0.04463},
{0.04003,0.18346,0.12218,0.03753},
{0.04245,0.17645,0.15814,0.03011},
{0.05011,0.16667,0.17074,0.04358},
{0.04058,0.17582,0.20943,0.02922},
{0.04001,0.17416,0.20986,0.05497},
{0.09685,0.14777,0.20981,0.04852},
{0.06667,0.17356,0.22710,0.05957},
{0.05118,0.16791,0.26700,0.06476},
{0.03204,0.18460,0.30764,0.05052},
{0.03866,0.17782,0.31329,0.06898},
{0.05455,0.16660,0.33208,0.06947},
{0.05942,0.17472,0.34423,0.06828},
{0.06049,0.16600,0.37302,0.07109},
{0.08034,0.15838,0.40116,0.05467},
{0.02948,0.19200,0.42222,0.07480},
{0.16157,0.32976,0.18964,0.06148},
{0.16184,0.35705,0.17618,0.07133},
{0.06190,0.18452,0.41600,0.12793},
{0.15913,0.41583,0.13385,0.10549},
{0.16514,0.41202,0.12900,0.13209},
{0.15798,0.41181,0.14254,0.14987},
{0.16535,0.44674,0.24245,0.03161},
{0.16039,0.44470,0.24661,0.05840},
{0.16619,0.44376,0.25613,0.06797},
{0.16794,0.44505,0.27188,0.07313},
{0.16552,0.45008,0.30474,0.06161},
{0.17327,0.44679,0.32441,0.06143},
{0.16424,0.45046,0.33749,0.07766},
{0.18750,0.44919,0.36323,0.05388},
{0.16081,0.45211,0.40343,0.06140},
{0.16599,0.43951,0.41478,0.08142},
{0.16547,0.44658,0.45401,0.05959},
{0.17154,0.43689,0.46392,0.07725},
{0.15752,0.44821,0.48186,0.08596},
{0.15732,0.44563,0.48507,0.10948},
{0.16971,0.42742,0.48779,0.13653},
{0.14927,0.43729,0.49444,0.16440},
{0.18053,0.44724,0.48163,0.15995},
{0.13141,0.43855,0.50035,0.22299},
{0.31397,0.55648,0.39828,0.04852},
{0.32756,0.53927,0.39830,0.07607},
{0.30887,0.53804,0.42265,0.09559},
{0.28398,0.53568,0.46662,0.10282},
{0.35160,0.56889,0.42010,0.07246},
{0.33810,0.58035,0.44442,0.07413},
{0.35449,0.59626,0.43868,0.07152},
{0.35559,0.60598,0.45165,0.07168},
{0.38379,0.64088,0.41710,0.06708},
{0.40352,0.64303,0.40488,0.08137},
{0.36838,0.64761,0.47222,0.06854},
{0.38514,0.68422,0.44359,0.06775},
{0.37280,0.67528,0.47337,0.08320},
{0.39335,0.70093,0.46774,0.06658},
{0.40587,0.71223,0.46598,0.06847},
{0.39728,0.73368,0.47795,0.06759},
{0.40697,0.73576,0.47481,0.08291},
{0.40122,0.78861,0.44658,0.08799},
{0.41127,0.76965,0.46563,0.10180},
{0.39978,0.77171,0.48541,0.11540},
{0.39130,0.80752,0.48702,0.11041},
{0.40436,0.80701,0.48445,0.12438},
{0.38816,0.80163,0.51922,0.13514},
{0.39551,0.80409,0.53365,0.13485},
{0.40850,0.83052,0.53325,0.11978},
{0.40092,0.85415,0.53346,0.12747},
{0.41872,0.88168,0.54551,0.09404},
{0.43358,0.88007,0.52966,0.12059},
{0.40858,0.87837,0.56392,0.13698},
{0.41637,0.85094,0.57749,0.16700},
{0.38951,0.83297,0.60557,0.20770},
{0.41677,0.88094,0.55170,0.21029},
{0.50089,1.00860,0.51420,0.05996},
{0.47470,0.99363,0.54721,0.09206},
{0.47810,0.98385,0.54905,0.12055},
{0.47903,0.97455,0.55883,0.14309},
{0.48351,0.98292,0.58877,0.12425},
{0.48664,0.98057,0.61483,0.12136},
{0.46078,0.97139,0.66506,0.13012},
{0.49148,0.98583,0.67674,0.09725},
{0.50865,0.98574,0.68109,0.09977},
{0.46259,0.97882,0.73056,0.12723},
{0.46221,0.95749,0.76259,0.14086},
{0.48500,0.95602,0.77234,0.13374}};
const double E_WK_BB[E_TYPE_NUMBER][4]={
{4.17218,16.05892,26.78365,69.45643},
{1.83008,7.20225,16.13585,18.75551},
{0.02620,2.00907,10.80597,130.49226},
{0.38968,1.99268,46.86913,108.84167},
{0.50627,3.68297,27.90586,74.98296},
{0.41335,10.98289,34.80286,177.19113},
{0.29792,7.84094,22.58809,72.59254},
{0.17964,2.60856,11.79972,38.02912},
{0.16403,3.96612,12.43903,40.05053},
{0.09101,0.41253,5.02463,17.52954},
{0.06008,2.07182,7.64444,146.00952},
{0.07424,2.87177,18.06729,97.00854},
{0.09086,2.53252,30.43883,98.26737},
{0.05396,1.86461,22.54263,72.43144},
{0.05564,1.62500,24.45354,64.38264},
{0.05214,1.40793,23.35691,53.59676},
{0.04643,1.15677,19.34091,52.88785},
{0.07991,1.01436,15.67109,39.60819},
{0.03352,0.82984,14.13679,200.97722},
{0.02289,0.71288,11.18914,135.02390},
{0.12527,1.34248,12.43524,131.71112},
{0.05198,0.86467,10.59984,103.56776},
{0.03786,0.57160,8.30305,91.78068},
{0.00240,0.44931,7.92251,86.64058},
{0.01836,0.41203,6.73736,76.30466},
{0.03947,0.43294,6.26864,71.29470},
{0.03962,0.43253,6.05175,68.72437},
{0.03558,0.39976,5.36660,62.46894},
{0.05475,0.45736,5.38252,60.43276},
{0.00137,0.26535,4.48040,54.26088},
{0.10455,2.18391,9.04125,75.16958},
{0.09890,2.06856,9.89926,68.13783},
{0.01642,0.32542,3.51888,44.50604},
{0.07669,1.89297,11.31554,46.32082},
{0.08199,1.76568,9.87254,38.10640},
{0.06939,1.53446,8.98025,33.04365},
{0.07044,1.59236,17.53592,215.26198},
{0.06199,1.41265,14.33812,152.80257},
{0.06364,1.34205,13.66551,125.72522},
{0.06565,1.25292,13.09355,109.50252},
{0.05921,1.15624,13.24924,98.69958},
{0.06162,1.11236,12.76149,90.92026},
{0.05081,0.99771,11.28925,84.28943},
{0.05120,1.08672,12.23172,85.27316},
{0.04662,0.85252,10.51121,74.53949},
{0.04933,0.79381,9.30944,41.17414},
{0.04481,0.75608,9.34354,67.91975},
{0.04867,0.71518,8.40595,64.24400},
{0.03672,0.64379,7.83687,73.37281},
{0.03308,0.60931,7.04977,64.83582},
{0.04023,0.58192,6.29247,55.57061},
{0.02842,0.50687,5.60835,48.28004},
{0.03830,0.58340,6.47550,47.08820},
{0.02097,0.41007,4.52105,37.18178},
{0.07813,1.45053,15.05933,199.48830},
{0.08444,1.40227,13.12939,160.56676},
{0.07206,1.19585,11.55866,127.31371},
{0.05717,0.98756,9.95556,117.31874},
{0.08249,1.43427,12.37363,150.55968},
{0.07081,1.31033,11.44403,144.17706},
{0.07442,1.38680,11.54391,143.72185},
{0.07155,1.34703,11.00432,140.09138},
{0.07794,1.55042,11.89283,142.79585},
{0.08508,1.60712,11.45367,116.64063},
{0.06520,1.32571,10.16884,134.69034},
{0.06850,1.43566,10.57719,131.88972},
{0.06264,1.26756,9.46411,107.50194},
{0.06750,1.35829,9.76480,127.40374},
{0.06958,1.38750,9.41888,122.10940},
{0.06574,1.31578,9.13448,120.98209},
{0.06517,1.29452,8.67569,100.34878},
{0.06213,1.30860,9.18871,91.20213},
{0.06292,1.23499,8.42904,77.59815},
{0.05693,1.15762,7.83077,67.14066},
{0.05145,1.11240,8.33441,65.71782},
{0.05573,1.11159,8.00221,57.35021},
{0.04855,0.99356,7.38693,51.75829},
{0.04981,0.97669,7.38024,44.52068},
{0.05151,1.00803,8.03707,45.01758},
{0.04693,0.98398,7.83562,46.51474},
{0.05161,1.02127,9.18455,64.88177},
{0.05154,1.03252,8.49678,58.79463},
{0.04200,0.90939,7.71158,57.79178},
{0.04661,0.87289,6.84038,51.36000},
{0.04168,0.73697,5.86112,43.78613},
{0.04488,0.83871,6.44020,43.51940},
{0.05786,1.20028,13.85073,172.15909},
{0.05239,1.03225,11.49796,143.12303},
{0.05167,0.98867,10.52682,112.18267},
{0.04931,0.95698,9.61135,95.44649},
{0.04748,0.93369,9.89867,102.06961},
{0.04660,0.89912,9.69785,100.23434},
{0.04323,0.78798,8.71624,92.30811},
{0.04641,0.85867,9.51157,111.02754},
{0.04918,0.87026,9.41105,104.98576},
{0.03904,0.72797,8.00506,86.41747},
{0.03969,0.68167,7.29607,75.72682},
{0.04291,0.69956,7.38554,77.18528}};

//0<S<2
const double E_GP_A1[E_TYPE_NUMBER][5]={
{0.0349, 0.1201, 0.197, 0.0573, 0.1195},
{0.0317, 0.0838, 0.1526, 0.1334, 0.0164},
{0.075, 0.2249, 0.5548, 1.4954, 0.9354},
{0.078, 0.221, 0.674, 1.3867, 0.6925},
{0.0909, 0.2551, 0.7738, 1.2136, 0.4606},
{0.0893, 0.2563, 0.757, 1.0487, 0.3575},
{0.1022, 0.3219, 0.7982, 0.8197, 0.1715},
{0.0974, 0.2921, 0.691, 0.699, 0.2039},
{0.1083, 0.3175, 0.6487, 0.5846, 0.1421},
{0.1269, 0.3535, 0.5582, 0.4674, 0.146},
{0.2142, 0.6853, 0.7692, 1.6589, 1.4482},
{0.2314, 0.6866, 0.9677, 2.1882, 1.1339},
{0.239, 0.6573, 1.2011, 2.5586, 1.2312},
{0.2519, 0.6372, 1.3795, 2.5082, 1.05},
{0.2548, 0.6106, 1.4541, 2.3204, 0.8477},
{0.2497, 0.5628, 1.3899, 2.1865, 0.7715},
{0.2443, 0.5397, 1.3919, 2.0197, 0.6621},
{0.2385, 0.5017, 1.3428, 1.8899, 0.6079},
{0.4115, 1.4031, 2.2784, 2.6742, 2.2162},
{0.4054, 1.388, 2.1602, 3.7532, 2.2063},
{0.3787, 1.2181, 2.0594, 3.2618, 2.387},
{0.3825, 1.2598, 2.0008, 3.0617, 2.0694},
{0.3876, 1.275, 1.9109, 2.8314, 1.8979},
{0.4046, 1.3696, 1.8941, 2.08, 1.2196},
{0.3796, 1.2094, 1.7815, 2.542, 1.5937},
{0.3946, 1.2725, 1.7031, 2.314, 1.4795},
{0.4118, 1.3161, 1.6493, 2.193, 1.283},
{0.386, 1.1765, 1.5451, 2.073, 1.3814},
{0.4314, 1.3208, 1.5236, 1.4671, 0.8562},
{0.4288, 1.2646, 1.4472, 1.8294, 1.0934},
{0.4818, 1.4032, 1.6561, 2.4605, 1.1054},
{0.4655, 1.3014, 1.6088, 2.6998, 1.3003},
{0.4517, 1.2229, 1.5852, 2.7958, 1.2638},
{0.4477, 1.1678, 1.5843, 2.8087, 1.1956},
{0.4798, 1.1948, 1.8695, 2.6953, 0.8203},
{0.4546, 1.0993, 1.7696, 2.7068, 0.8672},
{1.016, 2.8528, 3.5466, -7.7804, 12.1148},
{0.6703, 1.4926, 3.3368, 4.46, 3.1501},
{0.6894, 1.5474, 3.245, 4.2126, 2.9764},
{0.6719, 1.4684, 3.1668, 3.9557, 2.892},
{0.6123, 1.2677, 3.0348, 3.3841, 2.3683},
{0.6773, 1.4798, 3.1788, 3.0824, 1.8384},
{0.7082, 1.6392, 3.1993, 3.4327, 1.8711},
{0.6735, 1.4934, 3.0966, 2.7254, 1.5597},
{0.6413, 1.369, 2.9854, 2.6952, 1.5433},
{0.5904, 1.1775, 2.6519, 2.2875, 0.8689},
{0.6377, 1.379, 2.8294, 2.3631, 1.4553},
{0.6364, 1.4247, 2.7802, 2.5973, 1.7886},
{0.6768, 1.6589, 2.774, 3.1835, 2.1326},
{0.7224, 1.961, 2.7161, 3.5603, 1.8972},
{0.7106, 1.9247, 2.6149, 3.8322, 1.8899},
{0.6947, 1.869, 2.5356, 4.0013, 1.8955},
{0.7047, 1.9484, 2.594, 4.1526, 1.5057},
{0.6737, 1.7908, 2.4129, 4.21, 1.7058},
{1.2704, 3.8018, 5.6618, 0.9205, 4.8105},
{0.9049, 2.6076, 4.8498, 5.1603, 4.7388},
{0.8405, 2.3863, 4.6139, 5.1514, 4.7949},
{0.8551, 2.3915, 4.5772, 5.0278, 4.5118},
{0.9096, 2.5313, 4.5266, 4.6376, 4.369},
{0.8807, 2.4183, 4.4448, 4.6858, 4.1725},
{0.9471, 2.5463, 4.3523, 4.4789, 3.908},
{0.9699, 2.5837, 4.2778, 4.4575, 3.5985},
{0.8694, 2.2413, 3.9196, 3.9694, 4.5498},
{0.9673, 2.4702, 4.1148, 4.4972, 3.2099},
{0.9325, 2.3673, 3.8791, 3.9674, 3.7996},
{0.9505, 2.3705, 3.8218, 4.0471, 3.4451},
{0.9248, 2.2428, 3.6182, 3.791, 3.7912},
{1.0373, 2.4824, 3.6558, 3.8925, 3.0056},
{1.0075, 2.3787, 3.544, 3.6932, 3.1759},
{1.0347, 2.3911, 3.4619, 3.6556, 3.0052},
{0.9927, 2.2436, 3.3554, 3.7813, 3.0994},
{1.0295, 2.2911, 3.411, 3.9497, 2.4925},
{1.019, 2.2291, 3.4097, 3.9252, 2.2679},
{0.9853, 2.1167, 3.357, 3.7981, 2.2798},
{0.9914, 2.0858, 3.4531, 3.8812, 1.8526},
{0.9813, 2.0322, 3.3665, 3.6235, 1.9741},
{1.0194, 2.0645, 3.4425, 3.4914, 1.6976},
{0.9148, 1.8096, 3.2134, 3.2953, 1.5754},
{0.9674, 1.8916, 3.3993, 3.0524, 1.2607},
{1.0033, 1.9469, 3.4396, 3.1548, 1.418},
{1.0689, 2.1038, 3.6039, 3.4927, 1.8283},
{1.0891, 2.1867, 3.616, 3.8031, 1.8994},
{1.1007, 2.2306, 3.5689, 4.1549, 2.0382},
{1.1568, 2.4353, 3.6459, 4.4064, 1.7179},
{1.0909, 2.1976, 3.3831, 4.67, 2.1277},
{1.0756, 2.163, 3.3178, 4.8852, 2.0489},
{1.4282, 3.5081, 5.6767, 4.1964, 3.8946},
{1.3127, 3.1243, 5.2988, 5.3891, 5.4133},
{1.3128, 3.1021, 5.3385, 5.9611, 4.7562},
{1.2553, 2.9178, 5.0862, 6.1206, 4.7122},
{1.3218, 3.1444, 5.4371, 5.6444, 4.0107},
{1.3382, 3.2043, 5.4558, 5.4839, 3.6342},
{1.5193, 4.0053, 6.5327, -0.1402, 6.7489},
{1.3517, 3.2937, 5.3213, 4.6466, 3.5714},
{1.2135, 2.7962, 4.7545, 4.5731, 4.4786},
{1.2937, 3.11, 5.0393, 4.7546, 3.5031},
{1.2915, 3.1023, 4.9309, 4.6009, 3.4661},
{1.2089, 2.7391, 4.3482, 4.0047, 4.6497}
};
const double E_GP_B1[E_TYPE_NUMBER][5]={
{0.5347, 3.5867, 12.3471, 18.9525, 38.6269},
{0.2507, 1.4751, 4.4938, 12.6646, 31.1653},
{0.3864, 2.9383, 15.3829, 53.5545, 138.7337},
{0.3131, 2.2381, 10.1517, 30.9061, 78.3273},
{0.2995, 2.1155, 8.3816, 24.1292, 63.1314},
{0.2465, 1.71, 6.4094, 18.6113, 50.2523},
{0.2451, 1.7481, 6.1925, 17.3894, 48.1431},
{0.2067, 1.3815, 4.6943, 12.7105, 32.4726},
{0.2057, 1.3439, 4.2788, 11.3932, 28.7881},
{0.22, 1.3779, 4.0203, 9.4934, 23.1278},
{0.3334, 2.3446, 10.083, 48.3037, 138.27},
{0.3278, 2.272, 10.9241, 39.2898, 101.9748},
{0.3138, 2.1063, 10.4163, 34.4552, 98.5344},
{0.3075, 2.0174, 9.6746, 29.3744, 80.4732},
{0.2908, 1.874, 8.5176, 24.3434, 63.2996},
{0.2681, 1.6711, 7.0267, 19.5377, 50.3888},
{0.2468, 1.5242, 6.1537, 16.6687, 42.3086},
{0.2289, 1.3694, 5.2561, 14.0928, 35.5361},
{0.3703, 3.3874, 13.1029, 68.9592, 194.4329},
{0.3499, 3.0991, 11.9608, 53.9353, 142.3892},
{0.3133, 2.5856, 9.5813, 41.7688, 116.7282},
{0.304, 2.4863, 9.2783, 39.0751, 109.4583},
{0.2967, 2.378, 8.7981, 35.9528, 101.7201},
{0.2986, 2.3958, 9.1406, 37.4701, 113.7121},
{0.2699, 2.0455, 7.4726, 31.0604, 91.5622},
{0.2717, 2.0443, 7.6007, 29.9714, 86.2265},
{0.2742, 2.0372, 7.7205, 29.968, 84.9383},
{0.2478, 1.766, 6.3107, 25.2204, 74.3146},
{0.2694, 1.9223, 7.3474, 28.9892, 90.6246},
{0.2593, 1.7998, 6.75, 25.586, 73.5284},
{0.2825, 1.9785, 8.7546, 32.5238, 98.5523},
{0.2647, 1.7926, 7.6071, 26.5541, 77.5238},
{0.2493, 1.6436, 6.8154, 22.3681, 62.039},
{0.2405, 1.5442, 6.3231, 19.461, 52.0233},
{0.2504, 1.5963, 6.9653, 19.8492, 50.3233},
{0.2309, 1.4279, 5.9449, 16.6752, 42.2243},
{0.4853, 5.0925, 25.7851, 130.451, 138.6775},
{0.319, 2.2287, 10.3504, 52.3291, 151.2216},
{0.3189, 2.2904, 10.0062, 44.0771, 125.012},
{0.3036, 2.1249, 8.9236, 36.8458, 108.2049},
{0.2709, 1.7683, 7.2489, 27.9465, 98.5624},
{0.292, 2.0606, 8.1129, 30.5336, 100.0658},
{0.2976, 2.2106, 8.5246, 33.1456, 96.6377},
{0.2773, 1.9716, 7.3249, 26.6891, 90.5581},
{0.258, 1.7721, 6.3854, 23.2549, 85.1517},
{0.2324, 1.5019, 5.1591, 15.5428, 46.8213},
{0.2466, 1.6974, 5.7656, 20.0943, 76.7372},
{0.2407, 1.6823, 5.6588, 20.7219, 69.1109},
{0.2522, 1.8545, 6.2936, 25.1457, 84.5448},
{0.2651, 2.0604, 7.3011, 27.5493, 81.3349},
{0.2562, 1.9646, 6.8852, 24.7648, 68.9168},
{0.2459, 1.8542, 6.4411, 22.173, 59.2206},
{0.2455, 1.8638, 6.7639, 21.8007, 56.4395},
{0.2305, 1.689, 5.8218, 18.3928, 47.2496},
{0.4356, 4.2058, 23.4342, 136.778, 171.7561},
{0.3066, 2.4363, 12.1821, 54.6135, 161.9978},
{0.2791, 2.141, 10.34, 41.9148, 132.0204},
{0.2805, 2.12, 10.1808, 42.0633, 130.9893},
{0.2939, 2.2471, 10.8266, 48.8842, 147.602},
{0.2802, 2.0836, 10.0357, 47.4506, 146.9976},
{0.2977, 2.2276, 10.5762, 49.3619, 145.358},
{0.3003, 2.2447, 10.6487, 50.7994, 146.4179},
{0.2653, 1.859, 8.3998, 36.7397, 125.7089},
{0.2909, 2.1014, 9.7067, 43.427, 125.9474},
{0.2761, 1.9511, 8.9296, 41.5937, 131.0122},
{0.2773, 1.9469, 8.8862, 43.0938, 133.1396},
{0.266, 1.8183, 7.9655, 33.1129, 101.8139},
{0.2944, 2.0797, 9.4156, 45.8056, 132.772},
{0.2816, 1.9486, 8.7162, 41.842, 125.032},
{0.2855, 1.9679, 8.7619, 42.3304, 125.6499},
{0.2701, 1.8073, 7.8112, 34.4849, 103.3526},
{0.2761, 1.8625, 8.0961, 34.2712, 98.5295},
{0.2694, 1.7962, 7.6944, 31.0942, 91.1089},
{0.2569, 1.6745, 7.0098, 26.9234, 81.391},
{0.2548, 1.6518, 6.8845, 26.7234, 81.7215},
{0.2487, 1.5973, 6.4737, 23.2817, 70.9254},
{0.2554, 1.6475, 6.5966, 23.2269, 70.0272},
{0.2263, 1.3813, 5.3243, 17.5987, 60.0171},
{0.2358, 1.4712, 5.6758, 18.7119, 61.5286},
{0.2413, 1.5298, 5.8009, 19.452, 60.5753},
{0.254, 1.6715, 6.3509, 23.1531, 78.7099},
{0.2552, 1.7174, 6.5131, 23.917, 74.7039},
{0.2546, 1.7351, 6.4948, 23.6464, 70.378},
{0.2648, 1.8786, 7.1749, 25.1766, 69.2821},
{0.2466, 1.6707, 6.0197, 20.7657, 57.2663},
{0.2402, 1.6169, 5.7644, 19.4568, 52.5009},
{0.3183, 2.6889, 13.4816, 54.3866, 200.8321},
{0.2887, 2.2897, 10.8276, 43.5389, 145.6109},
{0.2861, 2.2509, 10.5287, 41.7796, 128.2973},
{0.2701, 2.0636, 9.3051, 34.5977, 107.92},
{0.2827, 2.225, 10.2454, 41.1162, 124.4449},
{0.2838, 2.2452, 10.2519, 41.7251, 124.9023},
{0.3213, 2.8206, 14.8878, 68.9103, 81.7257},
{0.2813, 2.2418, 9.9952, 42.7939, 132.1739},
{0.2483, 1.8437, 7.5421, 29.3841, 112.4579},
{0.2638, 2.0341, 8.7101, 35.2992, 109.4972},
{0.2611, 2.0023, 8.4377, 34.1559, 105.8911},
{0.2421, 1.7487, 6.7262, 23.2153, 80.3108}
};
//2<S<6
const double E_GP_A2[E_TYPE_NUMBER][5]={
{0.0088, 0.0449, 0.1481, 0.2356, 0.0914},
{0.0084, 0.0443, 0.1314, 0.1671, 0.0666},
{0.0478, 0.2048, 0.5253, 1.5225, 0.9853},
{0.0423, 0.1874, 0.6019, 1.4311, 0.7891},
{0.0436, 0.1898, 0.6788, 1.3273, 0.5544},
{0.0489, 0.2091, 0.7537, 1.142, 0.3555},
{0.0267, 0.1328, 0.5301, 1.102, 0.4215},
{0.0365, 0.1729, 0.5805, 0.8814, 0.3121},
{0.0382, 0.1822, 0.5972, 0.7707, 0.213},
{0.038, 0.1785, 0.5494, 0.6942, 0.1918},
{0.126, 0.6442, 0.8893, 1.8197, 1.2988},
{0.113, 0.5575, 0.9046, 2.158, 1.4735},
{0.1165, 0.5504, 1.0179, 2.6295, 1.5711},
{0.0567, 0.3365, 0.8104, 2.496, 2.1186},
{0.1005, 0.4615, 1.0663, 2.5854, 1.2725},
{0.0915, 0.4312, 1.0847, 2.4671, 1.0852},
{0.0799, 0.3891, 1.0037, 2.3332, 1.0507},
{0.1044, 0.4551, 1.4232, 2.1533, 0.4459},
{0.2149, 0.8703, 2.4999, 2.3591, 3.0318},
{0.2355, 0.9916, 2.3959, 3.7252, 2.5647},
{0.4636, 2.0802, 2.9003, 1.4193, 2.4323},
{0.2123, 0.896, 2.1765, 3.0436, 2.4439},
{0.2369, 1.0774, 2.1894, 3.0825, 1.719},
{0.197, 0.8228, 2.02, 2.1717, 1.7516},
{0.1943, 0.819, 1.9296, 2.4968, 2.0625},
{0.1929, 0.8239, 1.8689, 2.3694, 1.906},
{0.2186, 0.9861, 1.854, 2.3258, 1.4685},
{0.2313, 1.0657, 1.8229, 2.2609, 1.1883},
{0.3501, 1.6558, 1.9582, 0.2134, 1.4109},
{0.178, 0.8096, 1.6744, 1.9499, 1.4495},
{0.2135, 0.9768, 1.6669, 2.5662, 1.679},
{0.2135, 0.9761, 1.6555, 2.8938, 1.6356},
{0.2059, 0.9518, 1.6372, 3.049, 1.4756},
{0.1574, 0.7614, 1.4834, 3.0016, 1.7978},
{0.1899, 0.8983, 1.6358, 3.1845, 1.1518},
{0.1742, 0.8447, 1.5944, 3.1507, 1.1338},
{0.3781, 1.4904, 3.5753, 3.0031, 3.3272},
{0.3723, 1.4598, 3.5124, 4.4612, 3.3031},
{0.3234, 1.2737, 3.2115, 4.0563, 3.7962},
{0.2997, 1.1879, 3.1075, 3.974, 3.5769},
{0.168, 0.937, 2.73, 3.815, 3.0053},
{0.3069, 1.1714, 3.2293, 3.4254, 2.1224},
{0.2928, 1.1267, 3.1675, 3.6619, 2.5942},
{0.2604, 1.0442, 3.0761, 3.2175, 1.9448},
{0.2713, 1.0556, 3.1416, 3.0451, 1.7179},
{0.2003, 0.8779, 2.6135, 2.8594, 1.0258},
{0.2739, 1.0503, 3.1564, 2.7543, 1.4328},
{0.3072, 1.1303, 3.2046, 2.9329, 1.656},
{0.3564, 1.3011, 3.2424, 3.4839, 2.0459},
{0.2966, 1.1157, 3.0973, 3.8156, 2.5281},
{0.2725, 1.0651, 2.994, 4.0697, 2.5682},
{0.2422, 0.9692, 2.8114, 4.1509, 2.8161},
{0.2617, 1.0325, 2.8097, 4.4809, 2.319},
{0.2334, 0.9496, 2.6381, 4.468, 2.502},
{0.5713, 2.4866, 4.9795, 4.0198, 4.4403},
{0.5229, 2.2874, 4.7243, 5.0807, 5.6389},
{0.5461, 2.3856, 5.0653, 5.7601, 4.0463},
{0.2227, 1.076, 2.9482, 5.8496, 7.1834},
{0.5237, 2.2913, 4.6161, 4.7233, 4.8173},
{0.5368, 2.3301, 4.6058, 4.6621, 4.4622},
{0.5232, 2.2627, 4.4552, 4.4787, 4.5073},
{0.5162, 2.2302, 4.3449, 4.3598, 4.4292},
{0.5272, 2.2844, 4.3361, 4.3178, 4.0908},
{0.9664, 3.4052, 5.0803, 1.4991, 4.2528},
{0.511, 2.157, 4.0308, 3.9936, 4.2466},
{0.4974, 2.1097, 3.8906, 3.81, 4.3084},
{0.4679, 1.9693, 3.7191, 3.9632, 4.2432},
{0.5034, 2.1088, 3.8232, 3.7299, 3.8963},
{0.4839, 2.0262, 3.6851, 3.5874, 4.0037},
{0.5221, 2.1695, 3.7567, 3.6685, 3.4274},
{0.468, 1.9466, 3.5428, 3.849, 3.6594},
{0.4048, 1.737, 3.3399, 3.9448, 3.7293},
{0.3835, 1.6747, 3.2986, 4.0462, 3.4303},
{0.3661, 1.6191, 3.2455, 4.0856, 3.2064},
{0.3933, 1.6973, 3.4202, 4.1274, 2.6158},
{0.3854, 1.6555, 3.4129, 4.1111, 2.4106},
{0.351, 1.562, 3.2946, 4.0615, 2.4382},
{0.3083, 1.4158, 2.9662, 3.9349, 2.1709},
{0.3055, 1.3945, 2.9617, 3.899, 2.0026},
{0.3593, 1.5736, 3.5237, 3.8109, 1.6953},
{0.3511, 1.5489, 3.5676, 4.09, 2.5251},
{0.354, 1.5453, 3.5975, 4.3152, 2.7743},
{0.353, 1.5258, 3.5815, 4.5532, 3.0714},
{0.3673, 1.5772, 3.7079, 4.8582, 2.844},
{0.3547, 1.5206, 3.5621, 5.0184, 3.0075},
{0.4586, 1.7781, 3.9877, 5.7273, 1.546},
{0.8282, 2.9941, 5.6597, 4.9292, 4.2889},
{1.4129, 4.4269, 7.046, -1.0573, 8.643},
{0.7169, 2.571, 5.1791, 6.3484, 5.6474},
{0.6958, 2.4936, 5.1269, 6.6988, 5.0799},
{1.2502, 4.2284, 7.0489, 1.139, 5.8222},
{0.641, 2.2643, 4.8713, 5.9287, 5.3935},
{0.6938, 2.4652, 5.1227, 5.5965, 4.8543},
{0.6902, 2.4509, 5.1284, 5.0339, 4.8575},
{0.7577, 2.7264, 5.4184, 4.8198, 4.1013},
{0.7567, 2.7565, 5.4364, 5.1918, 3.5643},
{0.7492, 2.7267, 5.3521, 5.0369, 3.5321},
{0.81, 3.0001, 5.4635, 4.1756, 3.5066}
};
const double E_GP_B2[E_TYPE_NUMBER][5]={
{0.102, 1.0219, 4.6275, 22.8742, 80.1535},
{0.0989, 0.9845, 4.5527, 21.5563, 70.3903},
{0.0926, 0.9182, 4.3291, 19.2996, 58.9329},
{0.0686, 0.6808, 3.1163, 14.3458, 44.0455},
{0.081, 0.7957, 3.9054, 15.7701, 45.6124},
{0.0723, 0.7123, 3.5192, 13.7724, 39.1148},
{0.1557, 1.5347, 9.9947, 51.4251, 185.9828},
{0.148, 1.4643, 9.232, 49.8807, 148.0937},
{0.1244, 1.1948, 7.2756, 34.143, 111.2079},
{0.1121, 1.0638, 6.3891, 28.7081, 97.4289},
{0.0597, 0.6524, 4.4317, 19.554, 85.5011},
{0.1101, 1.0222, 5.9613, 25.1965, 93.5831},
{0.102, 0.9481, 5.4713, 23.8153, 82.8991},
{0.0887, 0.824, 4.8278, 19.8977, 80.4566},
{0.0907, 0.8324, 4.7702, 19.7862, 80.254},
{0.0659, 0.6111, 3.5563, 12.7638, 44.4283},
{0.0881, 0.8028, 4.4451, 18.7011, 79.2633},
{0.0966, 0.8856, 4.6273, 20.6789, 73.4723},
{0.1091, 1.0452, 5.09, 24.6578, 88.0513},
{0.0896, 0.8268, 4.2242, 20.69, 71.3399},
{0.0809, 0.7488, 3.871, 18.88, 60.6499},
{0.0708, 0.6472, 3.3609, 16.0752, 50.1724},
{0.0749, 0.6914, 3.4634, 16.3603, 48.2522},
{0.0655, 0.605, 3.0389, 14.0809, 41.0005},
{0.1626, 1.8213, 11.1049, 49.0568, 202.9987},
{0.1434, 1.6019, 9.4511, 42.7685, 148.4969},
{0.1479, 1.6552, 10.0059, 47.3245, 145.8464},
{0.0571, 0.5946, 3.2022, 16.4253, 95.703},
{0.136, 1.5068, 8.8213, 41.9536, 141.2424},
{0.1378, 1.514, 8.8719, 43.5967, 141.8065},
{0.1317, 1.4336, 8.3087, 40.601, 135.9196},
{0.1279, 1.3811, 7.9629, 39.1213, 132.7846},
{0.1285, 1.3943, 8.1081, 40.9631, 134.1233},
{0.2641, 2.6586, 16.2213, 80.206, 92.5359},
{0.121, 1.2704, 7.1368, 35.0354, 123.5062},
{0.1157, 1.2108, 6.7377, 32.415, 116.9225},
{0.1069, 1.0994, 5.9769, 27.1491, 96.3119},
{0.1141, 1.1769, 6.6087, 33.4332, 116.4913},
{0.1081, 1.1012, 6.1114, 30.3728, 110.5988},
{0.1148, 1.186, 6.752, 35.6807, 118.0692},
{0.1015, 1.0195, 5.6058, 27.4899, 95.2846},
{0.0868, 0.8585, 4.6378, 21.69, 80.2408},
{0.081, 0.802, 4.3545, 19.9644, 73.6337},
{0.0761, 0.7543, 4.0952, 18.2886, 68.0967},
{0.0806, 0.7972, 4.4237, 19.5692, 68.7477},
{0.0787, 0.7638, 4.2441, 18.37, 65.1071},
{0.0706, 0.6904, 3.8266, 16.0812, 58.7638},
{0.0609, 0.5993, 3.1921, 12.5285, 49.7675},
{0.0596, 0.5827, 3.1035, 11.9693, 47.9106},
{0.0694, 0.6758, 3.8457, 15.6203, 56.6614},
{0.0672, 0.6522, 3.742, 15.9791, 65.1354},
{0.0668, 0.6465, 3.6968, 16.2056, 61.4909},
{0.0661, 0.6324, 3.5906, 15.9962, 57.576},
{0.0678, 0.6527, 3.7396, 17.0668, 55.9789},
{0.0649, 0.6188, 3.4696, 15.609, 49.4818},
{0.0831, 0.784, 4.3599, 20.0128, 62.1535},
{0.1515, 1.6163, 9.7752, 42.848, 190.7366},
{0.2921, 3.1381, 19.6767, 102.043, 113.9798},
{0.1263, 1.29, 7.3686, 32.449, 118.0558},
{0.1211, 1.2247, 6.9398, 30.0991, 105.196},
{0.2415, 2.6442, 16.3313, 73.5757, 91.9401},
{0.1097, 1.0644, 5.7907, 25.0261, 101.3899},
{0.1171, 1.1757, 6.4053, 27.5217, 103.0482},
{0.1153, 1.1545, 6.2291, 27.0741, 111.315},
{0.1257, 1.3044, 7.1035, 32.4649, 118.8647},
{0.1239, 1.2979, 7.0798, 32.7871, 110.1512},
{0.1217, 1.2651, 6.8101, 31.6088, 106.4853},
{0.131, 1.4038, 7.6057, 34.0186, 90.5226}
};

#endif
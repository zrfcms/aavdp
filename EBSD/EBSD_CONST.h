#ifndef __EBSD_CONST_H__
#define __EBSD_CONST_H__
#define MAX_LENTH_IN_NAME 300
#define MAX_WORD_NUMBER 1000

#define MAX_MIN 1.0e8

#define PI 3.141592653589793
#define TWO_PI 6.283185307179586
#define FOUR_PI 12.566370614359172
#define PI_SQRT_HALF 0.886226925452758 //sqrt(PI)/2
#define PI_HALF_SQRT 1.253314137315500 //sqrt(PI/2)
#define PI_INVERSE 0.318309886183791 //1/PI
#define PI_PREG 1.55512030155621410 //2*sqrt(PI)/3^(3/4)
#define PI_PREC 0.9068996821171090 //PI/2*sqrt(3)
#define PI_PREB 1.0500751358086640 //3^(1/4)*sqrt(2/PI)
#define PI_PRED 2.0943951023931950 //2*PI/3
#define PI_PREA 0.5250375679043320 //3^(1/4)/sqrt(2*PI)
#define PI_PREE 0.7598356856515930 //3^(-1/4)
#define PI_PREF 0.7598356856515930 //3^(-1/4)
#define DEG_TO_RAD 0.017453292519943295 //PI/180
#define SQRT_3 1.732050807568877 //sqrt(3)
#define SQRT_3_INVERSE 0.577350269190 //1/sqrt(3)
#define SQRT_HALF_3 0.866025403780 //sqrt(3)/2

#define SCATTERING_EVENT_NUMBER 500 //Number of scattering events along one trajectory
#define CEN_TO_ANG 1.0E8
#define ZERO 1e-12
#define ONE 0.99999
#define CONST_IN_AREA_INV 1.5273987E19 //1.0/(5.21E-21*(4*PI))

#define G_TO_WK 0.6283185307179586 //0.1*TWO_PI
#define DW_TO_WK 1.2665147955292222 //100.0/(8.0*pow(PI, 2))
#define PRE_CONST_V 0.04787801
#define PRE_CONST_U 0.664840340614319 //2.0*m0*e/h**2*1.0E-18
#define PRE_CONST_RI1 0.5772157
#define ELECTRON_CHARGE 1.602176634e-19 //Coulomb
#define ELECTRON_REST_MASS 9.1093837090e-31 //kg
#define ELECTRON_REST_ENERGY 0.511e6 //eV
#define LIGHT_VELOCITY 299792458.0 //m/s
#define AVOGADRO_CONSTANT 6.02214076e23
#define PLANK_CONSTANT 6.62607015e-34 //J·s

#define SPACE_GROUP_NUMBER 237
#define SYMMORPHIC_SPACE_GROUP_NUMBER 73
#define POINT_SYMMETRY_NUMBER 48 //number of point group symmetry elements
#define POINT_GROUP_NUMBER 32 //number of point group
#define SAMPLING_TYPE_NUMBER 36
#define MAX_MULTIPLICITY 192
#define TYPE_NUMBER 98

const char SPACE_GROUP_SYMBOL[SPACE_GROUP_NUMBER][11]={
"P1","P-1",
"P2","P21","C2","Pm",
"Pc","Cm","Cc","P2/m",
"P21/m","C2/m","P2/c","P21/c",
"C2/c",
"P222","P2221","P21212","P212121",
"C2221","C222","F222","I222",
"I212121","Pmm2","Pmc21","Pcc2",
"Pma2","Pca21","Pnc2","Pmn21",
"Pba2","Pna21","Pnn2","Cmm2",
"Cmc21","Ccc2","Amm2","Abm2",
"Ama2","Aba2","Fmm2","Fdd2",
"Imm2","Iba2","Ima2","Pmmm",
"Pnnn","Pccm","Pban","Pmma",
"Pnna","Pmna","Pcca","Pbam",
"Pccn","Pbcm","Pnnm","Pmmn",
"Pbcn","Pbca","Pnma","Cmcm",
"Cmca","Cmmm","Cccm","Cmma",
"Ccca","Fmmm","Fddd","Immm",
"Ibam","Ibca","Imma",
"P4","P41","P42","P43",
"I4","I41","P-4","I-4",
"P4/m","P42/m","P4/n","P42/n",
"I4/m","I41/a","P422","P4212",
"P4122","P41212","P4222","P42212",
"P4322","P43212","I422","I4122",
"P4mm","P4bm","P42cm","P42nm",
"P4cc","P4nc","P42mc","P42bc",
"I4mm","I4cm","I41md","I41cd",
"P-42m","P-42c","P-421m","P-421c",
"P-4m2","P-4c2","P-4b2","P-4n2",
"I-4m2","I-4c2","I-42m","I-42d",
"P4/mmm","P4/mcc","P4/nbm","P4/nnc",
"P4/mbm","P4/mnc","P4/nmm","P4/ncc",
"P42/mmc","P42/mcm","P42/nbc","P42/nnm",
"P42/mbc","P42/mnm","P42/nmc","P42/ncm",
"I4/mmm","I4/mcm","I41/amd","I41/acd",
"P3","P31","P32","R3",
"P-3","R-3","P312","P321",
"P3112","P3121","P3212","P3221",
"R32","P3m1","P31m","P3c1",
"P31c","R3m","R3c","P-31m",
"P-31c","P-3m1","P-3c1","R-3m",
"R-3c",
"P6","P61","P65","P62",
"P64","P63","P-6","P6/m",
"P63/m","P622","P6122","P6522",
"P6222","P6422","P6322","P6mm",
"P6cc","P63cm","P63mc","P-6m2",
"P-6c2","P-62m","P-62c","P6/mmm",
"P6/mcc","P63/mcm","P63/mmc",
"P23","F23","I23","P213",
"I213","Pm3","Pn3","Fm3",
"Fd3","Im3","Pa3","Ia3",
"P432","P4232","F432","F4132",
"I432","P4332","P4132","I4132",
"P-43m","F-43m","I-43m","P-43n",
"F-43c","I-43d","Pm3m","Pn3n",
"Pm3n","Pn3m","Fm3m","Fm3c",
"Fd3m","Fd3c","Im3m","Ia3d",
"R3|146","R-3|148","R32|155","R3m|160",
"R3c|161","R-3m|166","R-3c|167"};
const char GENERATOR_STRING[SPACE_GROUP_NUMBER][40]={
"000","100","01cOOO0",
"01cODO0","02aDDOcOOO0","01jOOO0",
"01jOOD0","02aDDOjOOO0","02aDDOjOOD0",
"11cOOO0","11cODO0","12aDDOcOOO0",
"11cOOD0","11cODD0","12aDDOcOOD0",
"02bOOOcOOO0","02bOODcOOD0","02bOOOcDDO0",
"02bDODcODD0","03aDDObOODcOOD0","03aDDObOOOcOOO0",
"04aODDaDODbOOOcOOO0","03aDDDbOOOcOOO0","03aDDDbDODcODD0",
"02bOOOjOOO0","02bOODjOOD0","02bOOOjOOD0",
"02bOOOjDOO0","02bOODjDOO0","02bOOOjODD0",
"02bDODjDOD0","02bOOOjDDO0","02bOODjDDO0",
"02bOOOjDDD0","03aDDObOOOjOOO0","03aDDObOODjOOD0",
"03aDDObOOOjOOD0","03aODDbOOOjOOO0","03aODDbOOOjODO0",
"03aODDbOOOjDOO0","03aODDbOOOjDDO0","04aODDaDODbOOOjOOO0",
"04aODDaDODbOOOjBBB0","03aDDDbOOOjOOO0","03aDDDbOOOjDDO0",
"03aDDDbOOOjDOO0","12bOOOcOOO0","03bOOOcOOOhDDD1BBB",
"12bOOOcOOD0","03bOOOcOOOhDDO1BBO","12bDOOcOOO0",
"12bDOOcDDD0","12bDODcDOD0","12bDOOcOOD0",
"12bOOOcDDO0","12bDDOcODD0","12bOODcODD0",
"12bOOOcDDD0","03bOOOcDDOhDDO1BBO","12bDDDcOOD0",
"12bDODcODD0","12bDODcODO0","13aDDObOODcOOD0",
"13aDDObODDcODD0","13aDDObOOOcOOO0","13aDDObOOOcOOD0",
"13aDDObODOcODO0","04aDDObDDOcOOOhODD1OBB","14aODDaDODbOOOcOOO0",
"05aODDaDODbOOOcOOOhBBB1ZZZ","13aDDDbOOOcOOO0","13aDDDbOOOcDDO0",
"13aDDDbDODcODD0","13aDDDbODOcODO0","02bOOOgOOO0",
"02bOODgOOB0","02bOOOgOOD0","02bOODgOOF0",
"03aDDDbOOOgOOO0","03aDDDbDDDgODB0","02bOOOmOOO0",
"03aDDDbOOOmOOO0","12bOOOgOOO0","12bOOOgOOD0",
"03bOOOgDDOhDDO1YBO","03bOOOgDDDhDDD1YYY","13aDDDbOOOgOOO0",
"04aDDDbDDDgODBhODB1OYZ","03bOOOgOOOcOOO0","03bOOOgDDOcDDO0",
"03bOODgOOBcOOO0","03bOODgDDBcDDB0","03bOOOgOODcOOO0",
"03bOOOgDDDcDDD0","03bOODgOOFcOOO0","03bOODgDDFcDDF0",
"04aDDDbOOOgOOOcOOO0","04aDDDbDDDgODBcDOF0","03bOOOgOOOjOOO0",
"03bOOOgOOOjDDO0","03bOOOgOODjOOD0","03bOOOgDDDjDDD0",
"03bOOOgOOOjOOD0","03bOOOgOOOjDDD0","03bOOOgOODjOOO0",
"03bOOOgOODjDDO0","04aDDDbOOOgOOOjOOO0","04aDDDbOOOgOOOjOOD0",
"04aDDDbDDDgODBjOOO0","04aDDDbDDDgODBjOOD0","03bOOOmOOOcOOO0",
"03bOOOmOOOcOOD0","03bOOOmOOOcDDO0","03bOOOmOOOcDDD0",
"03bOOOmOOOjOOO0","03bOOOmOOOjOOD0","03bOOOmOOOjDDO0",
"03bOOOmOOOjDDD0","04aDDDbOOOmOOOjOOO0","04aDDDbOOOmOOOjOOD0",
"04aDDDbOOOmOOOcOOO0","04aDDDbOOOmOOOcDOF0","13bOOOgOOOcOOO0",
"13bOOOgOOOcOOD0","04bOOOgOOOcOOOhDDO1YYO","04bOOOgOOOcOOOhDDD1YYY",
"13bOOOgOOOcDDO0","13bOOOgOOOcDDD0","04bOOOgDDOcDDOhDDO1YBO",
"04bOOOgDDOcDDDhDDO1YBO","13bOOOgOODcOOO0","13bOOOgOODcOOD0",
"04bOOOgDDDcOODhDDD1YBY","04bOOOgDDDcOOOhDDD1YBY","13bOOOgOODcDDO0",
"13bOOOgDDDcDDD0","04bOOOgDDDcDDDhDDD1YBY","04bOOOgDDDcDDOhDDD1YBY",
"14aDDDbOOOgOOOcOOO0","14aDDDbOOOgOOOcOOD0","05aDDDbDDDgODBcDOFhODB1OBZ",
"05aDDDbDDDgODBcDOBhODB1OBZ","01nOOO0","01nOOC0",
"01nOOE0","02aECCnOOO0","11nOOO0",
"12aECCnOOO0","02nOOOfOOO0","02nOOOeOOO0",
"02nOOCfOOE0","02nOOCeOOO0","02nOOEfOOC0",
"02nOOEeOOO0","03aECCnOOOeOOO0","02nOOOkOOO0",
"02nOOOlOOO0","02nOOOkOOD0","02nOOOlOOD0",
"03aECCnOOOkOOO0","03aECCnOOOkOOD0","12nOOOfOOO0",
"12nOOOfOOD0","12nOOOeOOO0","12nOOOeOOD0",
"13aECCnOOOeOOO0","13aECCnOOOeOOD0","02nOOObOOO0",
"02nOOCbOOD0","02nOOEbOOD0","02nOOEbOOO0",
"02nOOCbOOO0","02nOOObOOD0","02nOOOiOOO0",
"12nOOObOOO0","12nOOObOOD0","03nOOObOOOeOOO0",
"03nOOCbOODeOOC0","03nOOEbOODeOOE0","03nOOEbOOOeOOE0",
"03nOOCbOOOeOOC0","03nOOObOODeOOO0","03nOOObOOOkOOO0",
"03nOOObOOOkOOD0","03nOOObOODkOOD0","03nOOObOODkOOO0",
"03nOOOiOOOkOOO0","03nOOOiOODkOOD0","03nOOOiOOOeOOO0",
"03nOOOiOODeOOO0","13nOOObOOOeOOO0","13nOOObOOOeOOD0",
"13nOOObOODeOOD0","13nOOObOODeOOO0","03bOOOcOOOdOOO0",
"05aODDaDODbOOOcOOOdOOO0","04aDDDbOOOcOOOdOOO0","03bDODcODDdOOO0",
"04aDDDbDODcODDdOOO0","13bOOOcOOOdOOO0","04bOOOcOOOdOOOhDDD1YYY",
"15aODDaDODbOOOcOOOdOOO0","06aODDaDODbOOOcOOOdOOOhBBB1ZZZ","14aDDDbOOOcOOOdOOO0",
"13bDODcODDdOOO0","14aDDDbDODcODDdOOO0","04bOOOcOOOdOOOeOOO0",
"04bOOOcOOOdOOOeDDD0","06aODDaDODbOOOcOOOdOOOeOOO0","06aODDaDODbODDcDDOdOOOeFBF0",
"05aDDDbOOOcOOOdOOOeOOO0","04bDODcODDdOOOeBFF0","04bDODcODDdOOOeFBB0",
"05aDDDbDODcODDdOOOeFBB0","04bOOOcOOOdOOOlOOO0","06aODDaDODbOOOcOOOdOOOlOOO0",
"05aDDDbOOOcOOOdOOOlOOO0","04bOOOcOOOdOOOlDDD0","06aODDaDODbOOOcOOOdOOOlDDD0",
"05aDDDbDODcODDdOOOlBBB0","14bOOOcOOOdOOOeOOO0","05bOOOcOOOdOOOeOOOhDDD1YYY",
"14bOOOcOOOdOOOeDDD0","05bOOOcOOOdOOOeDDDhDDD1YYY","16aODDaDODbOOOcOOOdOOOeOOO0",
"16aODDaDODbOOOcOOOdOOOeDDD0","07aODDaDODbODDcDDOdOOOeFBFhBBB1ZZZ","07aODDaDODbODDcDDOdOOOeFBFhFFF1XXX",
"15aDDDbOOOcOOOdOOOeOOO0","15aDDDbDODcODDdOOOeFBB0","01dOOO0",
"11dOOO0","02dOOOfOOO0","02dOOOlOOO0",
"02dOOOlDDD0","12dOOOfOOO0","12dOOOfDDD0"};
const int SYMMORPHIC_SPACE_GROUP[SYMMORPHIC_SPACE_GROUP_NUMBER]={
1,2,3,5,6,8,10,12,16,21,22,23,25,35,38,42,44,47,
65,69,71,75,79,81,82,83,87,89,97,99,107,111,115,
119,121,123,139,143,146,147,148,149,150,155,156,
157,160,162,164,166,168,174,175,177,183,187,189,
191,195,196,197,200,202,204,207,209,211,215,216,
217,221,225,229};
//the first space group for a given point group, to determine the point group for a given space group
const int PG_FIRST_SPACE_GROUP[POINT_GROUP_NUMBER]={
1,2,3,6,10,16,25,47,75,81,83,89,99,111,123,143,
147,149,156,162,168,174,175,177,183,187,191,195,200,207,215,221};
const int PG_SAMPLING_TYPE[SAMPLING_TYPE_NUMBER]={1,2,3,4,5,5,5,6,5,5,
6,6,7,-1,9,-1,-1,-1,-1,-1,
15,12,17,16,18,-1,19,3,6,6,
8,9,-1,-1,-1,-1};

const char ATOM_SYMBOL[TYPE_NUMBER][5]={
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
const int ATOM_NUMBER[TYPE_NUMBER]={
1,2,3,4,5,
6,7,8,9,10,
11,12,13,14,15,
16,17,18,19,20,
21,22,23,24,25,
26,27,28,29,30,
31,32,33,34,35,
36,37,38,39,40,
41,42,43,44,45,
46,47,48,49,50,
51,52,53,54,55,
56,57,58,59,60,
61,62,63,64,65,
66,67,68,69,70,
71,72,73,74,75,
76,77,78,79,80,
81,82,83,84,85,
86,87,88,89,90,
91,92,93,94,95,
96,97,98};
const double ATOM_MASS[TYPE_NUMBER]={
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

const double AA[TYPE_NUMBER][4]={
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
const double BB[TYPE_NUMBER][4]={
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

#define DURCH_NUMBER 21
const double DURCH_TABLE[DURCH_NUMBER]={
1.000000,1.005051,1.010206,1.015472,1.020852,
1.026355,1.031985,1.037751,1.043662,1.049726,
1.055956,1.062364,1.068965,1.075780,1.082830,
1.090140,1.097737,1.105647,1.113894,1.122497,
1.131470};
#endif
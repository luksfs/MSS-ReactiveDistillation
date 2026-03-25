$title Rective ditillation - MSS analysis

$onText
Just to show the application of our methodology in obtaining the best steady state solution
$offText
* ============================================================================ *
* =================================== SETS =================================== *
* ============================================================================ *
SET

 i "Number of components" /1*4/
 j "Number of stages"  /1*30/
 iter "Number of iteration" /1*315/ 
 ;


* ============================================================================ *
* ============================= Discrete Variables ============================= *
* ============================================================================ *
  
SCALAR
  Nc "Total number of components" /4/
  Ns "Total number of stages"
  NR "First reactive tray"
  NFE "Ethanol feed stage"
  NFB "Butenes feed stage";
  
* ============================================================================ *
* ============================= Solution Parameters ============================= *
* ============================================================================ *
ALIAS(j,jj,zz,xx,hh);


PARAMETERS
Z_cost(j,jj,zz,xx,hh)          "Compressibility factor for each solution"
v_mol_cost(j,jj,zz,xx,hh)   "Molar volume for each solution (m3/mol)"
Mw_cost(j,jj,zz,xx,hh)      "Molar weight for each soltution (kg/kmol)"
Qr_cost(j,jj,zz,xx,hh)        "Condenser duty for each solution (J/min)"
Qc_cost(j,jj,zz,xx,hh)        "Reboiler duty for each solution (J/min)"
T_cost(j,jj,zz,xx,hh)          "Temp profile (K)"
L_cost(j,jj,zz,xx,hh)          "Molar flowrate in liquid (mol/min)"
V_cost(j,jj,zz,xx,hh)          "Molar flowrate (mol/min)"
x_column(i,j,jj,zz,xx,hh)                "liquid molar frac"
y_column(i,j,jj,zz,xx,hh)                "vapor molar frac"
Tcond_cost(j,jj,zz,xx,hh)          "Temp condenser (K)"
Treb_cost(j,jj,zz,xx,hh)          "Temp reboiler (K)"
RR_array(j,jj,zz,xx,hh)        "Reflux Ratio"
x_ETBE_array(j,jj,zz,xx,hh)  "ETBE molar fraction at bottom product"
x_ETBE_D_array(j,jj,zz,xx,hh)  "ETBE molar fraction at top product"
Bottom(j,jj,zz,xx,hh)             "Bottom product (mol/min)"
Top(j,jj,zz,xx,hh)             "Top product (mol/min)"
count_s                             "Successful candidates counter"
count_f                             "Failed candidates counter"
count_s_conopt                 "Successful candidates counter for conopt"
Mw_mix_temp(j)
Tcond_t
profit_obj(j,jj,zz,xx,hh)
profit_obj_conopt(j,jj,zz,xx,hh) 'objective function conopt'
X_ETOH_copnopt(j,jj,zz,xx,hh) 'Ethanol conversion  conopt'
x_ETBE_conopt(j,jj,zz,xx,hh)  "ETBE molar fraction at bottom product - conopt"
Breb_cost(j,jj,zz,xx,hh)
model_stat(j,jj,zz,xx,hh)
model_stat_conopt(j,jj,zz,xx,hh)
X_EtOH_array(j,jj,zz,xx,hh)
D_colMax_arr(j,jj,zz,xx,hh)
D_col_arr(j,jj,zz,xx,hh)
results(iter)
it
;

it=0;
count_s = 0;
count_f = 0;
count_s_conopt = 0;

SCALAR R_cost "Gas constant J K-1 kmol-1" /8.314463e3/;
SCALAR Elapsed_time;

PARAMETER X_ETOH(j), x_D_ETBE(j), jmax, test_par(j);
SCALAR RR;

SCALAR
 Tmin /315/
 Tmax /410/
 jmin /1/;
* jmax /10/;
 
* jmax = Ns;

* ============================================================================ *
* ============================= GIVEN PARAMETERS ============================= *
* ============================================================================ *

*SRK  - components propeerties 1- n-butene | CAS RN: 106-98-9 ;  2 - Ethanol | CAS RN: 64-17-5 ;  3 - Isobutene/Isobutylene | CAS RN: 115-11-7 ; 4 - ETBE - Ethyl tert-butyl ether | CAS RN: 637-92-3

PARAMETER
 v_i(i)   "Stoichiometric coefficients" /1 0, 2 -1, 3 -1, 4 1/
 m_cat    "Catalyst mass (kg)" /0.4/  
 FE       "Ethanol molar flowrate (mol/min)"
 zE(i)    "Molar fraction of ethanol" /1 0, 2 1, 3 0, 4 0/
 FB       "Butenes molar flowrate (mol/min)" 
 zB(i)    "Molar fraction of butenes" /1 0.7, 2 0, 3 0.3, 4 0/  
 R        "Universal gas constant (J/mol-K)" /8.314/
 Omega_a  "SRK parameter" /0.42747/
 Omega_b  "SRK parameter" /0.08664/
 Mw(i)    "Molecular weights (kg/kmol)" / 1 56.10752, 2 46.06904, 3 56.10752, 4 102.17656 /
* Mw_mix(j)  "Molecular weights (kg/kmol)"
 Pc(i)    "Critical pressure (Pa)"      / 1 4.02E6, 2 6.137E6, 3 4.0E6, 4 2.934E6 /
 Tc(i)    "Critical temperature (K)"    / 1 419.5, 2 514, 3 417.9, 4 509.4 /
 w(i)     "Acentric factor"             / 1 0.184495, 2 0.643558, 3 0.19484, 4 0.316231 /
 Tb(i)    "Boiling temperature (K)"     / 1 342.6, 2 421.5, 3 341.8, 4 441.4 /
 Hform(i) "Heat of formation (J/mol)"   / 1 -500, 2 -234950, 3 -17100, 4 -313900 /
 kappa(i) "SRK equation parameter"
 ac(i)    "Component-specific parameter"
 bii(i)   "SRK b-parameter"
 F_factor "Scaling factor (mol/min)"
 FH_factor  "Scaling factor for duty (J/min)"
 H_factor "Scaling factor for enthalpy (J/mol)"
 TF_factor "Scaling Temperature (K)"
*Tcond "Condenser temperature (K)"
*Treb "Reboiler Temperature (K)"; 
;
 
* ============================================================================ *
* Enthalpy for two feeds
* ============================================================================ *

PARAMETERS
    HFB_d    "Enthalpy of butenos feed [J/mol]" /-21481.62933928521/
    HFE_d    "Enthalpy of ethanol feed [J/mol]" /-273241.2338487912/
    HFB      "Enthalpy of butenos feed [-]"
    HFE      "Enthalpy of ethanol feed [-]";
SCALAR NFE "Ethanol feed stage", NFB "Butenes feed stage";

*scaling factors
F_factor = (1.7118+5.774);
FH_factor = 1.7118*abs(HFE_d) + 5.774*abs(HFB_d);
H_factor = F_factor/FH_factor;
TF_factor = (342.38 + 325)/2;

HFE = HFE_d*H_factor;
HFB = HFB_d*H_factor;
*Feed
FE = 1.7118/F_factor;
FB = 5.774/F_factor;

* Calculate component-specific parameters
kappa(i) = 0.480 + 1.574*w(i) - 0.176*w(i)**2;
ac(i) = Omega_a * (R * Tc(i)) ** 2 / Pc(i);
bii(i) = Omega_b * R * Tc(i) / Pc(i);

* ============================================================================ *
* ================================= VARIABLES ================================ *
* ============================================================================ *

POSITIVE VARIABLES
                   L(j)        "Liquid molar flowrate (mol/min)"
                   V(j)        "Vapor molar flowrate (mol/min)"
                   x(i,j)      "Liquid molar fraction"
                   y(i,j)      "Vapor molar fraction"
                   T(j)        "Temperature (K)"
                   D           "Destillate product flowrate (mol/min)"
                   Qc          "Condenser Duty (J/min)"
                   Qr          "Reboiler Duty (J/min)";

* ============================================================================ *
* ================================= EQUATIONS ================================ *
* ============================================================================ *

$ontext

* Scalars
Parameter D  / 0.5748006462583984 /;
Parameter Qc / 0.5451343354854316 /;
Parameter Qr / 0.4592022879232921 /;


*POSITIVE VARIABLE D,Qc,Qr;

*D.l =  0.574800646258398;
*Qc.l = 0.5451343354854316;
*Qr.l = 0.4592022879232921;

* Liquid flow rates
Parameter L(j) /
1 1.7244, 2 1.9468, 3 1.7747, 4 1.5311, 5 1.2221, 
6 1.2554, 7 1.9097, 8 1.6147, 9 1.3435, 10 0.2330
/;

* Vapor flow rates
Parameter V(j) /
1 1.9888E-24, 2 2.2992, 3 2.2929, 4 2.1597, 5 1.9758, 
6 1.7604, 7 1.7937, 8 1.6767, 9 1.3817, 10 1.1105
/;


* Liquid compositions
Parameter x(i,j) /
1.1 0.8996, 1.2 0.7902, 1.3 0.7753, 1.4 0.7422, 1.5 0.6762,
1.6 0.6274, 1.7 0.6266, 1.8 0.5229, 1.9 0.2933, 1.10 0.0980,
2.1 0.0444, 2.2 0.1560, 2.3 0.1396, 2.4 0.0967, 2.5 0.0389,
2.6 0.0370, 2.7 0.0365, 2.8 0.0660, 2.9 0.0847, 2.10 0.0471,
3.1 0.0554, 3.2 0.0486, 3.3 0.0476, 3.4 0.0608, 3.5 0.0916,
3.6 0.1574, 3.7 0.2033, 3.8 0.1704, 3.9 0.0957, 3.10 0.0316,
4.1 0.0006, 4.2 0.0052, 4.3 0.0376, 4.4 0.1003, 4.5 0.1933,
4.6 0.1782, 4.7 0.1336, 4.8 0.2408, 4.9 0.5263, 4.10 0.8234
/;

* Vapor compositions
Parameter y(i,j) /
1.1 0.9219, 1.2 0.8996, 1.3 0.8964, 1.4 0.8765, 1.5 0.8369,
1.6 0.7631, 1.7 0.7274, 1.8 0.7000, 1.9 0.5945, 1.10 0.3343,
2.1 0.0209, 2.2 0.0444, 2.3 0.0439, 2.4 0.0386, 2.5 0.0220,
2.6 0.0208, 2.7 0.0198, 2.8 0.0351, 2.9 0.0691, 2.10 0.0926,
3.1 0.0571, 3.2 0.0554, 3.3 0.0551, 3.4 0.0718, 3.5 0.1131,
3.6 0.1909, 3.7 0.2351, 3.8 0.2271, 3.9 0.1938, 3.10 0.1092,
4.1 0.00007, 4.2 0.0006, 4.3 0.0046, 4.4 0.0130, 4.5 0.0279,
4.6 0.0252, 4.7 0.0178, 4.8 0.0377, 4.9 0.1425, 4.10 0.4640
/;


* Temperature profile
Parameter T(j) /
1 334.9684, 2 336.9010, 3 338.0852, 4 340.1857, 5 343.4608,
6 342.6741, 7 340.5410, 8 346.5390, 9 366.6270, 10 398.0755
/;

$offtext

* --- scaling parameter for pressure
Parameter P_scale 'Scaling factor to convert bar to Pa' /1e5/;

EQUATION psat_def(i,j);
set c   'Psat equation coefficients'   /j1*j7/
*  Saturation Pressure expanded Antoine (Aspen data)
POSITIVE VARIABLES logPsat(i,j) "Saturation Pressure (Pa)";

Table kk(i,c) 'Coefficients for the saturation pressure equation'
           j1         j2        j3       j4        j5       j6             j7
    1     51.836   -4019.2      0        0      -4.5229    4.8833E-17      6
    2     73.304   -7122.3      0        0      -7.1424    2.8853E-06      2
    3     78.01    -4634.1      0        0      -8.9575    1.3413E-05      2
    4     64.188   -5820.2      0        0      -6.1343    2.1405E-17      6;

*$macro psat_1 kk(i,'coeff1') + kk(i,'coeff2')/(kk(i,'coeff3') + T(j))
*$macro psat_2 kk(i,'coeff5')*log(T(j))
*$macro psat_3 kk(i,'coeff6')*(T(j)**kk(i,'coeff7'))
*
psat_def(i,j)$(ord(j) <= Ns) ..
     logPsat(i,j) =E= ( kk(i,'j1') + kk(i,'j2')/(kk(i,'j3') + T(j))
               + kk(i,'j4')*T(j) + kk(i,'j5')*log(T(j))
               + kk(i,'j6')*T(j)**kk(i,'j7') );
 
* ============================================================================ *
* Soave Redlich Kwong for Mixtures (phi_mix,HR_mix)
* ============================================================================ *

VARIABLES
    Z(j)          "Compressibility factor"
    v_mol(j)      "Molar volume"
    phi(i,j)      "Fugacity coefficient"
    HR(j)         "Residual enthalpy",
    alpha_ij(i,j) "Alpha parameter"
    aii(i,j)    "Component-wise SRK aii parameter"
    a(j)        "Mixture parameter a"
    b(j) "Mixture parameter b"
    dadT_ii(i,j) "Partial derivative of a to T"
    dadT(j)     "Total temperature derivative of a";

PARAMETER P(j) "Pressure (Pa)";
* Assign the same value to all trays
P(j) = 9.5e5;

EQUATIONS
    cubic_eos,
    v_definition,
    fugacity_eq,
    enthalpy_eq,
    alpha_definition,
    aii_definition,
    mixing_rule_b,
    mixing_rule_a,
    dadT_ii_eq,
    dadT_eq,
    cubic_eos_ideal,
    fugacity_eq_ideal,
    enthalpy_eq_ideal;

$macro Tr(i,j) ( T(j)/Tc(i) )
$macro Betav(j) ( b(j) * P(j) ) / ( R*T(j) )
$macro qv(j) ( a(j) / (b(j)*R*T(j) ) )

* Complete model
cubic_eos(j)$(ord(j) <= Ns)  ..
    -Z(j) + 1 + Betav(j) - qv(j) * Betav(j) * ( Z(j) - Betav(j) ) / ( Z(j)*( Z(j) + Betav(j) ) ) =E= 0;

alpha_definition(i,j)$(ord(j) <= Ns)  ..
    alpha_ij(i,j) =E= ( 1 + kappa(i) * ( 1 - sqrt(Tr(i,j)) )  )**2;

aii_definition(i,j)$(ord(j) <= Ns)  ..
    aii(i,j) =E= alpha_ij(i,j) * ac(i);
    
mixing_rule_b(j)$(ord(j) <= Ns)  ..
    b(j) =E= SUM(i, y(i,j) * bii(i));
    
* Create an alias for set i
ALIAS(i, k);
mixing_rule_a(j)$(ord(j) <= Ns)  .. 
    a(j) =E= SUM(i, SUM(k, y(i,j) * y(k,j) * sqrt(aii(i,j) * aii(k,j))));

v_definition(j)$(ord(j) <= Ns)  ..
    v_mol(j) =E= Z(j)*R*T(j)/P(j);

fugacity_eq(i,j)$(ord(j) <= Ns)  ..
    phi(i,j) =E= EXP(
        bii(i)*(Z(j)-1)/b(j) - log( P(j)*(v_mol(j)-b(j))/(R*T(j)) ) +
        (a(j)/(b(j)*R*T(j))) * ( bii(i)/b(j) - 2*SUM(k, y(k,j) * SQRT(aii(i,j) * aii(k,j)))/a(j)) * log((v_mol(j)+b(j))/v_mol(j))
    );

dadT_ii_eq(i,j)$(ord(j) <= Ns)  ..
    dadT_ii(i,j) =E= -(ac(i) * kappa(i)) * sqrt(alpha_ij(i,j) * Tr(i,j)) / T(j);

dadT_eq(j)$(ord(j) <= Ns)  ..
    dadT(j) =E= - SUM(i, SUM(k, y(i,j) * y(k,j) * sqrt(dadT_ii(i,j) * dadT_ii(k,j))));

enthalpy_eq(j)$(ord(j) <= Ns)  ..
    HR(j) =E= R*T(j)*(Z(j)-1 + (T(j)*dadT(j) - a(j)) / (b(j)*R*T(j)) * log((v_mol(j) + b(j))/v_mol(j)));

* Ideal gas model
cubic_eos_ideal(j)$(ord(j) <= Ns)  ..
    Z(j) =E= 1;
    
enthalpy_eq_ideal(j)$(ord(j) <= Ns)  ..
    HR(j) =E= 0;
    
fugacity_eq_ideal(i,j)$(ord(j) <= Ns)  ..
    phi(i,j) =E=1;

* ============================================================================ *
* NRTL model
* ============================================================================ *

Parameters
    b_nrtl(i,i) /
        1.2  595.529982, 1.3 110.239086, 1.4 226.373398,
        2.1 164.57256,  2.3 141.963213, 2.4 187.104064,
        3.1 -95.4252161,3.2 623.581001, 3.4 219.73407,
        4.1 -177.88565, 4.2 344.481315, 4.3 -172.59152
    /
    c_nrtl(i,i) /
    1.1  0.0,  1.2  0.3,  1.3  0.3,  1.4  0.3,
    2.1  0.3,  2.2  0.0,  2.3  0.3,  2.4  0.3,
    3.1  0.3,  3.2  0.3,  3.3  0.0,  3.4  0.3,
    4.1  0.3,  4.2  0.3,  4.3  0.3,  4.4  0.0
    /;


*Parameter a_nrtl(i,i);
*a_nrtl(i,i) = 0

Alias(i,k);
Alias(i,h);
Alias(i,m);

Variables
    gamma_nrtl(i,j);

Equations
    Compute_gamma(i,j), Compute_gamma_ideal(i,j);

*calc_tao(i,k,j) ..
*    tao_nrtl(i,k,j) =E= b_nrtl(i,k) / T(j);
    
*calc_g(i,k,j)..
*    g_nrtl(i,k,j) =E= exp(-c_nrtl(i,k) * tao_nrtl(i,k,j));

$macro tao_nrtl(i,k,j)  ( b_nrtl(i,k) / T(j) )
$macro g_nrtl(i,k,j) (exp(-c_nrtl(i,k) * ( b_nrtl(i,k) / T(j) ) ) )

*Equation 29
Compute_gamma(i,j)$(ord(j) <= Ns) ..
    gamma_nrtl(i,j) =E= exp(
        SUM(h, x(h,j) * tao_nrtl(h,i,j) * g_nrtl(h,i,j)) /
        SUM(k, x(k,j) * g_nrtl(k,i,j))
    +
        SUM(h, (x(h,j) * g_nrtl(i,h,j)) /
        SUM(k, x(k,j) * g_nrtl(k,h,j)) *
        (tao_nrtl(i,h,j) - 
        SUM(m, x(m,j) * tao_nrtl(m,h,j) * g_nrtl(m,h,j)) /
        SUM(k, x(k,j) * g_nrtl(k,h,j))))
        );

Compute_gamma_ideal(i,j)$(ord(j) <= Ns) ..
    gamma_nrtl(i,j) =E= 1;
    
* ============================================================================ *
* Partiion coeffient
* ============================================================================ *

Positive Variable K_part(i,j);

Equation K_equation(i,j);

K_equation(i,j)$(ord(j) <= Ns)  ..
    K_part(i,j)*(phi(i,j) * P(j)) =E= gamma_nrtl(i,j) * exp(logPsat(i,j));
    

* ============================================================================ *
* Enthalpy vapor phase
* ============================================================================ *

SET cc /1*5/;

PARAMETERS
    C_ig(i,cc) Heat capacity coefficients for ideal gas;

TABLE C_ig(i,cc)
       1              2              3              4              5
1  111.022862   -0.750243683  0.00374664805  -6.16980174E-06   3.65183052E-09
2   78.1195715  -0.495919524  0.00257066991  -4.30332398E-06   2.57806954E-09
3  150.319619   -0.966221476  0.00404410804  -6.00376885E-06   3.23582716E-09
4   47.4579586  -0.349114826  0.00180207897  -3.04885369E-06   1.83440907E-09;
    

SCALAR T0 "Reference temperature [K]" /298.15/; 

VARIABLES DH_ig(i,j) "Ideal gas enthalpy change [J/mol]";
VARIABLES Hvi(i,j) "Ideal gas enthlapy of vapor phase for component i + Hform[J/mol]";
VARIABLES Hv(j) "Total vapor phase enthalpy";
*VARIABLES DH_ig_2(i,j);

EQUATIONS calc_DH_ig, calc_Hvi, calc_Hv;

calc_DH_ig(i,j) $(ord(j) <= Ns) ..
    DH_ig(i,j) =E= (C_ig(i,'1')*(T(j) - T0)
               + C_ig(i,'2')*(T(j)**2 - T0**2)/2
               + C_ig(i,'3')*(T(j)**3 - T0**3)/3
               + C_ig(i,'4')*(T(j)**4 - T0**4)/4
               + C_ig(i,'5')*(T(j)**5 - T0**5)/5);

calc_Hvi(i,j) $(ord(j) <= Ns) ..
    Hvi(i,j) =E= DH_ig(i,j) + Hform(i);

calc_Hv(j)$(ord(j) <= Ns) ..
    Hv(j)/H_factor =E= (sum(i, y(i,j) * Hvi(i,j)) + HR(j));


* ============================================================================ *
* Enthalpy liquid phase
* ============================================================================ *
PARAMETERS
    C_liq(i,cc) Heat capacity coefficients for liquid phase,
    C_vap(i,cc) DIPPR coefficients for vaporization enthalpy;
    
TABLE C_liq(i,cc)
       1              2              3                  4                   5
1  2.42239762e+05  -2.74746461e+03   1.16801221e+01  -2.20506997e-02   1.56048214e-05
2  1.21280234e+03  -1.30223888e+01   5.71358761e-02  -1.12370042e-04   8.73258381e-08
3  3.08689696e+05  -3.49815711e+03   1.48572135e+01  -2.80226414e-02   1.98131733e-05
4  1.20736602e+03  -1.29812644e+01   5.93572571e-02  -1.16271343e-04   8.66983296e-08;

TABLE C_vap(i,cc)
       1          2      3          4         5
1   33774     0.5107    -0.17304   0.05181    0
2   65831     1.1905    -1.7666    1.0012     0
3   43172     1.5334    -1.9       0.83816    0
4   45290     0.27343    0.21645  -0.11756    0;
    
Parameter Tbr(i)
          alpha_Tb(i)
          ai_Tb
          e
          sigma
          ;
Tbr(i) = Tb(i)/Tc(i);

Parameter DH_vap(i) "Vaporization enthalpy at Tb [J/mol]";
DH_vap(i) = C_vap(i,'1')*(1 - Tbr(i))**(C_vap(i,'2')
            + C_vap(i,'3')*Tbr(i)
            + C_vap(i,'4')*Tbr(i)**2
            + C_vap(i,'5')*Tbr(i)**3);

Parameter DH_ig_Tb(i) "DH_ig from T0 to TB [J/mol]";
DH_ig_Tb(i) = C_ig(i,'1')*(Tb(i) - T0)
               + C_ig(i,'2')*(sqr(Tb(i)) - sqr(T0))/2
               + C_ig(i,'3')*(power(Tb(i),3) - power(T0,3))/3
               + C_ig(i,'4')*(power(Tb(i),4) - power(T0,4))/4
               + C_ig(i,'5')*(power(Tb(i),5) - power(T0,5))/5;

VARIABLES DH_liq(i,j) "Liquid enthalpy change from Tb to T [J/mol]";
VARIABLES Hv_Tbi_P(i,j) "Vapor phase enthalpy at Tb [J/mol]";
VARIABLES Hli(i,j) "Liquid enthalpy for each component [J/mol]";
VARIABLES Hl(j) "Total liquid phase enthalpy [J/mol]";
* SRK for pure components at Tb
* Define the equation for Z_vi

alpha_Tb(i) = (1 + kappa(i) * (1 - Tbr(i) ** 0.5)) ** 2;
ai_Tb(i) = alpha_Tb(i)*ac(i);

*Cubic EOS - SRK

e = 0;
sigma = 1;

Positive Variable Zvi(i,j), v_mol_i(i,j);
Equation Zvi_eq(i,j), v_eq(i,j), Zvi_eq_ideal(i,j);
Variable  HR_pure(i,j) "Residual enthalpy for pure substance";

$macro Betav_p(i,j) ( bii(i)*P(j)/(R*Tb(i)) )
$macro qv_p(i) ( ai_Tb(i) / (bii(i)*R*Tb(i)) )

Zvi_eq(i,j) $(ord(j) <= Ns) ..
    -Zvi(i,j) + 1 + Betav_p(i,j) -
    qv_p(i) * Betav_p(i,j) * ( Zvi(i,j) - Betav_p(i,j) ) / ( Zvi(i,j)*( Zvi(i,j) + Betav_p(i,j) ) ) =E= 0;

v_eq(i,j) $(ord(j) <= Ns) ..
    v_mol_i(i,j) =E= Zvi(i,j)*R*Tb(i)/P(j);    

Variable  HR_pure(i,j);
Parameter dadT_Tb(i);
Equation calc_HR_pure(i,j), calc_HR_pure_ideal(i,j);

dadT_Tb(i) = - ac(i)*kappa(i)*sqrt(alpha_Tb(i)*Tbr(i))/Tb(i);
    
calc_HR_pure(i,j) $(ord(j) <= Ns) ..
    HR_pure(i,j) =E= R * Tb(i) * (Zvi(i,j) - 1 + (Tb(i) * dadT_Tb(i) - ai_Tb(i)) / (bii(i) * R * Tb(i)) * log((v_mol_i(i,j) + bii(i)) / v_mol_i(i,j)));

Zvi_eq_ideal(i,j) $(ord(j) <= Ns) .. 
    Zvi(i,j) =E= 1;

calc_HR_pure_ideal(i,j) $(ord(j) <= Ns) ..
    HR_pure(i,j) =E= 0;


EQUATIONS calc_DH_liq, calc_Hli, calc_Hl, calc_Hv_Tbi_P(i,j);

calc_DH_liq(i,j) $(ord(j) <= Ns) ..
        DH_liq(i,j) =e= C_liq(i,'1')*(T(j) - Tb(i)) 
        + C_liq(i,'2')*(T(j)**2 - Tb(i)**2)/2 
        + C_liq(i,'3')*(T(j)**3 - Tb(i)**3)/3 
        + C_liq(i,'4')*(T(j)**4 - Tb(i)**4)/4 
        + C_liq(i,'5')*(T(j)**5 - Tb(i)**5)/5;
        
calc_Hv_Tbi_P(i,j) $(ord(j) <= Ns) ..
    Hv_Tbi_P(i,j) =E= DH_ig_Tb(i) + Hform(i) + HR_pure(i,j);

calc_Hli(i,j) $(ord(j) <= Ns) ..
    Hli(i,j) =E= DH_liq(i,j) + Hv_Tbi_P(i,j) - DH_vap(i);

calc_Hl(j)$(ord(j) <= Ns) ..
    Hl(j)/H_factor =E= sum(i, x(i,j) * Hli(i,j));
        
* ============================================================================ *
* Chemical reaction rate
* ============================================================================ *
*SETS
*   reactive_trays(j) "Reactive trays"  /3,4,5 /

PARAMETER  reactive_trays(j);

*reactive_trays(j)$(ord(j) = NR) = 1;
*reactive_trays(j)$(ord(j) = NR+1) = 1;
*reactive_trays(j)$(ord(j) = NR+2) = 1;


VARIABLES
    K_eq(j) "Reaction equilibrium constant"
    k_rate(j) "Specific reaction rate [mol/(kg_cat-min)]"
    k_A(j) "Adsorption rate"
    Rxn_Rate(j) "Reaction rate [mol/(kg_cat-min)]";


$macro Xg(i,j) ( x(i,j)*gamma_nrtl(i,j) ) 
* Ignore inert component (i=1)
*calc_Xg(i,j)$(ord(i) > 1) ..  
*    Xg(i,j) =E= x(i,j)*gamma_nrtl(i,j);
    
* Reaction equilibrium constant
EQUATION calc_K_eq(j);
calc_K_eq(j)$((ord(j) >= NR) AND (ord(j) <= NR+2)) ..
    K_eq(j) =E= EXP(10.387 + 4060.59/T(j) - 2.89055*LOG(T(j))
                - 0.01915144*T(j) + 5.28586E-5*T(j)**2 - 5.32977E-8*T(j)**3);

* Specific reaction rate
EQUATION calc_k_rate(j);
calc_k_rate(j)$((ord(j) >= NR) AND (ord(j) <= NR+2)) ..
    k_rate(j) =E= 7.41816E15 * EXP(-60.4E3 / (R * T(j))) / 60;

* Adsorption rate
EQUATION calc_k_A(j);
calc_k_A(j)$((ord(j) >= NR) AND (ord(j) <= NR+2)) ..
    k_A(j) =E= EXP(-1.0707 + 1323.1 / T(j));

* Reaction rate equation
EQUATION calc_Rx_Rate(j);
calc_Rx_Rate(j)$((ord(j) >= NR) AND (ord(j) <= NR+2)) ..
*    Rxn_Rate(j)*K_eq(j)*Xg('2',j) =E= k_rate(j)*(Xg('2',j))**2*(Xg('3',j)*K_eq(j)*Xg('2',j) - Xg('4',j))/(1+k_A(j)*Xg('2',j))**3;
    Rxn_Rate(j) =E= k_rate(j)*Xg('2',j) * ( Xg('2',j) * Xg('3',j) - Xg('4',j)/K_eq(j) ) / (1 + k_A(j)*Xg('2',j) )**3;
        
* Set reaction rate to zero for non-reactive trays
EQUATION assign_zero(j);
assign_zero(j)$( ((ord(j) < NR) OR (ord(j) > NR+2)) AND (ord(j) <= Ns) )  ..
   Rxn_Rate(j) =E= 0;



* ============================================================================ *
* MESHR model
* ============================================================================ *

EQUATIONS mol_balance_eq(j), mol_balance_condenser_eq, mol_balance_reboiler_eq,
          comp_balance_eq(i,j), comp_balance_condenser_eq(i), comp_balance_reboiler_eq(i,j),
          eq_rel_eq(i,j), summation_eq(j), energy_balance_eq(j), energy_balance_condenser_eq,
          energy_balance_reboiler_eq(j), spec_1_eq(j), spec_2_eq(j);
          
Parameters
*Spec_1 "Reflux ratio" /3/
*Spec_2 "Bottom to feed ratio" /0.233/;
 
Spec_1 "Ethanol conversion" /0.85/
Spec_2 "ETBE molar fraction on the bottom produtct" /0.83/;

* Total Molar Balance
mol_balance_eq(j)$((ord(j) GT 1) AND (ord(j) LT Ns)) ..
    L(j) + V(j) - V(j+1) - L(j-1)
    - FE$(ord(j) EQ NFE)
    - FB$(ord(j) EQ NFB)
    - SUM(i, v_i(i) * m_cat * Rxn_Rate(j))/F_factor =E= 0;

* Special case for the condenser (first stage)
mol_balance_condenser_eq ..
    V('1') + (L('1') + D) - V('2') =E= 0;

* Special case for the reboiler (last stage)
*mol_balance_reboiler_eq ..
*    V('10') + L('10') - L('9') =E= 0;
    
* Special case for the reboiler (last stage)
mol_balance_reboiler_eq(j)$(ord(j) EQ Ns) ..
    V(j) + L(j) - L(j-1) =E= 0;

* Component Molar Balance
comp_balance_eq(i,j)$((ord(j) GT 1) AND (ord(j) LT Ns)) ..
        L(j) * x(i,j) + V(j) * y(i,j) - V(j+1) * y(i,j+1) - L(j-1) * x(i,j-1) 
       - FE$(ord(j) EQ NFE) * zE(i) - FB$(ord(j) EQ NFB) * zB(i) - m_cat * v_i(i) * Rxn_Rate(j)/F_factor =E= 0;
    
* Special case for component balance in the condenser
comp_balance_condenser_eq(i) ..
    V('1') * y(i,'1') + (L('1') + D) * x(i,'1') - V('2') * y(i,'2') =E= 0;

* Special case for component balance in the reboiler
*comp_balance_reboiler_eq(i) ..
*    V('10') * y(i,'10') + L('10') * x(i,'10') - L('9') * x(i,'9') =E= 0;
    
* Special case for component balance in the reboiler
comp_balance_reboiler_eq(i,j)$(ord(j) EQ Ns) ..
    V(j) * y(i,j) + L(j) * x(i,j) - L(j-1) * x(i,j-1) =E= 0;

* Equilibrium Relations
eq_rel_eq(i,j) $(ord(j) <= Ns) ..
    K_part(i,j) * x(i,j) - y(i,j) =E= 0;

* Summation Constraint
summation_eq(j) $(ord(j) <= Ns) ..
    SUM(i, y(i,j) - x(i,j)) =E= 0;

* Energy Balance
energy_balance_eq(j)$((ord(j) GT 1) AND (ord(j) LT Ns)) ..
    V(j) * Hv(j) + L(j) * Hl(j) - V(j+1) * Hv(j+1) - L(j-1) * Hl(j-1) - FE$(ord(j) EQ NFE)*HFE - FB$(ord(j) EQ NFB)*HFB
    =E= 0;

* Special case for the condenser energy balance
energy_balance_condenser_eq ..
    V('1') * Hv('1') + (L('1') + D) * Hl('1') - V('2') * Hv('2') + Qc =E= 0;

* Special case for the reboiler energy balance
*energy_balance_reboiler_eq ..
*   V('10') * Hv('10') + L('10') * Hl('10') - L('9') * Hl('9') - Qr =E= 0;
    
* Special case for the reboiler energy balance
energy_balance_reboiler_eq(j)$(ord(j) EQ Ns) ..
    V(j) * Hv(j) + L(j) * Hl(j) - L(j-1) * Hl(j-1) - Qr =E= 0;

* Specifications
spec_1_eq(j)$(ord(j) EQ Ns) .. (FE - D*x['2','1'] - L[j]*x['2',j])/FE =e= 0.85;
spec_2_eq(j)$(ord(j) EQ Ns) .. 0.83 - x['4', j] =E= 0;

Equation spec_3_eq;
spec_3_eq ..  L['1']/D =g= 9;

Equation spec_4_eq;
spec_4_eq ..  L['1']/D =l= 15;
*spec_1_eq(j)$(ord(j) EQ Ns)  .. (FE - D*x['2','1'] - L[j]*x['2',j])/FE  =E= Spec_1;
*  molar flowarate spec 
*spec_2_eq(j)$(ord(j) EQ Ns)  .. x['4', j] =G= Spec_2;

* Specifications
*spec_1_eq .. (FE - D*x['2','1'] - L['10']*x['2','10'])/FE  =g= Spec_1;
*spec_2_eq .. x['4', '10'] =g= Spec_2;

*spec_1_eq .. D - L['1']/Spec_1 =E= 0;
*spec_2_eq .. (FE+FB)*Spec_2 - L['10'] =E= 0;

equation V1_eq;
V1_eq .. V('1') =E= 0;



*Psat


* Objective function cost
parameters
           CostETBE    "Price of ETBE bottom product [$/mol]"
           CostEth      "Price of Ethanol [$/mol]"
           CostBut      "Price of Butanes [$/mol]";

CostETBE=25.3e-3; 
CostEth=15e-3;
CostBut=8.25e-3;

* Calculate steam cost: $/J * J/s * 8150 h/year * 3600 s/h
*Steam_Cost = c_steam * Qr * 8150 * 3600;

* Cooling cost
*Cooling_Cost = 0.378e-9*Qc*8150*3600;

* Diameter calculation
VARIABLE
D_col(j)                    'Column diameter (ft)'
Dcol_max                "Maximum column diameter (ft)"
Mw_mix(j)  "Molecular weights (kg/kmol)"
;
EQUATION
def_D_col(j)            'Equation for column diameter calculation'
def_Mw_mix(j)      "Define mixture molar mass at stage j"
def_Dcol_max           'Equation for maximumcolumn diameter calculation';
SCALAR FF 'Parameter for flooding velocity estimation - Douglas book #kg^1/2 m^-1/2 s^-1 ';

def_Mw_mix(j)$(ord(j) <= Ns) ..
Mw_mix(j) =E= sum( i, y(i,j)*Mw(i) );
FF = 1.51*(0.45359237/0.3048)**0.5;

$macro u FF/sqrt (Mw_mix(j)*1e-3/v_mol(j) )
$macro  Q V(j)*F_factor*v_mol(j)/60 
*$macro Area Q/(u*0.88)
*$macro Vmass V(j)*F_factor/60*Mw_mix(j)*1e-3
*$macro rho_v Mw_mix(j)/ v_mol(j)

def_D_col(j)$( (ord(j) > 1) AND (ord(j) < Ns) ) ..
    D_col(j) =e=sqrt(4*(Q/(u*0.88))/Pi)*3.28084;
*def_D_col(j) ..
*    D_col(j) =e= sqrt(4 * (V(j)*F_factor*v_mol(j)/60) / ( ( FF/sqrt(Mw_mix(j)*1e-3/v_mol(j)) ) * 0.88 * Pi) );
*def_Dcol_max.l ..
*    Dcol_max.l =e= D_col('1');
def_Dcol_max(j)$( (ord(j) > 1) AND (ord(j) < Ns) ) ..
    Dcol_max =g= D_col(j);

    
* Column and tray cost
SCALAR
MS           "Marshall &  Swift Parameter"  /1773.4/;


*EQUATION def_Lc(j);
*VARIABLE Lc;

*def_Lc(j)$(ord(j) = Ns) ..
*   Lc =E=  0.7315 * Nt(j)$(ord(j) = Ns) * 3.28084
$macro Lc (0.7315*(Ns-2)*3.28084)
$macro ColCost  (MS/280*(101.9*Dcol_max**1.066*(0.7315*(Ns-2)*3.28084)**0.802*7.05))
$macro TrayCost (MS/280*(4.7*Dcol_max**1.55*(0.7315*(Ns-2)*3.28084)*2.7))



*Reboiler cost
EQUATION def_Treb(j);
EQUATION def_Tcond;
VARIABLE
 Tcond "Condenser temperature (K)"
 Treb "Reboiler Temperature (K)"; 

def_Treb(j)$(ord(j) = Ns) ..
    Treb =E= T(j)$(ord(j) = Ns);


*$macro Areb  Qr*FH_factor/((250/0.17611)*(433.15-Treb))*10.7639
*$macro RebCost  (MS/280*101.3*(Areb**.65*7.3525))
$macro RebCost  (MS/280*101.3*(( Qr*FH_factor/60/((250/0.17611)*(433.15-Treb))*10.7639)**.65*7.3525))

* Condenser cost
def_Tcond ..
    Tcond =E= T('1');
*$macro DTlm  (10)/log((Tcond-303.15)/(Tcond-313.15))
*$macro Acond  Qc*FH_factor/60/(150/0.17611 * DTlm)*10.7639
*$macro CondCost  ((MS*101.3*2.665*(Acond**0.65))/280)
$macro CondCost  (Qc*FH_factor/60/(150/0.17611 * ((10)/log((Tcond-303.15)/(Tcond-313.15))))*10.7639)**.65*1709.8394439285714


* Steam costs [$/J]
Parameter c_steam, T_steam, steam_price;
* Determine steam temperature and cost
c_steam = 4.54e-9;
T_steam = 433.15;

$macro CatCost  7.7*3*0.4
*VARIABLE
*CAP_costA  "Capital cost 10^6$/yr"
*OP_cost  "Operational cost 10^6$/yr"
*;
*
*Equations
*CAP_eq "Capital cost";
**OP_eq "Operational Cost";
*CAP_eq ..CAP_costA =E= ColCost;
*
$macro CAP_cost 1/3*(TrayCost + ColCost + CondCost + RebCost + CatCost)
$macro OP_cost (0.378e-9*Qc+  c_steam * Qr)* 8150 * 3600 * FH_factor + (FE*CostEth + FB*CostBut- CostETBE*Breb)*F_factor

Equation def_Breb(j), def_Profit;
Variable Breb, Profit;
def_Breb(j)$(ord(j) = Ns) .. Breb =E= L(j)$(ord(j) = Ns);
*def_Profit .. Profit*1e6=E= CAP_cost+OP_cost; 
def_Profit .. Profit*1e6=E=  1/3*(TrayCost + ColCost + CondCost + RebCost + CatCost)+(0.378e-9*Qc+  c_steam * Qr) * 8150 * 60 * FH_factor + (FE*CostEth + FB*CostBut- CostETBE*Breb)*F_factor* 8150 * 60; 



EQUATION obj_def;
VARIABLE obj;

obj_def ..
*        obj =E= Qr + Dcol_max;
    obj =E= Profit;



MODEL MESHR_simple /
psat_def,
cubic_eos_ideal,
v_definition,
fugacity_eq_ideal,
enthalpy_eq_ideal,
*alpha_definition,
*aii_definition,
*mixing_rule_b,
*mixing_rule_a,
*dadT_ii_eq,
*dadT_eq,
Compute_gamma_ideal,
K_equation,
calc_DH_ig,
calc_Hvi,
calc_Hv,
Zvi_eq_ideal,
v_eq,
calc_HR_pure_ideal,
calc_DH_liq,
calc_Hli,
calc_Hl,
calc_Hv_Tbi_P,
calc_K_eq,
calc_k_rate,
calc_k_A,
calc_Rx_Rate,
assign_zero,
mol_balance_eq,
mol_balance_condenser_eq,
mol_balance_reboiler_eq,
comp_balance_eq,
comp_balance_condenser_eq,
comp_balance_reboiler_eq,
eq_rel_eq,
summation_eq,
energy_balance_eq,
energy_balance_condenser_eq,
energy_balance_reboiler_eq,
spec_1_eq,
*spec_2_eq,
V1_eq,
def_D_col,
def_Mw_mix,
def_Dcol_max,
def_Profit,
def_Treb,
def_Tcond,
def_Breb,
*CAP_eq,
*OP_eq,
obj_def
/;

MODEL MESHR_Rigorous /
psat_def,
cubic_eos,
v_definition,
fugacity_eq,
enthalpy_eq,
alpha_definition,
aii_definition,
mixing_rule_b,
mixing_rule_a,
dadT_ii_eq,
dadT_eq,
Compute_gamma, 
K_equation,
calc_DH_ig,
calc_Hvi,
calc_Hv,
Zvi_eq,
v_eq,
*dadT_Tb_eq,
calc_HR_pure,
calc_DH_liq,
calc_Hli,
calc_Hl,
calc_Hv_Tbi_P,
calc_K_eq,
calc_k_rate,
calc_k_A,
calc_Rx_Rate,
assign_zero,
mol_balance_eq,
mol_balance_condenser_eq,
mol_balance_reboiler_eq,
comp_balance_eq,
comp_balance_condenser_eq,
comp_balance_reboiler_eq,
eq_rel_eq,
summation_eq,
energy_balance_eq,
energy_balance_condenser_eq,
energy_balance_reboiler_eq,
spec_1_eq,
spec_2_eq,
V1_eq,
def_D_col,
def_Mw_mix,
def_Dcol_max,
def_Profit,
def_Treb,
def_Tcond,
def_Breb,
*CAP_eq,
*OP_eq,
obj_def,
spec_3_eq,
*spec_4_eq,
/;



SCALAR nsmin, nsmax;
nsmin = 6;
nsmax = 30;

* Positive Variables
L.lo(j) = 1e-5 ;  L.up(j) = 50; L.l(j) = .5;
V.lo(j) = 1e-5 ;  V.up(j) = 50; V.l(j) = .5;
V.lo('1') = 0; V.up('1') = 0.001; V.l('1') = 0;
x.lo(i,j) = 0;  x.up(i,j) = 1; x.l(i,j) = 0.25;
y.lo(i,j) = 0;  y.up(i,j) = 1; y.l(i,j) = 0.25;
D.lo = 0.01; D.up = 1;
Qc.lo = 0.01;  Qc.up = 10; 
Qr.lo = 0.01;  Qr.up = 10;
Breb.lo = 0.001; Breb.up = 1;
T.lo(j) = Tmin; T.up(j) = Tmax; T.scale(j) = TF_factor;

logPsat.lo(i,j) = (
                kk(i,'j1') + kk(i,'j2')/(kk(i,'j3') + Tmin) 
               + kk(i,'j4')*Tmin + kk(i,"j5")*LOG(Tmin) 
               + kk(i,'j6')*(Tmin**kk(i,'j7')));
               

logPsat.up(i,j) = (
                kk(i,'j1') + kk(i,'j2')/(kk(i,'j3') + Tmax) 
               + kk(i,'j4')*Tmax + kk(i,"j5")*LOG(Tmax) 
               + kk(i,'j6')*(Tmax**kk(i,'j7')));
logPsat.scale(i,j) = 10; 

* Apply linear profile
*    jmax = Ns;
*    T.L(j) = Tmin + (ord(j) - jmin) / (jmax - jmin) * (Tmax - Tmin);
*                    T.L(j) = Tmin + (ord(j) - jmin) / (jmax - jmin) * (Tmax - Tmin);
    
* SRK parameters and properties
*dadT_Tb.lo(i) = -0.01; dadT_Tb.up(i) = -0.00001; 
    
alpha_ij.lo(i,j) = 0.5;  alpha_ij.up(i,j) = 2;
aii.lo(i,j) = 1;  aii.up(i,j) = 4;
a.lo(j) = 1;  a.up(j) = 3;
b.lo(j) = 1e-5;  b.up(j) = 1e-3; b.scale(j) = 1e-5;
dadT_ii.lo(i,j) = -0.009;  dadT_ii.up(i,j) = -0.001; dadT_ii.scale(i,j) = 1e-3;
dadT.lo(j) = -0.009;  dadT.up(j) = -0.001; dadT.scale(j) = 1e-3;
K_part.lo(i,j) = 0.001 ; K_part.up(i,j) = 10;
phi.lo(i,j) = 0.1; phi.up(i,j) = 1.2;
            
*    *Enthalpy Variables
DH_ig.lo(i,j) = (C_ig(i,'1')*(Tmin - T0)
               + C_ig(i,'2')*(Tmin**2 - T0**2)/2
               + C_ig(i,'3')*(Tmin**3 - T0**3)/3
               + C_ig(i,'4')*(Tmin**4 - T0**4)/4
               + C_ig(i,'5')*(Tmin**5 - T0**5)/5);
DH_ig.up(i,j) = (C_ig(i,'1')*(Tmax - T0)
               + C_ig(i,'2')*(Tmax**2 - T0**2)/2
               + C_ig(i,'3')*(Tmax**3 - T0**3)/3
               + C_ig(i,'4')*(Tmax**4 - T0**4)/4
               + C_ig(i,'5')*(Tmax**5 - T0**5)/5);
DH_ig.scale(i,j) = 1e3;
                
Hvi.lo(i,j) = DH_ig.lo(i,j) + Hform(i);
Hvi.up(i,j) = DH_ig.up(i,j) + Hform(i);
Hvi.scale(i,j) = 1e5;

HR.lo(j) = -2500;  HR.up(j) = -1100; HR.scale(j) = 1e3;

Hv.lo(j) = ( smin(i, Hvi.lo(i,j)) + HR.lo(j) ) * H_factor;
Hv.up(j) = ( smax(i, Hvi.up(i,j)) + HR.up(j) ) * H_factor;


DH_liq.lo(i,j) = C_liq(i,'1')*(Tmin - Tb(i)) 
                + C_liq(i,'2')*(Tmin**2 - Tb(i)**2)/2 
                + C_liq(i,'3')*(Tmin**3 - Tb(i)**3)/3 
                + C_liq(i,'4')*(Tmin**4 - Tb(i)**4)/4 
                + C_liq(i,'5')*(Tmin**5 - Tb(i)**5)/5;
DH_liq.up(i,j) = C_liq(i,'1')*(Tmax - Tb(i)) 
                + C_liq(i,'2')*(Tmax**2 - Tb(i)**2)/2 
                + C_liq(i,'3')*(Tmax**3 - Tb(i)**3)/3 
                + C_liq(i,'4')*(Tmax**4 - Tb(i)**4)/4 
                + C_liq(i,'5')*(Tmax**5 - Tb(i)**5)/5;
DH_liq.scale(i,j) = 1e4;
                
HR_pure.lo(i,j) = -3e3;  HR_pure.up(i,j) = -1e3; HR_pure.scale(i,j) = 1e3;
Hv_Tbi_P.lo(i,j) = DH_ig_Tb(i)*0.99 + Hform(i) + HR_pure.lo(i,j);
Hv_Tbi_P.up(i,j) = DH_ig_Tb(i)*1.01 + Hform(i) + HR_pure.up(i,j);
Hv_Tbi_P.scale(i,j) = 1e4;

Hli.lo(i,j) = DH_liq.lo(i,j) + Hv_Tbi_P.lo(i,j) - DH_vap(i);
Hli.up(i,j) = DH_liq.up(i,j) + Hv_Tbi_P.up(i,j) - DH_vap(i);
Hli.scale(i,j) = 1e5;

Hl.lo(j) = smin(i, Hvi.lo(i,j)) * H_factor;
Hl.up(j) = smax(i, Hvi.up(i,j)) * H_factor;
    
*    * molar volume
Z.lo(j) = 0.75; Z.up(j) = 0.9; Z.l(j) = 0.822;
v_mol.lo(j) = Z.lo(j)*R*Tmin/(9.5e5);
v_mol.up(j) = Z.up(j)*R*Tmax/(9.5e5);
v_mol.scale(j) = 1e-3;

Zvi.lo(i,j) = 0.75; Zvi.up(i,j) = 1;  Zvi.l(i,j) = 0.9;
v_mol_i.lo(i,j) = Zvi.lo(i,j)*R*Tmin/(9.5e5);
v_mol_i.up(i,j) = Zvi.up(i,j)*R*Tmax/(9.5e5);
v_mol_i.scale(i,j) = 1e-3;
    
*    * NRTL Parameters
gamma_nrtl.lo(i,j) =  0.01;  gamma_nrtl.up(i,j) = 20;
*    *tao_nrtl.lo(i,k,j) =-5; tao_nrtl.up(i,k,j) = 5; tao_nrtl.l(i,k,j) = 1;
*    *g_nrtl.lo(i,k,j) = 1e-10; g_nrtl.up(i,k,j) = 2; g_nrtl.l(i,k,j) = 1;
    
*    * ReactionVariables
*    *Xg.lo(i,j) = 1e-6;  Xg.up(i,j) = 1;
*    *K_eq.lo(j) = 11;  K_eq.up(j) = 26;
K_eq.lo(j) = EXP(10.387 + 4060.59/Tmax - 2.89055*LOG(Tmax) - 0.01915144*Tmax + 5.28586E-5*Tmax**2 - 5.32977E-8*Tmax**3);
K_eq.up(j) = EXP(10.387 + 4060.59/Tmin - 2.89055*LOG(Tmin) - 0.01915144*Tmin + 5.28586E-5*Tmin**2 - 5.32977E-8*Tmin**3);
K_eq.scale(j) = 10;
*    *k_rate.lo(j) = 4e4;  k_rate.up(j) = 4e5; k_rate.scale(j) = 1e4;
k_rate.lo(j) = 7.41816E15 * EXP(-60.4E3 / (R * Tmin)) / 60;
k_rate.up(j) = 7.41816E15 * EXP(-60.4E3 / (R * Tmax)) / 60;
k_rate.scale(j) = 1e4;
*    *k_A.lo(j) = 14;  k_A.up(j) = 20;
k_A.lo(j) = EXP( -1.0707 + 1323.1 / Tmax );
k_A.up(j) = EXP( -1.0707 + 1323.1 / Tmin );
k_A.scale(j) = 10;
Rxn_Rate.lo(j) = -10;  Rxn_Rate.up(j) = 10;

D.l  = 0.5748006462583984;
Qc.l = 0.5451343354854316;
Qr.l = 0.4592022879232921;
Breb.l = 0.24;
*               cost estimations
Mw_mix.lo(j) = 45; Mw_mix.up(j) = 103; Mw_mix.scale(j) = 10;
D_col.lo(j) = 0.001; D_col.up(j) = 2; D_col.l(j)=0.03; D_col.scale(j)=0.1;
Dcol_max.lo = 0.001; Dcol_max.up = 2; Dcol_max.l=0.03; Dcol_max.scale=0.1;
Tcond.lo = Tmin; Tcond.up = Tmin+50; Tcond.scale = TF_factor;
Treb.lo = Tmax-50; Treb.up = Tmax; Treb.scale = TF_factor;
*FE is a parameter for this case
*FE.lo = 0.00075; FE.up = 37.43; FE.l = 5.774/F_factor; # FE is a parameter for this case
    
    

*   
*        for(NFB = NFE+2 to (Ns - 1) by 1,
*        
*            for(NR = NFE to (NFB-2) by 1,
*            

*                option NLP = CONOPT4;
*                option threads = 8;
*                option optca = 1e-8;
*                option optcr = 1e-8;
*                SOLVE MESHR_Simple USING NLP MINIMIZING obj;
*                X_ETOH(j)$(ord(j) = Ns) = (FE - D.l*x.l['2','1'] - L.l[j]*x.l['2',j])/FE;
*                X_D_ETBE(j)$(ord(j) = Ns)  = x.l['4',j];
*                RR = L.l('1')/D.l;
*                DISPLAY X_ETOH, x_D_ETBE, RR, Qr.l;




$macro CONDITION ( (ord(j) <= Ns) AND (ord(jj) EQ Ns)  AND (ord(zz) EQ NFE) AND (ord(xx) EQ NFB) AND (ord(hh) EQ NR) )
$macro CONDITION_2 ( (ord(j) = Ns) AND (ord(jj) EQ Ns) AND (ord(zz) EQ NFE) AND (ord(xx) EQ NFB) AND (ord(hh) EQ NR) )
$macro CONDITION_3 ( (ord(j) = 1) AND (ord(jj) EQ Ns) AND (ord(zz) EQ NFE) AND (ord(xx) EQ NFB) AND (ord(hh) EQ NR) )
                
Ns = 9;
NFE = 3;
NFB = 7;
NR = 3;

option NLP = CONOPT;
option reslim = 1000;
option optcr = 1e-4;
*option optca = 1e-4;
option threads = 12;
MESHR_Rigorous.scaleopt = 1;
SOLVE MESHR_Rigorous USING NLP MINIMIZING obj;
display MESHR_Rigorous.resusd;
display MESHR_Rigorous.resusd, MESHR_Rigorous.etSolve;
X_ETOH(j)$(ord(j) = Ns) = (FE - D.l*x.l['2','1'] - L.l[j]*x.l['2',j])/FE;
X_D_ETBE(j)$(ord(j) = Ns)  = x.l['4',j];
RR = L.l('1')/D.l;
test_par(j)$(ord(j) EQ Ns) = (FE - D.l*x.l['2','1'] - L.l[j]*x.l['2',j])/FE - Spec_1;
DISPLAY Ns, NFE, NFB, NR, Dcol_max.l, X_ETOH, x_D_ETBE, RR, Qr.l, MESHR_Rigorous.modelStat,test_par;

option NLP = BARON;
option reslim = 1000;
option optcr = 1e-4;
*option optca = 1e-4;
option threads = 12;
MESHR_Rigorous.scaleopt = 1;
SOLVE MESHR_Rigorous USING NLP MINIMIZING obj;
display MESHR_Rigorous.resusd;
display MESHR_Rigorous.resusd, MESHR_Rigorous.etSolve;
X_ETOH(j)$(ord(j) = Ns) = (FE - D.l*x.l['2','1'] - L.l[j]*x.l['2',j])/FE;
X_D_ETBE(j)$(ord(j) = Ns)  = x.l['4',j];
RR = L.l('1')/D.l;
DISPLAY Ns, NFE, NFB, NR, Dcol_max.l, X_ETOH, x_D_ETBE, RR, Qr.l, MESHR_Rigorous.modelStat;


*option NLP = CONOPT;
*option reslim = 1000;
*option optcr = 1e-4;
**option optca = 1e-4;
*option threads = 12;
**MESHR_Rigorous.scaleopt = 1;
*SOLVE MESHR_Rigorous USING NLP MINIMIZING obj;

*option NLP = ANTIGONE;
*option reslim = 1000;
*option optcr = 1e-4;
*option threads = 12;
*SOLVE MESHR_Rigorous USING NLP MINIMIZING obj;
*display MESHR_Rigorous.resusd;
*display MESHR_Rigorous.resusd, MESHR_Rigorous.etSolve;
*X_ETOH(j)$(ord(j) = Ns) = (FE - D.l*x.l['2','1'] - L.l[j]*x.l['2',j])/FE;
*X_D_ETBE(j)$(ord(j) = Ns)  = x.l['4',j];
*RR = L.l('1')/D.l;
*test_par(j)$(ord(j) EQ Ns) = (FE - D.l*x.l['2','1'] - L.l[j]*x.l['2',j])/FE - Spec_1;
*DISPLAY Ns, NFE, NFB, NR, Dcol_max.l, X_ETOH, x_D_ETBE, RR, Qr.l, MESHR_Rigorous.modelStat,test_par;

*                * Mw_cost
Mw_mix_temp(j) = sum( i, y.l(i,j)*Mw(i) );
*                *     If the candidate is viable, save the results

Tcond_t = T.l('1');
model_stat(j,jj,zz,xx,hh)$(CONDITION_2) = MESHR_Rigorous.modelStat;
if ((MESHR_Rigorous.modelStat = 1) or (MESHR_Rigorous.modelStat = 2),

    Z_cost(j,jj,zz,xx,hh)$(CONDITION)          = Z.l(j);
    v_mol_cost(j,jj,zz,xx,hh)$(CONDITION)   = v_mol.l(j);
    Mw_cost(j,jj,zz,xx,hh)$(CONDITION)       = Mw_mix_temp(j);
    Qr_cost(j,jj,zz,xx,hh)$(CONDITION_2)     = Qr.l*FH_Factor;
    Qc_cost(j,jj,zz,xx,hh)$(CONDITION_2)     = Qc.l*FH_factor;
    T_cost(j,jj,zz,xx,hh)$(CONDITION)          = T.l(j);
    V_cost(j,jj,zz,xx,hh)$(CONDITION)          = V.l(j)*F_factor;
    L_cost(j,jj,zz,xx,hh)$(CONDITION)          = L.l(j)*F_factor;
    x_column(i,j,jj,zz,xx,hh)$(CONDITION)    = x.l(i,j);
    y_column(i,j,jj,zz,xx,hh)$(CONDITION)    = y.l(i,j);
    Treb_cost(j,jj,zz,xx,hh)$(CONDITION_2)  = T.l(j);
    Tcond_cost(j,jj,zz,xx,hh)$(CONDITION_2)  = Tcond_t;
    RR_array(j,jj,zz,xx,hh)$(CONDITION_2)         =  RR;
    x_ETBE_array(j,jj,zz,xx,hh) $(CONDITION_2) =  X_D_ETBE(j);
    x_ETBE_D_array(j,jj,zz,xx,hh)$(CONDITION_3) = x.l['4','1'];
    X_EtOH_array(j,jj,zz,xx,hh) $(CONDITION_2) =  X_ETOH(j);
    profit_obj(j,jj,zz,xx,hh)$(CONDITION_2)  = obj.l;
    Breb_cost(j,jj,zz,xx,hh)$(CONDITION_2) = Breb.l*F_Factor;
    D_colMax_arr(j,jj,zz,xx,hh)$(CONDITION_2) = Dcol_max.l;
    D_col_arr(j,jj,zz,xx,hh)$(CONDITION)          = D_col.l(j);
    count_s = count_s+1;
    
    else
    
        Z_cost(j,jj,zz,xx,hh)$(CONDITION)        = -1;
        v_mol_cost(j,jj,zz,xx,hh)$(CONDITION) = -1;
        Mw_cost(j,jj,zz,xx,hh)$(CONDITION)     = -1;
        Qr_cost(j,jj,zz,xx,hh)$(CONDITION_2)   = -1;
        Qc_cost(j,jj,zz,xx,hh)$(CONDITION_2)   = -1;
        V_cost(j,jj,zz,xx,hh)$(CONDITION)          = -1;
        L_cost(j,jj,zz,xx,hh)$(CONDITION)          = -1;
        T_cost(j,jj,zz,xx,hh)$(CONDITION)          = -1;
        x_column(i,j,jj,zz,xx,hh)$(CONDITION)    = -1;
        y_column(i,j,jj,zz,xx,hh)$(CONDITION)    = -1;
        Treb_cost(j,jj,zz,xx,hh)$(CONDITION_2)  = -1;
        Tcond_cost(j,jj,zz,xx,hh)$(CONDITION_2)  = -1;
        RR_array(j,jj,zz,xx,hh)$(CONDITION_2)    =  -1;
        x_ETBE_array(j,jj,zz,xx,hh) $(CONDITION_2) =  -1;
        x_ETBE_D_array(j,jj,zz,xx,hh)$(CONDITION_3) = -1;
        X_EtOH_array(j,jj,zz,xx,hh) $(CONDITION_2) =  -1;
        profit_obj(j,jj,zz,xx,hh)$(CONDITION_2)  = -1;
        Breb_cost(j,jj,zz,xx,hh)$(CONDITION_2) = -1;
        D_colMax_arr(j,jj,zz,xx,hh)$(CONDITION_2) = -1;
        D_col_arr(j,jj,zz,xx,hh)$(CONDITION)          = -1;
        count_f = count_f+1;
);
*                        it = it+1;
**                        Save results to a GDX file named differently in each iteration
*                        execute_unload 'result_iter' + it:0 + '.gdx', result;
*
**                        Store result for later
*                        results(ord(iter)=it) = result;
*            );
*        );

Elapsed_time = timeElapsed;
Display   Elapsed_time
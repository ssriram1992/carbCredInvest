$TITLE "Renewable investment encouragement with carbon credits"

$SETGLOBAL Limit_Listing "*"

%Limit_Listing%$ONTEXT
$offlisting
option dispwidth=60;
$offsymlist
$Offsymxref
$oninline
Option Solprint = off;
Option limrow = 0;
Option limcol = 0;
$ONTEXT
$OFFTEXT


Sets 
countries /c1, c2, c3/
prodTypes /coal, gas, renew/
producers /p1, p2/
;

Alias(countries, cc);
Alias(prodTypes, pt);
Alias(producers, pp);

Parameters
LinCost(countries, prodTypes) 										"Linear cost of production"
QuadCost(countries, prodTypes)										"Quadratic cost of production"
LinInvCost(countries)															"Cost of investment in renewables"
QuadInvCost(countries)														"Cost of investment in renewables"
InitCapacity(countries, producers, prodTypes)			"Initial capacity"
InitCarb(countries, producers)										"Initial carbon credit allocation"
CarbonCost(prodTypes)															"Number of carbon credits consumed per unit of production"
DemInt(countries)                                 "Intercept of the demand curve"
DemSlope(countries)                               "Slope of the demand curve"
;


$ontext
$offtext

Loop(cc,
DemInt(cc) = ORD(cc)*1000+2000;
DemSlope(cc) = ORD(cc)*0.5+1.5;

LinInvCost(cc) = 10+20*ord(cc);
QuadInvCost(cc) = 0.5+0.1*(2-ord(cc));

Loop(pt,
	LinCost(cc, pt) = 30 - 10*ord(pt);
	QuadCost(cc, pt) = 1.5 - 0.5*ord(pt);
	Loop(pp,
		InitCapacity(cc, pp, pt) = 30 + ord(cc)*10 - ord(pp) * 5 + 20*(3-ORD(pt));
		InitCapacity(cc, pp, "renew") = 0;
		InitCarb(cc, pp) = 500;
	);
);
);
InitCapacity("c2", "p1", "renew") = 500;
CarbonCost("coal") = 7;
CarbonCost("gas") = 5;
CarbonCost("renew") = 0.1;

$ontext

#######
#        #    #   #####  ######  #####
#        ##   #     #    #       #    #
#####    # #  #     #    #####   #    #
#        #  # #     #    #       #####
#        #   ##     #    #       #   #
#######  #    #     #    ######  #    #

######
#     #    ##     #####    ##
#     #   #  #      #     #  #
#     #  #    #     #    #    #
#     #  ######     #    ######
#     #  #    #     #    #    #
######   #    #     #    #    #

#     #
#     #  ######  #####   ######
#     #  #       #    #  #
#######  #####   #    #  #####
#     #  #       #####   #
#     #  #       #   #   #
#     #  ######  #    #  ######

$offtext

* Primal variables
Positive Variables
PRODUCTION(countries, producers, prodTypes)
INVESTMENT(countries, producers)
CARBON_BUY(countries, producers)
CARBON_SELL(countries, producers)
CARBON_BUY_C(countries)
CARBON_SELL_C(countries)
;
Variable
CARBPRICE(countries)
INTERNATIONAL_CARBPRICE
OBJ_CP(countries, producers)
OBJ
;

Equations
ep_CarbNeutr(countries, producers)							"Production is limited by carbon capacity"
ep_capacity(countries, producers, prodTypes) 		"Production capacity limits"
eo_producer(countries, producers)								"Producer's objective function"
em_carbontrade(countries)
em_internationalTrade
eo_total
;

Positive Variables 
d_CarbNeutr(countries, producers)								"Production is limited by carbon capacity"
d_capacity(countries, producers, prodTypes) 		"Production capacity limits"
;

Binary Variables
b_CarbNeutr(countries, producers)								"Production is limited by carbon capacity"
b_capacity(countries, producers, prodTypes) 		"Production capacity limits"
;

ep_CarbNeutr(countries, producers).. 
	CARBON_BUY(countries, producers)  - CARBON_SELL(countries, producers) + InitCarb(countries, producers) - sum(pt, CarbonCost(pt)*PRODUCTION(countries, producers, pt)) =g= 0;

ep_capacity(countries, producers, prodTypes).. 
	InitCapacity(countries, producers, prodTypes) + INVESTMENT(countries, producers)$SAMEAS(prodTypes, "renew") -  PRODUCTION(countries, producers, prodTypes) =g= 0;

em_carbontrade(countries)..
sum(pp, CARBON_BUY(countries, pp) - CARBON_SELL(countries, pp)) =e= 0;

em_internationalTrade..

eo_producer(countries, producers)..
OBJ_CP(countries, producers) =g= 
	sum(prodTypes, (LinCost(countries, prodTypes) + 0.5*QuadCost(countries, prodTypes)*PRODUCTION(countries, producers, prodTypes) )*PRODUCTION(countries, producers, prodTypes)) 
+ (LinInvCost(countries) + QuadInvCost(countries)*INVESTMENT(countries, producers))*INVESTMENT(countries, producers)
- (DemInt(countries) - DemSlope(countries)*sum((pp,pt), PRODUCTION(countries, pp, pt)))*(sum(pt, PRODUCTION(countries, producers, pt)))
+ (CARBON_BUY(countries, producers) - CARBON_SELL(countries, producers))*CARBPRICE(countries)
;

eo_total..
OBJ =g= sum((cc, pp), OBJ_CP(cc, pp));

CARBPRICE.lo(countries) = 600;
CARBPRICE.up(countries) = 600;

Model TotCostMin
/
/
;
Solve TotCostMin min OBJ use NLP;
Display LinCost, QuadCost, LinInvCost, QuadInvCost, InitCapacity, InitCarb, DemInt, DemSlope;
Display CARBON_BUY.L;
Display CARBON_SELL.L;
Display PRODUCTION.L;
Display INVESTMENT.L;

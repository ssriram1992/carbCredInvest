from EPECinterface import *
from gurobipy import *
import itertools
import csv

timeLimitPerJob = 60

dirty = ["coal", "gas"]
green = ["wind", "solar"]

countries = ["c1", "c2"]
producers = ["p1_1", "p1_2", "p2_1", "p2_2"]
domesticity = tuplelist([("c1","p1_1"), ("c1","p1_2"),
          ("c2","p2_1"), ("c2","p2_2")
         ])


scenario = ["s"+str(i+1) for i in range(2)]






def invSubsidy(data,scen):
    """
    Subsidy scenario
    
    PARAMETERS
    data: 
        inputData objecti
    scen:
        -1 or 0 or 1 implying low, medium or high subsidy
    """
    fact = (scen + 1)*1.0
    for xx in data.LinInvCost:
        data.LinInvCost[xx] = data.LinInvCost[xx]*(1 - fact*5/100)
        data.QuadInvCost[xx] = data.LinInvCost[xx]*(1 - fact*5/100)
    return data

def demSlope(data, scen):
    """
    Demand Slope altering scenario
    
    PARAMETERS
    data: 
        inputData objecti
    scen:
        -1 or 0 or 1 implying low, medium or high demand Slope
    """
    fact = 0.1*scen
    for xx in data.DemSlope:
        data.DemSlope[xx] = data.DemSlope[xx]*(1 + fact)
        data.QuadInvCost[xx] = data.DemSlope[xx]*(1 + fact)
    return data

def stdDevnCapFac(data, scen):
    """
    Standard deviation of the capacity factor in scenarios
    
    PARAMETERS
    data: 
        inputData objecti
    scen:
        -1 or 0 or 1 implying low, medium or high capacity factors
    
    WARNING
        Only works for 2 scenario programs
    """
    fact = 0.1*(2+scen)
    for pp in data.producers:
        data.CapacityFactor[pp, 'wind','s1'] = 0.5 - fact
        data.CapacityFactor[pp, 'solar','s1'] = 0.5 + fact
        data.CapacityFactor[pp, 'wind','s2'] = 0.5 + fact
        data.CapacityFactor[pp, 'solar','s2'] = 0.5 - fact
    return data

def emmInv(data, scen, cc):
    fact = scen + 1
    data.EmissionValue[cc] = 10 + fact*40
    data.EmissionValueQuad[cc] = 1 + fact*2
    data.InvestValue[cc, "solar"] =25 + 50*fact 
    data.InvestValue[cc, "wind"] = 5 + 20*fact
    return data

def emmInvc1(data, scen):
    return emmInv(data, scen, "c1")

def emmInvc2(data, scen):
    return emmInv(data, scen, "c2")

def crossCoeff(data, scen, cc):
    if scen == -1:
        data.EmissionCrossValue[cc] = 0
        for ee in data.green:
            data.InvestCrossValue[cc,ee] = 0
        return data
    data.EmissionCrossValue[cc] = data.EmissionValueQuad[cc]*(1+scen)
    data.InvestCrossValue[cc,"wind"] = data.InvestValue[cc,"wind"]*(1+scen)/200
    data.InvestCrossValue[cc,"solar"] = data.InvestValue[cc,"solar"]*(1+scen)/40
    return data
    
def crossCoeffc1(data, scen):
    return crossCoeff(data, scen, "c1")

def crossCoeffc2(data, scen):
    return crossCoeff(data, scen, "c2")


def carbCredits(data, val): 
    """
    Number of carbon credits allotted to each producer in the country
    """
    for pp in data.producers:
        data.InitCredits[pp] = val 
    for cc in data.countries:
        data.InitLeaderCredits[cc] = 0
    return data




dataBase = [15 ,25 ,15 ,25 ,15 ,25 ,15 ,25 ,1 ,4 ,1 ,4 ,1 ,4 ,1 ,4 ,4 ,6 ,0 ,0 ,4 ,6 ,0 ,0 ,4 ,6 ,0 ,0 ,4 ,6 ,0 ,0 ,0.2 ,0.2 ,0 ,0 ,0.2 ,0.2 ,0 ,0 ,0.2 ,0.2 ,0 ,0 ,0.2 ,0.2 ,0 ,0 ,0.4 ,0.6 ,0.6 ,0.4 ,0.4 ,0.6 ,0.6 ,0.4 ,0.4 ,0.6 ,0.6 ,0.4 ,0.4 ,0.6 ,0.6 ,0.4 ,200 ,200 ,20 ,20 ,200 ,200 ,20 ,20 ,200 ,200 ,20 ,20 ,200 ,200 ,20 ,20 ,100 ,100 ,100 ,100 ,3 ,2 ,0 ,0 ,800 ,1200 ,1200 ,800 ,1 ,1 ,1 ,1 ,200 ,200 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0]





# mylist = readInputFromFile("inputs.csv")


solution = []
with open("dat/log.txt", "a") as file:
	file.write( "  ".join( [vv for vv in  ("Investment_subsidy", "Demand_Slope", "capFac_std_dev", "emission-inv_val_for_c1", "emission-inv_val_for_c2", "crossterms_for_c1", "crossterms_for_c1") ]   ))
	file.write ("\n")
    
bad = []


# foundBad  = [21,47,60,99,101,114,118,119,127,128,129,133,145,161,164,182,208,214,228,236,254,260,269,276,290,296,312,323,339,354,361,370,371,375,377,390,393,399,411,417,424,425,435,438,449,457,458,465,484,488,503,528,533,546,554,555,556,561,627,646,659,668,682,683,700,708,710,714,717,722,728,740,744,770,786,795,797,798,799,807,820,821,851,852,856,860,884,887,892,902,910,924,928,935,948,965,970,971,991,992,1005,1024,1041,1043,1047,1050,1052,1086,1091,1100,1113,1118,1119,1130,1131,1133,1135,1141,1153,1176,1195,1209,1248,1258,1266,1275,1281,1282,1284,1286,1294,1295,1315,1329,1343,1348,1365,1366,1373,1375,1394,1397,1419,1443,1447,1448,1451,1453,1455,1478,1505,1508,1511,1518,1523,1532,1534,1536,1576,1580,1586,1590,1610,1613,1617,1640,1658,1661,1666,1669,1673,1686,1689,1699,1700,1712,1745,1748,1751,1768,1770,1774,1777,1779,1781,1793,1802,1823,1828,1829,1834,1842,1851,1856,1859,1862,1883,1888,1889,1901,1913,1914,1924,1934,1943,1956,1966,1977,1997,2021,2044,2045,2063,2066,2069,2072,2076,2085,2086,2126,2139,2156,2157,2159]


count = 0
for invS, demS, sdCF, eI1, eI2, cT1, cT2, carbCre in itertools.product([-1,0,1],[-1,0,1],[-1,0,1],[-1,0,1],[-1,0,1], [-1,0,1],[-1,0,1], [5*i for i in range(21)]):
    if count not in foundBad:
        # pass
        count = count+1
        continue
    print(count, end=" ", flush=True) 
    data = inputData(dirty, green, countries, producers, domesticity, scenario)
    baseData = dataBase.copy()
    data.buildFromList(baseData)
    # Varying the data  for each scenario
    data = invSubsidy(data, invS)
    data = demSlope(data, demS)
    data = stdDevnCapFac(data, sdCF)
    data = emmInvc1(data, eI1)
    data = emmInvc2(data, eI2)
    data = crossCoeffc1(data, cT1)
    data = crossCoeffc2(data, cT2)
    data = carbCredits(data, carbCre)
    # Make the MPEC
    MPEC = makeMPEC(data)
    # for cc in MPEC.data.countries: 
        # MPEC.CARB_IMP[cc].lb = 0
        # MPEC.CARB_IMP[cc].ub = 0
    MPEC.Model.Params.TimeLimit = timeLimitPerJob 
    MPEC.Model.Params.Threads = 4
    MPEC.Model.optimize()
    try:
        # eqn = MPEC.Model.getObjective()
        # MPEC.Model.addConstr(eqn <= MPEC.Model.ObjVal, name = "equatingObjective")
        # MPEC.Model.setObjective((MPEC.CARBON_PRICE['c1'] - MPEC.CARBON_PRICE['c2'])*(MPEC.CARBON_PRICE['c1'] - MPEC.CARBON_PRICE['c2']))
        # MPEC.Model.setObjective(quicksum(
                   # MPEC.CARB_IMP[cc]*MPEC.CARB_IMP[cc]*1e-3 for cc in MPEC.data.countries))
        # MPEC.Model.optimize()
        ans = MPEC.listify()
        # Write individual file
        # with open("dat/NoCoop_"+str(count)+".txt", "w") as file:
            # file.write(",".join([str (vv) for vv in ans]))
        
        # Write log file
        with open("dat/logCoop2.txt", "a") as file:
            file.write(str(count)+": ")
            file.write( "  ".join( [str(vv) for vv in  (invS, demS, sdCF, eI1, eI2, cT1, cT2, carbCre) ]   ))
            file.write ("\n")
        solution.append( [ans] )
    except :
        print ("\n Infeasibility for:" , count,": ", (invS, demS, sdCF, eI1, eI2, cT1, cT2, carbCre),"")
        bad.append(count)
    # if count == 2:
        # break
    count = count+1

with open("dat/badCoop2.txt","a") as file:
    towrite = ",".join ([str(bb) for bb in bad])
    file.write(towrite)

with open("dat/solutionsCoop2.txt","w") as file:
    for vals in solution:
        file.write(",".join([str(vv) for vv in vals]))
        file.write('\n')

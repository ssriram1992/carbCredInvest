from gurobipy import *
import numpy as np
import csv

def build_tupledict(array,iterator):
    """
    Assign a value to a tupledict from a list
    Keys are generated by the iterator
    """
    ans = tupledict()
    for val,it in zip(array,iterator):
        ans[it] = val
    return ans


                                        ######
    #    #    #  #####   #    #   ##### #     #    ##     #####    ##
    #    ##   #  #    #  #    #     #   #     #   #  #      #     #  #
    #    # #  #  #    #  #    #     #   #     #  #    #     #    #    #
    #    #  # #  #####   #    #     #   #     #  ######     #    ######
    #    #   ##  #       #    #     #   #     #  #    #     #    #    #
    #    #    #  #        ####      #   ######   #    #     #    #    #


    

class inputData:
    """
    Input date for EPEC
    """
    def __init__(self, dirty, green, countries, producers, domesticity, scenario):
        self.dirty = dirty
        self.green = green
        self.energy = dirty+green
        self.countries = countries
        self.producers = producers
        self.domesticity = domesticity
        self.scenario = scenario
    def listify(self):
        """
        Creates a list in the canonical order of 
        all input fields. 

        WARNING
        --------
        Works only if the object has all the necessary fields. 
        """
        ans = []

        # LinInvCost
        ans = ans + [(self.LinInvCost[pp,ee]) for pp in self.producers for ee in self.green]
        # QuadInvCost
        ans = ans + [(self.QuadInvCost[pp,ee]) for pp in self.producers for ee in self.green]
        # LinProdCost
        ans = ans + [self.LinProdCost[pp,ee] for pp in self.producers for ee in self.energy]
        # QuadProdCost
        ans = ans + [self.QuadProdCost[pp,ee] for pp in self.producers for ee in self.energy]
        # CapacityFactor
        CF = [self.CapacityFactor[pp,ee,ss] for  pp in self.producers for ee in self.green for ss in self.scenario]
        ans = ans + CF
        # InitCapacity
        ans = ans + [self.InitCapacity[pp,ee] for pp in self.producers for ee in self.energy]
        # InitCredits
        ans = ans + [self.InitCredits[pp] for pp in self.producers]
        # Emission
        ans = ans + [self.Emission[ee] for ee in self.energy]
        # DemInt
        ans = ans + [self.DemInt[cc,ss] for cc in self.countries for ss in self.scenario]
        # DemSlope 
        ans = ans + [self.DemSlope[cc,ss] for cc in self.countries for ss in self.scenario]
        # InitLeaderCredits
        ans = ans + [self.InitLeaderCredits[cc] for cc in self.countries]
        # minCons
        ans = ans + [self.minCons[cc] for cc in self.countries]
        # InvestValue
        ans = ans + [self.InvestValue[cc,ee] for cc in self.countries for ee in self.green]
        # InvestCrossValue
        ans = ans + [self.InvestCrossValue[cc,ee] for cc in self.countries for ee in self.green]
        # EmissionValue
        ans = ans + [self.EmissionValue[cc] for cc in self.countries]
        # EmissionValueQuad
        ans = ans + [self.EmissionValueQuad[cc] for cc in self.countries]
        # EmissionCrossValue
        ans = ans + [self.EmissionCrossValue[cc] for cc in self.countries]
        #
        ### RETURN ###
        #
        return ans
    def buildFromList(self,ans):
        """
        Build the fields of inputData using a given list
        """
        nco = len(self.countries)
        npr = len(self.producers)
        ndi = len(self.dirty)
        ngr = len(self.green)
        nen = len(self.energy)
        nsc = len(self.scenario)

        counter = 0
        # LinInvCost
        arr = ans [counter:counter + npr*ngr]
        counter += (npr*ngr)
        self.LinInvCost = build_tupledict(arr, itertools.product(self.producers, self.green))
        # QuadInvCost
        arr = ans [counter:counter + npr*ngr]
        counter += (npr*ngr)
        self.QuadInvCost = build_tupledict(arr, itertools.product(self.producers, self.green))
        # LinProdCost
        arr = ans [counter:counter + npr*nen]
        counter += (npr*nen)
        self.LinProdCost = build_tupledict(arr, itertools.product(self.producers, self.energy))
        # QuadProdCost
        arr = ans [counter:counter + npr*nen]
        counter += (npr*nen)
        self.QuadProdCost = build_tupledict(arr, itertools.product(self.producers, self.energy))
        # CapacityFactor
        arr = ans [counter:counter + npr*ngr*nsc]
        counter += (npr*ngr*nsc)
        self.CapacityFactor  = build_tupledict(arr, itertools.product(self.producers, self.green, self.scenario))
        # InitCapacity
        arr = ans [counter:counter + npr*nen]
        counter += (npr*nen)
        self.InitCapacity = build_tupledict(arr, itertools.product(self.producers, self.energy))
        # InitCredits
        arr = ans [counter:counter + npr]
        counter += (npr)
        self.InitCredits = build_tupledict(arr, self.producers)
        # Emission
        arr = ans [counter:counter + nen]
        counter += (nen)
        self.Emission = build_tupledict(arr, self.energy)
        # DemInt
        arr = ans [counter:counter + nco*nsc]
        counter += (nco*nsc)
        self.DemInt = build_tupledict(arr, itertools.product(self.countries, self.scenario))
        # DemSlope
        arr = ans [counter:counter + nco*nsc]
        counter += (nco*nsc)
        self.DemSlope = build_tupledict(arr, itertools.product(self.countries, self.scenario))
        # InitLeaderCredits
        arr = ans [counter:counter + nco]
        counter += (nco)
        self.InitLeaderCredits = build_tupledict(arr, self.countries)
        # minCons
        arr = ans [counter:counter + nco]
        counter += (nco)
        self.minCons = build_tupledict(arr, self.countries)
        # InvestValue
        arr = ans [counter:counter + nco*ngr]
        counter += (nco*ngr)
        self.InvestValue = build_tupledict(arr, itertools.product(self.countries, self.green))
        # InvestCrossValue
        arr = ans [counter:counter + nco*ngr]
        counter += (nco*ngr)
        self.InvestCrossValue = build_tupledict(arr, itertools.product(self.countries, self.green))
        # EmissionValue
        arr = ans [counter:counter + nco]
        counter += (nco)
        self.EmissionValue = build_tupledict(arr, self.countries)
        # EmissionValueQuad
        arr = ans [counter:counter + nco]
        counter += (nco)
        self.EmissionValueQuad = build_tupledict(arr, self.countries)
        # EmissionCrossValue
        arr = ans [counter:counter + nco]
        counter += (nco)
        self.EmissionCrossValue = build_tupledict(arr, self.countries)
        if counter >= len(ans):
            self.probability = {ss:1/len(self.scenario) for ss in self.scenario}
        else:
            arr = ans [counter:counter + nsc]
            counter += (nsc)
            self.probability = build_tupledict(arr, self.scenario)
        return 



    ####### ######  #######  #####  #     #
    #       #     # #       #     # ##   ##   ####   #####   ######  #
    #       #     # #       #       # # # #  #    #  #    #  #       #
    #####   ######  #####   #       #  #  #  #    #  #    #  #####   #
    #       #       #       #       #     #  #    #  #    #  #       #
    #       #       #       #     # #     #  #    #  #    #  #       #
    ####### #       #######  #####  #     #   ####   #####   ######  ######




        
class EPECModel:
    """
    Contains input data and the model to beautifully print the output
    """
    def __init__(self, data, gurobiModel, Production, Investment, CarbonPrice, 
                 FollCarbBuy, CarbImp, LeadCarbBuy):
        self.data = data
        self.Model = gurobiModel
        self.PRODUCTION = Production
        self.INVESTMENT = Investment
        self.CARBON_PRICE = CarbonPrice
        self.FOLL_CARB_BUY = FollCarbBuy
        self.CARB_IMP = CarbImp
        self.LEAD_CARB_BUY = LeadCarbBuy
    
    def printLeader(self,cc):
        print()
        print("********************************************************")
        print("********************************************************")
        print("COUNTRY ",cc)
        print("********************************************************")
        print("********************************************************")
        print()
        print("Carbon price\t:\t",self.CARBON_PRICE[cc].X)
        print("Carbon buyback\t:\t", self.LEAD_CARB_BUY[cc].X)
        print("Carbon import\t:\t", self.CARB_IMP[cc].X)
        print("")
        
    def printFollower(self, cc, pp):
        if (cc,pp) not in self.data.domesticity:
            raise KeyError("Invalid country-producer pair: " + str( (cc,pp)))
            return
        print("****************************")
        print("Producer",pp,"of country",cc)
        print("****************************")
        print("*** INVESTMENTS ***")
#         print(self.INVESTMENT)
        for ee in self.data.green:
            print(ee,"\t:\t",self.INVESTMENT[pp,ee].X)
        print("*** CARBON CREDITS ***")
        print("Credits bought\t:\t", self.FOLL_CARB_BUY[pp].X)
        print("*** PRODUCTION ***")
        expec = {ee:0 for ee in self.data.energy}
        for ss in self.data.scenario:
            self.printScenarioDetail(ss)
            print("TOTAL:",sum([self.PRODUCTION[pp,ee,ss].X for ee in self.data.energy]))
            for ee in self.data.energy:
                prod = self.PRODUCTION[pp,ee,ss].X
                print(ee,"\t:\t",prod)
                if hasattr(self.data, "probability"):
                    expec[ee] += prod*self.data.probability[ss]
        if hasattr(self.data, "probability"):
            print("EXPECTED PRODUCTION: ","\t:\t",sum(expec[ee] for ee in expec))
            for ee in self.data.energy:
                print(ee,"\t:\t", expec[ee])
        return 
    def printScenarioDetail(self,ss):
        if ss not in self.data.scenario:
            raise KeyError("Scenario "+str(ss)+" invalid!")
            return
        print("Scenario: ", ss, end=' ')
        if hasattr(self.data, "probability"):
                print("with probability",self.data.probability[ss])
        else:
            print()
        return
    def printFollowers(self, cc, withAggr = False):
        if withAggr:
            self.printFollowersAggregate(cc)
        for pp in self.data.producers:
            if(cc, pp) in self.data.domesticity:
                self.printFollower(cc, pp)
        return
    def printFollowersAggregate(self, cc, withDetail = False):
        TotalProd = {ss:0 for ss in self.data.scenario}
        TotalEmit = {ss:0 for ss in self.data.scenario}
        TotalInv = {ee:0 for ee in self.data.green}
        print("****************************")
        print("***PRODUCER AGGREGATED DETAILS***")
        print("****************************")
        for pp in self.data.producers:
            if (cc,pp) in self.data.domesticity:
                for ee in self.data.green:
                    TotalInv[ee] += self.INVESTMENT[pp, ee].X
        for ee in self.data.green:
            print("Total investments",ee,"\t:\t",TotalInv[ee])
        
        print()
        
        for ss in self.data.scenario:
            for pp in self.data.producers:
                if(cc, pp) in self.data.domesticity:
                    TotalProd[ss] += sum([self.PRODUCTION[pp,ee,ss].X for ee in self.data.energy])
            self.printScenarioDetail(ss)
            print("Production\t:\t", TotalProd[ss])
            TotalEmit[ss] = self.Model.getVarByName("TOTAL_EMIT["+cc+","+ss+"]").X
            print ("Emission\t:\t", TotalEmit[ss])
        if hasattr(self.data, "probability"):
            print()
            print("Expected Total Production\t:\t", 
                  sum([TotalProd[ss]*self.data.probability[ss] for ss in self.data.scenario])
                 )
            print("Expected Total Emission\t:\t", 
                  sum([TotalEmit[ss]*self.data.probability[ss] for ss in self.data.scenario])
                 )
        print()
        if withDetail:
            self.printFollowers(cc)
        return
    def listify(self, dataFrame = False):
        """
        Returns a list of answers
        """
        ans = []
        name = []
        coun = []
        prod = []
        ener = []
        scen = []

        # PRODUCTION
        ans = ans + [self.PRODUCTION[pp,ee,ss].X for pp,ee,ss in itertools.product(self.data.producers,self.data.energy,self.data.scenario)]
        name = name + ["PRODUCTION" for pp,ee,ss in itertools.product(self.data.producers,self.data.energy,self.data.scenario)]
        coun = coun + ["" for pp,ee,ss in itertools.product(self.data.producers,self.data.energy,self.data.scenario)]
        prod = prod + [pp for pp,ee,ss in itertools.product(self.data.producers,self.data.energy,self.data.scenario)]
        ener = ener + [ee for pp,ee,ss in itertools.product(self.data.producers,self.data.energy,self.data.scenario)]
        scen = scen + [ss for pp,ee,ss in itertools.product(self.data.producers,self.data.energy,self.data.scenario)]
        # INVESTMENT
        ans = ans + [self.INVESTMENT[pp,ee].X for pp,ee in itertools.product(self.data.producers, self.data.green)]
        name = name + ["INVESTMENT" for pp,ee in itertools.product(self.data.producers,self.data.green)]
        coun = coun + ["" for pp,ee in itertools.product(self.data.producers,self.data.green)]
        prod = prod + [pp for pp,ee in itertools.product(self.data.producers,self.data.green)]
        ener = ener + [ee for pp,ee in itertools.product(self.data.producers,self.data.green)]
        scen = scen + ["" for pp,ee in itertools.product(self.data.producers,self.data.green)]
        # FOLL_CARB_BUY
        ans = ans + [self.FOLL_CARB_BUY[pp].X for pp in self.data.producers]
        name = name + ["FOLL_CARB_BUY" for pp in itertools.product(self.data.producers)]
        coun = coun + ["" for pp in (self.data.producers)]
        prod = prod + [pp for pp in (self.data.producers)]
        ener = ener + ["" for pp in (self.data.producers)]
        scen = scen + ["" for pp in (self.data.producers)]
        # CARBON_PRICE
        ans = ans + [self.CARBON_PRICE[cc].X for cc in self.data.countries]
        name = name + ["CARBON_PRICE" for cc in itertools.product(self.data.countries)]
        coun = coun + [cc for cc in (self.data.countries)]
        prod = prod + ["" for cc in (self.data.countries)]
        ener = ener + ["" for cc in (self.data.countries)]
        scen = scen + ["" for cc in (self.data.countries)]
        # CARB_IMP
        ans = ans + [self.CARB_IMP[cc].X for cc in self.data.countries]    
        name = name + ["CARB_IMP" for cc in itertools.product(self.data.countries)]
        coun = coun + [cc for cc in (self.data.countries)]
        prod = prod + ["" for cc in (self.data.countries)]
        ener = ener + ["" for cc in (self.data.countries)]
        scen = scen + ["" for cc in (self.data.countries)]
        # LEAD_CARB_BUY
        ans = ans + [self.LEAD_CARB_BUY[cc].X for cc in self.data.countries]
        name = name + ["LEAD_CARB_BUY" for cc in itertools.product(self.data.countries)]
        coun = coun + [cc for cc in (self.data.countries)]
        prod = prod + ["" for cc in (self.data.countries)]
        ener = ener + ["" for cc in (self.data.countries)]
        scen = scen + ["" for cc in (self.data.countries)]
        # COUNTRY_OBJ
        ans = ans + [self.Model.getVarByName("COUNTRY_OBJ["+cc+"]").X for cc in self.data.countries]
        name = name + ["COUNTRY_OBJ" for cc in itertools.product(self.data.countries)]
        coun = coun + [cc for cc in (self.data.countries)]
        prod = prod + ["" for cc in (self.data.countries)]
        ener = ener + ["" for cc in (self.data.countries)]
        scen = scen + ["" for cc in (self.data.countries)]

        if dataFrame:
            return ans, name, coun, prod, ener, scen
        else:
            return ans
    def comprehensivePrint(details = True):
        """
        Comprehensively prints the solutions.
        details: Bool
            If this is true, prints a more detailed solution for each 
            follower, if not aggregates them a bit.
        """
        for cc in data.countries:
            self.printLeader(cc)
            self.printFollowersAggregate(cc, details)

                                #     # ######  #######  #####
 #    #    ##    #    #  ###### ##   ## #     # #       #     #
 ##  ##   #  #   #   #   #      # # # # #     # #       #
 # ## #  #    #  ####    #####  #  #  # ######  #####   #
 #    #  ######  #  #    #      #     # #       #       #
 #    #  #    #  #   #   #      #     # #       #       #     #
 #    #  #    #  #    #  ###### #     # #       #######  #####




def  makeMPEC(data):
    """
    data: Instance of class inputData. Should have the input data with all the parameter fields populated. 
    RETURNS
    -------
    An object of type EPECModel with unsolved model. The Model hence obtained has to be optimized first,
    before other read/write option can be performed
    EXAMPLE
    ------- 
    myModel = makeEPEC(data)
    myModel.Model.optimize()
    for cc in myModel.countries:
        myModel.printLeader(cc)
    ans = myModel.listify()

    """
    M = Model()
    # ## Variables
    # ### Follower Variables 
    PRODUCTION = M.addVars(data.producers, data.energy, data.scenario, 
                           vtype=GRB.CONTINUOUS, name = "PROD")
    INVESTMENT = M.addVars(data.producers, data.energy, 
                           vtype=GRB.CONTINUOUS, name = "INV") 
    FOLL_CARB_BUY = M.addVars(data.producers, lb = -GRB.INFINITY, 
                                  name = "FOLL_CARB_BUY")

    ### All data.dirty investments must be 0 
    for ee in data.dirty:
        for x in INVESTMENT.select("*",ee):
            x.ub = 0


    # ### Leader Variables 
    CARBON_PRICE = M.addVars(data.countries, name="CARB_PRI")
    C_PROD = M.addVars(data.countries, data.scenario, name="C_PROD")
    CARB_IMP = M.addVars(data.countries, lb=-GRB.INFINITY, name = "CARB_IMP")

    TOTAL_INV = M.addVars(data.countries, data.energy, name = "TOTAL_INV")
    TOTAL_EMIT = M.addVars(data.countries, data.scenario, name = "TOTAL_EMIT")
    LEAD_CARB_BUY = M.addVars(data.countries, lb = -GRB.INFINITY, name = "LEAD_CARB_BUY") 

    COUNTRY_OBJ = M.addVars(data.countries, lb= -10000+0*GRB.INFINITY, name = "COUNTRY_OBJ") 

    # ### Duals 
    D_INFRALIMIT = M.addVars(data.producers, data.energy, data.scenario,
                            name = "D_INFRA")
    D_EMITLIMIT = M.addVars(data.producers, data.scenario, name = "D_EMIT")

    # ### Binary variables 
    B_PRODUCTION= M.addVars(data.producers, data.energy, data.scenario, 
                           vtype= GRB.BINARY, name = "B_Prod")
    B_INVESTMENT= M.addVars(data.producers, data.energy, 
                           vtype= GRB.BINARY, name = "B_data.Invest")
    B_INFRALIM = M.addVars(data.producers, data.energy, data.scenario, 
                          vtype = GRB. BINARY, name = "B_InfraLim")
    B_EMITLIM = M.addVars(data.producers, data.scenario,  
                          vtype = GRB. BINARY, name = "B_EmitLim")


    # Constraints
    # ### Follower KKT conditions 

    # eq_infraLimit
    for pp,ee,ss in itertools.product(data.producers,data.energy,data.scenario):
        eqn = 0
        if ee in data.green:
            eqn = data.CapacityFactor[pp, ee, ss]*(INVESTMENT[pp,ee]+data.InitCapacity[pp, ee]) - PRODUCTION[pp, ee, ss]
        else:
            eqn = INVESTMENT[pp,ee]+data.InitCapacity[pp, ee]- PRODUCTION[pp, ee, ss]
        M.addConstr(eqn >= 0,name = "eq_infraLimit"+"_".join([pp,ee,ss]))
        M.addGenConstrIndicator(B_INFRALIM[pp,ee,ss],0,eqn ==0, name = "eq_infraLimit1"+"_".join([pp,ee,ss]))
        M.addGenConstrIndicator(B_INFRALIM[pp,ee,ss],1,D_INFRALIMIT[pp,ee,ss] ==0, name = "eq_infraLimit2"+"_".join([pp,ee,ss]))

    # eq_emitLimit
    for pp,ss in itertools.product(data.producers, data.scenario):
        eqn = data.InitCredits[pp] + FOLL_CARB_BUY[pp] -     quicksum(data.Emission[ee]*PRODUCTION[pp,ee,ss] for ee in data.dirty)
        M.addConstr(eqn >= 0,name = "eq_emitLimit"+"_".join([pp,ss]))
        M.addGenConstrIndicator(B_EMITLIM[pp,ss],0,eqn ==0, name = "eq_emitLimit1"+"_".join([pp,ss]))
        M.addGenConstrIndicator(B_EMITLIM[pp,ss],1,D_EMITLIMIT[pp,ss] ==0, name = "eq_emitLimit2"+"_".join([pp,ss]))
        
    # eq_countryProduction
    for ss in data.scenario:
        for cc in data.countries:
            eqn = 0
            for pp in data.producers:
                if (cc,pp) in data.domesticity:
                    eqn = eqn + quicksum(PRODUCTION[pp,ee,ss] for ee  in data.energy)
            M.addConstr(eqn == C_PROD[cc,ss], name="eq_countryProduction"+"_".join([cc,ss]))
                    
    # eq_production
    for pp,ee,ss in itertools.product(data.producers,data.energy,data.scenario):
        eqn = data.probability[ss]*(
            data.LinProdCost[pp,ee] + data.QuadProdCost[pp,ee]*PRODUCTION[pp,ee,ss] -
            quicksum(data.DemInt[cc[0],ss] - data.DemSlope[cc[0],ss]*C_PROD[cc[0],ss] for cc in data.domesticity.select("*",pp)) +
            quicksum(data.DemSlope[cc[0],ss]*PRODUCTION[pp,ee2,ss] for cc in data.domesticity.select("*",pp) for ee2 in data.energy) 
        )+  D_INFRALIMIT[pp,ee,ss] + data.Emission[ee]*D_EMITLIMIT[pp,ss]
        
        M.addConstr(eqn >= 0,name = "eq_PRODUCTION"+"_".join([pp,ee,ss]))
        M.addGenConstrIndicator(B_PRODUCTION[pp,ee,ss],0,eqn ==0, name = "eq_PRODUCTION1"+"_".join([pp,ee,ss]))
        M.addGenConstrIndicator(B_PRODUCTION[pp,ee,ss],1,PRODUCTION[pp,ee,ss] ==0, name = "eq_PRODUCTION2"+"_".join([pp,ee,ss]))
        
    # eq_investment
    for pp,ee in itertools.product(data.producers,data.green):
        eqn = data.LinInvCost[pp,ee] + data.QuadInvCost[pp,ee]*INVESTMENT[pp,ee] - quicksum(data.CapacityFactor[pp,ee,ss]*D_INFRALIMIT[pp,ee,ss] for ss in data.scenario)
        M.addConstr(eqn >= 0,name = "eq_INVESTMENT"+"_".join([pp,ee]))
        M.addGenConstrIndicator(B_INVESTMENT[pp,ee],0,eqn ==0, name = "eq_INVESTMENT1"+"_".join([pp,ee]))
        M.addGenConstrIndicator(B_INVESTMENT[pp,ee],1,INVESTMENT[pp,ee] ==0, name = "eq_INVESTMENT2"+"_".join([pp,ee]))
        
    # eq_carbonPurchase
    for pp in data.producers:
        cc = data.domesticity.select("*", pp)[0][0]
        M.addConstr(CARBON_PRICE[cc] -quicksum(D_EMITLIMIT[pp,ss] for ss in data.scenario) == 0, name = "eq_carbonPurchase"+pp
                   )


    # ### Leader Constraints 
    #     Mininmum consumption constraint
    for cc in data.countries:
        M.addConstrs((C_PROD[cc,ss] >= data.minCons[cc] for ss in data.scenario), name = "data.minCons"+cc) 

    # Total data.Emission
    for ss in data.scenario:
        eqn = {cc:0 for cc in data.countries}
        for (cc,pp) in data.domesticity:
            eqn[cc] = eqn[cc] + quicksum(data.Emission[ee]*PRODUCTION[pp,ee,ss] for ee in data.energy)
        M.addConstrs((TOTAL_EMIT[cc,ss] - eqn[cc]==0 for cc in data.countries), name="Totaldata.Emission"+str(ss))
        
    # Total data.Investment
    for ee in data.green:
        eqn = {cc:0 for cc in data.countries}
        for (cc,pp) in data.domesticity:
            eqn[cc] = eqn[cc] + INVESTMENT[pp,ee] 
        M.addConstrs((TOTAL_INV[cc,ee] - eqn[cc]==0 for cc in data.countries), name="TotalInv"+"_".join([cc,ee]))
        
    # Carbon Buy summing
    eqn = {cc:0 for cc in data.countries}
    for (cc,pp) in data.domesticity:
        eqn[cc] = eqn[cc] + FOLL_CARB_BUY[pp] 
    M.addConstrs((LEAD_CARB_BUY[cc] + eqn[cc] == 0 for cc in data.countries), name = "CarbBuySum")

    # Carbon Credits >= 0
    M.addConstrs((LEAD_CARB_BUY[cc] + data.InitLeaderCredits[cc] + CARB_IMP[cc]  >= 0 for cc in data.countries), 
                 name = "CarbCredPos")

    M.addConstr(quicksum(CARB_IMP[cc] for cc in data.countries) == 0)

    M.update()

    # ### Leader Objective 
    for cc in data.countries:
        M.addConstr(
           1000 * COUNTRY_OBJ[cc] == CARBON_PRICE[cc]*LEAD_CARB_BUY[cc] + # Buy credits from follower
                   quicksum(data.probability[ss]*TOTAL_EMIT[cc,ss]*data.EmissionValue[cc] for ss in data.scenario)- # self emission
                    quicksum(TOTAL_INV[cc,ee]*data.InvestValue[cc,ee] for ee in data.green) +  # self investment
                    quicksum(data.probability[ss]*TOTAL_EMIT[cc,ss]*TOTAL_EMIT[c2c,ss]*data.EmissionCrossValue[cc]
                             for c2c in data.countries for ss in data.scenario
                            )- # Emission cross
                   quicksum(TOTAL_INV[cc,ee]*TOTAL_INV[c2c,ee]*data.InvestCrossValue[cc,ee] 
                            for ee in data.green for c2c in data.countries), # cross investment    ,
            name = "eq_countryObj"+cc)
        
    activeCountry = {cc:1 for cc in data.countries}
    M.setObjective(quicksum(activeCountry[cc]*COUNTRY_OBJ[cc] for cc in data.countries))
    
    M.Params.LogToConsole = 0
    M.Params.NonConvex=2
    M.Params.Presolve = 2
    M.Params.CutPasses = 1
    M.Params.PreQLinearize = 1
    M.update()
    
    MPEC = EPECModel(data, M, PRODUCTION, INVESTMENT, CARBON_PRICE, FOLL_CARB_BUY, CARB_IMP, LEAD_CARB_BUY)
    return MPEC

def readInputFromFile(filename, delim = ','):
    """
    Reads inputs from a csv file
    and returns it in a list
    PARAMETERS
    ----------
    filename: str
        Complete filename with extension
    delim: str
        Delimiting character separting multiple values in the same row
    
    RETURNS
    -------
    mylist
    """
    mylist = []

    with open(filename,"r") as myfile:
        csv_read =  csv.reader(myfile,delimiter=delim)
        for lines in csv_read:
            numlines = [float(l) for l in lines]
            mylist = mylist+[numlines]
    return mylist

def solveMPECs(sets, mylist, timeLimitPerJob = 60):
    """
    Given a list of lists, solve the MPEC for each list.
    """
    solution = []
    count = 0
    length = len(mylist)
    dirty = sets.dirty
    green = sets.green
    countries = sets.countries
    producers = sets.producers
    domesticity = sets.domesticity
    scenario = sets.scenario
    for inp in mylist:
        count = count+1
        print("Processing",count,"out of", length)
        data = inputData(dirty, green, countries, producers, domesticity, scenario)
        data.buildFromList(inp)
        MPEC = makeMPEC(data)
        MPEC.Model.Params.TimeLimit = timeLimitPerJob 
        MPEC.Model.Params.Threads = 4
        MPEC.Model.optimize()
        eqn = MPEC.Model.getObjective()
        MPEC.Model.addConstr(eqn <= MPEC.Model.ObjVal, name = "equatingObjective")
        MPEC.Model.setObjective(
					(MPEC.CARBON_PRICE['c1'] - MPEC.CARBON_PRICE['c2'])*(MPEC.CARBON_PRICE['c1'] - MPEC.CARBON_PRICE['c2'])
					)
        # MPEC.Model.setObjective(quicksum(
                    # MPEC.CARB_IMP[cc]*MPEC.CARB_IMP[cc]*1e-3 for cc in MPEC.data.countries
                # ))
        MPEC.Model.optimize()
        ans = MPEC.listify()
        solution.append(ans)
    return solution

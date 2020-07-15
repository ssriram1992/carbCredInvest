import sys
import re
import itertools
import gurobipy as gp

# value = "x_"+str(11)
def getSolFromFile(value, filename):
    with open(filename, "r") as fo:
        lines = fo.readlines()
    for line in lines:
        ans = re.split(" +",line) # Regular expression split. Treat consective delimiters as one
        v = ans[0]
        s = ans[1]
        if v==value:
            return( float(s))
    return 0.0

def getSolsFromFile(values, filename):
    """
    values = ["x_11", "x_33"]
    ans = getSolsFromFile(values, "dat/dbgSol.sol")
    print(ans)
    
    The above code fetches the value of variables x_11 and x_33 
    in a gurobi or a scip written solution file dat/dbgSol.sol
    and returns the values in a dictionary
    """    
    ans = {k:0.0 for k in values}
    with open(filename, "r") as fo:
        lines = fo.readlines()
    for line in lines:
        ind = re.split(" +",line) # Regular expression split. Treat consective delimiters as one
        v = ind[0]
        s = ind[1]
        if v in values:
            print(v)
            ans[v] = float(s)
    return ans

def getParLocVar(filename):
    """
    Get the set of parameters, variable locations in the lcp model and values
    from a given file written using EPEC::appendSolution4XL
    If the file has multiple solution value rows, the first of them is chosen
    """
    with open(filename, "r") as fo:
        line1 = fo.readline().split(" ")
        line2 = fo.readline().split(" ")
        line3 = fo.readline().split(" ")

    param = dict()
    var = dict()
    locs = dict()
    for w1, w2, w3 in zip(line1, line2, line3):
        if w2.strip() == '':break
        key = w1.split("_")
        if len(key)>1: key = tuple(key)
        else: key = key[0]
        if float(w2) == -1: # i.e., a Parameter
            param[key] = float(w3)
        else:
            locs[key] = int(w2)
            var[key] =  float(w3)
    return param, locs, var

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def formatIt(a,b = ""):
    if is_number(b):
        return "{:>30}: {:>8}\n".format(a,round(b,2))
    else:
        return "{:>30}: {:>8}\n".format(a,b)
    
    
def getGurVar(model, locationDict, variable):
    return model.getVarByName("x_"+str(locationDict[variable]))

########################################################################
class EPECio:
    def __init__(self, countries, producers, dirty, clean, n_Scen=1):
        """
        countries should be a list of countries
        producers should be a dictionary, where the keys are elements of countries. 
                producers[cc] should contain the set of producers in country cc
        dirty is the set of dirty energy technologies
        clean is the set of clean energy technologies
        """
        self.countries = countries
        self.producers = producers
        self.dirty = dirty
        self.clean = clean
        self.n_Scen = n_Scen
        self.energy = dirty+clean
    def acquireData(self,filename):
        self.param, self.locs, self.var = getParLocVar (filename)
    def writeFoll(self,filename, cc, ff):
        with open(filename,"a") as fo:
            fo.write(cc+"---"+ff+"\n")
            fo.write(formatIt("Carbon Limit",self.var[(cc,ff,"CarbBuy")]))
            fo.write("\nInvestments by "+ff+"\n")
            for gg in self.clean:
                fo.write(formatIt("Investment in "+gg,self.var[(cc,ff,"Inv",gg)]))
            fo.write("\n")
            for ee in self.energy:
                ans = self.param[(cc,ff,"Cap",ee)] + (self.var[(cc,ff,"Inv",ee)] if ee in self.clean else 0)
                fo.write(formatIt("Total Capacity in "+ee, ans))
            fo.write("\n")
            # Scenario wise production
            expProd = {ee:0 for ee in self.energy}
            for xi in range(self.n_Scen):
                if self.n_Scen>1: fo.write("Scenario: "+str(xi+1)+"\n")
                if self.n_Scen>1:fo.write("  Production: "+"\n")
                for ee in self.energy:
                    val = self.var[(cc,ff,"Prod",ee,"xi"+str(xi))]
                    if self.n_Scen>1: fo.write(formatIt(ee+": ",val))
                    expProd[ee] += val/self.n_Scen
                if self.n_Scen>1: fo.write("\n")
            # Expected production
            fo.write(" EXPECTED Production: "+"\n")
            for ee in self.energy:
                fo.write(formatIt(ee+": ",expProd[ee]))
            fo.write(formatIt("Total: ",sum(expProd[ee] for ee in self.energy)))
            fo.write("\n")
            fo.write("-----------------------------------------")
            fo.write("\n")
            return expProd
                
    def writeCoun(self,filename, cc):
        with open(filename,"a") as fo:
            fo.write(cc + ":\n")
            fo.write(formatIt("Carbon Import:",self.var[(cc,"CarbImp")]))
            fo.write(formatIt("Total Emission:",self.var[(cc,"TotEmission")]))
            fo.write(formatIt("Carbon Tax:",self.var[(cc,"CarbTax")]))
        Prod = dict()
        for ff in self.producers[cc]:
            Prod[ff] = self.writeFoll(filename,cc,ff)
        with open(filename,"a") as fo:
            for ee in self.energy:
                fo.write(formatIt("Country total "+ee+": ",sum(Prod[ff][ee] for ff in self.producers[cc])))
            production = sum(Prod[ff][ee] for ff in self.producers[cc] for ee in self.energy)
            fo.write(formatIt("Country Total Prodn: ", production))
            demInt = sum([self.param[cc,'DemInt','xi'+str(xi)] for xi in range(self.n_Scen)])/self.n_Scen
            demSl = sum([self.param[cc,'DemSlope','xi'+str(xi)] for xi in range(self.n_Scen)])/self.n_Scen
            # fo.write(formatIt("demInt: ", demInt))
            # fo.write(formatIt("demSlp: ", demSl))
            fo.write(formatIt("Energy price Estimate: ", demInt - demSl*production))

    def writeAll(self,filename):
        with open(filename,"w") as fo:
            fo.write("###########################################\n")
            carbSuppl = sum( [self.var [(cc, "CarbImp")] for cc in self.countries ])
            fo.write(formatIt("Total Carbon Supplied:", carbSuppl))
            fo.write(formatIt("Carbon Price:", str(self.param["suppInt"]) + " + " + str(self.param["suppSlope"]) + "Q"))
            fo.write(formatIt(" ", self.param["suppInt"] + self.param["suppSlope"]*carbSuppl ))
            fo.write("###########################################\n")
        for cc in self.countries:
            self.writeCoun(filename, cc)
            with open(filename,"a") as fo:
                fo.write("###########################################\n")
                
        
#########################################################################################

if __name__ == '__main__':
    """
    Usage readWrite [prefix='Main'] [reoptimize=True] [verbose=True]
    """
    if len(sys.argv) > 1:
        prefix = sys.argv[1]
    else:
        print ("No prefix provided in commandline. Use prefix: 'Main'")
        prefix = 'Main'
        
    positions = "./dat/"+prefix+"solLog.dat"
    lpfile = "./dat/"+prefix+"lcpmodel.lp"
    unmod = "./dat/"+prefix+"unMod.txt"
    mod = "./dat/"+prefix+"mod"
    
    if len(sys.argv) > 2:
        reoptimize = int(sys.argv[2])
    else:
        print("No second argument. So no reoptimizing...")
        reoptimize = False
        
    if len(sys.argv) > 3:
        verbose = int(sys.argv[3])
    else:
        print("No third argument. So not verbose...")
        verbose = False
    
    # Preliminary data
    countries = ["c1","c2"]
    producers = {"c1":["F1","F2"],"c2":["F3","F4"]}
    dirty = ["coal", "gas"]
    clean = ["solar", "wind"]
    
    # Making the EPECio object
    epec = EPECio(countries, producers, dirty, clean)
    epec.acquireData(positions)
    
    # Writing the output
    epec.writeAll(unmod)
    
    if reoptimize:
        # Solving the gurobi model
        M = gp.read(lpfile)
        M.params.LogToConsole = verbose
        expr = 0 # M.getObjective()
        expr += gp.quicksum([getGurVar(M, epec.locs, (cc,"CarbImp")) for cc in countries])*0.1
        M.params.SolFiles = "./dat/tempSol"
        M.setObjective(expr)
        print("Using the fact that the data is symmetric")
        ci1 = getGurVar(M, epec.locs, ('c1','CarbImp'))
        ci2 = getGurVar(M, epec.locs, ('c2','CarbImp'))
        # M.addConstr(ci1==ci2)
        M.update()
        M.optimize()
        M.write("./dat/_tempSol.sol")

        # Getting the solution from the solved model
        import os
        for sol in range (M.getAttr("SolCount")):
          newVars = dict()
          for kk in epec.locs.keys():
              vv = epec.locs[kk]
              newVars[kk] = getSolFromFile("x_"+str(vv),"./dat/tempSol_"+str(sol)+".sol")
          epec.var = newVars
          # Writing the output
          # os.remove("./dat/tempSol_"+str(sol)+".sol") 
          epec.writeAll(mod+str(sol)+".txt")
        # os.remove("./dat/_tempSol.sol") 
        # Writing the best output
        epec.writeAll(mod+".txt")
        print((M.getAttr("SolCount")), "solutions found and written to files")
    print("Task completed successfully")



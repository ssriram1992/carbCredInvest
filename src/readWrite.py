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
            for gg in ["solar","wind"]:
                fo.write(formatIt("Investment in "+gg,self.var[(cc,ff,"Inv",gg)]))
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
            fo.write(formatIt("Carbon Tax:",self.var[(cc,"CarbTax")]))
        Prod = dict()
        for ff in self.producers[cc]:
            Prod[ff] = self.writeFoll(filename,cc,ff)
        with open(filename,"a") as fo:
            for ee in self.energy:
                fo.write(formatIt("Country total "+ee+": ",sum(Prod[ff][ee] for ff in self.producers[cc])))
            fo.write(formatIt("Country Total Prodn: ",sum(Prod[ff][ee] for ff in self.producers[cc]
                                                          for ee in self.energy
                                                         )))

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
#     argv = ["readWrite", "../dat/dbgsolLog.txt", "../dat/posSol.sol", "../dat/dbgSol.sol"]
    argv = sys.argv
    if len(argv) != 4:
        print("Error! \n Format: readWrite dat/solutionLog.txt dat/lcpmodel.lp dat/outputName.txt\n")
        exit(1)
    
    # Preliminary data
    countries = ["c1","c2"]
    producers = {"c1":["F1","F2"],"c2":["F3","F4"]}
    dirty = ["coal", "gas"]
    clean = ["solar", "wind"]
    
    # Making the EPECio object
    epec = EPECio(countries, producers, dirty, clean)
    epec.acquireData(argv[1])
    
    # Solving the gurobi model
    M = gp.read(argv[2])
    expr = M.getObjective()
    expr += gp.quicksum([getGurVar(M, epec.locs, (cc,"CarbImp")) for cc in countries])*0.1
    M.setObjective(expr)
    M.optimize()
    M.write("./dat/_tempSol.sol")
    
    # Getting the solution from the solved model
    newVars = dict()
    for kk in epec.locs.keys():
        vv = epec.locs[kk]
        newVars[kk] = getSolFromFile("x_"+str(vv), "./dat/_tempSol.sol")
    epec.var = newVars
    
    # Writing the output
    epec.writeAll(argv[3])











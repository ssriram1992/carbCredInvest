from gurobipy import *
import numpy as np


def solveBiMatrix(A, B, costs = False, verbose=False):
    """
    Identifies a Nash equilibrium for a bimatrix Game using Gurobi 9.
    Formulation based on Chapter 1 in "The Linear complementarity 
    Problem by Cottle, Pang and Stone. Published by SIAM."
    
    PARAMETERS:
    -----------
    A: np.array 
       A[i,j] should be the payoff (or costs) for the first player
       if the first player plays strategy i and the second player 
       plays strategy j
    
    B: np.array 
       B[i,j] should be the payoff (or costs) for the first player
       if the first player plays strategy i and the second player 
       plays strategy j
       
    costs: Bool (default=False)
       If this is true, then the entries of A and B are considered as 
       costs as opposed to payoffs. As a result, the game is played
       as if the players are minimizing the expected costs (as opposed
       to maximizing the expected pay offs.)
       
    verbose: Bool (default = False)
       If this is true, prints the Gurobi Logs during solve. If not, prints nothing. 
         
    ASSERTS:
    --------
    A.shape == B.shape
    
    RETURNS:
    --------
    A_MNE:
        np.array of shape (m,) giving the probabilities with which each
        strategy should be played by the first player at a Nash equilibrium
    B_MNE:
        np.array of shape (n,) giving the probabilities with which each
        strategy should be played by the second player at a Nash equilibrium
    """
    assert A.shape == B.shape
    if not costs: # Internally, always solve it as costs
        A = -A
        B = -B
    
    temp = min (np.min(A), np.min(B))
    if temp <= 0.1:
        A = A + (1 - temp)
        B = B + (1 - temp)

    m = A.shape[0]
    n = A.shape[1]
    M = Model()
    if not verbose:
        M.params.LogToConsole = 0

    x = M.addMVar(m, lb = np.zeros(m), name = "x")
    y = M.addMVar(n, lb = np.zeros(n), name = "y")
    u = M.addMVar(m, lb = np.zeros(m), name = "u")
    v = M.addMVar(n, lb = np.zeros(n), name = "v")

    M.addConstr(A@y - np.ones(m) == u )
    M.addConstr((B.T)@x - np.ones(n) == v )
    
    M.addConstr(u@x ==0, name = "uPerpx")
    M.addConstr(v@y ==0, name = "vPerpy")
    
    M.params.NonConvex = 2
    M.optimize()

    return (x.X/sum(x.X), y.X/sum(y.X))

if __name__  == "__main__":
	if len(sys.argv) < 6:
		print('Usage: Bimatrix.py payoff1 payoff2 output costs verbose\nFor example: \npython Bimatrix.py A.csv B.csv out 0 1\n');
		quit()
	A = np.genfromtxt(sys.argv[1])
	B = np.genfromtxt(sys.argv[2])
	ans1, ans2 = solveBiMatrix(A, B, int(sys.argv[4]), int(sys.argv[5]))
	print(ans1, ans2)
	out = sys.argv[3]
	np.savetxt(out+"A.csv",ans1)
	np.savetxt(out+"B.csv",ans2)


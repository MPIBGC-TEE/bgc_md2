from copy import deepcopy 
class MDP():
    def __init__(self, Rs, Ts):
        self.Rs=Rs
        self.Ts=Ts
        self.actions={k: list(v.keys()) for k,v in Ts.items()}

        self.Qs={
            k: {vl: 0  for vl in v}
            for k,v in self.actions.items()
        }

    def update_bellman(self,gamma):
        Ts=self.Ts
        Rs=self.Rs
        actions=self.actions
        Qs_old=self.Qs
        def q(s,a): 
            target_states=Ts[s][a].keys()
            return sum(
                [
                    Ts[s][a][sp]*(
                        Rs[s][a][sp] + gamma * (
                            max(
                                [   
                                    Qs_old[sp][ap]
                                    for ap in actions[sp]
                                ]
                            ) if len(actions[sp])>0 else 0
                        )
                    )
                    for sp in target_states
                ]
            )
    
        Qs_new={
            s :{
                a: q(s,a) 
                for a in actions[s]
                #for a in v.keys()
            }
            for s,v in Qs_old.items()
        }
        self.Qs=Qs_new

    def update(self,Sample,gamma,alpha):
        s,a,sp,r=Sample
        Qs=deepcopy(self.Qs)
        self.Qs[s][a] = (
            (1-alpha)*Qs[s][a]
            +
            alpha*(
                r +
                gamma * (
                    max( Qs[sp].values()) if len(Qs[sp])>0 else 0
                )    
            )
        )



Ts_1 = { 
    "S":{"r":{"A": 1.0}},
    "A":{"r":{"E10": 1.0}, "u":{"E1": 1.0}, "l":{"S": 1.0}},
    "E1":{"ex":{"X": 1.0}},
    "E10":{"ex":{"X": 1.0}},
    "X":{}
}    
Rs_1 = { 
    "S":{"r":{"A": 0}},
    "A":{"r":{"E10": 0}, "u":{"E1": 0}, "l":{"S": 0}},
    "E1":{"ex":{"X": 1}},
    "E10":{"ex":{"X": 10}},
}    
mdp_1 = MDP(Rs_1, Ts_1)
#mdp_1.update_bellman(gamma=1)

S1_mdp1=[
        ("S","r","A",0),
        ("A","u","E1",0),
        ("E1","ex","X",1),
        ("S","r","A",0),
        ("A","r","E10",0),
        ("E10","ex","X",10),
]


mdp_1 = MDP(Rs_1, Ts_1)
for i in range(100):
    for ex in S1_mdp1:
        mdp_1.update(ex,gamma=1,alpha=0.5)
print(mdp_1.Qs["A"])
print(mdp_1.Qs["S"])

from IPython import embed; embed()


#S1_mdp2=[
#        ("S","r","A",0),
#        ("A","esc","E1",0),
#        ("E1","ex","X",1),
#        ("S","r","A",0),
#        ("A","esc","E10",0),
#        ("E10","ex","X",10),
#]
#
#S2_mdp2=[
#        ("S","r","A",0),
#        ("A","esc","E1",0),
#        ("E1","ex","X",1),
#        ("S","r","A",0),
#        ("A","esc","E10",0),
#        ("E10","ex","X",10),
#        ("S","r","A",0),
#        ("A","esc","E10",0),
#        ("E10","ex","X",10),
#]
#
##Ts_2 = { 
##    "S":{"r":{"A": 1.0}},
##    "A":{"esc":{"S": 1.0/3, "E1": 1.0/3, "E10": 1.0/3}},
##    "E1":{"ex":{"X": 1.0}},
##    "E10":{"ex":{"X": 1.0}}
##}    
##Rs_2 = { 
##    "S":{"r":{"A": 0}},
##    "A":{"esc":{"S": 0, "E1": 0, "E10": 0}},
##    "E1":{"ex":{"X": 1}},
##    "E10":{"ex":{"X": 10}}
##}    
##

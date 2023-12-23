from typing import Dict
from class_city import City
import gurobipy as gp
from gurobipy import GRB

class USApHMP:
    """This class implements the uncapacitated, single allocation p-hub median problem from
    ANDREAS T. ERNST and MOHAN KRISHNAMOORTHY (1996)
    """

    def __init__(self, cities_data: Dict[int, City], max_nodes:int = 81, number_hubs:int = 1):
        assert 2 <= max_nodes <= 81
        assert 1 <= number_hubs
        self.cities_data = cities_data
        
        # Define set and parameters
        self.set_cities = self.__get_set_cities__(max_nodes)    # N
        self.par_flow = self.__get_flow__()                     # Wij
        self.par_supply = self.__get_flow_supply__()            # Oi
        self.par_demand = self.__get_flow_demand__()            # Di
        self.par_distance = self.__get_distance__()             # dij
        self.par_unit_cost_collection = 1.0                     # χ
        self.par_unit_cost_transfer = 0.5                       # α
        self.par_unit_cost_distribution = 1.0                   # δ
        self.par_number_hubs = number_hubs                      # p
        
        # Create optimization model
        self.model = gp.Model('USApHMP')
        
        # Create variables
        self.var_Z = self.__add_var_Z__()
        self.var_Y = self.__add_var_Y__()

        # Create constraints
        self._add_constraint_locate_hubs()
        self._add_constraint_assign_hub_to_node()
        self._add_constraint_hub_implication()
        self._add_constraint_multicommodity_flow_equation()

        self.model.update()

        print(f"\n\n\n\n\n\n\nConstraints:")
        for con in self.model.getConstrs():
            print(f"\n{con}:\n{self.model.getRow(con)} {con.Sense} {con.RHS}")
        


    def __get_set_cities__(self, max_nodes: int) -> list:
        return list(cities_data.keys())[0: max_nodes]

    def __get_flow__(self) -> dict():
        flow = dict() # (i,j) -> flow_ij
        for i in self.set_cities:
            city_i = self.cities_data[i]
            for j, value_ij in city_i.flow_goods_to_other_cities.items():
                if j in self.set_cities:
                    flow[(i, j)] = int(value_ij)
        return flow
    
    def __get_flow_supply__(self) -> dict:
        flow_origin = dict() # supply of node i
        for i in self.set_cities:
            city_i = self.cities_data[i]
            total_supply = 0
            for j, value_ij in city_i.flow_goods_to_other_cities.items():
                if j in self.cities_data:
                    total_supply += value_ij
            flow_origin[i] = total_supply
        return flow_origin
    
    def __get_flow_demand__(self) -> dict:
        demand = dict() # demand of node i
        for i in self.set_cities:
            demand_i = 0
            for j in self.set_cities:
                city_j = self.cities_data[j]
                demand_i += city_j.flow_goods_to_other_cities[i]
            demand[i] = demand_i
        return demand

    def __get_distance__(self) -> dict():
        distance = dict() # (i,j) -> flow_ij
        for i in self.set_cities:
            city_i = self.cities_data[i]
            for j, value_ij in city_i.distance_km_to_other_cities.items():
                if j in self.set_cities:
                    distance[(i, j)] = int(value_ij)
        return distance

    def __add_var_Z__(self):
        """
        Z_ik: 1 if node i is assigned to hub at node k, else 0 (binary)
        """
        N = self.set_cities
        return self.model.addVars(N, N, lb=0, ub=1, vtype = GRB.BINARY, name="Z")

    def __add_var_Y__(self):
        """
        Y_ikl : flow from i that goes between hubs k and l
        """
        N = self.set_cities
        return self.model.addVars(N, N, N, lb=0, vtype = GRB.CONTINUOUS, name="Y")
    
    def _add_constraint_locate_hubs(self):
        """
        We need to locate p hubs
        """
        N = self.set_cities
        self.model.addConstr(
            gp.quicksum(self.var_Z[i, i] for i in N) == self.par_number_hubs,
            name = "locate p hubs")

    def _add_constraint_assign_hub_to_node(self):
        """
        Every node has to be assigned to a single hub
        """
        N = self.set_cities
        self.model.addConstrs(
            (self.var_Z.sum(i, "*") == 1 for i in N), 
            name = "node to hub assignment")

    def _add_constraint_hub_implication(self):
        """
        Assigning a node to a hub implies the hub was located
        """
        N = self.set_cities
        self.model.addConstrs(
            (self.var_Z[i,k] <= self.var_Z[k,k] for i in N for k in N if i != k), 
            name = "node to hub implies hub")

    def _add_constraint_multicommodity_flow_equation(self):
        """
        Flow conservation from each node at allocated hubs
        """
        N = self.set_cities
        for i in N:
            for k in N:
                self.model.addConstr(
                    (gp.quicksum(self.var_Y[i,k,l] for l in N) - 
                     gp.quicksum(self.var_Y[i,l,k] for l in N) ==
                     self.par_supply[i]*self.var_Z[i,k] - 
                     gp.quicksum(self.par_flow[i,j] * self.var_Z[j,k] for j in N)
                    ),
                    name = f"flow_conservation[{i},{k}]"
                )
                

        
        
        

            






if __name__ == "__main__":
    from dataloader import load_data
    cities_data = load_data()
    problem = USApHMP(cities_data, max_nodes=4, number_hubs=2)

from typing import Dict
from class_city import City
import gurobipy as gp
from gurobipy import GRB

class USApSCT:
    """This class implements the uncapacitated, single allocation p-hub set covering with
    time service constraints, a mixture between:
    - Ernst, A.T., Krishnamoorthy, M. (1996)
    - Tan, P.Z., Kara, B. Y. (2007)
    
    References:
    - Tan, P.Z., Kara, B. Y. (2007). “A hub covering model for cargo delivery systems” 
    Networks 49(1): 28-39, doi: 10.1002/net.20139.
    - Ernst, A.T., Krishnamoorthy, M. (1996). “Efficient algorithms for the uncapacitated 
    ingle allocation p-hub median problem”. Location Science 4(3): 139-154.

    """

    def __init__(self, cities_data: Dict[int, City], max_nodes:int = 81, number_hubs:int = 1, max_time:int = 2000):
        assert 2 <= max_nodes <= 81
        assert 1 <= number_hubs
        assert 1 <= max_time
        self.cities_data = cities_data
        
        # Define sets starting with set_
        self.set_cities = self.__get_set_cities__(max_nodes)

        # Define parameters starting with par_
        self.par_flow = self.__get_city_to_city_parameter__("flow")
        self.par_fixed_hub_cost = self.__get_hub_cost__()
        self.par_distance_km = self.__get_city_to_city_parameter__("distance_km")
        self.par_travel_time_min = self.__get_city_to_city_parameter__("travel_time_min")
        self.par_fixed_link_cost = self.__get_city_to_city_parameter__("fixed_link_cost")
        self.par_supply = self.__get_node_supply__()
        self.par_cost_per_km = 10
        self.par_hub_to_hub_discount_factor = 0.5
        self.par_number_hubs = number_hubs
        self.par_max_arrival_time = max_time
        
        # Create optimization model
        self.instance_name = "USApSCT_n" + {max_nodes} + "_p" + {self.par_number_hubs} + "_t" + {self.par_max_arrival_time}
        self.model = gp.Model('USApSCT')
        
        # Create variables (var_)
        self.var_flow_orig_hub = self._add_var_flow_orig_hub()
        self.var_assign_orig_hub = self._add_var_assign_orig_hub()
        self.var_flow_orig_hubs = self._add_var_flow_orig_hubs()
        self.var_flow_orig_hub_dest = self._add_var_flow_orig_hub_dest()
        self.var_link_hubs = self._add_var_link_hubs()
        self.var_departure_from_hub_to_hub = self._add_var_departure_from_hub_to_hub()
        self.var_departure_from_hub_to_dest = self._add_var_departure_from_hub_to_dest()
        self.var_arrival_to_dest = self._add_var_arrival_to_dest()

        # Create constraints
        self._add_constraint_locate_hubs()
        self._add_constraint_assign_node_to_single_hub()
        self._add_constraint_hub_implication()
        self._add_constraint_flow_from_source()
        self._add_constraint_flow_between_hubs_implies_link()


        # Set objective function
        self.__set_objective_function__()

        self.model.update()

    
    def __get_set_cities__(self, max_nodes: int) -> list:
        return list(cities_data.keys())[0: max_nodes]

    def __get_city_to_city_parameter__(self, param_name: str) -> dict():
        param = dict() # (i,j) -> parameter_ij
        for i in self.set_cities:
            city_i = self.cities_data[i]
            
            city_i_to_j_dict = dict()
            if param_name == "distance_km":
                city_i_to_j_dict = city_i.distance_km_to_other_cities
            elif param_name == "travel_time_min":
                city_i_to_j_dict = city_i.travel_time_min_to_other_cities
            elif param_name == "flow":
                city_i_to_j_dict = city_i.flow_goods_to_other_cities
            elif param_name == "fixed_link_cost":
                city_i_to_j_dict = city_i.fixed_hub_cost
            else:
                raise ValueError(f"param_name {param_name} not recognized")

            for j, value_ij in city_i_to_j_dict.items():
                if j in self.set_cities:
                    param[(i, j)] = value_ij
        return param

    def __get_hub_cost__(self) -> dict():
        return {i: city.fixed_hub_cost for i, city in self.cities_data.items()}

    def __get_node_supply__(self) -> dict:
        flow_origin = dict() # supply of node i
        for i in self.set_cities:
            city_i = self.cities_data[i]
            total_supply = 0
            for j, value_ij in city_i.flow_goods_to_other_cities.items():
                if j in self.set_cities:
                    total_supply += value_ij
            flow_origin[i] = total_supply
        return flow_origin



    def _add_var_flow_orig_hub(self):
        """
        F_ik: Flow from node i to hub in node k (>= 0)
        """
        N = self.set_cities
        return self.model.addVars(N, N, lb=0, vtype = GRB.CONTINUOUS, name="flow_orig_hub")

    def _add_var_assign_orig_hub(self):
        """
        Z_ik: 1 if node i is assigned to hub at node k, else 0 (binary)
        """
        N = self.set_cities
        return self.model.addVars(N, N, lb=0, ub=1, vtype = GRB.BINARY, name="assign_orig_hub")

    def _add_var_flow_orig_hubs(self):
        """
        Y_ikl: flow from i that goes between hubs k and l
        """
        N = self.set_cities
        return self.model.addVars(N, N, N, lb=0, vtype = GRB.CONTINUOUS, name="flow_orig_hubs")
    
    def _add_var_flow_orig_hub_dest(self):
        """X_ilj: flow from i going through hub l to destination j
        """
        N = self.set_cities
        return self.model.addVars(N, N, N, lb=0, vtype = GRB.CONTINUOUS, name="flow_orig_hub_dest")

    def _add_var_link_hubs(self):
        """R_kl: 1 if hubs k and l are connected to transport goods
        """
        N = self.set_cities
        return self.model.addVars(N, N, lb=0, vtype = GRB.BINARY, name="link_hubs")

    def _add_var_departure_from_hub_to_hub(self):
        """D^hat_k: departure time from hub k for transportation to other hubs
        """
        N = self.set_cities
        return self.model.addVars(N, lb=0, ub=self.par_max_arrival_time, vtype = GRB.CONTINUOUS, name="departure_from_hub_to_hub")

    def _add_var_departure_from_hub_to_dest(self):
        """D_k: departure time from hub k for transportation to destinations
        """
        N = self.set_cities
        return self.model.addVars(N, lb=0, ub=self.par_max_arrival_time, vtype = GRB.CONTINUOUS, name="departure_from_hub_to_dest")

    def _add_var_arrival_to_dest(self):
        """A_j: arrival time to destination j
        """
        N = self.set_cities
        return self.model.addVars(N, lb=0, ub=self.par_max_arrival_time, vtype = GRB.CONTINUOUS, name="arrival_to_dest")








    def _add_constraint_locate_hubs(self):
        """
        We need to locate the predefined number of hubs

        Σ_i Z_ii = p
        """
        N = self.set_cities
        self.model.addConstr(
            gp.quicksum(self.var_assign_orig_hub[i, i] for i in N) == self.par_number_hubs,
            name = "locate p hubs")

    def _add_constraint_assign_node_to_single_hub(self):
        """
        Every node has to be assigned to a single hub

        Σ_k Z_ik = 1, Ɐ i ∈ N
        """
        N = self.set_cities
        self.model.addConstrs(
            (self.var_assign_orig_hub.sum(i, "*") == 1 for i in N), 
            name = "node to hub assignment")

    def _add_constraint_hub_implication(self):
        """
        Assigning a node to a hub implies the hub was located

        Z_ik <= Z_kk, Ɐ i ∈ N, k ∈ N: i != k
        """
        N = self.set_cities
        self.model.addConstrs(
            (self.var_assign_orig_hub[i,k] <= self.var_assign_orig_hub[k,k] 
            for i in N for k in N if i != k), 
            name = "node to hub implies hub")

    def _add_constraint_flow_from_source(self):
        """
        The flow from every node goes completely to its hub

        F_ik = Z_ik * s_i    Ɐ i ∈ N, k ∈ N: i != k
        """
        N = self.set_cities
        self.model.addConstrs(
            (self.var_flow_orig_hub[i,k] == self.var_assign_orig_hub[i,k] * self.par_supply[i] 
             for i in N for k in N if i != k),
             name = "flow from source goes to hub")

    def _add_constraint_flow_between_hubs_implies_link(self):
        """
        The flow from every node goes completely to its hub

        Σ_i Y_ikl <= R_kl, Ɐ k ∈ N, l ∈ N: k != l
        """
        N = self.set_cities
        self.model.addConstrs(
            (self.var_flow_orig_hubs.sum()





    
    def __set_objective_function__(self):
        N = self.set_cities
        self.model.setObjective(
            # collection and transfer costs
            gp.quicksum(
                self.par_distance[i,k] * self.var_Z[i,k] * (self.par_unit_cost_collection *self.par_supply[i] + 
                self.par_unit_cost_transfer * self.par_demand[i]) for i in N for k in N
            ) + 
            # distribution costs
            gp.quicksum(
                self.par_unit_cost_distribution * self.par_distance[k,l] * self.var_Y[i,k,l] for i in N for k in N for l in N
            ),
            sense=GRB.MINIMIZE
        )

    def solve(self):
        print("\n\n")
        self.model.optimize()    

    def get_solution(self):

        for index, val in self.model.getAttr('X', self.var_Z).items():
            if val > 0:
                print(f"Z {index}: {val}")

        for index, val in self.model.getAttr('X', self.var_Y).items():
            if val > 0:
                print(f"Y {index}: {val}")



if __name__ == "__main__":
    from dataloader import load_data
    cities_data = load_data()
    problem = USApHMP(cities_data, max_nodes=81, number_hubs=5)
    problem.solve()
    problem.get_solution()

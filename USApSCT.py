from typing import Dict
from class_city import City
import gurobipy as gp
from gurobipy import GRB
import plotter
import utils

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
        self.set_cities = self.__get_set_cities__(max_nodes)                                # N

        # Define parameters starting with par_
        self.par_flow = self.__get_city_to_city_parameter__("flow")                         # w_ij
        self.par_fixed_hub_cost = self.__get_hub_cost__()                                   # h_i
        self.par_distance_km = self.__get_city_to_city_parameter__("distance_km")           # d_ij
        self.par_travel_time_min = self.__get_city_to_city_parameter__("travel_time_min")   # t_ij
        self.par_fixed_link_cost = self.__get_city_to_city_parameter__("fixed_link_cost")   # f_ij
        self.par_supply = self.__get_node_supply__()                                        # s_i
        self.par_unit_cost_per_km_flow = 0.1                                                # c
        self.par_hub_to_hub_discount_factor = 0.5                                           # α
        self.par_number_hubs = number_hubs                                                  # p
        self.par_max_arrival_time = max_time                                                # β
        
        # Create optimization model
        self.instance_name = "USApSCT_n" + str(max_nodes)
        self.instance_name += "_p" + str(self.par_number_hubs)
        self.instance_name += "_t" + str(self.par_max_arrival_time)
        print(f"Instance name: {self.instance_name}")
        self.model = gp.Model('USApSCT')
        
        # Create variables (var_)
        self.var_flow_orig_hub = self._add_var_flow_orig_hub()                              # F_ik
        self.var_assign_orig_hub = self._add_var_assign_orig_hub()                          # Z_ik
        self.var_flow_orig_hub_hub = self._add_var_flow_orig_hub_hub()                      # Y_ikl
        self.var_flow_orig_hub_dest = self._add_var_flow_orig_hub_dest()                    # X_ilj
        self.var_link_hubs = self._add_var_link_hubs()                                      # R_kl
        self.var_departure_from_hub_to_hub = self._add_var_departure_from_hub_to_hub()      # DHH_k
        self.var_departure_from_hub_to_dest = self._add_var_departure_from_hub_to_dest()    # DHD_k
        self.var_arrival_to_dest = self._add_var_arrival_to_dest()                          # A_j

        # Create constraints
        self._add_constraint_locate_hubs()
        self._add_constraint_assign_node_to_single_hub()
        self._add_constraint_hub_implication()
        self._add_constraint_flow_from_source()
        self._add_constraint_flow_between_hubs_implies_link()
        self._add_constraint_conservation_of_flow_at_hub_from_origin()
        self._add_constraint_flow_node_in_hubs_implies_node_to_hub()
        self._add_constraint_flow_origin_destination_covered()
        self._add_constraint_link_between_hubs_implies_hubs()
        self._add_constraint_departure_from_hub_to_hub_waits_customers()
        self._add_constraint_departure_from_hub_to_destination_waits_for_hub_trucks()
        self._add_constraint_latest_arrival_time_to_destination()
        self._add_constraint_latest_arrival_to_dest_within_limit()


        # Set objective function
        # Variables for objective function (var_OF_)
        self.var_OF_location_cost = self._add_var_OF_location_cost()
        self.var_OF_transport_collection_cost = self._add_var_OF_transport_collection_cost()
        self.var_OF_transport_transfer_cost = self._add_var_OF_transport_transfer_cost()
        self.var_OF_transport_distribution_cost = self._add_var_OF_transport_distribution_cost()
        self.var_OF_links_cost = self._add_var_OF_links_cost()
        # Constraints for setting the variables of the OF
        self._add_constraint_OF_location_cost()
        self._add_constraint_OF_transport_collection_cost()
        self._add_constraint_OF_transport_transfer_cost()
        self._add_constraint_OF_transport_distribution_cost()
        self._add_constraint_OF_links_cost()
        # Setting the objecting function in the model
        self.__set_objective_function__()

        self.model.update()

        # print(f"\n\nCONSTRAINTS:")
        # for i, con in enumerate(self.model.getConstrs()):
        #     print(f"\n{con}:\n{self.model.getRow(con)} {con.Sense} {con.RHS}")            
        

    
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
                city_i_to_j_dict = city_i.fixed_link_cost_to_other_cities
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


    #########################################
    #                                       #
    #              VARIABLES                #
    #                                       #
    #########################################

    def _add_var_flow_orig_hub(self):
        """
        F_ik: Flow from node i to hub in node k (>= 0)
        """
        N = self.set_cities
        return self.model.addVars(N, N, lb=0, vtype = GRB.CONTINUOUS, name="F")

    def _add_var_assign_orig_hub(self):
        """
        Z_ik: 1 if node i is assigned to hub at node k, else 0 (binary)
        """
        N = self.set_cities
        return self.model.addVars(N, N, lb=0, ub=1, vtype = GRB.BINARY, name="Z")

    def _add_var_flow_orig_hub_hub(self):
        """
        * Y_ikl: flow from i that goes between hubs k and l
        * Ɐ i ∈ N, Ɐ k ∈ N, Ɐ l ∈ N: k != l & i != l
        """
        N = self.set_cities
        N3 = set()
        for i in N:
            for k in N:
                for l in N:
                    if k != l and i != l:
                        N3.add((i,k,l))
        return self.model.addVars(N3, lb=0, vtype = GRB.CONTINUOUS, name="Y")
    
    def _add_var_flow_orig_hub_dest(self):
        """
        * X_ilj: flow from i going through hub l to destination j
        * Ɐ i ∈ N, l ∈ N, j ∈ N: i != j
        """
        N = self.set_cities
        N3 = set()
        for i in N:
            for l in N:
                for j in N:
                    if i != j:
                        N3.add((i,l,j))
        return self.model.addVars(N3, lb=0, vtype = GRB.CONTINUOUS, name="X")

    def _add_var_link_hubs(self):
        """R_kl: 1 if hubs k and l are connected to transport goods
        """
        N = self.set_cities
        return self.model.addVars(N, N, lb=0, vtype = GRB.BINARY, name="R")

    def _add_var_departure_from_hub_to_hub(self):
        """DHH_k: departure time from hub k for transportation to other hubs
        """
        N = self.set_cities
        return self.model.addVars(N, lb=0, ub=self.par_max_arrival_time, vtype = GRB.CONTINUOUS, name="DHH")

    def _add_var_departure_from_hub_to_dest(self):
        """DHD_k: departure time from hub k for transportation to destinations
        """
        N = self.set_cities
        return self.model.addVars(N, lb=0, ub=self.par_max_arrival_time, vtype = GRB.CONTINUOUS, name="DHD")

    def _add_var_arrival_to_dest(self):
        """A_j: arrival time to destination j
        """
        N = self.set_cities
        return self.model.addVars(N, lb=0, ub=self.par_max_arrival_time, vtype = GRB.CONTINUOUS, name="A")


    #########################################
    #                                       #
    #            CONSTRAINTS                #
    #                                       #
    #########################################


    def _add_constraint_locate_hubs(self):
        """
        We need to locate the predefined number of hubs

        * Σ_i Z_ii = p
        """
        N = self.set_cities
        self.model.addConstr(
            gp.quicksum(self.var_assign_orig_hub[i, i] for i in N) == self.par_number_hubs,
            name = "locate p hubs")

    def _add_constraint_assign_node_to_single_hub(self):
        """
        Every node has to be assigned to a single hub

        * Σ_k Z_ik = 1
        * Ɐ i ∈ N
        """
        N = self.set_cities
        self.model.addConstrs(
            (self.var_assign_orig_hub.sum(i, "*") == 1 for i in N), 
            name = "node to hub assignment")

    def _add_constraint_hub_implication(self):
        """
        Assigning a node to a hub implies the hub was located

        * Z_ik <= Z_kk
        * Ɐ i ∈ N, k ∈ N: i != k
        """
        N = self.set_cities
        self.model.addConstrs(
            (self.var_assign_orig_hub[i,k] <= self.var_assign_orig_hub[k,k] 
            for i in N for k in N if i != k), 
            name = "node to hub implies hub")

    def _add_constraint_flow_from_source(self):
        """
        The flow from every node goes completely to its hub

        * F_ik = Z_ik * s_i
        * Ɐ i ∈ N, k ∈ N 
        """
        N = self.set_cities
        self.model.addConstrs(
            (self.var_flow_orig_hub[i,k] == self.var_assign_orig_hub[i,k] * self.par_supply[i] 
             for i in N for k in N),
             name = "flow from source goes to hub")

    def _add_constraint_flow_between_hubs_implies_link(self):
        """
        The flow from every node goes completely to its hub

        * Y_ikl <= R_kl * s_i
        * Ɐ i ∈ N, Ɐ k ∈ N, l ∈ N: k != l & i != l
        """
        N = self.set_cities
        self.model.addConstrs(
            (self.var_flow_orig_hub_hub[i,k,l]
             <= 
             self.var_link_hubs[k,l] * self.par_supply[i]
             for i in N for k in N for l in N if (k != l) and (l != i)),
             name = "flow between hubs implies link between them")

    def _add_constraint_conservation_of_flow_at_hub_from_origin(self):
        """
        Conservation flow from commodity (origin) i at hub k: flow in = flow out
        * F_ik + Σ_l (k != l and i != l) Y_ilk  = Σ_l (k != l and i != l) Y_ikl + Σ_j X_ikj
        * Ɐ i ∈ N, Ɐ k ∈ N
        """
        N = self.set_cities
        self.model.addConstrs(
            (self.var_flow_orig_hub[i,k] +
             gp.quicksum(self.var_flow_orig_hub_hub[i,l,k] for l in N if (i != k) and (l != k))
             == 
             gp.quicksum(self.var_flow_orig_hub_hub[i,k,l] for l in N if (k != l) and (l != i)) + 
             gp.quicksum(self.var_flow_orig_hub_dest[i,k,j] for j in N if i != j)
             for i in N for k in N),
             name = "multicomodity flow conservation at each hub")

    def _add_constraint_flow_node_in_hubs_implies_node_to_hub(self):
        """
        If there is flow from origin i between hubs (k,l) then origin i must have been
        assigned to hub k

        * Y_ikl <= Z_ik * s_i
        * Ɐ i ∈ N, Ɐ k ∈ N, l ∈ N: k != l & l != i
        """
        N = self.set_cities
        self.model.addConstrs(
            (self.var_flow_orig_hub_hub[i,k,l] <= 
             self.var_assign_orig_hub[i,k] * self.par_supply[i]
             for i in N for k in N for l in N if (k != l) and (l != i)
             ),
             name = "flow from origin in hubs implies origin is assigned to hub"
        )

    def _add_constraint_flow_origin_destination_covered(self):
        """
        The flow from origin i to destination j is covered

        * X_ilj = Z_jl * w_ij
        * Ɐ i ∈ N, Ɐ l ∈ N, j ∈ N: i != j
        """
        N = self.set_cities
        self.model.addConstrs(
            (self.var_flow_orig_hub_dest[i,l,j] == 
             self.var_assign_orig_hub[j,l] * self.par_flow[i,j]
             for i in N for l in N for j in N if i != j
             ),
             name = "flow from origin to destination is covered"
        )

    def _add_constraint_link_between_hubs_implies_hubs(self):
        """
        If there is a link between hubs k and l, then there is a hub in k and in l

        * 2*R_kl <= Z_kk + Z_ll
        * Ɐ k ∈ N, Ɐ l ∈ N: k != l
        """
        N = self.set_cities
        self.model.addConstrs(
            (2*self.var_link_hubs[k,l] <=
             self.var_assign_orig_hub[k,k] + self.var_assign_orig_hub[l,l]
             for k in N for l in N if k != l
             ),
             name = "link between hubs implies hubs"
        )

    def _add_constraint_departure_from_hub_to_hub_waits_customers(self):
        """
        A vehicle departuring from hub k to other hubs must wait cargo from
        its customers (nodes assigned to k)

        * 2*DHH_k >= t_ik * Z_ik
        * Ɐ i ∈ N, Ɐ k ∈ N: i != k
        """
        N = self.set_cities
        self.model.addConstrs(
            (2*self.var_departure_from_hub_to_hub[k] >=
             self.par_travel_time_min[i,k] * 
             self.var_assign_orig_hub[i,k]
             for i in N for k in N if i != k
             ),
             name = "departure from hub to hub waits for its customers"
        )
    
    def _add_constraint_departure_from_hub_to_destination_waits_for_hub_trucks(self):
        """
        The departuring time from hub l to its customers DHC_l must
        wait for incoming cargo from other hubs k DHH_k to consolidate cargo to customers

        * Original non-linear: DHD_l >= (DHH_k + t_kl) * R_kl
        * Linearized: DHD_l >= DHH_k + t_kl - M(1 - R_kl)
        * Where M is a sufficient big number, in this case: M = β (max delivery time)
        * Ɐ k ∈ N, Ɐ l ∈ N
        """
        N = self.set_cities
        self.model.addConstrs(
            (self.var_departure_from_hub_to_dest[l] >=
             self.var_departure_from_hub_to_hub[k] + 
             self.par_travel_time_min[k,l] - 
             self.par_max_arrival_time*(1 - self.var_link_hubs[k,l])
             for k in N for l in N
             ),
             name = "departure from hub to destination waits for hub to hub trucks"
        )
    
    def _add_constraint_latest_arrival_time_to_destination(self):
        """
        The latest arrival time to destination j is subject to its deliveries from hub l

        * Original non-linear: A_j >= (DHD_l + t_lj) * Z_jl
        * Linearized: A_j >= DHD_l + t_lj - M(1 - Z_jl)
        * Where M is a sufficient big number, in this case: M = β (max delivery time)
        * Ɐ j ∈ N, Ɐ l ∈ N
        """
        N = self.set_cities
        self.model.addConstrs(
            (self.var_arrival_to_dest[j] >=
             self.var_departure_from_hub_to_dest[l] + 
             self.par_travel_time_min[l,j] - 
             self.par_max_arrival_time*(1 - self.var_assign_orig_hub[j,l])
             for j in N for l in N
             ),
             name = "departure from hub to destination waits for hub to hub trucks"
        )

    def _add_constraint_latest_arrival_to_dest_within_limit(self):
        """
        The latest arrival time to destination j must be within the standard

        * A_j <= β
        * Ɐ j ∈ N
        """
        N = self.set_cities
        self.model.addConstrs(
            (self.var_arrival_to_dest[j] <= self.par_max_arrival_time
             for j in N
             ),
             name = "max arrival time within standard"
        )

    #########################################
    #                                       #
    #          OBJECTIVE FUNCTION           #
    #                                       #
    #########################################

    def _add_var_OF_location_cost(self):
        return self.model.addVar(lb=0, vtype = GRB.CONTINUOUS, name="OF_location_cost")
    
    def _add_constraint_OF_location_cost(self):
        """The location cost of the hubs

        * Σ_i Z_ii * h_i
        """
        N = self.set_cities
        self.model.addConstr(
            self.var_OF_location_cost == 
            gp.quicksum(
                self.var_assign_orig_hub[i,i] * self.par_fixed_hub_cost[i]
                for i in N
            ),
            name = "OF_location_cost"
        )

    def _add_var_OF_transport_collection_cost(self):
        return self.model.addVar(lb=0, vtype = GRB.CONTINUOUS, name="OF_transport_collection_cost")

    def _add_constraint_OF_transport_collection_cost(self):
        """The transportation cost of collection activities

        * Σ_i Σ_k F_ik * d_ik * c
        """
        N = self.set_cities
        self.model.addConstr(
            self.var_OF_transport_collection_cost == 
            gp.quicksum(
                self.var_flow_orig_hub[i,k] * self.par_distance_km[i,k] * self.par_unit_cost_per_km_flow
                for i in N for k in N
            ),
            name = "OF_transport_collection_cost"
        )

    def _add_var_OF_transport_transfer_cost(self):
        return self.model.addVar(lb=0, vtype = GRB.CONTINUOUS, name="OF_transport_transfer_cost")

    def _add_constraint_OF_transport_transfer_cost(self):
        """The transportation cost of transfer activities

        * Σ_i Σ_k Σ_l (k != l & l != i) Y_ikl * d_kl * c
        """
        N = self.set_cities
        self.model.addConstr(
            self.var_OF_transport_transfer_cost == 
            gp.quicksum(
                self.var_flow_orig_hub_hub[i,k,l] * self.par_distance_km[k,l] * self.par_unit_cost_per_km_flow
                for i in N for k in N for l in N if (k != l) and (l != i)
            ),
            name = "OF_transport_transfer_cost"
        )

    def _add_var_OF_transport_distribution_cost(self):
        return self.model.addVar(lb=0, vtype = GRB.CONTINUOUS, name="OF_transport_distribution_cost")
    
    def _add_constraint_OF_transport_distribution_cost(self):
        """The transportation cost of distribution activities

        * Σ_i Σ_l Σ_j X_ilj * d_lj * c
        """
        N = self.set_cities
        self.model.addConstr(
            self.var_OF_transport_distribution_cost == 
            gp.quicksum(
                self.var_flow_orig_hub_dest[i,l,j] * self.par_distance_km[l,j] * self.par_unit_cost_per_km_flow
                for i in N for l in N for j in N if i != j
            ),
            name = "OF_transport_distribution_cost"
        )

    def _add_var_OF_links_cost(self):
        return self.model.addVar(lb=0, vtype = GRB.CONTINUOUS, name="OF_links_cost")

    def _add_constraint_OF_links_cost(self):
        """The fixed cost of stablishing links between nodes

        * Σ_i Σ_k Z_ik * f_ik + Σ_k Σ_l R_kl * f_kl
        """
        N = self.set_cities
        self.model.addConstr(
            self.var_OF_links_cost == 
            gp.quicksum(
                self.var_assign_orig_hub[i,k] * self.par_fixed_link_cost[i,k]
                for i in N for k in N
            ) + 
            gp.quicksum(
                self.var_link_hubs[k,l] * self.par_fixed_link_cost[k,l]
                for k in N for l in N
            ),
            name = "OF_links_cost"
        )
    
    def __set_objective_function__(self):
        self.model.setObjective(
            self.var_OF_links_cost +
            self.var_OF_location_cost + 
            self.var_OF_transport_collection_cost + 
            self.var_OF_transport_distribution_cost + 
            self.var_OF_transport_transfer_cost,
            sense=GRB.MINIMIZE
        )

    def solve(self):
        print("\n\n\nSOLVING...\n")
        self.model.optimize()

        if self.model.status == GRB.INFEASIBLE:
            print(f"\n\nModel infeasable, computing IIS:")
            self.model.computeIIS()
            report = self.instance_name + "_IIS.ilp"
            self.model.write(report)
            print(f"IIS report written in {report}")

    def save_solution(self):
        """Saves the variables with values != 0 in a json file
        """
        
        if self.model.Status != GRB.OPTIMAL:
            return
        
        solution = dict()
        
        def remove_0_values(d: dict()) -> dict():
            return { str(k): v for (k,v) in d.items() if v != 0}
        
        # Objective function value
        solution["OF_val"] = self.model.ObjVal

        # Append objective function components
        solution['OF'] = {}
        solution['OF']['location_cost'] = self.var_OF_location_cost.X
        solution['OF']['transport_collection_cost'] = self.var_OF_transport_collection_cost.X
        solution['OF']['transport_transfer_cost'] = self.var_OF_transport_transfer_cost.X
        solution['OF']['transport_distribution_cost'] = self.var_OF_transport_distribution_cost.X
        solution['OF']['links_cost'] = self.var_OF_links_cost.X

        # Append variables
        solution['variables'] = {}
        solution['variables']['F'] = remove_0_values(self.model.getAttr('X', self.var_flow_orig_hub))
        solution['variables']['Z'] = remove_0_values(self.model.getAttr('X', self.var_assign_orig_hub))
        solution['variables']['Y'] = remove_0_values(self.model.getAttr('X', self.var_flow_orig_hub_hub))
        solution['variables']['X'] = remove_0_values(self.model.getAttr('X', self.var_flow_orig_hub_dest))
        solution['variables']['R'] = remove_0_values(self.model.getAttr('X', self.var_link_hubs))
        solution['variables']['DHH'] = remove_0_values(self.model.getAttr('X', self.var_departure_from_hub_to_hub))
        solution['variables']['DHD'] = remove_0_values(self.model.getAttr('X', self.var_departure_from_hub_to_dest))
        solution['variables']['A'] = remove_0_values(self.model.getAttr('X', self.var_arrival_to_dest))

        utils.save_solution(solution, file_name=self.instance_name)

    def load_solution(self):
        """Loads the values of the variables in a saved solution
        """        
        return utils.load_solution(self.instance_name)

    def plot_solution(self):
        solution = self.load_solution()
        
        # Input for plotting
        hubs_ids = set()
        origin_hub_flow_collection = solution["variables"]["F"]
        origin_hub_hub_flow_transfer = solution["variables"]["Y"]
        origin_hub_destination_flow_distribution = solution["variables"]["X"]

        # Getting hubs_ids
        Z = solution["variables"]["Z"]
        for (i,j), _ in Z.items():
            if i == j:
                hubs_ids.add(j)
        
        # Now we can plot
        plotter.plot_map(cities_data = self.cities_data,
                         cities_considered = self.set_cities,
                         hubs_ids=hubs_ids,
                         collection=origin_hub_flow_collection,
                         transfer=origin_hub_hub_flow_transfer,
                         distribution=origin_hub_destination_flow_distribution,
                         plot_name=self.instance_name
                         )




if __name__ == "__main__":
    from dataloader import load_data
    cities_data = load_data()
    problem = USApSCT(cities_data, max_nodes=2, number_hubs=2, max_time=3000)
    problem.solve()
    problem.save_solution()
    problem.plot_solution()
